"""Homology tier (miniprot) tests. The GFF3 parser is a pure function of the GFF text + genomic
sequence, so gene-model recovery is verified off captured fixtures with no WSL. Input guards and
the provenance seal (the reference protein IS the evidence and must seal the run) are checked at
the engine/provenance layer, also WSL-free."""
import hashlib
import pytest
import engine
from engine import BadRequest
from teagle_core import wsl, sequtil, provenance

_CODON = {'M':'ATG','K':'AAA','T':'ACA','A':'GCA','Y':'TAT','I':'ATT','Q':'CAA','R':'CGT','G':'GGT',
          'S':'TCA','L':'CTT','V':'GTT','P':'CCA','F':'TTT','N':'AAT','D':'GAT','E':'GAA','H':'CAT',
          'W':'TGG','C':'TGT'}
_PROT = "MKTAYIAKQRGSLVPFNDEHWCMKTAYIAK"
_DNA = "".join(_CODON[a] for a in _PROT)
_INTRON = "GT" + "GAGT" + "A"*18 + "CTGAC" + "T"*14 + "TTTTTTCCCCAG"       # canonical GT..AG, 55 bp
_G_PLUS = _DNA[:45] + _INTRON + _DNA[45:]
_rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
_G_MINUS = _rc(_G_PLUS)


def _row(*c):
    return "\t".join(str(x) for x in c)


def _gff(contig, prot, strand, cds, extra_mrna="", paf_tags=("fs:i:0", "st:i:0")):
    lines = ["##gff-version 3",
             _row("##PAF", prot, len(_PROT), 0, len(_PROT), strand, contig, 145, 0, 145, 90, 90, 0, *paf_tags)]
    mrna_attr = f"ID=MP1;Rank=1;Identity=1.0000;Positive=1.0000;{extra_mrna}Target={prot} 1 {len(_PROT)}"
    lines.append(_row(contig, "miniprot", "mRNA", 1, 145, 164, strand, ".", mrna_attr))
    for (s, e, tgt) in cds:
        lines.append(_row(contig, "miniprot", "CDS", s, e, 90, strand, 0,
                          f"Parent=MP1;Rank=1;Identity=1.0000;Target={prot} {tgt}"))
    return "\n".join(lines)


# ---------- parser: gene-model recovery ----------
def test_parser_plus_strand_intron():
    gff = _gff("c", "ref", "+", [(1, 45, "1 15"), (101, 145, "16 30")])
    h = wsl._parse_miniprot_gff(gff, _G_PLUS)
    assert len(h) == 1
    x = h[0]
    assert x["strand"] == "+" and x["counts"] == {"exons": 2, "introns": 1, "canonical_introns": 1}
    assert x["introns"][0]["donor"] == "GT" and x["introns"][0]["acceptor"] == "AG"
    assert x["introns"][0]["length"] == 55 and x["protein_coverage"] == 1.0


def test_parser_minus_strand_reads_motif_on_transcribed_strand():
    # on the minus strand miniprot lists CDS in reverse genomic order; the intron gap is still
    # found by genomic-coordinate sorting and the GT..AG motif read from the reverse strand.
    gff = _gff("c", "ref", "-", [(101, 145, "1 15"), (1, 45, "16 30")])
    x = wsl._parse_miniprot_gff(gff, _G_MINUS)[0]
    assert x["strand"] == "-"
    assert x["counts"]["introns"] == 1 and x["introns"][0]["canonical"] is True
    assert x["introns"][0]["donor"] == "GT" and x["introns"][0]["acceptor"] == "AG"


def test_parser_frameshift_is_not_an_intron():
    # miniprot keeps a frameshift inside ONE CDS (Frameshift= attr), so it must never be read as an intron.
    gff = "\n".join([
        "##gff-version 3",
        _row("##PAF", "p", 60, 0, 60, "+", "c", 181, 0, 181, 180, 181, 0, "fs:i:1", "st:i:0", "cg:Z:29M1F31M"),
        _row("c", "miniprot", "mRNA", 1, 181, 305, "+", ".", "ID=MP1;Rank=1;Identity=0.9945;Frameshift=1;Target=p 1 60"),
        _row("c", "miniprot", "CDS", 1, 181, 305, "+", 0, "Parent=MP1;Rank=1;Identity=0.9945;Frameshift=1;Target=p 1 60"),
    ])
    x = wsl._parse_miniprot_gff(gff, "A" * 181)[0]
    assert x["counts"]["introns"] == 0 and x["counts"]["exons"] == 1
    assert x["frameshifts"] == 1


def test_parser_inframe_stops_from_stopcodon_attr_without_paf():
    # no ##PAF line (e.g. a future format): the in-frame-stop count must fall back to the StopCodon= mRNA attribute
    gff = "\n".join([
        "##gff-version 3",
        _row("c", "miniprot", "mRNA", 1, 90, 200, "+", ".", "ID=MP1;Rank=1;Identity=0.9;StopCodon=2;Target=p 1 30"),
        _row("c", "miniprot", "CDS", 1, 90, 200, "+", 0, "Parent=MP1;Target=p 1 30"),
    ])
    assert wsl._parse_miniprot_gff(gff, "A" * 90)[0]["inframe_stops"] == 2


def test_parser_flags_noncanonical_intron():
    # same structure but the genomic gap boundaries are not GT..AG -> intron reported, canonical False.
    g = _DNA[:45] + ("CC" + _INTRON[2:-2] + "TT") + _DNA[45:]
    gff = _gff("c", "ref", "+", [(1, 45, "1 15"), (101, 145, "16 30")])
    x = wsl._parse_miniprot_gff(gff, g)[0]
    assert x["counts"]["introns"] == 1 and x["counts"]["canonical_introns"] == 0
    assert x["introns"][0]["canonical"] is False


def test_parser_ranks_multiple_hits_by_score():
    a = _row("c", "miniprot", "mRNA", 1, 90, 50, "+", ".", "ID=MP1;Rank=1;Identity=0.80;Target=pa 1 30")
    a_cds = _row("c", "miniprot", "CDS", 1, 90, 50, "+", 0, "Parent=MP1;Target=pa 1 30")
    b = _row("c", "miniprot", "mRNA", 1, 90, 300, "+", ".", "ID=MP2;Rank=1;Identity=0.99;Target=pb 1 30")
    b_cds = _row("c", "miniprot", "CDS", 1, 90, 300, "+", 0, "Parent=MP2;Target=pb 1 30")
    gff = "\n".join(["##gff-version 3", a, a_cds, b, b_cds])
    h = wsl._parse_miniprot_gff(gff, "A" * 90, max_hits=5)
    assert [x["protein"] for x in h] == ["pb", "pa"]        # higher score first


def test_parser_respects_max_hits():
    rows = ["##gff-version 3"]
    for i in range(6):
        rows.append(_row("c", "miniprot", "mRNA", 1, 90, 100 + i, "+", ".", f"ID=MP{i};Target=p{i} 1 30"))
        rows.append(_row("c", "miniprot", "CDS", 1, 90, 100 + i, "+", 0, f"Parent=MP{i};Target=p{i} 1 30"))
    assert len(wsl._parse_miniprot_gff("\n".join(rows), "A" * 90, max_hits=3)) == 3


# ---------- engine input guards (raise before touching WSL) ----------
def test_run_miniprot_rejects_missing_sequence():
    with pytest.raises(BadRequest):
        engine.run_miniprot({"protein": ">p\nMKTAYIAKQR"})


def test_run_miniprot_rejects_missing_protein():
    with pytest.raises(BadRequest):
        engine.run_miniprot({"sequence": ">g\n" + _DNA})


def test_run_miniprot_rejects_nucleotide_as_protein():
    with pytest.raises(BadRequest):
        engine.run_miniprot({"sequence": ">g\n" + _DNA, "protein": ">bad\n" + _DNA})


def test_run_miniprot_rejects_nonstring_protein():
    with pytest.raises(BadRequest):
        engine.run_miniprot({"sequence": ">g\n" + _DNA, "protein": 123})


def test_require_aa_accepts_real_protein():
    engine._require_aa("MKTAYIAKQRGSLVPFNDEHW")            # no raise


# ---------- seal: the reference protein must seal the run ----------
def _seal(protein_sha):
    return provenance.build_manifest("homology", _G_PLUS, "contig",
                                     {"tool": "miniprot", "protein_sha256": protein_sha})["manifestSha256"]


def test_homology_seal_distinguishes_reference_protein():
    # different proteins on the same genomic sequence must NOT collide (the transcript-unhashed
    # class of bug present in the splice tier must not be inherited here).
    assert _seal(hashlib.sha256(b"PROT_A").hexdigest()) != _seal(hashlib.sha256(b"PROT_B").hexdigest())


def test_homology_seal_deterministic():
    s = hashlib.sha256(b"PROT_A").hexdigest()
    assert _seal(s) == _seal(s)


def _homology_protein_sha(precs):
    ptext = "\n".join(f">{pid}\n{pseq}" for pid, pseq in precs)     # matches engine.run_miniprot
    return hashlib.sha256(ptext.encode()).hexdigest()


def test_homology_seal_distinguishes_multiprotein_split():
    # two different record splits that concatenate to the same residues must NOT collide (headers seal them)
    a = [("gag", "MK"), ("pol", "TAYIAKQR")]
    b = [("gag", "MKTAY"), ("pol", "IAKQR")]
    assert _seal(_homology_protein_sha(a)) != _seal(_homology_protein_sha(b))


# ---------- protein parser keeps selenocysteine ----------
def test_parse_protein_keeps_U_not_converted_to_T():
    recs = sequtil.parse_protein(">p\nMKUAYW")             # U = selenocysteine, must survive
    assert recs[0][1] == "MKUAYW"


def test_parse_fasta_still_normalizes_rna_u_to_t():
    recs = sequtil.parse_fasta(">n\nACGUACGU")             # nucleotide path unchanged
    assert recs[0][1] == "ACGTACGT"
