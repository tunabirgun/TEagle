"""Coordinate display, primer-table annotation, oligo order sheets and on-screen-only table columns.

Pure helpers come first so they run without a window; the two widget subclasses at the end only change what a
table or viewer SHOWS or COPIES, never what a file receives."""
from __future__ import annotations
import re

from PySide6.QtGui import QKeySequence
from PySide6.QtWidgets import QApplication

import figures
from widgets import DataTable, GenomePanel
from teagle_core import primers as primers_mod

# ---------------------------------------------------------------------------------------------- coordinates
_REC = {0: "Location in the sequence, 0-based half-open [start, end): the first base is 0, the end is not included, "
           "length = end − start. Add 1 to the start (and keep the end) for a 1-based, inclusive position, as "
           "GenBank, SnapGene and IGV count.",
        1: "Location in the sequence, 1-based, inclusive [start, end]: the GenBank, GFF3 and SnapGene convention. "
           "Display only; exported tables stay 0-based, half-open."}
# Per-table tooltip overrides. The 'Coords (1-based)' key is also the isPcr genome-scan header in the shared
# glossary, so every record table passes its own meaning for it.
COORD_GLOSS = {
    0: {"Coords (0-based)": _REC[0],
        "Start (0-based)": "First base of the ORF, 0-based (the first base of the record is 0).",
        "End (exclusive)": "Position just past the ORF's last base (0-based, half-open), so End − Start is its length in bp.",
        "nt (0-based)": "Where the domain hit lies on the record, 0-based half-open [start, end).",
        "Intron span (0-based)": "Intron location on the aligned sequence, 0-based half-open [start, end)."},
    1: {"Coords (1-based)": _REC[1],
        "Start (1-based)": "First base of the ORF, 1-based (the first base of the record is 1).",
        "End (inclusive)": "Last base of the ORF, 1-based, inclusive.",
        "nt (1-based)": "Where the domain hit lies on the record, 1-based, inclusive [start, end]. Display only.",
        "Intron span (1-based)": "Intron location on the aligned sequence, 1-based, inclusive. Display only."},
}
COORD_GLOSS_ALL = {
    "aa (1-based)": "Residue range of the domain within its translated ORF, as HMMER reports it: 1-based, inclusive. "
                    "Not changed by the Coords setting.",
    "Consensus": "Aligned span on the Dfam consensus, as RepeatMasker reports it: 1-based, inclusive, then "
                           "the consensus length. Not changed by the Coords setting.",
}
BASES = {0: "0-based, half-open [start, end) — TEagle tables, BED, Python",
         1: "1-based, inclusive [start, end] — GenBank, GFF3, SnapGene, IGV"}


def span(s, e, base=0):
    """'s–e' in the chosen display convention. The stored interval is 0-based half-open: a 1-based inclusive start is
    s + 1 and the end is unchanged. An empty interval reads 's^s+1' (between two bases) in 1-based mode."""
    if base and e <= s:
        return f"{s}^{s + 1}"
    return f"{s + base}–{e}"


def col(name, base):
    return f"{name} ({base}-based)"


def gloss(base):
    return {**COORD_GLOSS[base], **COORD_GLOSS_ALL}


# ----------------------------------------------------------------------------------------- element geometry
def _pairs(rec):
    return [e for e in rec.get("structural") or []
            if not e.get("advisory") and str(e.get("type", "")).startswith(("LTR", "TIR"))
            and e.get("five_prime") and e.get("three_prime")]


def element_segments(rec, word=None):
    """([(start, end, name)], pair row) for the accepted terminal-repeat pair with the widest extent, in half-open
    record coordinates: 5′ flank, 5′ repeat, internal, 3′ repeat, 3′ flank; empty parts are dropped. ([], None) when
    the record has no accepted pair."""
    rows = _pairs(rec)
    if not rows:
        return [], None
    e = max(rows, key=lambda r: r["three_prime"][1] - r["five_prime"][0])
    word = word or str(e["type"]).split(" ")[0]
    n = (rec.get("composition") or {}).get("length") or len(rec.get("seq") or "")
    (a0, a1), (b0, b1) = e["five_prime"], e["three_prime"]
    parts = [(0, a0, "5′ flank"), (a0, a1, f"5′ {word}"), (a1, b0, "internal"), (b0, b1, f"3′ {word}"), (b1, n, "3′ flank")]
    return [p for p in parts if p[1] > p[0]], e


def insertion_bounds(rec):
    """{start, end, up, down, kind, identity} of the element the accepted terminal-repeat pair delimits, or None."""
    rows = [e for e in _pairs(rec) if e.get("element_span")]
    if not rows:
        return None
    e = max(rows, key=lambda r: r["element_span"][1] - r["element_span"][0])
    s, t = e["element_span"]
    n = (rec.get("composition") or {}).get("length") or len(rec.get("seq") or "")
    return {"start": s, "end": t, "up": s, "down": n - t, "length": n,
            "kind": str(e["type"]).split(" ")[0], "identity": e.get("identity")}


def interior(rec):
    """[start, end) between the two terminal repeats, or None without an accepted pair."""
    _seg, e = element_segments(rec)
    return [e["five_prime"][1], e["three_prime"][0]] if e else None


def domain_choices(rec):
    """Domain hits for the assay picker, strongest HMMER score first: (index in rec['domains'], hit)."""
    doms = list(enumerate(rec.get("domains") or []))
    return sorted(doms, key=lambda x: -(x[1].get("score") or 0))


def _inside_domain(a, b, doms):
    return next((d.get("domain") for d in doms or [] if d.get("nt") and d["nt"][0] <= a and b <= d["nt"][1]), None)


def footprint_place(a, b, segments, doms=None):
    """Where [a, b) lies: the segments it touches joined by '/', then '(DOM)' when it sits wholly inside a domain hit."""
    names = [n for s, e, n in segments if a < e and s < b]
    dom = _inside_domain(a, b, doms)
    return "/".join(names) + (f" ({dom})" if dom and names else "")


def primer_core(c, side):
    """The 3′-terminal MIN_PRIMER_SIZE bases of a primer, half-open on the designed sequence (the shortest 3′ core an
    in-silico PCR accepts as a binding site)."""
    k = primers_mod.MIN_PRIMER_SIZE
    if side == "F":
        a, b = c["left_pos"]
        return [max(a, b - k), b]
    a, b = c["right_pos"]
    return [a, min(b, a + k)]


def pair_location(c, rec, bound, base=0, template_label=None):
    """The Location cell: product span, then where each primer sits. A pair designed on another sequence reads
    'a–b on <label>' and names no segment, because its coordinates index that sequence, not this record."""
    span_ = c.get("intended_product_span") or ([c["left_pos"][0], c["right_pos"][1]]
                                               if c.get("left_pos") and c.get("right_pos") else None)
    if not span_:
        return ""
    text = span(span_[0], span_[1], base)
    if not bound:
        return f"{text} on {template_label or 'a different record'}"
    segments, _pair = element_segments(rec)
    doms = rec.get("domains") or []
    if not (c.get("left_pos") and c.get("right_pos")):
        return text
    if segments:
        f_, r_ = (footprint_place(*c["left_pos"], segments, doms), footprint_place(*c["right_pos"], segments, doms))
        return f"{text} · {f_} → {r_}" if f_ != r_ else f"{text} · {f_}"
    dom = _inside_domain(*c["left_pos"], doms)
    return f"{text} · ({dom})" if dom and dom == _inside_domain(*c["right_pos"], doms) else text


def _hits(a, b):
    return a[0] < b[1] and b[0] < a[1]


def pair_notes(cands, rec, bound, word=None):
    """[(short, tooltip)] per candidate. Three flags, none of which changes a design: a primer whose 3′ core lies
    inside one terminal repeat (it may also bind the other copy and solo copies), two pairs whose sites overlap at
    both ends (nearly the same product), and a primer sequence shared by two pairs (order it once)."""
    segments, pair = element_segments(rec, word) if (bound and rec) else ([], None)
    k = primers_mod.MIN_PRIMER_SIZE
    out = []
    for i, c in enumerate(cands):
        short, tips = [], []
        if pair and c.get("left_pos") and c.get("right_pos"):
            wd = str(pair["type"]).split(" ")[0] if not word else word
            tir = wd == "TIR"
            ident = f" (copies {pair['identity']}% identical)" if pair.get("identity") is not None else ""
            for side, nm in (("F", "forward"), ("R", "reverse")):
                a, b = primer_core(c, side)
                hit = next((n for s, e, n in segments if "flank" not in n and n != "internal" and s <= a and b <= e), None)
                if not hit:
                    continue
                other = "3′" if hit.startswith("5′") else "5′"
                if tir:
                    short.append(f"{side} in {hit}: may prime from the inverted {other} {wd} (single-primer product)")
                    tips.append(f"The {nm} primer's 3′-terminal {k} nt lie inside the {hit}{ident}, which is the inverted "
                                f"counterpart of the {other} {wd}, so it may also prime from that arm and give a single-primer product.")
                else:
                    short.append(f"{side} in {hit}: may also bind {other} {wd} / solo {wd}s")
                    tips.append(f"The {nm} primer's 3′-terminal {k} nt lie inside the {hit}{ident}, so it may also bind the "
                                f"{other} {wd} and solo {wd}s elsewhere in a genome.")
        near, same = [], {"F": [], "R": []}
        for j, d in enumerate(cands):
            if j == i or not (c.get("left_pos") and d.get("left_pos") and c.get("right_pos") and d.get("right_pos")):
                continue
            if _hits(c["left_pos"], d["left_pos"]) and _hits(c["right_pos"], d["right_pos"]):
                near.append(d.get("id"))
            else:
                for side, key in (("F", "left_seq"), ("R", "right_seq")):
                    if c.get(key) and str(c[key]).upper() == str(d.get(key, "")).upper():
                        same[side].append(d.get("id"))
        if near:
            short.append("≈" + ", ".join(map(str, near)))
            tips.append("Both primer sites overlap those of " + ", ".join(map(str, near)) + ": the pairs amplify nearly the same product.")
        for side in "FR":
            if same[side] and not near:
                short.append(f"same {side} as " + ", ".join(map(str, same[side])))
                tips.append(f"The {'forward' if side == 'F' else 'reverse'} primer is the same sequence as in "
                            + ", ".join(map(str, same[side])) + ": order it once.")
        out.append((" · ".join(short) or "—", " ".join(tips)))
    return out


NOTES_GLOSS = (f"A primer whose 3′-terminal {primers_mod.MIN_PRIMER_SIZE} nt lie inside one copy of a terminal repeat may also bind "
               f"the other copy and solo copies elsewhere ({primers_mod.MIN_PRIMER_SIZE} nt is the shortest 3′ core in-silico PCR "
               "accepts as a binding site). ≈ marks pairs whose forward and reverse sites both overlap: they amplify nearly the "
               "same product. “same F/R as” marks a primer shared by two pairs; order it once.")
LOCATION_GLOSS = ("Product span on the analysed record, then where each primer sits: 5′ or 3′ flank, 5′ or 3′ terminal repeat, or "
                  "internal (from the accepted terminal-repeat pair in panel 02), with the protein domain when the primer lies "
                  "inside one. A → B: forward in A, reverse in B. A pair designed on a sub-sequence shows its span on that "
                  "sequence.")


# --------------------------------------------------------------------------------------------- order sheet
def pair_token(text):
    """A pair label as a filename / FASTA-header token ('D2·P3' -> 'D2-P3')."""
    return re.sub(r"[^A-Za-z0-9_.]+", "-", str(text)).strip("-") or "pair"


def oligo_name(prefix, pair, side):
    """Filename- and vendor-safe oligo name: '<record id>_<pair>_<F|R>'."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", f"{prefix}_{pair}_{side}" if prefix else f"{pair}_{side}")


def conditions_mix(c, temp=False):
    """'50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM primer' from a conditions dict, None when a field is missing."""
    try:
        return ((f"{c['temp_c']:g} °C, " if temp else "") + f"{c['mv_conc']:g} mM Na⁺, {c['dv_conc']:g} mM Mg²⁺, "
                f"{c['dntp_conc']:g} mM dNTP, {c['dna_conc']:g} nM primer")
    except (KeyError, TypeError, ValueError):
        return None


def tm_header(cands):
    """Tm column header naming the design conditions behind the numbers; the generic form when a pair carries none."""
    found = {conditions_mix(((c.get("design_provenance") or {}).get("parameters") or {}).get("tm_conditions") or {})
             for c in cands}
    if len(found) == 1 and None not in found:
        return f"Tm (°C; Primer3: {found.pop()})"
    return ("Tm (°C; Primer3 design conditions — differ between pairs, see each design manifest)" if len(found) > 1 else
            "Tm (°C; Primer3 design conditions)")


def order_sheet(entries):
    """(headers, rows, paste text) for [(pair label, candidate)]: one row per unique oligo sequence, forward then
    reverse of each pair. A primer shared by two pairs is listed once and names both under 'Used in'. Tm and GC are the
    Primer3 design values. Scale and Purification are left empty for the user to fill in the vendor form."""
    uniq = {}
    for label, c in entries:
        for side, sk, tk, gk in (("F", "left_seq", "left_tm", "left_gc"), ("R", "right_seq", "right_tm", "right_gc")):
            seq = str(c.get(sk) or "").upper()
            if not seq:
                continue
            used = f"{label} {side}" + (f" ({c['product_size']} bp)" if c.get("product_size") is not None else "")
            if seq in uniq:
                uniq[seq]["used"].append(used)
                continue
            uniq[seq] = {"name": oligo_name(c.get("_name_prefix"), pair_token(label), side), "seq": seq, "len": len(seq),
                         "tm": c.get(tk), "gc": c.get(gk), "used": [used]}
    headers = ["Name", "Sequence (5'-3')", "Scale", "Purification", "Length (nt)",
               tm_header([c for _l, c in entries]), "GC (%)", "Used in"]
    rows = [[u["name"], u["seq"], "", "", u["len"], "" if u["tm"] is None else u["tm"],
             "" if u["gc"] is None else u["gc"], "; ".join(u["used"])] for u in uniq.values()]
    return headers, rows, "\n".join(f"{r[0]}\t{r[1]}" for r in rows)


# ----------------------------------------------------------------------------------------------- widgets
class CoordTable(DataTable):
    """A DataTable whose FILE and CLIPBOARD output can differ from its screen.

    canon_rows (rows by original index) are what an export writes when the screen shows a converted coordinate
    convention, including a capped table whose sort runs over the complete row set. display_only_columns never reach
    an export, a row copy or Ctrl+C. Cell copy stays what the cell shows."""

    def __init__(self, headers, tooltips=None, parent=None):
        super().__init__(headers, tooltips, parent)
        self.canon_rows = None
        self.display_only_columns = frozenset()

    def _file_cols(self):
        return [j for j in range(self.columnCount()) if j not in self.display_only_columns]

    def export_headers(self):
        heads = super().export_headers()
        return [h for j, h in enumerate(heads) if j not in self.display_only_columns]

    def rows_data(self):
        if self.canon_rows is None:
            rows = super().rows_data()
        elif self._complete_rows is not None:
            rows = [list(self.canon_rows[i]) for i, _ in self._sorted_complete_pairs()]
        else:
            rows = [["" if v is None else str(v) for v in self.canon_rows[self._orig(i)]] for i in range(self.rowCount())]
        keep = self._file_cols()
        return rows if len(keep) == self.columnCount() else [[r[j] for j in keep if j < len(r)] for r in rows]

    def _copy_row(self, row):
        if row >= 0:
            QApplication.clipboard().setText("\t".join(self._value(row, j) for j in self._file_cols()))

    def keyPressEvent(self, e):
        if self.display_only_columns and e.matches(QKeySequence.Copy):
            rows = sorted({i.row() for i in self.selectedIndexes()}) or ([self.currentRow()] if self.currentRow() >= 0 else [])
            if rows:
                QApplication.clipboard().setText("\n".join("\t".join(self._value(r, j) for j in self._file_cols()) for r in rows))
            e.accept()
            return
        super().keyPressEvent(e)


class CoordGenomePanel(GenomePanel):
    """Genome viewer whose position label and hover text also give the 1-based reading when coord_base is 1. The figure
    itself, and every SVG, PNG or PDF export, keep their 0-based half-open text."""
    coord_base = 0

    def _render(self):
        super()._render()
        if not self.coord_base:
            return
        s, e = int(self.view["start"]), int(self.view["end"])
        self.pos.setText(f"{s + 1}–{e} bp · figure text 0-based, half-open · {(self.view['end'] - self.view['start']) / 1000:.2f} kb")
        for r in getattr(self.canvas, "regions", None) or []:
            if r.get("start") is not None and r.get("end") is not None:
                r["tip"] = (r.get("tip") or "") + (f"\n1-based, inclusive: {figures.fmt_interval(r['start'] + 1, r['end'])} "
                                                   "(figure text is 0-based, half-open)")
