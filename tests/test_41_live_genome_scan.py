"""Live check of the whole-genome scan's coordinate contract against real isPcr output.

Needs the WSL backend and network access (yeast R64 download, NCBI coordinate fetch), so it is excluded
from the hermetic suite: run with  pytest -m "wsl and network" tests/test_41_live_genome_scan.py"""
import pytest

pytestmark = [pytest.mark.wsl, pytest.mark.network]

ORG, ACC, START, END = "Saccharomyces cerevisiae", "NC_001136.10", 100001, 102000


@pytest.mark.parametrize("strand", ["+", "-"])
def test_designed_product_is_the_single_on_target_hit(strand):
    import engine
    from teagle_core import sequtil
    assert engine.run_genome_prepare({"organism": ORG}).get("ok")
    fetched = engine.run_fetch_coords({"regions": f"chrIV:{START}-{END}", "organism": ORG, "strand": strand})
    seq = sequtil.parse_fasta(fetched["fasta"])[0][1]
    cand = engine.run_primers({"sequence": seq, "params": {"prod_min": 150, "prod_max": 500}})["candidates"][0]
    a, b = cand["intended_product_span"]
    res = engine.run_genome_pcr({
        "fwd": cand["left_seq"], "rev": cand["right_seq"], "organism": ORG,
        "design_locus": [{"accession": ACC, "start": START, "end": END, "strand": 1 if strand == "+" else 2,
                          "intended_product_span": [a, b]}]})
    lo, hi = (START + a, START + b - 1) if strand == "+" else (END - b + 1, END - a)
    on = [x for x in res["amplicons"] if x.get("on_target")]
    assert res["n_on_target"] == 1 and len(on) == 1
    assert (on[0]["start"], on[0]["end"]) == (lo, hi)
    assert on[0]["length"] == hi - lo + 1 == b - a
