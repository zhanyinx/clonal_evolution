"""Synthetic checks for bin/create_input4pyclone.py. Run: python tests/test_create_input4pyclone.py"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import create_input4pyclone as m  # noqa: E402


def test_count_bases():
    # 3 ref (. . ,), read-start mark ^] then A, read-end $, a, +2AC skipped, a, -1T skipped, t (other), * (deletion)
    assert m.count_bases("..,^]A$a+2ACa-1Tt*", "C", "A") == (3, 3)
    # the mapping-quality char after ^ is never counted as a base
    assert m.count_bases("^AA^aa", "C", "A") == (0, 2)
    # multi-digit indel lengths
    assert m.count_bases("+12ACGTACGTACGTa", "C", "A") == (0, 1)


def test_assign_cn():
    seg = pd.DataFrame(dict(chr=["chr1", "chr1", "chr1", "chr2"], startpos=[100, 300, 5_000_000, 10],
                            endpos=[200, 400, 6_000_000, 20], major_cn=[1, 1, 2, 1], minor_cn=[1, 1, 0, 1]))
    f = lambda c, p: m.assign_cn(seg, c, p, 1_000_000)  # noqa: E731
    assert f("chr1", 100) == (1, 1, "segment")  # 1-based inclusive start
    assert f("chr1", 200) == (1, 1, "segment")  # 1-based inclusive end
    assert f("chr1", 250) == (1, 1, "flanks_agree")
    assert f("chr1", 1000) == (1, 1, "nearest_segment")  # flanks disagree, left is 600 bp away
    assert f("chr1", 3_000_000)[2] == "no_segment_within_max_gap"
    assert f("chr1", 7_500_000)[2] == "no_segment_within_max_gap"
    assert f("chr3", 5)[2] == "no_segment_on_chromosome"


if __name__ == "__main__":
    test_count_bases()
    test_assign_cn()
    print("ok")
