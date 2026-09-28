"""Tests for transplant.py's own logic: which deletions, how drawn, how judged.

Run: python3 -m pytest scripts/transplant/test_transplant.py
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import transplant as t  # noqa: E402


def test_a_pure_deletion_is_its_deleted_bases_zero_based_half_open():
    # POS 1001 (1-based) is the anchor; the 60 bases after it are deleted.
    assert t.deletion_of("chr1", 1001, "A" + "C" * 60, "A") == ("chr1", 1001, 1061)


def test_only_pure_deletions_of_50_bp_or_more_are_taken():
    assert t.deletion_of("chr1", 1001, "A" + "C" * 49, "A") is None      # 49 bp
    assert t.deletion_of("chr1", 1001, "A" + "C" * 60, "G") is None      # ALT is not the anchor
    assert t.deletion_of("chr1", 1001, "A" + "C" * 60, "AT") is None     # a deletion with inserted bases


def test_size_bins_follow_the_plan():
    assert [t.size_bin(n) for n in (50, 299, 300, 999, 1000, 9999, 10000)] == [
        "50-299", "50-299", "300-999", "300-999", "1k-10k", "1k-10k", "10k+"]


def test_het_means_one_copy_phased_or_not():
    assert [t.is_het(g) for g in ("0/1", "1/0", "0|1", "1|0", "1/1", "0/0", "./1")] == [
        True, True, True, True, False, False, False]


def test_a_deletion_with_another_sv_within_1_kb_is_dropped():
    events = [("chr1", 10_000, 10_500), ("chr1", 50_000, 50_500), ("chr2", 10_000, 10_500)]
    others = [("chr1", 10_000, 10_500),   # itself: ignored
              ("chr1", 11_400, 11_450),   # 900 bp past the first one's end
              ("chr1", 51_600, 51_700),   # 1,100 bp past the second one's end
              ("chr2", 10_000, 10_500),   # the third event itself: ignored
              ("chr3", 10_000, 10_500)]   # another chromosome
    assert t.isolated(events, others) == [("chr1", 50_000, 50_500), ("chr2", 10_000, 10_500)]


def test_the_draw_is_seeded_and_takes_the_first_that_pass_in_shuffled_order():
    cands = [("chr1", i * 1000, i * 1000 + 400) for i in range(50)]
    a = t.draw(cands, 10, seed=1, keep=lambda e: True)
    assert a == t.draw(list(reversed(cands)), 10, seed=1, keep=lambda e: True)  # input order does not matter
    assert a != t.draw(cands, 10, seed=2, keep=lambda e: True)
    assert len(a) == 10 and len(set(a)) == 10
    odd = t.draw(cands, 10, seed=1, keep=lambda e: e[1] // 1000 % 2 == 1)
    assert all(e[1] // 1000 % 2 == 1 for e in odd) and len(odd) == 10
    assert len(t.draw(cands, 100, seed=1, keep=lambda e: True)) == 50  # fewer than asked: all of them


def test_events_in_one_spike_run_are_at_least_100_kb_apart():
    events = [("chr1", 0, 500), ("chr1", 50_000, 50_500), ("chr1", 200_000, 200_500), ("chr2", 0, 500)]
    runs = t.batches(events, gap=100_000)
    assert runs == [[("chr1", 0, 500), ("chr1", 200_000, 200_500), ("chr2", 0, 500)],
                    [("chr1", 50_000, 50_500)]]


def test_refused_events_are_read_back_from_spike_s_message():
    stderr = (
        "Error: The reads over 2 events would carry less than half of what truth.vcf would claim, "
        "so spike stops rather than write it (RF8):\n"
        "  DEL  chr20:27100001-27110000 (10000bp): 5731 of 5749 reads over it (99.7%) are ones spike cannot edit\n"
        "  DEL  chr1:501-1000 (500bp): 30 of 50 reads over it (60.0%) are ones spike cannot edit\n")
    assert t.parse_refused(stderr) == [("chr20", 27_100_000, 27_110_000), ("chr1", 500, 1_000)]


def test_percentiles_interpolate_linearly():
    assert t.percentile([1, 2, 3, 4], 50) == 2.5
    assert t.percentile([1, 2, 3, 4], 25) == pytest.approx(1.75)
    assert t.percentile([7], 90) == 7


def test_the_pilot_events_are_left_out():
    pool = [("chr1", 0, 400), ("chr1", 5_000, 5_400), ("chr2", 0, 400)]
    assert t.without(pool, {("chr1", 5_000, 5_400)}) == [("chr1", 0, 400), ("chr2", 0, 400)]


PASS, FAIL = {"pass": True}, {"pass": False}


def test_a_metric_is_judged_only_when_the_controls_that_must_fail_it_fail():
    assert t.verdict(PASS, {"B1": FAIL, "B2": FAIL}, "J") == "pass"
    assert t.verdict(FAIL, {"B1": FAIL, "B2": FAIL}, "J") == "fail"
    assert t.verdict(PASS, {"B1": FAIL, "B2": PASS}, "J") == "inconclusive"   # B2 blind to J
    assert t.verdict(PASS, {"B1": PASS, "B2": FAIL}, "J") == "inconclusive"
    assert t.verdict(PASS, {"B1": FAIL, "B2": PASS}, "E1") == "pass"          # B2 is not asked about E1
    assert t.verdict(PASS, {"B1": PASS, "B2": FAIL}, "E1") == "inconclusive"
    assert t.verdict(PASS, {"B2": FAIL}, "E1") == "inconclusive"              # B1 missing


def test_a_bin_is_supported_only_when_every_metric_passes():
    assert t.overall(["pass", "pass"]) == "supported"
    assert t.overall(["pass"]) == "supported"
    assert t.overall(["pass", "fail"]) == "refuted"
    assert t.overall(["fail", "inconclusive"]) == "refuted"
    assert t.overall(["pass", "inconclusive"]) == "inconclusive"
    assert t.overall([]) == "inconclusive"


def test_the_pass_rule_needs_the_median_inside_the_middle_half_and_no_wider_spread():
    c = list(range(-10, 11))                 # 25th..75th: -5..5; 10th..90th width 16
    assert t.judge([0, 1, -1, 2, -2], c)["pass"]
    assert not t.judge([6, 7, 8, 9, 10], c)["pass"]           # median 8: biased
    wide = [-30, -20, 0, 20, 30]                              # median 0, width 52 > 1.5 x 16
    j = t.judge(wide, c)
    assert j["bias_ok"] and not j["spread_ok"] and not j["pass"]
    edge = t.judge([5, 5, 5], c)                              # median exactly on the 75th: inside
    assert edge["bias_ok"]
