"""Tests for round3.py: split-half counting noise at each site, and the round 3 rule.

Run: python3 -m pytest scripts/transplant/test_round3.py
"""
import hashlib
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import round3 as r3  # noqa: E402
import small_evidence as se  # noqa: E402
import test_small_evidence as tse  # noqa: E402


def test_a_pair_s_half_is_the_plan_s_hash_of_salt_and_name():
    names = [f"read{i}" for i in range(2_000)]
    halves = [r3.half_of(n, 0) for n in names]
    assert halves == [hashlib.blake2b(f"0:{n}".encode()).digest()[0] & 1 for n in names]
    assert 900 < sum(halves) < 1_100
    assert halves != [r3.half_of(n, 1) for n in names]


DUP_POS = 11_000  # a 100 bp duplication of [11000, 11100) in tse.REF, anchored at POS 11000
DUP = ("DUP", tse.CHROM, DUP_POS, tse.REF[DUP_POS - 1], tse.REF[DUP_POS - 1] + tse.REF[DUP_POS:DUP_POS + 100])


def _dup_bams(tmp_path):
    s, e = DUP_POS, DUP_POS + 100
    main = tse.background((s, e))
    # 30 pairs clipped at the copy's end (a mate each side), 5 of them replaced by spike
    for i in range(30):
        main += [tse.segment(f"c{i}", e - 90, "90M40S"), tse.segment(f"c{i}", e + 300, "100M", flag=16)]
    sim = []
    for i in range(20):  # spike's pairs, clipped at the copy's start
        sim += [tse.segment(f"s{i}", s - 300, "100M"), tse.segment(f"s{i}", s, "40S90M", flag=16)]
    skip = frozenset(f"c{i}" for i in range(5))
    return (tse.write_bam(tmp_path / "main.bam", main), tse.write_bam(tmp_path / "sim.bam", sim), skip)


def test_measure_on_all_reads_is_round_2_s_measure(tmp_path):
    main, sim, skip = _dup_bams(tmp_path)
    got = r3.measure_on(r3.paths_for(main, skip, [sim]), tse.fetch, DUP)
    assert got == se.measure(main, tse.fetch, DUP, skip, [sim])
    assert got["n_any"] == 25 + 20
    snv = ("SNV", tse.CHROM, 9_001, "C", "A")
    assert r3.measure_on(r3.paths_for(main), tse.fetch, snv) == se.measure(main, tse.fetch, snv)


def test_the_two_halves_share_out_every_counted_pair_and_both_honour_the_skip_set(tmp_path):
    main, sim, skip = _dup_bams(tmp_path)
    full = r3.measure_on(r3.paths_for(main, skip, [sim]), tse.fetch, DUP)
    for salt in range(3):
        h0, h1 = (r3.measure_on(r3.paths_for(main, skip, [sim], half=(salt, h)), tse.fetch, DUP) for h in (0, 1))
        assert h0["n_any"] + h1["n_any"] == full["n_any"]
        assert 0 < h0["n_any"] < full["n_any"]
        assert h0["flank_depth"] + h1["flank_depth"] == pytest.approx(full["flank_depth"])


def test_split_noise_is_the_mean_squared_half_gap_over_4_skipping_undefined_halves():
    assert r3.split_noise([(0.5, 0.3), (0.4, 0.4), (None, 0.2), (0.1, None)]) == pytest.approx((0.2 ** 2 / 4 + 0) / 2)
    assert r3.split_noise([(None, 0.1)]) is None
    assert r3.split_noise([]) is None


def test_mismatch_ratio_is_squared_differences_over_summed_noise():
    assert r3.mismatch_ratio([(0.2, 0.01), (-0.1, 0.03)]) == pytest.approx((0.04 + 0.01) / 0.04)
    assert r3.mismatch_ratio([(0.0, 0.0)]) is None


def test_items_pair_each_fake_with_its_real_twin_and_add_their_noise():
    e1, e2, e3 = (("SNV", "chr1", p, "A", "G") for p in (1, 2, 3))
    fake = {e1: {"A": (0.4, 0.01)}, e2: {"A": (0.5, None)}, e3: {"A": (0.3, 0.02)}}
    real = {e1: {"A": (0.5, 0.02)}, e2: {"A": (0.5, 0.01)}}
    assert r3.items_of(fake, real, "A") == [pytest.approx((-0.1, 0.03))]
    real[e3] = {"A": (None, 0.01)}
    assert r3.items_of(fake, real, "A") == [pytest.approx((-0.1, 0.03))]
    real[e1] = {"A": (0.5, None)}
    assert r3.items_of(fake, real, "A") == []


def test_rho_interval_is_seeded_and_resamples_each_side_apart():
    same = [((-1) ** i * 0.3, 0.01) for i in range(40)]
    assert r3.rho_interval(same, [(x / 3, v) for x, v in same], draws=50, seed=3) == pytest.approx((9, 9))
    f = [(0.1, 0.01)] * 10                              # R_f is 1 in every draw
    c = [(0.1, 0.01)] * 45 + [(1.0, 0.01)] * 5          # R_c is about 10.9 on average
    lo, hi = r3.rho_interval(f, c, draws=500, seed=3)
    assert hi < 0.5                                     # c drawn from all 50, not from f's 10 slots
    varied = [(0.05 * i, 0.01) for i in range(1, 51)]
    got = r3.rho_interval(f, varied, draws=500, seed=3)
    assert got == r3.rho_interval(f, varied, draws=500, seed=3)
    assert got != r3.rho_interval(f, varied, draws=500, seed=4)
    lo, hi = r3.rho_interval([(0.1, 0.0), (0.2, 0.01)], c, draws=200, seed=3)   # draws with no noise are skipped
    assert lo is not None and hi < float("inf")
    assert r3.rho_interval([(0.1, 0.0)], c, draws=20, seed=3) == (None, None)


def test_the_interval_is_the_5th_and_95th_percentiles():
    assert r3.interval(list(range(1, 101))) == pytest.approx((5.95, 95.05))
    assert r3.interval([]) == (None, None)


def _result(key, which, e, value, v=0.0025):
    return (key, which, e, {}, {"A": (value, v, 10)})


def test_verdicts_pair_fakes_with_real_twins_and_lean_on_the_forward_b1(tmp_path):
    events = [("SNV", "chr1", 1_000 * i, "C", "T") for i in range(1, 31)]
    sign = lambda i: (-1) ** i  # noqa: E731
    res = []
    for i, e in enumerate(events):          # c = +-0.1 over noise 0.005: R_c = 2
        res += [_result("shared", "HG001", e, 0.5 + 0.05 * sign(i)), _result("shared", "HG002", e, 0.5 - 0.05 * sign(i))]
        res += [_result("forward", "real", e, 0.5), _result("reverse", "real", e, 0.5)]
        res += [_result("forward", "fake_normal", e, 0.5 + 0.05 * sign(i)),    # R_f 0.5: rho 0.25
                _result("forward", "fake_B1", e, 0.25),                        # rho 6.25
                _result("reverse", "fake_normal", e, 0.5 + 0.3 * sign(i))]     # rho 9
    r3.judge_all(str(tmp_path), res, draws=200)
    rows = [line.rstrip("\n").split("\t") for line in open(tmp_path / "verdicts3.tsv")][1:]
    got = {(s, g, m): v for s, g, m, v in rows if g == "SNV"}
    assert got == {("forward", "SNV", "A"): "pass", ("forward", "SNV", "overall"): "supported",
                   ("reverse", "SNV", "A"): "fail", ("reverse", "SNV", "overall"): "refuted"}
    judge = {(r[0], r[1]): r for r in (line.rstrip("\n").split("\t") for line in open(tmp_path / "judge3.tsv"))}
    assert judge[("forward", "normal")][8:11] == ["0.5000", "2.0000", "0.2500"]   # R_f, R_c, rho
    assert judge[("forward", "B1")][-1] == "fail"


def test_a_metric_passes_below_2_25_fails_above_and_is_unsure_across():
    assert r3.judge3(0.5, 2.25) == "pass"
    assert r3.judge3(0.5, 2.26) == "unsure"
    assert r3.judge3(2.25, 3.0) == "unsure"
    assert r3.judge3(2.26, 3.0) == "fail"
    assert r3.judge3(None, None) == "unsure"


def test_a_metric_counts_only_when_each_must_fail_control_fails():
    assert r3.verdict3("pass", {"B1": "fail"}, "A") == "pass"
    assert r3.verdict3("pass", {"B1": "unsure"}, "A") == "inconclusive"
    assert r3.verdict3("pass", {"B1": "pass"}, "A") == "inconclusive"
    assert r3.verdict3("fail", {"B1": "fail"}, "E") == "fail"
    assert r3.verdict3("unsure", {"B1": "fail"}, "A") == "unsure"
    assert r3.verdict3("pass", {"B1": "fail"}, "J") == "inconclusive"          # B2 missing
    assert r3.verdict3("pass", {"B1": "fail", "B2": "fail"}, "J") == "pass"
    assert r3.verdict3("pass", {"B1": "fail", "B2": "unsure"}, "J") == "inconclusive"


def test_calibration_compares_split_noise_with_binomial_noise_on_snvs():
    rows = [{"A": 0.5, "N": 40, "v": 0.00625}, {"A": 0.2, "N": 20, "v": 0.016},
            {"A": None, "N": 0, "v": None}, {"A": 0.3, "N": 30, "v": None}]
    assert r3.calibration(rows) == pytest.approx((0.00625 + 0.016) / (0.25 / 40 + 0.16 / 20))
    assert r3.calibrated(0.8) and r3.calibrated(1.25)
    assert not r3.calibrated(0.79) and not r3.calibrated(1.26)


def test_a_value_matches_round_2b_only_as_written():
    assert r3.same_as_written(0.16892, "0.1689")
    assert not r3.same_as_written(0.16896, "0.1689")
    assert r3.same_as_written(None, "")
    assert not r3.same_as_written(None, "0.0000")
    assert r3.same_as_written(12, "12")


def test_events_get_the_run_round2_gave_them_and_must_lie_in_its_bed():
    events = [("SNV", "chr1", p, "C", "T") for p in (1_000, 50_000, 200_000)] + [("SNV", "chr2", 1_000, "C", "T")]
    span = lambda e: (e[1], e[2] - 1, e[2])  # noqa: E731
    assert r3.runs_of(events, span) == {events[0]: 0, events[1]: 1, events[2]: 0, events[3]: 0}
    bed = [("chr1", 900, 1_100), ("chr2", 0, 10)]
    assert r3.in_run_bed(("SNV", "chr1", 1_000, "C", "T"), bed)
    assert r3.in_run_bed(("SNV", "chr1", 901, "C", "T"), bed)
    assert not r3.in_run_bed(("SNV", "chr1", 900, "C", "T"), bed)            # POS 900 is base 899
    assert not r3.in_run_bed(("SNV", "chr3", 5, "C", "T"), bed)
