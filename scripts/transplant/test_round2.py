"""Tests for round2.py's own logic: which variants, how spiked, how judged.

Run: python3 -m pytest scripts/transplant/test_round2.py
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import round2 as r2  # noqa: E402


def test_a_record_is_an_snv_or_a_1_to_49_bp_pure_indel():
    assert r2.kind_of("C", "T") == ("SNV", 1)
    assert r2.kind_of("GCA", "G") == ("DEL", 2)
    assert r2.kind_of("G", "GCA") == ("INS", 2)
    assert r2.kind_of("G" + "A" * 49, "G") == ("DEL", 49)
    assert r2.kind_of("G" + "A" * 50, "G") is None         # 50 bp: an SV
    assert r2.kind_of("CA", "TG") is None                   # MNV
    assert r2.kind_of("GCA", "T") is None                   # complex: the anchor changes
    assert r2.kind_of("G", "TCA") is None


def test_groups_follow_the_plan():
    assert [r2.group_of(k, n) for k, n in (("SNV", 1), ("DEL", 1), ("DEL", 4), ("DEL", 5), ("DEL", 19),
                                            ("DEL", 20), ("INS", 49), ("DUP", 50), ("DUP", 299))] == [
        "SNV", "DEL1-4", "DEL1-4", "DEL5-19", "DEL5-19", "DEL20-49", "INS20-49", "DUP50-299", "DUP50-299"]
    assert r2.group_of("DUP", 300) is None
    assert r2.metrics_of("SNV") == ["A"]
    assert r2.metrics_of("INS5-19") == ["A", "E"]
    assert r2.metrics_of("DUP50-299") == ["J"]


HET, HOM = "0|1", "1|1"


def test_transplants_are_het_in_one_sample_and_absent_from_the_other_and_shared_are_het_in_both():
    region = {"chr1": [(0, 100_000)]}
    h1 = [("chr1", 1_000, "C", "T", HET),         # HG001 only: forward
          ("chr1", 5_000, "C", "T", "1|0"),       # in both, het in both: shared
          ("chr1", 9_000, "C", "T", HET),         # HG002 hom: neither
          ("chr1", 13_000, "C", "T", HOM),        # HG001 hom: neither
          ("chr1", 17_000, "CA", "TG", HET),      # MNV: neither
          ("chr1", 200_000, "C", "T", HET)]       # outside the region
    h2 = [("chr1", 5_000, "C", "T", "0/1"), ("chr1", 9_000, "C", "T", HOM),
          ("chr1", 21_000, "GCA", "G", HET),      # HG002 only: reverse
          ("chr1", 25_000, "C", "T", ".|."),      # no genotype: neither
          ("chr1", 29_000, "C", "T", HET)]        # HG002 only,
    h1.append(("chr1", 29_000, "C", "A", HET))    #   but HG001 has another allele there: neither
    sets = r2.classify(h1, h2, region, margin=150)
    assert sets["forward"] == [("SNV", "chr1", 1_000, "C", "T")]
    assert sets["reverse"] == [("DEL", "chr1", 21_000, "GCA", "G")]
    assert sets["shared"] == [("SNV", "chr1", 5_000, "C", "T")]


def test_a_small_variant_with_another_record_within_150_bp_is_dropped():
    region = {"chr1": [(0, 100_000)]}
    h1 = [("chr1", 1_000, "C", "T", HET), ("chr1", 1_140, "A", "G", HOM),     # 140 bp on: dropped
          ("chr1", 5_000, "C", "T", HET),
          ("chr1", 9_000, "GCA", "G", HET)]
    h2 = [("chr1", 5_160, "A", "G", HOM),                                       # 160 bp on: kept
          ("chr1", 8_850, "A", "G", HOM)]    # ends 150 bp before 9,000
    sets = r2.classify(h1, h2, region, margin=150)
    assert sets["forward"] == [("SNV", "chr1", 5_000, "C", "T")]
    others = [("chr1", 5_100, 5_300)]                                          # an SV of either sample
    assert r2.classify(h1, h2, region, margin=150, svs=others)["forward"] == []


def test_a_variant_must_lie_inside_one_region_interval():
    region = {"chr1": [(0, 1_000), (1_000, 2_000)]}
    h1 = [("chr1", 500, "C", "T", HET), ("chr1", 999, "GCAT", "G", HET)]     # [999, 1002] crosses 1,000
    assert r2.classify(h1, [], region, margin=0)["forward"] == [("SNV", "chr1", 500, "C", "T")]


def test_a_duplication_s_isolation_ignores_its_own_records():
    dups = [(("DUP", "chr1", 10_000, "A", "A" + "C" * 100), [("chr1", 10_000, 10_001)]),
            (("DUP", "chr1", 50_000, "A", "A" + "C" * 100), [("chr1", 50_000, 50_001), ("chr1", 50_004, 50_005)]),
            (("DUP", "chr1", 90_000, "A", "A" + "C" * 100), [("chr1", 90_000, 90_001)])]
    svs = [("chr1", 10_000, 10_001), ("chr1", 50_000, 50_001), ("chr1", 50_004, 50_005),   # their own
           ("chr1", 90_900, 90_950)]                                                        # 1 kb from the third
    kept = r2.isolated_dups(dups, svs, margin=1_000, segment=lambda e: (e[2], e[2] + 100))
    assert kept == [dups[0][0], dups[1][0]]


def test_spike_specs():
    snv, dele, dup = ("SNV", "chr1", 1_001, "C", "T"), ("DEL", "chr1", 5_000, "GCA", "G"), ("DUP", "chr1", 10_000, "A", "AC")
    assert r2.event_spec(snv, "normal") == "snp:chr1:1001:C:T;af=0.5"
    assert r2.event_spec(dele, "B1") == "snp:chr1:5000:GCA:G;af=0.25"
    assert r2.event_spec(dup, "normal", segment=(10_000, 10_100)) == "dup:chr1:10000-10100;af=0.5"
    assert r2.event_spec(dup, "B2", segment=(10_000, 10_100)) == "dup:chr1:10200-10300;af=0.5"
    with pytest.raises(ValueError):
        r2.event_spec(snv, "B2")


def test_refused_events_are_found_from_spike_s_labels():
    stderr = (
        "Error: The reads over 3 events would carry less than half ... (RF8):\n"
        "  SNV  chr1:1001 C>T: 30 of 50 reads over it (60.0%) are ones spike cannot edit\n"
        "  SNV  chr2:5000 GCA>G: 30 of 50 reads over it (60.0%) are ones spike cannot edit\n"
        "  DUP  chr3:10201-10300 (100bp): 60 of 70 reads over it (85.7%) are ones spike cannot edit\n")
    listed = r2.parse_refused(stderr)
    events = [("SNV", "chr1", 1_001, "C", "T"), ("DEL", "chr2", 5_000, "GCA", "G"),
              ("DUP", "chr3", 10_000, "A", "AC"), ("SNV", "chr4", 7, "A", "G")]
    seg = {events[2]: (10_000, 10_100)}
    assert r2.refused_of(events, listed, "B2", seg.get) == events[:3]


def test_a_refusal_that_names_no_event_of_the_run_stops_the_run():
    # Unshifted, the DUP spike names is not this run's: rather stop than loop.
    stderr = "  DUP  chr3:10201-10300 (100bp): 60 of 70 reads over it (85.7%) are ones spike cannot edit\n"
    dup = ("DUP", "chr3", 10_000, "A", "AC")
    with pytest.raises(RuntimeError):
        r2.refused_of([dup], r2.parse_refused(stderr), "normal", {dup: (10_000, 10_100)}.get)
    listed = r2.parse_refused("  SNV  chr9:1 A>G: 3 of 5 reads over it (60.0%) are ones spike cannot edit\n")
    with pytest.raises(RuntimeError):
        r2.refused_of([("SNV", "chr1", 1_001, "C", "T")], listed, "normal", lambda e: None)


PASS, FAIL = {"pass": True}, {"pass": False}


def test_each_metric_needs_its_own_broken_controls_to_fail():
    assert r2.verdict(PASS, {"B1": FAIL}, "A") == "pass"
    assert r2.verdict(PASS, {"B1": PASS}, "A") == "inconclusive"
    assert r2.verdict(FAIL, {"B1": FAIL}, "E") == "fail"
    assert r2.verdict(PASS, {"B1": FAIL}, "J") == "inconclusive"             # B2 missing
    assert r2.verdict(PASS, {"B1": FAIL, "B2": FAIL}, "J") == "pass"
    assert r2.verdict(PASS, {"B1": FAIL, "B2": PASS}, "J") == "inconclusive"


def test_the_pilot_stops_when_b1_passes_any_metric():
    ok = [["forward", "B1", "SNV", "A", False], ["forward", "B1", "DUP50-299", "J", False],
          ["forward", "B2", "DUP50-299", "J", True],       # the plan's pilot stops on B1 only
          ["forward", "normal", "SNV", "A", True]]
    assert r2.pilot_stops(ok) == []
    bad = ok + [["forward", "B1", "INS20-49", "E", True]]
    assert r2.pilot_stops(bad) == ["B1 passes E in INS20-49"]
