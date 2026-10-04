"""Tests for the duplicates check (docs/superpowers/plans/2026-10-04-duplicates.md), counted by hand.

Run: python3 -B -m pytest scripts/duplicates/test_duplicates.py
"""
import os
import random
import sys

import pysam
import pytest

sys.path.insert(0, os.path.dirname(__file__))
import control  # noqa: E402
import draw_events  # noqa: E402
import measure  # noqa: E402

CHROM, LENGTH = "chrT", 200_000


def make_ref():
    rng = random.Random(3)
    seq = [rng.choice("ACGT") for _ in range(LENGTH)]
    seq[999] = "C"                   # 1-based 1000: the counted site
    seq[150_000:150_010] = "N" * 10  # an N run at 1-based 150001-150010
    return "".join(seq)


REF = make_ref()
HEADER = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                          "SQ": [{"SN": CHROM, "LN": LENGTH}]})


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / "ref.fa"
    path.write_text(f">{CHROM}\n" + "\n".join(REF[i:i + 60] for i in range(0, LENGTH, 60)) + "\n")
    pysam.faidx(str(path))
    with pysam.FastaFile(str(path)) as f:
        yield f


def segment(name, pos0, cigar="100M", flag=0x1 | 0x2 | 0x40, mapq=60, sub=None, qual=40, tags=()):
    """A read at 0-based `pos0`, its bases taken from REF except `sub` {0-based pos: base}."""
    s = pysam.AlignedSegment(HEADER)
    s.query_name, s.flag, s.reference_id, s.reference_start, s.mapping_quality = name, flag, 0, pos0, mapq
    s.next_reference_id, s.next_reference_start = 0, pos0
    s.cigarstring = cigar
    bases, rpos = [], pos0
    for op, n in s.cigartuples:
        if op == 0:
            bases += [(sub or {}).get(p, REF[p]) for p in range(rpos, rpos + n)]
            rpos += n
        elif op == 2:
            rpos += n
    s.query_sequence = "".join(bases)
    quals = [qual] * len(bases)
    s.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in quals))
    for k, v in tags:
        s.set_tag(k, v)
    return s


def write_bam(path, segments):
    path = str(path)
    with pysam.AlignmentFile(path, "wb", header=HEADER) as out:
        for s in sorted(segments, key=lambda s: s.reference_start):
            out.write(s)
    pysam.index(path)
    return pysam.AlignmentFile(path)


# --- the allele counter ----------------------------------------------------------------

def test_the_counter_reads_the_base_at_the_one_based_position(tmp_path, fasta):
    # Site 1-based 1000 = 0-based 999, reference C. Alt reads carry T there.
    alt = {999: "T"}
    bam = write_bam(tmp_path / "a.bam", [
        segment("ref1", 950), segment("ref2", 999),           # 999 is the first base of ref2
        segment("alt1", 900, sub=alt), segment("alt2", 990, sub=alt),
        segment("left", 899),                                 # covers 899..998: ends one base short
        segment("right", 1000),                               # starts one base past
    ])
    ref_names, total = measure.allele_counts(bam, fasta, CHROM, 1000)
    assert sorted(ref_names) == ["ref1", "ref2"]
    assert total == 4


def test_the_counter_skips_what_deepvariant_skips(tmp_path, fasta):
    bam = write_bam(tmp_path / "a.bam", [
        segment("kept", 950),
        segment("mapq5", 950, mapq=5),                        # exactly the floor: counted
        segment("mapq4", 950, mapq=4),
        segment("dup", 950, flag=0x1 | 0x2 | 0x40 | 0x400),
        segment("qcfail", 950, flag=0x1 | 0x2 | 0x40 | 0x200),
        segment("secondary", 950, flag=0x1 | 0x2 | 0x40 | 0x100),
        segment("supplementary", 950, flag=0x1 | 0x2 | 0x40 | 0x800),
        segment("bq9", 950, qual=9),
        segment("bq10", 950, qual=10),                        # exactly the floor: counted
        segment("deleted", 950, cigar="40M20D60M"),           # 990..1009 deleted: no base at 999
    ])
    ref_names, total = measure.allele_counts(bam, fasta, CHROM, 1000)
    assert sorted(ref_names) == ["bq10", "kept", "mapq5"]
    assert total == 3


def test_pooled_share_is_pooled_reads_not_a_mean_of_sites():
    # 1 of 2 at one site, 0 of 8 at another: pooled 1/10, not the mean 0.25.
    assert measure.pooled([(["a"], 2), ([], 8)]) == (1, 10)


# --- the verdict -----------------------------------------------------------------------

@pytest.mark.parametrize("r_spk, r_real, leak, want", [
    (0.050, 0.010, 0.80, "matters"),
    (0.050, 0.010, 0.50, "matters"),            # L exactly 0.5 is enough
    (0.050, 0.010, 0.49, "something else"),
    (0.020, 0.000, 0.90, "does not matter"),    # exactly 0.02 over is not more than 0.02
    (0.030, 0.000, 0.90, "matters"),
    (0.010, 0.010, 0.00, "does not matter"),
])
def test_the_verdict_follows_the_locked_rule(r_spk, r_real, leak, want):
    assert measure.verdict(r_spk, r_real, leak) == want


def test_old_allele_reads_split_four_ways():
    split = measure.split_ref_reads(
        ["SPIKE_ev0001_hap_000001:46:FC:2:1101:1:1", "srcdup", "untaken", "kept"],
        source_dups={"srcdup"}, replaced={"kept"})
    assert split == {"source duplicate": 1, "not taken": 1, "replaced": 1, "SPIKE_": 1}


# --- the coverage test -----------------------------------------------------------------

def covering(n, start=0, **kw):
    return [segment(f"r{start + i}", 950, **kw) for i in range(n)]


@pytest.mark.parametrize("n, ok", [(19, False), (20, True), (45, True), (46, False)])
def test_coverage_is_20_to_45_records(tmp_path, n, ok):
    bam = write_bam(tmp_path / "a.bam", covering(n))
    assert measure.coverage_ok(bam, CHROM, 1000) is ok


def test_coverage_counts_no_duplicate_and_needs_95_percent_good(tmp_path):
    dups = covering(10, start=100, flag=0x1 | 0x2 | 0x40 | 0x400)
    bam = write_bam(tmp_path / "a.bam", covering(19) + dups)          # 19 counted: too few
    assert measure.coverage_ok(bam, CHROM, 1000) is False
    # 40 records, 2 of them MAPQ 19: 38/40 = 0.95, enough. A third (37/40) is not.
    two = covering(38) + covering(2, start=50, mapq=19)
    assert measure.coverage_ok(write_bam(tmp_path / "b.bam", two), CHROM, 1000) is True
    three = covering(37) + covering(3, start=50, mapq=19)
    assert measure.coverage_ok(write_bam(tmp_path / "c.bam", three), CHROM, 1000) is False
    improper = covering(37) + covering(3, start=50, flag=0x1 | 0x40)
    assert measure.coverage_ok(write_bam(tmp_path / "d.bam", improper), CHROM, 1000) is False


# --- the event draw --------------------------------------------------------------------

def even_bam(tmp_path):
    """25 good records over every base of chrT (100 bp reads every 4 bp), plus 10 MAPQ 0 ones
    over 1-based 60001-60100 (25 of 35 there is under 0.95)."""
    reads = [segment(f"e{p}", p) for p in range(0, LENGTH - 100, 4)]
    bad = [segment(f"low{k}", 60_000, mapq=0) for k in range(10)]
    return write_bam(tmp_path / "even.bam", reads + bad)


def test_the_draw_keeps_every_locked_rule(tmp_path, fasta):
    bam = even_bam(tmp_path)
    calls = [100_000]
    events = draw_events.draw(bam, fasta, calls, CHROM, 5_000, 195_000, seed=7,
                              n_het=3, n_hom=3, n_del=1, max_candidates=20_000)
    kinds = [e["kind"] for e in events]
    assert kinds == ["het"] * 3 + ["hom"] * 3 + ["del"]
    for e in events:
        lo, hi = e["start"], e["end"]
        assert 5_000 <= lo and hi <= 195_000
        assert all(abs(c - lo) > 1_000 and abs(c - hi) > 1_000 for c in calls)
        assert not (lo - 1_000 <= 150_010 and hi + 1_000 >= 150_001), "no N within 1 kb"
        assert not (lo <= 60_100 and hi >= 60_001), "the low-MAPQ stretch fails coverage"
        if e["kind"] != "del":
            assert lo == hi
            assert (e["ref"], e["alt"]) in {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
            assert e["ref"] == REF[lo - 1]
        else:
            assert hi - lo + 1 == 300
    for i, a in enumerate(events):
        for b in events[i + 1:]:
            assert max(a["start"], b["start"]) - min(a["end"], b["end"]) >= 20_000
    again = draw_events.draw(bam, fasta, calls, CHROM, 5_000, 195_000, seed=7,
                             n_het=3, n_hom=3, n_del=1, max_candidates=20_000)
    assert again == events


def test_the_draw_refuses_a_site_within_1_kb_of_a_call_or_an_n(tmp_path, fasta):
    bam = even_bam(tmp_path)
    one = dict(seed=7, n_het=1, n_hom=0, n_del=0, max_candidates=500)
    # Every position in 99,500-100,500 is within 1 kb of the call at 100,000.
    with pytest.raises(SystemExit):
        draw_events.draw(bam, fasta, [100_000], CHROM, 99_500, 100_500, **one)
    # Every position in 149,500-151,000 is within 1 kb of the N run at 150,001-150,010.
    with pytest.raises(SystemExit):
        draw_events.draw(bam, fasta, [], CHROM, 149_500, 151_000, **one)


def test_a_deletion_needs_coverage_along_its_length(tmp_path, fasta):
    # Every 300 bp deletion starting at 59,801-59,851 has a checked base in the MAPQ 0 stretch at
    # 60,001-60,100, although its first base is fine.
    bam = even_bam(tmp_path)
    assert measure.coverage_ok(bam, CHROM, 59_801) and not measure.coverage_ok(bam, CHROM, 60_001)
    with pytest.raises(SystemExit):
        draw_events.draw(bam, fasta, [], CHROM, 59_801, 60_150, seed=7, n_het=0, n_hom=0, n_del=1,
                         max_candidates=500)


def test_the_draw_stops_when_too_few_are_accepted(tmp_path, fasta):
    bam = even_bam(tmp_path)
    with pytest.raises(SystemExit):
        draw_events.draw(bam, fasta, [], CHROM, 5_000, 195_000, seed=7,
                         n_het=20, n_hom=20, n_del=3, max_candidates=20_000)


def test_the_vcf_states_what_spike_reads(tmp_path, fasta):
    events = [{"kind": "het", "chrom": CHROM, "start": 1000, "end": 1000, "ref": "C", "alt": "T"},
              {"kind": "hom", "chrom": CHROM, "start": 30_000, "end": 30_000, "ref": "A", "alt": "G"},
              {"kind": "del", "chrom": CHROM, "start": 60_001, "end": 60_300}]
    path = tmp_path / "e.vcf"
    draw_events.write_vcf(events, fasta, str(path))
    rows = [l.rstrip("\n").split("\t") for l in path.read_text().splitlines() if not l.startswith("#")]
    assert rows[0][:5] == [CHROM, "1000", ".", "C", "T"] and rows[0][7] == "SIM_VAF=0.5"
    assert rows[1][:5] == [CHROM, "30000", ".", "A", "G"] and rows[1][7] == "SIM_VAF=1.0"
    # spike reads a DEL with END as 0-based (POS, END): POS is the base before the deletion.
    assert rows[2][:5] == [CHROM, "60000", ".", REF[59_999], "<DEL>"]
    assert rows[2][7] == "SVTYPE=DEL;END=60300;SIM_VAF=1.0"


# --- reported rows ---------------------------------------------------------------------

def test_reads_wholly_inside_a_deletion(tmp_path):
    bam = write_bam(tmp_path / "a.bam", [
        segment("inside", 2_000),                                       # 2000..2099
        segment("edge", 1_999),                                         # starts one base early
        segment("dupinside", 2_050, flag=0x1 | 0x2 | 0x40 | 0x400),
    ])
    # The deletion is 1-based 2001-2300 = 0-based [2000, 2300).
    assert measure.inside(bam, CHROM, 2_001, 2_300) == ["inside"]


def test_duplicate_share_near_a_site(tmp_path):
    bam = write_bam(tmp_path / "a.bam", [
        segment("a", 950), segment("b", 950, flag=0x1 | 0x2 | 0x40 | 0x400),
        segment("c", 950, flag=0x1 | 0x2 | 0x40 | 0x100),               # secondary: not counted
        segment("far", 5_000),                                          # outside 700..1300
    ])
    assert measure.dup_counts(bam, CHROM, 1000, flank=300) == (1, 2)


# --- C1 --------------------------------------------------------------------------------

def dup_set(di, rep, dups, pos0, improper=False):
    flag = 0x1 | (0 if improper else 0x2) | 0x40
    return ([segment(rep, pos0, flag=flag, tags=[("DI", di)])]
            + [segment(d, pos0, flag=flag | 0x400, tags=[("DI", di)]) for d in dups])


def test_c1_picks_sets_with_one_kept_pair_and_some_duplicates(tmp_path):
    segs = (dup_set(1, "k1", ["d1a", "d1b"], 1_000)
            + dup_set(2, "k2", ["d2a"], 2_000)
            + dup_set(3, "k3", ["d3a"], 3_000, improper=True)                    # not proper: out
            + [segment("k4", 4_000, tags=[("DI", 4)]), segment("k4b", 4_000, tags=[("DI", 4)])]  # no duplicate
            + [segment("plain", 5_000)])
    bam = write_bam(tmp_path / "a.bam", segs)
    sets = control.pick_sets(bam, CHROM, 0, LENGTH, seed=7, n=2)
    assert sorted((s["rep"], tuple(sorted(s["dups"]))) for s in sets) == [("k1", ("d1a", "d1b")), ("k2", ("d2a",))]
    with pytest.raises(SystemExit):
        control.pick_sets(bam, CHROM, 0, LENGTH, seed=7, n=3)


def test_c1_counts_a_set_freed_when_a_member_lost_its_flag(tmp_path):
    sets = [{"rep": "k1", "dups": ["d1a", "d1b"]}, {"rep": "k2", "dups": ["d2a"]}]
    after = write_bam(tmp_path / "a.bam", [
        segment("d1a", 1_000), segment("d1b", 1_000, flag=0x1 | 0x2 | 0x40 | 0x400),  # set 1 freed
        segment("d2a", 2_000, flag=0x1 | 0x2 | 0x40 | 0x400),                          # set 2 not
        segment("d2a", 2_000, flag=0x1 | 0x2 | 0x80 | 0x100),                          # secondary: not counted
    ])
    assert control.freed(after, CHROM, 0, LENGTH, sets) == 1
