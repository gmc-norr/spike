"""Tests for small_evidence.py: the read evidence of one SNV, indel or duplication, counted by hand.

Run: python3 -m pytest scripts/transplant/test_small_evidence.py
"""
import os
import random
import sys

import pysam
import pytest

sys.path.insert(0, os.path.dirname(__file__))
import small_evidence as se  # noqa: E402

CHROM, LENGTH = "chrT", 20_000


def make_ref():
    rng = random.Random(7)
    seq = [rng.choice("ACGT") for _ in range(LENGTH)]
    seq[4999:5013] = list("G" + "CA" * 6 + "T")     # a CA repeat at [5000, 5012)
    seq[7999:8004] = list("ACGTA")                   # CGT at [8000, 8003), no repeat
    seq[9000] = "C"                                  # the SNV's base
    seq[14999:16501] = list("C" + "A" * 1500 + "C")  # a 1.5 kb run of A at [15000, 16500)
    return "".join(seq)


REF = make_ref()


def fetch(chrom, start, end):
    assert chrom == CHROM
    return REF[start:end]


HEADER = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                          "SQ": [{"SN": CHROM, "LN": LENGTH}]})


def segment(name, pos, cigar, flag=0, ins="", sub=None, sa=None):
    """A read at `pos` with `cigar`, its bases taken from REF: `ins` fills I ops in turn,
    `sub` replaces reference bases {pos: base}, soft clips read T."""
    s = pysam.AlignedSegment(HEADER)
    s.query_name, s.flag, s.reference_id, s.reference_start, s.mapping_quality = name, flag, 0, pos, 60
    s.cigarstring = cigar
    bases, rpos, ins_left = [], pos, ins
    for op, n in s.cigartuples:
        if op in (0, 7, 8):
            bases += [(sub or {}).get(p, REF[p]) for p in range(rpos, rpos + n)]
            rpos += n
        elif op == 1:
            bases += list(ins_left[:n])
            ins_left = ins_left[n:]
        elif op == 4:
            bases += ["T"] * n
        elif op in (2, 3):
            rpos += n
    s.query_sequence = "".join(bases)
    s.query_qualities = pysam.qualitystring_to_array("I" * len(bases))
    if sa:
        s.set_tag("SA", sa)
    return s


def background(around):
    """100 bp reads every 10 bp over the 1.2 kb on each side of `around` = (lo, hi), none overlapping it."""
    lo, hi = around
    return ([segment(f"bl{i}", p, "100M") for i, p in enumerate(range(lo - 1_200, lo - 100 - 20, 10))]
            + [segment(f"br{i}", p, "100M") for i, p in enumerate(range(hi + 20, hi + 1_200, 10))])


def write_bam(path, segments):
    path = str(path)
    with pysam.AlignmentFile(path, "wb", header=HEADER) as out:
        for s in sorted(segments, key=lambda s: s.reference_start):
            out.write(s)
    pysam.index(path)
    return str(path)


def reference_depth(paths, skip, start, end):
    """Mean depth over [start, end) by pysam's own counter, same read filter."""
    total = 0
    for path, use_skip in paths:
        with pysam.AlignmentFile(path) as bam:
            cov = bam.count_coverage(
                CHROM, start, end, quality_threshold=0,
                read_callback=lambda r, u=use_skip: not (r.flag & 0xF04) and not (u and r.query_name in skip))
        total += sum(sum(base) for base in cov)
    return total / (end - start)


def flank(paths, skip, lo, hi):
    return (reference_depth(paths, skip, lo - 1_000, lo) + reference_depth(paths, skip, hi, hi + 1_000)) / 2


# --- the repeat region ---------------------------------------------------------------

def test_a_repeat_region_is_the_unit_extended_while_the_reference_repeats_it():
    assert se.repeat_region(fetch, CHROM, 5000, 5002, "CA") == (5000, 5012)   # deleting one CA
    assert se.repeat_region(fetch, CHROM, 5000, 5000, "CA") == (5000, 5012)   # inserting one CA
    assert se.repeat_region(fetch, CHROM, 5004, 5004, "CA") == (5000, 5012)   # from inside it
    assert se.repeat_region(fetch, CHROM, 8000, 8003, "CGT") == (8000, 8003)  # no repeat: the event itself
    assert se.repeat_region(fetch, CHROM, 8000, 8000, "GGT") == (8000, 8000)  # an insertion in no repeat: empty
    assert se.repeat_region(fetch, CHROM, 15000, 15001, "A") == (15000, 16500)  # longer than one fetch, rightward
    assert se.repeat_region(fetch, CHROM, 16400, 16400, "A") == (15000, 16500)  # and leftward


# --- SNV -------------------------------------------------------------------------------

def test_an_snv_s_allele_fraction_is_alt_reads_over_reads_with_a_base_there(tmp_path):
    alt = {9000: "T"}
    reads = background((9000, 9001)) + [
        segment("a1", 8950, "100M", sub=alt), segment("a2", 8960, "100M", sub=alt),
        segment("a3", 8990, "100M", sub=alt),
        segment("r1", 8950, "100M"), segment("r2", 8999, "100M"),
        segment("t1", 8950, "100M", sub={9000: "G"}),          # a third base: in the denominator
        segment("d1", 8950, "50M1D49M"),                        # no base at 9000: not counted
        segment("n1", 8900, "100M"),                            # ends at 9000: not over it
        segment("dup", 8950, "100M", flag=0x400, sub=alt),
        segment("skip1", 8950, "100M", sub=alt),
    ]
    main = write_bam(tmp_path / "main.bam", reads)
    extra = write_bam(tmp_path / "extra.bam", [segment("x1", 8950, "100M", sub=alt)])
    row = se.measure(main, fetch, ("SNV", CHROM, 9001, "C", "T"), {"skip1"}, [extra])
    assert (row["n_carry"], row["n_ref"]) == (4, 3)
    assert row["A"] == pytest.approx(4 / 7)
    assert row["E"] is None and row["J"] is None


# --- DEL -------------------------------------------------------------------------------

def test_a_deletion_s_exact_form_fraction_and_any_form_evidence(tmp_path):
    # Deleting one CA (VCF POS 5000, anchor G at 0-based 4999): repeat region [5000, 5012).
    reads = background((5000, 5012)) + [
        segment("k1", 4950, "50M2D50M"),                  # D2 at 5000: carrier
        segment("k2", 4960, "50M2D50M"),                  # D2 at 5010, in the repeat: carrier
        segment("k3", 4955, "45M2D55M"),                  # a pair, both mates carriers:
        segment("k3", 4965, "35M2D65M"),                  #   two reads for A, one pair for E
        segment("r1", 4950, "100M"),                      # over [4990, 5022), no indel: reference
        segment("r3", 4940, "100M2D20M"),                 # its D2 at 5040 is past 5022: reference
        segment("r2", 4920, "100M"),                      # ends at 5020, short of 5022: neither
        segment("o1", 4950, "50M4D50M"),                  # wrong length: E only
        segment("o2", 4950, "50M2I50M", ins="CA"),        # wrong type: neither
        segment("c1", 4960, "40M60S"),                    # clip at 5000: E only
        segment("c2", 4960, "40M3S"),                     # clip under 5 bp: neither
        segment("f2", 4900, "118M2D30M"),                 # D2 at 5018: past the region +-1, inside +-10: E only
        segment("k5", 4949, "50M2D50M"),                  # D2 at 4999, the region's start - 1: carrier
        segment("f3", 4948, "50M2D50M"),                  # D2 at 4998, - 2: E only
        segment("h1", 4950, "50M1D50M"),                  # D1: half the length is enough for E
        segment("e1", 4990, "3S100M"),                    # starts at 4990 with a 3 bp clip: neither
        segment("dup", 4950, "50M2D50M", flag=0x400),
        segment("sec", 4950, "50M2D50M", flag=0x100),
        segment("qc", 4950, "50M2D50M", flag=0x200),
        segment("skip1", 4950, "50M2D50M"),
    ]
    main = write_bam(tmp_path / "main.bam", reads)
    extra = write_bam(tmp_path / "extra.bam", [segment("x1", 4950, "50M2D50M")])
    row = se.measure(main, fetch, ("DEL", CHROM, 5000, "GCA", "G"), {"skip1"}, [extra])
    assert (row["n_carry"], row["n_ref"]) == (6, 2)     # k1 k2 k3 k3 x1 k5 ; r1 r3
    assert row["A"] == pytest.approx(6 / 8)
    assert row["n_any"] == 10                           # k1 k2 k3 x1 k5 o1 c1 f2 f3 h1
    paths = [(main, True), (extra, False)]
    assert row["flank_depth"] == pytest.approx(flank(paths, {"skip1"}, 5000, 5012))
    assert row["E"] == pytest.approx(10 / row["flank_depth"])
    assert row["J"] is None


def test_without_the_skip_list_and_extra_bam_the_deletion_counts_change(tmp_path):
    reads = background((5000, 5012)) + [segment("k1", 4950, "50M2D50M"), segment("r1", 4950, "100M"),
                                        segment("skip1", 4950, "50M2D50M")]
    main = write_bam(tmp_path / "main.bam", reads)
    row = se.measure(main, fetch, ("DEL", CHROM, 5000, "GCA", "G"))
    assert (row["n_carry"], row["n_ref"], row["n_any"]) == (2, 1, 2)


# --- INS -------------------------------------------------------------------------------

def test_an_insertion_carrier_must_spell_the_truth_s_sequence(tmp_path):
    # Inserting one CA after the G at 0-based 4999: repeat region [5000, 5012).
    reads = background((5000, 5012)) + [
        segment("k1", 4950, "50M2I50M", ins="CA"),        # CA at 5000: carrier
        segment("k2", 4951, "49M2I50M", ins="AC"),        # AC at 5000 reads G AC CA..: another sequence, E only
        segment("k2b", 4951, "50M2I50M", ins="AC"),       # AC at 5001 reads G C AC A.., the same: carrier
        segment("k4", 4950, "50M2I50M", ins="GG"),        # wrong bases: E only
        segment("o1", 4950, "50M4I50M", ins="CACA"),      # wrong length: E only
        segment("r1", 4950, "100M"),                      # reference
    ]
    main = write_bam(tmp_path / "main.bam", reads)
    row = se.measure(main, fetch, ("INS", CHROM, 5000, "G", "GCA"))
    assert (row["n_carry"], row["n_ref"]) == (2, 1)     # k1 k2b ; r1
    assert row["A"] == pytest.approx(2 / 3)
    assert row["n_any"] == 5                            # k1 k2 k2b k4 o1


# --- DUP -------------------------------------------------------------------------------

def test_a_duplication_s_segment_is_the_copy_the_insertion_repeats():
    seg = REF[12000:12100]
    assert se.segment_of(fetch, CHROM, 12000, REF[11999], REF[11999] + seg) == (12000, 12100)  # copy to the right
    assert se.segment_of(fetch, CHROM, 12100, REF[12099], REF[12099] + seg) == (12000, 12100)  # copy to the left
    assert se.segment_of(fetch, CHROM, 12000, REF[11999], REF[11999] + "A" * 100) is None


def test_a_duplication_s_junction_evidence_counts_each_kind_once_per_pair(tmp_path):
    s, e = 12000, 12100
    seg = REF[s:e]
    reads = background((s, e)) + [
        segment("i1", 11950, "100M100I50M", ins=seg),     # a 100 bp I at 12050
        segment("i2", 11950, "100M40I50M", ins=seg),      # 40 bp: under half, not counted
        segment("i3", 11950, "100M60I50M", ins=seg),      # 60 bp: over half, counted
        segment("i4", 11950, "175M100I10M", ins=seg),     # at 12125, 25 bp past e: not counted
        segment("c1", 12000, "100M50S"),                  # clip at e
        segment("c2", 12000, "30S100M"),                  # clip at s
        segment("c3", 12200, "30S100M"),                  # clip far off: not counted
        segment("s1", 11960, "140M50H", sa="chrT,12001,+,140H50M,60,0;"),   # ends at e, rest at s
        segment("s2", 12000, "50H100M", sa="chrT,11961,+,140M50H,60,0;"),   # starts at s, rest ends at e
        segment("s3", 11960, "140M50H", sa="chrT,15001,+,140H50M,60,0;"),   # rest far off
    ]
    main = write_bam(tmp_path / "main.bam", reads)
    row = se.measure(main, fetch, ("DUP", CHROM, 12000, REF[11999], REF[11999] + seg))
    assert row["n_any"] == 6                            # i1 i3 c1 c2 s1 s2
    assert row["flank_depth"] == pytest.approx(flank([(main, False)], set(), s, e))
    assert row["J"] == pytest.approx(6 / row["flank_depth"])
    assert row["A"] is None and row["E"] is None


def test_the_count_n1_reads_is_carriers_for_small_variants_and_pairs_for_duplications():
    assert se.recipient_count({"n_carry": 3, "n_any": 9}, "SNV") == 3
    assert se.recipient_count({"n_carry": 3, "n_any": 9}, "INS") == 3
    assert se.recipient_count({"n_carry": None, "n_any": 9}, "DUP") == 9
