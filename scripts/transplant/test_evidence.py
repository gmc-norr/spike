"""Tests for evidence.py: the read evidence of one deletion, counted by hand.

Run: python3 -m pytest scripts/transplant/test_evidence.py
"""
import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.dirname(__file__))
import evidence  # noqa: E402

CHROM, LENGTH = "chrT", 20_000
START, END = 5_000, 5_500  # the deletion: 0-based, half-open
BOUND = 878


def segment(header, name, pos, cigar, flag=0, mate_pos=None, tlen=0, sa=None):
    s = pysam.AlignedSegment(header)
    s.query_name = name
    s.flag = flag
    s.reference_id = 0
    s.reference_start = pos
    s.mapping_quality = 60
    s.cigarstring = cigar
    qlen = sum(n for op, n in s.cigartuples if op in (0, 1, 4, 7, 8))
    s.query_sequence = "A" * qlen
    s.query_qualities = pysam.qualitystring_to_array("I" * qlen)
    if mate_pos is not None:
        s.next_reference_id = 0
        s.next_reference_start = mate_pos
        s.template_length = tlen
    if sa:
        s.set_tag("SA", sa)
    return s


def write_bam(path, segments, header):
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for s in sorted(segments, key=lambda s: s.reference_start):
            out.write(s)
    pysam.index(path)


@pytest.fixture
def bams(tmp_path):
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                              "SQ": [{"SN": CHROM, "LN": LENGTH}]})
    h = header
    reads = []
    # Background: 100 bp reads every 10 bp over both flanks, every 20 bp inside
    # (the other haplotype), none carrying evidence.
    for i, pos in enumerate(list(range(3_800, 4_901, 10)) + list(range(5_500, 6_601, 10))):
        reads.append(segment(h, f"t{i}", pos, "100M"))
    for i, pos in enumerate(range(5_000, 5_401, 20)):
        reads.append(segment(h, f"in{i}", pos, "100M"))
    # Evidence.
    reads += [
        segment(h, "g1", 4_950, "50M500D50M"),                  # gap at START, 500 >= 250
        segment(h, "c1", 4_900, "100M50S"),                     # clip at START
        segment(h, "c2", 5_500, "30S100M"),                     # clip at END
        segment(h, "s1", 4_900, "100M50S", sa="chrT,5501,+,100S50M,60,0;"),  # split
        segment(h, "d1", 4_500, "100M", flag=0x1 | 0x20 | 0x40, mate_pos=5_800, tlen=1_400),
        segment(h, "d1", 5_800, "100M", flag=0x1 | 0x10 | 0x80, mate_pos=4_500, tlen=-1_400),
    ]
    # Not evidence, or not counted.
    reads += [
        segment(h, "g2", 4_950, "50M100D50M"),                  # gap too short (100 < 250)
        segment(h, "c3", 4_800, "100M50S"),                     # clip 100 bp before START
        segment(h, "c4", 4_900, "100M5S"),                      # clip shorter than 10
        segment(h, "d2", 4_200, "100M", flag=0x1 | 0x20 | 0x40, mate_pos=4_500, tlen=400),
        segment(h, "d2", 4_500, "100M", flag=0x1 | 0x10 | 0x80, mate_pos=4_200, tlen=-400),
        segment(h, "dup", 4_950, "50M500D50M", flag=0x400),
        segment(h, "sec", 4_950, "50M500D50M", flag=0x100),
        segment(h, "qc", 4_950, "50M500D50M", flag=0x200),
        segment(h, "skip1", 4_950, "50M500D50M"),                # named in the skip list
    ]
    main = str(tmp_path / "main.bam")
    write_bam(main, reads, header)
    extra = str(tmp_path / "extra.bam")
    write_bam(extra, [segment(h, "x1", 4_950, "50M500D50M")], header)
    skip = str(tmp_path / "skip.txt")
    with open(skip, "w") as fh:
        fh.write("skip1\n")
    return main, extra, skip


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


def test_one_deletion_counts_each_kind_of_evidence_once_per_pair(bams):
    main, extra, skip = bams
    row = evidence.measure(main, CHROM, START, END, BOUND, skip_names={"skip1"}, extra_bams=[extra])
    assert (row["n_gap"], row["n_clip"], row["n_split"], row["n_disc"]) == (2, 3, 1, 1)
    # g1, x1, c1, c2, s1 (clip and split, one pair), d1.
    assert row["n_any"] == 6
    flanks = [(main, True), (extra, False)]
    left = reference_depth(flanks, {"skip1"}, START - 1_000, START)
    right = reference_depth(flanks, {"skip1"}, END, END + 1_000)
    assert row["flank_depth"] == pytest.approx((left + right) / 2)
    assert row["J"] == pytest.approx(6 / row["flank_depth"])
    inside = reference_depth(flanks, {"skip1"}, START, END)
    assert row["E1"] == pytest.approx(inside / row["flank_depth"])


def test_without_the_skip_list_and_extra_bam_the_counts_change(bams):
    main, _, _ = bams
    row = evidence.measure(main, CHROM, START, END, BOUND)
    # skip1 now counts, x1 is gone.
    assert (row["n_gap"], row["n_any"]) == (2, 6)
    assert row["n_clip"] == 3


def test_e1_is_left_out_under_300_bp(bams):
    main, _, _ = bams
    row = evidence.measure(main, CHROM, START, START + 299, BOUND)
    assert row["E1"] is None
