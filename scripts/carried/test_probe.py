"""Tests for probe.py's site and read rules, on hand-made reads."""
import pysam
import pytest

from probe import classify, decide, trim

#            0123456789012345678901234567890
REF = 'ACGTACGTACGTACGTACGTACGTACGTACGT'


def ref(p):
    return REF[p]


def read(start, cigar, seq):
    r = pysam.AlignedSegment()
    r.reference_start = start
    r.cigarstring = cigar
    r.query_sequence = seq
    return r


@pytest.mark.parametrize('pos, r, a, want', [
    (10, 'G', 'T', (10, 11)),            # SNV
    (10, 'GT', 'CA', (10, 12)),          # MNV
    (9, 'CGT', 'C', (10, 12)),           # deletion: the anchor goes
    (9, 'C', 'CAA', (10, 10)),           # insertion: between 9 and 10
    (9, 'CGTA', 'CTA', (10, 11)),        # shared suffix too
])
def test_trim_keeps_only_the_changed_bases(pos, r, a, want):
    assert trim(pos, r, a) == want


def test_a_matching_read_covers_and_does_not_carry():
    assert classify(read(0, '20M', REF[:20]), 10, 11, ref) == (True, False)


def test_a_mismatch_counts_inside_the_site_only():
    inside = REF[:10] + 'T' + REF[11:20]
    outside = REF[:12] + 'C' + REF[13:20]
    assert outside != REF[:20]
    assert classify(read(0, '20M', inside), 10, 11, ref) == (True, True)
    assert classify(read(0, '20M', outside), 10, 11, ref) == (True, False)


def test_an_n_base_is_not_another_allele():
    assert classify(read(0, '20M', REF[:10] + 'N' + REF[11:20]), 10, 11, ref) == (True, False)


def test_a_deletion_over_the_site_carries():
    assert classify(read(0, '10M2D8M', REF[:10] + REF[12:20]), 10, 12, ref) == (True, True)
    # one that ends before the site does not
    assert classify(read(0, '7M2D11M', REF[:7] + REF[9:20]), 10, 12, ref) == (True, False)


@pytest.mark.parametrize('cigar, want', [
    ('10M2I8M', True),   # boundary 10 = s
    ('11M2I7M', True),   # boundary 11 = e
    ('12M2I6M', False),  # boundary 12: past the site
])
def test_an_insertion_carries_at_the_site_boundaries_only(cigar, want):
    m = int(cigar.split('M')[0])
    seq = REF[:m] + 'AA' + REF[m:18]
    assert classify(read(0, cigar, seq), 10, 11, ref) == (True, want)


def test_a_read_that_stops_inside_the_site_does_not_cover():
    assert classify(read(0, '11M', REF[:11]), 10, 12, ref) == (False, False)
    assert classify(read(11, '10M', REF[11:21]), 10, 12, ref) == (False, False)


def test_decide_needs_ten_reads_and_a_fifth():
    assert decide(9, 9) == 'unchecked'
    assert decide(10, 2) == 'refuse'
    assert decide(10, 1) == 'pass'
