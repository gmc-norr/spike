# `--edit-model origin` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `--edit-model origin`, which removes each original read by its chance of having come from an event's edited copy (from MAPQ and `XA`), at the event and at its look-alikes, and then tests it on a made-up genome with an exact twin.

**Architecture:** A new module `src/origin.rs` holds the model. It covers placements and chances, fragments, look-alike regions, origin depth, removal decisions and the BAM/CRAM scan. `simulate.rs` gains an origin path: it keeps every pool pair, scales the tiling by origin depth, and returns per-fragment removal chances. `main.rs` gathers a site per event and, after all events, draws once per fragment family. It then removes those pairs from the pools' outputs and lists every removed name for `merge.sh`. `clean` stays the default and byte-identical.

**Tech Stack:** Rust 2021, noodles 0.88 (noodles-sam 0.69, noodles-bam 0.73, noodles-cram 0.74), rand 0.8. The physics test uses Python 3, wgsim, bwa-mem2 and samtools.

**Spec:** `docs/superpowers/specs/2026-09-26-edit-model-origin-design.md` (branch `edit-model`; fixes R1–R5). Numbers stated there live there; this plan points at them.

## Global Constraints

- The default is `--edit-model clean`. Its output must stay byte-identical to master `94f3f20`: the same `R1.fq.gz`/`R2.fq.gz` contents, `replaced_reads.txt` and `truth.vcf` body.
- `origin` applies only to events that remove reads: DEL, INV, INS, small variants and the full DUP model. Additive events (fusion, `--dup-model junction`) run as `clean`.
- Never push. Never merge without the user's word. Never touch the user's ` M .gitignore` or `?? .ignore` in `/home/parlar_ai/dev/spike`. Always `git add` named files, never `-A` or `.`.
- TDD: every new test is seen failing for the right reason before its code is written.
- The judgment-gate order is: the plan commit (this file), then `code:` commits, then one `result:` commit. The plan and the result are never in one commit.
- Build in its own `CARGO_TARGET_DIR` per commit under the scratchpad, and md5 each binary that is run.
- Use at most 16 threads (`-j 16`, `--test-threads=16`, align/merge threads ≤ 16).
- Commit trailer: `Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>`.
- Use `LC_ALL=C` for any awk that prints decimals. Label any number not yet run as `PREDICTED (not run)`.
- No new clippy warnings. The baseline, measured on `94f3f20`, is `bin "spike" generated 13 warnings` and `bin "spike" test generated 14 warnings (12 duplicates)`.
- `S=/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad` below. Use a new scratch dir if the session changed.

## Measured before planning

- **The 35x HG002 BAM keeps `XA` tags** (`data/giab_hg38/HG002/HG002.novaseq.pcr-free.35x.bwamem2.dedup.grch38_no_alt.bam`):
  - Its first 100,000 records hold 15,255 with `XA`, the first at record 64.
  - `chr20:10000000-11000000` holds 2,841 of 298,950 with `XA`.
  - The Task 9 slice `chr20:14500000-14600000` holds 328 of 29,427, the first at record 498.
  - Some spots have MAPQ 0 reads without `XA`. At `chr20:7117236-7121236`, 822 are MAPQ 0 and 1 carries `XA`. In Task 9's footprint `chr20:14546421-14550735`, 4 are MAPQ 0 and none carries it, while 1 MAPQ 40 read does. Under bwa-mem's `-h 5` rule those reads have more than 5 hits.
  - This plan's first revision (`84ca4c8`) said the BAM had no `XA`, from one window. That was wrong. The input check is now on the file (spec R8).
- **Tools found:** `wgsim`, `bwa-mem2`, `bwa` and `samtools` are in `/home/parlar_ai/.pixi/bin`, and python has `pysam 0.23.3`.
- **The tests pass on `94f3f20`:** `cargo test`: 554 passed, 0 failed, 1 ignored.

## Revision 2 (2026-09-27, before any run)

A review of the first revision (`84ca4c8`) found five gaps. Each was checked against the code or the BAM before this revision. Nothing had been built or run, so the gate's order holds: this revision is itself a `plan:` commit.

| | Gap | Fix | Where |
| --- | --- | --- | --- |
| P1 | origin depth counted fragments R4 keeps (spec R6) | depth over removable fragments only | Task 5 |
| P2 | phasing skips duplicates, so a family could split on its shared draw (spec R7) | the family's phase call, one fate at its highest member's total | Task 6 |
| P3 | the gate counted a crashed control as evidence | NO VERDICT unless `origin` and `--min-mapq 0` finish and `clean` finishes or stops with its expected refusal | Task 12 |
| P4 | the XA check passed on any record with `XA`; a per-spot check cannot tell a stripped tag from more than 5 hits (spec R8) | `require_xa` on the file's first 100,000 records | Tasks 7, 9 |
| P5 | M12's test put every pool pair outside the footprint, so it could not catch early suppression | a test with pool pairs inside it | Tasks 8, 11 |

## File structure

| File | Change | Responsibility |
| --- | --- | --- |
| `src/origin.rs` | create | the origin model: `Span`, `Placement`, chances, `OriginRecord`, `Fragment`, look-alike regions, `OriginSite` (depth, `f`, removal chances), `decide`, and the BAM/CRAM `gather` |
| `src/types.rs` | modify `SplicedOutput` | carry `origin_chances` |
| `src/simulate.rs` | modify | `is_additive`, `simulate_event_origin`, the origin path in the event body, `origin_coverage_for_tiling`, `depth_fold_by`, `apply_removals` |
| `src/main.rs` | modify | `mod origin`, `--edit-model`, `validate_edit_model`, `origin_footprint`, the per-event site, one `decide` after all events |
| `src/extract.rs` | modify `test_fixtures` | `write_one_contig_bam` for origin's BAM tests |
| `scripts/origin_physics.py` | create | the physics test (Task 12) |
| `README.md`, `docs/review/REVIEW.md` | modify | user docs; the physics result |

---

### Task 1: The `--edit-model` flag

**Files:**
- Modify: `src/main.rs`. Add the `Args` field after `dup_model` (`src/main.rs:190-195`), validate it after the `--dup-model` check (`src/main.rs:426-432`), and add tests in `mod tests`.

**Interfaces:**
- Produces: `Args::edit_model: String` (`"clean"` or `"origin"`), `fn validate_edit_model(model: &str) -> Result<()>`.

- [ ] **Step 1: Write the failing tests** (in `src/main.rs` `mod tests`)

```rust
    #[test]
    fn test_validate_edit_model_accepts_clean_and_origin() {
        assert!(validate_edit_model("clean").is_ok());
        assert!(validate_edit_model("origin").is_ok());
    }

    #[test]
    fn test_validate_edit_model_rejects_anything_else_and_names_both() {
        let err = validate_edit_model("Origin").unwrap_err().to_string();
        assert!(err.contains("'Origin'") && err.contains("clean") && err.contains("origin"), "{}", err);
    }

    #[test]
    fn test_edit_model_defaults_to_clean() {
        let args = Args::try_parse_from(["spike", "--bam", "x.bam", "--reference", "x.fa"]).unwrap();
        assert_eq!(args.edit_model, "clean");
    }
```

- [ ] **Step 2: Run them and see them fail**

Run: `CARGO_TARGET_DIR=$S/target-t1 cargo test -j 16 edit_model -- --test-threads=16`
Expected: a compile error, `cannot find function validate_edit_model` and `no field edit_model`.

- [ ] **Step 3: Implement**

Add to `Args`, after `dup_model`:

```rust
    /// Which original reads an event replaces. "clean" (default): only the
    /// donor pool's pairs (both mates at --min-mapq or above, a proper pair,
    /// no duplicate, secondary, supplementary or QC-fail flag) inside the
    /// event's footprint. "origin" (experimental): every primary read at the
    /// event and at its look-alikes, each removed by its chance of having
    /// come from the edited copy, read from its MAPQ and its XA tag. It
    /// needs the aligner's XA tags (bwa-mem and bwa-mem2 write them).
    #[arg(long, default_value = "clean")]
    edit_model: String,
```

Add beside `validate_allele_fraction`:

```rust
/// Check `--edit-model`: "clean" or "origin".
fn validate_edit_model(model: &str) -> Result<()> {
    if model != "clean" && model != "origin" {
        bail!("invalid --edit-model '{}', expected 'clean' or 'origin'", model);
    }
    Ok(())
}
```

In `main`, right after the `--dup-model` check:

```rust
    validate_edit_model(&args.edit_model)?;
```

- [ ] **Step 4: Run the tests and see them pass**

Run: `CARGO_TARGET_DIR=$S/target-t1 cargo test -j 16 edit_model -- --test-threads=16`
Expected: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add src/main.rs
git commit -m "code: origin -- the --edit-model flag (clean default, origin accepted, no behavior yet)

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Placements and chances

**Files:**
- Create: `src/origin.rs`
- Modify: `src/main.rs`. Add `mod origin;` after `mod loh;` (`src/main.rs:12`).

**Interfaces:**
- Produces: `pub struct Span { chrom: String, start: u64, end: u64 }` with `new`, `holds`, `overlaps` and `Display`; `pub struct Placement { span: Span, chance: f64 }`; `pub fn p_here(mapq: u8, n_xa: usize) -> f64`; `pub fn parse_xa(xa: &str) -> Vec<Span>`; `pub fn placements(primary: Span, mapq: u8, alternatives: &[Span]) -> Vec<Placement>`; `pub fn chance_within(placements: &[Placement], footprint: &Span) -> f64`; `fn reference_length(ops: &[Op]) -> u64`.

- [ ] **Step 1: Write the failing tests.** Create `src/origin.rs` holding only the module doc, the imports and this test module:

```rust
//! `--edit-model origin`: which original reads came from an event's edited
//! copy, wherever the aligner put them.
//!
//! The default (`clean`) edits only donor-pool pairs: both mates at
//! `--min-mapq` or above, a proper pair, no duplicate, secondary,
//! supplementary or QC-fail flag. Where the aligner cannot tell a spot from a
//! look-alike most reads fail that, so the event is barely planted (RF8).
//! `origin` gives every primary read a chance of having come from the
//! event's footprint, read from its MAPQ and its `XA` hits, and removes it
//! by that chance. The design is
//! `docs/superpowers/specs/2026-09-26-edit-model-origin-design.md`; comments
//! name its fixes R1-R5.

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn test_p_here_follows_mapq_and_the_xa_count() {
        assert!(close(p_here(60, 0), 1.0 - 1e-6));
        assert!(close(p_here(3, 0), 1.0 - 10f64.powf(-0.3)));
        assert!(close(p_here(0, 1), 0.5));
        // No XA at MAPQ 0: bwa-mem lists at most 5, so there are more.
        assert!(close(p_here(0, 0), 1.0 / 6.0));
    }

    #[test]
    fn test_parse_xa_reads_chrom_start_strand_and_cigar_length() {
        let hits = parse_xa("cluster,+2002,151M,0;cluster,-1002,100M1D51M,1;");
        assert_eq!(
            hits,
            vec![Span::new("cluster", 2001, 2152), Span::new("cluster", 1001, 1153)]
        );
    }

    #[test]
    fn test_parse_xa_skips_a_malformed_entry() {
        assert_eq!(parse_xa("chr1,+10,5Q,0;chr1,+20,5M,0;"), vec![Span::new("chr1", 19, 24)]);
    }

    #[test]
    fn test_a_read_on_a_remote_copy_with_three_hits_in_the_footprint_gets_three_quarters() {
        // R2: bwa-mem put `cluster_1` on the remote copy at MAPQ 0 with all
        // three XA hits in one cluster (the review's probe).
        let alternatives = parse_xa("cluster,+2002,151M,0;cluster,+1002,151M,0;cluster,+1502,151M,0;");
        let placed = placements(Span::new("remote", 1001, 1152), 0, &alternatives);
        assert!(close(chance_within(&placed, &Span::new("cluster", 1000, 2400)), 0.75));
    }

    #[test]
    fn test_a_mapq0_read_whose_one_hit_is_also_in_the_footprint_gets_one() {
        let placed = placements(Span::new("chr1", 100, 250), 0, &[Span::new("chr1", 600, 750)]);
        assert!(close(chance_within(&placed, &Span::new("chr1", 0, 1000)), 1.0));
    }

    #[test]
    fn test_a_placement_sticking_out_of_the_footprint_does_not_count() {
        let placed = placements(Span::new("chr1", 950, 1100), 60, &[]);
        assert!(close(chance_within(&placed, &Span::new("chr1", 0, 1000)), 0.0));
    }
}
```

Add `mod origin;` to `src/main.rs`.

- [ ] **Step 2: Run them and see them fail**

Run: `CARGO_TARGET_DIR=$S/target-t2 cargo test -j 16 origin:: -- --test-threads=16`
Expected: a compile error, `cannot find function p_here` (and `Span`, `parse_xa`, `placements`, `chance_within`).

- [ ] **Step 3: Implement.** Put this above `#[cfg(test)]` in `src/origin.rs`:

```rust
use noodles::sam::alignment::record::cigar::Op;

/// A reference interval, 0-based half-open.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Span {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl Span {
    pub fn new(chrom: &str, start: u64, end: u64) -> Self {
        Span {
            chrom: chrom.to_string(),
            start,
            end,
        }
    }

    /// Whether `other` lies wholly inside this span.
    pub fn holds(&self, other: &Span) -> bool {
        self.chrom == other.chrom && other.start >= self.start && other.end <= self.end
    }

    /// Whether `other` shares at least one base with this span.
    pub fn overlaps(&self, other: &Span) -> bool {
        self.chrom == other.chrom && other.start < self.end && self.start < other.end
    }
}

impl std::fmt::Display for Span {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
    }
}

/// One place a read may have come from, and the chance that it did.
#[derive(Debug, Clone, PartialEq)]
pub struct Placement {
    pub span: Span,
    pub chance: f64,
}

/// bwa-mem lists alternative hits in `XA` only when there are at most this
/// many (its `-h 5` default), so a MAPQ 0 read without `XA` has more.
const MAX_LISTED_HITS: usize = 5;

/// The chance a read came from where the aligner put it.
///
/// MAPQ 1 or more is the aligner's own estimate, `1 - 10^(-MAPQ/10)`. At
/// MAPQ 0 the hits are equally good: with `n_xa` listed alternatives the
/// primary is one of `1 + n_xa`; with none listed there are more than
/// [`MAX_LISTED_HITS`], so it is one of at least 6.
pub fn p_here(mapq: u8, n_xa: usize) -> f64 {
    if mapq > 0 {
        1.0 - 10f64.powf(-(mapq as f64) / 10.0)
    } else if n_xa > 0 {
        1.0 / (1 + n_xa) as f64
    } else {
        1.0 / (1 + MAX_LISTED_HITS) as f64
    }
}

/// The reference bases a CIGAR covers.
fn reference_length(ops: &[Op]) -> u64 {
    ops.iter()
        .filter(|op| op.kind().consumes_reference())
        .map(|op| op.len() as u64)
        .sum()
}

/// The reference bases a CIGAR string such as `100M1D51M` covers, or `None`
/// when it is malformed.
fn cigar_reference_length(cigar: &str) -> Option<u64> {
    let (mut len, mut n, mut digits) = (0u64, 0u64, false);
    for c in cigar.bytes() {
        if c.is_ascii_digit() {
            n = n * 10 + u64::from(c - b'0');
            digits = true;
            continue;
        }
        if !digits {
            return None;
        }
        match c {
            b'M' | b'D' | b'N' | b'=' | b'X' => len += n,
            b'I' | b'S' | b'H' | b'P' => {}
            _ => return None,
        }
        n = 0;
        digits = false;
    }
    (!digits).then_some(len)
}

/// The alternative hits in an `XA:Z` value: `chrom,±pos,CIGAR,NM;` each, with
/// `pos` 1-based and its sign the strand. Malformed entries are skipped.
pub fn parse_xa(xa: &str) -> Vec<Span> {
    xa.split(';')
        .filter(|entry| !entry.is_empty())
        .filter_map(|entry| {
            let mut fields = entry.split(',');
            let chrom = fields.next()?;
            let pos: i64 = fields.next()?.parse().ok()?;
            let len = cigar_reference_length(fields.next()?)?;
            let start = pos.unsigned_abs().checked_sub(1)?;
            Some(Span::new(chrom, start, start + len))
        })
        .collect()
}

/// A read's placements: its primary at [`p_here`], then each `XA` hit at an
/// equal share of the rest, `(1 - p_here) / k`.
pub fn placements(primary: Span, mapq: u8, alternatives: &[Span]) -> Vec<Placement> {
    let here = p_here(mapq, alternatives.len());
    let each = if alternatives.is_empty() {
        0.0
    } else {
        (1.0 - here) / alternatives.len() as f64
    };
    std::iter::once(Placement {
        span: primary,
        chance: here,
    })
    .chain(alternatives.iter().map(|span| Placement {
        span: span.clone(),
        chance: each,
    }))
    .collect()
}

/// The chance a read came from inside `footprint`: the sum over its
/// placements lying wholly inside it, the primary included (R2).
pub fn chance_within(placements: &[Placement], footprint: &Span) -> f64 {
    placements
        .iter()
        .filter(|p| footprint.holds(&p.span))
        .map(|p| p.chance)
        .sum()
}
```

Until Task 9 wires the module into `main`, the non-test build warns that these items are unused (`dead_code`). That is expected. The clippy baseline is checked in Task 9, when every item has a caller.

- [ ] **Step 4: Run the tests and see them pass**

Run: `CARGO_TARGET_DIR=$S/target-t2 cargo test -j 16 origin:: -- --test-threads=16`
Expected: 6 passed.

- [ ] **Step 5: Commit**

```bash
git add src/origin.rs src/main.rs
git commit -m "code: origin -- placements and their chances (MAPQ, XA; R2 sums a footprint's placements)

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Records and fragments

**Files:**
- Modify: `src/origin.rs`

**Interfaces:**
- Consumes: `Span`, `Placement`, `chance_within`, `reference_length` (Task 2).
- Produces: `pub struct FivePrime { chrom, pos: u64, reverse: bool }` (Ord); `pub fn five_prime(chrom: &str, start: u64, ops: &[Op], reverse: bool) -> FivePrime`; `pub struct OriginRecord { name, first, placements, duplicate, qc_fail, mate_unmapped, five_prime }` with `fn primary(&self) -> &Placement`; `pub struct Fragment<'a> { name: &'a str, mates: Vec<&'a OriginRecord> }` with `chance(&self, &Span) -> f64`, `removable(&self, &Span) -> bool`, `at_spot(&self, &Span) -> bool`, `family(&self) -> Vec<FivePrime>`.

- [ ] **Step 1: Write the failing tests** (add to `origin::tests`)

```rust
    use noodles::sam::alignment::record::cigar::op::Kind;

    /// A record placed at `start..start+100` with `alternatives` as its XA.
    pub(super) fn record(name: &str, first: bool, start: u64, mapq: u8, alternatives: &[Span]) -> OriginRecord {
        OriginRecord {
            name: name.to_string(),
            first,
            placements: placements(Span::new("chr1", start, start + 100), mapq, alternatives),
            duplicate: false,
            qc_fail: false,
            mate_unmapped: false,
            five_prime: FivePrime { chrom: "chr1".to_string(), pos: start, reverse: !first },
        }
    }

    fn fp() -> Span {
        Span::new("chr1", 0, 1000)
    }

    #[test]
    fn test_five_prime_steps_back_over_a_leading_clip_or_past_a_trailing_one() {
        let fwd = [Op::new(Kind::SoftClip, 5), Op::new(Kind::Match, 95)];
        assert_eq!(five_prime("chr1", 100, &fwd, false).pos, 95);
        let rev = [Op::new(Kind::Match, 95), Op::new(Kind::SoftClip, 5)];
        assert_eq!(five_prime("chr1", 300, &rev, true).pos, 400);
    }

    #[test]
    fn test_a_pair_takes_its_surest_mates_chance() {
        let a = record("p", true, 100, 60, &[]);
        let b = record("p", false, 300, 0, &[Span::new("chr1", 5000, 5100)]);
        let f = Fragment { name: "p", mates: vec![&a, &b] };
        assert!(close(f.chance(&fp()), 1.0 - 1e-6));
    }

    #[test]
    fn test_a_tie_on_the_primary_goes_to_the_higher_footprint_chance() {
        // Both MAPQ 0 with one hit (1/2 each). a's hit is outside (1/2 in the
        // footprint); b's is inside too (1).
        let a = record("p", true, 100, 0, &[Span::new("chr1", 5000, 5100)]);
        let b = record("p", false, 300, 0, &[Span::new("chr1", 600, 700)]);
        let f = Fragment { name: "p", mates: vec![&a, &b] };
        assert!(close(f.chance(&fp()), 1.0));
    }

    #[test]
    fn test_a_pair_is_removable_only_when_every_mate_could_come_from_the_footprint() {
        // R4.
        let a = record("p", true, 100, 60, &[]);
        let inside = record("p", false, 300, 60, &[]);
        let outside = record("p", false, 1200, 60, &[]);
        let outside_with_hit_inside = record("p", false, 1200, 0, &[Span::new("chr1", 700, 800)]);
        assert!(Fragment { name: "p", mates: vec![&a, &inside] }.removable(&fp()));
        assert!(!Fragment { name: "p", mates: vec![&a, &outside] }.removable(&fp()));
        assert!(Fragment { name: "p", mates: vec![&a, &outside_with_hit_inside] }.removable(&fp()));
    }

    #[test]
    fn test_an_unseen_mate_blocks_removal_unless_it_is_unmapped() {
        let alone = record("p", true, 100, 60, &[]);
        assert!(!Fragment { name: "p", mates: vec![&alone] }.removable(&fp()));
        let orphan = OriginRecord { mate_unmapped: true, ..record("p", true, 100, 60, &[]) };
        assert!(Fragment { name: "p", mates: vec![&orphan] }.removable(&fp()));
    }

    #[test]
    fn test_at_spot_means_a_mate_the_aligner_put_inside_the_footprint() {
        let here = record("p", true, 100, 0, &[Span::new("chr1", 5000, 5100)]);
        let there = record("q", true, 5000, 0, &[Span::new("chr1", 100, 200)]);
        assert!(Fragment { name: "p", mates: vec![&here] }.at_spot(&fp()));
        assert!(!Fragment { name: "q", mates: vec![&there] }.at_spot(&fp()));
    }

    #[test]
    fn test_duplicates_share_a_family_and_other_fragments_do_not() {
        let (a1, a2) = (record("a", true, 100, 60, &[]), record("a", false, 300, 60, &[]));
        let (d1, d2) = (record("d", true, 100, 60, &[]), record("d", false, 300, 60, &[]));
        let (o1, o2) = (record("o", true, 110, 60, &[]), record("o", false, 300, 60, &[]));
        let fam = |x: &OriginRecord, y: &OriginRecord| Fragment { name: "x", mates: vec![x, y] }.family();
        assert_eq!(fam(&a1, &a2), fam(&d2, &d1));
        assert_ne!(fam(&a1, &a2), fam(&o1, &o2));
    }
```

- [ ] **Step 2: Run them and see them fail**

Run: `CARGO_TARGET_DIR=$S/target-t3 cargo test -j 16 origin:: -- --test-threads=16`
Expected: a compile error, `cannot find struct OriginRecord` (and `FivePrime`, `five_prime`, `Fragment`).

- [ ] **Step 3: Implement** (in `src/origin.rs`, after `chance_within`)

```rust
use noodles::sam::alignment::record::cigar::op::Kind;

/// A read's unclipped 5' end: where its first sequenced base would align.
/// Duplicate marking keys a fragment on its two mates' 5' ends.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FivePrime {
    pub chrom: String,
    pub pos: u64,
    pub reverse: bool,
}

/// Clipped bases at the start and at the end of a CIGAR.
fn clips(ops: &[Op]) -> (u64, u64) {
    let clipped = |op: &&Op| matches!(op.kind(), Kind::SoftClip | Kind::HardClip);
    let lead = ops.iter().take_while(clipped).map(|op| op.len() as u64).sum();
    let trail = ops.iter().rev().take_while(clipped).map(|op| op.len() as u64).sum();
    (lead, trail)
}

/// The unclipped 5' end of a read aligned at `start` (0-based) with `ops`.
pub fn five_prime(chrom: &str, start: u64, ops: &[Op], reverse: bool) -> FivePrime {
    let (lead, trail) = clips(ops);
    let pos = if reverse {
        start + reference_length(ops) + trail
    } else {
        start.saturating_sub(lead)
    };
    FivePrime {
        chrom: chrom.to_string(),
        pos,
        reverse,
    }
}

/// One primary record, as far as `origin` needs it.
#[derive(Debug, Clone)]
pub struct OriginRecord {
    pub name: String,
    /// Read 1 of its pair; the two records of one name differ here.
    pub first: bool,
    /// The primary placement first, then the `XA` hits.
    pub placements: Vec<Placement>,
    pub duplicate: bool,
    pub qc_fail: bool,
    pub mate_unmapped: bool,
    pub five_prime: FivePrime,
}

impl OriginRecord {
    /// Where the aligner put the read.
    pub fn primary(&self) -> &Placement {
        &self.placements[0]
    }
}

/// The records of one read name that `origin` read: one or both mates.
#[derive(Debug, Clone)]
pub struct Fragment<'a> {
    pub name: &'a str,
    pub mates: Vec<&'a OriginRecord>,
}

impl Fragment<'_> {
    /// The mate whose own placement is surest: the highest primary chance,
    /// ties to the higher chance of having come from `footprint`. One mate
    /// pinned uniquely somewhere pins the fragment there.
    fn surest(&self, footprint: &Span) -> &OriginRecord {
        self.mates
            .iter()
            .copied()
            .max_by(|a, b| {
                a.primary().chance.total_cmp(&b.primary().chance).then(
                    chance_within(&a.placements, footprint)
                        .total_cmp(&chance_within(&b.placements, footprint)),
                )
            })
            .expect("a fragment has at least one mate")
    }

    /// `p_origin`: the chance this fragment came from `footprint`, its
    /// surest mate's.
    pub fn chance(&self, footprint: &Span) -> f64 {
        chance_within(&self.surest(footprint).placements, footprint)
    }

    /// Whether spike's new reads can replace this fragment (R4): every mate
    /// has a placement wholly inside `footprint`, and a mate `origin` did not
    /// read is unmapped. The new fragments never reach past the footprint,
    /// so removing one that sticks out would leave a depth dip.
    pub fn removable(&self, footprint: &Span) -> bool {
        let every_mate_read = self.mates.len() == 2 || self.mates.iter().all(|m| m.mate_unmapped);
        every_mate_read
            && self
                .mates
                .iter()
                .all(|m| m.placements.iter().any(|p| footprint.holds(&p.span)))
    }

    /// Whether the aligner put a mate inside `footprint`, where the sample's
    /// phase call applies. Elsewhere the copy is unknown.
    pub fn at_spot(&self, footprint: &Span) -> bool {
        self.mates.iter().any(|m| footprint.holds(&m.primary().span))
    }

    /// The duplicate family: the mates' 5' ends, sorted. Duplicates of one
    /// molecule share it (R3).
    pub fn family(&self) -> Vec<FivePrime> {
        let mut ends: Vec<FivePrime> = self.mates.iter().map(|m| m.five_prime.clone()).collect();
        ends.sort();
        ends
    }
}
```

- [ ] **Step 4: Run the tests and see them pass**

Run: `CARGO_TARGET_DIR=$S/target-t3 cargo test -j 16 origin:: -- --test-threads=16`
Expected: 13 passed.

- [ ] **Step 5: Commit**

```bash
git add src/origin.rs
git commit -m "code: origin -- records and fragments (surest mate, R4 removable, duplicate family)

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Look-alike regions

**Files:**
- Modify: `src/origin.rs`

**Interfaces:**
- Consumes: `OriginRecord`, `Span` (Tasks 2–3); the test helper `record` (Task 3).
- Produces: `pub fn lookalike_regions(spot: &[OriginRecord], footprint: &Span, read_length: u64) -> Vec<Span>`.

- [ ] **Step 1: Write the failing tests**

```rust
    #[test]
    fn test_hits_within_a_kilobase_form_one_region_grown_by_the_read_length() {
        let spot = [
            record("b", true, 200, 0, &[Span::new("chr1", 40_000, 40_100)]),
            record("b", false, 400, 0, &[Span::new("chr1", 40_900, 41_000)]),
        ];
        assert_eq!(
            lookalike_regions(&spot, &fp(), 150),
            vec![Span::new("chr1", 39_850, 41_150)]
        );
    }

    #[test]
    fn test_a_region_one_read_points_into_is_not_a_lookalike() {
        let spot = [
            record("b", true, 200, 0, &[Span::new("chr1", 40_000, 40_100)]),
            record("c", true, 300, 0, &[Span::new("chr1", 50_000, 50_100)]),
        ];
        assert!(lookalike_regions(&spot, &fp(), 150).is_empty());
    }

    #[test]
    fn test_hits_inside_the_footprint_are_not_lookalikes() {
        let spot = [
            record("b", true, 200, 0, &[Span::new("chr1", 600, 700)]),
            record("b", false, 400, 0, &[Span::new("chr1", 650, 750)]),
        ];
        assert!(lookalike_regions(&spot, &fp(), 150).is_empty());
    }

    #[test]
    fn test_hits_more_than_a_kilobase_apart_are_two_regions() {
        let spot = [
            record("b", true, 200, 0, &[Span::new("chr1", 40_000, 40_100), Span::new("chr1", 60_000, 60_100)]),
            record("c", true, 300, 0, &[Span::new("chr1", 40_050, 40_150), Span::new("chr1", 60_050, 60_150)]),
        ];
        assert_eq!(
            lookalike_regions(&spot, &fp(), 100),
            vec![Span::new("chr1", 39_900, 40_250), Span::new("chr1", 59_900, 60_250)]
        );
    }
```

- [ ] **Step 2: Run them and see them fail**

Run: `CARGO_TARGET_DIR=$S/target-t4 cargo test -j 16 origin:: -- --test-threads=16`
Expected: a compile error, `cannot find function lookalike_regions`.

- [ ] **Step 3: Implement** (add `use std::collections::BTreeSet;` at the top of `src/origin.rs`)

```rust
/// `XA` hits within this many bases of each other form one look-alike region.
const LOOKALIKE_GAP: u64 = 1000;
/// A region is a look-alike only when this many reads point into it.
const LOOKALIKE_MIN_READS: usize = 2;

/// The look-alike regions of `footprint`: the `XA` hits of the reads placed
/// over it that lie outside it, grouped when within [`LOOKALIKE_GAP`] of each
/// other and grown by `read_length` on each side. A region counts only when
/// at least [`LOOKALIKE_MIN_READS`] reads point into it.
pub fn lookalike_regions(spot: &[OriginRecord], footprint: &Span, read_length: u64) -> Vec<Span> {
    let mut hits: Vec<(&Span, (&str, bool))> = spot
        .iter()
        .flat_map(|r| {
            r.placements[1..]
                .iter()
                .map(move |p| (&p.span, (r.name.as_str(), r.first)))
        })
        .filter(|(span, _)| !footprint.overlaps(span))
        .collect();
    hits.sort();

    let mut regions = Vec::new();
    let mut i = 0;
    while i < hits.len() {
        let (first, _) = hits[i];
        let mut end = first.end;
        let mut reads: BTreeSet<(&str, bool)> = BTreeSet::new();
        while i < hits.len()
            && hits[i].0.chrom == first.chrom
            && hits[i].0.start <= end + LOOKALIKE_GAP
        {
            end = end.max(hits[i].0.end);
            reads.insert(hits[i].1);
            i += 1;
        }
        if reads.len() >= LOOKALIKE_MIN_READS {
            regions.push(Span::new(
                &first.chrom,
                first.start.saturating_sub(read_length),
                end + read_length,
            ));
        }
    }
    regions
}
```

- [ ] **Step 4: Run the tests and see them pass**

Run: `CARGO_TARGET_DIR=$S/target-t4 cargo test -j 16 origin:: -- --test-threads=16`
Expected: 17 passed.

- [ ] **Step 5: Commit**

```bash
git add src/origin.rs
git commit -m "code: origin -- look-alike regions from the spot's XA hits

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: The site, origin depth and `f`

**Files:**
- Modify: `src/origin.rs`

**Interfaces:**
- Consumes: Tasks 2–4; `crate::types::{ReadPair, ReadPool}`; `crate::stats::FragmentDist` in tests.
- Produces: `pub struct OriginSite { footprint: Span, lookalikes: Vec<Span>, records: Vec<OriginRecord>, f: f64 }` with `fragments(&self) -> Vec<Fragment<'_>>`, `removable_names(&self) -> Vec<String>`, `read_coverage_at(&self, chrom: &str, pos: u64, window: u64) -> f64` and `fragment_coverage_at(&self, chrom, pos, window) -> f64`; `pub fn fragment_to_read_ratio(pool: &ReadPool) -> anyhow::Result<f64>`.

- [ ] **Step 1: Write the failing tests**

```rust
    use crate::types::{ReadPair, ReadPool};

    /// Twin spots L = chr1:2000-2150 and P = chr1:12000-12150: 20 reads at
    /// each, all MAPQ 0 with their one XA hit on the other. Each read's mate
    /// is unmapped, so a fragment is one record, and each has its own 5' end
    /// so no two share a duplicate family.
    pub(super) fn twin_site() -> OriginSite {
        let (l, p) = (Span::new("chr1", 2000, 2150), Span::new("chr1", 12_000, 12_150));
        let at = |name: String, end5: u64, here: &Span, there: &Span| OriginRecord {
            placements: placements(here.clone(), 0, std::slice::from_ref(there)),
            mate_unmapped: true,
            ..record(&name, true, end5, 0, &[])
        };
        let records = (0..20u64)
            .map(|i| at(format!("l{}", i), i, &l, &p))
            .chain((0..20u64).map(|i| at(format!("p{}", i), 100 + i, &p, &l)))
            .collect();
        OriginSite { footprint: Span::new("chr1", 0, 5000), lookalikes: vec![Span::new("chr1", 11_850, 12_300)], records, f: 1.5 }
    }

    #[test]
    fn test_origin_depth_at_a_twin_is_the_true_read_depth() {
        // 20 reads at L count 1/2 each and so do their 20 twins at P.
        let site = twin_site();
        assert!(close(site.read_coverage_at("chr1", 2075, 100), 20.0));
        assert!(close(site.fragment_coverage_at("chr1", 2075, 100), 30.0));
    }

    #[test]
    fn test_duplicate_and_qc_fail_reads_add_no_depth() {
        // R3.
        let mut site = twin_site();
        let copy = site.records[0].clone();
        site.records.push(OriginRecord { name: "dup".into(), duplicate: true, ..copy.clone() });
        site.records.push(OriginRecord { name: "qc".into(), qc_fail: true, ..copy });
        assert!(close(site.read_coverage_at("chr1", 2075, 100), 20.0));
    }

    #[test]
    fn test_reads_spike_cannot_remove_add_no_depth() {
        // R6: a read at L whose mapped mate spike never read is kept (R4),
        // so its depth must not pay for new reads on top of it.
        let mut site = twin_site();
        let copy = site.records[0].clone();
        site.records.push(OriginRecord { name: "stray".into(), mate_unmapped: false, ..copy });
        assert!(!site.removable_names().contains(&"stray".to_string()));
        assert!(close(site.read_coverage_at("chr1", 2075, 100), 20.0));
    }

    #[test]
    fn test_fragments_are_the_names_with_a_placement_in_the_footprint() {
        let mut site = twin_site();
        site.records.push(record("unique_at_p", true, 12_500, 60, &[]));
        let names: Vec<&str> = site.fragments().iter().map(|f| f.name).collect();
        assert_eq!(names.len(), 40);
        assert!(!names.contains(&"unique_at_p"));
    }

    fn pair(name: &str, start: u64, span: u64, bases: usize) -> ReadPair {
        ReadPair {
            name: name.to_string(),
            seq1: vec![b'A'; bases],
            qual1: vec![b'I'; bases],
            seq2: vec![b'A'; bases],
            qual2: vec![b'I'; bases],
            ref_start: start,
            ref_end: start + span,
            insert_size: span as i64,
            chrom: "chr1".to_string(),
        }
    }

    pub(super) fn pool(pairs: Vec<ReadPair>) -> ReadPool {
        ReadPool { pairs, frag_dist: crate::stats::FragmentDist::from_stats(400.0, 80.0) }
    }

    #[test]
    fn test_f_is_fragment_span_over_read_bases_across_the_whole_pool() {
        // R1: 400 bp fragments of two 150 bp reads, far from any footprint.
        let p = pool((0..30).map(|i| pair(&format!("r{}", i), 500_000 + 10 * i, 400, 150)).collect());
        assert!(close(fragment_to_read_ratio(&p).unwrap(), 400.0 / 300.0));
    }

    #[test]
    fn test_f_refuses_a_pool_without_bases() {
        let p = pool(vec![pair("r", 0, 400, 0)]);
        assert!(fragment_to_read_ratio(&p).is_err());
    }
```

- [ ] **Step 2: Run them and see them fail**

Run: `CARGO_TARGET_DIR=$S/target-t5 cargo test -j 16 origin:: -- --test-threads=16`
Expected: a compile error, `cannot find struct OriginSite` and `cannot find function fragment_to_read_ratio`.

- [ ] **Step 3: Implement.** Add `use std::collections::BTreeMap;`, `use anyhow::{bail, Result};` and `use crate::types::ReadPool;` at the top, then:

```rust
/// Everything `origin` read for one event: the footprint the event's new
/// reads cover, the look-alike regions read besides it, every primary record
/// read in either, and the pool's fragment-to-read ratio.
#[derive(Debug, Clone)]
pub struct OriginSite {
    pub footprint: Span,
    pub lookalikes: Vec<Span>,
    pub records: Vec<OriginRecord>,
    /// See [`fragment_to_read_ratio`] (R1).
    pub f: f64,
}

impl OriginSite {
    /// The fragments with a placement inside the footprint, in name order.
    pub fn fragments(&self) -> Vec<Fragment<'_>> {
        let mut by_name: BTreeMap<&str, Vec<&OriginRecord>> = BTreeMap::new();
        for r in &self.records {
            by_name.entry(r.name.as_str()).or_default().push(r);
        }
        by_name
            .into_iter()
            .map(|(name, mates)| Fragment { name, mates })
            .filter(|f| {
                f.mates
                    .iter()
                    .any(|m| m.placements.iter().any(|p| self.footprint.holds(&p.span)))
            })
            .collect()
    }

    /// The names of the fragments spike may remove (R4).
    pub fn removable_names(&self) -> Vec<String> {
        self.fragments()
            .into_iter()
            .filter(|f| f.removable(&self.footprint))
            .map(|f| f.name.to_string())
            .collect()
    }

    /// Read depth that came from around `pos` on `chrom`. At up to 50 points
    /// over `window` (the points `simulate::estimate_coverage_at` samples),
    /// sum the chances of every placement covering the point, then average.
    /// Only fragments spike can remove count (R6): one kept by R4 stays in
    /// the BAM, so counting it would add new reads on top of it. Duplicate
    /// and QC-fail reads add nothing either, since spike's new reads are
    /// never flagged (R3).
    pub fn read_coverage_at(&self, chrom: &str, pos: u64, window: u64) -> f64 {
        let start = pos.saturating_sub(window / 2);
        let end = pos.saturating_add(window / 2);
        let range = end - start;
        let n = range.min(50).max(1);
        let step = if n > 1 { range / n } else { 1 };
        let removable: BTreeSet<String> = self.removable_names().into_iter().collect();
        let placed: Vec<&Placement> = self
            .records
            .iter()
            .filter(|r| !r.duplicate && !r.qc_fail && removable.contains(&r.name))
            .flat_map(|r| r.placements.iter())
            .filter(|p| p.span.chrom == chrom && p.span.start < end && start < p.span.end)
            .collect();
        let total: f64 = (0..n)
            .map(|i| {
                let at = start + i * step;
                placed
                    .iter()
                    .filter(|p| p.span.start <= at && at < p.span.end)
                    .map(|p| p.chance)
                    .sum::<f64>()
            })
            .sum();
        total / n as f64
    }

    /// [`read_coverage_at`](Self::read_coverage_at) in the fragment units the
    /// tiling count is in (R1).
    pub fn fragment_coverage_at(&self, chrom: &str, pos: u64, window: u64) -> f64 {
        self.read_coverage_at(chrom, pos, window) * self.f
    }
}

/// `f`: the pool's summed fragment spans over its summed read lengths (R1).
///
/// It depends on the library's fragment and read lengths, not on the spot,
/// so it is taken over the whole pool: an event inside a perfect twin has no
/// pool read in its footprint. Pool pairs keep no CIGAR, so read lengths are
/// whole reads. A soft-clipped read covers fewer bases than its length, so
/// where clips are common `f` runs a little low.
pub fn fragment_to_read_ratio(pool: &ReadPool) -> Result<f64> {
    let spans: u64 = pool.pairs.iter().map(|p| p.ref_end.saturating_sub(p.ref_start)).sum();
    let bases: u64 = pool.pairs.iter().map(|p| (p.seq1.len() + p.seq2.len()) as u64).sum();
    if bases == 0 {
        bail!(
            "the donor pool's {} read pair(s) hold no bases, so --edit-model origin cannot \
             convert read depth into fragment depth",
            pool.pairs.len()
        );
    }
    Ok(spans as f64 / bases as f64)
}
```

- [ ] **Step 4: Run the tests and see them pass**

Run: `CARGO_TARGET_DIR=$S/target-t5 cargo test -j 16 origin:: -- --test-threads=16`
Expected: 23 passed.

- [ ] **Step 5: Commit**

```bash
git add src/origin.rs
git commit -m "code: origin -- the site, origin depth (R3 no flagged depth, R6 only removable fragments) and f over the whole pool (R1)

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: Removal chances and one draw per family

**Files:**
- Modify: `src/origin.rs`

**Interfaces:**
- Consumes: Tasks 2–5; `crate::synth::copy_rate(copy: Option<bool>, vaf: f64) -> f64` (`src/synth.rs:843`).
- Produces: `pub struct Chance { name: String, family: Vec<FivePrime>, chance: f64 }`; `OriginSite::removal_chances(&self, read_copy: &HashMap<String, bool>, vaf: f64) -> Vec<Chance>`; `pub fn decide(chances: &[Chance], rng: &mut StdRng) -> BTreeSet<String>`.

- [ ] **Step 1: Write the failing tests**

```rust
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use std::collections::HashMap;

    fn chance(name: &str, family: u64, p: f64) -> Chance {
        Chance { name: name.to_string(), family: vec![FivePrime { chrom: "chr1".into(), pos: family, reverse: false }], chance: p }
    }

    #[test]
    fn test_removal_chance_is_p_origin_times_the_copy_rate() {
        // At the spot the phase call applies; at a look-alike the rate is vaf.
        let site = twin_site();
        let read_copy: HashMap<String, bool> = [("l0".to_string(), true)].into();
        let chances = site.removal_chances(&read_copy, 0.5);
        let of = |n: &str| chances.iter().find(|c| c.name == n).unwrap().chance;
        assert!(close(of("l0"), 0.5 * 1.0));
        assert!(close(of("l1"), 0.5 * 0.5));
        assert!(close(of("p0"), 0.5 * 0.5));
    }

    #[test]
    fn test_two_events_at_both_twins_remove_half_not_seven_sixteenths() {
        // R5: each event gives 1/2 x 1/2 = 1/4; together 1/2. Separate draws
        // would give 1 - (3/4)^2 = 7/16.
        let chances: Vec<Chance> = (0..20_000u64)
            .flat_map(|i| [chance(&format!("f{}", i), i, 0.25), chance(&format!("f{}", i), i, 0.25)])
            .collect();
        let removed = decide(&chances, &mut StdRng::seed_from_u64(1));
        let share = removed.len() as f64 / 20_000.0;
        assert!((0.48..0.52).contains(&share), "removed share {}", share);
    }

    #[test]
    fn test_a_total_above_one_always_removes() {
        let chances = [chance("a", 1, 0.7), chance("a", 1, 0.7)];
        for seed in 0..20 {
            assert!(decide(&chances, &mut StdRng::seed_from_u64(seed)).contains("a"));
        }
    }

    #[test]
    fn test_a_duplicate_shares_its_originals_fate() {
        // R3: 1000 families of two, each at 1/2.
        let chances: Vec<Chance> = (0..1000u64)
            .flat_map(|i| [chance(&format!("o{}", i), i, 0.5), chance(&format!("d{}", i), i, 0.5)])
            .collect();
        let removed = decide(&chances, &mut StdRng::seed_from_u64(2));
        let split = (0..1000).filter(|i| removed.contains(&format!("o{}", i)) != removed.contains(&format!("d{}", i))).count();
        assert_eq!(split, 0);
        assert!((400..600).contains(&(removed.len() / 2)), "{} removed", removed.len());
    }

    #[test]
    fn test_a_duplicate_takes_its_familys_phase_call() {
        // R7: phasing skips duplicates (src/loh.rs:628), so only the original
        // is in read_copy. At VAF 0.5 it gets rate 1, and so must its duplicate.
        let orig = OriginRecord { mate_unmapped: true, ..record("orig", true, 100, 60, &[]) };
        let dup = OriginRecord { name: "dup".into(), duplicate: true, ..orig.clone() };
        let site = OriginSite { footprint: fp(), lookalikes: vec![], records: vec![orig, dup], f: 1.0 };
        let on_event: HashMap<String, bool> = [("orig".to_string(), true)].into();
        let chances = site.removal_chances(&on_event, 0.5);
        assert_eq!(chances.len(), 2);
        assert!(close(chances[0].chance, chances[1].chance), "{:?}", chances);
        // On the other copy the original's rate is 0, and so is its duplicate's.
        let on_other: HashMap<String, bool> = [("orig".to_string(), false)].into();
        assert!(site.removal_chances(&on_other, 0.5).is_empty());
    }

    #[test]
    fn test_a_family_shares_one_fate_even_when_its_chances_differ() {
        // R7: one draw per family, against its highest member's total. Member
        // by member, a draw between 0.5 and 0.99 would remove only "orig".
        let chances = [chance("orig", 7, 0.99), chance("dup", 7, 0.5)];
        for seed in 0..200 {
            let removed = decide(&chances, &mut StdRng::seed_from_u64(seed));
            assert_eq!(removed.contains("orig"), removed.contains("dup"), "seed {}", seed);
        }
    }
```

- [ ] **Step 2: Run them and see them fail**

Run: `CARGO_TARGET_DIR=$S/target-t6 cargo test -j 16 origin:: -- --test-threads=16`
Expected: a compile error, `cannot find struct Chance`, `no method removal_chances` and `cannot find function decide`.

- [ ] **Step 3: Implement.** Add `use std::collections::HashMap;`, `use rand::rngs::StdRng;`, `use rand::Rng;` and `use crate::synth::copy_rate;` at the top, then:

```rust
/// One event's chance of removing one fragment.
#[derive(Debug, Clone, PartialEq)]
pub struct Chance {
    pub name: String,
    pub family: Vec<FivePrime>,
    pub chance: f64,
}

impl OriginSite {
    /// This event's chance of removing each fragment it can remove:
    /// `p_origin x copy_rate(copy, vaf)`. At the spot `copy` is the sample's
    /// phase call for the fragment's duplicate family (R7): phasing skips
    /// duplicates (`src/loh.rs:628`), so a family takes the call of whichever
    /// member has one, and members that disagree get none. A fragment the
    /// aligner put only at a look-alike has no call, so its rate is `vaf`.
    pub fn removal_chances(&self, read_copy: &HashMap<String, bool>, vaf: f64) -> Vec<Chance> {
        let fragments: Vec<Fragment<'_>> = self
            .fragments()
            .into_iter()
            .filter(|f| f.removable(&self.footprint))
            .collect();
        let mut family_copy: BTreeMap<Vec<FivePrime>, Option<bool>> = BTreeMap::new();
        for f in &fragments {
            if let Some(&copy) = read_copy.get(f.name) {
                family_copy
                    .entry(f.family())
                    .and_modify(|call| {
                        if *call != Some(copy) {
                            *call = None;
                        }
                    })
                    .or_insert(Some(copy));
            }
        }
        fragments
            .into_iter()
            .filter_map(|f| {
                let family = f.family();
                let copy = if f.at_spot(&self.footprint) {
                    family_copy.get(&family).copied().flatten()
                } else {
                    None
                };
                let chance = f.chance(&self.footprint) * copy_rate(copy, vaf);
                (chance > 0.0).then(|| Chance {
                    name: f.name.to_string(),
                    family,
                    chance,
                })
            })
            .collect()
    }
}

/// Which fragments to remove, over every event's chances at once (R5).
///
/// A fragment came from one copy, so "it came from event 1's edited copy"
/// and "it came from event 2's" cannot both be true: its chances add,
/// capped at 1. One molecule has one origin, so a duplicate family is
/// removed or kept whole (R3, R7). There is one draw per family, in family
/// order, against its highest member's total.
pub fn decide(chances: &[Chance], rng: &mut StdRng) -> BTreeSet<String> {
    let mut totals: BTreeMap<&str, (&[FivePrime], f64)> = BTreeMap::new();
    for c in chances {
        totals
            .entry(c.name.as_str())
            .or_insert((c.family.as_slice(), 0.0))
            .1 += c.chance;
    }
    let mut families: BTreeMap<&[FivePrime], f64> = BTreeMap::new();
    for (family, total) in totals.values() {
        let highest = families.entry(*family).or_insert(0.0);
        *highest = highest.max(total.min(1.0));
    }
    let removed_families: BTreeSet<&[FivePrime]> = families
        .into_iter()
        .filter(|(_, highest)| rng.gen::<f64>() < *highest)
        .map(|(family, _)| family)
        .collect();
    totals
        .into_iter()
        .filter(|(_, (family, _))| removed_families.contains(family))
        .map(|(name, _)| name.to_string())
        .collect()
}
```

- [ ] **Step 4: Run the tests and see them pass**

Run: `CARGO_TARGET_DIR=$S/target-t6 cargo test -j 16 origin:: -- --test-threads=16`
Expected: 29 passed.

- [ ] **Step 5: Commit**

```bash
git add src/origin.rs
git commit -m "code: origin -- removal chances; chances add over events, one draw and one phase call per family (R5, R3, R7)

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: Reading the site from a BAM or CRAM

**Files:**
- Modify: `src/extract.rs`. Add `write_one_contig_bam` to `pub(crate) mod test_fixtures` (`src/extract.rs:910-1073`).
- Modify: `src/origin.rs`

**Interfaces:**
- Consumes: Tasks 2–6. It also uses `crate::extract::{is_cram, build_fasta_repository, open_cram_reader_for_region, record_is_on_queried_reference, safe_noodles_position}`, the same calls `validate::scan_region` makes (`src/validate.rs:3279-3355`).
- Produces: `pub fn gather(bam_path: &str, ref_path: &str, footprint: &Span, read_length: usize, pool: &ReadPool) -> Result<OriginSite>`; `pub const XA_PROBE_RECORDS: usize = 100_000`; `pub fn first_xa_record(bam_path: &str, ref_path: &str, limit: usize) -> Result<Option<usize>>`; `pub fn require_xa(bam_path: &str, ref_path: &str) -> Result<usize>`; `pub(crate) fn test_fixtures::write_one_contig_bam(path: &Path, contig: &str, contig_len: usize, records: &[RecordBuf]) -> String`.

- [ ] **Step 1: Add the fixture** (to `extract::test_fixtures`; this is test code, so no red step)

```rust
    /// Write `records`, sorted by position, as a BAM on one contig. Add a
    /// hand-written `.bai` of one bin that holds them all, the way
    /// `census::tests::write_pairs_bam` does: a query anywhere on the
    /// contig reads the whole file. Returns the BAM's path.
    pub(crate) fn write_one_contig_bam(
        path: &std::path::Path,
        contig: &str,
        contig_len: usize,
        records: &[noodles::sam::alignment::RecordBuf],
    ) -> String {
        use noodles::sam::alignment::io::Write as _;
        let header = noodles::sam::Header::builder()
            .add_reference_sequence(
                contig,
                noodles::sam::header::record::value::Map::<
                    noodles::sam::header::record::value::map::ReferenceSequence,
                >::new(std::num::NonZeroUsize::try_from(contig_len).unwrap()),
            )
            .build();
        {
            let mut writer = noodles::bam::io::writer::Builder
                .build_from_path(path)
                .unwrap();
            writer.write_header(&header).unwrap();
            for r in records {
                writer.write_alignment_record(&header, r).unwrap();
            }
            writer.try_finish().unwrap();
        }
        let (first_record, end_of_file) = {
            let mut reader = noodles::bam::io::Reader::new(std::fs::File::open(path).unwrap());
            reader.read_header().unwrap();
            let first_record = reader.get_ref().virtual_position();
            let mut record = noodles::bam::Record::default();
            while reader.read_record(&mut record).unwrap() != 0 {}
            (first_record, reader.get_ref().virtual_position())
        };
        let mut bai: Vec<u8> = Vec::new();
        bai.extend_from_slice(b"BAI\x01");
        bai.extend_from_slice(&1u32.to_le_bytes()); // one reference
        bai.extend_from_slice(&1u32.to_le_bytes()); // one bin
        bai.extend_from_slice(&0u32.to_le_bytes()); // bin 0 spans the contig
        bai.extend_from_slice(&1u32.to_le_bytes()); // one chunk
        bai.extend_from_slice(&u64::from(first_record).to_le_bytes());
        bai.extend_from_slice(&u64::from(end_of_file).to_le_bytes());
        bai.extend_from_slice(&0u32.to_le_bytes()); // no linear index
        std::fs::write(format!("{}.bai", path.display()), bai).unwrap();
        path.to_str().unwrap().to_string()
    }
```

- [ ] **Step 2: Write the failing tests** (in `origin::tests`)

```rust
    use noodles::sam::alignment::record::{Flags, MappingQuality};
    use noodles::sam::alignment::RecordBuf;

    fn bam_record(name: &str, flags: u16, pos: usize, mapq: u8, xa: Option<&str>, mate_pos: usize) -> RecordBuf {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::{Data, QualityScores, Sequence};
        let data: Data = xa.map(|xa| (Tag::new(b'X', b'A'), Value::from(xa))).into_iter().collect();
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::from(flags))
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles::core::Position::new(pos).unwrap())
            .set_mapping_quality(MappingQuality::new(mapq).unwrap())
            .set_cigar([Op::new(Kind::Match, 100)].into_iter().collect())
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(noodles::core::Position::new(mate_pos).unwrap())
            .set_sequence(Sequence::from(vec![b'A'; 100]))
            .set_quality_scores(QualityScores::from(vec![30u8; 100]))
            .set_data(data)
            .build()
    }

    fn scratch(label: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("spike_origin_{}_{}", label, std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn bases_pool() -> ReadPool {
        pool((0..30).map(|i| pair(&format!("r{}", i), 500_000 + 10 * i, 400, 150)).collect())
    }

    /// chrT, 60 kb; footprint chrT:20000-25000, twin P at 42000.
    /// a: unique at the spot. b: MAPQ 0 at the spot, hits at P. c: MAPQ 0 at
    /// P, hits at the spot. d: unique at P. e: read 1 at the spot, read 2
    /// past the footprint's end.
    fn twin_bam(dir: &std::path::Path) -> String {
        let (r1, r2) = (0x63u16, 0x93u16);
        let records = [
            bam_record("a", r1, 21_001, 60, None, 21_201),
            bam_record("a", r2, 21_201, 60, None, 21_001),
            bam_record("b", r1, 22_001, 0, Some("chrT,+42001,100M,0;"), 22_201),
            bam_record("b", r2, 22_201, 0, Some("chrT,-42201,100M,0;"), 22_001),
            bam_record("e", r1, 24_801, 60, None, 25_101),
            bam_record("e", r2, 25_101, 60, None, 24_801),
            bam_record("c", r1, 42_001, 0, Some("chrT,+22001,100M,0;"), 42_201),
            bam_record("c", r2, 42_201, 0, Some("chrT,-22201,100M,0;"), 42_001),
            bam_record("d", r1, 42_301, 60, None, 42_501),
            bam_record("d", r2, 42_501, 60, None, 42_301),
        ];
        crate::extract::test_fixtures::write_one_contig_bam(&dir.join("twin.bam"), "chrT", 60_000, &records)
    }

    #[test]
    fn test_gather_reads_the_spot_and_its_lookalike() {
        let dir = scratch("gather");
        let bam = twin_bam(&dir);
        let site = gather(&bam, "", &Span::new("chrT", 20_000, 25_000), 100, &bases_pool()).unwrap();

        assert_eq!(site.lookalikes, vec![Span::new("chrT", 41_900, 42_400)]);
        let names: Vec<&str> = site.fragments().iter().map(|f| f.name).collect();
        assert_eq!(names, ["a", "b", "c", "e"]);
        assert_eq!(site.removable_names(), ["a", "b", "c"]);
        let c = site.fragments().into_iter().find(|f| f.name == "c").unwrap();
        assert!(!c.at_spot(&site.footprint));
        assert!(close(c.chance(&site.footprint), 0.5));
        assert!(close(site.f, 400.0 / 300.0));
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Two MAPQ 0 reads at the spot, with no XA anywhere in the file.
    fn noxa_bam(dir: &std::path::Path) -> String {
        let records = [
            bam_record("z", 0x63, 21_001, 0, None, 21_201),
            bam_record("z", 0x93, 21_201, 0, None, 21_001),
        ];
        crate::extract::test_fixtures::write_one_contig_bam(&dir.join("noxa.bam"), "chrT", 60_000, &records)
    }

    #[test]
    fn test_gather_does_not_stop_at_a_spot_whose_mapq0_reads_lack_xa() {
        // R8: under bwa-mem's -h 5 rule such reads have more than 5 hits, so
        // each gets 1/6. Whether the file kept XA at all is `require_xa`'s question.
        let dir = scratch("spot_noxa");
        let bam = noxa_bam(&dir);
        let site = gather(&bam, "", &Span::new("chrT", 20_000, 25_000), 100, &bases_pool()).unwrap();
        let z = site.fragments().into_iter().find(|f| f.name == "z").unwrap();
        assert!(close(z.chance(&site.footprint), 1.0 / 6.0));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_first_xa_record_counts_records_up_to_the_first_xa() {
        // In `twin_bam` the first record with XA is b's read 1, the third.
        let dir = scratch("first_xa");
        let bam = twin_bam(&dir);
        assert_eq!(first_xa_record(&bam, "", 100).unwrap(), Some(3));
        assert_eq!(first_xa_record(&bam, "", 2).unwrap(), None);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_require_xa_stops_a_file_without_xa() {
        // R8.
        let dir = scratch("require_xa");
        let err = require_xa(&noxa_bam(&dir), "").unwrap_err().to_string();
        assert!(err.contains("XA") && err.contains("--edit-model clean"), "{}", err);
        assert_eq!(require_xa(&twin_bam(&dir), "").unwrap(), 3);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_first_xa_record_reads_a_cram() {
        // The two-contig CRAM carries SA:Z on every record and no XA.
        let dir = scratch("cram_xa");
        let (fasta, cram) = crate::extract::test_fixtures::write_two_contig_cram(&dir);
        assert_eq!(first_xa_record(&cram, &fasta, 100).unwrap(), None);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_gather_reads_only_the_queried_contig_of_a_cram() {
        // A container holding two contigs is decoded whole (L2, N4).
        let dir = scratch("cram");
        let (fasta, cram) = crate::extract::test_fixtures::write_two_contig_cram(&dir);
        let site = gather(&cram, &fasta, &Span::new("chrA", 0, 1000), 100, &bases_pool()).unwrap();
        let names: Vec<&str> = site.fragments().iter().map(|f| f.name).collect();
        assert_eq!(names, ["chrA_pair0", "chrA_pair1", "chrA_pair2"]);
        let _ = std::fs::remove_dir_all(&dir);
    }
```

- [ ] **Step 3: Run them and see them fail**

Run: `CARGO_TARGET_DIR=$S/target-t7 cargo test -j 16 origin:: -- --test-threads=16`
Expected: a compile error, `cannot find function gather` (and `first_xa_record`, `require_xa`).

- [ ] **Step 4: Implement.** Change the top-level anyhow import to `use anyhow::{bail, Context, Result};`, then:

```rust
/// One [`OriginRecord`] from an alignment record, or `None` for a record
/// `origin` does not judge on its own. That is an unmapped, secondary or
/// supplementary record (those share their primary's name and go with it),
/// or one without a name. MAPQ 255, "unavailable", counts as 0.
fn origin_record(
    header: &noodles::sam::Header,
    buf: &noodles::sam::alignment::RecordBuf,
) -> Option<OriginRecord> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    let flags = buf.flags();
    if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
        return None;
    }
    let name = buf.name()?.to_string();
    let (chrom, _) = header
        .reference_sequences()
        .get_index(buf.reference_sequence_id()?)?;
    let chrom = chrom.to_string();
    let start = usize::from(buf.alignment_start()?) as u64 - 1;
    let ops: &[Op] = buf.cigar().as_ref();
    let xa = match buf.data().get(&Tag::new(b'X', b'A')) {
        Some(Value::String(s)) => Some(s.to_string()),
        _ => None,
    };
    let alternatives = xa.as_deref().map(parse_xa).unwrap_or_default();
    let mapq = buf.mapping_quality().map(|q| q.get()).unwrap_or(0);
    Some(OriginRecord {
        name,
        first: flags.is_first_segment(),
        placements: placements(
            Span::new(&chrom, start, start + reference_length(ops)),
            mapq,
            &alternatives,
        ),
        duplicate: flags.is_duplicate(),
        qc_fail: flags.is_qc_fail(),
        mate_unmapped: flags.is_mate_unmapped(),
        five_prime: five_prime(&chrom, start, ops, flags.is_reverse_complemented()),
    })
}

/// Every primary record overlapping `span`, from a BAM or a CRAM.
fn scan(bam_path: &str, ref_path: &str, span: &Span) -> Result<Vec<OriginRecord>> {
    let region = noodles::core::Region::new(
        span.chrom.as_str(),
        crate::extract::safe_noodles_position(span.start + 1)
            ..=crate::extract::safe_noodles_position(span.end),
    );
    let mut out = Vec::new();
    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let (mut reader, header) =
            crate::extract::open_cram_reader_for_region(bam_path, &repository, &region)
                .context("failed to open CRAM for the origin scan")?;
        let queried = header
            .reference_sequences()
            .get_index_of(span.chrom.as_bytes());
        for result in reader.query(&header, &region)? {
            let buf = result?.try_into_alignment_record(&header)?;
            // A container holding several contigs is decoded whole (L2, N4).
            if !crate::extract::record_is_on_queried_reference(&buf, queried) {
                continue;
            }
            out.extend(origin_record(&header, &buf));
        }
    } else {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .with_context(|| format!("failed to open BAM for the origin scan: {}", bam_path))?;
        let header = reader.read_header()?;
        for result in reader.query(&header, &region)? {
            let record = result?;
            let buf =
                noodles::sam::alignment::RecordBuf::try_from_alignment_record(&header, &record)?;
            out.extend(origin_record(&header, &buf));
        }
    }
    Ok(out)
}

/// How many records, from the start of the file, [`require_xa`] reads
/// looking for an `XA` tag. The 35x HG002 BAM's first `XA` is on record 64.
pub const XA_PROBE_RECORDS: usize = 100_000;

/// The number of records read up to and including the first that carries
/// `XA`, or `None` when none of the first `limit` does.
pub fn first_xa_record(bam_path: &str, ref_path: &str, limit: usize) -> Result<Option<usize>> {
    use noodles::sam::alignment::record::data::field::Tag;
    let xa = Tag::new(b'X', b'A');
    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let mut reader = noodles::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(bam_path)
            .with_context(|| format!("failed to open CRAM: {}", bam_path))?;
        let header = reader.read_header()?;
        for (i, result) in reader.records(&header).take(limit).enumerate() {
            let buf = result?.try_into_alignment_record(&header)?;
            if buf.data().get(&xa).is_some() {
                return Ok(Some(i + 1));
            }
        }
    } else {
        let mut reader = noodles::bam::io::reader::Builder
            .build_from_path(bam_path)
            .with_context(|| format!("failed to open BAM: {}", bam_path))?;
        reader.read_header()?;
        for (i, result) in reader.records().take(limit).enumerate() {
            if result?.data().get(&xa).is_some() {
                return Ok(Some(i + 1));
            }
        }
    }
    Ok(None)
}

/// Stop unless the file keeps the aligner's `XA` tags (R8); otherwise every
/// MAPQ 0 read would silently get 1/6 and no look-alike would be found.
///
/// The check is on the file, not the spot: a spot whose MAPQ 0 reads carry
/// no `XA` is normal, since under bwa-mem's `-h 5` rule their hits number
/// more than 5. Returns the record the first `XA` is on.
pub fn require_xa(bam_path: &str, ref_path: &str) -> Result<usize> {
    match first_xa_record(bam_path, ref_path, XA_PROBE_RECORDS)? {
        Some(n) => Ok(n),
        None => bail!(
            "--edit-model origin needs the aligner's XA tags (its alternative hits), but none \
             of the first {} records of {} carries one. bwa-mem and bwa-mem2 write XA by \
             default, and a later step can strip it. Re-align with one of them, or use \
             --edit-model clean.",
            XA_PROBE_RECORDS,
            bam_path,
        ),
    }
}

/// Read everything `origin` needs for one event: the footprint, its
/// look-alike regions and `f` (see [`OriginSite`]). Whether the file keeps
/// `XA` at all is [`require_xa`]'s check, made once per run.
pub fn gather(
    bam_path: &str,
    ref_path: &str,
    footprint: &Span,
    read_length: usize,
    pool: &ReadPool,
) -> Result<OriginSite> {
    let spot = scan(bam_path, ref_path, footprint)?;
    let lookalikes = lookalike_regions(&spot, footprint, read_length as u64);
    let mut seen: BTreeSet<(String, bool)> =
        spot.iter().map(|r| (r.name.clone(), r.first)).collect();
    let mut records = spot;
    for region in &lookalikes {
        for r in scan(bam_path, ref_path, region)? {
            if seen.insert((r.name.clone(), r.first)) {
                records.push(r);
            }
        }
    }
    Ok(OriginSite {
        footprint: footprint.clone(),
        lookalikes,
        records,
        f: fragment_to_read_ratio(pool)?,
    })
}
```

- [ ] **Step 5: Run the tests and see them pass**

Run: `CARGO_TARGET_DIR=$S/target-t7 cargo test -j 16 origin:: -- --test-threads=16`
Expected: 35 passed.

- [ ] **Step 6: Commit**

```bash
git add src/origin.rs src/extract.rs
git commit -m "code: origin -- gather a site from BAM or CRAM; the XA check reads the file, not the spot (R8)

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 8: The origin path in `simulate`

**Files:**
- Modify: `src/types.rs`. Add a `SplicedOutput` field (`src/types.rs:232-256`).
- Modify: `src/simulate.rs`: the event body at `src/simulate.rs:114-283`, `depth_fold` at `:925-952`, and the `spliced` test helper at `:3056-3066`. Add new functions beside `donor_coverage_for_tiling` (`:453`).

**Interfaces:**
- Consumes: `OriginSite`, `Span`, `Chance` (Tasks 5–6).
- Produces:
  - `SplicedOutput::origin_chances: Vec<crate::origin::Chance>`;
  - `pub fn is_additive(event: &SimEvent, dup_model: &str) -> bool`;
  - `pub fn simulate_event_origin(event_index, event, pool, haplotype, config, synth_gen, vaf, origin: Option<&OriginSite>, rng) -> Result<SplicedOutput>`;
  - `pub fn apply_removals(outputs: &mut [SplicedOutput], removed: &BTreeSet<String>)`;
  - `fn depth_fold_by(haplotype, cov, coverage_at: &dyn Fn(&str, u64, u64) -> f64) -> DepthFold`;
  - `fn origin_coverage_for_tiling(site, haplotype, pool, fallback_bp) -> Result<(f64, Vec<String>)>`.

- [ ] **Step 1: Write the failing tests** (in `simulate::tests`)

```rust
    use crate::origin::{self, OriginRecord, OriginSite, Span};

    /// A MAPQ 60 read of 150 bp at `start` on chr1, with no XA.
    fn origin_read(name: &str, first: bool, start: u64) -> OriginRecord {
        OriginRecord {
            name: name.to_string(),
            first,
            placements: origin::placements(Span::new("chr1", start, start + 150), 60, &[]),
            duplicate: false,
            qc_fail: false,
            mate_unmapped: false,
            five_prime: origin::FivePrime { chrom: "chr1".into(), pos: start, reverse: !first },
        }
    }

    /// 80 unique pairs over chr1:0-5000, the footprint of `del_haplotype(2000, 1000)`.
    fn origin_site(footprint: Span) -> OriginSite {
        let records = (0..80u64)
            .flat_map(|i| {
                let s = 100 + 50 * i;
                [origin_read(&format!("o{}", i), true, s), origin_read(&format!("o{}", i), false, s + 250)]
            })
            .collect();
        OriginSite { footprint, lookalikes: vec![], records, f: 1.0 }
    }

    fn origin_del_event() -> SimEvent {
        SimEvent::Deletion {
            chrom: "chr1".to_string(),
            del_start: 2000,
            del_end: 3000,
            gene: "TEST".to_string(),
            exons: vec![],
            allele_fraction: Some(0.5),
        }
    }

    #[test]
    fn test_origin_keeps_an_event_whose_pool_has_no_read_in_the_footprint() {
        // R1: the pool is 500 kb away, so `clean` refuses this event
        // (test_simulate_event_refuses_pool_with_no_coverage_at_breakpoint).
        // Under origin the depth comes from the site.
        let mut hap = del_haplotype(2000, 1000);
        let pool = make_covering_pool(500_000, 540_000, 200);
        let site = origin_site(Span::new("chr1", 0, 5000));
        let mut rng = StdRng::seed_from_u64(42);
        let out = simulate_event_origin(
            1, &origin_del_event(), &pool, &mut hap, &make_config(), &mock_synth_gen(150), 0.5, Some(&site), &mut rng,
        )
        .unwrap();

        // Nothing is suppressed here: `origin::decide` draws after all events.
        assert_eq!(out.kept_originals.len(), 200);
        assert_eq!(out.suppressed_count, 0);
        assert!(!out.chimeric_pairs.is_empty());
        assert_eq!(out.origin_chances.len(), 80);
        assert!(out.origin_chances.iter().all(|c| (c.chance - 0.5 * (1.0 - 1e-6)).abs() < 1e-9));
    }

    #[test]
    fn test_origin_leaves_the_pool_pairs_inside_the_footprint_for_the_final_draw() {
        // These are the pairs `clean` suppresses here, at rate 1/2. Under
        // origin every one comes back kept: `origin::decide` removes them only
        // after every event's chances are summed (R5).
        let mut hap = del_haplotype(2000, 1000);
        let pool = make_covering_pool(0, 4000, 100);
        assert!(pool.pairs.iter().all(|p| p.ref_end <= 5000), "every pair lies inside chr1:0-5000");
        let site = origin_site(Span::new("chr1", 0, 5000));
        let mut rng = StdRng::seed_from_u64(42);
        let out = simulate_event_origin(
            1, &origin_del_event(), &pool, &mut hap, &make_config(), &mock_synth_gen(150), 0.5, Some(&site), &mut rng,
        )
        .unwrap();
        assert_eq!(out.suppressed_count, 0);
        assert_eq!(out.kept_originals.len(), 100);
    }

    #[test]
    fn test_origin_site_that_misses_the_haplotype_footprint_is_a_bug() {
        let mut hap = del_haplotype(2000, 1000);
        let pool = make_covering_pool(0, 5000, 200);
        let site = origin_site(Span::new("chr1", 0, 4000));
        let mut rng = StdRng::seed_from_u64(42);
        let err = simulate_event_origin(
            1, &origin_del_event(), &pool, &mut hap, &make_config(), &mock_synth_gen(150), 0.5, Some(&site), &mut rng,
        )
        .err()
        .expect("a mismatched site must be refused")
        .to_string();
        assert!(err.contains("bug"), "{}", err);
    }

    #[test]
    fn test_depth_fold_by_measures_bins_with_the_given_estimator() {
        let hap = tandem_dup_haplotype(11_000, 15_000, 1_000);
        let fold = depth_fold_by(&hap, 40.0, &|_, _, _| 10.0);
        assert!((fold.fold - 41.0 / 11.0).abs() < 1e-9, "{:?}", fold);
    }

    #[test]
    fn test_apply_removals_moves_removed_pool_pairs_to_suppressed() {
        let mut outputs = vec![spliced(&["a", "b"], &[], &["c"])];
        let removed: BTreeSet<String> = ["a".to_string(), "not_in_a_pool".to_string()].into();
        apply_removals(&mut outputs, &removed);
        assert_eq!(sorted_names(&outputs[0].kept_originals), ["b"]);
        assert_eq!(outputs[0].suppressed_names, ["c", "a"]);
        assert_eq!(outputs[0].suppressed_count, 2);
        // A removed read in no pool is neither re-emitted nor listed here;
        // main lists it in replaced_reads.txt itself.
        assert!(!consumed_original_names(&outputs).contains("not_in_a_pool"));
    }

    #[test]
    fn test_is_additive_names_fusions_and_junction_dups() {
        assert!(!is_additive(&origin_del_event(), "full"));
        let dup = SimEvent::Duplication { chrom: "chr1".into(), dup_start: 10, dup_end: 20, gene: "G".into(), allele_fraction: None };
        assert!(is_additive(&dup, "junction"));
        assert!(!is_additive(&dup, "full"));
    }
```

- [ ] **Step 2: Run them and see them fail**

Run: `CARGO_TARGET_DIR=$S/target-t8 cargo test -j 16 simulate:: -- --test-threads=16`
Expected: a compile error, `cannot find function simulate_event_origin` (and `depth_fold_by`, `apply_removals`, `is_additive`, `no field origin_chances`).

- [ ] **Step 3: Implement**

(a) `src/types.rs`: add to `SplicedOutput`, after `depth_fold`:

```rust
    /// Under `--edit-model origin`, this event's chance of removing each
    /// fragment it can remove. `origin::decide` sums them over every event
    /// and draws once per fragment. Empty under `clean`.
    pub origin_chances: Vec<crate::origin::Chance>,
```

In `simulate.rs`, add `origin_chances: Vec::new(),` to the `spliced` test helper.

(b) `src/simulate.rs`: add `is_additive` and use it in the event body. Replace the `let is_additive = match event { ... };` block (`src/simulate.rs:146-150`) with `let is_additive = is_additive(event, &config.dup_model);`, keeping the comment above it. Add:

```rust
/// Whether `event` only adds reads: a fusion, or a DUP under the legacy
/// junction model. It removes no original, so `--edit-model` does not apply.
pub fn is_additive(event: &SimEvent, dup_model: &str) -> bool {
    match event {
        SimEvent::Fusion { .. } => true,
        SimEvent::Duplication { .. } => dup_model == "junction",
        _ => false,
    }
}
```

(c) Entry points. `simulate_event` becomes a wrapper:

```rust
pub fn simulate_event(
    event_index: usize,
    event: &SimEvent,
    pool: &ReadPool,
    haplotype: &mut VariantHaplotype,
    config: &SimConfig,
    synth_gen: &SynthReadGenerator,
    vaf: f64,
    rng: &mut StdRng,
) -> Result<SplicedOutput> {
    simulate_event_origin(event_index, event, pool, haplotype, config, synth_gen, vaf, None, rng)
}

/// [`simulate_event`] under `--edit-model origin` when `origin` is `Some`.
/// The site's fragments are judged by where they came from, and the pool's
/// pairs are not suppressed here. `None` is exactly `simulate_event`.
#[allow(clippy::too_many_arguments)]
pub fn simulate_event_origin(
    event_index: usize,
    event: &SimEvent,
    pool: &ReadPool,
    haplotype: &mut VariantHaplotype,
    config: &SimConfig,
    synth_gen: &SynthReadGenerator,
    vaf: f64,
    origin: Option<&crate::origin::OriginSite>,
    rng: &mut StdRng,
) -> Result<SplicedOutput> {
    let copies = sample_copies_for_event(event, haplotype, config, synth_gen.reference(), rng)?;
    simulate_event_inner(
        event_index, event, pool, haplotype, config, synth_gen, vaf, &copies, origin, rng,
    )
}
```

Rename the body of `simulate_event_with_copies` to `fn simulate_event_inner(..., copies: &[(String, loh::SampleCopies)], origin: Option<&crate::origin::OriginSite>, rng: &mut StdRng)`, with the same `#[allow(clippy::too_many_arguments)]`. Keep `simulate_event_with_copies` as a wrapper that passes `None`:

```rust
#[allow(clippy::too_many_arguments)]
fn simulate_event_with_copies(
    event_index: usize,
    event: &SimEvent,
    pool: &ReadPool,
    haplotype: &mut VariantHaplotype,
    config: &SimConfig,
    synth_gen: &SynthReadGenerator,
    vaf: f64,
    copies: &[(String, loh::SampleCopies)],
    rng: &mut StdRng,
) -> Result<SplicedOutput> {
    simulate_event_inner(
        event_index, event, pool, haplotype, config, synth_gen, vaf, copies, None, rng,
    )
}
```

(d) In `simulate_event_inner`, right after `let (hap_ref_start, hap_ref_end) = ...;`:

```rust
    if let Some(site) = origin {
        let footprint =
            crate::origin::Span::new(haplotype.primary_chrom(), hap_ref_start, hap_ref_end);
        anyhow::ensure!(
            !is_additive && site.footprint == footprint,
            "origin site {} does not match the haplotype's footprint {} -- this is a bug; \
             please report it",
            site.footprint,
            footprint,
        );
    }
```

In the suppression loop, change the `replaceable` line to:

```rust
        // Under origin nothing is suppressed here: every event's chances are
        // summed first and drawn once, in `origin::decide` (R5).
        let replaceable = origin.is_none()
            && !is_additive
            && classify_pair_relation(pair, hap_ref_start, hap_ref_end) == PairRelation::Inside;
```

`rng.gen` stays behind `replaceable &&`, so `clean` draws exactly as before.

Replace the `donor_coverage_for_tiling` call and the `depth_fold` line with:

```rust
    let (cov, uncovered_breakpoint_sides) = match origin {
        Some(site) => {
            origin_coverage_for_tiling(site, haplotype, pool, (&first_bp_chrom, first_bp_ref))?
        }
        None => donor_coverage_for_tiling(
            event,
            haplotype,
            pool,
            (&sv_chrom, sv_start, sv_end),
            (&first_bp_chrom, first_bp_ref),
        )?,
    };

    // CR2: measured only. Under origin each bin is measured with origin
    // depth, the estimator `cov` came from (the T3 rule).
    let depth_fold = match origin {
        Some(site) => depth_fold_by(haplotype, cov, &|chrom, pos, window| {
            site.fragment_coverage_at(chrom, pos, window)
        }),
        None => depth_fold(haplotype, pool, cov),
    };
```

In the returned `SplicedOutput`, add:

```rust
        origin_chances: origin
            .map(|site| site.removal_chances(&read_copy, vaf))
            .unwrap_or_default(),
```

(e) `depth_fold_by`. Move `depth_fold`'s body into it, with the one depth line changed. `depth_fold` then calls it:

```rust
/// [`depth_fold`] with each bin's depth measured by `coverage_at(chrom,
/// center, width)` instead of the pool.
fn depth_fold_by(
    haplotype: &VariantHaplotype,
    cov: f64,
    coverage_at: &dyn Fn(&str, u64, u64) -> f64,
) -> DepthFold {
    const BIN: u64 = 1000;
    let mut worst = DepthFold {
        fold: 1.0,
        scaled_by: cov,
        worst_bin: String::new(),
        worst_depth: cov,
    };
    for origin in haplotype.segments.iter().filter_map(|seg| seg.origin.as_ref()) {
        let len = origin.ref_end.saturating_sub(origin.ref_start);
        if len == 0 {
            continue;
        }
        let n_bins = ((len as f64 / BIN as f64).round() as u64).max(1);
        for b in 0..n_bins {
            let start = origin.ref_start + len * b / n_bins;
            let end = origin.ref_start + len * (b + 1) / n_bins;
            let depth = coverage_at(&origin.chrom, (start + end) / 2, end - start);
            let fold = ((depth + 1.0) / (cov + 1.0)).max((cov + 1.0) / (depth + 1.0));
            if fold > worst.fold {
                worst.fold = fold;
                worst.worst_bin = format!("{}:{}-{}", origin.chrom, start, end);
                worst.worst_depth = depth;
            }
        }
    }
    worst
}

fn depth_fold(haplotype: &VariantHaplotype, pool: &ReadPool, cov: f64) -> DepthFold {
    depth_fold_by(haplotype, cov, &|chrom, pos, window| {
        estimate_coverage_at(pool, chrom, pos, window)
    })
}
```

Keep `depth_fold`'s existing doc comment on `depth_fold`.

(f) `origin_coverage_for_tiling` (next to `donor_coverage_for_tiling`):

```rust
/// [`donor_coverage_for_tiling`] under `--edit-model origin` (R1).
///
/// The depth is the site's origin depth in fragment units, not the pool's,
/// which is 0 inside a perfect twin. The event is refused only when no
/// breakpoint side has origin depth above 0; the depth returned is the first
/// covered side's. Only events that remove reads get a site, and those come
/// from one locus, so the fusion rule does not arise.
fn origin_coverage_for_tiling(
    site: &crate::origin::OriginSite,
    haplotype: &VariantHaplotype,
    pool: &ReadPool,
    fallback_bp: (&str, u64),
) -> Result<(f64, Vec<String>)> {
    let mut sides = breakpoint_sides(haplotype);
    if sides.is_empty() {
        sides.push((fallback_bp.0.to_string(), fallback_bp.1));
    }
    let covs: Vec<f64> = sides
        .iter()
        .map(|(chrom, pos)| site.fragment_coverage_at(chrom, *pos, 2000))
        .collect();
    let covered = |cov: f64| !cov.is_nan() && cov > 0.0;
    let uncovered: Vec<String> = sides
        .iter()
        .zip(&covs)
        .filter(|(_, &c)| !covered(c))
        .map(|((chrom, pos), _)| format!("{}:{}", chrom, pos))
        .collect();
    let Some(((chrom, pos), cov)) = sides
        .iter()
        .zip(covs.iter().copied())
        .find(|(_, c)| covered(*c))
    else {
        anyhow::bail!(
            "event over {} has no origin depth at any of its breakpoints ({}): no read at \
             the spot or at its {} look-alike region(s) could have come from there, so \
             spike would invent the reads it plants and still write a truth VCF beside them.",
            site.footprint,
            uncovered.join(", "),
            site.lookalikes.len(),
        );
    };
    log::info!(
        "  origin depth at {}:{}: {:.1}x (the donor pool's there: {:.1}x)",
        chrom,
        pos,
        cov,
        estimate_coverage_at(pool, chrom, *pos, 2000),
    );
    Ok((cov, uncovered))
}
```

(g) `apply_removals` (after `combine_event_outputs`):

```rust
/// Under `--edit-model origin`, move each event's pool pairs that
/// `origin::decide` removed from `kept_originals` to `suppressed_names`, so
/// [`combine_event_outputs`] drops them and [`consumed_original_names`] still
/// lists them. Removed names in no pool are the caller's to list.
pub fn apply_removals(outputs: &mut [SplicedOutput], removed: &BTreeSet<String>) {
    for output in outputs.iter_mut() {
        let (gone, kept): (Vec<ReadPair>, Vec<ReadPair>) =
            std::mem::take(&mut output.kept_originals)
                .into_iter()
                .partition(|p| removed.contains(&p.name));
        output.kept_originals = kept;
        output.suppressed_names.extend(gone.into_iter().map(|p| p.name));
        output.suppressed_count = output.suppressed_names.len();
    }
}
```

- [ ] **Step 4: Run the tests and see them pass**

Run: `CARGO_TARGET_DIR=$S/target-t8 cargo test -j 16 -- --test-threads=16`
Expected: every test passes, including the 6 new ones and every existing `simulate::` test. If `test_merge_script_aborts_when_original_bam_does_not_match_replaced_reads` fails, it is the known flaky test: re-run it alone and report both outputs.

- [ ] **Step 5: Commit**

```bash
git add src/types.rs src/simulate.rs
git commit -m "code: origin -- simulate's origin path (no suppression here, origin depth and refusal R1, depth fold on origin depth)

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 9: Wiring it into `main`

**Files:**
- Modify: `src/main.rs`: the event loop (`src/main.rs:537-640`), the block after the refusal bail (`:642-660`), and a new `origin_footprint` next to `build_haplotype` (`:1439`). Add tests in `mod tests`.

**Interfaces:**
- Consumes: `origin::{require_xa, gather, decide, Span, OriginSite, Chance}`, `simulate::{is_additive, simulate_event_origin, apply_removals}`, `HAP_FLANK` (`src/main.rs:201`), `SharedReference::chromosome_length(&str) -> Option<u64>`.
- Produces: `fn origin_footprint(event: &SimEvent, contig_len: u64) -> Option<origin::Span>`.

- [ ] **Step 1: Write the failing test** (in `main.rs` `mod tests`)

```rust
    #[test]
    fn test_origin_footprint_is_the_haplotypes_reference_range() {
        let del = |start, end| SimEvent::Deletion {
            chrom: "chr1".into(), del_start: start, del_end: end, gene: "G".into(), exons: vec![], allele_fraction: None,
        };
        let span = |s, e| Some(origin::Span::new("chr1", s, e));
        assert_eq!(origin_footprint(&del(5000, 6000), 100_000), span(3000, 8000));
        // The flanks stop at the contig's ends, as `VariantHaplotype::from_deletion` does.
        assert_eq!(origin_footprint(&del(99_000, 99_500), 100_000), span(97_000, 100_000));
        assert_eq!(origin_footprint(&del(500, 600), 100_000), span(0, 2600));
    }
```

- [ ] **Step 2: Run it and see it fail**

Run: `CARGO_TARGET_DIR=$S/target-t9 cargo test -j 16 origin_footprint -- --test-threads=16`
Expected: a compile error, `cannot find function origin_footprint`.

- [ ] **Step 3: Implement**

(a) The footprint helper, next to `build_haplotype`:

```rust
/// The reference range an event's haplotype covers: its region grown by
/// `HAP_FLANK` on each side, stopped at the contig's ends, as
/// `VariantHaplotype::ref_range` reports it once the haplotype is built.
/// `simulate` checks the two agree. `None` for a fusion, which has no single
/// region.
fn origin_footprint(event: &SimEvent, contig_len: u64) -> Option<origin::Span> {
    let (chrom, start, end) = event.primary_region()?;
    Some(origin::Span::new(
        chrom,
        start.saturating_sub(HAP_FLANK),
        (end + HAP_FLANK).min(contig_len),
    ))
}
```

(b) Before the event loop:

```rust
    let edit_origin = args.edit_model == "origin";
    if edit_origin {
        // R8: once per run, on the file. A spot whose MAPQ 0 reads lack XA is
        // normal (more than 5 hits); a file with no XA at all was stripped.
        let first = origin::require_xa(&args.bam, &args.reference)?;
        log::info!(
            "--edit-model origin (experimental): reads are removed by their chance of having \
             come from each event's edited copy, at the event and at its look-alikes. The BAM \
             keeps XA tags (the first is on record {}).",
            first,
        );
    }
    // Under origin: each event's label and site, for the removal log after
    // the one draw over every event.
    let mut origin_sites: Vec<(String, origin::OriginSite)> = Vec::new();
```

(c) In the loop, right after `extract_pool_for_event(...)`:

```rust
        // --edit-model origin: every primary read at the event and at its
        // look-alikes. Additive events remove nothing, so they get no site.
        let site = if edit_origin && !simulate::is_additive(event, &config.dup_model) {
            let chrom = event
                .primary_region()
                .map(|(chrom, _, _)| chrom.to_string())
                .expect("an event that removes reads has one region");
            let contig_len = shared_ref
                .chromosome_length(&chrom)
                .ok_or_else(|| anyhow::anyhow!("{} is not in the loaded reference", chrom))?;
            let footprint = origin_footprint(event, contig_len)
                .expect("an event that removes reads has one region");
            let site = origin::gather(
                &config.bam_path,
                &config.ref_path,
                &footprint,
                config.read_length,
                &pool,
            )?;
            log::info!(
                "  origin: {} fragment(s) could have come from {}; look-alike region(s): {}",
                site.fragments().len(),
                footprint,
                if site.lookalikes.is_empty() {
                    "none".to_string()
                } else {
                    site.lookalikes.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(", ")
                },
            );
            Some(site)
        } else {
            None
        };
```

Extend `editable` so the census counts origin's removable reads as editable:

```rust
        let editable: std::collections::HashSet<String> = pool
            .pairs
            .iter()
            .map(|p| p.name.clone())
            .chain(unusable_qual_names.iter().cloned())
            .chain(site.iter().flat_map(|s| s.removable_names()))
            .collect();
```

Replace the `simulate::simulate_event(...)` call with:

```rust
        let output = simulate::simulate_event_origin(
            i + 1,
            event,
            &pool,
            &mut haplotype,
            &config,
            &synth_gen,
            vaf,
            site.as_ref(),
            &mut rng,
        )?;
```

At the end of the loop body, after `event_outputs.push(output);`:

```rust
        if let Some(site) = site {
            origin_sites.push((event_label(event), site));
        }
```

(d) After the `refusal_message` bail and before `consumed_original_names`:

```rust
    // --edit-model origin: one draw per fragment family over every event's
    // chances (R5). Removed pool pairs move to their event's suppressed
    // names; every removed name goes to replaced_reads.txt.
    let mut origin_removed: BTreeSet<String> = BTreeSet::new();
    if edit_origin {
        let chances: Vec<origin::Chance> = event_outputs
            .iter()
            .flat_map(|o| o.origin_chances.iter().cloned())
            .collect();
        origin_removed = origin::decide(&chances, &mut rng);
        simulate::apply_removals(&mut event_outputs, &origin_removed);
        for (stat, output) in event_stats.iter_mut().zip(&event_outputs) {
            stat.suppressed = output.suppressed_count;
        }
        for (label, site) in &origin_sites {
            let (mut at_spot, mut elsewhere) = (0usize, 0usize);
            for f in site.fragments() {
                if origin_removed.contains(f.name) {
                    if f.at_spot(&site.footprint) {
                        at_spot += 1;
                    } else {
                        elsewhere += 1;
                    }
                }
            }
            log::info!(
                "{}: origin removed {} fragment(s) at the spot and {} at its look-alikes",
                label,
                at_spot,
                elsewhere,
            );
        }
    }
```

Then change `let mut replaced_names = simulate::consumed_original_names(&event_outputs);` to be followed by:

```rust
    replaced_names.extend(origin_removed.iter().cloned());
```

- [ ] **Step 4: Run the whole suite, clippy and fmt**

Run: `CARGO_TARGET_DIR=$S/target-t9 cargo test -j 16 -- --test-threads=16`
Expected: 0 failed. Record the passed count; it should be 554 + the new tests, counted from the output and not predicted.
Run: `CARGO_TARGET_DIR=$S/target-t9 cargo test --release -j 16 -- --test-threads=16`
Expected: 0 failed.
Run: `CARGO_TARGET_DIR=$S/target-t9 cargo clippy --all-targets -j 16 2>&1 | grep -E "generated [0-9]+ warnings"`
Expected: `13 warnings` (bin) and `14 warnings (12 duplicates)` (test), as on `94f3f20`.
Run: `cargo fmt --check` and report its output (fix only lines this plan added).

- [ ] **Step 5: Check that `clean` is byte-identical to master, on real data**

```bash
S=/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad
REF=data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
BAM=data/giab_hg38/HG002/HG002.novaseq.pcr-free.35x.bwamem2.dedup.grch38_no_alt.bam
mkdir -p $S/origin-c1 && cd /home/parlar_ai/dev/spike
git worktree add --detach $S/origin-c1/base 94f3f20
(cd $S/origin-c1/base && CARGO_TARGET_DIR=$S/target-base cargo build --release -j 16)
CARGO_TARGET_DIR=$S/target-t9 cargo build --release -j 16
md5sum $S/target-base/release/spike $S/target-t9/release/spike
samtools view -b -o $S/origin-c1/slice.bam $BAM chr20:14500000-14600000 && samtools index $S/origin-c1/slice.bam
E="--bam $S/origin-c1/slice.bam --reference $REF --event del:chr20:14548422-14548735 --seed 1 --allow-resistant"
$S/target-base/release/spike $E -o $S/origin-c1/base_run 2> $S/origin-c1/base.log
$S/target-t9/release/spike $E -o $S/origin-c1/new_default 2> $S/origin-c1/new_default.log
$S/target-t9/release/spike $E --edit-model clean -o $S/origin-c1/new_clean 2> $S/origin-c1/new_clean.log
for d in base_run new_default new_clean; do
  echo "$d $(zcat $S/origin-c1/$d/R1.fq.gz | md5sum | cut -c1-8) $(zcat $S/origin-c1/$d/R2.fq.gz | md5sum | cut -c1-8) $(md5sum < $S/origin-c1/$d/replaced_reads.txt | cut -c1-8) $(grep -v '^##' $S/origin-c1/$d/truth.vcf | md5sum | cut -c1-8)"
done
git worktree remove $S/origin-c1/base
```

Expected: the three rows hold the same four md5s. A row that differs fails Global Constraint 1; stop and report it.

- [ ] **Step 6: Run origin on real data, with the XA tags kept and with them stripped**

```bash
# Kept: the slice holds 328 records with XA ("Measured before planning").
$S/target-t9/release/spike $E --edit-model origin -o $S/origin-c1/new_origin 2> $S/origin-c1/new_origin.log; echo "exit $?"
grep -E "keeps XA tags|origin:|origin depth|origin removed" $S/origin-c1/new_origin.log
# Stripped: the same slice without any XA.
samtools view -b -x XA -o $S/origin-c1/slice_noxa.bam $S/origin-c1/slice.bam && samtools index $S/origin-c1/slice_noxa.bam
samtools view $S/origin-c1/slice_noxa.bam | grep -c "XA:Z:"
$S/target-t9/release/spike --bam $S/origin-c1/slice_noxa.bam --reference $REF --event del:chr20:14548422-14548735 --seed 1 --allow-resistant --edit-model origin -o $S/origin-c1/noxa_origin 2> $S/origin-c1/noxa_origin.log; echo "exit $?"
grep -c "needs the aligner's XA tags" $S/origin-c1/noxa_origin.log
```

Expected (PREDICTED, not run):
- The kept run exits 0 and logs a `keeps XA tags` line, an `origin:` line, an `origin depth` line and an `origin removed` line.
- The stripped slice counts 0 `XA:Z:`.
- The stripped run exits non-zero, with at least 1 matching error line.

Record the actual lines.

- [ ] **Step 7: Commit**

```bash
git add src/main.rs
git commit -m "code: origin -- main gathers a site per event and draws once over every event

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 10: User docs

**Files:**
- Modify: `README.md`. Add a subsection after the `--min-mapq 0` bullet list under "Reads spike cannot edit" (after `README.md:475`), and add the `--edit-model` entry to the CLI reference block (after `--dup-model`, `README.md:1390-1393`).

- [ ] **Step 1: Add the subsection**

```markdown
#### Editing hard spots: `--edit-model origin` (experimental)

`--min-mapq 0` makes the MAPQ 0 reads editable, but it removes each one as if it certainly came from the event. It also never touches the look-alike copy the aligner split those reads with. `--edit-model origin` models what a real variant does there instead:

- Every primary read over the event's footprint (the event ± 2 kb), and over its look-alike regions, gets a chance of having come from the footprint. The chance comes from its MAPQ and its `XA` alternative hits.
- It is removed by that chance times its copy's rate, but only if spike's new reads can replace it: every mapped mate could have come from inside the footprint.
- Chances from several events add up. A duplicate shares its original's fate.
- The number of new reads comes from where reads came from (the origin depth). It does not come from the donor pool, which can be empty inside a perfect twin.

It needs the aligner's `XA` tags. bwa-mem and bwa-mem2 write them by default, and spike stops if none of the first 100,000 records of the BAM carries one. A spot whose MAPQ 0 reads lack `XA` is normal: under bwa-mem's `-h 5` rule their hits number more than 5, so each gets a chance of 1/6. At `chr20:7117236-7121236` in the 35x HG002 BAM, 822 reads are MAPQ 0 and 1 carries `XA`, yet the file's first 100,000 records hold 15,255 with it. The default stays `clean` until `origin` is tested against real data. The design is in `docs/superpowers/specs/2026-09-26-edit-model-origin-design.md`.
```

- [ ] **Step 2: Paste the `--edit-model` entry.** Build the release binary, run `$S/target-t9/release/spike --help`, and copy its `--edit-model <EDIT_MODEL>` entry verbatim into the CLI reference block after `--dup-model`.

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "docs: origin -- README section and CLI reference for --edit-model

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 11: Break it on purpose (mutation check)

Apply each mutation alone and run `CARGO_TARGET_DIR=$S/target-mut cargo test -j 16 -- --test-threads=16`. Record which tests went red, then revert with `git checkout src/`. **Pass:** every mutation reddens at least one test through a FAIL, not a compile error. Record the actual red tests; do not copy the "should redden" column.

| # | Mutation | Should redden (PREDICTED, not run) |
| --- | --- | --- |
| M1 | `MAX_LISTED_HITS` 5 → 4 | `test_p_here_follows_mapq_and_the_xa_count` |
| M2 | `chance_within` counts only `placements[..1]` | the 3/4 and the 1 tests |
| M3 | `surest` uses `min_by` | `test_a_pair_takes_its_surest_mates_chance` |
| M4 | `removable` drops `every_mate_read` | `test_an_unseen_mate_blocks_removal_unless_it_is_unmapped` |
| M5 | `removable` uses `.any(` for the mates | `test_a_pair_is_removable_only_when_every_mate_could_come_from_the_footprint` |
| M6 | `read_coverage_at` drops the duplicate/QC-fail filter | `test_duplicate_and_qc_fail_reads_add_no_depth` |
| M7 | `decide` draws once per chance and removes on any hit (a union) | `test_two_events_at_both_twins_remove_half_not_seven_sixteenths` |
| M8 | `decide` keys draws by name, not family | `test_a_duplicate_shares_its_originals_fate` |
| M9 | `LOOKALIKE_MIN_READS` 2 → 1 | `test_a_region_one_read_points_into_is_not_a_lookalike` |
| M10 | `require_xa` returns `Ok(0)` without reading the file | `test_require_xa_stops_a_file_without_xa` |
| M11 | `scan` (CRAM) drops `record_is_on_queried_reference` | `test_gather_reads_only_the_queried_contig_of_a_cram` |
| M12 | `replaceable` drops `origin.is_none() &&` | `test_origin_leaves_the_pool_pairs_inside_the_footprint_for_the_final_draw` |
| M13 | the origin arm calls `donor_coverage_for_tiling` | `test_origin_keeps_an_event_whose_pool_has_no_read_in_the_footprint` |
| M14 | `apply_removals` leaves `kept_originals` alone | `test_apply_removals_moves_removed_pool_pairs_to_suppressed` |
| M15 | `removal_chances` uses `copy = None` always | `test_removal_chance_is_p_origin_times_the_copy_rate` |
| M16 | `fragment_to_read_ratio` divides by the pair count | `test_f_is_fragment_span_over_read_bases_across_the_whole_pool` |
| M17 | `read_coverage_at` drops `removable.contains(&r.name)` (R6) | `test_reads_spike_cannot_remove_add_no_depth` |
| M18 | `removal_chances` takes `read_copy.get(f.name)` instead of the family's call (R7) | `test_a_duplicate_takes_its_familys_phase_call` |
| M19 | `decide` removes a name when its family's draw is below the name's own total (R7) | `test_a_family_shares_one_fate_even_when_its_chances_differ` |
| M20 | `first_xa_record`'s CRAM branch returns `Ok(Some(1))` | `test_first_xa_record_reads_a_cram` |
| M21 | `gather` stops when the spot's MAPQ 0 reads carry no `XA` (the per-spot check R8 replaced) | `test_gather_does_not_stop_at_a_spot_whose_mapq0_reads_lack_xa` |

- [ ] **Step 1: Run M1–M21 and record the red tests per mutation** in `$S/origin-mut.txt`.
- [ ] **Step 2: Confirm `git status --short` shows only ` M .gitignore` and `?? .ignore`.** Nothing is committed in this task; its record goes into the result commit (Task 12).

---

### Task 12: The physics test (the judgment gate)

**Files:**
- Create: `scripts/origin_physics.py`
- Modify: `docs/review/REVIEW.md`. Add a new section at the end, "`--edit-model origin`: physics test".

**Locked in this plan, before any run.** Commit 1 of the gate is the commit of this file.

- **Genome.** `chrT`, 80 kb, from `random.Random(20260927)`: U1 `[0,20000)`, S `[20000,30000)`, U2 `[30000,50000)`, an exact copy of S at `[50000,60000)`, U3 `[60000,80000)`.
- **Event.** `del:chrT:24500-25500`, het (`--allele-fraction` default 0.5). Its 5 kb footprint lies inside S.
- **Reads.** wgsim `-1 150 -2 150 -d 400 -s 50 -e 0.001 -r 0 -R 0 -X 0`, at 30x per copy (`N = 30 × total_haplotype_length / 300` pairs). Qualities are rewritten to Q30 and `/1` `/2` suffixes stripped. Aligned with `bwa-mem2 mem -R '@RG\tID:toy\tSM:TOY'`, then sorted and indexed.
- **Samples.** Truth seeds 1 and 2: copy A with `[24500,25500)` removed, plus copy B. Donor seed 3: both copies intact. spike `--seed 7 --allow-resistant`, three ways: `clean`, `--min-mapq 0` and `--edit-model origin`. Each goes through `align.sh` and `merge.sh` with 8 threads.
- **Measured.** Mean depth from `samtools depth -a` (default flags; any MAPQ) over L = `chrT:24600-25400` and P = `chrT:54600-55400`. Each is divided by the length-weighted mean over `chrT:5000-15000`, `35000-45000` and `65000-75000`.
- **Band.** For X in {L, P}: `[min(truth1_X, truth2_X) − 0.05, max(truth1_X, truth2_X) + 0.05]`.
- **Verdict**, checked in this order. A control counts as evidence only when it ran as meant:
  - **NO VERDICT** in three cases: `origin` exits non-zero; `--min-mapq 0` exits non-zero; or `clean` exits non-zero without its expected refusal, whose log line holds `has no donor coverage`. Nothing was measured. Fix the cause, commit it as `code:`, run again, and report that it happened.
  - **INCONCLUSIVE** if `--min-mapq 0` lies inside both bands, or `clean` exits 0 and lies inside both. The test then cannot tell the models apart. Stop and report.
  - **SUPPORTED** if `origin` lies inside both bands.
  - **REFUTED** otherwise. The default stays `clean` and no tuning follows. Report and stop.
  - So the gate can reject only a `--min-mapq 0` run that finished and landed outside a band, and a `clean` run that finished outside a band or stopped with its expected refusal.
- **Harness check (not a verdict).** In the donor BAM, at `chrT:22500-27500`, count primary records, MAPQ 0 records, and records with `XA`. If the MAPQ 0 records carry no `XA`, the harness is broken: fix it and re-run. This is not a verdict.
- **Prediction on paper** (spec table, PREDICTED, not run): truth ≈ 0.75 at L and P; `--min-mapq 0` ≈ 0.5 at L and ≈ 1 at P; `origin` ≈ 0.75 at both; `clean` refused.

- [ ] **Step 1: Write the harness** `scripts/origin_physics.py`

```python
#!/usr/bin/env python3
"""--edit-model origin, the physics test: a made-up genome with an exact twin.

chrT (80 kb): unique U1 [0,20000), S [20000,30000), unique U2 [30000,50000), an
exact copy of S [50000,60000), unique U3 [60000,80000). The event is a het 1 kb
deletion in the middle of the first S, so its 5 kb footprint lies inside S.

Truth: reads from the diploid genome carrying the deletion on one copy (two
seeds). Donor: reads from it without. spike runs on the donor three ways --
clean, --min-mapq 0, --edit-model origin -- each through align.sh and merge.sh.
Measured: depth over L and P (any MAPQ) over the unique windows' mean. The rule
is locked in docs/superpowers/plans/2026-09-27-edit-model-origin.md, Task 12.

Usage: origin_physics.py OUT_DIR SPIKE [THREADS]
"""
import os
import random
import subprocess
import sys

GENOME_SEED = 20260927
TRUTH_SEEDS = (1, 2)
DONOR_SEED = 3
SPIKE_SEED = 7
PER_COPY = 30
READ_LEN = 150
MARGIN = 0.05
DELETED = (24500, 25500)
L = (24600, 25400)
P = (54600, 55400)
BASE = [(5000, 15000), (35000, 45000), (65000, 75000)]
EVENT = "del:chrT:24500-25500"


def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def genome():
    rng = random.Random(GENOME_SEED)
    seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    u1, s, u2, u3 = seq(20000), seq(10000), seq(20000), seq(20000)
    return u1 + s + u2 + s + u3


def reads(out, name, copies, seed, ref, threads):
    """wgsim from `copies`, qualities to Q30, aligned; returns the BAM."""
    haps = f"{out}/{name}.haps.fa"
    write_fasta(haps, [(f"copy{i}", c) for i, c in enumerate(copies)])
    pairs = PER_COPY * sum(len(c) for c in copies) // (2 * READ_LEN)
    raw = [f"{out}/{name}.raw{i}.fq" for i in (1, 2)]
    sh(["wgsim", "-N", str(pairs), "-1", str(READ_LEN), "-2", str(READ_LEN), "-d", "400",
        "-s", "50", "-e", "0.001", "-r", "0", "-R", "0", "-X", "0", "-S", str(seed),
        haps, *raw], stdout=subprocess.DEVNULL)
    fq = [f"{out}/{name}.R{i}.fq" for i in (1, 2)]
    for src, dst in zip(raw, fq):
        with open(src) as fin, open(dst, "w") as fout:
            for i, line in enumerate(fin):
                line = line.rstrip("\n")
                if i % 4 == 0 and line[-2:] in ("/1", "/2"):
                    line = line[:-2]
                if i % 4 == 3:
                    line = "?" * len(line)
                fout.write(line + "\n")
    bam = f"{out}/{name}.bam"
    align = subprocess.Popen(["bwa-mem2", "mem", "-t", str(threads), "-R", r"@RG\tID:toy\tSM:TOY",
                              ref, *fq], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    sh(["samtools", "sort", "-@", str(threads), "-o", bam, "-"], stdin=align.stdout)
    if align.wait() != 0:
        raise SystemExit(f"bwa-mem2 failed for {name}")
    sh(["samtools", "index", bam])
    return bam


def spike(out, name, spike_bin, donor, ref, extra, threads):
    """spike -> align.sh -> merge.sh; returns (spike exit, merged BAM or None)."""
    d = f"{out}/{name}"
    with open(f"{d}.log", "w") as log:
        rc = subprocess.run([spike_bin, "--bam", donor, "--reference", ref, "--event", EVENT,
                             "--seed", str(SPIKE_SEED), "--allow-resistant", "-o", d, *extra],
                            stderr=log, stdout=log).returncode
        if rc != 0:
            return rc, None
        sh(["bash", f"{d}/align.sh", ref, str(threads)], stdout=log, stderr=log)
        sh(["bash", f"{d}/merge.sh", donor, ref, str(threads)], stdout=log, stderr=log)
    return 0, f"{d}/merged.bam"


def mean_depth(bam, start, end):
    out = subprocess.run(["samtools", "depth", "-a", "-r", f"chrT:{start + 1}-{end}", bam],
                         capture_output=True, text=True, check=True).stdout
    vals = [int(line.split("\t")[2]) for line in out.splitlines()]
    return sum(vals) / len(vals)


def ratios(bam):
    base = sum(mean_depth(bam, a, b) * (b - a) for a, b in BASE) / sum(b - a for a, b in BASE)
    return mean_depth(bam, *L) / base, mean_depth(bam, *P) / base


def harness_check(donor):
    out = subprocess.run(["samtools", "view", "-F", "0x904", donor, "chrT:22501-27500"],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    mapq0 = [r for r in out if r.split("\t")[4] == "0"]
    with_xa = [r for r in mapq0 if "\tXA:Z:" in r]
    print(f"harness: {len(out)} primary records over the footprint, {len(mapq0)} MAPQ 0, "
          f"{len(with_xa)} of those with XA")
    if mapq0 and not with_xa:
        raise SystemExit("harness broken: the donor's MAPQ 0 reads carry no XA")


def main(out, spike_bin, threads="8"):
    os.makedirs(out, exist_ok=True)
    g = genome()
    ref = f"{out}/ref.fa"
    write_fasta(ref, [("chrT", g)])
    sh(["samtools", "faidx", ref])
    sh(["bwa-mem2", "index", ref], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    deleted = g[:DELETED[0]] + g[DELETED[1]:]

    truths = [ratios(reads(out, f"truth{s}", [deleted, g], s, ref, threads)) for s in TRUTH_SEEDS]
    donor = reads(out, "donor", [g, g], DONOR_SEED, ref, threads)
    harness_check(donor)

    runs = {}
    for name, extra in [("clean", []), ("mapq0", ["--min-mapq", "0"]),
                        ("origin", ["--edit-model", "origin"])]:
        rc, merged = spike(out, name, spike_bin, donor, ref, extra, threads)
        runs[name] = (rc, ratios(merged) if merged else None)

    band = [(min(t[i] for t in truths) - MARGIN, max(t[i] for t in truths) + MARGIN) for i in (0, 1)]
    inside = lambda r: r is not None and all(band[i][0] <= r[i] <= band[i][1] for i in (0, 1))
    for s, t in zip(TRUTH_SEEDS, truths):
        print(f"truth{s}\texit -\tL {t[0]:.3f}\tP {t[1]:.3f}")
    for name, (rc, r) in runs.items():
        cells = f"L {r[0]:.3f}\tP {r[1]:.3f}" if r else "refused"
        print(f"{name}\texit {rc}\t{cells}\tinside {inside(r)}")
    print(f"band L [{band[0][0]:.3f}, {band[0][1]:.3f}]  P [{band[1][0]:.3f}, {band[1][1]:.3f}]")

    mapq0_rc, mapq0 = runs["mapq0"]
    clean_rc, clean = runs["clean"]
    origin_rc, origin = runs["origin"]
    with open(f"{out}/clean.log") as fh:
        clean_refused_as_expected = clean_rc != 0 and "has no donor coverage" in fh.read()
    # A control is evidence only when it ran as meant.
    if origin_rc != 0:
        print("NO VERDICT: origin failed to run -- fix it and run again")
    elif mapq0_rc != 0:
        print("NO VERDICT: the --min-mapq 0 control failed to run -- fix it and run again")
    elif clean_rc != 0 and not clean_refused_as_expected:
        print("NO VERDICT: the clean control failed, but not with its expected refusal")
    elif inside(mapq0) or (clean_rc == 0 and inside(clean)):
        print("VERDICT: INCONCLUSIVE -- the test cannot tell the models apart")
    elif inside(origin):
        print("VERDICT: SUPPORTED")
    else:
        print("VERDICT: REFUTED")


if __name__ == "__main__":
    main(*sys.argv[1:])
```

- [ ] **Step 2: Commit the harness before running it**

```bash
git add scripts/origin_physics.py
git commit -m "code: origin -- physics test harness (made-up twin genome; rule locked in the plan)

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 3: Run it once**

```bash
CARGO_TARGET_DIR=$S/target-t12 cargo build --release -j 16 && md5sum $S/target-t12/release/spike
python3 scripts/origin_physics.py $S/origin-physics $S/target-t12/release/spike 8 2>&1 | tee $S/origin-physics.txt
grep -h "origin" $S/origin-physics/origin.log | head -20
```

Record the full output. Do not re-run with other seeds, windows or margins. Only a harness-check failure or a NO VERDICT allows a second run. Fix the cause, commit the fix as `code:`, run again, and report that it happened and why.

- [ ] **Step 4: Write the result into `docs/review/REVIEW.md`.** It holds the locked rule (pointing to this plan), the harness check line, the table of truth1, truth2, clean, mapq0 and origin (exit, L, P), the bands and the verdict. Add the Task 9 byte-identity md5 rows and the Task 11 mutation record. For each measured number, give the command that produced it. Update the spec's Status line to say whether the physics test passed. Update the README section's last paragraph with the verdict in one sentence.

- [ ] **Step 5: Commit the result** (never together with the plan commit)

```bash
git add docs/review/REVIEW.md README.md docs/superpowers/specs/2026-09-26-edit-model-origin-design.md
git commit -m "result: --edit-model origin physics test -- <SUPPORTED / REFUTED / INCONCLUSIVE>

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 6: If REFUTED or INCONCLUSIVE**, add a case to `.claude/judgment-gate-cases.md`: what was claimed, what it rested on, and which gate caught it. Commit it separately as `docs:`.

---

## Not in this plan

- **The real test on the user's GIAB BAMs (arriving 2026-09-28).** It gets its own gate plan once the BAMs are in hand. First count their `XA` tags and read their `@PG` lines (duplicate marking); the spec's assumptions table asks for both. Then pick sites, and lock the plan before any comparison is run.
- Making `origin` the default. Only the real test can do that.
- `spike validate`'s coverage expectations at look-alike spots (the spec's known side effect).
- Merging `edit-model` into master. That waits for the user's word.
