//! BAM validation subcommand for spike.
//!
//! Reads a simulated BAM + truth VCF and runs automated checks to verify
//! that the spike-in reads look realistic.

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};

use anyhow::{bail, Context, Result};
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::Cigar as CigarTrait;

use crate::census;
use crate::extract::record_is_on_queried_reference;

/// Parsed command-line arguments for `spike validate`.
struct ValidateArgs {
    bam_path: String,
    truth_path: String,
    ref_path: String,
    min_mapq: u8,
    flank_bp: u64,
    json_output: bool,
    /// Whether the advisory checks count towards the exit status.
    strict: bool,
}

/// The event label the three whole-sample checks are reported under.
const GLOBAL_LABEL: &str = "[global]";

/// A truth VCF event with its expected properties.
#[derive(Debug)]
struct TruthEvent {
    chrom: String,
    start: u64,      // 0-based
    end: u64,        // 0-based, exclusive
    sv_type: String, // DEL, DUP, INV, INS, BND, SNP
    expected_vaf: f64,
    gene: String,
    /// The other breakpoint (chrom, 0-based position): a split read at this
    /// event should have its other part there. None for INS and SNP.
    partner: Option<(String, u64)>,
    /// For SmallVariant: explicit REF/ALT alleles.
    ref_allele: Option<Vec<u8>>,
    alt_allele: Option<Vec<u8>>,
    /// For INS: the inserted length (`SVLEN`). An insertion has no reference
    /// span, so `end == start` and the length cannot be read off the
    /// coordinates the way every other type's can. None for every other type.
    ins_len: Option<u64>,
    /// For INS: the record's ALT column exactly as the truth VCF wrote it --
    /// the anchor base and then the inserted sequence, or a symbolic `<INS>`
    /// from a spike old enough not to have recorded the bases (CR7). None for
    /// every other type.
    ///
    /// Deliberately **not** `alt_allele`. [`NearbyRecords::new`] selects the
    /// records the small-variant path treats as edits by
    /// `ref_allele.is_some() && alt_allele.is_some()`, and
    /// [`NearbyRecords::inside`] applies each to the reference at its own
    /// `start`; an INS's `start` is the VCF POS under the SV convention where a
    /// small variant's is POS - 1, so an insertion put there would enter the
    /// haplotype enumeration `allele_freq` votes against, one base off, and
    /// move a non-advisory row. This field is read by [`ins_alt_sequence`],
    /// which is the only reader of it anywhere.
    ins_alt: Option<Vec<u8>>,
    /// `N` when the record's ID is spike's `sim_ins_N` (an INS) or `sim_del_N`
    /// (a DEL): the number spike's own reads for the event carry as
    /// `evNNNN_hap_` (RF13, RF14). None for any other ID, and for every other
    /// type.
    sim_number: Option<u32>,
    /// The census numbers spike recorded for this event, as INFO holds them.
    census: CensusInfo,
}

/// The census numbers spike wrote into one truth record's INFO, as text.
///
/// `SIM_RESIST` is the share of the reads over the event spike could not edit
/// (CR4); `SIM_DEPTH_FOLD` is how far the donor's depth departs from the one
/// depth the event was scaled by (CR2). The text is kept unparsed so that a
/// field that is not there -- an older spike's truth VCF -- can be told from
/// one that is there and unreadable: the first is no row at all, the second an
/// advisory failure.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
struct CensusInfo {
    resist: Option<String>,
    depth_fold: Option<String>,
}

impl CensusInfo {
    /// The two fields read off one record's INFO column. A field that is
    /// absent and VCF's own `.` are both `None`: neither is a number spike
    /// measured, and neither is a malformed one.
    fn from_info(info: &str) -> Self {
        let field = |key: &str| match parse_info_field(info, key) {
            None | Some(".") => None,
            Some(value) => Some(value.to_string()),
        };
        CensusInfo {
            resist: field("SIM_RESIST"),
            depth_fold: field("SIM_DEPTH_FOLD"),
        }
    }
}

/// Result of a single validation check.
struct CheckResult {
    event_label: String,
    check_name: String,
    expected: String,
    observed: String,
    pass: bool,
    /// An advisory row reports something spike measured rather than something
    /// this run verified. It is printed and counted with the rest, but it
    /// stays out of the exit status unless `--strict` is given, and it never
    /// stands in for a check of the event (M11).
    advisory: bool,
}

/// A check's result; a check that could not run is a failed result, so an
/// unevaluable run never reports all-PASS.
///
/// `advisory` is the standing of the check being wrapped, and it governs the
/// row either way -- the result the check produced and the "check runs"
/// failure it leaves when it could not run. An advisory check that errored
/// must stay advisory: a non-advisory error row would enter the exit status
/// (a run that exits 0 today would start exiting 1) and would satisfy
/// `check_event`'s "a check applies" fallback on an event only advisory rows
/// cover (M11). Every check that verifies something itself passes `false`.
fn check_outcome(
    label: &str,
    check_name: &str,
    result: Result<CheckResult>,
    advisory: bool,
) -> CheckResult {
    match result {
        Ok(mut outcome) => {
            outcome.advisory = advisory;
            outcome
        }
        Err(e) => {
            log::warn!("{} check failed to run for {}: {:#}", check_name, label, e);
            CheckResult {
                event_label: label.to_string(),
                check_name: check_name.to_string(),
                expected: "check runs".to_string(),
                observed: format!("error: {:#}", e),
                pass: false,
                advisory,
            }
        }
    }
}

/// Entry point for `spike validate`.
pub fn run() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args = parse_validate_args()?;

    log::info!("spike validate");
    log::info!("  BAM:       {}", args.bam_path);
    log::info!("  Truth VCF: {}", args.truth_path);
    log::info!("  Reference: {}", args.ref_path);

    let truth_events = load_truth_events(&args.truth_path)?;
    log::info!("Loaded {} truth events", truth_events.len());

    let mut results: Vec<CheckResult> = Vec::new();

    // Per-event checks.
    let nearby = NearbyRecords::new(&truth_events);
    for event in &truth_events {
        check_event(&args, event, &nearby, &mut results);
    }

    // Global checks. All three read one sample, taken from the truth events'
    // own regions rather than from the head of the file (L15).
    match sample_event_regions(
        &args.bam_path,
        &args.ref_path,
        &truth_events,
        args.flank_bp,
        GLOBAL_SAMPLE_MAX,
    ) {
        Ok(sample) => {
            results.push(check_insert_size(&sample));
            results.push(check_dup_rate(&sample));
            results.push(check_mapq(&sample));
        }
        Err(e) => {
            // A sample that could not be taken is three failed checks, not
            // three checks that quietly pass on a default (M10, M11).
            for check_name in ["insert_size", "dup_rate", "mean_mapq"] {
                results.push(check_outcome(
                    GLOBAL_LABEL,
                    check_name,
                    Err(anyhow::anyhow!("{:#}", e)),
                    false,
                ));
            }
        }
    }

    // Print results.
    print_results(&results, args.json_output, args.strict)?;

    if let Some(message) = failure_message(&results, args.strict) {
        bail!("{}", message);
    }

    Ok(())
}

/// Report what spike recorded for one event: `SIM_RESIST` as `resistant` and
/// `SIM_DEPTH_FOLD` as `depth_fold`, each against the threshold spike warns
/// at.
///
/// Both numbers are read back from the truth VCF, not measured from the BAM,
/// so both rows are advisory: they say what spike saw at simulation time, and
/// they inherit whatever that measurement could not see. A record carrying
/// neither field -- an older spike's truth VCF -- gets neither row.
fn push_census_rows(label: &str, event: &TruthEvent, results: &mut Vec<CheckResult>) {
    if let Some(raw) = event.census.resist.as_deref() {
        results.push(census_row(label, RESISTANT, raw, census::WARN_ABOVE, 3));
    }
    if let Some(raw) = event.census.depth_fold.as_deref() {
        results.push(census_row(
            label,
            DEPTH_FOLD,
            raw,
            census::DEPTH_FOLD_WARN_ABOVE,
            2,
        ));
    }
}

/// One advisory row: `raw` read as a number, printed to `decimals` places and
/// passed iff it is at or below `warn_above`. The expectation is formatted
/// from `warn_above` itself, so moving the constant spike warns at moves the
/// printed expectation with it.
///
/// A value that cannot be read is a failed row, never a dropped one: a truth
/// VCF whose census is unreadable has not been checked. Its `observed` reads
/// `bad: <raw>` -- short on purpose, because the Observed column is 14
/// characters wide and a longer prefix would truncate away the very value the
/// reader needs to see.
fn census_row(
    label: &str,
    check_name: &str,
    raw: &str,
    warn_above: f64,
    decimals: usize,
) -> CheckResult {
    let (observed, pass) = match raw.parse::<f64>() {
        Ok(value) => (format!("{:.*}", decimals, value), value <= warn_above),
        Err(_) => (format!("bad: {}", raw), false),
    };
    CheckResult {
        event_label: label.to_string(),
        check_name: check_name.to_string(),
        expected: format!("<={:.*}", decimals, warn_above),
        observed,
        pass,
        advisory: true,
    }
}

/// The failures the run exits on, as the message it exits with, or `None` when
/// none of the rows it counts failed.
///
/// The advisory rows stay out of the count by default, against a total that
/// leaves them out too: a run that fails three of five real checks says so
/// whether or not advisory rows were printed beside them. `--strict` counts
/// every row instead.
fn failure_message(results: &[CheckResult], strict: bool) -> Option<String> {
    let n_total = counted_rows(results, strict).count();
    let n_fail = counted_rows(results, strict).filter(|r| !r.pass).count();
    if n_fail == 0 {
        return None;
    }
    Some(format!("{}/{} validation checks failed", n_fail, n_total))
}

/// The rows the exit status is computed over: the non-advisory ones by default,
/// every row when `strict`.
///
/// [`failure_message`] and [`print_results_json`]'s `counted_*` trio both read
/// this one filter, so `--json`'s summary cannot disagree with the exit status
/// printed beside it. A second copy of `strict || !r.advisory` is exactly how
/// the two would drift apart.
fn counted_rows<'a>(
    results: &'a [CheckResult],
    strict: bool,
) -> impl Iterator<Item = &'a CheckResult> + 'a {
    results.iter().filter(move |r| strict || !r.advisory)
}

/// The depth-ratio row's name: event depth over flanking depth, at the run's
/// own `--min-mapq`.
const COVERAGE_RATIO: &str = "coverage_ratio";

/// The advisory coverage row's name: `coverage_ratio` recomputed with no MAPQ
/// floor. 17 characters, so it fits the Check column's 18 without moving a
/// thing; `test_the_any_mapq_row_does_not_move_the_status_column` pins that.
const COVERAGE_ANY_MAPQ: &str = "coverage_any_mapq";

/// The pooled split-read row's name: the minimum required over both
/// breakpoint windows together.
const SPLIT_READS: &str = "split_reads";

/// The advisory split-read row's name: the same evidence, the same minimum,
/// required at **each** breakpoint (NF5).
///
/// 20 characters against the Check column's `{:<18}`, so this row's own later
/// columns sit two characters right of every other row's. That was weighed and
/// accepted when the row was specified: renaming it was not this run's call,
/// and widening the column would move the non-advisory rows' spacing, which is
/// forbidden. `{:<18}` pads short names and only overflows long ones, so no
/// other row is affected; `test_the_each_end_row_is_the_only_row_its_width_moves`
/// pins both halves of that.
const SPLIT_READS_EACH_END: &str = "split_reads_each_end";

/// At least this many split reads must join an event's two breakpoints --
/// pooled over both ends for `split_reads`, and at each end on its own for
/// `split_reads_each_end`. Reads with an SA tag occur anywhere; reads joining
/// these two points do not.
const MIN_SPLIT_READS: usize = 2;

/// The INS row's name: reads carrying the inserted sequence.
const INS_READS: &str = "ins_reads";

/// The INS row that decides: whether the reads spike made for the event are
/// in the BAM at POS, carrying its bases (RF13).
const INS_PLANTED: &str = "ins_planted";

/// Fewest of spike's own reads that must carry an insertion for
/// [`INS_PLANTED`] to pass. One: the sample's reads are left out, so there is no
/// background count to rise above (RF13's plan).
const MIN_PLANTED_READS: usize = 1;

/// Reference bases a junction probe takes on its flank side; the other 16 are
/// the haplotype's past the junction.
const PLANTED_PROBE_FLANK: usize = 15;

/// Length of a junction probe.
const PLANTED_PROBE_LEN: usize = 31;

/// Flank positions of a probe that may differ in a carrying read: the sample's
/// own SNPs, which spike writes onto the event copy's reference bases, and
/// sequencing errors. An inserted base may not differ at all: RF13's first
/// attempt let 2 substitutions fall anywhere, and a truth with the wrong 4
/// inserted bases matched.
const PLANTED_MAX_FLANK_MISMATCH: usize = 2;

/// The DEL row that decides: whether the reads spike made for the event are
/// in the BAM at its breakpoints, carrying its join (RF14).
const DEL_PLANTED: &str = "del_planted";

/// How far either side of each breakpoint [`DEL_PLANTED`] looks for spike's
/// reads: the windows `split_reads` reads, so the row sees at least the records
/// that one did.
const DEL_PLANTED_PAD: u64 = 500;

/// A deletion's junction probe, in both orientations and without a duplicate:
/// the haplotype's bases across the join, `left` (the last
/// [`PLANTED_PROBE_FLANK`] reference bases before START) then `right` (the
/// first bases from END on). No position is inserted, so the mask is all
/// false and [`carries_junction_probe`]'s tolerance covers every base.
fn del_junction_probes(left: &[u8], right: &[u8]) -> Vec<(Vec<u8>, Vec<bool>)> {
    let mut probe = left.to_ascii_uppercase();
    probe.extend_from_slice(&right.to_ascii_uppercase());
    let mask = vec![false; probe.len()];
    let mut reverse = probe.clone();
    crate::extract::reverse_complement(&mut reverse);
    let mut probes = vec![(probe, mask.clone())];
    if reverse != probes[0].0 {
        probes.push((reverse, mask));
    }
    probes
}

/// The read-name prefix spike gives the tiled reads of event `n`:
/// `simulate_event`'s `format!("ev{:04}", n)` and `tile_haplotype_reads`'s
/// `_hap_`. `simulate.rs` pins that the two agree.
pub(crate) fn planted_read_prefix(n: u32) -> String {
    format!("ev{:04}_hap_", n)
}

/// `N` if `id` is spike's truth ID `<prefix>N` (`sim_ins_N`, `sim_del_N`);
/// truth.rs numbers it with the same `i + 1` that names the event's reads.
fn sim_number(id: &str, prefix: &str) -> Option<u32> {
    let digits = id.strip_prefix(prefix)?;
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    digits.parse().ok()
}

/// An insertion's two junction probes, each with a mask marking the positions
/// that hold an inserted base, in both orientations and without duplicates.
///
/// `left` is the reference right up to the insertion point and `right` the
/// reference from it on, as the event's haplotype `left + inserted + right`
/// joins them. The probes are that haplotype's [`PLANTED_PROBE_LEN`] bases
/// across each junction, [`PLANTED_PROBE_FLANK`] of them on the reference
/// side: `[o-15, o+16)` and `[o+L-16, o+L+15)` for `o = left.len()`. For a
/// short insertion each probe reaches through it into the other flank.
fn ins_junction_probes(
    left: &[u8],
    inserted: &[u8],
    right: &[u8],
) -> Vec<(Vec<u8>, Vec<bool>)> {
    let mut hap: Vec<u8> = Vec::with_capacity(left.len() + inserted.len() + right.len());
    hap.extend_from_slice(left);
    hap.extend_from_slice(inserted);
    hap.extend_from_slice(right);
    let hap = hap.to_ascii_uppercase();
    let o = left.len();
    let l = inserted.len();
    let is_inserted = |i: usize| i >= o && i < o + l;
    let inside = PLANTED_PROBE_LEN - PLANTED_PROBE_FLANK;

    let mut probes: Vec<(Vec<u8>, Vec<bool>)> = Vec::new();
    for (a, b) in [
        (o.saturating_sub(PLANTED_PROBE_FLANK), o + inside),
        ((o + l).saturating_sub(inside), o + l + PLANTED_PROBE_FLANK),
    ] {
        let (a, b) = (a, b.min(hap.len()));
        let probe = hap[a..b].to_vec();
        let mask: Vec<bool> = (a..b).map(is_inserted).collect();
        let mut reverse = probe.clone();
        crate::extract::reverse_complement(&mut reverse);
        let reverse_mask: Vec<bool> = mask.iter().rev().copied().collect();
        for candidate in [(probe, mask), (reverse, reverse_mask)] {
            if !probes.contains(&candidate) {
                probes.push(candidate);
            }
        }
    }
    probes
}

/// True if some window of `seq` carries one of `probes`: every inserted
/// position equal, and at most [`PLANTED_MAX_FLANK_MISMATCH`] flank positions
/// not, case-insensitively. A deletion's probe marks no position inserted, so
/// the tolerance covers all of it.
fn carries_junction_probe(seq: &[u8], probes: &[(Vec<u8>, Vec<bool>)]) -> bool {
    let seq = seq.to_ascii_uppercase();
    probes.iter().any(|(probe, mask)| {
        seq.windows(probe.len()).any(|window| {
            let mut flank_mismatches = 0;
            for ((&got, &want), &inserted) in window.iter().zip(probe).zip(mask) {
                if got != want {
                    if inserted {
                        return false;
                    }
                    flank_mismatches += 1;
                    if flank_mismatches > PLANTED_MAX_FLANK_MISMATCH {
                        return false;
                    }
                }
            }
            true
        })
    })
}

/// At least this many reads must carry an insertion's evidence -- an alignment
/// that leaves the reference for [`INS_READS`], the inserted bases themselves
/// for [`INS_SEQUENCE`]. Two matches [`MIN_SPLIT_READS`] and for the same
/// reason: one clipped read is background anywhere, two at the same point are
/// not. Shared rather than written twice, so the two rows cannot drift apart.
const MIN_INS_READS: usize = 2;

/// The advisory INS row's name: reads carrying the bases the truth record's
/// own ALT names, not just an alignment that leaves the reference (CR9).
///
/// 12 characters against the Check column's `{:<18}`, so nothing moves;
/// `test_the_ins_sequence_row_does_not_move_the_status_column` pins that.
const INS_SEQUENCE: &str = "ins_sequence";

/// Longest probe k-mer [`INS_SEQUENCE`] takes out of an insertion's ALT. Long
/// enough to be specific in a read and short enough to sit inside one.
const INS_KMER_LEN: usize = 31;

/// Shortest probe k-mer [`INS_SEQUENCE`] will grade on. Below it the row is
/// not evaluable rather than graded: a given random 12-mer is expected in
/// about one 151 bp read in 10^5, an 8-mer in about one in 400, so a shorter
/// k-mer would be found in unedited reads and the count would mean nothing.
const MIN_INS_KMER_LEN: usize = 12;

/// How far either side of POS [`INS_SEQUENCE`] reads the reference looking for
/// its own probe k-mers. Finding one there makes the row not evaluable: an
/// unedited read would match it.
const INS_KMER_REF_PAD: u64 = 1_000;

/// How far either side of POS [`INS_SEQUENCE`] collects reads: one read length
/// on each side, so every read whose bases can reach the insertion point is
/// queried and nothing further is.
const INS_SEQUENCE_PAD: u64 = 150;

/// The small-variant row's name: the allele fraction read off the pileup.
const ALLELE_FREQ: &str = "allele_freq";

/// The "no check applies to this event" fallback row's name (M11).
const EVENT_CHECKED: &str = "event_checked";

/// The advisory census row's name for `SIM_RESIST`: the share of the reads over
/// the event spike could not edit, as spike recorded it (CR4).
const RESISTANT: &str = "resistant";

/// The advisory census row's name for `SIM_DEPTH_FOLD`: how far the donor's
/// depth departed from the depth the scaling assumed, as spike recorded it
/// (CR2).
const DEPTH_FOLD: &str = "depth_fold";

// The twelve check-name constants above name every row `check_event` pushes,
// and `print_usage` prints all twelve from these same constants: the six
// non-advisory ones
// (`COVERAGE_RATIO`, `DEL_PLANTED`, `SPLIT_READS` for DUP, INV and BND,
// `INS_PLANTED`, `ALLELE_FREQ`, `EVENT_CHECKED`) and the six advisory ones
// (`COVERAGE_ANY_MAPQ`, `SPLIT_READS_EACH_END`, `INS_READS` since RF13,
// `INS_SEQUENCE`, `RESISTANT`, `DEPTH_FOLD`), with `SPLIT_READS` advisory for a
// DEL since RF14. Renaming any of them moves the printed table and `--help`
// together, so neither can print one name while the other prints the old one.
// README.md is still edited by hand: it spells every row name out in prose (the
// table at "Which check covers which type", and a section each), and nothing in
// the build notices when a rename leaves it behind.
//
// The boundary of the shared-constant idiom is deliberate
// and stops there: `insert_size`, `dup_rate` and `mean_mapq` are `run()`'s
// global rows rather than `check_event`'s, and `check_ins_reads` and the
// allele-frequency path still spell their own names out where they build their
// `CheckResult`s. Finishing those too would edit the allele-frequency path,
// `check_ins_reads`, the three global checks and `run()`'s fallback list, for a
// cosmetic idiom, with five non-advisory rows' printed text as the stake.

/// Run every check that applies to one truth event, pushing one result per
/// check. An event no check applies to gets one failed "not evaluable"
/// result, so a truth VCF of such events cannot report all-PASS (M11).
fn check_event(
    args: &ValidateArgs,
    event: &TruthEvent,
    nearby: &NearbyRecords,
    results: &mut Vec<CheckResult>,
) {
    let label = format_event_label(event);
    let n_before = results.len();

    // Coverage ratio check (meaningful for DEL, DUP).
    if event.sv_type == "DEL" || event.sv_type == "DUP" {
        let r = check_coverage_ratio(
            &args.bam_path,
            &args.ref_path,
            event,
            args.flank_bp,
            args.min_mapq,
            COVERAGE_RATIO,
        );
        results.push(check_outcome(&label, COVERAGE_RATIO, r, false));

        // The same check, same window, no MAPQ floor. `coverage_ratio` cannot
        // see a read its own `--min-mapq` rejects, and those are exactly the
        // reads spike could not edit: its donor pool has the same floor, so an
        // AF=1 deletion over a window of low-MAPQ reads leaves them where they
        // were and still reports `observed 0.00, pass` (CR4). This row counts
        // every primary, non-duplicate, non-QC-fail read instead. It is
        // advisory and unconditional: at `--min-mapq 0` the two rows are
        // identical by construction, and both are still printed.
        let r = check_coverage_ratio(
            &args.bam_path,
            &args.ref_path,
            event,
            args.flank_bp,
            0,
            COVERAGE_ANY_MAPQ,
        );
        results.push(check_outcome(&label, COVERAGE_ANY_MAPQ, r, true));
    }

    // The counted DEL row: are spike's own reads for the event there, carrying
    // its join (RF14). Whether the aligner split them and named the partner is
    // the realism question `split_reads` below asks, and for a DEL that row is
    // advisory now. DUP, INV and BND were not measured and are unchanged.
    if event.sv_type == "DEL" {
        let r = check_del_planted(&args.bam_path, &args.ref_path, event);
        results.push(check_outcome(&label, DEL_PLANTED, r, false));
    }

    // Split reads joining the two breakpoints (not INS: its inserted
    // sequence has no second reference breakpoint).
    if matches!(event.sv_type.as_str(), "DEL" | "DUP" | "INV" | "BND") {
        // One pair of queries, two rows. `split_reads` requires its minimum
        // over the **union** of the two breakpoint windows, so an event whose
        // every split read sits at one end still passes while its `expected`
        // string reads as though each end contributed (NF5). The advisory
        // `split_reads_each_end` requires the same minimum at each end.
        //
        // Both rows are built from the same two name sets, so they can never
        // disagree about what evidence exists -- including when the BAM could
        // not be read at all, where the one error becomes both rows. The error
        // is re-wrapped rather than cloned (`anyhow::Error` is not `Clone`);
        // `{:#}` of the re-wrap is the same flattened chain the single row
        // printed before, so the pooled row's Observed did not move.
        let (pooled, each_end) = match check_split_reads(
            &args.bam_path,
            &args.ref_path,
            event,
            args.min_mapq,
            &label,
        ) {
            Ok((pooled, each_end)) => (Ok(pooled), Ok(each_end)),
            Err(e) => {
                // One error, two rows: flatten the chain once and hand each
                // row its own `anyhow::Error` over that one string.
                let flat = format!("{:#}", e);
                (
                    Err(anyhow::Error::msg(flat.clone())),
                    Err(anyhow::Error::msg(flat)),
                )
            }
        };
        results.push(check_outcome(&label, SPLIT_READS, pooled, event.sv_type == "DEL"));
        results.push(check_outcome(&label, SPLIT_READS_EACH_END, each_end, true));
    }

    // Reads carrying the inserted sequence (INS only: it is the one event
    // type with no second breakpoint and no reference span of its own).
    if event.sv_type == "INS" {
        // The counted row: are spike's own reads for the event there, carrying
        // its bases (RF13). Whether the aligner wrote them as a neat `I` is the
        // realism question, which a fixed CIGAR rule cannot answer: `ins_reads`
        // below failed 47 of 112 real HG002 insertions of 20-39 bp and passed
        // absent 50 bp ones at 12 of 200 empty sites, so it is advisory now.
        let r = check_ins_planted(&args.bam_path, &args.ref_path, event);
        results.push(check_outcome(&label, INS_PLANTED, r, false));

        let r = check_ins_reads(&args.bam_path, &args.ref_path, event, args.min_mapq);
        results.push(check_outcome(&label, INS_READS, r, true));

        // The same insertion, read rather than counted. `ins_reads` works from
        // the CIGAR alone and never looks at a base, so an insertion of
        // roughly the right length in roughly the right place passes it
        // whatever it spells (CR9). This row takes the inserted bases out of
        // the truth record's own ALT and looks for them in the reads.
        //
        // No row at all when the ALT recorded no sequence: a symbolic `<INS>`
        // is an older spike's truth VCF, written before CR7's fix put the
        // bases there, and that is not a failure of this run -- the same rule
        // T1 applied to a truth record carrying no census field.
        if let Some(inserted) = ins_alt_sequence(event) {
            let r = check_ins_sequence(
                &args.bam_path,
                &args.ref_path,
                event,
                inserted,
                args.min_mapq,
            );
            results.push(check_outcome(&label, INS_SEQUENCE, r, true));
        }
    }

    // Allele frequency (meaningful for SNPs/small variants).
    if event.sv_type == "SNP" && event.ref_allele.is_some() && event.alt_allele.is_some() {
        let r = check_allele_freq(&args.bam_path, &args.ref_path, event, nearby, args.min_mapq);
        results.push(check_outcome(&label, ALLELE_FREQ, r, false));
    }

    push_census_rows(&label, event, results);

    // An advisory row is not a check of the event, so it cannot stand in for
    // one: an INS-only truth VCF carrying SIM_RESIST must still report
    // `event_checked FAIL` (M11).
    let n_checked = results[n_before..].iter().filter(|r| !r.advisory).count();
    if n_checked == 0 {
        // No check applies to this event type -- an INS has no coverage
        // ratio, no split reads and no allele frequency -- so nothing about
        // the event was verified. An unverified event is a failed result, not
        // an absent one: an INS-only truth VCF used to score 3/3 PASS on the
        // three global checks alone (M11).
        log::warn!("no check applies to {}, so it was not evaluated", label);
        results.push(CheckResult {
            event_label: label.clone(),
            check_name: EVENT_CHECKED.to_string(),
            expected: "a check applies".to_string(),
            observed: format!("none for {}", event.sv_type),
            pass: false,
            advisory: false,
        });
    }

    // Log progress.
    log::info!("Checked: {}", label);
}

/// Parse validate-specific arguments from std::env::args().
fn parse_validate_args() -> Result<ValidateArgs> {
    let raw: Vec<String> = std::env::args().collect();
    parse_validate_args_from(&raw)
}

/// The flags of one `spike validate` command line, so that the parsing can be
/// tested without a process to run. `raw[0]` is the binary and `raw[1]` is
/// `validate`; the rest are flags.
///
/// **Not fully testable: `--help` prints the usage and exits the process with
/// status 0.** A test that passes `--help` here would end the test binary
/// mid-run and be read as a pass. Test every other flag, never that one.
fn parse_validate_args_from(raw: &[String]) -> Result<ValidateArgs> {
    if raw.len() < 3 {
        print_usage();
        bail!("missing required arguments");
    }

    let mut bam_path: Option<String> = None;
    let mut truth_path: Option<String> = None;
    let mut ref_path: Option<String> = None;
    let mut min_mapq: u8 = 20;
    let mut flank_bp: u64 = 5000;
    let mut json_output = false;
    let mut strict = false;

    let args = &raw[2..];
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--bam" | "-b" => {
                i += 1;
                bam_path = Some(
                    args.get(i)
                        .cloned()
                        .ok_or_else(|| anyhow::anyhow!("--bam requires a value"))?,
                );
            }
            "--truth" | "-t" => {
                i += 1;
                truth_path = Some(
                    args.get(i)
                        .cloned()
                        .ok_or_else(|| anyhow::anyhow!("--truth requires a value"))?,
                );
            }
            "--reference" | "-r" => {
                i += 1;
                ref_path = Some(
                    args.get(i)
                        .cloned()
                        .ok_or_else(|| anyhow::anyhow!("--reference requires a value"))?,
                );
            }
            "--min-mapq" => {
                i += 1;
                min_mapq = args
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("--min-mapq requires a value"))?
                    .parse()
                    .context("--min-mapq must be a number")?;
            }
            "--flank" => {
                i += 1;
                flank_bp = args
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("--flank requires a value"))?
                    .parse()
                    .context("--flank must be a number")?;
            }
            "--json" => {
                json_output = true;
            }
            "--strict" => {
                strict = true;
            }
            "--help" | "-h" => {
                print_usage();
                std::process::exit(0);
            }
            other => {
                bail!("unknown argument: {}", other);
            }
        }
        i += 1;
    }

    Ok(ValidateArgs {
        bam_path: bam_path.ok_or_else(|| anyhow::anyhow!("--bam is required"))?,
        truth_path: truth_path.ok_or_else(|| anyhow::anyhow!("--truth is required"))?,
        ref_path: ref_path.ok_or_else(|| anyhow::anyhow!("--reference is required"))?,
        min_mapq,
        flank_bp,
        json_output,
        strict,
    })
}

fn print_usage() {
    eprintln!("Usage: spike validate --bam <BAM> --truth <VCF> --reference <FASTA> [OPTIONS]");
    eprintln!();
    eprintln!("Options:");
    eprintln!("  --bam, -b        Simulated BAM/CRAM file (required)");
    eprintln!("  --truth, -t      Truth VCF from spike (required)");
    eprintln!("  --reference, -r  Reference FASTA with .fai index (required)");
    eprintln!("  --min-mapq       Minimum MAPQ for counting reads (default: 20)");
    eprintln!("  --flank          Flanking bp for coverage comparison (default: 5000)");
    eprintln!("  --json           Output JSON instead of text table");
    eprintln!("  --strict         Count the advisory checks in the exit status");
    eprintln!("  --help, -h       Show this help");
    eprintln!();
    eprintln!("Checks, by truth-event type:");
    eprintln!("  DEL              {}, {} (spike's own reads for the", COVERAGE_RATIO, DEL_PLANTED);
    eprintln!("                   event, named after the truth ID sim_del_N, carrying");
    eprintln!("                   its join at the breakpoints)");
    eprintln!("  DUP              {}, {}", COVERAGE_RATIO, SPLIT_READS);
    eprintln!("  INV, BND         {}", SPLIT_READS);
    eprintln!("  INS              {} (spike's own reads for the event,", INS_PLANTED);
    eprintln!("                   named after the truth ID sim_ins_N, carrying its");
    eprintln!("                   inserted bases across a junction at POS)");
    eprintln!(
        "  SNP, small indel {} (a substitution from the pileup; a small",
        ALLELE_FREQ
    );
    eprintln!("  and MNV          indel from each read's bases over it: nearer the");
    eprintln!("                   reference with the indel made carries it, nearer");
    eprintln!("                   the reference spans it, a tie does not vote; other");
    eprintln!("                   truth records nearby are tried in with it -- only");
    eprintln!("                   reads reaching 10bp past the indel's repeat on both");
    eprintln!("                   sides vote; an MNV");
    eprintln!("                   from the whole alt run, read by read). Each read");
    eprintln!("                   pair votes once; a pair whose mates disagree");
    eprintln!("                   does not vote");
    eprintln!("  every event      insert_size, dup_rate, mean_mapq, over the whole sample");
    eprintln!();
    eprintln!("Advisory checks, printed beside the checks above:");
    eprintln!("  DEL, DUP         {} (the depth ratio again with no", COVERAGE_ANY_MAPQ);
    eprintln!("                   MAPQ floor)");
    eprintln!("  DEL              {} (split reads joining the two", SPLIT_READS);
    eprintln!("                   breakpoints: how the aligner wrote the junction)");
    eprintln!(
        "  DEL, DUP, INV    {} (the same evidence and the same",
        SPLIT_READS_EACH_END
    );
    eprintln!("  and BND          minimum required at each breakpoint, not pooled)");
    eprintln!("  INS              {} (reads whose alignment leaves the", INS_READS);
    eprintln!("                   reference at POS -- an I operation, or a soft clip");
    eprintln!("                   once the insertion is 50 bp or longer)");
    eprintln!("  INS              {} (reads carrying the bases the truth", INS_SEQUENCE);
    eprintln!("                   record's own ALT names)");
    eprintln!(
        "  every event      {} and {}, when the truth record carries",
        RESISTANT, DEPTH_FOLD
    );
    eprintln!("                   SIM_RESIST or SIM_DEPTH_FOLD: what spike recorded at");
    eprintln!("                   simulation time, read back from the truth VCF rather");
    eprintln!("                   than measured from the BAM");
    eprintln!();
    eprintln!("An advisory row prints its own PASS or FAIL, marked `(advisory)` in the");
    eprintln!("Status column, and is left out of the exit status unless --strict is");
    eprintln!("given.");
    eprintln!();
    eprintln!("A check that cannot run is a FAILED check, never a silent pass. Exit 0");
    eprintln!("means every check the exit status counts ran and passed: every");
    eprintln!("non-advisory check by default, and every check including the advisory");
    eprintln!("ones under --strict. An event type no check covers (e.g. SVTYPE=CNV) is");
    eprintln!("reported as `{} FAIL`.", EVENT_CHECKED);
}

// ---------------------------------------------------------------------------
// Truth VCF parsing
// ---------------------------------------------------------------------------

/// A symbolic SV record's `END`, checked against its own POS.
///
/// An `END` at or before POS leaves the event region empty, and
/// `count_depth_in_region` answers a zero-length region with a hard-coded
/// `0.0`. A DEL at a high `SIM_VAF` then PASSes `coverage_ratio` -- expected
/// 0.10, observed 0.00 -- on a region no query ever read. Refuse the record
/// rather than grade the run on it.
fn sv_end(info: &str, chrom: &str, id: &str, start: u64) -> Result<u64> {
    let end = parse_info_u64(info, "END").unwrap_or(start + 1);
    if end <= start {
        bail!(
            "truth record {} at {}:{} has END={}, at or before its own POS: the event \
             region is empty, so no check can read it and a coverage ratio taken over \
             it means nothing. Fix the record's END.",
            id,
            chrom,
            start,
            end,
        );
    }
    Ok(end)
}

/// Load truth events from a spike-produced truth VCF.
fn load_truth_events(path: &str) -> Result<Vec<TruthEvent>> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read truth VCF: {}", path))?;

    let reader = BufReader::new(content.as_bytes());
    let mut events = Vec::new();
    let mut seen_bnd_ids: std::collections::HashSet<String> = std::collections::HashSet::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            continue;
        }

        let chrom = fields[0].to_string();
        let vcf_pos: u64 = match fields[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let id = fields[2].to_string();
        let ref_col = fields[3];
        let alt_col = fields[4];
        let info = fields[7];

        let sv_type_str = parse_info_field(info, "SVTYPE");
        let sim_vaf = parse_info_field(info, "SIM_VAF")
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(0.5);
        let gene = parse_info_field(info, "SIM_GENE")
            .unwrap_or("unknown")
            .to_string();
        // One record takes one arm below, so each arm may move this.
        let census_info = CensusInfo::from_info(info);

        match sv_type_str {
            Some("DEL") => {
                let end = sv_end(info, &chrom, &id, vcf_pos)?;
                let partner = Some((chrom.clone(), end));
                events.push(TruthEvent {
                    chrom,
                    start: vcf_pos, // VCF POS for SV = 0-based start
                    end,
                    sv_type: "DEL".to_string(),
                    partner,
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: None,
                    alt_allele: None,
                    ins_len: None,
                    ins_alt: None,
                    sim_number: sim_number(&id, "sim_del_"),
                    census: census_info,
                });
            }
            Some("DUP") => {
                let end = sv_end(info, &chrom, &id, vcf_pos)?;
                let partner = Some((chrom.clone(), end));
                events.push(TruthEvent {
                    chrom,
                    start: vcf_pos,
                    end,
                    sv_type: "DUP".to_string(),
                    partner,
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: None,
                    alt_allele: None,
                    ins_len: None,
                    ins_alt: None,
                    sim_number: None,
                    census: census_info,
                });
            }
            Some("INV") => {
                let end = sv_end(info, &chrom, &id, vcf_pos)?;
                let partner = Some((chrom.clone(), end));
                events.push(TruthEvent {
                    chrom,
                    start: vcf_pos,
                    end,
                    sv_type: "INV".to_string(),
                    partner,
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: None,
                    alt_allele: None,
                    ins_len: None,
                    ins_alt: None,
                    sim_number: None,
                    census: census_info,
                });
            }
            Some("INS") => {
                // SVLEN is the inserted length; spike writes it positive, but
                // other producers sign it, so take the magnitude.
                let ins_len = parse_info_field(info, "SVLEN")
                    .and_then(|v| v.parse::<i64>().ok())
                    .map(|v| v.unsigned_abs());
                events.push(TruthEvent {
                    chrom,
                    start: vcf_pos,
                    end: vcf_pos,
                    sv_type: "INS".to_string(),
                    partner: None,
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: None,
                    alt_allele: None,
                    ins_len,
                    ins_alt: Some(alt_col.as_bytes().to_vec()),
                    sim_number: sim_number(&id, "sim_ins_"),
                    census: census_info,
                });
            }
            Some("BND") => {
                if seen_bnd_ids.contains(&id) {
                    continue;
                }
                seen_bnd_ids.insert(id.clone());
                if let Some(mate_id) = parse_info_field(info, "MATEID") {
                    seen_bnd_ids.insert(mate_id.to_string());
                }
                // BND POS is 1-based breakpoint in truth VCF.
                events.push(TruthEvent {
                    chrom,
                    start: vcf_pos.saturating_sub(1), // 1-based → 0-based
                    end: vcf_pos,
                    sv_type: "BND".to_string(),
                    partner: bnd_partner(alt_col),
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: None,
                    alt_allele: None,
                    ins_len: None,
                    ins_alt: None,
                    sim_number: None,
                    census: census_info,
                });
            }
            Some(other) => {
                // An SVTYPE no check covers. It used to fall into the small
                // variant arm below, where a symbolic ALT like `<CNV>` is
                // longer than one base and took `check_allele_freq`'s indel
                // exit -- one silent PASS per unrecognised type. Keep the
                // type: `check_event` then reports it as `event_checked`,
                // which is a failed result, so an unknown type is visible.
                let end = sv_end(info, &chrom, &id, vcf_pos)?;
                log::warn!(
                    "truth record {} at {}:{} has SVTYPE={}, which no check covers",
                    id,
                    chrom,
                    vcf_pos,
                    other
                );
                events.push(TruthEvent {
                    chrom,
                    start: vcf_pos,
                    end,
                    sv_type: other.to_string(),
                    partner: None,
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: None,
                    alt_allele: None,
                    ins_len: None,
                    ins_alt: None,
                    sim_number: None,
                    census: census_info,
                });
            }
            None => {
                // No SVTYPE: small variant (SNP/indel).
                // VCF POS is 1-based for small variants in truth VCF.
                let pos_0based = vcf_pos.saturating_sub(1);
                let ref_allele = ref_col.as_bytes().to_vec();
                let alt_allele = alt_col.as_bytes().to_vec();
                let end = pos_0based + ref_allele.len() as u64;
                events.push(TruthEvent {
                    chrom,
                    start: pos_0based,
                    end,
                    sv_type: "SNP".to_string(),
                    partner: None,
                    expected_vaf: sim_vaf,
                    gene,
                    ref_allele: Some(ref_allele),
                    alt_allele: Some(alt_allele),
                    ins_len: None,
                    ins_alt: None,
                    sim_number: None,
                    census: census_info,
                });
            }
        }
    }

    Ok(events)
}

// ---------------------------------------------------------------------------
// Validation checks
// ---------------------------------------------------------------------------

/// Check coverage ratio: event depth vs flanking depth.
///
/// `min_mapq` is the floor every depth here is counted at, and `check_name`
/// names the row it produces, so one body serves both coverage rows: the real
/// `coverage_ratio` at the run's own `--min-mapq`, and the advisory
/// [`COVERAGE_ANY_MAPQ`] at a floor of 0. Everything else -- the window, the
/// flank-averaging rule, the expected ratio and its tolerance -- is shared by
/// construction rather than by agreement, so the two rows can only ever differ
/// in the reads the floor admits.
///
/// The row's advisory standing is deliberately **not** a parameter here.
/// [`check_outcome`] stamps it on every row it wraps -- the row this function
/// returns and the "check runs" row an errored check leaves alike -- so one
/// flag at the call site governs both. A second copy threaded through here
/// could only ever agree with it or be silently overwritten by it, and nothing
/// would notice which (F1).
fn check_coverage_ratio(
    bam_path: &str,
    ref_path: &str,
    event: &TruthEvent,
    flank_bp: u64,
    min_mapq: u8,
    check_name: &str,
) -> Result<CheckResult> {
    let label = format_event_label(event);

    // Compute depth in the event region.
    let event_depth = count_depth_in_region(
        bam_path,
        ref_path,
        &event.chrom,
        event.start,
        event.end,
        min_mapq,
    )?;

    // Compute depth in left and right flanking regions.
    let left_start = event.start.saturating_sub(flank_bp);
    let left_end = event.start;
    let right_start = event.end;
    let right_end = event.end.saturating_add(flank_bp);

    let left_depth = count_depth_in_region(
        bam_path,
        ref_path,
        &event.chrom,
        left_start,
        left_end,
        min_mapq,
    )?;
    let right_depth = count_depth_in_region(
        bam_path,
        ref_path,
        &event.chrom,
        right_start,
        right_end,
        min_mapq,
    )?;

    // Average only flanking regions that have non-zero length (handles events
    // near chromosome start/end where one flank may be empty).
    let flank_depth = match (left_start < left_end, right_start < right_end) {
        (true, true) => (left_depth + right_depth) / 2.0,
        (true, false) => left_depth,
        (false, true) => right_depth,
        (false, false) => 0.0,
    };

    Ok(coverage_ratio_result(
        label,
        check_name,
        &event.sv_type,
        event.expected_vaf,
        event_depth,
        flank_depth,
    ))
}

/// Judge an event's depth against its flanks, as the row `check_name`.
///
/// The name is the caller's, because the same judgement serves the two
/// coverage rows (see [`check_coverage_ratio`]); the expectations and the
/// tolerance below are not, because the two rows compare their depths against
/// exactly the same thing. The rows are built non-advisory and every
/// production caller reaches them through [`check_outcome`], which is the one
/// place that flag is set (F1).
fn coverage_ratio_result(
    label: String,
    check_name: &str,
    sv_type: &str,
    expected_vaf: f64,
    event_depth: f64,
    flank_depth: f64,
) -> CheckResult {
    if flank_depth < 1.0 {
        return CheckResult {
            event_label: label,
            check_name: check_name.to_string(),
            expected: "N/A".to_string(),
            observed: "no flanking coverage".to_string(),
            pass: false, // can't evaluate: don't report it as a pass
            advisory: false, // stamped by check_outcome
        };
    }

    let ratio = event_depth / flank_depth;

    // Expected ratio depends on SV type and VAF.
    let (expected_str, pass) = match sv_type {
        "DEL" => {
            // Spike suppresses reads at rate VAF → expected ratio = 1 - VAF.
            let expected_ratio = 1.0 - expected_vaf;
            let tolerance = 0.3;
            let pass = (ratio - expected_ratio).abs() < tolerance;
            (format!("{:.2}", expected_ratio), pass)
        }
        "DUP" => {
            // Spike adds depth copies at rate VAF → expected ratio = 1 + VAF.
            let expected_ratio = 1.0 + expected_vaf;
            let tolerance = 0.3;
            let pass = (ratio - expected_ratio).abs() < tolerance;
            (format!("{:.2}", expected_ratio), pass)
        }
        // Only DEL and DUP are routed here. Any other type has no expected
        // ratio, so it gets no verdict rather than a loose pass (NF6).
        _ => ("N/A".to_string(), false),
    };

    CheckResult {
        event_label: label,
        check_name: check_name.to_string(),
        expected: expected_str,
        observed: format!("{:.2}", ratio),
        pass,
        advisory: false, // stamped by check_outcome
    }
}

/// Check for split reads joining the event's two breakpoints: reads at one
/// breakpoint whose SA:Z alignment lands at the other.
///
/// Returns **both** split-read rows, `(split_reads, split_reads_each_end)`,
/// out of the same two queries. One pair of queries rather than two is not an
/// optimisation: it is what stops the two rows disagreeing about what evidence
/// exists, since neither can see a read the other did not (NF5). Neither row
/// is stamped advisory here -- `check_outcome` at the call site is the only
/// place that flag is set.
///
/// `label` is the caller's own, not recomputed: `check_event` already holds it
/// for its `check_outcome` calls, and one label now feeds two rows.
fn check_split_reads(
    bam_path: &str,
    ref_path: &str,
    event: &TruthEvent,
    min_mapq: u8,
    label: &str,
) -> Result<(CheckResult, CheckResult)> {
    let pad = 500u64;
    let partner = event
        .partner
        .clone()
        .with_context(|| format!("{}: no partner breakpoint for split reads", label))?;
    let here = (event.chrom.clone(), event.start);

    let at_here = split_reads_to_partner(
        bam_path,
        ref_path,
        &here.0,
        here.1.saturating_sub(pad),
        here.1.saturating_add(pad),
        min_mapq,
        &partner,
        pad,
    )?;
    let at_partner = split_reads_to_partner(
        bam_path,
        ref_path,
        &partner.0,
        partner.1.saturating_sub(pad),
        partner.1.saturating_add(pad),
        min_mapq,
        &here,
        pad,
    )?;

    Ok(split_read_rows(label, &partner, at_here, at_partner))
}

/// The two split-read rows for one event, from the two sets of read names its
/// breakpoint windows gave.
///
/// The pooled row counts the **union** of the two sets, exactly as it always
/// has: a read name seen at both ends is one read, and it votes once. The
/// per-end row takes each set's own size and requires the minimum of each,
/// which is why a junction whose every split read sits at one breakpoint fails
/// it while the pooled row still passes.
fn split_read_rows(
    label: &str,
    partner: &(String, u64),
    at_here: HashSet<String>,
    at_partner: HashSet<String>,
) -> (CheckResult, CheckResult) {
    let (n_here, n_partner) = (at_here.len(), at_partner.len());
    let mut pooled = at_here;
    pooled.extend(at_partner);

    (
        CheckResult {
            event_label: label.to_string(),
            check_name: SPLIT_READS.to_string(),
            expected: format!(">={} joining {}:{}", MIN_SPLIT_READS, partner.0, partner.1 + 1),
            observed: format!("{}", pooled.len()),
            pass: pooled.len() >= MIN_SPLIT_READS,
            advisory: false, // stamped by check_outcome
        },
        CheckResult {
            event_label: label.to_string(),
            check_name: SPLIT_READS_EACH_END.to_string(),
            expected: format!(">={} at each end", MIN_SPLIT_READS),
            // `print_results_text` truncates Observed at 14 characters, which
            // this pair reaches only once the two counts run to 14 digits
            // between them. Measured against `truncate` itself, not reasoned
            // about: `100000/1000000` (13 digits, 14 characters) is printed
            // whole, `1000000/1000000` (14 digits, 15 characters) prints as
            // `1000000/100...` and drops the partner count silently, and so
            // does the lopsided `9/1000000000000`. That is a million distinct
            // read names at each end of one 1 kb window, so it is unreachable,
            // and `--json` prints the pair untruncated regardless. Nobody need
            // re-derive this.
            observed: format!("{}/{}", n_here, n_partner),
            pass: n_here >= MIN_SPLIT_READS && n_partner >= MIN_SPLIT_READS,
            advisory: false, // stamped by check_outcome
        },
    )
}

/// Longest clip an insertion of any size is required to produce. A read is at
/// most a few hundred bases, so a 1 kb insertion still only clips part of one;
/// past this the threshold would reject the event's own reads. It is also the
/// length at which a soft clip starts counting as evidence at all: see
/// [`cigar_shows_insertion_near`].
const INS_MAX_EVIDENCE_LEN: u64 = 50;

/// Check that reads at an INS breakpoint carry the inserted sequence.
///
/// An insertion has no second reference breakpoint and no reference span, so
/// neither the coverage-ratio nor the split-read check applies to it. What it
/// does leave behind is sequence that does not fit the reference at one point,
/// which the aligner records as an `I` operation or a soft clip. spike writes
/// `SVTYPE=INS` itself (`truth.rs`), so without this check `spike validate`
/// could never exit 0 on spike's own output.
fn check_ins_reads(
    bam_path: &str,
    ref_path: &str,
    event: &TruthEvent,
    min_mapq: u8,
) -> Result<CheckResult> {
    // How far from POS the alignment may leave the reference. An aligner
    // places the boundary within a few bases of the insertion point; 100 bp
    // covers that without reaching the next feature.
    const PAD: u64 = 100;
    let label = format_event_label(event);
    let Some(ins_len) = event.ins_len.filter(|&n| n > 0) else {
        bail!("{}: no SVLEN, so there is no insertion length to look for", label);
    };
    // A short insertion fits inside a read as an `I` operation of its own
    // length; a long one is clipped, and only part of it is in any one read.
    let min_len = ins_len.min(INS_MAX_EVIDENCE_LEN);

    let names = reads_with_inserted_sequence(
        bam_path, ref_path, &event.chrom, event.start, PAD, min_len, min_mapq,
    )?;

    Ok(CheckResult {
        event_label: label,
        check_name: "ins_reads".to_string(),
        expected: format!(
            ">={} reads with >={}bp inserted at {}:{}",
            MIN_INS_READS,
            min_len,
            event.chrom,
            // `load_truth_events` reads an INS as `start: vcf_pos`, so `start`
            // is already the POS the truth record names.
            event.start
        ),
        observed: format!("{}", names.len()),
        pass: names.len() >= MIN_INS_READS,
        advisory: false,
    })
}

/// The bases one INS truth record says were inserted, or `None` when it
/// recorded none.
///
/// The ALT column is a literal sequence when it does not begin with `<` and is
/// longer than one base; its first base is the anchor VCF requires, so the
/// inserted sequence is everything after it. A symbolic ALT (`<INS>`) and a
/// bare anchor both mean the bases were not recorded, and the caller pushes no
/// row for them.
fn ins_alt_sequence(event: &TruthEvent) -> Option<&[u8]> {
    let alt = event.ins_alt.as_deref()?;
    if alt.len() <= 1 || alt.starts_with(b"<") {
        return None;
    }
    Some(&alt[1..])
}

/// The probe k-mers of an insertion, in both orientations and without
/// duplicates.
///
/// The **first** and **last** `k` bases, never the middle: a read anchored
/// left of POS carries the insertion's beginning and one anchored right of it
/// carries its end, while for an insertion longer than a read no read contains
/// the middle at all -- a 2000 bp insertion's middle k-mer sits 1000 bases in,
/// far past the reach of a 151 bp read. Both orientations, because a read's
/// stored sequence is in reference orientation but the insertion may be read
/// from either side.
///
/// One set serves both uses -- the reads and the reference guard -- so the
/// guard covers exactly the k-mers a read is tested against and the two cannot
/// disagree. For `inserted.len() == k` the first and last k-mers are the same
/// sequence, and a palindromic k-mer is its own reverse complement, so the
/// duplicates are dropped rather than probed twice.
fn alt_probe_kmers(inserted: &[u8], k: usize) -> Vec<Vec<u8>> {
    let upper = inserted.to_ascii_uppercase();
    let mut probes: Vec<Vec<u8>> = Vec::new();
    for kmer in [&upper[..k], &upper[upper.len() - k..]] {
        let mut reverse = kmer.to_vec();
        crate::extract::reverse_complement(&mut reverse);
        for candidate in [kmer.to_vec(), reverse] {
            if !probes.contains(&candidate) {
                probes.push(candidate);
            }
        }
    }
    probes
}

/// True if `bases` holds any of `probes`, case-insensitively.
fn holds_a_probe(bases: &[u8], probes: &[Vec<u8>]) -> bool {
    let upper = bases.to_ascii_uppercase();
    probes
        .iter()
        .any(|probe| upper.windows(probe.len()).any(|w| w == probe.as_slice()))
}

/// Check that reads at an INS breakpoint carry the bases the truth record's
/// own ALT names (CR9).
///
/// [`check_ins_reads`] beside this row counts reads whose *alignment* leaves
/// the reference near POS and reads no base at all, so an insertion of roughly
/// the right length in roughly the right place passes it whatever it spells.
/// CR7's fix put the inserted sequence into the truth VCF's ALT, so it can now
/// be looked for: this is the first check here that verifies inserted sequence
/// rather than a CIGAR.
///
/// `inserted` is the ALT past its anchor base, from [`ins_alt_sequence`]. Two
/// probe k-mers of `k = min(inserted.len(), INS_KMER_LEN)` bases are taken from
/// its two ends and a read supports the insertion if its bases hold either, in
/// either orientation. The reads are the distinct names over
/// `POS +/- INS_SEQUENCE_PAD` that [`usable_alignment`] admits -- the same
/// filter `ins_reads` applies, reused rather than restated -- and the floor is
/// the shared [`MIN_INS_READS`].
///
/// Two cases the row cannot answer, each a **failed** not-evaluable row rather
/// than a silent pass (M10):
///
/// - `k < MIN_INS_KMER_LEN`: a k-mer that short is not specific enough inside a
///   read, so a count of matches would say nothing about this insertion.
/// - either probe k-mer already in the reference within [`INS_KMER_REF_PAD`] of
///   POS: an unedited read would match it, so the count would mean nothing. A
///   random insertion makes that vanishingly unlikely; an explicit one copied
///   from nearby sequence does not.
///
/// The row is built non-advisory here and stamped by [`check_outcome`], which
/// is the single source of that flag (T2's F1).
fn check_ins_sequence(
    bam_path: &str,
    ref_path: &str,
    event: &TruthEvent,
    inserted: &[u8],
    min_mapq: u8,
) -> Result<CheckResult> {
    let label = format_event_label(event);
    // `load_truth_events` reads an INS as `start: vcf_pos`, so `start` is
    // already the POS the truth record names -- the same position
    // `check_ins_reads` scans around.
    let pos = event.start;
    let k = inserted.len().min(INS_KMER_LEN);
    if k < MIN_INS_KMER_LEN {
        return Ok(event_not_evaluable(
            &label,
            INS_SEQUENCE,
            &format!(">={} with a >={}bp kmer", MIN_INS_READS, MIN_INS_KMER_LEN),
            &format!("alt is {}bp", inserted.len()),
            &format!(
                "the ALT records {} inserted bases, fewer than the {} a probe k-mer \
                 needs to be specific inside a read",
                inserted.len(),
                MIN_INS_KMER_LEN
            ),
        ));
    }
    let probes = alt_probe_kmers(inserted, k);
    let expected = format!(">={} with a {}bp alt kmer", MIN_INS_READS, k);

    // The reference first: a probe the reference already holds near POS would
    // be matched by reads spike never touched.
    let (_, window) = crate::reference::fetch_window(
        ref_path,
        &event.chrom,
        pos.saturating_sub(INS_KMER_REF_PAD),
        pos + INS_KMER_REF_PAD,
    )?;
    if holds_a_probe(&window, &probes) {
        return Ok(event_not_evaluable(
            &label,
            INS_SEQUENCE,
            &expected,
            "kmer in ref",
            &format!(
                "the reference within {}bp of {}:{} already holds one of the ALT's own \
                 {}bp k-mers, so an unedited read would match it",
                INS_KMER_REF_PAD, event.chrom, pos, k
            ),
        ));
    }

    let mut names: HashSet<String> = HashSet::new();
    for_each_alignment(
        bam_path,
        ref_path,
        &event.chrom,
        pos.saturating_sub(INS_SEQUENCE_PAD),
        pos + INS_SEQUENCE_PAD,
        min_mapq,
        &mut |name, _align_start, _ops, seq| {
            if holds_a_probe(seq, &probes) {
                names.insert(String::from_utf8_lossy(name).into_owned());
            }
        },
    )?;

    Ok(CheckResult {
        event_label: label,
        check_name: INS_SEQUENCE.to_string(),
        expected,
        observed: format!("{}", names.len()),
        pass: names.len() >= MIN_INS_READS,
        advisory: false, // stamped by check_outcome
    })
}

/// Check that the reads spike made for an insertion are in the BAM at its
/// position, carrying its bases (RF13).
///
/// The question is whether the event was planted, not how the aligner wrote
/// it down: [`check_ins_reads`] asks the second, and the realism probe found it
/// failing 47 of 112 real HG002 insertions of 20-39 bp and passing absent ones
/// at empty sites. So this row looks only at spike's own reads for the event,
/// named [`planted_read_prefix`] of the truth ID's `N`, over POS +/-
/// [`INS_SEQUENCE_PAD`]. Every such record counts but a secondary or
/// supplementary one: its MAPQ, and duplicate, QC-fail and unmapped flags, are
/// the aligner's verdict, not the question. A read carries the event when its
/// bases hold one of [`ins_junction_probes`] as [`carries_junction_probe`] reads it.
///
/// No read of the sample's own can count, whatever it spells, so there is no
/// background to rise above and [`MIN_PLANTED_READS`] is 1. Two cases are a
/// failed not-evaluable row rather than a silent pass (M10): an ID that is not
/// `sim_ins_N`, which leaves spike's reads unknown, and an ALT with no bases.
fn check_ins_planted(bam_path: &str, ref_path: &str, event: &TruthEvent) -> Result<CheckResult> {
    let label = format_event_label(event);
    let Some(n) = event.sim_number else {
        return Ok(event_not_evaluable(
            &label,
            INS_PLANTED,
            &format!(">={} of spike's reads carrying", MIN_PLANTED_READS),
            "no sim_ins_N id",
            "the truth record's ID is not spike's sim_ins_N, so which reads are spike's \
             for this event is unknown",
        ));
    };
    let expected = format!(">={} ev{:04} read carrying", MIN_PLANTED_READS, n);
    let Some(inserted) = ins_alt_sequence(event) else {
        return Ok(event_not_evaluable(
            &label,
            INS_PLANTED,
            &expected,
            "alt has no bases",
            "the ALT records no inserted bases, so there is nothing to look for",
        ));
    };

    // `load_truth_events` reads an INS as `start: vcf_pos`, which is the 0-based
    // position spike inserts before: the haplotype is `ref[..pos] + INS +
    // ref[pos..]`.
    let pos = event.start;
    let reach = PLANTED_PROBE_LEN as u64;
    let (_, left) =
        crate::reference::fetch_window(ref_path, &event.chrom, pos.saturating_sub(reach), pos)?;
    let (_, right) = crate::reference::fetch_window(ref_path, &event.chrom, pos, pos + reach)?;
    let probes = ins_junction_probes(&left, inserted, &right);
    let prefix = planted_read_prefix(n);

    let mut names: HashSet<String> = HashSet::new();
    scan_region(
        bam_path,
        ref_path,
        &event.chrom,
        pos.saturating_sub(INS_SEQUENCE_PAD),
        pos + INS_SEQUENCE_PAD,
        &|flags, _mapq| !(flags.is_secondary() || flags.is_supplementary()),
        &mut |name, _align_start, _ops, seq| {
            if name.starts_with(prefix.as_bytes()) && carries_junction_probe(seq, &probes) {
                names.insert(String::from_utf8_lossy(name).into_owned());
            }
        },
    )?;

    Ok(CheckResult {
        event_label: label,
        check_name: INS_PLANTED.to_string(),
        expected,
        observed: format!("{}", names.len()),
        pass: names.len() >= MIN_PLANTED_READS,
        advisory: false, // stamped by check_outcome
    })
}

/// Check that the reads spike made for a deletion are in the BAM at its
/// breakpoints, carrying its join (RF14).
///
/// [`check_ins_planted`]'s question for a DEL. `split_reads` asks whether the
/// aligner split the junction reads and named the partner; on real SV sites it
/// often does not (the realism probe: 12 of the pipeline's 20 real HG002
/// deletions fail it in HG002's own BAM). So this row looks only at spike's own
/// reads for the event, named [`planted_read_prefix`] of the truth ID's `N`,
/// over START and END +/- [`DEL_PLANTED_PAD`]. Every such record counts but a
/// secondary or supplementary one. A read carries the event when its bases hold
/// [`del_junction_probes`] within [`PLANTED_MAX_FLANK_MISMATCH`] substitutions.
///
/// Measured on RF14's fresh sites: 33 of 33 correct deletions at VAF 0.5 carry
/// (8 to 35 reads); at VAF 0.1, 2 of 30 correct ones have no carrier by chance,
/// where `split_reads` fails 14. Where a real deletion's join is already
/// reference sequence (repeat units), unedited reads spell it too; the name
/// filter keeps them out, but the row then shows spike's reads are there, not
/// that the bases were removed -- `coverage_ratio` measures that.
///
/// Not evaluable, which fails (M10): an ID that is not `sim_del_N`, and a START
/// with fewer than [`PLANTED_PROBE_FLANK`] bases before it.
fn check_del_planted(bam_path: &str, ref_path: &str, event: &TruthEvent) -> Result<CheckResult> {
    let label = format_event_label(event);
    let Some(n) = event.sim_number else {
        return Ok(event_not_evaluable(
            &label,
            DEL_PLANTED,
            &format!(">={} of spike's reads carrying", MIN_PLANTED_READS),
            "no sim_del_N id",
            "the truth record's ID is not spike's sim_del_N, so which reads are spike's \
             for this event is unknown",
        ));
    };
    let expected = format!(">={} ev{:04} read carrying", MIN_PLANTED_READS, n);
    let flank = PLANTED_PROBE_FLANK as u64;
    if event.start < flank {
        return Ok(event_not_evaluable(
            &label,
            DEL_PLANTED,
            &expected,
            "start < 15",
            "fewer than 15 bases precede the deletion, so its join cannot be spelled",
        ));
    }

    // `load_truth_events` reads a DEL as `start: vcf_pos`, the 0-based first
    // deleted base, and `end` as INFO END: the haplotype is
    // `ref[..start] + ref[end..]`.
    let right_len = (PLANTED_PROBE_LEN - PLANTED_PROBE_FLANK) as u64;
    let (_, left) =
        crate::reference::fetch_window(ref_path, &event.chrom, event.start - flank, event.start)?;
    let (_, right) =
        crate::reference::fetch_window(ref_path, &event.chrom, event.end, event.end + right_len)?;
    let probes = del_junction_probes(&left, &right);
    let prefix = planted_read_prefix(n);

    let mut names: HashSet<String> = HashSet::new();
    for breakpoint in [event.start, event.end] {
        scan_region(
            bam_path,
            ref_path,
            &event.chrom,
            breakpoint.saturating_sub(DEL_PLANTED_PAD),
            breakpoint + DEL_PLANTED_PAD,
            &|flags, _mapq| !(flags.is_secondary() || flags.is_supplementary()),
            &mut |name, _align_start, _ops, seq| {
                if name.starts_with(prefix.as_bytes()) && carries_junction_probe(seq, &probes) {
                    names.insert(String::from_utf8_lossy(name).into_owned());
                }
            },
        )?;
    }

    Ok(CheckResult {
        event_label: label,
        check_name: DEL_PLANTED.to_string(),
        expected,
        observed: format!("{}", names.len()),
        pass: names.len() >= MIN_PLANTED_READS,
        advisory: false, // stamped by check_outcome
    })
}

/// Column of an `allele_counts` entry for one base, or None for anything
/// that is not A, C, G or T.
fn pileup_base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Below this the fraction is noise: at 4 reads the only fractions that exist
/// are 0, 0.25, 0.5, 0.75 and 1. Not enough data is not a pass. Every
/// allele-fraction count is of fragments (`fragment_vote`), so this is five
/// molecules whatever the variant's shape (N16).
const MIN_PILEUP_DEPTH: u32 = 5;

/// The one vote a fragment casts, from the votes of its reads.
///
/// Both mates of a pair arrive under one read name, and where they overlap
/// they read the same molecule. Counted apart, one fragment was two draws in
/// the binomial `allele_freq_result` grades against, and three pairs cleared
/// `MIN_PILEUP_DEPTH` and `MIN_ALT_READS` meant for five and three
/// molecules (N16). Mates that disagree give the fragment no vote, the rule
/// `mnv_allele_freq` applies per offset.
fn fragment_vote<T: Copy + PartialEq>(votes: &[T]) -> Option<T> {
    let (&first, rest) = votes.split_first()?;
    rest.iter().all(|&v| v == first).then_some(first)
}

/// Chance that a read of the reference shows the alt allele anyway
/// (sequencing error, mismapping). Sets the floor of alt reads that errors
/// alone can't explain.
const ALT_ERROR_RATE: f64 = 0.001;

/// Chance that a read of the alt allele shows the reference instead.
const REF_ERROR_RATE: f64 = 0.01;

/// "Outside what a correct run gives": below 0.5% on either side.
const AF_TAIL: f64 = 0.005;

/// A correct run must clear the error floor this often, or the depth is too
/// shallow to tell a spike-in from none.
const AF_POWER: f64 = 0.99;

/// Fewest alt reads that count as evidence at any depth.
const MIN_ALT_READS: u32 = 3;

/// Deepest pileup the "depth it needs" hint searches up to.
const AF_MAX_DEPTH_HINT: u32 = 1_000_000;

/// Aligned bases a read needs past an indel's repeat region, on each side,
/// before it may vote on the indel either way (N18). Near its end a read's
/// indel is written as mismatches or a clip rather than a gap, and a read that
/// stops inside the repeat cannot show an extra or missing unit at all, so
/// such a read aligns as the reference whatever it carries. Chosen on HG002
/// 35x chr20 and confirmed on held-out chr21/chr22: it takes the het-indel
/// out-of-range rate from 9.3% to 4.1% and keeps 81% of fragments.
const INDEL_FLANK: u64 = 10;

/// Reference bases fetched on each side of a small indel to find its repeat
/// region. A repeat longer than this is cut off here, which only makes the
/// region -- and the span a read needs -- shorter than it should be.
const INDEL_REF_WINDOW: u64 = 2_000;

/// What a truth record's REF/ALT pair describes, and with it how the reads
/// carrying it have to be counted.
enum SmallVariantShape {
    /// One base swapped for another: a single-position pileup of A/C/G/T.
    Substitution,
    /// Equal-length runs longer than one base (`AC` > `GT`): the whole alt run
    /// against the whole ref run, read by read.
    Mnv,
    /// `ACG` > `A`: this many reference bases after the anchor base are gone.
    Deletion(u64),
    /// `A` > `ACG`: this many bases are inserted after the anchor base.
    Insertion(u64),
}

/// The shape of a small variant, or `None` for a complex allele -- one that
/// changes length *and* rewrites the anchor base (`AC` > `GTT`, `A` > `CG`).
/// A complex allele leaves neither one CIGAR operation nor one allele run to
/// count, so it is not measured here.
fn small_variant_shape(reference: &[u8], alt: &[u8]) -> Option<SmallVariantShape> {
    let shares_anchor = matches!(
        (reference.first(), alt.first()),
        (Some(r), Some(a)) if r.eq_ignore_ascii_case(a)
    );
    match (reference.len(), alt.len()) {
        (0, _) | (_, 0) => None,
        (1, 1) => Some(SmallVariantShape::Substitution),
        (r, a) if r == a => Some(SmallVariantShape::Mnv),
        (r, 1) if shares_anchor => Some(SmallVariantShape::Deletion(r as u64 - 1)),
        (1, a) if shares_anchor => Some(SmallVariantShape::Insertion(a as u64 - 1)),
        _ => None,
    }
}

/// P(X = k) for every k in 0..=n, X ~ Bin(n, p), for 0 < p < 1. Summed in
/// log space: `(1 - p)^n` underflows long before depths validate sees.
fn binomial_pmf(n: u32, p: f64) -> Vec<f64> {
    let step = p.ln() - (1.0 - p).ln();
    let mut ln = n as f64 * (1.0 - p).ln();
    let mut out = Vec::with_capacity(n as usize + 1);
    for k in 0..=n {
        out.push(ln.exp());
        ln += ((n - k) as f64).ln() - ((k + 1) as f64).ln() + step;
    }
    out
}

/// The fewest alt reads (never under `MIN_ALT_READS`) that reference reads
/// misread at `ALT_ERROR_RATE` reach less than `AF_TAIL` of the time.
fn alt_error_floor(n: u32) -> u32 {
    let mut at_least = 1.0; // P(X >= k), from k = 0
    for (k, prob) in binomial_pmf(n, ALT_ERROR_RATE).iter().enumerate() {
        if at_least < AF_TAIL {
            return (k as u32).max(MIN_ALT_READS);
        }
        at_least -= prob;
    }
    n + 1
}

/// The alt fraction a correct run shows once error reads are counted in.
fn observed_fraction_model(expected_vaf: f64) -> f64 {
    expected_vaf * (1.0 - REF_ERROR_RATE) + (1.0 - expected_vaf) * ALT_ERROR_RATE
}

/// Whether a correct run at depth `n` clears the error floor `AF_POWER` of
/// the time, i.e. whether the reads can tell the spike-in from none.
fn af_evaluable(n: u32, p_obs: f64) -> bool {
    let floor = alt_error_floor(n) as usize;
    binomial_pmf(n, p_obs).iter().skip(floor).sum::<f64>() >= AF_POWER
}

/// Roughly the shallowest depth at which `af_evaluable` holds, for the
/// "too shallow" message; None past `AF_MAX_DEPTH_HINT`.
fn af_depth_needed(from: u32, p_obs: f64) -> Option<u32> {
    let mut hi = from.max(MIN_PILEUP_DEPTH);
    while !af_evaluable(hi, p_obs) {
        if hi >= AF_MAX_DEPTH_HINT {
            return None;
        }
        hi = hi.saturating_mul(2).min(AF_MAX_DEPTH_HINT);
    }
    let mut lo = from.max(MIN_PILEUP_DEPTH);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if af_evaluable(mid, p_obs) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    Some(hi)
}

/// One allele-fraction verdict from an alt count and a total. The same
/// grading applies whichever counting rule produced the two numbers, so
/// every shape of small variant is graded here.
///
/// It passes only when the reads could tell the spike-in from none (the
/// alt count clears what errors alone reach, at a depth where a correct run
/// does that 99% of the time) and the count sits inside the central 99% of
/// what a correct run gives. A fixed tolerance can do neither: +-0.15
/// passed every SIM_VAF under 0.15 on zero alt reads.
fn allele_freq_result(event: &TruthEvent, alt: u32, total: u32) -> CheckResult {
    let label = format_event_label(event);
    let expected = format!("{:.2}", event.expected_vaf);

    if total == 0 {
        return CheckResult {
            event_label: label,
            check_name: "allele_freq".to_string(),
            expected,
            observed: "no coverage".to_string(),
            pass: false,
            advisory: false,
        };
    }
    if total < MIN_PILEUP_DEPTH {
        return event_not_evaluable(
            &label,
            "allele_freq",
            &expected,
            &format!("low depth ({})", total),
            "fewer reads than the check needs to measure a fraction",
        );
    }

    let p_obs = observed_fraction_model(event.expected_vaf);
    let observed_vaf = alt as f64 / total as f64;
    let floor = alt_error_floor(total);
    let pmf = binomial_pmf(total, p_obs);
    if pmf.iter().skip(floor as usize).sum::<f64>() < AF_POWER {
        let needed = match af_depth_needed(total, p_obs) {
            Some(n) => format!("about {} reads", n),
            None => format!("more than {} reads", AF_MAX_DEPTH_HINT),
        };
        return event_not_evaluable(
            &label,
            "allele_freq",
            &expected,
            &format!(
                "too shallow ({:.2} at {} reads; needs {})",
                observed_vaf, total, needed
            ),
            "a correct spike-in would too often show no more alt reads than errors do",
        );
    }

    let x = alt.min(total) as usize;
    let at_most: f64 = pmf[..=x].iter().sum();
    let at_least: f64 = pmf[x..].iter().sum();
    // A hom truth has no "too many alt reads".
    let in_range = at_most >= AF_TAIL && (event.expected_vaf >= 1.0 || at_least >= AF_TAIL);
    CheckResult {
        event_label: label,
        check_name: "allele_freq".to_string(),
        expected,
        observed: format!("{:.2}", observed_vaf),
        pass: alt >= floor && in_range,
        advisory: false,
    }
}

/// Check allele frequency for small variants via pileup.
fn check_allele_freq(
    bam_path: &str,
    ref_path: &str,
    event: &TruthEvent,
    nearby: &NearbyRecords,
    min_mapq: u8,
) -> Result<CheckResult> {
    let label = format_event_label(event);
    let expected = format!("{:.2}", event.expected_vaf);

    let ref_allele = event.ref_allele.as_ref().unwrap();
    let alt_allele = event.alt_allele.as_ref().unwrap();

    // Each shape of small variant leaves its own mark in an alignment, so each
    // is counted its own way. A shape none of them fits is not measured, and
    // an unmeasured allele fraction is a failed check, not a pass (N9): a
    // result row is still pushed, so the event counts as covered and the row
    // says out loud that nothing was measured.
    let Some(shape) = small_variant_shape(ref_allele, alt_allele) else {
        return Ok(event_not_evaluable(
            &label,
            "allele_freq",
            &expected,
            "N/A (complex allele)",
            "REF and ALT are neither a substitution, an MNV, nor a small indel \
             sharing an anchor base",
        ));
    };
    match shape {
        SmallVariantShape::Deletion(len) => {
            return indel_allele_freq(bam_path, ref_path, event, nearby, min_mapq, Kind::Deletion, len)
        }
        SmallVariantShape::Insertion(len) => {
            return indel_allele_freq(bam_path, ref_path, event, nearby, min_mapq, Kind::Insertion, len)
        }
        SmallVariantShape::Mnv => return mnv_allele_freq(bam_path, ref_path, event, min_mapq),
        SmallVariantShape::Substitution => {}
    }

    // A non-ACGT alt has no column in the pileup. Checked before the BAM is
    // read: it is a property of the truth record, not of the data.
    let alt_base = alt_allele[0].to_ascii_uppercase();
    let Some(alt_idx) = pileup_base_index(alt_base) else {
        return Ok(event_not_evaluable(
            &label,
            "allele_freq",
            &expected,
            "unknown alt base",
            "the alt allele is not one of A, C, G, T",
        ));
    };

    // Pileup at the variant position, one vote per fragment: the per-read
    // map rather than the per-base counts, which see an overlapping pair
    // twice (N16).
    let mut dummy_allele_counts: HashMap<u64, [u32; 4]> = HashMap::new();
    let mut read_alleles: HashMap<String, Vec<(u64, u8)>> = HashMap::new();

    pileup_region(
        bam_path,
        ref_path,
        &event.chrom,
        event.start,
        event.start + 1,
        min_mapq,
        &mut dummy_allele_counts,
        &mut read_alleles,
    )?;

    let mut counts = [0u32; 4];
    for bases in read_alleles.values() {
        let at_site: Vec<u8> = bases
            .iter()
            .filter(|&&(rp, _)| rp == event.start)
            .map(|&(_, base)| base)
            .collect();
        if let Some(idx) = fragment_vote(&at_site).and_then(pileup_base_index) {
            counts[idx] += 1;
        }
    }
    Ok(allele_freq_result(
        event,
        counts[alt_idx],
        counts.iter().sum(),
    ))
}

/// Allele fraction of a small indel, counted off the alignments' CIGARs.
///
/// An indel is not a column in a pileup, and where an aligner writes its gap
/// -- or whether it writes one at all -- depends on the read. So each read's
/// own bases over the site are compared with the reference with and without
/// the site's edit, and with every other truth record nearby tried in too,
/// and the read counts for the nearer kind (N15). This is the MNV rule (N10)
/// with a nearest match instead of an exact one.
fn indel_allele_freq(
    bam_path: &str,
    ref_path: &str,
    event: &TruthEvent,
    nearby: &NearbyRecords,
    min_mapq: u8,
    kind: Kind,
    indel_len: u64,
) -> Result<CheckResult> {
    let (carries, spans) =
        count_indel_reads(bam_path, ref_path, event, nearby, kind, indel_len, min_mapq)?;
    log::debug!(
        "allele_freq {}: {} carry, {} span",
        format_event_label(event),
        carries,
        spans
    );
    Ok(allele_freq_result(event, carries, carries + spans))
}

/// Fragments whose reads' bases over the site are nearest haplotypes holding
/// the site, and fragments whose reads' bases are nearest haplotypes without
/// it -- one vote per fragment (`fragment_vote`), not per read.
fn count_indel_reads(
    bam_path: &str,
    ref_path: &str,
    event: &TruthEvent,
    nearby: &NearbyRecords,
    kind: Kind,
    indel_len: u64,
    min_mapq: u8,
) -> Result<(u32, u32)> {
    let ref_len = event.ref_allele.as_ref().map_or(1, |r| r.len() as u64);
    let inserted: Vec<u8> = match kind {
        Kind::Insertion => event
            .alt_allele
            .as_ref()
            .map_or(Vec::new(), |a| a[1..].to_ascii_uppercase()),
        _ => Vec::new(),
    };
    let (window_start, window) = crate::reference::fetch_window(
        ref_path,
        &event.chrom,
        event.start.saturating_sub(INDEL_REF_WINDOW),
        event.start + ref_len + INDEL_REF_WINDOW,
    )?;
    let (region_start, region_end) =
        indel_repeat_region(&window, window_start, event.start, kind, indel_len, &inserted);
    // The site: the repeat region, the base on each side, and INDEL_FLANK
    // more (N18). A read has to cover all of it to vote, and its bases over
    // it are what it votes with.
    let window_end = window_start + window.len() as u64;
    let first = region_start.saturating_sub(1 + INDEL_FLANK).max(window_start);
    let last = (region_end + INDEL_FLANK).min(window_end.saturating_sub(1));
    let site: Edit = (
        event.start,
        event.ref_allele.as_deref().unwrap_or_default(),
        event.alt_allele.as_deref().unwrap_or_default(),
    );
    // The other truth records in the window: the reads may carry any of them
    // with the site, as one local change the truth set wrote in parts (N19).
    let mut others = nearby.inside(event, first, last);
    if others.len() > MAX_NEARBY_EDITS {
        log::debug!(
            "allele_freq {}: {} truth records nearby, more than {}; comparing with the site alone",
            format_event_label(event),
            others.len(),
            MAX_NEARBY_EDITS
        );
        others.clear();
    }
    let Some(haplotypes) = site_haplotypes(&window, window_start, first, last, site, &others) else {
        return Ok((0, 0));
    };
    let mut votes: HashMap<Vec<u8>, Vec<IndelVote>> = HashMap::new();
    // Reads often share their bases over the site; each distinct run is
    // compared with the haplotypes once.
    let mut vote_of: HashMap<Vec<u8>, Option<IndelVote>> = HashMap::new();

    // Every read that votes either way covers the anchor base and the base
    // just past the REF allele, so that span is also the query window.
    for_each_alignment(
        bam_path,
        ref_path,
        &event.chrom,
        event.start,
        event.start + ref_len + 1,
        min_mapq,
        &mut |name, align_start, ops, seq| {
            let Some(bases) = read_bases_over(ops, seq, align_start, first, last) else {
                return;
            };
            let vote = *vote_of
                .entry(bases.to_vec())
                .or_insert_with(|| haplotype_vote(bases, &haplotypes));
            if let Some(vote) = vote {
                votes.entry(name.to_vec()).or_default().push(vote);
            }
        },
    )?;

    let (mut carries, mut spans) = (0u32, 0u32);
    for fragment in votes.values() {
        match fragment_vote(fragment) {
            Some(IndelVote::Carries) => carries += 1,
            Some(IndelVote::Spans) => spans += 1,
            None => {}
        }
    }
    Ok((carries, spans))
}

/// The truth set's small variants on each chromosome, sorted by start, so an
/// indel can find the records near it (N15).
#[derive(Default)]
struct NearbyRecords<'a> {
    by_chrom: HashMap<&'a str, Vec<&'a TruthEvent>>,
}

impl<'a> NearbyRecords<'a> {
    fn new(events: &'a [TruthEvent]) -> Self {
        let mut by_chrom: HashMap<&str, Vec<&TruthEvent>> = HashMap::new();
        for e in events.iter().filter(|e| e.ref_allele.is_some() && e.alt_allele.is_some()) {
            by_chrom.entry(e.chrom.as_str()).or_default().push(e);
        }
        for records in by_chrom.values_mut() {
            records.sort_by_key(|e| e.start);
        }
        NearbyRecords { by_chrom }
    }

    /// The edits of the records other than `event` whose REF lies inside
    /// 0-based `[first, last]`, one per record.
    fn inside(&self, event: &TruthEvent, first: u64, last: u64) -> Vec<Edit<'a>> {
        let Some(records) = self.by_chrom.get(event.chrom.as_str()) else {
            return Vec::new();
        };
        let same = |r: &[u8], a: &[u8]| {
            event.ref_allele.as_deref().is_some_and(|x| x.eq_ignore_ascii_case(r))
                && event.alt_allele.as_deref().is_some_and(|x| x.eq_ignore_ascii_case(a))
        };
        let from = records.partition_point(|e| e.start < first);
        records[from..]
            .iter()
            .take_while(|e| e.start <= last)
            .filter_map(|e| Some((e.start, e.ref_allele.as_deref()?, e.alt_allele.as_deref()?)))
            .filter(|&(p, r, a)| p + (r.len().max(1) as u64) - 1 <= last && !(p == event.start && same(r, a)))
            .collect()
    }
}

/// One truth record's change to the reference: its 0-based start, REF and
/// ALT.
type Edit<'a> = (u64, &'a [u8], &'a [u8]);

/// With more truth records than this near an indel, it is compared with its
/// own record alone: the combinations double with each one.
const MAX_NEARBY_EDITS: usize = 10;

/// How one read votes on a small indel at a position.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum IndelVote {
    /// The haplotypes nearest its bases over the site all hold the site.
    Carries,
    /// The haplotypes nearest its bases over the site all lack it.
    Spans,
}

/// The reference over 0-based `[first, last]` with every combination of the
/// site's edit and the `others` applied whose REF spans do not overlap, each
/// marked with whether it holds the site; `None` if the site does not fit.
/// No genotype is used: which records a read carries, its bases say.
fn site_haplotypes(
    window: &[u8],
    window_start: u64,
    first: u64,
    last: u64,
    site: Edit,
    others: &[Edit],
) -> Option<Vec<(bool, Vec<u8>)>> {
    let mut haplotypes = Vec::new();
    for mask in 0u32..1 << others.len() {
        let chosen = others.iter().enumerate().filter(|(i, _)| mask >> i & 1 == 1);
        let mut edits: Vec<Edit> = chosen.map(|(_, &e)| e).collect();
        if let Some(h) = apply_edits(window, window_start, first, last, &mut edits) {
            haplotypes.push((false, h));
        }
        edits.push(site);
        if let Some(h) = apply_edits(window, window_start, first, last, &mut edits) {
            haplotypes.push((true, h));
        }
    }
    haplotypes.iter().any(|&(holds, _)| holds).then_some(haplotypes)
}

/// The reference over 0-based `[first, last]` with `edits` applied, upper
/// case; `None` if two of them overlap or one does not lie inside.
fn apply_edits(
    window: &[u8],
    window_start: u64,
    first: u64,
    last: u64,
    edits: &mut [Edit],
) -> Option<Vec<u8>> {
    let at = |p: u64| usize::try_from(p.checked_sub(window_start)?).ok();
    let (f, l) = (at(first)?, at(last)?);
    if l >= window.len() {
        return None;
    }
    edits.sort_by_key(|e| e.0);
    let mut h = Vec::new();
    let mut from = f;
    for &mut (pos, reference, alt) in edits {
        let p = at(pos)?;
        let past = p + reference.len();
        if p < from || past > l + 1 {
            return None;
        }
        h.extend_from_slice(&window[from..p]);
        h.extend_from_slice(alt);
        from = past;
    }
    h.extend_from_slice(&window[from..=l]);
    h.make_ascii_uppercase();
    Some(h)
}

/// A read's bases from the one aligned to reference position `first` to the
/// one aligned to `last`, inclusive, inserted bases between them included; or
/// `None` unless both ends sit on an aligned (`M`, `=`, `X`) base. A read that
/// stops short of either end, or has a gap or clip there, cannot show the
/// site whole (N18).
fn read_bases_over<'a>(
    ops: &[noodles::sam::alignment::record::cigar::Op],
    seq: &'a [u8],
    align_start: u64,
    first: u64,
    last: u64,
) -> Option<&'a [u8]> {
    let (mut ref_pos, mut q) = (align_start, 0usize);
    let (mut q_first, mut q_last) = (None, None);
    for op in ops {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let span = ref_pos..ref_pos + len as u64;
                if span.contains(&first) {
                    q_first = Some(q + (first - ref_pos) as usize);
                }
                if span.contains(&last) {
                    q_last = Some(q + (last - ref_pos) as usize);
                }
                ref_pos += len as u64;
                q += len;
            }
            Kind::Insertion | Kind::SoftClip => q += len,
            Kind::Deletion | Kind::Skip => ref_pos += len as u64,
            Kind::HardClip | Kind::Pad => {}
        }
    }
    seq.get(q_first?..=q_last?)
}

/// Levenshtein distance, ignoring case.
fn levenshtein(a: &[u8], b: &[u8]) -> usize {
    let mut row: Vec<usize> = (0..=b.len()).collect();
    for (i, x) in a.iter().enumerate() {
        let mut diag = row[0];
        row[0] = i + 1;
        for (j, y) in b.iter().enumerate() {
            let next = if x.eq_ignore_ascii_case(y) { diag } else { diag + 1 };
            diag = row[j + 1];
            row[j + 1] = next.min(row[j] + 1).min(row[j + 1] + 1);
        }
    }
    row[b.len()]
}

/// A read's vote from its bases over the site: `Carries` if every haplotype
/// nearest them holds the site, `Spans` if none does, none otherwise.
fn haplotype_vote(bases: &[u8], haplotypes: &[(bool, Vec<u8>)]) -> Option<IndelVote> {
    let distances: Vec<(usize, bool)> = haplotypes
        .iter()
        .map(|(holds, h)| (levenshtein(bases, h), *holds))
        .collect();
    let nearest = distances.iter().map(|&(d, _)| d).min()?;
    let mut holds = distances.iter().filter(|&&(d, _)| d == nearest).map(|&(_, h)| h);
    let first = holds.next()?;
    if holds.all(|h| h == first) {
        Some(if first { IndelVote::Carries } else { IndelVote::Spans })
    } else {
        None
    }
}

/// The stretch of reference a small indel can slide along, 0-based
/// `[start, end)`: the deleted bases, or the empty junction an insertion
/// goes into, widened both ways for as long as the reference repeats the
/// deleted or inserted unit. `pos` is the anchor base; `window` holds the
/// reference from `window_start` on, and the region stops at its edges.
fn indel_repeat_region(
    window: &[u8],
    window_start: u64,
    pos: u64,
    kind: Kind,
    len: u64,
    inserted: &[u8],
) -> (u64, u64) {
    let base = |p: u64| -> Option<u8> {
        let i = usize::try_from(p.checked_sub(window_start)?).ok()?;
        window.get(i).map(u8::to_ascii_uppercase)
    };
    let same = |a: Option<u8>, b: Option<u8>| matches!((a, b), (Some(x), Some(y)) if x == y);
    match kind {
        Kind::Deletion => {
            let (mut s, mut e) = (pos + 1, pos + 1 + len);
            while s > 0 && same(base(s - 1), base(s - 1 + len)) {
                s -= 1;
            }
            while same(base(e), base(e - len)) {
                e += 1;
            }
            (s, e)
        }
        _ if inserted.is_empty() => (pos + 1, pos + 1),
        _ => {
            let unit = |k: usize| Some(inserted[k % inserted.len()].to_ascii_uppercase());
            let (mut s, mut e) = (pos + 1, pos + 1);
            let mut k = 0;
            while same(base(e), unit(k)) {
                e += 1;
                k += 1;
            }
            k = 0;
            while s > 0 && same(base(s - 1), unit(inserted.len() - 1 - k % inserted.len())) {
                s -= 1;
                k += 1;
            }
            (s, e)
        }
    }
}

/// Allele fraction of an MNV: reads whose bases are the whole alt run against
/// reads whose bases are the whole ref run.
///
/// The two runs are the same length, so this is a pileup like a substitution's
/// -- but a fraction per base answers a different question at every offset, and
/// a read that carries only one of the two substitutions is not this variant.
/// The alleles are therefore read jointly, one read at a time, which is what
/// `pileup_region`'s per-read map is for. A read matching neither run whole is
/// evidence for neither and enters neither count.
fn mnv_allele_freq(
    bam_path: &str,
    ref_path: &str,
    event: &TruthEvent,
    min_mapq: u8,
) -> Result<CheckResult> {
    let ref_allele = event.ref_allele.as_ref().unwrap();
    let alt_allele = event.alt_allele.as_ref().unwrap();

    // A run holding a base outside A/C/G/T has no column in the pileup, so it
    // could never be matched. Checked before the BAM is read: it is a property
    // of the truth record, not of the data.
    if ref_allele
        .iter()
        .chain(alt_allele.iter())
        .any(|b| pileup_base_index(*b).is_none())
    {
        return Ok(event_not_evaluable(
            &format_event_label(event),
            "allele_freq",
            &format!("{:.2}", event.expected_vaf),
            "unknown allele base",
            "an allele holds a base that is not one of A, C, G, T",
        ));
    }

    let mut allele_counts: HashMap<u64, [u32; 4]> = HashMap::new();
    let mut read_alleles: HashMap<String, Vec<(u64, u8)>> = HashMap::new();

    pileup_region(
        bam_path,
        ref_path,
        &event.chrom,
        event.start,
        event.start + ref_allele.len() as u64,
        min_mapq,
        &mut allele_counts,
        &mut read_alleles,
    )?;

    let (mut alt_reads, mut ref_reads) = (0u32, 0u32);
    for bases in read_alleles.values() {
        // Both mates of a pair arrive under one name; where they overlap and
        // disagree, the offset is left unreadable so the read matches neither
        // run.
        let mut seen: HashMap<u64, u8> = HashMap::new();
        for &(rp, base) in bases {
            seen.entry(rp)
                .and_modify(|b| {
                    if *b != base {
                        *b = b'.';
                    }
                })
                .or_insert(base);
        }
        let run: Option<Vec<u8>> = (0..ref_allele.len() as u64)
            .map(|i| seen.get(&(event.start + i)).copied())
            .collect();
        let Some(run) = run else { continue };
        if run.eq_ignore_ascii_case(alt_allele) {
            alt_reads += 1;
        } else if run.eq_ignore_ascii_case(ref_allele) {
            ref_reads += 1;
        }
    }

    Ok(allele_freq_result(event, alt_reads, alt_reads + ref_reads))
}

// ---------------------------------------------------------------------------
// Global checks
// ---------------------------------------------------------------------------

/// Records to sample for the global checks, over all event regions together.
const GLOBAL_SAMPLE_MAX: u64 = 200_000;

/// The smallest sample a missing duplicate flag says anything about. At a 1%
/// duplicate rate a 41-record window holds no duplicate about two times in
/// three, while 1 000 records hold none with probability 4e-5, so below this
/// "no record carries the flag" is a small sample, not an unmarked file.
const DUP_RATE_MIN_SAMPLE: u64 = 1_000;

/// What the three global checks are computed from: primary, mapped records
/// sampled from the truth events' own regions.
#[derive(Default)]
struct GlobalSample {
    /// Records sampled.
    total: u64,
    /// ...of which carry the duplicate flag.
    dups: u64,
    /// Sum of their mapping qualities.
    mapq_sum: u64,
    /// Template lengths of the properly-paired, non-duplicate records.
    insert_sizes: Vec<f64>,
    /// Whether the alignment file's header records a duplicate-marking step.
    duplicates_marked: bool,
}

impl GlobalSample {
    /// Add one record. Unmapped, secondary and supplementary records are the
    /// caller's to skip.
    fn add(
        &mut self,
        flags: noodles::sam::alignment::record::Flags,
        mapq: u8,
        template_length: i32,
    ) {
        self.total += 1;
        self.mapq_sum += mapq as u64;
        if flags.is_duplicate() {
            self.dups += 1;
        }
        // The same filter `bam_stats` uses: a duplicate's or a QC-failed
        // record's template length is not an insert-size observation.
        if !flags.is_duplicate()
            && !flags.is_qc_fail()
            && flags.is_properly_segmented()
            && !flags.is_mate_unmapped()
            && template_length > 0
        {
            self.insert_sizes.push(template_length as f64);
        }
    }

    /// Mean and standard deviation of the sampled insert sizes, or `None`
    /// when nothing properly paired was sampled.
    fn insert_stats(&self) -> Option<(f64, f64)> {
        if self.insert_sizes.is_empty() {
            return None;
        }
        let mean = self.insert_sizes.iter().sum::<f64>() / self.insert_sizes.len() as f64;
        let variance = self
            .insert_sizes
            .iter()
            .map(|x| (x - mean).powi(2))
            .sum::<f64>()
            / self.insert_sizes.len() as f64;
        Some((mean, variance.sqrt()))
    }
}

/// Sample up to `budget` records for the global checks from every truth
/// event's own window (event +/- `flank_bp`).
///
/// The head of a file is not a sample of it: on whole-genome HG002 the first
/// 100k records are chr1's telomere, mean MAPQ 10.0, which fails a check the
/// rest of the file passes (L15). The event windows are both representative of
/// the reads `validate` is judging and cheap to read, because they are indexed
/// queries like every other check here rather than a walk from the top.
///
/// Every event gets an equal share of the budget and no event is skipped. The
/// sample is still the *head* of each window, in coordinate order, so an event
/// longer than its share covers (above about 1 Mb at 35x with the default
/// budget and one event) is represented by its first records rather than by a
/// spread over it. The size of the sample, and of each region's contribution,
/// is logged, because a check computed over 41 records reads exactly like one
/// computed over 200 000.
fn sample_event_regions(
    bam_path: &str,
    ref_path: &str,
    events: &[TruthEvent],
    flank_bp: u64,
    budget: u64,
) -> Result<GlobalSample> {
    if events.is_empty() {
        bail!("truth VCF holds no events, so there is no region to sample");
    }

    // An equal share each, with no floor. A floor above the fair share
    // silently truncates the event list: at 1 000 records each, the first 200
    // events of a 250-event panel spend the whole 200 000 and events 201..250
    // are never read, with nothing in the log or the output to say so.
    let per_region = (budget / events.len() as u64).max(1);

    let mut sample = GlobalSample::default();
    let mut visited = 0usize;
    let mut failed = 0usize;
    for event in events {
        if sample.total >= budget {
            break;
        }
        let start = event.start.saturating_sub(flank_bp);
        let end = event.end + flank_bp;
        visited += 1;
        match sample_region(
            bam_path,
            ref_path,
            &event.chrom,
            start,
            end,
            per_region,
            &mut sample,
        ) {
            Ok(taken) => log::info!(
                "{} sampled {} records from {}:{}-{}",
                GLOBAL_LABEL,
                taken,
                event.chrom,
                start + 1,
                end
            ),
            // One region that cannot be queried -- a truth event on a contig
            // the alignment file does not have -- must not fail all three
            // global checks, the way the per-event checks isolate an error per
            // check through `check_outcome`.
            Err(e) => {
                failed += 1;
                log::warn!(
                    "{} skipping {}:{}-{}: {:#}",
                    GLOBAL_LABEL,
                    event.chrom,
                    start + 1,
                    end,
                    e
                );
            }
        }
    }

    if visited < events.len() {
        log::warn!(
            "{} the {}-record budget ran out after {} of {} event regions; the rest went unsampled",
            GLOBAL_LABEL,
            budget,
            visited,
            events.len()
        );
    }

    if sample.total == 0 && failed > 0 {
        bail!(
            "no records could be sampled: all {} truth event regions failed to query",
            failed
        );
    }

    log::info!(
        "{} sample: {} records over {} of {} event regions (up to {} each)",
        GLOBAL_LABEL,
        sample.total,
        visited - failed,
        events.len(),
        per_region
    );

    Ok(sample)
}

/// Add up to `max_records` of one region's primary alignments to `sample`,
/// returning how many were taken. The records are the window's first in
/// coordinate order, not a spread over it.
fn sample_region(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    start: u64,
    end: u64,
    max_records: u64,
    sample: &mut GlobalSample,
) -> Result<u64> {
    let start_pos = crate::extract::safe_noodles_position(start + 1);
    let end_pos = crate::extract::safe_noodles_position(end);
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
    let mut taken: u64 = 0;

    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let (mut reader, header) =
            crate::extract::open_cram_reader_for_region(bam_path, &repository, &region)
                .context("failed to open CRAM for the global sample")?;
        sample.duplicates_marked |= header_marks_duplicates(&header);
        let query = reader.query(&header, &region)?;
        // `query` has already rejected an unknown contig, so this is `Some`.
        let queried_reference_sequence_id =
            header.reference_sequences().get_index_of(chrom.as_bytes());

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            // A multi-contig container is decoded whole and `Query` filters on
            // coordinates alone, so another contig's reads would enter the
            // sample (L2, N4). The BAM arm needs no such guard.
            if !record_is_on_queried_reference(&buf, queried_reference_sequence_id) {
                continue;
            }
            let flags = buf.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            let mq: u8 = buf.mapping_quality().map(u8::from).unwrap_or(0);
            sample.add(flags, mq, buf.template_length());
            taken += 1;
            if taken >= max_records {
                break;
            }
        }
    } else {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .context("failed to open BAM for the global sample")?;
        let header = reader.read_header()?;
        sample.duplicates_marked |= header_marks_duplicates(&header);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let record = rec_result?;
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            let mq: u8 = record.mapping_quality().map(u8::from).unwrap_or(0);
            sample.add(flags, mq, record.template_length());
            taken += 1;
            if taken >= max_records {
                break;
            }
        }
    }

    Ok(taken)
}

/// Whether the header records a duplicate-marking step.
///
/// A `@PG` record from a duplicate marker (`samtools markdup`, Picard
/// `MarkDuplicates`, `sambamba markdup`, biobambam's `bammarkduplicates`, or a
/// UMI pipeline's `dedup`) means a record without the duplicate flag is a
/// record that step decided is not a duplicate -- so a sample carrying no
/// duplicate flag is a measured 0%, not a file nothing ever marked. The tool's
/// name is usually only in the command line: HG002's marker is `ID:samtools.4`
/// `PN:samtools`, `CL:... samtools markdup -@ 4 - out.bam`.
fn header_marks_duplicates(header: &noodles::sam::Header) -> bool {
    header.programs().as_ref().iter().any(|(id, program)| {
        let mut fields = vec![String::from_utf8_lossy(id.as_ref()).into_owned()];
        for value in program.other_fields().values() {
            fields.push(String::from_utf8_lossy(value.as_ref()).into_owned());
        }
        fields
            .iter()
            .any(|f| f.split_whitespace().any(names_a_duplicate_marker))
    })
}

/// Whether one `@PG` token names a duplicate marker.
///
/// A token, not a substring of the whole record: HG002's `samtools markdup`
/// writes to a path with `dedup` in its name, and a file name is not a
/// program. A false positive here reports a fabricated 0% duplicate rate on a
/// file nothing marked -- the silent pass this check exists to avoid -- so the
/// match is deliberately narrow.
fn names_a_duplicate_marker(token: &str) -> bool {
    let name = token
        .rsplit('/')
        .next()
        .unwrap_or(token)
        .to_ascii_lowercase();
    name == "markdup"
        || name == "dedup"
        || name.starts_with("markduplicates")
        || name.starts_with("bammarkduplicates")
}

/// True if a read aligned at `align_start` carries at least `min_len` bases
/// that do not fit the reference within `pad` of `pos`.
///
/// An aligner represents an insertion two ways depending on its size: an `I`
/// operation when the inserted sequence is short enough to fit inside a read
/// that still anchors on both sides, and a soft clip at the insertion point
/// when it is not. Both are counted, at the reference position where the
/// alignment leaves the reference -- for a leading clip that is the alignment
/// start, for a trailing clip the alignment end. A soft clip counts only once
/// `min_len` reaches [`INS_MAX_EVIDENCE_LEN`], which is exactly the point at
/// which a read can no longer anchor both sides of the insertion: below it the
/// aligner writes an `I` operation, and a short clip is background anywhere.
fn cigar_shows_insertion_near(
    ops: impl Iterator<Item = std::io::Result<noodles::sam::alignment::record::cigar::Op>>,
    align_start: u64,
    pos: u64,
    pad: u64,
    min_len: u64,
) -> bool {
    let mut ref_pos = align_start;
    for op in ops.flatten() {
        match op.kind() {
            Kind::Insertion => {
                if op.len() as u64 >= min_len && ref_pos.abs_diff(pos) <= pad {
                    return true;
                }
            }
            Kind::SoftClip => {
                // A clip only counts once the insertion is long enough that a
                // read cannot anchor both sides of it -- which is exactly the
                // `INS_MAX_EVIDENCE_LEN` cap, since `min_len` reaches it only
                // for an insertion at least that long. Below it the aligner
                // writes an `I` operation, and a short clip is background: at
                // a 3 bp threshold five positions on the HG002 chr20 slice
                // where nothing was planted give 1, 3, 0, 0 and 1 clipped
                // reads, and 2 is a pass.
                if min_len >= INS_MAX_EVIDENCE_LEN
                    && op.len() as u64 >= min_len
                    && ref_pos.abs_diff(pos) <= pad
                {
                    return true;
                }
            }
            Kind::Match
            | Kind::Deletion
            | Kind::Skip
            | Kind::SequenceMatch
            | Kind::SequenceMismatch => {
                ref_pos += op.len() as u64;
            }
            // HardClip and Pad consume neither reference nor stored sequence.
            _ => {}
        }
    }
    false
}

/// A global check its sample cannot answer: a failed result, never a silent
/// pass (M10).
fn not_evaluable(check_name: &str, expected: &str, observed: &str, why: &str) -> CheckResult {
    event_not_evaluable(GLOBAL_LABEL, check_name, expected, observed, why)
}

/// A check on one event that its data cannot answer: a failed result, never a
/// silent pass.
///
/// The rule is the same one [`not_evaluable`] applies to the three global
/// checks and `check_event` applies to an event type no check covers -- a
/// verdict nothing was measured for is a failure, not a pass. A result row is
/// pushed either way, so `check_event`'s "a check applies" gate counts the
/// event as covered; only `pass` says whether anything was actually measured.
fn event_not_evaluable(
    label: &str,
    check_name: &str,
    expected: &str,
    observed: &str,
    why: &str,
) -> CheckResult {
    log::warn!("{} check is not evaluable for {}: {}", check_name, label, why);
    CheckResult {
        event_label: label.to_string(),
        check_name: check_name.to_string(),
        expected: expected.to_string(),
        observed: observed.to_string(),
        pass: false,
        advisory: false,
    }
}

/// Check the insert size distribution of the sampled records.
fn check_insert_size(sample: &GlobalSample) -> CheckResult {
    let expected = "mean 50-1000, sd 5-300";
    let Some((mean, stddev)) = sample.insert_stats() else {
        return not_evaluable(
            "insert_size",
            expected,
            "no pairs",
            "no properly-paired record in the event regions",
        );
    };

    // Reasonable Illumina ranges.
    let pass = (50.0..=1000.0).contains(&mean) && (5.0..=300.0).contains(&stddev);

    CheckResult {
        event_label: GLOBAL_LABEL.to_string(),
        check_name: "insert_size".to_string(),
        expected: expected.to_string(),
        observed: format!("{:.0}+/-{:.0}", mean, stddev),
        pass,
        advisory: false,
    }
}

/// Check the duplicate rate of the sampled records.
fn check_dup_rate(sample: &GlobalSample) -> CheckResult {
    if sample.total == 0 {
        return not_evaluable("dup_rate", "<50%", "no reads", "no record sampled");
    }
    if sample.dups == 0 && !sample.duplicates_marked {
        // No duplicate flag in a file whose header names no duplicate marker.
        // Either nothing marked them -- in which case 0% would be a default,
        // not a measurement -- or the window is too small to hold one. Say
        // which, rather than telling the owner of a marked BAM to mark it.
        if sample.total < DUP_RATE_MIN_SAMPLE {
            return not_evaluable(
                "dup_rate",
                "<50%",
                "too few reads",
                &format!(
                    "{} sampled records carry no duplicate flag, too few to tell an unmarked \
                     file from a window that holds no duplicate",
                    sample.total
                ),
            );
        }
        return not_evaluable(
            "dup_rate",
            "<50%",
            "no dup flags",
            "no sampled record carries the duplicate flag and no @PG record marked duplicates \
             -- mark duplicates to evaluate it",
        );
    }

    let rate = sample.dups as f64 / sample.total as f64;
    let pass = rate < 0.50;

    CheckResult {
        event_label: GLOBAL_LABEL.to_string(),
        check_name: "dup_rate".to_string(),
        expected: "<50%".to_string(),
        observed: format!("{:.1}%", rate * 100.0),
        pass,
        advisory: false,
    }
}

/// Check the mean mapping quality of the sampled records.
fn check_mapq(sample: &GlobalSample) -> CheckResult {
    if sample.total == 0 {
        return not_evaluable("mean_mapq", ">20", "no reads", "no record sampled");
    }

    let mean_mapq = sample.mapq_sum as f64 / sample.total as f64;
    let pass = mean_mapq > 20.0;

    CheckResult {
        event_label: GLOBAL_LABEL.to_string(),
        check_name: "mean_mapq".to_string(),
        expected: ">20".to_string(),
        observed: format!("{:.1}", mean_mapq),
        pass,
        advisory: false,
    }
}

// ---------------------------------------------------------------------------
// BAM reading helpers
// ---------------------------------------------------------------------------

/// Mean read depth (aligned bases per base) in a region.
fn count_depth_in_region(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    start: u64,
    end: u64,
    min_mapq: u8,
) -> Result<f64> {
    if start >= end {
        return Ok(0.0);
    }

    let mut spans: Vec<(u64, u64)> = Vec::new();

    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let start_pos = crate::extract::safe_noodles_position(start + 1);
        let end_pos = crate::extract::safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let (mut reader, header) =
            crate::extract::open_cram_reader_for_region(bam_path, &repository, &region)
                .context("failed to open CRAM for depth counting")?;
        let query = reader.query(&header, &region)?;
        // `query` has already rejected an unknown contig, so this is `Some`.
        let queried_reference_sequence_id =
            header.reference_sequences().get_index_of(chrom.as_bytes());

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            // A container holding several contigs is decoded whole and `Query`
            // filters on coordinates alone, so another contig's reads would be
            // counted as this region's (L2, N4). The BAM arm below needs no
            // such guard: noodles-bam's `Query` compares the reference id.
            if !record_is_on_queried_reference(&buf, queried_reference_sequence_id) {
                continue;
            }
            let flags = buf.flags();
            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_duplicate()
                || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = buf.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            if let Some(p) = buf.alignment_start() {
                let s = usize::from(p) as u64 - 1;
                spans.push((s, s + ref_span(CigarTrait::iter(&buf.cigar()))));
            }
        }
    } else {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .context("failed to open BAM for depth counting")?;
        let header = reader.read_header()?;

        let start_pos = crate::extract::safe_noodles_position(start + 1);
        let end_pos = crate::extract::safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let record = rec_result?;
            let flags = record.flags();
            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_duplicate()
                || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = record.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            if let Some(Ok(p)) = record.alignment_start() {
                let s = usize::from(p) as u64 - 1;
                spans.push((s, s + ref_span(record.cigar().iter())));
            }
        }
    }

    Ok(mean_depth(&spans, start, end))
}

/// Mean depth over [start, end) from read reference spans [s, e).
fn mean_depth(spans: &[(u64, u64)], start: u64, end: u64) -> f64 {
    if start >= end {
        return 0.0;
    }
    let bases: u64 = spans
        .iter()
        .map(|&(s, e)| e.min(end).saturating_sub(s.max(start)))
        .sum();
    bases as f64 / (end - start) as f64
}

/// Reference bases covered by an alignment's CIGAR.
fn ref_span(
    ops: impl Iterator<Item = std::io::Result<noodles::sam::alignment::record::cigar::Op>>,
) -> u64 {
    ops.filter_map(|op| op.ok())
        .filter(|op| {
            matches!(
                op.kind(),
                Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch
            )
        })
        .map(|op| op.len() as u64)
        .sum()
}

/// Names of reads in [start, end) whose SA:Z tag has an alignment within
/// `pad` of `partner` (chrom, 0-based position).
#[allow(clippy::too_many_arguments)]
fn split_reads_to_partner(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    start: u64,
    end: u64,
    min_mapq: u8,
    partner: &(String, u64),
    pad: u64,
) -> Result<HashSet<String>> {
    let mut names = HashSet::new();
    let points_to_partner =
        |sa: Option<String>| sa.is_some_and(|sa| sa_points_near(&sa, &partner.0, partner.1, pad));

    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let start_pos = crate::extract::safe_noodles_position(start + 1);
        let end_pos = crate::extract::safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let (mut reader, header) =
            crate::extract::open_cram_reader_for_region(bam_path, &repository, &region)
                .context("failed to open CRAM for SA tag counting")?;
        let query = reader.query(&header, &region)?;
        // `query` has already rejected an unknown contig, so this is `Some`.
        let queried_reference_sequence_id =
            header.reference_sequences().get_index_of(chrom.as_bytes());

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            // A container holding several contigs is decoded whole and `Query`
            // filters on coordinates alone, so another contig's reads would be
            // counted as this region's (L2, N4). The BAM arm below needs no
            // such guard: noodles-bam's `Query` compares the reference id.
            if !record_is_on_queried_reference(&buf, queried_reference_sequence_id) {
                continue;
            }
            let flags = buf.flags();
            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_duplicate()
                || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = buf.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            if points_to_partner(sa_value_buf(&buf)) {
                if let Some(n) = buf.name() {
                    names.insert(String::from_utf8_lossy(n.as_ref()).into_owned());
                }
            }
        }
    } else {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .context("failed to open BAM for SA tag counting")?;
        let header = reader.read_header()?;

        let start_pos = crate::extract::safe_noodles_position(start + 1);
        let end_pos = crate::extract::safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let record = rec_result?;
            let flags = record.flags();
            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_duplicate()
                || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = record.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            if points_to_partner(sa_value_bam(&record)) {
                if let Some(n) = record.name() {
                    names.insert(String::from_utf8_lossy(n.as_ref()).into_owned());
                }
            }
        }
    }

    Ok(names)
}

/// Names of reads within `pad` of 0-based `pos` on `chrom` whose alignment
/// carries at least `min_len` bases of sequence the reference has not got
/// there (see [`cigar_shows_insertion_near`]).
fn reads_with_inserted_sequence(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    pos: u64,
    pad: u64,
    min_len: u64,
    min_mapq: u8,
) -> Result<HashSet<String>> {
    let mut names = HashSet::new();
    let start = pos.saturating_sub(pad);
    let end = pos.saturating_add(pad);

    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let start_pos = crate::extract::safe_noodles_position(start + 1);
        let end_pos = crate::extract::safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let (mut reader, header) =
            crate::extract::open_cram_reader_for_region(bam_path, &repository, &region)
                .context("failed to open CRAM for the insertion check")?;
        let query = reader.query(&header, &region)?;
        // `query` has already rejected an unknown contig, so this is `Some`.
        let queried_reference_sequence_id =
            header.reference_sequences().get_index_of(chrom.as_bytes());

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            // A container holding several contigs is decoded whole and `Query`
            // filters on coordinates alone (L2, N4).
            if !record_is_on_queried_reference(&buf, queried_reference_sequence_id) {
                continue;
            }
            if !usable_alignment(buf.flags(), buf.mapping_quality().map(u8::from), min_mapq) {
                continue;
            }
            let Some(align_start) = buf.alignment_start() else {
                continue;
            };
            let align_start = usize::from(align_start).saturating_sub(1) as u64;
            let cigar = buf.cigar();
            if cigar_shows_insertion_near(
                CigarTrait::iter(&cigar),
                align_start,
                pos,
                pad,
                min_len,
            ) {
                if let Some(n) = buf.name() {
                    names.insert(String::from_utf8_lossy(n.as_ref()).into_owned());
                }
            }
        }
    } else {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .context("failed to open BAM for the insertion check")?;
        let header = reader.read_header()?;

        let start_pos = crate::extract::safe_noodles_position(start + 1);
        let end_pos = crate::extract::safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let record = rec_result?;
            if !usable_alignment(
                record.flags(),
                record.mapping_quality().map(u8::from),
                min_mapq,
            ) {
                continue;
            }
            let align_start = match record.alignment_start() {
                Some(Ok(p)) => usize::from(p).saturating_sub(1) as u64,
                _ => continue,
            };
            let cigar = record.cigar();
            if cigar_shows_insertion_near(cigar.iter(), align_start, pos, pad, min_len) {
                if let Some(n) = record.name() {
                    names.insert(String::from_utf8_lossy(n.as_ref()).into_owned());
                }
            }
        }
    }

    Ok(names)
}

/// What `for_each_alignment` hands each record to: its read name, 0-based
/// alignment start, CIGAR operations and bases.
pub(crate) type AlignmentVisitor<'a> =
    dyn FnMut(&[u8], u64, &[noodles::sam::alignment::record::cigar::Op], &[u8]) + 'a;

/// Walk every usable, named alignment overlapping 0-based `[start, end)`,
/// handing the visitor each record's read name, 0-based alignment start, CIGAR
/// operations and bases.
///
/// The BAM and the CRAM reader hand out different record types, so the scans
/// above each carry their own copy of this query; a check that needs nothing
/// from a record but its name, CIGAR and bases can share one.
pub(crate) fn for_each_alignment(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    start: u64,
    end: u64,
    min_mapq: u8,
    visit: &mut AlignmentVisitor<'_>,
) -> Result<()> {
    scan_region(
        bam_path,
        ref_path,
        chrom,
        start,
        end,
        &|flags, mapq| usable_alignment(flags, mapq, min_mapq),
        visit,
    )
}

/// Which records a [`scan_region`] hands on, from their flags and MAPQ (`None`
/// when the record says 255, "unavailable").
type RecordFilter<'a> = dyn Fn(noodles::sam::alignment::record::Flags, Option<u8>) -> bool + 'a;

/// [`for_each_alignment`] with the record filter given rather than fixed to
/// [`usable_alignment`]: `ins_planted` counts spike's reads however the aligner
/// flagged or scored them.
fn scan_region(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    start: u64,
    end: u64,
    keep: &RecordFilter<'_>,
    visit: &mut AlignmentVisitor<'_>,
) -> Result<()> {
    let start_pos = crate::extract::safe_noodles_position(start + 1);
    let end_pos = crate::extract::safe_noodles_position(end);
    let region = noodles::core::Region::new(chrom, start_pos..=end_pos);

    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let (mut reader, header) =
            crate::extract::open_cram_reader_for_region(bam_path, &repository, &region)
                .context("failed to open CRAM for a read scan")?;
        let query = reader.query(&header, &region)?;
        // `query` has already rejected an unknown contig, so this is `Some`.
        let queried_reference_sequence_id =
            header.reference_sequences().get_index_of(chrom.as_bytes());

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            // A container holding several contigs is decoded whole and `Query`
            // filters on coordinates alone (L2, N4).
            if !record_is_on_queried_reference(&buf, queried_reference_sequence_id) {
                continue;
            }
            if !keep(buf.flags(), buf.mapping_quality().map(u8::from)) {
                continue;
            }
            let Some(name) = buf.name() else {
                continue;
            };
            let Some(align_start) = buf.alignment_start() else {
                continue;
            };
            let cigar = buf.cigar();
            let ops: Vec<_> = CigarTrait::iter(&cigar).collect::<std::io::Result<_>>()?;
            visit(
                name.as_ref(),
                usize::from(align_start).saturating_sub(1) as u64,
                &ops,
                buf.sequence().as_ref(),
            );
        }
    } else {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .context("failed to open BAM for a read scan")?;
        let header = reader.read_header()?;
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let record = rec_result?;
            if !keep(record.flags(), record.mapping_quality().map(u8::from)) {
                continue;
            }
            let Some(name) = record.name() else {
                continue;
            };
            let align_start = match record.alignment_start() {
                Some(Ok(p)) => usize::from(p).saturating_sub(1) as u64,
                _ => continue,
            };
            let cigar = record.cigar();
            let ops: Vec<_> = cigar.iter().collect::<std::io::Result<_>>()?;
            let seq: Vec<u8> = record.sequence().iter().collect();
            visit(name.as_ref(), align_start, &ops, &seq);
        }
    }

    Ok(())
}

/// The record filter every scan here shares: a primary, mapped, non-duplicate
/// alignment that clears `min_mapq`.
fn usable_alignment(
    flags: noodles::sam::alignment::record::Flags,
    mapping_quality: Option<u8>,
    min_mapq: u8,
) -> bool {
    !(flags.is_unmapped()
        || flags.is_secondary()
        || flags.is_supplementary()
        || flags.is_duplicate()
        || flags.is_qc_fail())
        && mapping_quality.unwrap_or(0) >= min_mapq
}

/// True if an SA:Z value (`chrom,pos,strand,CIGAR,mapQ,NM;...`, 1-based pos)
/// has an alignment on `chrom` within `pad` of 0-based `pos`.
fn sa_points_near(sa: &str, chrom: &str, pos: u64, pad: u64) -> bool {
    sa.split(';').filter(|e| !e.is_empty()).any(|entry| {
        let mut fields = entry.split(',');
        let (Some(c), Some(p)) = (fields.next(), fields.next()) else {
            return false;
        };
        match p.parse::<u64>() {
            Ok(p1) => c == chrom && p1.saturating_sub(1).abs_diff(pos) <= pad,
            Err(_) => false,
        }
    })
}

/// Partner breakpoint (chrom, 0-based) from a BND ALT such as `A]chr2:500]`.
fn bnd_partner(alt: &str) -> Option<(String, u64)> {
    let open = alt.find(['[', ']'])?;
    let close = alt.rfind(['[', ']'])?;
    let (chrom, pos) = alt.get(open + 1..close)?.rsplit_once(':')?;
    let pos: u64 = pos.parse().ok()?;
    Some((chrom.to_string(), pos.checked_sub(1)?))
}

/// The SA:Z value of a BAM record, if any.
fn sa_value_bam(record: &noodles::bam::Record) -> Option<String> {
    use noodles::sam::alignment::record::data::field::{Tag, Value};
    match record.data().get(&Tag::new(b'S', b'A'))? {
        Ok(Value::String(s)) => Some(String::from_utf8_lossy(s).into_owned()),
        _ => None,
    }
}

/// The SA:Z value of a RecordBuf (CRAM), if any.
fn sa_value_buf(buf: &noodles::sam::alignment::RecordBuf) -> Option<String> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    match buf.data().get(&Tag::new(b'S', b'A'))? {
        Value::String(s) => Some(String::from_utf8_lossy(s).into_owned()),
        _ => None,
    }
}

/// Pileup a region to count alleles at each position.
/// Reuses the same CIGAR walking pattern as loh.rs.
#[allow(clippy::too_many_arguments)]
fn pileup_region(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    region_start: u64,
    region_end: u64,
    min_mapq: u8,
    allele_counts: &mut HashMap<u64, [u32; 4]>,
    read_alleles: &mut HashMap<String, Vec<(u64, u8)>>,
) -> Result<()> {
    if crate::extract::is_cram(bam_path) {
        let repository = crate::extract::build_fasta_repository(ref_path)?;
        let start_pos = crate::extract::safe_noodles_position(region_start + 1);
        let end_pos = crate::extract::safe_noodles_position(region_end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let (mut reader, header) =
            crate::extract::open_cram_reader_for_region(bam_path, &repository, &region)
                .context("failed to open CRAM for pileup")?;
        let query = reader.query(&header, &region)?;
        // `query` has already rejected an unknown contig, so this is `Some`.
        let queried_reference_sequence_id =
            header.reference_sequences().get_index_of(chrom.as_bytes());

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record.try_into_alignment_record(&header)?;
            // A container holding several contigs is decoded whole and `Query`
            // filters on coordinates alone, so another contig's reads would be
            // counted as this region's (L2, N4). The BAM arm below needs no
            // such guard: noodles-bam's `Query` compares the reference id.
            if !record_is_on_queried_reference(&buf, queried_reference_sequence_id) {
                continue;
            }
            let flags = buf.flags();
            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_duplicate()
                || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = buf.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            let name = match buf.name() {
                Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
                None => continue,
            };
            let align_start = match buf.alignment_start() {
                Some(p) => usize::from(p).saturating_sub(1) as u64,
                None => continue,
            };
            let seq: Vec<u8> = buf.sequence().as_ref().to_vec();
            let cigar = buf.cigar();
            walk_cigar_pileup(
                &seq,
                CigarTrait::iter(&cigar),
                align_start,
                region_start,
                region_end,
                &name,
                allele_counts,
                read_alleles,
            );
        }
    } else {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .context("failed to open BAM for pileup")?;
        let header = reader.read_header()?;

        let start_pos = crate::extract::safe_noodles_position(region_start + 1);
        let end_pos = crate::extract::safe_noodles_position(region_end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let record = rec_result?;
            let flags = record.flags();
            if flags.is_unmapped()
                || flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_duplicate()
                || flags.is_qc_fail()
            {
                continue;
            }
            let mq: u8 = record.mapping_quality().map(u8::from).unwrap_or(0);
            if mq < min_mapq {
                continue;
            }
            let name = match record.name() {
                Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
                None => continue,
            };
            let align_start = match record.alignment_start() {
                Some(Ok(p)) => usize::from(p).saturating_sub(1) as u64,
                _ => continue,
            };
            let seq: Vec<u8> = record.sequence().iter().collect();
            let cigar = record.cigar();
            walk_cigar_pileup(
                &seq,
                Box::new(cigar.iter()),
                align_start,
                region_start,
                region_end,
                &name,
                allele_counts,
                read_alleles,
            );
        }
    }

    Ok(())
}

/// Walk CIGAR operations and collect allele counts + per-read alleles.
/// Local copy of the same logic from loh.rs.
#[allow(clippy::too_many_arguments)]
fn walk_cigar_pileup(
    seq: &[u8],
    cigar_ops: Box<
        dyn Iterator<Item = std::io::Result<noodles::sam::alignment::record::cigar::Op>> + '_,
    >,
    align_start: u64,
    region_start: u64,
    region_end: u64,
    name: &str,
    allele_counts: &mut HashMap<u64, [u32; 4]>,
    read_alleles: &mut HashMap<String, Vec<(u64, u8)>>,
) {
    let mut ref_pos = align_start;
    let mut seq_pos = 0usize;

    for op_result in cigar_ops {
        let op = match op_result {
            Ok(o) => o,
            Err(_) => break,
        };
        let len = op.len();

        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..len {
                    let rp = ref_pos + i as u64;
                    if rp >= region_start && rp < region_end {
                        if let Some(&base) = seq.get(seq_pos + i) {
                            let base = base.to_ascii_uppercase();
                            if let Some(idx) = pileup_base_index(base) {
                                allele_counts.entry(rp).or_insert([0; 4])[idx] += 1;
                                read_alleles
                                    .entry(name.to_string())
                                    .or_default()
                                    .push((rp, base));
                            }
                        }
                    }
                }
                ref_pos += len as u64;
                seq_pos += len;
            }
            Kind::Insertion | Kind::SoftClip => {
                seq_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos += len as u64;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
}

// ---------------------------------------------------------------------------
// VCF INFO helpers (local copies from vcf_input.rs)
// ---------------------------------------------------------------------------

fn parse_info_field<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    for field in info.split(';') {
        if let Some(value) = field
            .strip_prefix(key)
            .and_then(|rest| rest.strip_prefix('='))
        {
            return Some(value);
        }
    }
    None
}

fn parse_info_u64(info: &str, key: &str) -> Option<u64> {
    parse_info_field(info, key)?.parse().ok()
}

// ---------------------------------------------------------------------------
// Output formatting
// ---------------------------------------------------------------------------

fn format_event_label(event: &TruthEvent) -> String {
    let region = if event.start == event.end {
        format!("{}:{}", event.chrom, event.start)
    } else {
        format!("{}:{}-{}", event.chrom, event.start, event.end)
    };
    format!("{} {} ({})", event.sv_type, region, event.gene)
}

fn print_results(results: &[CheckResult], json: bool, strict: bool) -> Result<()> {
    let n_pass = results.iter().filter(|r| r.pass).count();
    let n_total = results.len();
    let stdout = std::io::stdout();
    let mut out = stdout.lock();

    if json {
        print_results_json(&mut out, results, n_total, n_pass, strict)?;
    } else {
        print_results_text(&mut out, results, n_total, n_pass, strict)?;
    }

    Ok(())
}

fn print_results_text(
    out: &mut impl Write,
    results: &[CheckResult],
    n_total: usize,
    n_pass: usize,
    strict: bool,
) -> Result<()> {
    writeln!(out)?;
    writeln!(out, "spike validate -- {} checks", n_total)?;
    writeln!(out)?;
    writeln!(
        out,
        "{:<35} {:<18} {:<25} {:<15} Status",
        "Event", "Check", "Expected", "Observed"
    )?;
    writeln!(out, "{}", "-".repeat(100))?;

    for r in results {
        // Only the last column changes: an advisory row is marked where it
        // stands, so every other column keeps the width it had.
        let status = match (r.pass, r.advisory) {
            (true, false) => "PASS",
            (false, false) => "FAIL",
            (true, true) => "PASS (advisory)",
            (false, true) => "FAIL (advisory)",
        };
        writeln!(
            out,
            "{:<35} {:<18} {:<25} {:<15} {}",
            truncate(&r.event_label, 34),
            r.check_name,
            truncate(&r.expected, 24),
            truncate(&r.observed, 14),
            status,
        )?;
    }

    writeln!(out)?;
    writeln!(out, "Result: {}/{} PASS", n_pass, n_total,)?;

    // Only when there is one: a truth VCF with no census carries no advisory
    // row, and its report reads as it always did.
    let n_advisory = results.iter().filter(|r| r.advisory).count();
    if n_advisory > 0 {
        let n_advisory_pass = results.iter().filter(|r| r.advisory && r.pass).count();
        writeln!(
            out,
            "Advisory: {} checks, {} PASS, {} FAIL {}",
            n_advisory,
            n_advisory_pass,
            n_advisory - n_advisory_pass,
            if strict {
                "(in the exit status: --strict)"
            } else {
                "(not in the exit status; --strict includes them)"
            },
        )?;
    }

    Ok(())
}

fn print_results_json(
    out: &mut impl Write,
    results: &[CheckResult],
    n_total: usize,
    n_pass: usize,
    strict: bool,
) -> Result<()> {
    // `total`, `pass` and `fail` are over every row printed, and keep exactly
    // the meaning they had before the advisory rows existed -- a parser written
    // against them still reads what it always read. The `counted_*` trio beside
    // them is over the rows the *exit status* is computed from: the non-advisory
    // rows, or every row under `--strict`, which `strict` names. They come from
    // `counted_rows`, the same filter `failure_message` exits on, so a consumer
    // reading `counted_pass == 0` is reading the number the exit status agrees
    // with. `scripts/validate_pipeline.sh`'s step-5 guard is that consumer: an
    // advisory row is read back from the truth VCF and passes whatever the BAM
    // holds, so `pass` alone cannot tell a working run from a header-only BAM.
    let n_counted_total = counted_rows(results, strict).count();
    let n_counted_pass = counted_rows(results, strict).filter(|r| r.pass).count();

    // Manual JSON to avoid serde dependency.
    writeln!(out, "{{")?;
    writeln!(
        out,
        "  \"summary\": {{ \"total\": {}, \"pass\": {}, \"fail\": {}, \"counted_total\": {}, \"counted_pass\": {}, \"counted_fail\": {}, \"strict\": {} }},",
        n_total,
        n_pass,
        n_total - n_pass,
        n_counted_total,
        n_counted_pass,
        n_counted_total - n_counted_pass,
        strict,
    )?;
    writeln!(out, "  \"checks\": [")?;

    for (i, r) in results.iter().enumerate() {
        let comma = if i + 1 < results.len() { "," } else { "" };
        writeln!(
            out,
            "    {{ \"event\": \"{}\", \"check\": \"{}\", \"expected\": \"{}\", \"observed\": \"{}\", \"pass\": {}, \"advisory\": {} }}{}",
            escape_json(&r.event_label),
            escape_json(&r.check_name),
            escape_json(&r.expected),
            escape_json(&r.observed),
            r.pass,
            r.advisory,
            comma,
        )?;
    }

    writeln!(out, "  ]")?;
    writeln!(out, "}}")?;

    Ok(())
}

fn truncate(s: &str, max_len: usize) -> String {
    if s.len() <= max_len {
        s.to_string()
    } else {
        // Back off to the nearest char boundary at or before the cut point:
        // slicing mid-character (e.g. a non-ASCII gene name) would panic.
        let mut end = max_len.saturating_sub(3).min(s.len());
        while end > 0 && !s.is_char_boundary(end) {
            end -= 1;
        }
        format!("{}...", &s[..end])
    }
}

fn escape_json(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::Op;

    #[test]
    fn test_sa_points_near_partner_only() {
        let sa = "chr20,45000101,-,100S50M,60,0;chr20,41000001,+,100M50S,60,1;";
        assert!(sa_points_near(sa, "chr20", 45000000, 500));
        assert!(sa_points_near(sa, "chr20", 41000000, 500)); // second entry
        assert!(!sa_points_near(sa, "chr20", 43000000, 500)); // elsewhere
        assert!(!sa_points_near(sa, "chr9", 45000000, 500)); // other chrom
        assert!(!sa_points_near("garbage", "chr20", 45000000, 500));
    }

    #[test]
    fn test_truth_events_know_their_partner_breakpoint() {
        let vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr20\t41000000\tsim_del_1\tA\t<DEL>\t999\tPASS\tSVTYPE=DEL;END=41010000;SVLEN=-10000\tGT\t0/1
chr20\t41000000\tsim_fus_2\tA\tA]chr20:45000000]\t999\tPASS\tSVTYPE=BND;MATEID=sim_fus_2_mate\tGT\t0/1
chr20\t42000000\tsim_ins_3\tA\t<INS>\t999\tPASS\tSVTYPE=INS;SVLEN=500\tGT\t0/1
";
        let path = std::env::temp_dir().join(format!("spike_partner_{}.vcf", std::process::id()));
        std::fs::write(&path, vcf).unwrap();
        let events = load_truth_events(path.to_str().unwrap()).unwrap();
        std::fs::remove_file(&path).ok();

        assert_eq!(events[0].partner, Some(("chr20".to_string(), 41010000)));
        assert_eq!(events[1].partner, Some(("chr20".to_string(), 44999999)));
        assert_eq!(events[2].partner, None);
    }

    #[test]
    fn test_mean_depth_counts_aligned_bases() {
        // Two reads over [0,100): 100 + 50 bases inside -> depth 1.5.
        assert_eq!(mean_depth(&[(0, 100), (50, 150)], 0, 100), 1.5);
        // 35x of 150 bp reads, one start every 150/35 bp: depth ~35, not ~0.23.
        let spans: Vec<(u64, u64)> = (0..3500u64).map(|i| {
            let s = i * 150 / 35;
            (s, s + 150)
        }).collect();
        let d = mean_depth(&spans, 1000, 2000);
        assert!((d - 35.0).abs() < 0.5, "depth {}", d);
    }

    #[test]
    fn test_coverage_ratio_without_flank_coverage_fails() {
        let r = coverage_ratio_result("DEL".to_string(), "coverage_ratio", "DEL", 0.5, 0.0, 0.0);
        assert!(!r.pass, "an unevaluable coverage check must not pass");
    }

    #[test]
    fn test_an_unevaluable_coverage_row_keeps_the_name_it_was_given() {
        // The `flank_depth < 1.0` branch builds a row of its own, so it has to
        // carry the caller's name like the graded one does. A DEL whose flanks
        // are under 1x otherwise prints two rows both called `coverage_ratio`,
        // with different `pass` values, in the table and in `--json`, where a
        // consumer selecting by check name takes whichever it hits first (F2).
        let r = coverage_ratio_result("DEL".to_string(), COVERAGE_ANY_MAPQ, "DEL", 1.0, 0.0, 0.0);

        assert_eq!(r.check_name, "coverage_any_mapq");
        assert_eq!(r.observed, "no flanking coverage");
        assert!(!r.pass, "an unevaluable coverage check must not pass");
    }

    #[test]
    fn test_coverage_ratio_has_no_verdict_for_a_type_without_an_expected_ratio() {
        // Only DEL and DUP reach this check. Any other type has no expected
        // ratio, so an untouched region must not pass it by default (NF6).
        let r = coverage_ratio_result("INV".to_string(), "coverage_ratio", "INV", 0.5, 35.0, 35.0);
        assert!(!r.pass, "observed {}, expected {}", r.observed, r.expected);
        assert_eq!(r.expected, "N/A");
    }

    #[test]
    fn test_coverage_ratio_tells_deletion_from_untouched() {
        assert!(coverage_ratio_result("DEL".to_string(), "coverage_ratio", "DEL", 0.5, 17.5, 35.0).pass);
        assert!(!coverage_ratio_result("DEL".to_string(), "coverage_ratio", "DEL", 0.5, 35.0, 35.0).pass);
    }

    #[test]
    fn test_check_that_cannot_run_is_a_failure() {
        // e.g. truth VCF uses "20" but the BAM uses "chr20".
        let r = check_outcome(
            "DEL 20:100-200",
            "split_reads",
            Err(anyhow::anyhow!("reference sequence not found: 20")),
            false,
        );
        assert!(!r.pass);
        assert_eq!(r.check_name, "split_reads");
        assert!(r.observed.contains("reference sequence not found: 20"), "{}", r.observed);
    }

    #[test]
    fn test_parse_info_field() {
        assert_eq!(
            parse_info_field("SVTYPE=DEL;END=100;SIM_VAF=0.500", "END"),
            Some("100")
        );
        assert_eq!(
            parse_info_field("SVTYPE=DEL;END=100;SIM_VAF=0.500", "SIM_VAF"),
            Some("0.500")
        );
        assert_eq!(parse_info_field("SVTYPE=DEL;END=100", "MISSING"), None);
    }

    #[test]
    fn test_format_event_label() {
        let event = TruthEvent {
            chrom: "chr7".to_string(),
            start: 55000,
            end: 56000,
            sv_type: "DEL".to_string(),
            expected_vaf: 0.5,
            gene: "EGFR".to_string(),
            partner: None,
            ref_allele: None,
            alt_allele: None,
            ins_len: None,
            ins_alt: None,
            sim_number: None,
            census: CensusInfo::default(),
        };
        assert_eq!(format_event_label(&event), "DEL chr7:55000-56000 (EGFR)");
    }

    #[test]
    fn test_format_event_label_point() {
        let event = TruthEvent {
            chrom: "chr7".to_string(),
            start: 55200,
            end: 55200,
            sv_type: "INS".to_string(),
            expected_vaf: 0.3,
            gene: "EGFR".to_string(),
            partner: None,
            ref_allele: None,
            alt_allele: None,
            ins_len: None,
            ins_alt: None,
            sim_number: None,
            census: CensusInfo::default(),
        };
        assert_eq!(format_event_label(&event), "INS chr7:55200 (EGFR)");
    }

    #[test]
    fn test_load_truth_events_del() {
        let vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr7\t55000\tsim_del_1\tN\t<DEL>\t999\tPASS\tSVTYPE=DEL;END=56000;SVLEN=-1000;SIM_VAF=0.500;SIM_GENE=EGFR;SIM_EXONS=.\tGT\t0/1
";
        // Write to temp file.
        let dir = std::env::temp_dir().join("spike_test_validate");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("truth_del.vcf");
        std::fs::write(&path, vcf).unwrap();

        let events = load_truth_events(path.to_str().unwrap()).unwrap();
        assert_eq!(events.len(), 1);
        assert_eq!(events[0].sv_type, "DEL");
        assert_eq!(events[0].chrom, "chr7");
        assert_eq!(events[0].start, 55000);
        assert_eq!(events[0].end, 56000);
        assert!((events[0].expected_vaf - 0.5).abs() < f64::EPSILON);
        assert_eq!(events[0].gene, "EGFR");
    }

    #[test]
    fn test_load_truth_events_snp() {
        let vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr7\t55201\tsim_var_1\tA\tT\t999\tPASS\tSIM_VAF=0.500;SIM_GENE=EGFR\tGT\t0/1
";
        let dir = std::env::temp_dir().join("spike_test_validate");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("truth_snp.vcf");
        std::fs::write(&path, vcf).unwrap();

        let events = load_truth_events(path.to_str().unwrap()).unwrap();
        assert_eq!(events.len(), 1);
        assert_eq!(events[0].sv_type, "SNP");
        assert_eq!(events[0].start, 55200); // 1-based 55201 → 0-based 55200
        assert_eq!(events[0].ref_allele.as_deref(), Some(b"A".as_slice()));
        assert_eq!(events[0].alt_allele.as_deref(), Some(b"T".as_slice()));
    }

    #[test]
    fn test_load_truth_events_bnd_pair() {
        let vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr2\t29416089\tsim_fus_1\tN\tN[chr2:42522656[\t999\tPASS\tSVTYPE=BND;MATEID=sim_fus_1_mate;SIM_VAF=0.050;SIM_GENE=ALK\tGT\t0/1
chr2\t42522656\tsim_fus_1_mate\tN\t]chr2:29416089]N\t999\tPASS\tSVTYPE=BND;MATEID=sim_fus_1;SIM_VAF=0.050;SIM_GENE=EML4\tGT\t0/1
";
        let dir = std::env::temp_dir().join("spike_test_validate");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("truth_bnd.vcf");
        std::fs::write(&path, vcf).unwrap();

        let events = load_truth_events(path.to_str().unwrap()).unwrap();
        // BND pairs should produce only 1 event (mate is deduplicated).
        assert_eq!(events.len(), 1);
        assert_eq!(events[0].sv_type, "BND");
    }

    #[test]
    fn test_truncate() {
        assert_eq!(truncate("short", 10), "short");
        assert_eq!(truncate("this is a longer string", 10), "this is...");
    }

    #[test]
    fn test_escape_json() {
        assert_eq!(escape_json("hello \"world\""), "hello \\\"world\\\"");
        assert_eq!(escape_json("a\\b"), "a\\\\b");
    }

    #[test]
    fn test_truncate_does_not_panic_on_multibyte_char_boundary() {
        // "IFN-\u{3b3}" (interferon gamma) is a real gene alias whose Greek
        // letter is a 2-byte UTF-8 character. The naive byte slice used to
        // land inside that character and panic ("byte index N is not a char
        // boundary"); it must instead back off to the previous boundary.
        let name = "IFN-\u{3b3}-associated-deletion-event";
        assert_eq!(truncate(name, 8), "IFN-...");
    }

    #[test]
    fn test_escape_json_escapes_all_control_chars() {
        assert_eq!(escape_json("a\tb\rc\nd"), "a\\tb\\rc\\nd");
        // Every other C0 control character below 0x20 must become \u00XX,
        // as bare JSON requires -- \u{1} (SOH) has no short-form escape.
        assert_eq!(escape_json("a\u{1}b"), "a\\u0001b");
    }

    #[test]
    fn test_escape_json_uses_the_short_escapes_for_backspace_and_form_feed() {
        // \b (0x08) and \f (0x0c) have short forms in JSON, and the escaper
        // emits them rather than \u0008 / \u000c. Nothing exercised those
        // two arms: the tests above cover \t, \r, \n and the \u00XX
        // fallback only.
        assert_eq!(escape_json("a\u{8}b\u{c}c"), "a\\bb\\fc");
    }

    #[test]
    fn test_truncate_backs_off_across_a_four_byte_character() {
        // The cut point can land 1, 2 or 3 bytes into a character. An emoji
        // in a gene name is far-fetched, but a 4-byte codepoint is the widest
        // UTF-8 gets and the deepest the walk-back loop has to go; the test
        // above only ever backs off one byte inside a 2-byte character.
        let name = "AB\u{1f600}CDEFGHIJ";
        assert_eq!(name.len(), 14, "2 + 4 + 8 bytes");
        // max_len 8 cuts at byte 5, which is inside the emoji's bytes 2..6.
        assert_eq!(truncate(name, 8), "AB...");
    }

    // --- N4: `validate`'s CRAM queries must not read another contig's reads ---

    /// The shared two-contig CRAM fixture in a scratch directory of its own.
    /// Returns `(dir, fasta_path, cram_path)`; the caller removes `dir`.
    fn two_contig_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let dir = std::env::temp_dir().join(format!(
            "spike_test_validate_{}_{}",
            tag,
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let (fasta, cram) = crate::extract::test_fixtures::write_two_contig_cram(&dir);
        (dir, fasta, cram)
    }

    #[test]
    fn test_count_depth_in_region_skips_another_contigs_reads() {
        // chrA's six reads cover 600 of the region's 2000 bases, so its mean
        // depth is 0.3. chrB's four reads cover another 400 bases of the same
        // window from the shared container; counted, they lift the answer to
        // 0.5 and with it every coverage-ratio verdict (L2).
        let (dir, fasta, cram) = two_contig_cram("depth");

        let depth = count_depth_in_region(&cram, &fasta, "chrA", 0, 2000, 20).unwrap();

        assert_eq!(
            depth, 0.3,
            "only chrA's reads may count toward chrA's depth"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_split_reads_to_partner_skips_another_contigs_reads() {
        // Every fixture record carries the same SA:Z, so what comes back is
        // exactly "the records this query saw" — which is what makes a
        // foreign contig's read names visible here.
        let (dir, fasta, cram) = two_contig_cram("split_reads");

        let found = split_reads_to_partner(
            &cram,
            &fasta,
            "chrA",
            0,
            2000,
            20,
            &("chrA".to_string(), 1500),
            10,
        )
        .unwrap();

        let mut names: Vec<String> = found.into_iter().collect();
        names.sort();
        assert_eq!(
            names,
            vec!["chrA_pair0", "chrA_pair1", "chrA_pair2"],
            "a chrB read's name must not be counted as a chrA split read"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_pileup_region_skips_another_contigs_reads() {
        // The allele-fraction check's pileup. chrA and chrB cover disjoint
        // intervals, so every count over chrB's 301-400, 501-600 and 701-800
        // is a base that arrived from the other contig.
        let (dir, fasta, cram) = two_contig_cram("pileup");

        let mut allele_counts: HashMap<u64, [u32; 4]> = HashMap::new();
        let mut read_alleles: HashMap<String, Vec<(u64, u8)>> = HashMap::new();
        pileup_region(
            &cram,
            &fasta,
            "chrA",
            0,
            2000,
            20,
            &mut allele_counts,
            &mut read_alleles,
        )
        .unwrap();

        // The three intervals only chrB covers hold all 400 of its bases.
        let foreign: u32 = (300..400u64)
            .chain(500..600)
            .chain(700..800)
            .filter_map(|pos| allele_counts.get(&pos))
            .map(|c| c.iter().sum::<u32>())
            .sum();
        assert_eq!(foreign, 0, "chrB's bases must not enter the chrA pileup");

        let mut names: Vec<String> = read_alleles.keys().cloned().collect();
        names.sort();
        assert_eq!(
            names,
            vec!["chrA_pair0", "chrA_pair1", "chrA_pair2"],
            "only chrA reads may carry an allele in a chrA pileup"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    // --- L15: the global checks read the event regions, not the file's head ---

    /// A one-contig CRAM shaped like a whole-genome BAM: a head nobody asked
    /// about, and an event window far from it. 15 pairs at chrA:101-1800 carry
    /// MAPQ 0 (a telomere's worth of multi-mapping reads) and 3 pairs inside
    /// chrA:10001-10500 carry MAPQ 60, so the mean MAPQ of all 36 records is
    /// exactly 10.0 -- the number REVIEW.md measured on whole-genome HG002 --
    /// while the event's own window reads 60. No record carries the duplicate
    /// flag. Returns `(dir, fasta_path, cram_path)`; the caller removes `dir`.
    fn head_and_event_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let seq = cycling_contig();
        let dir = fixture_dir(tag);
        let fasta_path = write_one_contig_fasta(&dir, "one_contig", &seq);
        let header = one_contig_header(seq.len());

        // One pair = two records, read1 forward (0x63) and read2 reverse
        // (0x93), both properly segmented, 300 bp apart.
        let record = |name: &str, start: usize, first: bool, mapq: u8| {
            let (pos, mate_pos) = if first {
                (start, start + 200)
            } else {
                (start + 200, start)
            };
            let span = 200 + TEST_READ_LEN;
            noodles::cram::Record::builder()
                .set_bam_flags(noodles::sam::alignment::record::Flags::from(if first {
                    0x63u16
                } else {
                    0x93u16
                }))
                .set_flags(noodles::cram::record::Flags::QUALITY_SCORES_STORED_AS_ARRAY)
                .set_reference_sequence_id(0)
                .set_read_length(TEST_READ_LEN)
                .set_alignment_start(noodles::core::Position::new(pos).unwrap())
                .set_name(name)
                .set_next_fragment_reference_sequence_id(0)
                .set_next_mate_alignment_start(noodles::core::Position::new(mate_pos).unwrap())
                .set_template_size(if first { span as i32 } else { -(span as i32) })
                .set_mapping_quality(
                    noodles::sam::alignment::record::MappingQuality::new(mapq).unwrap(),
                )
                .set_bases(noodles::sam::alignment::record_buf::Sequence::from(
                    seq[pos - 1..pos - 1 + TEST_READ_LEN].to_vec(),
                ))
                .set_quality_scores(noodles::sam::alignment::record_buf::QualityScores::from(
                    vec![40u8; TEST_READ_LEN],
                ))
                .build()
        };

        // 15 pairs of MAPQ 0 at the contig head, then 3 pairs of MAPQ 60 over
        // the event.
        let mut records: Vec<noodles::cram::Record> = Vec::new();
        for i in 0..15usize {
            let name = format!("head_pair{}", i);
            let start = 101 + i * 100;
            records.push(record(&name, start, true, 0));
            records.push(record(&name, start, false, 0));
        }
        for i in 0..3usize {
            let name = format!("event_pair{}", i);
            let start = 10_001 + i * 100;
            records.push(record(&name, start, true, 60));
            records.push(record(&name, start, false, 60));
        }
        let cram_path = write_indexed_cram(&dir, "head_and_event", &fasta_path, &header, &records);

        (
            dir,
            fasta_path.to_str().unwrap().to_string(),
            cram_path.to_str().unwrap().to_string(),
        )
    }

    /// A het DEL truth event over [start, end).
    fn del_event(chrom: &str, start: u64, end: u64) -> TruthEvent {
        TruthEvent {
            chrom: chrom.to_string(),
            start,
            end,
            sv_type: "DEL".to_string(),
            expected_vaf: 0.5,
            gene: "unknown".to_string(),
            partner: Some((chrom.to_string(), end)),
            ref_allele: None,
            alt_allele: None,
            ins_len: None,
            ins_alt: None,
            sim_number: None,
            census: CensusInfo::default(),
        }
    }

    #[test]
    fn test_global_sample_reads_the_event_region_not_the_head_of_the_file() {
        // The file's first 30 records are the head's MAPQ 0 pairs and only 6
        // lie in the event's window, so a sample taken from the top of the
        // file reads 10.0 and FAILs the >20 check on a file whose event
        // region is MAPQ 60 throughout (L15).
        let (dir, fasta, cram) = head_and_event_cram("global_sample");
        let event = del_event("chrA", 10_000, 10_200);

        let sample = sample_event_regions(
            &cram,
            &fasta,
            std::slice::from_ref(&event),
            5_000,
            GLOBAL_SAMPLE_MAX,
        )
        .unwrap();
        let result = check_mapq(&sample);

        assert_eq!(
            result.observed, "60.0",
            "mean MAPQ must come from the event's own region"
        );
        assert!(result.pass, "MAPQ 60 must pass the >20 check");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_dup_rate_without_duplicate_flags_is_not_evaluable() {
        // A file with no duplicate flags is not a file without duplicates:
        // nothing marked them, so there is no rate to report and none to
        // pass on (M10).
        let sample = GlobalSample {
            total: 10_000,
            dups: 0,
            mapq_sum: 600_000,
            insert_sizes: vec![300.0; 10_000],
            duplicates_marked: false,
        };

        let result = check_dup_rate(&sample);

        assert_eq!(result.observed, "no dup flags");
        assert!(!result.pass, "an unevaluable duplicate rate may not pass");
    }

    #[test]
    fn test_dup_rate_is_measured_when_records_carry_duplicate_flags() {
        // The other half of the rule: a marked file still gets a number.
        let sample = GlobalSample {
            total: 1_000,
            dups: 71,
            mapq_sum: 60_000,
            insert_sizes: vec![300.0; 1_000],
            duplicates_marked: false,
        };

        let result = check_dup_rate(&sample);

        assert_eq!(result.observed, "7.1%");
        assert!(result.pass, "7.1% is under the 50% limit");
    }

    #[test]
    fn test_insert_size_without_paired_records_is_not_evaluable() {
        // The fallback this check used to inherit from `bam_stats` was
        // 350 +/- 50, which passes: a default reported as a measurement.
        let sample = GlobalSample {
            total: 100,
            dups: 0,
            mapq_sum: 6_000,
            insert_sizes: Vec::new(),
            duplicates_marked: false,
        };

        let result = check_insert_size(&sample);

        assert_eq!(result.observed, "no pairs");
        assert!(!result.pass, "an unevaluable insert size may not pass");
    }

    // --- L15 fix pass 1: every event region is sampled, and the sample's own
    // size and provenance are reported ---

    /// An INS truth event: no coverage ratio, no split reads, no allele
    /// frequency -- no check applies to it at all.
    fn ins_event(chrom: &str, pos: u64) -> TruthEvent {
        TruthEvent {
            chrom: chrom.to_string(),
            start: pos,
            end: pos + 1,
            sv_type: "INS".to_string(),
            expected_vaf: 0.5,
            gene: "unknown".to_string(),
            partner: None,
            ref_allele: None,
            alt_allele: None,
            ins_len: Some(300),
            ins_alt: None,
            sim_number: None,
            census: CensusInfo::default(),
        }
    }

    fn args_for(bam: &str, reference: &str) -> ValidateArgs {
        ValidateArgs {
            bam_path: bam.to_string(),
            truth_path: "unused.vcf".to_string(),
            ref_path: reference.to_string(),
            min_mapq: 20,
            flank_bp: 5_000,
            json_output: false,
            strict: false,
        }
    }

    #[test]
    fn test_every_event_region_gets_a_share_of_the_sample_budget() {
        // Two regions and a budget of 4: with a per-region floor larger than
        // the fair share the first region spends the whole budget and the
        // second is never read (MAPQ 0.0 instead of the 30.0 the two regions
        // average to). A truth VCF with more than 200 events hit exactly this
        // against the real 200 000/1 000 numbers.
        let (dir, fasta, cram) = head_and_event_cram("budget_share");
        let events = vec![
            del_event("chrA", 100, 1_800),     // the head: 30 records at MAPQ 0
            del_event("chrA", 10_000, 10_200), // the event: MAPQ 60
        ];

        let sample = sample_event_regions(&cram, &fasta, &events, 0, 4).unwrap();
        let result = check_mapq(&sample);

        assert_eq!(sample.total, 4, "both regions must contribute their share");
        assert_eq!(
            result.observed, "30.0",
            "the last event region must reach the sample"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_sample_region_reports_how_many_records_it_took() {
        // The count is what says a global check was computed over 41 records
        // rather than over 200 000.
        let (dir, fasta, cram) = head_and_event_cram("region_count");
        let mut sample = GlobalSample::default();

        let n = sample_region(&cram, &fasta, "chrA", 10_000, 10_500, 100, &mut sample).unwrap();

        assert_eq!(n, 6, "3 pairs lie in chrA:10001-10500");
        assert_eq!(n, sample.total);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_one_unqueryable_region_does_not_fail_the_whole_sample() {
        // A truth event on a contig the alignment file does not have used to
        // propagate out of the sample and fail all three global checks, unlike
        // the per-event checks, which isolate an error per check.
        let (dir, fasta, cram) = head_and_event_cram("bad_contig");
        let events = vec![
            del_event("chrZ", 100, 200), // not in the file
            del_event("chrA", 10_000, 10_200),
        ];

        let sampled = sample_event_regions(&cram, &fasta, &events, 5_000, GLOBAL_SAMPLE_MAX);

        assert!(
            sampled.is_ok(),
            "one unqueryable region must not fail the whole sample: {:?}",
            sampled.as_ref().err()
        );
        let sample = sampled.unwrap();
        assert_eq!(check_mapq(&sample).observed, "60.0");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_global_sample_of_no_events_is_an_error() {
        // The empty-truth-VCF bail: without it the three global checks report
        // whole-file statistics for a run that validated nothing (M11).
        let (dir, fasta, cram) = head_and_event_cram("no_events");

        let sampled = sample_event_regions(&cram, &fasta, &[], 5_000, GLOBAL_SAMPLE_MAX);

        assert!(
            sampled.is_err(),
            "a truth VCF with no events has no region to sample"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_event_no_check_applies_to_is_not_evaluable() {
        // A truth VCF whose events no check covers used to score 3/3 PASS:
        // the events were never checked and the three globals passed on
        // background reads (M11's door, one step over). INS is checked now
        // (`ins_reads`), so this uses a type that still reaches no branch.
        // No check reads the BAM here, so the paths need not exist.
        let args = args_for("/nonexistent/no.bam", "/nonexistent/no.fa");
        let event = TruthEvent {
            sv_type: "CNV".to_string(),
            ..ins_event("chrA", 10_000)
        };
        let mut results: Vec<CheckResult> = Vec::new();

        check_event(&args, &event, &NearbyRecords::default(), &mut results);

        assert_eq!(results.len(), 1, "an unchecked event must leave a result");
        assert!(!results[0].pass, "an event no check covers may not pass");
    }

    /// A one-base small-variant truth event with the given REF/ALT.
    fn small_variant_event(reference: &[u8], alt: &[u8]) -> TruthEvent {
        TruthEvent {
            chrom: "chrA".to_string(),
            start: 250,
            end: 250 + reference.len() as u64,
            sv_type: "SNP".to_string(),
            expected_vaf: 0.5,
            gene: "unknown".to_string(),
            partner: None,
            ref_allele: Some(reference.to_vec()),
            alt_allele: Some(alt.to_vec()),
            ins_len: None,
            ins_alt: None,
            sim_number: None,
            census: CensusInfo::default(),
        }
    }

    #[test]
    fn test_allele_freq_it_cannot_measure_is_not_a_pass() {
        // Three returns inside check_allele_freq answered `pass: true` on
        // questions it had not asked: an indel, an MNV and an alt base that
        // is not one of ACGT. A result row was still pushed, so check_event's
        // "a check applies" gate counted all three as covered. N10 gave the
        // indel and the MNV a counting rule of their own, so what is left here
        // is every shape no rule fits -- the invariant is unchanged and the
        // cases it is asserted over moved. No BAM is read on any of these
        // paths, so the paths need not exist.
        for (reference, alt, what) in [
            (&b"AC"[..], &b"GTT"[..], "a complex allele"),
            (&b"A"[..], &b"CG"[..], "an insertion that rewrites its anchor"),
            (&b"AC"[..], &b"AGT"[..], "a deletion and an insertion at once"),
            (&b"T"[..], &b"N"[..], "a non-ACGT alt"),
            (&b"TG"[..], &b"AN"[..], "an MNV with a non-ACGT base"),
        ] {
            let event = small_variant_event(reference, alt);
            let r = check_allele_freq("/nonexistent/no.bam", "/nonexistent/no.fa", &event, &NearbyRecords::default(), 20)
                .expect("an unmeasurable allele fraction is a result, not an error");
            assert!(
                !r.pass,
                "{} is not an allele fraction this check measured, so it may not PASS \
                 (observed {:?})",
                what, r.observed
            );
        }
    }

    #[test]
    fn test_indel_and_mnv_allele_fractions_are_read_from_the_bam() {
        // The three shapes N10 gave a counting rule reached their verdict
        // before the BAM was opened, which is what made it a verdict about
        // nothing. An unreadable BAM must now be an error -- check_outcome
        // turns that into a failed row -- not a result reached without one.
        for (reference, alt, what) in [
            (&b"AC"[..], &b"A"[..], "a deletion"),
            (&b"A"[..], &b"ACGT"[..], "an insertion"),
            (&b"TG"[..], &b"AC"[..], "an MNV"),
        ] {
            let event = small_variant_event(reference, alt);
            let outcome = check_allele_freq("/nonexistent/no.bam", "/nonexistent/no.fa", &event, &NearbyRecords::default(), 20);
            assert!(
                outcome.is_err(),
                "{} must be measured from the BAM, not decided without opening one",
                what
            );
        }
    }

    #[test]
    fn test_allele_freq_below_the_depth_floor_is_not_a_pass() {
        // chrA is covered at depth 1-2 everywhere, below the 5-read floor the
        // check needs. "not enough data" was reported as a PASS.
        let (dir, fasta, cram) = two_contig_cram("low_depth_af");
        let event = small_variant_event(b"A", b"C");

        let r = check_allele_freq(&cram, &fasta, &event, &NearbyRecords::default(), 20).unwrap();

        assert!(
            !r.pass,
            "a depth the check calls too low to measure may not PASS (observed {:?})",
            r.observed
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    fn af_event(expected_vaf: f64) -> TruthEvent {
        TruthEvent {
            expected_vaf,
            ..small_variant_event(b"G", b"A")
        }
    }

    /// P(X = k) for k = 0..=n, X ~ Bin(n, p), summed in log space. Written
    /// here rather than borrowed from the code under test.
    fn binomial_pmf(n: u32, p: f64) -> Vec<f64> {
        if p >= 1.0 {
            let mut v = vec![0.0; n as usize + 1];
            v[n as usize] = 1.0;
            return v;
        }
        let mut ln = n as f64 * (1.0 - p).ln();
        let mut out = Vec::with_capacity(n as usize + 1);
        for k in 0..=n {
            out.push(ln.exp());
            ln += ((n - k) as f64).ln() - ((k + 1) as f64).ln() + p.ln() - (1.0 - p).ln();
        }
        out
    }

    #[test]
    fn test_allele_freq_never_passes_on_zero_alt_reads() {
        // REVIEW N14: 44 reads, not one carrying the alt, passed for every
        // SIM_VAF below 0.15 -- a correct spike-in and no spike-in at all
        // looked the same to the check.
        for vaf in [0.02, 0.05, 0.10, 0.1499, 0.15, 0.20, 0.50, 1.0] {
            let r = allele_freq_result(&af_event(vaf), 0, 44);
            assert!(!r.pass, "SIM_VAF {} passed on 0 of 44 alt reads ({:?})", vaf, r.observed);
        }
    }

    #[test]
    fn test_allele_freq_too_shallow_to_tell_says_so() {
        // 44 reads at VAF 0.05 expect ~2 alt reads: fewer than sequencing
        // errors can explain, so a correct run can't be told from none.
        let r = allele_freq_result(&af_event(0.05), 2, 44);
        assert!(!r.pass);
        assert!(r.observed.contains("too shallow"), "observed {:?}", r.observed);
    }

    #[test]
    fn test_allele_freq_grades_against_the_binomial_range() {
        // 36 of 100 at VAF 0.5 is outside the central 99% (P(X <= 36) is
        // 0.0033), though it sits within the old +-0.15 of 0.5.
        assert!(!allele_freq_result(&af_event(0.5), 36, 100).pass);
        assert!(allele_freq_result(&af_event(0.5), 50, 100).pass);
        // A hom truth takes a stray reference read, not many; and a perfect
        // hom run passes at any depth.
        assert!(allele_freq_result(&af_event(1.0), 39, 40).pass);
        assert!(!allele_freq_result(&af_event(1.0), 30, 40).pass);
        assert!(allele_freq_result(&af_event(1.0), 3000, 3000).pass);
    }

    #[test]
    fn test_allele_freq_needs_three_alt_reads_even_inside_the_range() {
        // 44 reads at VAF 0.19: P(X <= 2) is 0.0063, inside the central 99%,
        // so only the three-read floor keeps 2 alt reads from passing.
        assert!(!allele_freq_result(&af_event(0.19), 2, 44).pass);
        assert!(allele_freq_result(&af_event(0.19), 8, 44).pass);
    }

    #[test]
    fn test_allele_freq_passes_correct_runs_and_not_empty_ones() {
        // The N14 pass criteria locked in REVIEW.md, computed exactly:
        // wherever the check is evaluable, a correct run (x ~ Bin(n, p'))
        // passes >= 98% of the time and a run that planted nothing (only
        // 0.1% error reads) passes <= 0.5% of the time.
        let mut evaluable = Vec::new();
        for n in [20u32, 44, 100, 300, 1000, 3000] {
            for p in [0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0] {
                let event = af_event(p);
                let p_eff = p * (1.0 - 0.01) + (1.0 - p) * 0.001;
                let (correct_pmf, empty_pmf) = (binomial_pmf(n, p_eff), binomial_pmf(n, 0.001));
                // Counts neither run gives more than 1e-12 of the time move
                // either probability by < 3001e-12: skip grading them.
                let passes: Vec<bool> = (0..=n)
                    .map(|x| {
                        let i = x as usize;
                        (correct_pmf[i] > 1e-12 || empty_pmf[i] > 1e-12)
                            && allele_freq_result(&event, x, n).pass
                    })
                    .collect();
                if !passes.contains(&true) {
                    continue; // too shallow at this depth: nothing passes
                }
                evaluable.push((n, p));
                let pass_prob = |pmf: &[f64]| -> f64 {
                    pmf.iter().zip(&passes).filter(|(_, &ok)| ok).map(|(q, _)| q).sum()
                };
                let correct = pass_prob(&correct_pmf);
                let empty = pass_prob(&empty_pmf);
                assert!(correct >= 0.98, "n={} p={}: a correct run passes only {:.4}", n, p, correct);
                assert!(empty <= 0.005, "n={} p={}: an empty run passes {:.4}", n, p, empty);
            }
        }
        // A rule that never passes would satisfy both bounds vacuously.
        for cell in [(44, 0.5), (1000, 0.02), (3000, 0.01), (20, 1.0)] {
            assert!(evaluable.contains(&cell), "{:?} should be evaluable", cell);
        }
    }

    #[test]
    fn test_unrecognised_svtype_is_not_read_as_a_small_variant() {
        // `<CNV>` has no REF/ALT alleles to compare, but the SVTYPE fallback
        // routed it into the SNP arm, where ALT `<CNV>` is longer than one
        // base and took the indel exit -- one silent PASS per unknown type.
        let vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr20\t39100000\tcnv_1\tN\t<CNV>\t999\tPASS\tSVTYPE=CNV;END=39110000;SIM_VAF=0.500\tGT\t0/1
";
        let path =
            std::env::temp_dir().join(format!("spike_unknown_svtype_{}.vcf", std::process::id()));
        std::fs::write(&path, vcf).unwrap();
        let events = load_truth_events(path.to_str().unwrap()).unwrap();
        std::fs::remove_file(&path).ok();

        assert_eq!(events.len(), 1);
        assert_eq!(
            events[0].sv_type, "CNV",
            "an unrecognised SVTYPE must keep its own type, not be called a SNP"
        );
        assert!(events[0].ref_allele.is_none());
    }

    #[test]
    fn test_truth_record_ending_before_it_starts_is_rejected() {
        // END <= POS makes the event region empty; count_depth_in_region
        // returns a hard-coded 0.0 for it, so a DEL with SIM_VAF 0.9 PASSed
        // coverage_ratio (expected 0.10, observed 0.00) from a region no
        // query ever read.
        let vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr20\t39200000\tbad_1\tN\t<DEL>\t999\tPASS\tSVTYPE=DEL;END=39199000;SVLEN=-1000;SIM_VAF=0.900\tGT\t0/1
";
        let path =
            std::env::temp_dir().join(format!("spike_bad_end_{}.vcf", std::process::id()));
        std::fs::write(&path, vcf).unwrap();
        let err = load_truth_events(path.to_str().unwrap())
            .err()
            .map(|e| e.to_string())
            .unwrap_or_default();
        std::fs::remove_file(&path).ok();

        assert!(
            err.contains("END"),
            "a truth record ending before it starts must be refused; got {:?}",
            err
        );
    }

    #[test]
    fn test_ins_event_gets_a_real_check() {
        // spike writes SVTYPE=INS itself (truth.rs), so "no check applies"
        // made `spike validate` unable to exit 0 on spike's own output --
        // README step 4 broken for insertions. An INS must get a check that
        // reads the BAM, not a verdict reached without opening it.
        //
        // Since RF13 the counted row is `ins_planted`; `ins_reads` is advisory.
        let args = args_for("/nonexistent/no.bam", "/nonexistent/no.fa");
        let event = TruthEvent {
            sim_number: Some(1),
            ins_len: Some(4),
            ins_alt: Some(b"AGGGG".to_vec()),
            ..ins_event("chrA", 10_000)
        };
        let mut results: Vec<CheckResult> = Vec::new();

        check_event(&args, &event, &NearbyRecords::default(), &mut results);

        let counted: Vec<&CheckResult> = results.iter().filter(|r| !r.advisory).collect();
        assert_eq!(counted.len(), 1, "an INS must leave exactly one counted result");
        assert_eq!(
            counted[0].check_name, INS_PLANTED,
            "an INS must be checked for spike's reads carrying the inserted sequence"
        );
        assert!(
            counted[0].observed.starts_with("error"),
            "the row went to its data rather than decide without it: {}",
            counted[0].observed
        );
    }

    #[test]
    fn test_insertion_evidence_is_read_from_the_cigar() {
        // `std::io::Error` is not `Clone`, so each case rebuilds its ops.
        fn ops(spec: &[(Kind, usize)]) -> Vec<std::io::Result<Op>> {
            spec.iter().map(|&(k, n)| Ok(Op::new(k, n))).collect()
        }
        // 100M at 1000, then a 300 bp soft clip: the clip boundary is 1100.
        let clip_at_1100 = [(Kind::Match, 100), (Kind::SoftClip, 300)];
        assert!(cigar_shows_insertion_near(
            ops(&clip_at_1100).into_iter(),
            1000,
            1100,
            100,
            50
        ));
        // Same read, a breakpoint 500 bp away: not this event's evidence.
        assert!(!cigar_shows_insertion_near(
            ops(&clip_at_1100).into_iter(),
            1000,
            1600,
            100,
            50
        ));
        // Same read, clip too short to be the insertion.
        assert!(!cigar_shows_insertion_near(
            ops(&clip_at_1100).into_iter(),
            1000,
            1100,
            100,
            400
        ));

        // A leading clip's boundary is the alignment start.
        let leading = [(Kind::SoftClip, 120), (Kind::Match, 31)];
        assert!(cigar_shows_insertion_near(
            ops(&leading).into_iter(),
            1100,
            1100,
            100,
            50
        ));

        // An I operation after 50M, a 200 bp D and 10M from 1000 sits at
        // reference 1260: the deletion must advance the reference position.
        let with_del = [
            (Kind::Match, 50),
            (Kind::Deletion, 200),
            (Kind::Match, 10),
            (Kind::Insertion, 60),
            (Kind::Match, 91),
        ];
        assert!(cigar_shows_insertion_near(
            ops(&with_del).into_iter(),
            1000,
            1260,
            10,
            50
        ));
        assert!(!cigar_shows_insertion_near(
            ops(&with_del).into_iter(),
            1000,
            1050,
            10,
            50
        ));

        // A plain 151M read is evidence of nothing.
        assert!(!cigar_shows_insertion_near(
            ops(&[(Kind::Match, 151)]).into_iter(),
            1000,
            1100,
            100,
            50
        ));
    }

    #[test]
    fn test_a_short_insertions_evidence_must_be_an_insertion_not_a_clip() {
        // `ins_reads` takes min(SVLEN, 50) bases of inserted *or clipped*
        // sequence within 100 bp of POS as evidence, which is background once
        // SVLEN is small: on the merged HG002 chr20 slice a 3 bp threshold
        // finds 1, 3, 0, 0 and 1 reads at five positions where nothing was
        // planted, and 3 is over the 2-read pass mark -- an INS with SVLEN 3
        // PASSed on a position holding no insertion. A read that anchors both
        // sides of a short insertion writes an `I` operation; a 3 bp soft clip
        // is not evidence of one. The same five positions give 0, 0, 0, 0, 0
        // `I` operations.
        fn ops(spec: &[(Kind, usize)]) -> Vec<std::io::Result<Op>> {
            spec.iter().map(|&(k, n)| Ok(Op::new(k, n))).collect()
        }

        let clipped = [(Kind::Match, 100), (Kind::SoftClip, 3)];
        assert!(
            !cigar_shows_insertion_near(ops(&clipped).into_iter(), 1000, 1100, 100, 3),
            "a 3 bp soft clip is not evidence of a 3 bp insertion"
        );
        let inserted = [(Kind::Match, 100), (Kind::Insertion, 3), (Kind::Match, 48)];
        assert!(
            cigar_shows_insertion_near(ops(&inserted).into_iter(), 1000, 1100, 100, 3),
            "an I operation of the insertion's own length is"
        );
        // A clip is the only mark an insertion longer than a read can leave,
        // so at the 50 bp cap it still counts -- N8's 300 bp INS is read from
        // clipped reads and must keep passing.
        let long_clip = [(Kind::Match, 100), (Kind::SoftClip, 300)];
        assert!(
            cigar_shows_insertion_near(ops(&long_clip).into_iter(), 1000, 1100, 100, 50),
            "an insertion at or over the evidence cap is clipped, not inserted"
        );
    }

    #[test]
    fn test_header_marks_duplicates_from_the_pg_command_line() {
        use noodles::sam::header::record::value::{
            map::{program::tag, Program},
            Map,
        };

        // HG002's marker is ID:samtools.4 PN:samtools, and only its CL says
        // markdup -- so the ID and the name are not enough.
        let marked = noodles::sam::Header::builder()
            .add_program(
                "samtools.4",
                Map::<Program>::builder()
                    .insert(tag::NAME, "samtools")
                    .insert(tag::COMMAND_LINE, "/bin/samtools markdup -@ 4 - out.bam")
                    .build()
                    .unwrap(),
            )
            .build();
        let aligned_only = noodles::sam::Header::builder()
            .add_program(
                "bwa-mem2",
                Map::<Program>::builder()
                    .insert(tag::NAME, "bwa-mem2")
                    .insert(tag::COMMAND_LINE, "bwa-mem2 mem -t 8 ref.fa r1.fq r2.fq")
                    .build()
                    .unwrap(),
            )
            .build();

        // The same pipeline with markdup taken out: its command line still
        // mentions a path with "dedup" in the name, which is a file, not a
        // duplicate marker.
        let dedup_in_a_path = noodles::sam::Header::builder()
            .add_program(
                "samtools.4",
                Map::<Program>::builder()
                    .insert(tag::NAME, "samtools")
                    .insert(
                        tag::COMMAND_LINE,
                        "/bin/samtools sort -@ 4 - HG002.35x.bwamem2.dedup.grch38.bam",
                    )
                    .build()
                    .unwrap(),
            )
            .build();

        assert!(header_marks_duplicates(&marked));
        assert!(!header_marks_duplicates(&aligned_only));
        assert!(
            !header_marks_duplicates(&dedup_in_a_path),
            "a file name is not a duplicate marker"
        );
    }

    #[test]
    fn test_dup_rate_is_zero_when_the_header_marked_duplicates() {
        // A marked file whose sampled window happens to hold no duplicate is
        // not an unmarked file: 0% is what the pipeline decided, and telling
        // the user to "mark duplicates" is wrong advice about their BAM.
        let sample = GlobalSample {
            total: 41,
            dups: 0,
            mapq_sum: 2_460,
            insert_sizes: vec![300.0; 41],
            duplicates_marked: true,
        };

        let result = check_dup_rate(&sample);

        assert_eq!(result.observed, "0.0%");
        assert!(
            result.pass,
            "a marked file with no duplicate is 0%, and 0% passes"
        );
    }

    #[test]
    fn test_dup_rate_of_a_small_unmarked_sample_says_the_sample_is_too_small() {
        // With no markdup @PG to go on, 41 records without a duplicate flag
        // are no evidence of an unmarked file: at a 1% duplicate rate a
        // 41-record window holds none about two times in three.
        let sample = GlobalSample {
            total: 41,
            dups: 0,
            mapq_sum: 2_460,
            insert_sizes: vec![300.0; 41],
            duplicates_marked: false,
        };

        let result = check_dup_rate(&sample);

        assert_eq!(result.observed, "too few reads");
        assert!(!result.pass, "an unevaluable duplicate rate may not pass");
    }

    // --- N10: a small indel's and an MNV's allele fraction ---

    /// A one-contig CRAM holding a small deletion, a small insertion and an
    /// MNV, plus a site where nothing was planted.
    ///
    /// chrA's bases cycle `ACGT`, so the reference reads `ACG` at 5001-5003,
    /// `A` at 8001 and `AC` at 11001-11002 (1-based): the truth records are
    /// `ACG>A`, `A>ATTTT` and `AC>GT`. Six read1 records carry each variant in
    /// their CIGAR (`50M2D50M`, `50M4I46M`, and 100M with two substituted
    /// bases) and six span the same junction the way the reference does, so
    /// every site sits at an allele fraction of exactly 0.50. chrA:14001 has
    /// twelve reference reads and nothing planted. Each read1 has a plain
    /// read2 200 bp downstream, clear of the site it belongs to. Returns
    /// `(dir, fasta_path, cram_path)`; the caller removes `dir`.
    fn small_variant_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let seq = cycling_contig();

        // A full-length match over the reference's own bases.
        let plain = |start0: usize| -> TestRead {
            (
                start0,
                vec![(Kind::Match, TEST_READ_LEN)],
                seq[start0..start0 + TEST_READ_LEN].to_vec(),
            )
        };

        let mut pairs: Vec<(String, TestRead, TestRead)> = Vec::new();
        let mut push_pair = |name: String, start0: usize, ops: Vec<(Kind, usize)>, bases: Vec<u8>| {
            pairs.push((name, (start0, ops, bases), plain(start0 + 200)));
        };

        // Eight alt and eight reference pairs per site: at VAF 0.5 the
        // allele_freq check needs ~14 reads before it can tell a spike-in
        // from none (N14), so twelve would grade "too shallow".
        for i in 0..8usize {
            // `ACG` > `A` at chrA:5001: the two reference bases after the
            // anchor are gone, so the `D` operation starts at 0-based 5001.
            let s = 4951 + i;
            let (first_m, second_m) = (50 - i, 50 + i);
            let mut bases = seq[s..s + first_m].to_vec();
            bases.extend_from_slice(&seq[5003..5003 + second_m]);
            push_pair(
                format!("del_alt{}", i),
                s,
                vec![
                    (Kind::Match, first_m),
                    (Kind::Deletion, 2),
                    (Kind::Match, second_m),
                ],
                bases,
            );
            let (_, ops, bases) = plain(4941 + 3 * i);
            push_pair(format!("del_ref{}", i), 4941 + 3 * i, ops, bases);

            // `A` > `ATTTT` at chrA:8001: four bases inserted after the anchor.
            let s = 7951 + i;
            let (first_m, second_m) = (50 - i, 46 + i);
            let mut bases = seq[s..s + first_m].to_vec();
            bases.extend_from_slice(b"TTTT");
            bases.extend_from_slice(&seq[8001..8001 + second_m]);
            push_pair(
                format!("ins_alt{}", i),
                s,
                vec![
                    (Kind::Match, first_m),
                    (Kind::Insertion, 4),
                    (Kind::Match, second_m),
                ],
                bases,
            );
            let (_, ops, bases) = plain(7941 + 3 * i);
            push_pair(format!("ins_ref{}", i), 7941 + 3 * i, ops, bases);

            // `AC` > `GT` at chrA:11001: both bases swapped, no length change.
            let s = 10941 + 3 * i;
            let (_, ops, mut bases) = plain(s);
            bases[11_000 - s] = b'G';
            bases[11_001 - s] = b'T';
            push_pair(format!("mnv_alt{}", i), s, ops, bases);
            let (_, ops, bases) = plain(10945 + 3 * i);
            push_pair(format!("mnv_ref{}", i), 10945 + 3 * i, ops, bases);
        }
        // chrA:14001: sixteen reference reads and nothing planted -- deep
        // enough that "no alt read" is a measured FAIL, not "too shallow".
        for i in 0..16usize {
            let (_, ops, bases) = plain(13941 + 3 * i);
            push_pair(format!("clean_ref{}", i), 13941 + 3 * i, ops, bases);
        }

        pairs_cram(tag, &seq, &pairs)
    }

    /// Length of every read in the test CRAMs below.
    const TEST_READ_LEN: usize = 100;

    /// One read of a test pair: its 0-based start, CIGAR and bases.
    type TestRead = (usize, Vec<(Kind, usize)>, Vec<u8>);

    /// chrA for the test CRAMs: 20 kb of bases cycling `ACGT`, so 0-based
    /// position `p` holds `ACGT[p % 4]`.
    fn cycling_contig() -> Vec<u8> {
        (0..20_000).map(|i| b"ACGT"[i % 4]).collect()
    }

    // The five helpers below are the scaffolding every CRAM fixture in this
    // module needs and none of them is about: a scratch directory, the FASTA
    // and its `.fai`, the one-contig header, the writer plus its `.crai`, and
    // one 100 bp record on chrA. Only the read plan differs between fixtures,
    // so only that stays written out (F3, F4).

    /// A scratch directory of its own for one fixture, named after `tag` and
    /// this process, emptied first. The caller removes it.
    fn fixture_dir(tag: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "spike_test_validate_{}_{}",
            tag,
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    /// Write `seq` as a one-contig chrA FASTA at `dir/<basename>.fa`, with the
    /// `.fai` the CRAM reader needs beside it. Returns the FASTA's path.
    fn write_one_contig_fasta(
        dir: &std::path::Path,
        basename: &str,
        seq: &[u8],
    ) -> std::path::PathBuf {
        let mut fasta = String::from(">chrA\n");
        let offset = fasta.len();
        for chunk in seq.chunks(60) {
            fasta.push_str(std::str::from_utf8(chunk).unwrap());
            fasta.push('\n');
        }
        let fasta_path = dir.join(format!("{}.fa", basename));
        std::fs::write(&fasta_path, &fasta).unwrap();
        std::fs::write(
            dir.join(format!("{}.fa.fai", basename)),
            format!("chrA\t{}\t{}\t60\t61\n", seq.len(), offset),
        )
        .unwrap();
        fasta_path
    }

    /// A SAM header naming chrA alone, `len` bases long.
    fn one_contig_header(len: usize) -> noodles::sam::Header {
        noodles::sam::Header::builder()
            .add_reference_sequence(
                "chrA",
                noodles::sam::header::record::value::Map::<
                    noodles::sam::header::record::value::map::ReferenceSequence,
                >::new(std::num::NonZeroUsize::try_from(len).unwrap()),
            )
            .build()
    }

    /// Write `records` as `dir/<basename>.cram` against `fasta_path`, with its
    /// `.crai` beside it, and return the CRAM's path. The index is written
    /// from the CRAM itself, so a fixture cannot end up indexed as something
    /// it did not write.
    fn write_indexed_cram(
        dir: &std::path::Path,
        basename: &str,
        fasta_path: &std::path::Path,
        header: &noodles::sam::Header,
        records: &[noodles::cram::Record],
    ) -> std::path::PathBuf {
        let cram_path = dir.join(format!("{}.cram", basename));
        let repository =
            crate::extract::build_fasta_repository(fasta_path.to_str().unwrap()).unwrap();
        {
            let mut writer = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repository)
                .build_from_path(&cram_path)
                .unwrap();
            writer.write_header(header).unwrap();
            for rec in records {
                writer.write_record(header, rec.clone()).unwrap();
            }
            writer.try_finish(header).unwrap();
        }

        let index = noodles::cram::index(&cram_path).unwrap();
        let mut index_writer = noodles::cram::crai::io::Writer::new(
            std::fs::File::create(dir.join(format!("{}.cram.crai", basename))).unwrap(),
        );
        index_writer.write_index(&index).unwrap();
        index_writer.finish().unwrap();
        cram_path
    }

    /// One `TEST_READ_LEN`-long read of a one-contig chrA fixture: primary,
    /// mapped, unmarked, not a duplicate, single-end, its bases taken from
    /// `seq` at 0-based `start0` and its quality scores flat at 40. `tags`
    /// carries whatever else the fixture needs -- an `SA:Z` entry, or
    /// `Data::default()` for none.
    ///
    /// Two fixtures built this record independently and differed only in the
    /// MAPQ being a constant rather than a parameter and in the tags, so it is
    /// one function (F4).
    fn one_contig_record(
        seq: &[u8],
        name: &str,
        start0: usize,
        mapq: u8,
        tags: noodles::sam::alignment::record_buf::Data,
    ) -> noodles::cram::Record {
        noodles::cram::Record::builder()
            .set_bam_flags(noodles::sam::alignment::record::Flags::from(0u16))
            .set_flags(noodles::cram::record::Flags::QUALITY_SCORES_STORED_AS_ARRAY)
            .set_reference_sequence_id(0)
            .set_read_length(TEST_READ_LEN)
            .set_alignment_start(noodles::core::Position::new(start0 + 1).unwrap())
            .set_name(name)
            .set_mapping_quality(
                noodles::sam::alignment::record::MappingQuality::new(mapq).unwrap(),
            )
            .set_bases(noodles::sam::alignment::record_buf::Sequence::from(
                seq[start0..start0 + TEST_READ_LEN].to_vec(),
            ))
            .set_quality_scores(
                noodles::sam::alignment::record_buf::QualityScores::from(vec![
                    40u8;
                    TEST_READ_LEN
                ]),
            )
            .set_tags(tags)
            .build()
    }

    /// Write `pairs` -- `(name, read1, read2)` -- onto chrA as an indexed CRAM
    /// in a scratch directory of its own. Returns `(dir, fasta_path,
    /// cram_path)`; the caller removes `dir`.
    fn pairs_cram(
        tag: &str,
        seq: &[u8],
        pairs: &[(String, TestRead, TestRead)],
    ) -> (std::path::PathBuf, String, String) {
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

        let dir = fixture_dir(tag);
        let fasta_path = write_one_contig_fasta(&dir, "small_variant", seq);
        let header = one_contig_header(seq.len());

        // The CIGAR and the read bases are given explicitly: a small indel's
        // allele fraction is read off the CIGAR, so the fixture has to carry
        // real `I` and `D` operations rather than a full-length match.
        let record = |name: &str, read: &TestRead, mate: &TestRead, first: bool| {
            let (start0, ops, bases) = read;
            let flags = noodles::cram::record::Flags::QUALITY_SCORES_STORED_AS_ARRAY;
            let cigar: Cigar = ops.iter().map(|&(k, n)| Op::new(k, n)).collect();
            let sequence = Sequence::from(bases.clone());
            let quality_scores = QualityScores::from(vec![40u8; TEST_READ_LEN]);
            let features = noodles::cram::record::Features::from_cigar(
                flags,
                &cigar,
                &sequence,
                &quality_scores,
            );
            // Read 1 is the leftmost of every pair here, and a mate's
            // reference span is its M and D operations.
            let (r1_start, r2) = if first { (*start0, mate) } else { (mate.0, read) };
            let r2_span: usize = r2
                .1
                .iter()
                .filter(|(k, _)| matches!(k, Kind::Match | Kind::Deletion))
                .map(|&(_, n)| n)
                .sum();
            let template = (r2.0 + r2_span - r1_start) as i32;
            noodles::cram::Record::builder()
                .set_bam_flags(noodles::sam::alignment::record::Flags::from(if first {
                    0x63u16
                } else {
                    0x93u16
                }))
                .set_flags(flags)
                .set_reference_sequence_id(0)
                .set_read_length(TEST_READ_LEN)
                .set_alignment_start(noodles::core::Position::new(start0 + 1).unwrap())
                .set_name(name)
                .set_next_fragment_reference_sequence_id(0)
                .set_next_mate_alignment_start(noodles::core::Position::new(mate.0 + 1).unwrap())
                .set_template_size(if first { template } else { -template })
                .set_mapping_quality(
                    noodles::sam::alignment::record::MappingQuality::new(60).unwrap(),
                )
                .set_bases(sequence)
                .set_quality_scores(quality_scores)
                .set_features(features)
                .build()
        };

        let mut records: Vec<noodles::cram::Record> = Vec::new();
        for (name, read1, read2) in pairs {
            records.push(record(name, read1, read2, true));
            records.push(record(name, read2, read1, false));
        }

        let cram_path = write_indexed_cram(&dir, "small_variant", &fasta_path, &header, &records);

        (
            dir,
            fasta_path.to_str().unwrap().to_string(),
            cram_path.to_str().unwrap().to_string(),
        )
    }

    /// `small_variant_event` at a position of its own.
    fn small_variant_event_at(pos: u64, reference: &[u8], alt: &[u8]) -> TruthEvent {
        TruthEvent {
            start: pos,
            end: pos + reference.len() as u64,
            ..small_variant_event(reference, alt)
        }
    }

    #[test]
    fn test_small_indel_and_mnv_allele_fractions_are_measured() {
        // N9 stopped check_allele_freq calling an indel or an MNV a pass, but
        // it did not measure one, so a truth VCF from
        // `--event "snp:chr20:30000000:ACG:A"` carried a permanent
        // `allele_freq FAIL` -- spike's own round trip broken for a second
        // variant class, the same shape as N8's INS. Each site below is half
        // alt reads and half reference reads.
        let (dir, fasta, cram) = small_variant_cram("small_variant_af");

        for (pos, reference, alt, what) in [
            (5_000u64, &b"ACG"[..], &b"A"[..], "a 2 bp deletion"),
            (8_000, &b"A"[..], &b"ATTTT"[..], "a 4 bp insertion"),
            (11_000, &b"AC"[..], &b"GT"[..], "an MNV"),
        ] {
            let event = small_variant_event_at(pos, reference, alt);
            let r = check_allele_freq(&cram, &fasta, &event, &NearbyRecords::default(), 20).unwrap();
            assert_eq!(
                r.observed,
                "0.50",
                "{} is carried by half the reads at chrA:{}",
                what,
                pos + 1
            );
            assert!(r.pass, "{} at its own allele fraction must PASS", what);

            // The same record at a site where nothing was planted.
            let clean = small_variant_event_at(14_000, reference, alt);
            let r = check_allele_freq(&cram, &fasta, &clean, &NearbyRecords::default(), 20).unwrap();
            assert_eq!(
                r.observed, "0.00",
                "nothing was planted at chrA:14001, so no read carries {}",
                what
            );
            assert!(
                !r.pass,
                "{} must not pass at a site where nothing was planted",
                what
            );
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    // --- N16: one vote per fragment ---

    /// A CRAM whose mates overlap at four sites, read 2 starting 20 bp after
    /// read 1, so both reads of every pair cover the site.
    ///
    /// - chrA:16001 `A>T`: three pairs, both mates alt.
    /// - chrA:17001 `A>T`: five pairs with both mates alt, two whose read 1
    ///   is alt and read 2 is reference, and one with both mates reference.
    /// - chrA:18001 `ACG>A`: three pairs, both mates carrying the deletion.
    /// - chrA:19001 `ACG>A`: the same split as chrA:17001 -- five pairs
    ///   carrying the deletion in both mates, two carrying it in read 1 only
    ///   (read 2 spans the junction without it), one spanning it in both.
    fn overlapping_pairs_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let seq = cycling_contig();
        let snv = |start0: usize, site: usize, alt: bool| -> TestRead {
            let mut bases = seq[start0..start0 + TEST_READ_LEN].to_vec();
            if alt {
                bases[site - start0] = b'T';
            }
            (start0, vec![(Kind::Match, TEST_READ_LEN)], bases)
        };
        // `ACG` > `A` with the anchor at 0-based `site`: the `D` starts one
        // base after it.
        let del = |start0: usize, site: usize| -> TestRead {
            let first_m = site + 1 - start0;
            let second_m = TEST_READ_LEN - first_m;
            let mut bases = seq[start0..site + 1].to_vec();
            bases.extend_from_slice(&seq[site + 3..site + 3 + second_m]);
            (
                start0,
                vec![
                    (Kind::Match, first_m),
                    (Kind::Deletion, 2),
                    (Kind::Match, second_m),
                ],
                bases,
            )
        };

        let plain = |start0: usize| -> TestRead {
            (start0, vec![(Kind::Match, TEST_READ_LEN)], seq[start0..start0 + TEST_READ_LEN].to_vec())
        };

        let mut pairs: Vec<(String, TestRead, TestRead)> = Vec::new();
        for i in 0..8usize {
            let (s, m) = (15_960 + i, 15_980 + i);
            if i < 3 {
                pairs.push((format!("snv_once{}", i), snv(s, 16_000, true), snv(m, 16_000, true)));
                let (s, m) = (17_950 + i, 17_970 + i);
                pairs.push((format!("del_once{}", i), del(s, 18_000), del(m, 18_000)));
            }
            // Pairs 0-4 agree on alt, 5-6 split, 7 agrees on reference.
            let (read1_alt, read2_alt) = (i < 7, i < 5);
            let (s, m) = (s + 1000, m + 1000);
            pairs.push((
                format!("snv_split{}", i),
                snv(s, 17_000, read1_alt),
                snv(m, 17_000, read2_alt),
            ));
            let (s, m) = (18_950 + i, 18_970 + i);
            pairs.push((
                format!("del_split{}", i),
                if read1_alt { del(s, 19_000) } else { plain(s) },
                if read2_alt { del(m, 19_000) } else { plain(m) },
            ));
        }
        pairs_cram(tag, &seq, &pairs)
    }

    /// A hom truth record at `pos`, so every read that votes is expected alt.
    fn hom_event_at(pos: u64, reference: &[u8], alt: &[u8]) -> TruthEvent {
        TruthEvent {
            expected_vaf: 1.0,
            ..small_variant_event_at(pos, reference, alt)
        }
    }

    #[test]
    fn test_allele_freq_counts_an_overlapping_pair_once() {
        // The two mates of a pair read one molecule. Counted as two reads,
        // three pairs clear the five-read depth floor and the three-alt-read
        // floor that are meant to need five and three molecules (N16).
        let (dir, fasta, cram) = overlapping_pairs_cram("n16_once");
        for (pos, reference, what) in [(16_000u64, &b"A"[..], "an SNV"), (18_000, &b"ACG"[..], "a deletion")] {
            let alt = if reference.len() == 1 { &b"T"[..] } else { &reference[..1] };
            let r = check_allele_freq(&cram, &fasta, &hom_event_at(pos, reference, alt), &NearbyRecords::default(), 20).unwrap();
            assert_eq!(
                r.observed, "low depth (3)",
                "{}: three pairs are three votes, not six",
                what
            );
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_mates_that_disagree_give_their_pair_no_vote() {
        // Five pairs agree on alt, one on reference, and two have one mate
        // each way. A pair that says both is evidence for neither, the rule
        // MNVs already use: 5/6. Counting reads gave 12/16; letting a split
        // pair vote with its read 1 or its read 2 gives 7/8 or 5/8.
        let (dir, fasta, cram) = overlapping_pairs_cram("n16_split");
        for (pos, reference, what) in [(17_000u64, &b"A"[..], "an SNV"), (19_000, &b"ACG"[..], "a deletion")] {
            let alt = if reference.len() == 1 { &b"T"[..] } else { &reference[..1] };
            let r = check_allele_freq(&cram, &fasta, &hom_event_at(pos, reference, alt), &NearbyRecords::default(), 20).unwrap();
            assert_eq!(r.observed, "0.83", "{}: the two split pairs must not vote", what);
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    // --- N18: a read that stops at an indel cannot show it ---

    /// A read over the reference from `s` with one `kind` operation at
    /// reference position `q` (inserting `inserted`).
    fn read_with_indel(seq: &[u8], s: usize, q: usize, kind: Kind, len: usize, inserted: &[u8]) -> TestRead {
        let first = q - s;
        let mut bases = seq[s..q].to_vec();
        let second = match kind {
            Kind::Deletion => {
                let second = TEST_READ_LEN - first;
                bases.extend_from_slice(&seq[q + len..q + len + second]);
                second
            }
            _ => {
                let second = TEST_READ_LEN - first - len;
                bases.extend_from_slice(inserted);
                bases.extend_from_slice(&seq[q..q + second]);
                second
            }
        };
        (s, vec![(Kind::Match, first), (kind, len), (Kind::Match, second)], bases)
    }

    /// Two sites on `cycling_contig`, each with 8 pairs carrying the indel,
    /// 8 reference pairs reaching well past it, and 10 reference-aligned
    /// pairs that stop at it -- 6 ending 3 bases past the site, 4 starting 2
    /// bases before the anchor -- the way an aligner writes a carrier whose
    /// read ends there. Read 2 is always 200 bp downstream, clear of both.
    ///
    /// - chrA:5001 `ACG>A`, deleting 0-based 5001..5003.
    /// - chrA:8001 `A>ATTTT`, inserting `TTTT` in front of 8001.
    fn read_end_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let seq = cycling_contig();
        let plain = |s: usize| -> TestRead {
            (s, vec![(Kind::Match, TEST_READ_LEN)], seq[s..s + TEST_READ_LEN].to_vec())
        };
        let mut pairs: Vec<(String, TestRead, TestRead)> = Vec::new();
        let mut push = |name: String, read1: TestRead| {
            let mate = plain(read1.0 + 200);
            pairs.push((name, read1, mate));
        };
        for (site, far, kind, len, inserted, tag) in [
            (5_000usize, 5_003usize, Kind::Deletion, 2usize, &b""[..], "d"),
            (8_000, 8_001, Kind::Insertion, 4, &b"TTTT"[..], "i"),
        ] {
            for i in 0..8usize {
                push(format!("{}_alt{}", tag, i), read_with_indel(&seq, site - 49 + i, site + 1, kind, len, inserted));
                push(format!("{}_ref{}", tag, i), plain(site - 59 + i));
            }
            for i in 0..6usize {
                push(format!("{}_ends{}", tag, i), plain(far + 4 - TEST_READ_LEN));
            }
            for i in 0..4usize {
                push(format!("{}_starts{}", tag, i), plain(site - 2));
            }
        }
        pairs_cram(tag, &seq, &pairs)
    }

    #[test]
    fn test_a_read_that_stops_at_an_indel_does_not_vote_on_it() {
        // Near a read's end an aligner writes an indel as mismatches or a
        // clip, not a gap, so a carrier whose read stops there aligns as the
        // reference and voted `Spans` (N18). Reads that stop at the site are
        // dropped from both counts: 8 of 16, not 8 of 26.
        let (dir, fasta, cram) = read_end_cram("n18_ends");
        let seq = cycling_contig();
        for (pos, reference, alt, what) in [
            (5_000u64, seq[5_000..5_003].to_vec(), seq[5_000..5_001].to_vec(), "a deletion"),
            (8_000, seq[8_000..8_001].to_vec(), [&seq[8_000..8_001], &b"TTTT"[..]].concat(), "an insertion"),
        ] {
            let r = check_allele_freq(&cram, &fasta, &small_variant_event_at(pos, &reference, &alt), &NearbyRecords::default(), 20).unwrap();
            assert_eq!(r.observed, "0.50", "{}: reads stopping at it must not vote", what);
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    // --- N15: a read votes by the haplotype its bases spell ---

    /// Two sites on `cycling_contig`, 16 pairs each, read 2 always a plain
    /// 100M 200 bp downstream:
    ///
    /// - chrA:5001 `ACG>A` (deleting 0-based 5001..5003): 8 pairs carry it,
    ///   8 carry a 2 bp deletion at 5009 instead -- another sequence, since
    ///   no 2 bp deletion slides in `ACGT` repeated.
    /// - chrA:8002 `CG>C` (deleting 8002): 8 pairs carry the whole local
    ///   change a split truth record is one part of, a 3 bp deletion of
    ///   8002..8005; 8 are plain reference.
    fn haplotype_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let seq = cycling_contig();
        let plain = |s: usize| -> TestRead {
            (s, vec![(Kind::Match, TEST_READ_LEN)], seq[s..s + TEST_READ_LEN].to_vec())
        };
        let mut pairs: Vec<(String, TestRead, TestRead)> = Vec::new();
        let mut push = |name: String, read1: TestRead| {
            let mate = plain(read1.0 + 200);
            pairs.push((name, read1, mate));
        };
        for i in 0..8usize {
            push(format!("a_here{}", i), read_with_indel(&seq, 4_951 + i, 5_001, Kind::Deletion, 2, b""));
            push(format!("a_other{}", i), read_with_indel(&seq, 4_959 + i, 5_009, Kind::Deletion, 2, b""));
            push(format!("b_whole{}", i), read_with_indel(&seq, 7_952 + i, 8_002, Kind::Deletion, 3, b""));
            push(format!("b_ref{}", i), plain(7_941 + i));
        }
        pairs_cram(tag, &seq, &pairs)
    }

    #[test]
    fn test_a_read_votes_by_the_haplotype_its_bases_spell() {
        // Site A: a same-size gap 8 bp on is another variant, but the pad
        // rule counted it (N15): 16 of 16. Site B: the reads carry the whole
        // change the truth record is part of, written as one 3 bp gap, which
        // the pad rule could not see (N19): 0 carriers. Each read's bases are
        // closer to one haplotype, so each site is 8 of 16.
        let (dir, fasta, cram) = haplotype_cram("n15_haplotype");
        let seq = cycling_contig();
        for (pos, reference, alt, what) in [
            (5_000u64, &seq[5_000..5_003], &seq[5_000..5_001], "another deletion 8 bp on"),
            (8_001, &seq[8_001..8_003], &seq[8_001..8_002], "the whole change as one gap"),
        ] {
            let r = check_allele_freq(&cram, &fasta, &small_variant_event_at(pos, reference, alt), &NearbyRecords::default(), 20).unwrap();
            assert_eq!(r.observed, "0.50", "{}", what);
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// chrA:12001 `ACG>A` and chrA:12004 `T>TCG` on `cycling_contig`: one
    /// local change, `CGT` to `TCG` at 12001..12004, written by the truth set
    /// as two records the way GIAB often does (N19). 8 pairs carry it as the
    /// aligner writes it, three mismatches and no gap; 8 are plain reference.
    fn split_record_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let seq = cycling_contig();
        let plain = |s: usize| -> TestRead {
            (s, vec![(Kind::Match, TEST_READ_LEN)], seq[s..s + TEST_READ_LEN].to_vec())
        };
        let mut pairs: Vec<(String, TestRead, TestRead)> = Vec::new();
        for i in 0..8usize {
            let (s, ops, mut bases) = plain(11_951 + i);
            bases[12_001 - s..12_004 - s].copy_from_slice(b"TCG");
            pairs.push((format!("both{}", i), (s, ops, bases), plain(s + 200)));
            pairs.push((format!("ref{}", i), plain(11_941 + i), plain(12_141 + i)));
        }
        pairs_cram(tag, &seq, &pairs)
    }

    #[test]
    fn test_a_read_votes_by_the_nearby_truth_records_it_carries() {
        // Compared with each record alone, a read carrying both is two
        // edits from the reference and two from the record: a tie, no vote,
        // so 0 of 8. Compared with every combination of the two, it is the
        // pair, which holds each record: 8 of 16 for both.
        let (dir, fasta, cram) = split_record_cram("n15_split_record");
        let seq = cycling_contig();
        let events = vec![
            small_variant_event_at(12_000, &seq[12_000..12_003], &seq[12_000..12_001]),
            small_variant_event_at(12_003, &seq[12_003..12_004], b"TCG"),
        ];
        let nearby = NearbyRecords::new(&events);
        for event in &events {
            let r = check_allele_freq(&cram, &fasta, event, &nearby, 20).unwrap();
            assert_eq!(r.observed, "0.50", "{}", format_event_label(event));
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_a_crowded_site_is_compared_with_its_own_record_alone() {
        // Eleven other records in the window, one past the cap: the pair's
        // second half is left out with the rest, and its carriers tie again.
        let (dir, fasta, cram) = split_record_cram("n15_crowded");
        let seq = cycling_contig();
        let mut events = vec![
            small_variant_event_at(12_000, &seq[12_000..12_003], &seq[12_000..12_001]),
            small_variant_event_at(12_003, &seq[12_003..12_004], b"TCG"),
        ];
        for p in (11_991..12_000).chain(12_005..12_006) {
            events.push(small_variant_event_at(p as u64, &seq[p..p + 1], b"N"));
        }
        assert_eq!(events.len() - 1, MAX_NEARBY_EDITS + 1);
        let r = check_allele_freq(&cram, &fasta, &events[0], &NearbyRecords::new(&events), 20).unwrap();
        assert_eq!(r.observed, "too shallow (0.00 at 8 reads; needs about 14 reads)");
        let r = check_allele_freq(&cram, &fasta, &events[0], &NearbyRecords::new(&events[..11]), 20).unwrap();
        assert_eq!(r.observed, "0.50", "at the cap, every record is used");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_nearby_records_are_the_others_wholly_inside_the_window() {
        let seq = cycling_contig();
        let at = |p: usize, r: usize, a: &[u8]| small_variant_event_at(p as u64, &seq[p..p + r], a);
        let mut other_chrom = at(1_010, 1, b"T");
        other_chrom.chrom = "chrB".to_string();
        let events = vec![
            at(1_010, 3, &seq[1_010..1_011]), // the site
            at(1_000, 1, b"T"),               // first base: inside
            at(1_030, 1, b"T"),               // last base: inside
            at(999, 2, b"T"),                 // starts before the window
            at(1_030, 2, b"T"),               // runs past it
            at(1_010, 1, b"T"),               // same start, another record
            other_chrom,
            TruthEvent { ref_allele: None, alt_allele: None, ..at(1_020, 1, b"T") }, // no alleles
        ];
        let nearby = NearbyRecords::new(&events);
        let starts: Vec<u64> = nearby.inside(&events[0], 1_000, 1_030).iter().map(|e| e.0).collect();
        assert_eq!(starts, vec![1_000, 1_010, 1_030]);
    }

    #[test]
    fn test_overlapping_records_never_share_a_haplotype() {
        // Two alleles at one place are two haplotypes, never one.
        let window = cycling_contig();
        let site = (1_000, &window[1_000..1_001], &b"ATT"[..]);
        let same_place = (1_000, &window[1_000..1_001], &b"AGG"[..]);
        let elsewhere = (1_020, &window[1_020..1_021], &b"T"[..]);
        let count = |others: &[Edit]| site_haplotypes(&window, 0, 980, 1_030, site, others).unwrap().len();
        assert_eq!(count(&[]), 2);
        assert_eq!(count(&[elsewhere]), 4);
        assert_eq!(count(&[same_place]), 3, "the reference, the site, the other allele");
        assert_eq!(count(&[same_place, elsewhere]), 6);
    }

    #[test]
    fn test_another_truth_allele_at_the_site_counts_as_spanning_it() {
        // Alone, a same-size insertion of other bases is no nearer the
        // reference than the site's own. When the truth set lists it, a read
        // carrying it is nearest that allele, which does not hold the site.
        let window = cycling_contig();
        let site = (1_000, &window[1_000..1_001], &b"ATTTT"[..]);
        let other = (1_000, &window[1_000..1_001], &b"AGGCA"[..]);
        let mut read = window[980..1_001].to_vec();
        read.extend_from_slice(b"GGCA");
        read.extend_from_slice(&window[1_001..=1_030]);
        let alone = site_haplotypes(&window, 0, 980, 1_030, site, &[]).unwrap();
        let listed = site_haplotypes(&window, 0, 980, 1_030, site, &[other]).unwrap();
        assert_eq!(haplotype_vote(&read, &alone), None);
        assert_eq!(haplotype_vote(&read, &listed), Some(IndelVote::Spans));
    }

    #[test]
    fn test_indel_repeat_region_covers_the_whole_repeat() {
        // A window from 0: `GATG` + (AC)x10 + `GATC`, the repeat at 4..24
        // (a `C` in front of it would belong to it too, read as (CA)n). A
        // 2 bp deletion anchored at 5 (deleting 6..8) can slide over all of
        // it; so can an `AC` insertion anchored at 9.
        let mut w = b"GATG".to_vec();
        w.extend(std::iter::repeat_n(&b"AC"[..], 10).flatten());
        w.extend_from_slice(b"GATC");
        assert_eq!(indel_repeat_region(&w, 0, 5, Kind::Deletion, 2, b""), (4, 24));
        assert_eq!(indel_repeat_region(&w, 0, 9, Kind::Insertion, 2, b"AC"), (4, 24));
        // Off the repeat nothing slides: the deleted bases, the bare junction.
        assert_eq!(indel_repeat_region(&w, 0, 0, Kind::Deletion, 1, b""), (1, 2));
        assert_eq!(indel_repeat_region(&w, 0, 0, Kind::Insertion, 1, b"T"), (1, 1));
        // A poly-T insertion slides along the whole run, and the window
        // offset is honoured.
        let t = b"GGATTTTTTCGG";
        assert_eq!(indel_repeat_region(t, 100, 103, Kind::Insertion, 1, b"T"), (103, 109));
    }

    #[test]
    fn test_read_bases_over_needs_an_aligned_base_at_both_ends() {
        // A read that stops one base short of either end, or has a gap
        // there, cannot show the site whole and gives no bases (N18).
        let seq: Vec<u8> = (0..100).map(|i| b"ACGT"[i % 4]).collect();
        let plain = [Op::new(Kind::Match, 100)];
        assert_eq!(read_bases_over(&plain, &seq, 1_000, 1_010, 1_020), Some(&seq[10..=20]));
        assert_eq!(read_bases_over(&plain, &seq, 1_000, 999, 1_020), None, "starts one base late");
        assert_eq!(read_bases_over(&plain, &seq, 1_000, 1_010, 1_100), None, "ends one base early");
        // A deletion over the last position: no base is aligned to it.
        let gap = [Op::new(Kind::Match, 20), Op::new(Kind::Deletion, 2), Op::new(Kind::Match, 80)];
        assert_eq!(read_bases_over(&gap, &seq, 1_000, 1_010, 1_021), None);
        // An insertion inside the span is part of what the read shows.
        let ins = [Op::new(Kind::Match, 15), Op::new(Kind::Insertion, 3), Op::new(Kind::Match, 82)];
        assert_eq!(read_bases_over(&ins, &seq, 1_000, 1_010, 1_020), Some(&seq[10..=23]));
    }

    /// `ACG` > `A` at 0-based 1000 over `window` (from 0): the reference and
    /// truth haplotypes over 980..=1030, and the vote of a read from 950
    /// whose 2 bp `D` sits `shift` bases from 1001.
    fn shifted_deletion_vote(window: &[u8], shift: i64) -> Option<IndelVote> {
        const POS: u64 = 1_000;
        let site = (POS, &window[1_000..1_003], &window[1_000..1_001]);
        let haplotypes = site_haplotypes(window, 0, 980, 1_030, site, &[]).unwrap();
        let q = (POS as i64 + 1 + shift) as usize;
        let mut bases = window[950..q].to_vec();
        bases.extend_from_slice(&window[q + 2..1_052]);
        let ops = [Op::new(Kind::Match, q - 950), Op::new(Kind::Deletion, 2), Op::new(Kind::Match, 1_050 - q)];
        read_bases_over(&ops, &bases, 950, 980, 1_030).and_then(|b| haplotype_vote(b, &haplotypes))
    }

    #[test]
    fn test_a_deletion_written_anywhere_along_its_repeat_carries() {
        // Along a repeat the aligner may write one deletion anywhere (N13),
        // further than any fixed pad (N15). The read's bases are the same
        // whatever the spelling, so the vote is too.
        let window: Vec<u8> = (0..2_000).map(|i| b"AC"[i % 2]).collect();
        for shift in -18i64..=18 {
            assert_eq!(shifted_deletion_vote(&window, shift), Some(IndelVote::Carries), "shift {}", shift);
        }
    }

    #[test]
    fn test_a_deletion_off_a_repeat_votes_by_how_far_its_bases_are() {
        // In `ACGT` repeated no 2 bp deletion slides, so one written
        // elsewhere deletes other bases. One base off it is one edit from
        // the truth and two from the reference; two off is two from each,
        // a tie; three or more off it is another variant (N15).
        let window = cycling_contig();
        for shift in -12i64..=12 {
            let want = match shift.abs() {
                0 | 1 => Some(IndelVote::Carries),
                2 => None,
                _ => Some(IndelVote::Spans),
            };
            assert_eq!(shifted_deletion_vote(&window, shift), want, "shift {}", shift);
        }
    }

    #[test]
    fn test_an_insertion_votes_by_its_own_bases() {
        // `A` > `ATTTT` at 1000 in `ACGT` repeated.
        let window = cycling_contig();
        let site = (1_000, &window[1_000..1_001], &b"ATTTT"[..]);
        let haplotypes = site_haplotypes(&window, 0, 980, 1_030, site, &[]).unwrap();
        let read = |inserted: &[u8]| {
            let mut b = window[980..1_001].to_vec();
            b.extend_from_slice(inserted);
            b.extend_from_slice(&window[1_001..=1_030]);
            haplotype_vote(&b, &haplotypes)
        };
        assert_eq!(read(b"TTTT"), Some(IndelVote::Carries));
        assert_eq!(read(b"TTAT"), Some(IndelVote::Carries), "one base misread");
        assert_eq!(read(b""), Some(IndelVote::Spans), "the reference itself");
        // Another insertion of the same size is as many edits from the truth
        // as from the reference, or fewer, so it never counts as spanning --
        // the pad rule before this counted any same-size insertion.
        assert_eq!(read(b"GGCA"), None);
        assert_eq!(read(b"CGTA"), Some(IndelVote::Carries));
    }

    #[test]
    fn test_a_read_equally_far_from_both_does_not_vote() {
        // Two bases of the truth's four: two edits from either haplotype.
        let window = cycling_contig();
        let site = (1_000, &window[1_000..1_001], &b"ATTTT"[..]);
        let haplotypes = site_haplotypes(&window, 0, 980, 1_030, site, &[]).unwrap();
        let mut b = window[980..1_001].to_vec();
        b.extend_from_slice(b"TT");
        b.extend_from_slice(&window[1_001..=1_030]);
        assert_eq!(haplotype_vote(&b, &haplotypes), None);
        assert_eq!(levenshtein(b"kitten", b"SITTING"), 3, "the ruler, case aside");
    }

    #[test]
    fn test_ins_check_names_the_truth_records_own_pos() {
        // `load_truth_events` reads an INS as `start: vcf_pos` -- `start` is
        // already the insertion point the truth record names -- so
        // `event.start + 1` labelled the check with a base one past the POS
        // it was measured at.
        let (dir, fasta, cram) = small_variant_cram("ins_pos_label");
        let event = ins_event("chrA", 10_000);

        let r = check_ins_reads(&cram, &fasta, &event, 20).unwrap();

        assert!(
            r.expected.ends_with("chrA:10000"),
            "the check must name the truth record's own POS; got {:?}",
            r.expected
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
    // --- T1: the census numbers spike recorded, reported as advisory rows ---

    /// The advisory rows one truth record's INFO column leaves behind.
    fn census_rows(info: &str) -> Vec<CheckResult> {
        let event = TruthEvent {
            census: CensusInfo::from_info(info),
            ..ins_event("chrA", 10_000)
        };
        let mut results: Vec<CheckResult> = Vec::new();
        push_census_rows(&format_event_label(&event), &event, &mut results);
        results
    }

    fn row_names(rows: &[CheckResult]) -> Vec<&str> {
        rows.iter().map(|r| r.check_name.as_str()).collect()
    }

    /// One result with the pass/advisory shape a printing test needs.
    fn row(check_name: &str, pass: bool, advisory: bool) -> CheckResult {
        CheckResult {
            event_label: "DEL chrA:100-200 (unknown)".to_string(),
            check_name: check_name.to_string(),
            expected: "expected".to_string(),
            observed: "observed".to_string(),
            pass,
            advisory,
        }
    }

    fn text_report(results: &[CheckResult], strict: bool) -> String {
        let n_pass = results.iter().filter(|r| r.pass).count();
        let mut out: Vec<u8> = Vec::new();
        print_results_text(&mut out, results, results.len(), n_pass, strict).unwrap();
        String::from_utf8(out).unwrap()
    }

    /// The table line for one check, by name.
    ///
    /// **This is a substring match, not a field match.** `split_reads` is a
    /// prefix of `split_reads_each_end`, so on a report carrying both rows a
    /// bare `"split_reads"` returns whichever of the two comes first. The safe
    /// form for the pooled row is `"split_reads "` **with a trailing space**:
    /// the Check column is `{:<18}` and the format string adds one more, so
    /// every name shorter than 18 is followed by at least one space. A
    /// hand-built row list passed here is safe with the bare string only while
    /// it holds no `SPLIT_READS_EACH_END` row -- add one above the pooled row
    /// and this silently asserts about the wrong line.
    fn line_for<'a>(report: &'a str, check_name: &str) -> &'a str {
        report
            .lines()
            .find(|l| l.contains(check_name))
            .unwrap_or_else(|| panic!("no line for {} in:\n{}", check_name, report))
    }

    #[test]
    fn test_a_resistant_share_above_the_warn_threshold_is_an_advisory_fail() {
        let rows = census_rows("SVTYPE=INS;SIM_VAF=0.500;SIM_RESIST=0.500");

        assert_eq!(row_names(&rows), vec!["resistant"]);
        assert_eq!(rows[0].expected, "<=0.100");
        assert_eq!(rows[0].observed, "0.500");
        assert!(!rows[0].pass, "0.500 is above census::WARN_ABOVE");
        assert!(rows[0].advisory, "the census rows are advisory");
    }

    #[test]
    fn test_a_resistant_share_at_the_warn_threshold_is_an_advisory_pass() {
        // spike warns above the threshold, not at it, so the row passes at it.
        let rows = census_rows("SVTYPE=INS;SIM_RESIST=0.100");

        assert_eq!(row_names(&rows), vec!["resistant"]);
        assert_eq!(rows[0].observed, "0.100");
        assert!(rows[0].pass, "0.100 is census::WARN_ABOVE itself");
        assert!(rows[0].advisory);
    }

    #[test]
    fn test_a_depth_fold_above_the_warn_threshold_is_an_advisory_fail() {
        // 3.88 is the fold spike measured on the review's CR2 probe.
        let rows = census_rows("SVTYPE=DUP;SIM_DEPTH_FOLD=3.88");

        assert_eq!(row_names(&rows), vec!["depth_fold"]);
        assert_eq!(rows[0].expected, "<=1.50");
        assert_eq!(rows[0].observed, "3.88");
        assert!(!rows[0].pass, "3.88 is above census::DEPTH_FOLD_WARN_ABOVE");
        assert!(rows[0].advisory);
    }

    #[test]
    fn test_a_depth_fold_at_the_warn_threshold_is_an_advisory_pass() {
        let rows = census_rows("SVTYPE=DUP;SIM_DEPTH_FOLD=1.50");

        assert_eq!(row_names(&rows), vec!["depth_fold"]);
        assert_eq!(rows[0].observed, "1.50");
        assert!(rows[0].pass, "1.50 is census::DEPTH_FOLD_WARN_ABOVE itself");
    }

    #[test]
    fn test_the_two_expectations_are_the_thresholds_spike_warns_at() {
        // The printed expectation is derived from the constant, so moving the
        // constant moves the expectation with it.
        let rows = census_rows("SIM_RESIST=0.010;SIM_DEPTH_FOLD=1.00");

        assert_eq!(row_names(&rows), vec!["resistant", "depth_fold"]);
        assert_eq!(rows[0].expected, format!("<={:.3}", census::WARN_ABOVE));
        assert_eq!(
            rows[1].expected,
            format!("<={:.2}", census::DEPTH_FOLD_WARN_ABOVE)
        );
    }

    #[test]
    fn test_a_truth_record_without_the_census_fields_gets_no_row() {
        // An older spike's truth VCF carries neither field. That is nothing to
        // report, not a failure.
        assert!(census_rows("SVTYPE=INS;SVLEN=300;SIM_VAF=0.500").is_empty());
        // `.` is VCF's own missing value -- spike writes it when it has no
        // number -- so it is absent too, not a malformed number.
        assert!(census_rows("SIM_RESIST=.;SIM_DEPTH_FOLD=.").is_empty());
    }

    #[test]
    fn test_a_malformed_census_field_is_an_advisory_fail_not_a_dropped_row() {
        let rows = census_rows("SIM_RESIST=lots");

        assert_eq!(row_names(&rows), vec!["resistant"]);
        assert!(!rows[0].pass, "a number that cannot be read is not a pass");
        assert!(rows[0].advisory);
        assert!(
            rows[0].observed.contains("lots"),
            "the row must quote what it found; got {:?}",
            rows[0].observed
        );
    }

    #[test]
    fn test_load_truth_events_carries_the_census_numbers() {
        let vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr7\t55000\tsim_del_1\tN\t<DEL>\t999\tPASS\tSVTYPE=DEL;END=56000;SIM_VAF=0.500;SIM_RESIST=0.010;SIM_DEPTH_FOLD=3.88;SIM_GENE=EGFR\tGT\t0/1
chr7\t55201\tsim_var_1\tA\tT\t999\tPASS\tSIM_VAF=0.500;SIM_GENE=EGFR\tGT\t0/1
";
        let dir = std::env::temp_dir().join("spike_test_validate");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("truth_census.vcf");
        std::fs::write(&path, vcf).unwrap();

        let events = load_truth_events(path.to_str().unwrap()).unwrap();

        assert_eq!(events.len(), 2);
        assert_eq!(events[0].census.resist.as_deref(), Some("0.010"));
        assert_eq!(events[0].census.depth_fold.as_deref(), Some("3.88"));
        assert_eq!(
            events[1].census,
            CensusInfo::default(),
            "a record without the fields carries neither"
        );
    }

    #[test]
    fn test_an_advisory_row_alone_is_not_a_check_of_the_event() {
        // M11: an event no check covers must still report `event_checked
        // FAIL`, whether or not spike recorded a census for it. An advisory
        // row says what spike measured, not that this run verified anything.
        let args = args_for("/nonexistent/no.bam", "/nonexistent/no.fa");
        let event = TruthEvent {
            sv_type: "CNV".to_string(),
            census: CensusInfo::from_info("SIM_RESIST=0.010;SIM_DEPTH_FOLD=1.00"),
            ..ins_event("chrA", 10_000)
        };
        let mut results: Vec<CheckResult> = Vec::new();

        check_event(&args, &event, &NearbyRecords::default(), &mut results);

        assert_eq!(
            row_names(&results),
            vec!["resistant", "depth_fold", "event_checked"]
        );
        let fallback = results.last().unwrap();
        assert!(!fallback.pass, "an event no check covers may not pass");
        assert!(!fallback.advisory, "the fallback is a real failure");
    }

    #[test]
    fn test_advisory_failures_stay_out_of_the_exit_status_until_strict() {
        let results = vec![
            row("coverage_ratio", true, false),
            row("split_reads", false, false),
            row("resistant", false, true),
        ];

        assert_eq!(
            failure_message(&results, false).as_deref(),
            Some("1/2 validation checks failed"),
            "by default the exit status counts the real checks alone"
        );
        assert_eq!(
            failure_message(&results, true).as_deref(),
            Some("2/3 validation checks failed"),
            "--strict counts every row"
        );

        let only_advisory_fails = vec![row("coverage_ratio", true, false), row("resistant", false, true)];
        assert_eq!(
            failure_message(&only_advisory_fails, false),
            None,
            "an advisory failure alone leaves the run passing"
        );
        assert_eq!(
            failure_message(&only_advisory_fails, true).as_deref(),
            Some("1/2 validation checks failed")
        );
    }

    #[test]
    fn test_the_text_table_marks_the_advisory_rows_and_counts_them_apart() {
        let results = vec![
            row("coverage_ratio", true, false),
            row("resistant", false, true),
            row("depth_fold", true, true),
        ];

        let report = text_report(&results, false);

        assert!(
            line_for(&report, "coverage_ratio").ends_with(" PASS"),
            "a real check keeps its bare status; got {:?}",
            line_for(&report, "coverage_ratio")
        );
        assert!(
            line_for(&report, "resistant").ends_with(" FAIL (advisory)"),
            "got {:?}",
            line_for(&report, "resistant")
        );
        assert!(
            line_for(&report, "depth_fold").ends_with(" PASS (advisory)"),
            "got {:?}",
            line_for(&report, "depth_fold")
        );
        assert!(
            report.contains("\nResult: 2/3 PASS\n"),
            "the result line counts every row; got:\n{}",
            report
        );
        assert!(
            report.contains(
                "\nAdvisory: 2 checks, 1 PASS, 1 FAIL (not in the exit status; \
                 --strict includes them)\n"
            ),
            "got:\n{}",
            report
        );
    }

    #[test]
    fn test_the_any_mapq_row_does_not_move_the_status_column() {
        // The Check column is `{:<18}`, which *grows* past a longer name
        // rather than truncating it: measured, a 19-character name puts Status
        // at offset 98 where every other row has it at 97, with nothing else
        // to notice. `coverage_any_mapq` is 17, so it neither moves the column
        // nor eats its padding; this is what says so (F4).
        let report = text_report(
            &[
                row("coverage_ratio", true, false),
                row(COVERAGE_ANY_MAPQ, true, true),
            ],
            false,
        );
        let status_at = |name: &str| {
            let line = line_for(&report, name);
            line.find("PASS")
                .unwrap_or_else(|| panic!("no status on {:?}", line))
        };

        assert_eq!(
            status_at("coverage_any_mapq"),
            status_at("coverage_ratio"),
            "the any-MAPQ row must not move the Status column:\n{}",
            report
        );
        assert!(
            COVERAGE_ANY_MAPQ.len() < 18,
            "the name must still fit the 18-wide Check column with a pad left; \
             it is {} characters",
            COVERAGE_ANY_MAPQ.len()
        );
    }

    /// Every row a one-DEL truth VCF leaves on a **header-only** BAM: the five
    /// non-advisory checks all fail, and the only rows that pass are the two
    /// census rows, which are read back from the truth VCF and never touch the
    /// BAM at all.
    ///
    /// This is the report `scripts/validate_pipeline.sh`'s step-5 guard sees
    /// when `merge.sh` silently produces nothing, and the reason that guard
    /// cannot count `pass`: 2 of 9 rows passed and not one of them looked at a
    /// read.
    fn only_advisory_rows_pass() -> Vec<CheckResult> {
        vec![
            row(COVERAGE_RATIO, false, false),
            row(COVERAGE_ANY_MAPQ, false, true),
            row(SPLIT_READS, false, false),
            row(SPLIT_READS_EACH_END, false, true),
            row(RESISTANT, true, true),
            row(DEPTH_FOLD, true, true),
            row("insert_size", false, false),
            row("dup_rate", false, false),
            row("mean_mapq", false, false),
        ]
    }

    fn json_report(results: &[CheckResult], strict: bool) -> String {
        let n_pass = results.iter().filter(|r| r.pass).count();
        let mut out: Vec<u8> = Vec::new();
        print_results_json(&mut out, results, results.len(), n_pass, strict).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn test_a_json_report_whose_only_passing_rows_are_advisory_counts_no_pass() {
        // The Critical this pair of keys exists for: `passed == 0` is the only
        // signal `scripts/validate_pipeline.sh` has that the run measured
        // anything, and on this report `"pass"` is 2 while every check that
        // read the BAM failed. `counted_pass` is 0, so the guard fires again.
        let results = only_advisory_rows_pass();

        let report = json_report(&results, false);
        let summary = line_for(&report, "\"summary\"");

        assert!(
            summary.contains("\"counted_total\": 5")
                && summary.contains("\"counted_pass\": 0")
                && summary.contains("\"counted_fail\": 5")
                && summary.contains("\"strict\": false"),
            "the exit status counts the five non-advisory rows, none of which \
             passed; got:\n{}",
            summary
        );
        // The three original keys are over every row printed and must not have
        // moved: a parser written before the advisory rows existed still reads
        // what it always read.
        assert!(
            summary.contains("\"total\": 9")
                && summary.contains("\"pass\": 2")
                && summary.contains("\"fail\": 7"),
            "got:\n{}",
            summary
        );
        // And the trio is the exit status's own arithmetic, not a second copy
        // of it.
        assert_eq!(
            failure_message(&results, false).as_deref(),
            Some("5/5 validation checks failed"),
            "the run exits on the same five rows `counted_*` counts"
        );
    }

    #[test]
    fn test_strict_moves_the_counted_trio_and_nothing_else() {
        // Under --strict the exit status counts every row, so the trio counts
        // every row too -- the same filter, the same flag.
        let results = only_advisory_rows_pass();

        let report = json_report(&results, true);
        let summary = line_for(&report, "\"summary\"");

        assert!(
            summary.contains("\"counted_total\": 9")
                && summary.contains("\"counted_pass\": 2")
                && summary.contains("\"counted_fail\": 7")
                && summary.contains("\"strict\": true"),
            "got:\n{}",
            summary
        );
        assert!(
            summary.contains("\"total\": 9")
                && summary.contains("\"pass\": 2")
                && summary.contains("\"fail\": 7"),
            "the first three keys do not depend on --strict; got:\n{}",
            summary
        );
        assert_eq!(
            failure_message(&results, true).as_deref(),
            Some("7/9 validation checks failed"),
        );
    }

    #[test]
    fn test_a_report_with_no_advisory_row_counts_every_row_in_the_trio() {
        // The other end of it: nothing advisory, so `counted_*` and the three
        // original keys say the same thing and a consumer reading either gets
        // the same answer an older spike gave it.
        let results = vec![
            row(COVERAGE_RATIO, true, false),
            row(SPLIT_READS, false, false),
        ];

        let summary = {
            let report = json_report(&results, false);
            line_for(&report, "\"summary\"").to_string()
        };

        assert!(
            summary.contains("\"total\": 2")
                && summary.contains("\"pass\": 1")
                && summary.contains("\"fail\": 1")
                && summary.contains("\"counted_total\": 2")
                && summary.contains("\"counted_pass\": 1")
                && summary.contains("\"counted_fail\": 1"),
            "got:\n{}",
            summary
        );
    }

    #[test]
    fn test_the_advisory_summary_says_when_strict_counts_them() {
        let results = vec![row("resistant", false, true)];

        let report = text_report(&results, true);

        assert!(
            report.contains("\nAdvisory: 1 checks, 0 PASS, 1 FAIL (in the exit status: --strict)\n"),
            "got:\n{}",
            report
        );
    }

    #[test]
    fn test_a_report_with_no_advisory_row_prints_no_advisory_line() {
        // An older spike's truth VCF must print what master printed, to the
        // byte.
        let results = vec![row("coverage_ratio", true, false), row("split_reads", false, false)];

        for strict in [false, true] {
            let report = text_report(&results, strict);
            assert!(!report.contains("Advisory"), "got:\n{}", report);
            assert!(!report.contains("(advisory)"), "got:\n{}", report);
            assert!(line_for(&report, "split_reads").ends_with(" FAIL"));
        }
    }

    #[test]
    fn test_json_says_of_every_check_whether_it_is_advisory() {
        let results = vec![row("coverage_ratio", true, false), row("resistant", false, true)];
        let mut out: Vec<u8> = Vec::new();

        print_results_json(&mut out, &results, results.len(), 1, false).unwrap();

        let json = String::from_utf8(out).unwrap();
        assert!(
            json.contains("\"pass\": true, \"advisory\": false }"),
            "advisory follows pass; got:\n{}",
            json
        );
        assert!(
            json.contains("\"pass\": false, \"advisory\": true }"),
            "got:\n{}",
            json
        );
    }

    #[test]
    fn test_strict_is_off_unless_it_is_asked_for() {
        let argv = |extra: &[&str]| -> Vec<String> {
            let mut v: Vec<String> = [
                "spike",
                "validate",
                "--bam",
                "b.bam",
                "--truth",
                "t.vcf",
                "--reference",
                "r.fa",
            ]
            .iter()
            .map(|s| s.to_string())
            .collect();
            v.extend(extra.iter().map(|s| s.to_string()));
            v
        };

        assert!(
            !parse_validate_args_from(&argv(&[])).unwrap().strict,
            "--strict is off by default"
        );
        assert!(parse_validate_args_from(&argv(&["--strict"])).unwrap().strict);
        assert!(
            parse_validate_args_from(&argv(&["--strict=yes"])).is_err(),
            "--strict takes no value"
        );
    }

    #[test]
    fn test_an_advisory_check_that_cannot_run_stays_advisory() {
        // F1: `check_outcome` is the only factory for the "check runs" row,
        // and the advisory checks queued behind these two (coverage_ratio at
        // any MAPQ, split reads per breakpoint, INS sequence identity) all
        // read the BAM, so all of them will error through here. A
        // non-advisory error row would put such a failure into the exit
        // status by default, and would stand in for a check of the event.
        let r = check_outcome(
            "DEL chrA:100-200 (unknown)",
            "coverage_ratio",
            Err(anyhow::anyhow!("reference sequence not found: 20")),
            true,
        );

        assert!(r.advisory, "an advisory check that errored is still advisory");
        assert!(!r.pass, "a check that could not run is not a pass");
        assert_eq!(r.expected, "check runs");
        assert_eq!(
            failure_message(&[r], false),
            None,
            "an errored advisory check stays out of the default exit status"
        );
    }

    #[test]
    fn test_the_table_shows_what_a_malformed_census_field_held() {
        // F2: the Observed column is 14 characters, so the diagnostic has to
        // fit in it -- a row reading `unreadable:...` tells the reader
        // nothing. What the CheckResult holds is not enough; this pins the
        // table.
        let report = text_report(&census_rows("SIM_RESIST=lots"), false);
        let line = line_for(&report, "resistant");

        assert!(
            line.contains("bad: lots"),
            "the table must show the value it could not read; got {:?}",
            line
        );
        assert!(
            !line.contains("..."),
            "the diagnostic must not be truncated away; got {:?}",
            line
        );
        assert!(line.contains("FAIL (advisory)"), "got {:?}", line);
    }

    #[test]
    fn test_a_report_of_only_advisory_rows_never_counts_out_of_zero() {
        // F5: by default `failure_message` counts the non-advisory rows
        // against the non-advisory total, which is zero here. The `n_fail ==
        // 0` early return is what keeps `x/0` off the screen.
        let results = vec![row("resistant", false, true)];

        assert_eq!(
            failure_message(&results, false),
            None,
            "nothing the default counts failed, so there is no message at all"
        );
        assert_eq!(
            failure_message(&results, true).as_deref(),
            Some("1/1 validation checks failed"),
            "--strict counts the advisory row, against the advisory row"
        );
    }

    /// A one-contig CRAM whose deletion window holds nothing but MAPQ 0
    /// reads, and a second window that holds nothing at all.
    ///
    /// chrA is 20 kb of cycling bases. Every read is 100 bp, primary, mapped,
    /// unmarked and unpaired, so the depth over a window is exactly the bases
    /// its reads put there.
    ///
    /// - `chrA:[10000,11000)` -- the *hidden* deletion: 10 reads at **MAPQ 0**
    ///   (depth 1.0), the reads a spike-in could not touch because they never
    ///   entered its donor pool. Its flanks `[9000,10000)` and `[11000,12000)`
    ///   carry 20 reads each at MAPQ 60 (depth 2.0).
    /// - `chrA:[14000,15000)` -- the *clean* deletion: no read at any MAPQ,
    ///   with the same two flanks at depth 2.0 (`[13000,14000)` and
    ///   `[15000,16000)`).
    ///
    /// So with `--flank 1000`, a floor of 20 reads the hidden window as ratio
    /// 0.00 and a floor of 0 reads it as 0.50, while the clean window reads
    /// 0.00 at either floor.
    ///
    /// Returns `(dir, fasta_path, cram_path)`; the caller removes `dir`.
    fn mapq_hidden_deletion_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let seq = cycling_contig();
        let dir = fixture_dir(tag);
        let fasta_path = write_one_contig_fasta(&dir, "mapq_hidden", &seq);
        let header = one_contig_header(seq.len());

        // 20 reads at MAPQ 60 over a 1000 bp flank is depth 2.0; 10 inside a
        // 1000 bp window is depth 1.0. Every read lies wholly inside the
        // window it is meant to cover, so no depth leaks across a boundary.
        let mut plan: Vec<(usize, u8)> = Vec::new();
        for flank_start in [9_000usize, 11_100, 13_000, 15_100] {
            for i in 0..20usize {
                plan.push((flank_start + i * 40, 60));
            }
        }
        for i in 0..10usize {
            plan.push((10_050 + i * 80, 0));
        }
        plan.sort_unstable();

        // No tags: this fixture is about MAPQ and depth, and the split-read
        // scan is the only thing that reads a tag.
        let records: Vec<noodles::cram::Record> = plan
            .iter()
            .enumerate()
            .map(|(i, &(start0, mapq))| {
                one_contig_record(
                    &seq,
                    &format!("read{}", i),
                    start0,
                    mapq,
                    noodles::sam::alignment::record_buf::Data::default(),
                )
            })
            .collect();
        let cram_path = write_indexed_cram(&dir, "mapq_hidden", &fasta_path, &header, &records);

        (
            dir,
            fasta_path.to_str().unwrap().to_string(),
            cram_path.to_str().unwrap().to_string(),
        )
    }

    /// An AF=1 DEL over [start, end) on chrA, with `--flank 1000`.
    fn hom_del_args(cram: &str, fasta: &str) -> ValidateArgs {
        ValidateArgs {
            flank_bp: 1_000,
            ..args_for(cram, fasta)
        }
    }

    fn hom_del(start: u64, end: u64) -> TruthEvent {
        TruthEvent {
            expected_vaf: 1.0,
            ..del_event("chrA", start, end)
        }
    }

    /// The row with this name, or a panic naming what was there instead.
    fn named<'a>(rows: &'a [CheckResult], check_name: &str) -> &'a CheckResult {
        rows.iter()
            .find(|r| r.check_name == check_name)
            .unwrap_or_else(|| panic!("no {} row in {:?}", check_name, row_names(rows)))
    }

    #[test]
    fn test_the_any_mapq_row_sees_the_depth_the_mapq_floor_hides() {
        // CR4 on one input: half the donor pairs at MAPQ 0 leave half the
        // original depth inside an AF=1 deletion, and `coverage_ratio` cannot
        // see it, because the check applies the same `--min-mapq 20` floor
        // that kept those reads out of the donor pool. The two rows must
        // disagree here -- that disagreement is the row's whole purpose.
        let (dir, fasta, cram) = mapq_hidden_deletion_cram("any_mapq_hidden");
        let args = hom_del_args(&cram, &fasta);
        let mut results: Vec<CheckResult> = Vec::new();

        check_event(&args, &hom_del(10_000, 11_000), &NearbyRecords::default(), &mut results);
        let _ = std::fs::remove_dir_all(&dir);

        let floored = named(&results, "coverage_ratio");
        assert_eq!(floored.expected, "0.00");
        assert_eq!(floored.observed, "0.00", "the MAPQ 0 reads are invisible here");
        assert!(floored.pass, "the deletion looks perfect at --min-mapq 20");
        assert!(!floored.advisory, "coverage_ratio stays a real check");

        let any_mapq = named(&results, "coverage_any_mapq");
        assert_eq!(any_mapq.expected, "0.00", "the same expectation as the floored row");
        assert_eq!(any_mapq.observed, "0.50", "half the flank depth is still there");
        assert!(!any_mapq.pass, "0.50 is outside the 0.30 tolerance of 0.00");
        assert!(any_mapq.advisory, "the any-MAPQ row is always advisory");
    }

    #[test]
    fn test_a_deletion_with_nothing_left_inside_it_passes_at_every_mapq() {
        // The other side of the same fixture: a window emptied at every MAPQ
        // must not be reported as a problem by the new row. Without this, a
        // row that always failed would look just as "right" as one that
        // reads the depth.
        let (dir, fasta, cram) = mapq_hidden_deletion_cram("any_mapq_clean");
        let args = hom_del_args(&cram, &fasta);
        let mut results: Vec<CheckResult> = Vec::new();

        check_event(&args, &hom_del(14_000, 15_000), &NearbyRecords::default(), &mut results);
        let _ = std::fs::remove_dir_all(&dir);

        for name in ["coverage_ratio", "coverage_any_mapq"] {
            let r = named(&results, name);
            assert_eq!(r.observed, "0.00", "{} read {}", name, r.observed);
            assert!(r.pass, "{} must pass a deletion that really is empty", name);
        }
    }

    #[test]
    fn test_min_mapq_zero_still_prints_both_coverage_rows() {
        // At `--min-mapq 0` the two rows are identical by construction. Both
        // are still produced: a row that disappeared when the floor reached
        // its own would make the report's shape depend on a flag.
        let (dir, fasta, cram) = mapq_hidden_deletion_cram("any_mapq_floor_zero");
        let args = ValidateArgs {
            min_mapq: 0,
            ..hom_del_args(&cram, &fasta)
        };
        let mut results: Vec<CheckResult> = Vec::new();

        check_event(&args, &hom_del(10_000, 11_000), &NearbyRecords::default(), &mut results);
        let _ = std::fs::remove_dir_all(&dir);

        let floored = named(&results, "coverage_ratio");
        let any_mapq = named(&results, "coverage_any_mapq");
        assert_eq!(floored.observed, "0.50", "the floor is gone from the real check too");
        assert_eq!(any_mapq.observed, floored.observed);
        assert_eq!(any_mapq.expected, floored.expected);
        assert_eq!(any_mapq.pass, floored.pass);
        assert!(!floored.advisory, "the real check is still in the exit status");
        assert!(any_mapq.advisory, "its twin is still out of it");
    }

    #[test]
    fn test_the_any_mapq_row_is_advisory_in_the_table_and_in_the_json() {
        // No BAM here: a check that cannot run must stay advisory too, or a
        // missing file would start failing runs that exit 0 today (T1's F1).
        let args = args_for("/nonexistent/no.bam", "/nonexistent/no.fa");
        let mut results: Vec<CheckResult> = Vec::new();
        check_event(&args, &hom_del(10_000, 11_000), &NearbyRecords::default(), &mut results);

        let any_mapq = named(&results, "coverage_any_mapq");
        assert!(any_mapq.advisory, "an errored any-MAPQ row is still advisory");
        assert_eq!(any_mapq.expected, "check runs");
        assert!(!any_mapq.pass);
        assert_eq!(
            failure_message(&results, false).as_deref(),
            Some("2/2 validation checks failed"),
            "the advisory row is in neither the count nor the total"
        );
        assert_eq!(
            failure_message(&results, true).as_deref(),
            // Five rows on a DEL since RF14: this row, `coverage_ratio`,
            // `del_planted`, and the advisory `split_reads` and
            // `split_reads_each_end`. The default count above is the half that
            // matters: still 2/2, now `coverage_ratio` and `del_planted`.
            Some("5/5 validation checks failed"),
            "--strict counts it with the rest"
        );

        let report = text_report(&results, false);
        assert!(
            line_for(&report, "coverage_any_mapq").ends_with("FAIL (advisory)"),
            "the table must mark the row; got {:?}",
            line_for(&report, "coverage_any_mapq")
        );
        assert!(
            line_for(&report, "coverage_ratio").ends_with(" FAIL"),
            "the real check's row is untouched; got {:?}",
            line_for(&report, "coverage_ratio")
        );

        let mut json: Vec<u8> = Vec::new();
        print_results_json(&mut json, &results, results.len(), 0, false).unwrap();
        let json = String::from_utf8(json).unwrap();
        assert!(
            line_for(&json, "coverage_any_mapq").contains("\"advisory\": true"),
            "{}",
            json
        );
        assert!(
            line_for(&json, "\"check\": \"coverage_ratio\"").contains("\"advisory\": false"),
            "{}",
            json
        );
    }

    #[test]
    fn test_the_any_mapq_row_covers_the_types_coverage_ratio_covers_and_no_others() {
        // The row is the coverage check at another floor, so it applies
        // exactly where a ratio of event depth to flank depth means
        // something: a DEL and a DUP, and nothing else.
        let args = args_for("/nonexistent/no.bam", "/nonexistent/no.fa");
        let event_of = |sv_type: &str| TruthEvent {
            sv_type: sv_type.to_string(),
            ..hom_del(10_000, 11_000)
        };

        for sv_type in ["DEL", "DUP"] {
            let mut results: Vec<CheckResult> = Vec::new();
            check_event(&args, &event_of(sv_type), &NearbyRecords::default(), &mut results);
            let names = row_names(&results);
            assert!(
                names.contains(&"coverage_ratio") && names.contains(&"coverage_any_mapq"),
                "{} must get both coverage rows; got {:?}",
                sv_type,
                names
            );
        }

        for sv_type in ["INV", "BND", "INS", "SNP", "CNV"] {
            let mut results: Vec<CheckResult> = Vec::new();
            check_event(&args, &event_of(sv_type), &NearbyRecords::default(), &mut results);
            let names = row_names(&results);
            assert!(
                !names.contains(&"coverage_any_mapq"),
                "{} has no expected coverage ratio, so it gets no any-MAPQ row; got {:?}",
                sv_type,
                names
            );
            assert!(
                !names.contains(&"coverage_ratio"),
                "{} gets no coverage_ratio either -- the two rows must cover the same types; got {:?}",
                sv_type,
                names
            );
        }
    }

    #[test]
    fn test_an_event_whose_rows_are_all_advisory_still_reports_event_checked() {
        // M11 again, now with three advisory rows in play: the "a check
        // applies" fallback counts the non-advisory rows only, so an event no
        // real check covers must still report `event_checked FAIL` however
        // many advisory rows sit beside it.
        let args = args_for("/nonexistent/no.bam", "/nonexistent/no.fa");
        let event = TruthEvent {
            sv_type: "CNV".to_string(),
            census: CensusInfo::from_info("SIM_RESIST=0.500;SIM_DEPTH_FOLD=3.88"),
            ..hom_del(10_000, 11_000)
        };
        let mut results: Vec<CheckResult> = Vec::new();

        check_event(&args, &event, &NearbyRecords::default(), &mut results);

        assert_eq!(row_names(&results), vec!["resistant", "depth_fold", "event_checked"]);
        let checked = named(&results, "event_checked");
        assert!(!checked.pass, "an event no check covers may not pass");
        assert!(!checked.advisory, "it is the row that fails the run");
    }

    // --- NF5: the split-read minimum at each breakpoint, not pooled ---

    /// A one-contig chrA CRAM holding three deletion junctions that differ
    /// only in **where** the split reads sit. Every record is 100 bp, primary,
    /// mapped, unmarked, MAPQ 60, and carries an `SA:Z` entry; `split_reads`
    /// keys on that tag, so a record without one is invisible to both rows.
    /// The three junctions, with the 500 bp windows the check reads
    /// (1-based, inclusive):
    ///
    /// * **one-sided**, `chrA:10000-12000`: three reads in `[9501,10500]`
    ///   carrying `SA:Z:chrA,12001,...` and three in `[11501,12500]` carrying
    ///   `SA:Z:chrA,1501,...` -- present at the partner breakpoint, but
    ///   pointing somewhere else entirely, so nothing joins back. 3 here, 0
    ///   there.
    /// * **two-sided**, `chrA:5000-7000`: three reads in `[4501,5500]`
    ///   pointing at chrA:7001 and three in `[6501,7500]` pointing back at
    ///   chrA:5001. 3 here, 3 there.
    /// * **one name at both ends**, `chrA:15000-17000`: a single read *name*
    ///   with one record in each window, each pointing at the other -- a pair
    ///   whose two mates straddle the junction, which is what a real BAM looks
    ///   like. 1 here, 1 there, and **1** pooled, not 2.
    /// * **within the pad**, `chrA:2000-2500` (span 500 == `pad`): three reads
    ///   straddling the left breakpoint, naming the right one, and **nothing
    ///   at the right breakpoint at all**. The two windows overlap and each
    ///   window's SA test reaches the other breakpoint, so the same three
    ///   records answer both: 3 here, 3 there. The documented blind zone.
    /// * **past the pad**, `chrA:18000-18501` (span 501, one base more): the
    ///   same layout one base past the boundary, where the SA test at the
    ///   partner window no longer reaches `here`. 3 here, 0 there.
    ///
    /// The five pairs of windows are disjoint, so no junction's reads are
    /// visible to another's query. Returns `(dir, fasta_path, cram_path)`; the
    /// caller removes `dir`.
    fn one_sided_split_read_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let seq = cycling_contig();
        let dir = fixture_dir(tag);
        let fasta_path = write_one_contig_fasta(&dir, "split_ends", &seq);
        let header = one_contig_header(seq.len());

        // (name, 0-based start, the 1-based position its SA:Z names).
        let mut plan: Vec<(String, usize, u64)> = Vec::new();
        for i in 0..3usize {
            // One-sided: the left window points at the right breakpoint, the
            // right window points at chrA:1501 -- neither near chrA:10001.
            plan.push((format!("one_sided_left{}", i), 9_900 + i * 50, 12_001));
            plan.push((format!("one_sided_right{}", i), 11_900 + i * 50, 1_501));
            // Two-sided: each end points at the other.
            plan.push((format!("two_sided_left{}", i), 4_900 + i * 50, 7_001));
            plan.push((format!("two_sided_right{}", i), 6_900 + i * 50, 5_001));
        }
        // The same name in both windows of the third junction.
        plan.push(("both_ends".to_string(), 14_950, 17_001));
        plan.push(("both_ends".to_string(), 16_950, 15_001));
        // The blind zone and its boundary: three reads straddling the left
        // breakpoint, naming the right one, and nothing at the right
        // breakpoint. The two junctions differ only in their span -- 500 and
        // 501 -- so what separates them is `pad` alone.
        for i in 0..3usize {
            plan.push((format!("within_pad{}", i), 1_950 + i * 25, 2_501));
            plan.push((format!("past_pad{}", i), 17_950 + i * 25, 18_502));
        }
        plan.sort_by_key(|&(_, start0, _)| start0);

        // One `SA:Z` entry naming a 1-based position on chrA. The scan reads
        // its first two fields and nothing else, so the rest is filler.
        let sa_tag = |sa_pos: u64| -> noodles::sam::alignment::record_buf::Data {
            [(
                noodles::sam::alignment::record::data::field::Tag::new(b'S', b'A'),
                noodles::sam::alignment::record_buf::data::field::Value::from(
                    format!("chrA,{},+,50M50S,60,0;", sa_pos).as_str(),
                ),
            )]
            .into_iter()
            .collect()
        };

        let records: Vec<noodles::cram::Record> = plan
            .iter()
            .map(|(name, start0, sa_pos)| {
                one_contig_record(&seq, name, *start0, 60, sa_tag(*sa_pos))
            })
            .collect();
        let cram_path = write_indexed_cram(&dir, "split_ends", &fasta_path, &header, &records);

        (
            dir,
            fasta_path.to_str().unwrap().to_string(),
            cram_path.to_str().unwrap().to_string(),
        )
    }

    /// The two split-read rows for one junction of the fixture above.
    fn split_rows(cram: &str, fasta: &str, start: u64, end: u64) -> Vec<CheckResult> {
        let args = args_for(cram, fasta);
        let mut results: Vec<CheckResult> = Vec::new();
        check_event(&args, &del_event("chrA", start, end), &NearbyRecords::default(), &mut results);
        results
    }

    #[test]
    fn test_the_each_end_row_fails_evidence_that_all_sits_at_one_breakpoint() {
        // NF5 on one input: three reads join the left breakpoint to the right
        // and nothing joins back, so the pooled row's two-name minimum is met
        // at one end alone and it PASSes -- while its own `expected` reads as
        // though each end contributed. The per-end row is the one that says
        // so. The two rows disagreeing here is the whole point of the row.
        let (dir, fasta, cram) = one_sided_split_read_cram("split_one_sided");
        let results = split_rows(&cram, &fasta, 10_000, 12_000);
        let _ = std::fs::remove_dir_all(&dir);

        let pooled = named(&results, "split_reads");
        assert_eq!(pooled.expected, ">=2 joining chrA:12001");
        assert_eq!(pooled.observed, "3", "three distinct names, all at one end");
        assert!(pooled.pass, "the pooled row passes one-sided evidence");
        // A DEL here: since RF14 its `split_reads` is advisory and `del_planted`
        // decides. `test_the_del_rows_are_del_planted_counted_and_split_reads_advisory`
        // pins that a DUP's still counts.
        assert!(pooled.advisory, "a DEL's split_reads is advisory since RF14");

        let each_end = named(&results, "split_reads_each_end");
        assert_eq!(each_end.observed, "3/0", "nothing joins back from the partner");
        assert!(!each_end.pass, "0 at the partner breakpoint is below the minimum");
        assert!(each_end.advisory, "the per-end row is always advisory");

        // Both rows carry `check_event`'s own label, not one recomputed inside
        // the check: `coverage_ratio`'s row is built from that same local, so
        // the three must agree. Nothing pinned this before -- a wrong label
        // threaded into `split_read_rows` passed the whole suite (measured).
        let label = &named(&results, "coverage_ratio").event_label;
        assert_eq!(&pooled.event_label, label, "the pooled row's Event column");
        assert_eq!(&each_end.event_label, label, "the per-end row's Event column");
    }

    #[test]
    fn test_the_each_end_row_passes_a_junction_with_reads_at_both_breakpoints() {
        // The other side of the same fixture: a junction whose reads really do
        // sit at both ends must pass. Without this, a row that always failed
        // would look as "right" as one that reads each end.
        let (dir, fasta, cram) = one_sided_split_read_cram("split_two_sided");
        let results = split_rows(&cram, &fasta, 5_000, 7_000);
        let _ = std::fs::remove_dir_all(&dir);

        let pooled = named(&results, "split_reads");
        assert_eq!(pooled.observed, "6", "three names at each end, all distinct");
        assert!(pooled.pass);

        let each_end = named(&results, "split_reads_each_end");
        assert_eq!(each_end.observed, "3/3");
        assert!(each_end.pass, "three at each end clears a minimum of two");
        assert!(each_end.advisory, "a passing per-end row is advisory too");
    }

    #[test]
    fn test_a_read_name_at_both_breakpoints_counts_once_in_the_pooled_row() {
        // The pooled count is the **union** of the two name sets, not their
        // sum: one read name with a record in each window is one read, and it
        // votes once. Summing would read 2 and PASS on the evidence of a
        // single read.
        let (dir, fasta, cram) = one_sided_split_read_cram("split_both_ends");
        let results = split_rows(&cram, &fasta, 15_000, 17_000);
        let _ = std::fs::remove_dir_all(&dir);

        let pooled = named(&results, "split_reads");
        assert_eq!(pooled.observed, "1", "one name at both ends is one read");
        assert!(!pooled.pass, "one read is below the pooled minimum of two");

        let each_end = named(&results, "split_reads_each_end");
        assert_eq!(each_end.observed, "1/1", "each end sees that one name");
        assert!(!each_end.pass, "one at each end is below the minimum");
    }

    #[test]
    fn test_the_each_end_expectation_is_the_minimum_split_reads_asks_for() {
        // Both rows say what they ask for by naming the constant, not by
        // repeating "2" in a string: a changed minimum must change both
        // printed expectations with it.
        let (dir, fasta, cram) = one_sided_split_read_cram("split_expectation");
        let results = split_rows(&cram, &fasta, 10_000, 12_000);
        let _ = std::fs::remove_dir_all(&dir);

        assert_eq!(
            named(&results, "split_reads_each_end").expected,
            format!(">={} at each end", MIN_SPLIT_READS)
        );
        assert_eq!(
            named(&results, "split_reads").expected,
            format!(">={} joining chrA:12001", MIN_SPLIT_READS),
            "the pooled row's expected string did not move"
        );
    }

    #[test]
    fn test_the_each_end_row_is_advisory_in_the_table_and_in_the_json() {
        // No BAM here: a check that cannot run must stay advisory too, or a
        // missing file would start failing runs that exit 0 today (T1's F1).
        // Both rows come out of one pair of queries, so one unreadable BAM is
        // one error reported twice -- once in the exit status, once beside it.
        // A DUP, whose pooled row still counts: a DEL's is advisory since RF14.
        let args = args_for("/nonexistent/no.bam", "/nonexistent/no.fa");
        let mut results: Vec<CheckResult> = Vec::new();
        let dup = TruthEvent {
            sv_type: "DUP".to_string(),
            ..del_event("chrA", 10_000, 12_000)
        };
        check_event(&args, &dup, &NearbyRecords::default(), &mut results);

        let each_end = named(&results, "split_reads_each_end");
        assert!(each_end.advisory, "an errored per-end row is still advisory");
        assert_eq!(each_end.expected, "check runs");
        assert!(!each_end.pass);
        assert!(
            !named(&results, "split_reads").advisory,
            "the pooled row's error stays in the exit status"
        );
        assert_eq!(
            failure_message(&results, false).as_deref(),
            Some("2/2 validation checks failed"),
            "neither advisory row is in the count or the total"
        );
        assert_eq!(
            failure_message(&results, true).as_deref(),
            Some("4/4 validation checks failed"),
            "--strict counts them with the rest"
        );

        let report = text_report(&results, false);
        assert!(
            line_for(&report, "split_reads_each_end").ends_with("FAIL (advisory)"),
            "the table must mark the row; got {:?}",
            line_for(&report, "split_reads_each_end")
        );
        // A trailing space, so the name does not match the longer row's line:
        // `{:<18}` pads `split_reads` to 18 and the format adds one more.
        assert!(
            line_for(&report, "split_reads ").ends_with(" FAIL"),
            "the real check's row is untouched; got {:?}",
            line_for(&report, "split_reads ")
        );

        let mut json: Vec<u8> = Vec::new();
        print_results_json(&mut json, &results, results.len(), 0, false).unwrap();
        let json = String::from_utf8(json).unwrap();
        assert!(
            line_for(&json, "\"check\": \"split_reads_each_end\"").contains("\"advisory\": true"),
            "{}",
            json
        );
        assert!(
            line_for(&json, "\"check\": \"split_reads\",").contains("\"advisory\": false"),
            "{}",
            json
        );
    }

    #[test]
    fn test_the_each_end_row_covers_the_types_split_reads_covers_and_no_others() {
        // The row is the same evidence at a stricter bar, so it applies
        // exactly where a junction between two breakpoints means something:
        // DEL, DUP, INV and BND, and nothing else. INS is out for the reason
        // `split_reads` is: an inserted sequence has no second reference
        // breakpoint.
        let args = args_for("/nonexistent/no.bam", "/nonexistent/no.fa");
        let event_of = |sv_type: &str| TruthEvent {
            sv_type: sv_type.to_string(),
            ..del_event("chrA", 10_000, 12_000)
        };

        for sv_type in ["DEL", "DUP", "INV", "BND"] {
            let mut results: Vec<CheckResult> = Vec::new();
            check_event(&args, &event_of(sv_type), &NearbyRecords::default(), &mut results);
            let names = row_names(&results);
            assert!(
                names.contains(&"split_reads") && names.contains(&"split_reads_each_end"),
                "{} must get both split-read rows; got {:?}",
                sv_type,
                names
            );
        }

        for sv_type in ["INS", "SNP", "CNV"] {
            let mut results: Vec<CheckResult> = Vec::new();
            check_event(&args, &event_of(sv_type), &NearbyRecords::default(), &mut results);
            let names = row_names(&results);
            assert!(
                !names.contains(&"split_reads_each_end"),
                "{} has no second breakpoint, so it gets no per-end row; got {:?}",
                sv_type,
                names
            );
            assert!(
                !names.contains(&"split_reads"),
                "{} gets no split_reads either -- the two rows must cover the same types; got {:?}",
                sv_type,
                names
            );
        }
    }

    #[test]
    fn test_the_each_end_row_is_the_only_row_its_width_moves() {
        // The Check column is `{:<18}`, which *grows* past a longer name
        // rather than truncating it: T2's review measured a 19-character name
        // putting Status at 98 where every other row has it at 97.
        // `split_reads_each_end` is 20 characters, so this row's own later
        // columns sit two characters right -- the cost the plan weighed and
        // accepted rather than rename the row or widen the column. What has to
        // hold is that it costs nothing anywhere else: `{:<18}` pads short
        // names, so the rows either side keep the offset they had.
        let report = text_report(
            &[
                row(COVERAGE_RATIO, true, false),
                row(SPLIT_READS, true, false),
                row(SPLIT_READS_EACH_END, true, true),
                row("ins_reads", true, false),
            ],
            false,
        );
        let status_at = |name: &str| {
            let line = line_for(&report, name);
            line.find("PASS")
                .unwrap_or_else(|| panic!("no status on {:?}", line))
        };

        assert_eq!(
            SPLIT_READS_EACH_END.len(),
            20,
            "the accepted overflow of the 18-wide Check column is two characters"
        );
        assert_eq!(
            status_at(SPLIT_READS_EACH_END),
            status_at(COVERAGE_RATIO) + 2,
            "the per-end row's Status sits exactly two characters right:\n{}",
            report
        );
        // A trailing space on the pooled name, so it cannot match the longer
        // row's line.
        for name in [COVERAGE_RATIO, "split_reads ", "ins_reads"] {
            assert_eq!(
                status_at(name),
                status_at(COVERAGE_RATIO),
                "{} must keep the offset every other row has:\n{}",
                name,
                report
            );
        }
    }

    #[test]
    fn test_the_each_end_row_is_never_a_check_of_the_event() {
        // M11 with the per-end row in play. Two halves, because no event type
        // can receive the per-end row *without* also receiving the pooled
        // non-advisory row -- so "only advisory rows" cannot be built out of
        // this row, and what is checkable is that the fallback still counts
        // the non-advisory rows only, and that the row stays out of the
        // default exit status on its own.
        let args = args_for("/nonexistent/no.bam", "/nonexistent/no.fa");

        // An event no real check covers still reports `event_checked FAIL`,
        // and gets neither split-read row.
        let cnv = TruthEvent {
            sv_type: "CNV".to_string(),
            census: CensusInfo::from_info("SIM_RESIST=0.500;SIM_DEPTH_FOLD=3.88"),
            ..del_event("chrA", 10_000, 12_000)
        };
        let mut results: Vec<CheckResult> = Vec::new();
        check_event(&args, &cnv, &NearbyRecords::default(), &mut results);
        assert_eq!(row_names(&results), vec!["resistant", "depth_fold", "event_checked"]);
        let checked = named(&results, "event_checked");
        assert!(!checked.pass, "an event no check covers may not pass");
        assert!(!checked.advisory, "it is the row that fails the run");

        // A DEL is covered, so the fallback must not fire beside four rows of
        // which two are advisory.
        let mut results: Vec<CheckResult> = Vec::new();
        check_event(&args, &del_event("chrA", 10_000, 12_000), &NearbyRecords::default(), &mut results);
        assert!(
            !row_names(&results).contains(&"event_checked"),
            "a DEL is checked; got {:?}",
            row_names(&results)
        );

        // And the per-end row alone never fails a run by default.
        let i = results
            .iter()
            .position(|r| r.check_name == "split_reads_each_end")
            .unwrap();
        assert_eq!(failure_message(&results[i..=i], false), None);
    }

    #[test]
    fn test_the_each_end_rows_documented_blind_zone_for_spans_within_the_pad() {
        // A **documented blind zone**, pinned as the behaviour it is and not
        // as behaviour anyone wants. Both breakpoint windows are `pad` = 500
        // bp wide, and the SA test at each window asks whether the entry lands
        // within 500 bp of the *other* breakpoint. So for an event whose two
        // breakpoints are 500 bp apart or less, one read sitting only at the
        // left breakpoint whose `SA:Z` names the right one satisfies both
        // halves: it is inside both windows, and its SA is within 500 bp of
        // both. `at_here` and `at_partner` become the same set, the two counts
        // are the same number, and the row degenerates into the pooled row --
        // for that whole size class it cannot tell one-sided evidence from
        // two-sided. `scripts/validate_pipeline.sh` sets MIN_DEL_SIZE=500, so
        // the pipeline's smallest admissible deletion sits on the boundary.
        //
        // The locked plan mandated reusing `pad`, so this is **pinned, not
        // endorsed**: README.md's "what it does not" list for the row says the
        // same thing in prose. Both junctions below have three reads at the
        // left breakpoint naming the right one and nothing at the right
        // breakpoint at all, and differ only in their span -- END - POS of 500
        // against 501 -- so the boundary itself is pinned and a change to
        // `pad` reddens this test.
        let (dir, fasta, cram) = one_sided_split_read_cram("split_blind_zone");
        let within_pad = split_rows(&cram, &fasta, 2_000, 2_500);
        let past_pad = split_rows(&cram, &fasta, 18_000, 18_501);
        let _ = std::fs::remove_dir_all(&dir);

        // END - POS == 500: the two windows and the two SA tests coincide, so
        // one-sided evidence reads as two-sided and the advisory row PASSes
        // exactly where it exists to FAIL.
        let each_end = named(&within_pad, "split_reads_each_end");
        assert_eq!(
            each_end.observed, "3/3",
            "the blind zone: the same three reads are counted at both ends"
        );
        assert!(
            each_end.pass,
            "the blind zone: a span of 500 passes on one-sided evidence"
        );
        let pooled = named(&within_pad, "split_reads");
        assert_eq!(pooled.observed, "3", "the union is those same three names");
        assert!(pooled.pass, "the pooled row is what this row degenerates to");

        // END - POS == 501, one base past the boundary: the same read layout
        // now reads as the one-sided evidence it is.
        let each_end = named(&past_pad, "split_reads_each_end");
        assert_eq!(
            each_end.observed, "3/0",
            "one base past the pad, nothing joins back"
        );
        assert!(!each_end.pass, "past the blind zone the row FAILs");
        let pooled = named(&past_pad, "split_reads");
        assert_eq!(pooled.observed, "3", "the pooled row did not move either");
        assert!(pooled.pass, "the pooled row passes on either side of it");
    }

    // --- T6: the inserted bases the truth record names, found in the reads ---

    /// Deterministic pseudo-random ACGT, `len` bases from `seed`.
    ///
    /// chrA cycles `ACGT` with period 4, so every 31-mer in it is a rotation of
    /// `ACGT` repeated; an insertion has to bring a sequence of its own or the
    /// reference guard would fire on it. A fixed seed keeps the fixture
    /// identical run to run, and each test that depends on a k-mer being absent
    /// asserts that rather than trusting the generator.
    fn pseudo_random_bases(seed: u64, len: usize) -> Vec<u8> {
        let mut x = seed | 1;
        (0..len)
            .map(|_| {
                x = x
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                b"ACGT"[((x >> 33) % 4) as usize]
            })
            .collect()
    }

    /// The reverse complement of `bases`, as a value.
    fn revcomp_of(bases: &[u8]) -> Vec<u8> {
        let mut out = bases.to_vec();
        crate::extract::reverse_complement(&mut out);
        out
    }

    /// The 60 bp insertion the reads at [`INS_CARRIED_POS`] carry whole.
    fn carried_insertion() -> Vec<u8> {
        pseudo_random_bases(11, 60)
    }

    /// A 2 kb insertion: longer than any read, so only its first and last
    /// bases are ever inside one.
    fn long_insertion() -> Vec<u8> {
        pseudo_random_bases(13, 2_000)
    }

    /// An 8 bp insertion -- below [`MIN_INS_KMER_LEN`], so the row cannot
    /// grade it.
    fn short_insertion() -> Vec<u8> {
        pseudo_random_bases(17, 8)
    }

    /// The bases of `long_insertion` a read at [`INS_LONG_MIDDLE_POS`] carries:
    /// 60 out of the middle, reachable by no read of a real run.
    fn long_middle() -> Vec<u8> {
        long_insertion()[985..1_045].to_vec()
    }

    // The five insertion loci of `inserted_sequence_cram`, 3 kb apart, so that
    // no locus's +/-150 read window or `ins_reads`'s +/-100 CIGAR window
    // reaches another's reads.
    const INS_CARRIED_POS: u64 = 5_000;
    const INS_REVCOMP_POS: u64 = 8_000;
    const INS_LONG_ENDS_POS: u64 = 11_000;
    const INS_LONG_MIDDLE_POS: u64 = 14_000;
    const INS_SHORT_POS: u64 = 17_000;

    /// How far before the insertion point the inverted copy of the ALT's first
    /// probe k-mer sits in [`inverted_repeat_contig`]: inside
    /// [`INS_KMER_REF_PAD`], so the reference guard reads it, and inside
    /// [`INS_SEQUENCE_PAD`] as well, so an *unedited* read over the insertion
    /// point carries it in its own bases.
    const INVERTED_REPEAT_OFFSET: usize = 60;

    /// One single-end read of a one-contig chrA fixture with its CIGAR and its
    /// bases given explicitly.
    ///
    /// [`one_contig_record`] slices its bases out of the reference and sets no
    /// features, which is right for a read that matches; a read carrying an
    /// insertion cannot be built that way, because CRAM stores a match as a
    /// reference-relative feature and would hand back the reference's bases.
    /// The bases therefore arrive through `Features::from_cigar` together with
    /// the CIGAR that explains them, as `pairs_cram`'s paired records do.
    fn one_contig_spliced_record(
        name: &str,
        start0: usize,
        ops: &[(Kind, usize)],
        bases: &[u8],
    ) -> noodles::cram::Record {
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

        let flags = noodles::cram::record::Flags::QUALITY_SCORES_STORED_AS_ARRAY;
        let cigar: Cigar = ops.iter().map(|&(k, n)| Op::new(k, n)).collect();
        let sequence = Sequence::from(bases.to_vec());
        let quality_scores = QualityScores::from(vec![40u8; bases.len()]);
        let features = noodles::cram::record::Features::from_cigar(
            flags,
            &cigar,
            &sequence,
            &quality_scores,
        );
        noodles::cram::Record::builder()
            .set_bam_flags(noodles::sam::alignment::record::Flags::from(0u16))
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_read_length(bases.len())
            .set_alignment_start(noodles::core::Position::new(start0 + 1).unwrap())
            .set_name(name)
            .set_mapping_quality(
                noodles::sam::alignment::record::MappingQuality::new(60).unwrap(),
            )
            .set_bases(sequence)
            .set_quality_scores(quality_scores)
            .set_features(features)
            .build()
    }

    /// One indexed CRAM on the 20 kb cycling chrA holding five insertion loci
    /// that differ only in *which bases* their reads carry, so that
    /// `ins_reads` passes at every one of them and `ins_sequence` can only
    /// differ by reading the bases:
    ///
    /// - [`INS_CARRIED_POS`]: three reads carrying the whole 60 bp
    ///   `carried_insertion` as an `I` operation.
    /// - [`INS_REVCOMP_POS`]: two reads carrying its **reverse complement**.
    /// - [`INS_LONG_ENDS_POS`]: one read soft-clipping the **first** 60 bases
    ///   of the 2 kb `long_insertion`, one the **last** 60.
    /// - [`INS_LONG_MIDDLE_POS`]: two reads soft-clipping 60 bases out of that
    ///   insertion's **middle**.
    /// - [`INS_SHORT_POS`]: two reads carrying the 8 bp `short_insertion`.
    ///
    /// Every record is 100 bp, primary, mapped, unmarked, single-end and MAPQ
    /// 60, so `usable_alignment` admits all of them at the default
    /// `--min-mapq 20` and only the bases separate the loci.
    fn inserted_sequence_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let seq = cycling_contig();
        let dir = fixture_dir(tag);
        let fasta_path = write_one_contig_fasta(&dir, "ins_sequence", &seq);
        let header = one_contig_header(seq.len());

        // A read carrying `inserted` whole, anchored on both sides: an `I`
        // operation at `pos`, which is what an aligner writes for an insertion
        // shorter than a read.
        let spanning = |name: &str, pos: u64, inserted: &[u8]| {
            let flank = (TEST_READ_LEN - inserted.len()) / 2;
            let start0 = (pos as usize) - flank;
            let mut bases = seq[start0..start0 + flank].to_vec();
            bases.extend_from_slice(inserted);
            bases.extend_from_slice(&seq[start0 + flank..start0 + 2 * flank]);
            one_contig_spliced_record(
                name,
                start0,
                &[
                    (Kind::Match, flank),
                    (Kind::Insertion, inserted.len()),
                    (Kind::Match, flank),
                ],
                &bases,
            )
        };
        // A read anchored to the *left* of `pos` whose trailing 60 bases are
        // soft-clipped: what an aligner writes for an insertion too long to
        // anchor both sides of.
        let clipped_right = |name: &str, pos: u64, carried: &[u8]| {
            let start0 = (pos as usize) - 40;
            let mut bases = seq[start0..start0 + 40].to_vec();
            bases.extend_from_slice(carried);
            one_contig_spliced_record(
                name,
                start0,
                &[(Kind::Match, 40), (Kind::SoftClip, carried.len())],
                &bases,
            )
        };
        // The mirror image: anchored to the *right* of `pos`, leading bases
        // clipped.
        let clipped_left = |name: &str, pos: u64, carried: &[u8]| {
            let start0 = pos as usize;
            let mut bases = carried.to_vec();
            bases.extend_from_slice(&seq[start0..start0 + 40]);
            one_contig_spliced_record(
                name,
                start0,
                &[(Kind::SoftClip, carried.len()), (Kind::Match, 40)],
                &bases,
            )
        };

        let carried = carried_insertion();
        let long = long_insertion();
        let mut records: Vec<noodles::cram::Record> = Vec::new();
        // In alignment order, so the CRAM is coordinate-sorted.
        for i in 0..3u64 {
            records.push(spanning(
                &format!("carried{}", i),
                INS_CARRIED_POS + i,
                &carried,
            ));
        }
        for i in 0..2u64 {
            records.push(spanning(
                &format!("revcomp{}", i),
                INS_REVCOMP_POS + i,
                &revcomp_of(&carried),
            ));
        }
        records.push(clipped_right("long_first", INS_LONG_ENDS_POS, &long[..60]));
        records.push(clipped_left(
            "long_last",
            INS_LONG_ENDS_POS,
            &long[long.len() - 60..],
        ));
        for i in 0..2u64 {
            records.push(clipped_right(
                &format!("long_middle{}", i),
                INS_LONG_MIDDLE_POS + i,
                &long_middle(),
            ));
        }
        for i in 0..2u64 {
            records.push(spanning(
                &format!("short{}", i),
                INS_SHORT_POS + i,
                &short_insertion(),
            ));
        }

        let cram_path = write_indexed_cram(&dir, "ins_sequence", &fasta_path, &header, &records);
        (
            dir,
            fasta_path.to_str().unwrap().to_string(),
            cram_path.to_str().unwrap().to_string(),
        )
    }

    /// chrA with the **reverse complement** of `carried_insertion()`'s first
    /// [`INS_KMER_LEN`] bases written in [`INVERTED_REPEAT_OFFSET`] bases
    /// before [`INS_CARRIED_POS`] -- an inverted repeat of the insertion
    /// itself, which is the case the reference guard's reverse-complement half
    /// exists to catch. The probe's **forward** orientation is nowhere in the
    /// contig; only its reverse complement is, and only within
    /// [`INS_KMER_REF_PAD`] of POS.
    fn inverted_repeat_contig() -> Vec<u8> {
        let mut seq = cycling_contig();
        let inverted = revcomp_of(&carried_insertion()[..INS_KMER_LEN]);
        let at = INS_CARRIED_POS as usize - INVERTED_REPEAT_OFFSET;
        seq[at..at + inverted.len()].copy_from_slice(&inverted);
        seq
    }

    /// One indexed CRAM on [`inverted_repeat_contig`] holding, at
    /// [`INS_CARRIED_POS`], three reads that carry `carried_insertion()` whole
    /// as an `I` operation and two **unedited** reads whose plain
    /// [`TEST_READ_LEN`] match spans the inverted repeat. The second pair is
    /// the hazard in the flesh: reads spike never touched whose own bases hold
    /// a probe k-mer, because the reference under them does.
    fn inverted_repeat_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let seq = inverted_repeat_contig();
        let dir = fixture_dir(tag);
        let fasta_path = write_one_contig_fasta(&dir, "ins_inverted", &seq);
        let header = one_contig_header(seq.len());

        let carried = carried_insertion();
        let flank = (TEST_READ_LEN - carried.len()) / 2;
        let mut records: Vec<noodles::cram::Record> = Vec::new();
        // The unedited reads first: they start further left and the CRAM has to
        // be coordinate-sorted.
        for i in 0..2usize {
            records.push(one_contig_record(
                &seq,
                &format!("unedited{}", i),
                INS_CARRIED_POS as usize + i - TEST_READ_LEN,
                60,
                noodles::sam::alignment::record_buf::Data::default(),
            ));
        }
        for i in 0..3usize {
            let start0 = INS_CARRIED_POS as usize + i - flank;
            let mut bases = seq[start0..start0 + flank].to_vec();
            bases.extend_from_slice(&carried);
            bases.extend_from_slice(&seq[start0 + flank..start0 + 2 * flank]);
            records.push(one_contig_spliced_record(
                &format!("carried{}", i),
                start0,
                &[
                    (Kind::Match, flank),
                    (Kind::Insertion, carried.len()),
                    (Kind::Match, flank),
                ],
                &bases,
            ));
        }

        let cram_path = write_indexed_cram(&dir, "ins_inverted", &fasta_path, &header, &records);
        (
            dir,
            fasta_path.to_str().unwrap().to_string(),
            cram_path.to_str().unwrap().to_string(),
        )
    }

    /// One indexed CRAM on the cycling chrA holding two reads at
    /// [`INS_CARRIED_POS`] that carry `carried_insertion()` as an `I` operation
    /// spelled **lowercase** -- soft-masked bases, which a read's stored
    /// sequence may hold and which must still count against an uppercase ALT.
    fn lowercase_read_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let seq = cycling_contig();
        let dir = fixture_dir(tag);
        let fasta_path = write_one_contig_fasta(&dir, "ins_lower", &seq);
        let header = one_contig_header(seq.len());

        let carried = carried_insertion().to_ascii_lowercase();
        let flank = (TEST_READ_LEN - carried.len()) / 2;
        let mut records: Vec<noodles::cram::Record> = Vec::new();
        for i in 0..2usize {
            let start0 = INS_CARRIED_POS as usize + i - flank;
            let mut bases = seq[start0..start0 + flank].to_vec();
            bases.extend_from_slice(&carried);
            bases.extend_from_slice(&seq[start0 + flank..start0 + 2 * flank]);
            records.push(one_contig_spliced_record(
                &format!("lower{}", i),
                start0,
                &[
                    (Kind::Match, flank),
                    (Kind::Insertion, carried.len()),
                    (Kind::Match, flank),
                ],
                &bases,
            ));
        }

        let cram_path = write_indexed_cram(&dir, "ins_lower", &fasta_path, &header, &records);
        (
            dir,
            fasta_path.to_str().unwrap().to_string(),
            cram_path.to_str().unwrap().to_string(),
        )
    }

    /// An INS truth event at 0-based `pos` whose ALT is the reference base
    /// there followed by `inserted` -- the shape CR7's fix writes and
    /// `load_truth_events` reads back.
    fn ins_event_with_alt(pos: u64, inserted: &[u8]) -> TruthEvent {
        let mut alt = vec![cycling_contig()[pos as usize]];
        alt.extend_from_slice(inserted);
        TruthEvent {
            ins_len: Some(inserted.len() as u64),
            ins_alt: Some(alt),
            ..ins_event("chrA", pos)
        }
    }

    /// The rows one INS event leaves behind on a fixture.
    fn ins_rows(cram: &str, fasta: &str, event: &TruthEvent) -> Vec<CheckResult> {
        let args = args_for(cram, fasta);
        let mut results: Vec<CheckResult> = Vec::new();
        check_event(&args, event, &NearbyRecords::default(), &mut results);
        results
    }

    /// `expected` as the row builds it from the floor and `k`.
    fn ins_sequence_expected(k: usize) -> String {
        format!(">={} with a {}bp alt kmer", MIN_INS_READS, k)
    }

    #[test]
    fn test_the_ins_sequence_row_passes_reads_carrying_the_truth_records_own_bases() {
        let (dir, fasta, cram) = inserted_sequence_cram("ins_seq_pass");
        let rows = ins_rows(
            &cram,
            &fasta,
            &ins_event_with_alt(INS_CARRIED_POS, &carried_insertion()),
        );

        assert_eq!(row_names(&rows), vec![INS_PLANTED, INS_READS, INS_SEQUENCE]);
        let r = named(&rows, INS_SEQUENCE);
        assert_eq!(r.expected, ins_sequence_expected(INS_KMER_LEN));
        assert_eq!(r.observed, "3", "all three reads carry the inserted bases");
        assert!(r.pass, "three reads is above the floor of {}", MIN_INS_READS);
        assert!(r.advisory, "the row is always advisory");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_different_inserted_bases_fail_the_sequence_row_that_ins_reads_passes() {
        // The point of the row, pinned on one input: the same reads, the same
        // CIGARs, the same locus, two truth ALTs of the same length. This is
        // C1's separation -- run A's BAM against run B's truth VCF -- in a
        // unit test.
        let (dir, fasta, cram) = inserted_sequence_cram("ins_seq_wrong");
        let right = carried_insertion();
        let wrong = pseudo_random_bases(29, right.len());
        assert_eq!(wrong.len(), right.len(), "the two ALTs are the same length");
        assert_ne!(wrong, right, "and they are different sequences");

        let ok = ins_rows(&cram, &fasta, &ins_event_with_alt(INS_CARRIED_POS, &right));
        let bad = ins_rows(&cram, &fasta, &ins_event_with_alt(INS_CARRIED_POS, &wrong));

        // `ins_reads` cannot tell the two apart: it reads the CIGAR alone
        // (CR9), and the CIGAR is the same.
        let (a, b) = (named(&ok, INS_READS), named(&bad, INS_READS));
        assert_eq!(a.expected, b.expected);
        assert_eq!(a.observed, "3");
        assert_eq!(b.observed, "3");
        assert!(a.pass && b.pass, "ins_reads PASSes on both");

        // `ins_sequence` does.
        let a = named(&ok, INS_SEQUENCE);
        let b = named(&bad, INS_SEQUENCE);
        assert_eq!(a.observed, "3");
        assert!(a.pass, "the truth record's own bases are in the reads");
        assert_eq!(b.observed, "0", "no read carries the other sequence");
        assert!(!b.pass, "so the row FAILs where ins_reads passed");
        assert!(b.advisory, "and it is still advisory");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_a_reverse_complemented_alt_kmer_counts() {
        // A read's stored sequence is in reference orientation, but the
        // insertion may be read from either side, so both orientations count.
        let inserted = carried_insertion();
        let k = INS_KMER_LEN;
        let carried_by_the_reads = revcomp_of(&inserted);
        for probe in [&inserted[..k], &inserted[inserted.len() - k..]] {
            assert!(
                !carried_by_the_reads.windows(k).any(|w| w == probe),
                "the reads at this locus must carry the probe in no other way \
                 than reverse-complemented"
            );
        }

        let (dir, fasta, cram) = inserted_sequence_cram("ins_seq_revcomp");
        let rows = ins_rows(&cram, &fasta, &ins_event_with_alt(INS_REVCOMP_POS, &inserted));

        let r = named(&rows, INS_SEQUENCE);
        assert_eq!(r.observed, "2", "both reverse-complemented reads count");
        assert!(r.pass);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_only_the_first_and_last_bases_of_a_long_insertion_are_reachable() {
        // Why the two probes are the first and last k bases and not the
        // middle: for an insertion longer than a read, no read contains the
        // middle at all. A 2000 bp insertion's middle k-mer sits 1000 bases
        // in, ten times past the reach of a 100 bp fixture read (or a real
        // 151 bp one).
        let long = long_insertion();
        let k = INS_KMER_LEN;
        let middle = long_middle();
        assert!(
            middle.windows(k).any(|w| w == &long[985..985 + k]),
            "these reads do carry a middle k-mer -- the row must still not \
             count them"
        );
        for probe in [&long[..k], &long[long.len() - k..]] {
            for form in [probe.to_vec(), revcomp_of(probe)] {
                assert!(
                    !middle.windows(k).any(|w| w == form.as_slice()),
                    "and neither end's k-mer, in either orientation, is in the \
                     middle 60 bases"
                );
            }
        }

        let (dir, fasta, cram) = inserted_sequence_cram("ins_seq_long");
        let ends = ins_rows(&cram, &fasta, &ins_event_with_alt(INS_LONG_ENDS_POS, &long));
        let mid = ins_rows(&cram, &fasta, &ins_event_with_alt(INS_LONG_MIDDLE_POS, &long));

        let r = named(&ends, INS_SEQUENCE);
        assert_eq!(r.expected, ins_sequence_expected(k), "k is capped at {}", k);
        assert_eq!(
            r.observed, "2",
            "one read carries the first {} bases, one the last", k
        );
        assert!(r.pass);

        let r = named(&mid, INS_SEQUENCE);
        assert_eq!(r.observed, "0", "no read carries either end");
        assert!(!r.pass, "so the row FAILs on the insertion's middle");
        // And `ins_reads` passes at that same locus, so the FAIL is about the
        // bases and not about there being no evidence of an insertion.
        let r = named(&mid, INS_READS);
        assert_eq!(r.observed, "2");
        assert!(r.pass, "ins_reads sees two clipped reads either way");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_an_alt_that_records_no_sequence_gets_no_ins_sequence_row() {
        // A symbolic `<INS>` did not record the bases, and an older spike's
        // truth VCF must not become a FAIL *on this advisory row* for that --
        // the same rule T1 applied to a missing census field. So: no
        // `ins_sequence` row at all, not a failing one. (`ins_planted`, the
        // counted row since RF13, is still there: with no bases it has nothing
        // to look for and is a failed not-evaluable row, as its plan locked.)
        let args = args_for("/nonexistent/no.bam", "/nonexistent/no.fa");
        let carried = carried_insertion();
        for (what, alt) in [
            ("a symbolic <INS>", Some(b"<INS>".to_vec())),
            ("a symbolic <DUP:TANDEM>", Some(b"<DUP:TANDEM>".to_vec())),
            ("an anchor base alone", Some(b"A".to_vec())),
            ("an empty ALT column", Some(Vec::new())),
            ("no ALT carried at all", None),
        ] {
            let event = TruthEvent {
                ins_alt: alt,
                ..ins_event_with_alt(INS_CARRIED_POS, &carried)
            };
            let mut results: Vec<CheckResult> = Vec::new();
            check_event(&args, &event, &NearbyRecords::default(), &mut results);
            assert_eq!(
                row_names(&results),
                vec![INS_PLANTED, INS_READS],
                "{} records no sequence, so it gets no {} row",
                what,
                INS_SEQUENCE
            );
        }
    }

    #[test]
    fn test_an_insertion_below_the_kmer_bound_is_not_evaluable_and_advisory() {
        let (dir, fasta, cram) = inserted_sequence_cram("ins_seq_short");
        let short = short_insertion();
        assert!(short.len() < MIN_INS_KMER_LEN, "8 bases, below the bound");
        let rows = ins_rows(&cram, &fasta, &ins_event_with_alt(INS_SHORT_POS, &short));

        let r = named(&rows, INS_SEQUENCE);
        assert_eq!(
            r.expected,
            format!(">={} with a >={}bp kmer", MIN_INS_READS, MIN_INS_KMER_LEN)
        );
        assert_eq!(r.observed, format!("alt is {}bp", short.len()));
        assert!(!r.pass, "a check that cannot run is a FAILED check (M10)");
        assert!(r.advisory, "and being advisory it costs nothing in the exit status");
        // The reads are there and `ins_reads` counts them: what the row cannot
        // do is grade an 8-mer, not find evidence.
        let r = named(&rows, INS_READS);
        assert_eq!(r.observed, "2");
        assert!(r.pass);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_an_alt_kmer_already_in_the_reference_makes_the_row_not_evaluable() {
        // The guard: an unedited read would match a probe the reference
        // already holds, so the count would mean nothing. chrA cycles `ACGT`,
        // so an insertion of `ACGT` repeated is exactly that case -- and an
        // explicit insertion copied from nearby sequence is the real one.
        let repeated: Vec<u8> = b"ACGT".repeat(15);
        let k = INS_KMER_LEN;
        assert!(
            cycling_contig().windows(k).any(|w| w == &repeated[..k]),
            "the reference really does hold this probe"
        );

        let (dir, fasta, cram) = inserted_sequence_cram("ins_seq_in_ref");
        let guarded = ins_rows(&cram, &fasta, &ins_event_with_alt(INS_CARRIED_POS, &repeated));
        let r = named(&guarded, INS_SEQUENCE);
        assert_eq!(r.expected, ins_sequence_expected(k));
        assert_eq!(r.observed, "kmer in ref");
        assert!(!r.pass, "not evaluable is a FAILED row (M10)");
        assert!(r.advisory);

        // The same locus with a sequence the reference has not got is graded,
        // so the guard is about this ALT and not about the locus.
        let graded = ins_rows(
            &cram,
            &fasta,
            &ins_event_with_alt(INS_CARRIED_POS, &carried_insertion()),
        );
        assert!(named(&graded, INS_SEQUENCE).pass);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_an_inverted_alt_kmer_in_the_reference_makes_the_row_not_evaluable() {
        // The guard's reverse-complement half, on its own. A read's stored
        // sequence is in reference-forward orientation and `holds_a_probe`
        // matches `probe` **or** `revcomp(probe)`, so a reference holding
        // `revcomp(probe)` near POS -- an inverted repeat of the insertion,
        // within 1 kb -- is matched by reads spike never touched exactly as a
        // forward copy would be. The row must decline it.
        //
        // `test_an_alt_kmer_already_in_the_reference_makes_the_row_not_evaluable`
        // cannot see this: its `ACGT`-repeat ALT is in the `ACGT`-cycling contig
        // forwards as well, so a forward-only guard would fire there too. Here
        // the forward probe is absent from the whole contig, so only a guard
        // that searches both orientations declines the row.
        let carried = carried_insertion();
        let k = INS_KMER_LEN;
        let seq = inverted_repeat_contig();
        for probe in [&carried[..k], &carried[carried.len() - k..]] {
            assert!(
                !seq.windows(k).any(|w| w == probe),
                "the contig must hold neither probe k-mer in its forward \
                 orientation, or a forward-only guard would fire too"
            );
        }
        assert!(
            seq.windows(k)
                .any(|w| w == revcomp_of(&carried[..k]).as_slice()),
            "and it must hold the first probe's reverse complement"
        );
        // The hazard spelled out on the very bases an unedited read over POS
        // reads: they match a probe, so counting them would mean nothing.
        let under_an_unedited_read =
            &seq[INS_CARRIED_POS as usize - TEST_READ_LEN..INS_CARRIED_POS as usize];
        assert!(
            holds_a_probe(under_an_unedited_read, &alt_probe_kmers(&carried, k)),
            "an unedited read over the inverted repeat would be counted as \
             support, which is exactly why the row cannot be graded here"
        );

        let (dir, fasta, cram) = inverted_repeat_cram("ins_seq_inverted");
        let rows = ins_rows(&cram, &fasta, &ins_event_with_alt(INS_CARRIED_POS, &carried));

        let r = named(&rows, INS_SEQUENCE);
        assert_eq!(r.expected, ins_sequence_expected(k));
        assert_eq!(
            r.observed, "kmer in ref",
            "a reverse-complemented copy in the reference is as disqualifying \
             as a forward one"
        );
        assert!(!r.pass, "not evaluable is a FAILED row (M10)");
        assert!(r.advisory);
        // And the evidence really is there: `ins_reads` counts the three reads
        // that carry the insertion, so the declined row is about the reference
        // and not about there being nothing to see.
        let r = named(&rows, INS_READS);
        assert_eq!(r.observed, "3");
        assert!(r.pass);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_a_lowercase_alt_grades_as_its_uppercase_spelling_does() {
        // A soft-masked ALT is legal VCF -- `bcftools` and some callers emit
        // one -- and `alt_probe_kmers` uppercases before it cuts its probes, so
        // the reads that pass on an uppercase ALT pass on the lowercase
        // spelling of the same bases. Without that uppercasing this is an
        // advisory FAIL reading `observed 0` on a correct run, with nothing in
        // the output to say why.
        let (dir, fasta, cram) = inserted_sequence_cram("ins_seq_lower_alt");
        let upper = carried_insertion();
        let lower = upper.to_ascii_lowercase();
        assert_ne!(lower, upper, "the two spellings differ as bytes");

        let rows = ins_rows(&cram, &fasta, &ins_event_with_alt(INS_CARRIED_POS, &lower));
        let r = named(&rows, INS_SEQUENCE);
        assert_eq!(r.expected, ins_sequence_expected(INS_KMER_LEN));
        assert_eq!(
            r.observed, "3",
            "all three reads carry the bases the lowercase ALT names"
        );
        assert!(r.pass, "so a soft-masked ALT grades, and PASSes");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_lowercase_read_bases_hold_an_uppercase_probe() {
        // The other half of the same rule and the other `to_ascii_uppercase`:
        // a read whose stored sequence is soft-masked counts against an
        // uppercase ALT. Pinned end to end on a CRAM whose records really do
        // carry lowercase bases, and then on `holds_a_probe` itself, which is
        // where the uppercasing that makes it work lives.
        let inserted = carried_insertion();
        let (dir, fasta, cram) = lowercase_read_cram("ins_seq_lower_read");
        // The fixture is only worth anything if the lowercase bases survive the
        // CRAM round trip, so read them back and check before grading on them.
        let mut soft_masked = 0usize;
        for_each_alignment(
            &cram,
            &fasta,
            "chrA",
            INS_CARRIED_POS - INS_SEQUENCE_PAD,
            INS_CARRIED_POS + INS_SEQUENCE_PAD,
            0,
            &mut |_name, _start, _ops, seq| {
                if seq.iter().any(|b| b.is_ascii_lowercase()) {
                    soft_masked += 1;
                }
            },
        )
        .unwrap();
        assert_eq!(
            soft_masked, 2,
            "both records really do hand back lowercase bases"
        );

        let rows = ins_rows(&cram, &fasta, &ins_event_with_alt(INS_CARRIED_POS, &inserted));
        let r = named(&rows, INS_SEQUENCE);
        assert_eq!(r.observed, "2", "both soft-masked reads count");
        assert!(r.pass, "so soft-masked reads grade, and PASS");
        let _ = std::fs::remove_dir_all(&dir);

        // And the same rule at the one line that implements it, so a reader of
        // `holds_a_probe` sees why its uppercasing is not redundant with the
        // one `fetch_window` already did for the reference.
        let probes = alt_probe_kmers(&inserted, INS_KMER_LEN);
        assert!(
            holds_a_probe(&inserted.to_ascii_lowercase(), &probes),
            "lowercase bases hold an uppercase probe"
        );
    }

    #[test]
    fn test_the_ins_sequence_row_is_advisory_in_the_table_and_in_the_json() {
        // An errored row is advisory too: `check_outcome` stamps both arms, so
        // a BAM that cannot be read cannot leak a non-advisory FAIL into the
        // exit status.
        let rows = ins_rows(
            "/nonexistent/no.bam",
            "/nonexistent/no.fa",
            &ins_event_with_alt(INS_CARRIED_POS, &carried_insertion()),
        );

        let r = named(&rows, INS_SEQUENCE);
        assert_eq!(r.expected, "check runs");
        assert!(!r.pass);
        assert!(r.advisory, "an errored ins_sequence row is still advisory");
        assert!(
            !named(&rows, INS_PLANTED).advisory,
            "while the ins_planted row beside it is not"
        );
        assert!(named(&rows, INS_READS).advisory, "and since RF13 ins_reads is advisory too");

        assert_eq!(
            failure_message(&rows, false).as_deref(),
            Some("1/1 validation checks failed"),
            "the advisory rows are in neither the count nor the total by default"
        );
        assert_eq!(
            failure_message(&rows, true).as_deref(),
            Some("3/3 validation checks failed"),
            "--strict counts them"
        );

        let report = text_report(&rows, false);
        assert!(
            line_for(&report, INS_SEQUENCE).ends_with("FAIL (advisory)"),
            "in:\n{}",
            report
        );
        assert!(
            line_for(&report, INS_PLANTED).ends_with("FAIL")
                && !line_for(&report, INS_PLANTED).ends_with("(advisory)"),
            "in:\n{}",
            report
        );
        assert!(
            line_for(&report, INS_READS).ends_with("FAIL (advisory)"),
            "in:\n{}",
            report
        );

        let mut json: Vec<u8> = Vec::new();
        print_results_json(&mut json, &rows, rows.len(), 0, false).unwrap();
        let json = String::from_utf8(json).unwrap();
        assert!(
            line_for(&json, "\"ins_sequence\"").contains("\"advisory\": true"),
            "in:\n{}",
            json
        );
        assert!(
            line_for(&json, "\"ins_planted\"").contains("\"advisory\": false"),
            "in:\n{}",
            json
        );
        assert!(
            line_for(&json, "\"ins_reads\"").contains("\"advisory\": true"),
            "in:\n{}",
            json
        );
    }

    #[test]
    fn test_the_ins_sequence_row_is_produced_for_ins_and_for_no_other_type() {
        // The gate is the event's type, not whether an ALT happens to be
        // carried: every event below carries one.
        let args = args_for("/nonexistent/no.bam", "/nonexistent/no.fa");
        let alt = Some({
            let mut a = vec![b'A'];
            a.extend_from_slice(&carried_insertion());
            a
        });
        for sv_type in ["DEL", "DUP", "INV", "BND", "SNP", "CNV"] {
            let event = TruthEvent {
                sv_type: sv_type.to_string(),
                ins_alt: alt.clone(),
                ..del_event("chrA", 5_000, 7_000)
            };
            let mut results: Vec<CheckResult> = Vec::new();
            check_event(&args, &event, &NearbyRecords::default(), &mut results);
            assert!(
                !row_names(&results).contains(&INS_SEQUENCE),
                "{} has no insertion, so it gets no {} row; got {:?}",
                sv_type,
                INS_SEQUENCE,
                row_names(&results)
            );
        }

        let event = TruthEvent {
            ins_alt: alt,
            ..ins_event_with_alt(INS_CARRIED_POS, &carried_insertion())
        };
        let mut results: Vec<CheckResult> = Vec::new();
        check_event(&args, &event, &NearbyRecords::default(), &mut results);
        assert_eq!(row_names(&results), vec![INS_PLANTED, INS_READS, INS_SEQUENCE]);
    }

    #[test]
    fn test_the_ins_sequence_expectation_is_built_from_the_floor_and_the_kmer() {
        // `expected` is derived, so raising the floor or the k-mer length
        // cannot leave the string behind claiming the old one -- and it has to
        // survive the 24-wide Expected column at both ends of k's range.
        let (dir, fasta, cram) = inserted_sequence_cram("ins_seq_expected");

        // An insertion longer than the cap: k is the cap.
        let rows = ins_rows(&cram, &fasta, &ins_event_with_alt(INS_CARRIED_POS, &long_insertion()));
        assert_eq!(
            named(&rows, INS_SEQUENCE).expected,
            ins_sequence_expected(INS_KMER_LEN)
        );

        // A shorter one: k is the insertion's own length.
        for len in [MIN_INS_KMER_LEN, 20, INS_KMER_LEN] {
            let rows = ins_rows(
                &cram,
                &fasta,
                &ins_event_with_alt(INS_CARRIED_POS, &pseudo_random_bases(31, len)),
            );
            assert_eq!(
                named(&rows, INS_SEQUENCE).expected,
                ins_sequence_expected(len),
                "k follows the inserted length below the cap"
            );
        }

        for k in [MIN_INS_KMER_LEN, INS_KMER_LEN] {
            assert!(
                ins_sequence_expected(k).len() <= 24,
                "{:?} must survive the 24-wide Expected column",
                ins_sequence_expected(k)
            );
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_the_ins_reads_row_is_unchanged_by_sharing_its_minimum() {
        // Constraint 4: lifting `MIN_INS_READS` out of `check_ins_reads` to
        // module scope must not change one character the non-advisory row
        // prints, so the strings are pinned as literals here rather than
        // against the constant they are built from.
        let (dir, fasta, cram) = inserted_sequence_cram("ins_seq_unchanged");
        let rows = ins_rows(
            &cram,
            &fasta,
            &ins_event_with_alt(INS_CARRIED_POS, &carried_insertion()),
        );

        let r = named(&rows, INS_READS);
        assert_eq!(r.check_name, "ins_reads");
        assert_eq!(r.expected, ">=2 reads with >=50bp inserted at chrA:5000");
        assert_eq!(r.observed, "3");
        assert!(r.pass);
        // The one thing that moved, on purpose: since RF13 the row says how
        // the aligner wrote the insertion and no longer decides.
        assert!(r.advisory);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_load_truth_events_carries_an_ins_records_alt_and_nothing_elses() {
        let vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chrA\t5001\tsim_ins_1\tA\tAGGTT\t999\tPASS\tSVTYPE=INS;SVLEN=4;SIM_VAF=0.500\tGT\t0/1
chrA\t6001\tsim_ins_2\tA\t<INS>\t999\tPASS\tSVTYPE=INS;SVLEN=200;SIM_VAF=0.500\tGT\t0/1
chrA\t7001\tsim_del_1\tN\t<DEL>\t999\tPASS\tSVTYPE=DEL;END=8000;SIM_VAF=0.500\tGT\t0/1
chrA\t9001\tsim_var_1\tA\tT\t999\tPASS\tSIM_VAF=0.500\tGT\t0/1
";
        let dir = fixture_dir("ins_seq_truth");
        let path = dir.join("truth.vcf");
        std::fs::write(&path, vcf).unwrap();
        let events = load_truth_events(path.to_str().unwrap()).unwrap();

        assert_eq!(
            events.iter().map(|e| e.sv_type.as_str()).collect::<Vec<_>>(),
            vec!["INS", "INS", "DEL", "SNP"]
        );
        // The INS records carry their ALT column verbatim, symbolic or not.
        assert_eq!(events[0].ins_alt.as_deref(), Some(b"AGGTT".as_slice()));
        assert_eq!(events[1].ins_alt.as_deref(), Some(b"<INS>".as_slice()));
        // And nothing else does.
        assert_eq!(events[2].ins_alt, None);
        assert_eq!(events[3].ins_alt, None);

        // No other event type's parse moved: the small variant still carries
        // its own REF and ALT, and the INS records still carry neither -- which
        // is what keeps an insertion out of `NearbyRecords`, where it would be
        // applied to the reference at its own `start` and move the
        // non-advisory `allele_freq` row.
        assert_eq!(events[3].ref_allele.as_deref(), Some(b"A".as_slice()));
        assert_eq!(events[3].alt_allele.as_deref(), Some(b"T".as_slice()));
        for e in &events[..3] {
            assert_eq!(e.ref_allele, None, "{} carries no REF allele", e.sv_type);
            assert_eq!(e.alt_allele, None, "{} carries no ALT allele", e.sv_type);
        }
        let nearby = NearbyRecords::new(&events);
        assert_eq!(
            nearby.inside(&events[3], 0, 20_000).len(),
            0,
            "only the small variant is a nearby record, and it is the event itself"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_the_ins_sequence_row_does_not_move_the_status_column() {
        // 12 characters in an 18-wide Check column, so unlike
        // `split_reads_each_end` this row costs the table nothing.
        assert_eq!(INS_SEQUENCE, "ins_sequence");
        assert_eq!(INS_SEQUENCE.len(), 12);
        let report = text_report(
            &[
                row(COVERAGE_RATIO, true, false),
                row(INS_READS, true, false),
                row(INS_SEQUENCE, true, true),
            ],
            false,
        );
        let status_at = |name: &str| {
            let line = line_for(&report, name);
            line.find("PASS")
                .unwrap_or_else(|| panic!("no status on {:?}", line))
        };
        assert_eq!(
            status_at(INS_SEQUENCE),
            status_at(COVERAGE_RATIO),
            "the row must leave the Status column where every other row has it:\n{}",
            report
        );
        assert_eq!(status_at(INS_READS), status_at(COVERAGE_RATIO));
    }

    // ── ins_planted (RF13) ─────────────────────────────────────────────────

    /// `one_contig_spliced_record` with the BAM flags and MAPQ given: the
    /// `ins_planted` fixture needs reads a usable-alignment filter would drop.
    fn spliced_record_with(
        name: &str,
        start0: usize,
        ops: &[(Kind, usize)],
        bases: &[u8],
        bam_flags: u16,
        mapq: u8,
    ) -> noodles::cram::Record {
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};

        let flags = noodles::cram::record::Flags::QUALITY_SCORES_STORED_AS_ARRAY;
        let cigar: Cigar = ops.iter().map(|&(k, n)| Op::new(k, n)).collect();
        let sequence = Sequence::from(bases.to_vec());
        let quality_scores = QualityScores::from(vec![40u8; bases.len()]);
        let features = noodles::cram::record::Features::from_cigar(
            flags,
            &cigar,
            &sequence,
            &quality_scores,
        );
        noodles::cram::Record::builder()
            .set_bam_flags(noodles::sam::alignment::record::Flags::from(bam_flags))
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_read_length(bases.len())
            .set_alignment_start(noodles::core::Position::new(start0 + 1).unwrap())
            .set_name(name)
            .set_mapping_quality(
                noodles::sam::alignment::record::MappingQuality::new(mapq).unwrap(),
            )
            .set_bases(sequence)
            .set_quality_scores(quality_scores)
            .set_features(features)
            .build()
    }

    /// One indexed CRAM with six reads at [`INS_CARRIED_POS`], every one
    /// carrying `carried_insertion()` whole as an `I` operation, so only their
    /// names and flags tell them apart:
    ///
    /// - `ev0001_hap_000001` and `_000002`: event 1's own, MAPQ 60, unflagged;
    /// - `ev0001_hap_000003`: event 1's own, **MAPQ 0 and duplicate-flagged**;
    /// - `ev0001_hap_000004`: event 1's own, but **secondary**;
    /// - `A00744:46:HV3C3DSXX:2:1104:13160:1000`: a read of the sample's own;
    /// - `ev0002_hap_000001`: **another event's** read.
    ///
    /// For truth record `sim_ins_1`, `ins_planted` must count exactly the first
    /// three.
    fn planted_reads_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let seq = cycling_contig();
        let dir = fixture_dir(tag);
        let fasta_path = write_one_contig_fasta(&dir, "ins_planted", &seq);
        let header = one_contig_header(seq.len());

        let carried = carried_insertion();
        let flank = (TEST_READ_LEN - carried.len()) / 2;
        let reads: [(&str, u16, u8); 6] = [
            ("ev0001_hap_000001", 0, 60),
            ("ev0001_hap_000002", 0, 60),
            ("ev0001_hap_000003", 0x400, 0),
            ("ev0001_hap_000004", 0x100, 60),
            ("A00744:46:HV3C3DSXX:2:1104:13160:1000", 0, 60),
            ("ev0002_hap_000001", 0, 60),
        ];
        // All six at one start, so every read carries the insertion at POS
        // itself: `inserted_sequence_cram`'s `+ i` staggering would put read i's
        // insertion i bases away, which on the cycling contig is a different
        // junction and rightly carries nothing.
        let start0 = INS_CARRIED_POS as usize - flank;
        let mut records: Vec<noodles::cram::Record> = Vec::new();
        for (name, flags, mapq) in reads.iter() {
            let mut bases = seq[start0..start0 + flank].to_vec();
            bases.extend_from_slice(&carried);
            bases.extend_from_slice(&seq[start0 + flank..start0 + 2 * flank]);
            records.push(spliced_record_with(
                name,
                start0,
                &[
                    (Kind::Match, flank),
                    (Kind::Insertion, carried.len()),
                    (Kind::Match, flank),
                ],
                &bases,
                *flags,
                *mapq,
            ));
        }
        let cram_path = write_indexed_cram(&dir, "ins_planted", &fasta_path, &header, &records);
        (
            dir,
            fasta_path.to_str().unwrap().to_string(),
            cram_path.to_str().unwrap().to_string(),
        )
    }

    /// Truth record `sim_ins_1` at [`INS_CARRIED_POS`], carrying `inserted`.
    fn planted_event(inserted: &[u8]) -> TruthEvent {
        TruthEvent {
            sim_number: Some(1),
            ..ins_event_with_alt(INS_CARRIED_POS, inserted)
        }
    }

    #[test]
    fn test_ins_planted_counts_only_the_events_own_reads() {
        let (dir, fasta, cram) = planted_reads_cram("planted_own");
        let rows = ins_rows(&cram, &fasta, &planted_event(&carried_insertion()));
        let r = named(&rows, INS_PLANTED);
        // Event 1's three primary reads, the MAPQ 0 duplicate among them: how
        // the aligner scored or flagged a read spike made is not the question.
        // Not the secondary record, not the sample's read carrying the same
        // bases, and not event 2's read.
        assert_eq!(r.observed, "3", "{} / {}", r.expected, r.observed);
        assert!(r.pass);
        assert!(!r.advisory, "ins_planted is the INS event's counted row");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_ins_planted_passes_on_one_carrying_read() {
        // Event 2 has one read in the fixture. With the sample's reads left
        // out there is no background, so one of spike's own reads carrying the
        // insertion is the planted event (MIN_PLANTED_READS).
        let (dir, fasta, cram) = planted_reads_cram("planted_one");
        let event = TruthEvent {
            sim_number: Some(2),
            ..planted_event(&carried_insertion())
        };
        let rows = ins_rows(&cram, &fasta, &event);
        let r = named(&rows, INS_PLANTED);
        assert_eq!(r.observed, "1", "{} / {}", r.expected, r.observed);
        assert!(r.pass, "one carrying read of spike's own is enough");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_ins_planted_fails_a_truth_whose_inserted_letters_are_wrong() {
        let (dir, fasta, cram) = planted_reads_cram("planted_wrong");
        // Every inserted base swapped for a different one, as RF13's N2.
        let wrong: Vec<u8> = carried_insertion()
            .iter()
            .map(|b| match b {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            })
            .collect();
        let rows = ins_rows(&cram, &fasta, &planted_event(&wrong));
        let r = named(&rows, INS_PLANTED);
        assert_eq!(r.observed, "0", "{} / {}", r.expected, r.observed);
        assert!(!r.pass);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_ins_planted_is_not_evaluable_without_a_spike_id() {
        let (dir, fasta, cram) = planted_reads_cram("planted_noid");
        let event = TruthEvent {
            sim_number: None,
            ..planted_event(&carried_insertion())
        };
        let rows = ins_rows(&cram, &fasta, &event);
        let r = named(&rows, INS_PLANTED);
        assert!(!r.pass, "no sim_ins_N ID, no way to know spike's reads");
        assert!(r.observed.contains("id"), "{} / {}", r.expected, r.observed);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_ins_planted_is_not_evaluable_without_inserted_bases() {
        let (dir, fasta, cram) = planted_reads_cram("planted_noalt");
        let event = TruthEvent {
            ins_alt: Some(b"<INS>".to_vec()),
            ..planted_event(&carried_insertion())
        };
        let rows = ins_rows(&cram, &fasta, &event);
        let r = named(&rows, INS_PLANTED);
        assert!(!r.pass, "no bases, nothing to look for");
        assert!(r.observed.contains("alt"), "{} / {}", r.expected, r.observed);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_an_ins_probe_allows_two_flank_substitutions_and_no_inserted_one() {
        let left = pseudo_random_bases(41, 31);
        let right = pseudo_random_bases(43, 31);
        let inserted = pseudo_random_bases(47, 60);
        let probes = ins_junction_probes(&left, &inserted, &right);
        // A read over the left junction only: the last 31 flank bases, then
        // 40 inserted ones (not reaching the right junction).
        let read = |edit: &dyn Fn(&mut Vec<u8>)| {
            let mut r = left.clone();
            r.extend_from_slice(&inserted[..40]);
            edit(&mut r);
            r
        };
        // Probe positions: flank bases left[16..31], inserted bases at 31..47.
        let swap = |b: u8| if b == b'A' { b'C' } else { b'A' };
        assert!(carries_junction_probe(&read(&|_| {}), &probes), "exact");
        assert!(
            carries_junction_probe(&read(&|r| { r[20] = swap(r[20]); r[25] = swap(r[25]); }), &probes),
            "two flank substitutions: the sample's SNPs and errors"
        );
        assert!(
            !carries_junction_probe(
                &read(&|r| { r[18] = swap(r[18]); r[22] = swap(r[22]); r[27] = swap(r[27]); }),
                &probes
            ),
            "three flank substitutions"
        );
        assert!(
            !carries_junction_probe(&read(&|r| { r[35] = swap(r[35]); }), &probes),
            "one substitution in an inserted base"
        );
        assert!(
            carries_junction_probe(&revcomp_of(&read(&|_| {})), &probes),
            "either orientation"
        );
    }

    #[test]
    fn test_a_short_insertions_wrong_letters_do_not_carry() {
        // RF13's first attempt: a 4 bp insertion with wrong letters matched
        // within a tolerance that spanned the inserted bases.
        let left = pseudo_random_bases(53, 31);
        let right = pseudo_random_bases(59, 31);
        let inserted = b"GATC".to_vec();
        let wrong = b"TCGA".to_vec();
        let mut read = left.clone();
        read.extend_from_slice(&inserted);
        read.extend_from_slice(&right);
        assert!(carries_junction_probe(&read, &ins_junction_probes(&left, &inserted, &right)));
        assert!(!carries_junction_probe(&read, &ins_junction_probes(&left, &wrong, &right)));
    }

    #[test]
    fn test_the_ins_rows_are_ins_planted_then_ins_reads_and_ins_sequence_advisory() {
        let (dir, fasta, cram) = planted_reads_cram("planted_rows");
        let rows = ins_rows(&cram, &fasta, &planted_event(&carried_insertion()));
        let ins: Vec<(&str, bool)> = rows
            .iter()
            .filter(|r| r.check_name.starts_with("ins_"))
            .map(|r| (r.check_name.as_str(), r.advisory))
            .collect();
        assert_eq!(
            ins,
            vec![(INS_PLANTED, false), (INS_READS, true), (INS_SEQUENCE, true)],
            "ins_reads says how the aligner wrote it; it no longer decides (RF13)"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_load_truth_events_reads_the_spike_event_number_from_the_id() {
        let dir = fixture_dir("planted_ids");
        let path = dir.join("truth_ins_ids.vcf");
        std::fs::write(
            &path,
            "##fileformat=VCFv4.3\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
             chr1\t100\tsim_ins_7\tA\tAGGGG\t999\tPASS\tSVTYPE=INS;SVLEN=4\tGT\t0/1\n\
             chr1\t500\treal_ins_5\tA\tAGGGG\t999\tPASS\tSVTYPE=INS;SVLEN=4\tGT\t0/1\n\
             chr1\t900\tsim_ins_7x\tA\tAGGGG\t999\tPASS\tSVTYPE=INS;SVLEN=4\tGT\t0/1\n",
        )
        .unwrap();
        let events = load_truth_events(path.to_str().unwrap()).unwrap();
        let numbers: Vec<Option<u32>> = events.iter().map(|e| e.sim_number).collect();
        assert_eq!(numbers, vec![Some(7), None, None]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    // ── del_planted (RF14) ─────────────────────────────────────────────────

    /// chrA for the `del_planted` fixture: 20 kb of seeded random bases, so a
    /// 16-base stretch is unique and a probe shifted by one base is a different
    /// probe (the cycling contig repeats every 4 bases).
    fn del_planted_contig() -> Vec<u8> {
        pseudo_random_bases(61, 20_000)
    }

    /// One read of the `del_planted` fixture: name, 0-based start, CIGAR,
    /// bases, BAM flags and MAPQ.
    type DelPlantedRead = (&'static str, usize, Vec<(Kind, usize)>, Vec<u8>, u16, u8);

    /// The fixture deletion: `[DEL_PLANTED_START, DEL_PLANTED_END)`, 3 kb, so
    /// the two 500 bp breakpoint windows do not overlap.
    const DEL_PLANTED_START: u64 = 5_000;
    const DEL_PLANTED_END: u64 = 8_000;

    /// One indexed CRAM of 100 bp reads spelling the deletion's join
    /// `seq[..5000] + seq[8000..]`, told apart by name, placement, flags and
    /// bases:
    ///
    /// - `ev0001_hap_000001`: event 1's, aligned at the START side (60M40S);
    /// - `ev0001_hap_000002`: event 1's, aligned at the END side only (40S60M),
    ///   inside END's window and outside START's;
    /// - `ev0001_hap_000003`: event 1's, **MAPQ 0 and duplicate-flagged**;
    /// - `ev0001_hap_000004`: event 1's, but **secondary**;
    /// - `ev0001_hap_000005`: event 1's, **one** base of the probe changed;
    /// - `ev0001_hap_000006`: event 1's, **three** bases of the probe changed;
    /// - `A00744:46:HV3C3DSXX:2:1104:13160:1000`: a read of the sample's own;
    /// - `ev0002_hap_000001`: **another event's** read.
    ///
    /// For truth record `sim_del_1`, `del_planted` must count exactly 1, 2, 3
    /// and 5.
    fn del_planted_cram(tag: &str) -> (std::path::PathBuf, String, String) {
        let seq = del_planted_contig();
        let dir = fixture_dir(tag);
        let fasta_path = write_one_contig_fasta(&dir, "del_planted", &seq);
        let header = one_contig_header(seq.len());
        let (s, e) = (DEL_PLANTED_START as usize, DEL_PLANTED_END as usize);

        // Left-placed carrier: 60 bases before START, 40 from END on. Read
        // index 45..60 is the probe's left 15, 60..76 its right 16.
        let left_placed = || {
            let mut b = seq[s - 60..s].to_vec();
            b.extend_from_slice(&seq[e..e + 40]);
            b
        };
        let swap = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let with_changes = |at: &[usize]| {
            let mut b = left_placed();
            for &i in at {
                b[i] = swap(b[i]);
            }
            b
        };
        let left_ops = [(Kind::Match, 60), (Kind::SoftClip, 40)];
        let mut right_placed = seq[s - 40..s].to_vec();
        right_placed.extend_from_slice(&seq[e..e + 60]);

        let mut reads: Vec<DelPlantedRead> = vec![
            ("ev0001_hap_000001", s - 60, left_ops.to_vec(), left_placed(), 0, 60),
            ("ev0001_hap_000003", s - 60, left_ops.to_vec(), left_placed(), 0x400, 0),
            ("ev0001_hap_000004", s - 60, left_ops.to_vec(), left_placed(), 0x100, 60),
            ("ev0001_hap_000005", s - 60, left_ops.to_vec(), with_changes(&[55]), 0, 60),
            ("ev0001_hap_000006", s - 60, left_ops.to_vec(), with_changes(&[50, 55, 65]), 0, 60),
            ("A00744:46:HV3C3DSXX:2:1104:13160:1000", s - 60, left_ops.to_vec(), left_placed(), 0, 60),
            ("ev0002_hap_000001", s - 60, left_ops.to_vec(), left_placed(), 0, 60),
            (
                "ev0001_hap_000002",
                e,
                vec![(Kind::SoftClip, 40), (Kind::Match, 60)],
                right_placed,
                0,
                60,
            ),
        ];
        reads.sort_by_key(|r| r.1);
        let records: Vec<noodles::cram::Record> = reads
            .iter()
            .map(|(name, start0, ops, bases, flags, mapq)| {
                spliced_record_with(name, *start0, ops, bases, *flags, *mapq)
            })
            .collect();
        let cram_path = write_indexed_cram(&dir, "del_planted", &fasta_path, &header, &records);
        (
            dir,
            fasta_path.to_str().unwrap().to_string(),
            cram_path.to_str().unwrap().to_string(),
        )
    }

    /// Truth record `sim_del_1` over the fixture's deletion, ending at `end`.
    fn del_planted_event(end: u64) -> TruthEvent {
        TruthEvent {
            sim_number: Some(1),
            ..del_event("chrA", DEL_PLANTED_START, end)
        }
    }

    /// The rows one event leaves behind on a fixture.
    fn event_rows(cram: &str, fasta: &str, event: &TruthEvent) -> Vec<CheckResult> {
        let args = args_for(cram, fasta);
        let mut results: Vec<CheckResult> = Vec::new();
        check_event(&args, event, &NearbyRecords::default(), &mut results);
        results
    }

    #[test]
    fn test_del_planted_counts_only_the_events_own_reads() {
        let (dir, fasta, cram) = del_planted_cram("del_planted_own");
        let rows = event_rows(&cram, &fasta, &del_planted_event(DEL_PLANTED_END));
        let r = named(&rows, DEL_PLANTED);
        // Event 1's reads 1, 2, 3 and 5: the one placed at END alone, the MAPQ 0
        // duplicate (how the aligner scored or flagged spike's read is not the
        // question) and the one with a single changed base (a sample SNP or an
        // error). Not the secondary record, not the one three bases off, not
        // the sample's read and not event 2's.
        assert_eq!(r.observed, "4", "{} / {}", r.expected, r.observed);
        assert!(r.pass);
        assert!(!r.advisory, "del_planted is the DEL event's counted row");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_del_planted_passes_on_one_carrying_read() {
        // Event 2 has one read in the fixture; with the sample's reads left out
        // one carrying read of spike's own is the planted event.
        let (dir, fasta, cram) = del_planted_cram("del_planted_one");
        let event = TruthEvent {
            sim_number: Some(2),
            ..del_planted_event(DEL_PLANTED_END)
        };
        let rows = event_rows(&cram, &fasta, &event);
        let r = named(&rows, DEL_PLANTED);
        assert_eq!(r.observed, "1", "{} / {}", r.expected, r.observed);
        assert!(r.pass, "one carrying read of spike's own is enough");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_del_planted_fails_a_truth_whose_end_is_50_bp_off() {
        // RF14's N2a: the same reads, a truth record whose deletion runs 50 bp
        // too far. Its join is not what the reads spell.
        let (dir, fasta, cram) = del_planted_cram("del_planted_off");
        let rows = event_rows(&cram, &fasta, &del_planted_event(DEL_PLANTED_END + 50));
        let r = named(&rows, DEL_PLANTED);
        assert_eq!(r.observed, "0", "{} / {}", r.expected, r.observed);
        assert!(!r.pass);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_del_planted_is_not_evaluable_without_a_spike_id() {
        let (dir, fasta, cram) = del_planted_cram("del_planted_noid");
        let event = TruthEvent {
            sim_number: None,
            ..del_planted_event(DEL_PLANTED_END)
        };
        let rows = event_rows(&cram, &fasta, &event);
        let r = named(&rows, DEL_PLANTED);
        assert!(!r.pass, "no sim_del_N ID, no way to know spike's reads");
        assert!(r.observed.contains("id"), "{} / {}", r.expected, r.observed);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_a_del_probe_allows_two_substitutions_and_not_three() {
        let left = pseudo_random_bases(67, 40);
        let right = pseudo_random_bases(71, 40);
        let probes = del_junction_probes(&left[25..], &right[..16]);
        // A read over the join: 40 bases before it, 40 after. The probe is read
        // positions 25..56 (15 before the join, 16 after).
        let read = |at: &[usize]| {
            let mut r = left.clone();
            r.extend_from_slice(&right);
            for &i in at {
                r[i] = if r[i] == b'A' { b'C' } else { b'A' };
            }
            r
        };
        assert!(carries_junction_probe(&read(&[]), &probes), "exact");
        assert!(carries_junction_probe(&read(&[30, 50]), &probes), "two substitutions");
        assert!(
            !carries_junction_probe(&read(&[28, 38, 50]), &probes),
            "three substitutions"
        );
        assert!(carries_junction_probe(&revcomp_of(&read(&[])), &probes), "either orientation");
        // A read of the reference with no deletion does not carry it.
        let unjoined = pseudo_random_bases(67, 80);
        assert!(!carries_junction_probe(&unjoined, &probes), "no join, no carrier");
    }

    #[test]
    fn test_the_del_rows_are_del_planted_counted_and_split_reads_advisory() {
        let (dir, fasta, cram) = del_planted_cram("del_planted_rows");
        let flags = |event: &TruthEvent| -> Vec<(String, bool)> {
            event_rows(&cram, &fasta, event)
                .iter()
                .filter(|r| {
                    [COVERAGE_RATIO, COVERAGE_ANY_MAPQ, DEL_PLANTED, SPLIT_READS, SPLIT_READS_EACH_END]
                        .contains(&r.check_name.as_str())
                })
                .map(|r| (r.check_name.clone(), r.advisory))
                .collect()
        };
        let own = |v: &[(&str, bool)]| -> Vec<(String, bool)> {
            v.iter().map(|&(n, a)| (n.to_string(), a)).collect()
        };
        assert_eq!(
            flags(&del_planted_event(DEL_PLANTED_END)),
            own(&[
                (COVERAGE_RATIO, false),
                (COVERAGE_ANY_MAPQ, true),
                (DEL_PLANTED, false),
                (SPLIT_READS, true),
                (SPLIT_READS_EACH_END, true),
            ]),
            "split_reads says how the aligner split it; it no longer decides a DEL (RF14)"
        );
        let dup = TruthEvent {
            sv_type: "DUP".to_string(),
            ..del_planted_event(DEL_PLANTED_END)
        };
        assert_eq!(
            flags(&dup),
            own(&[
                (COVERAGE_RATIO, false),
                (COVERAGE_ANY_MAPQ, true),
                (SPLIT_READS, false),
                (SPLIT_READS_EACH_END, true),
            ]),
            "a DUP was not measured: its split_reads still decides, and it gets no del_planted"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_load_truth_events_reads_the_del_event_number_from_the_id() {
        let dir = fixture_dir("del_planted_ids");
        let path = dir.join("truth_del_ids.vcf");
        std::fs::write(
            &path,
            "##fileformat=VCFv4.3\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
             chr1\t100\tsim_del_4\tN\t<DEL>\t999\tPASS\tSVTYPE=DEL;END=500;SVLEN=-400\tGT\t0/1\n\
             chr1\t1000\treal_del_1000\tN\t<DEL>\t999\tPASS\tSVTYPE=DEL;END=1400;SVLEN=-400\tGT\t0/1\n\
             chr1\t2000\tsim_ins_4\tN\t<DEL>\t999\tPASS\tSVTYPE=DEL;END=2400;SVLEN=-400\tGT\t0/1\n\
             chr1\t3000\tsim_ins_7\tA\tAGGGG\t999\tPASS\tSVTYPE=INS;SVLEN=4\tGT\t0/1\n\
             chr1\t4000\tsim_del_7\tA\tAGGGG\t999\tPASS\tSVTYPE=INS;SVLEN=4\tGT\t0/1\n",
        )
        .unwrap();
        let events = load_truth_events(path.to_str().unwrap()).unwrap();
        let numbers: Vec<Option<u32>> = events.iter().map(|e| e.sim_number).collect();
        // Each type reads only its own spike ID.
        assert_eq!(numbers, vec![Some(4), None, None, Some(7), None]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_planted_read_prefix_spells_spikes_read_names() {
        // `simulate.rs` names event N's tiled reads `format!("ev{:04}", N)` +
        // `_hap_` + a counter; simulate's own test pins that the two agree.
        assert_eq!(planted_read_prefix(1), "ev0001_hap_");
        assert_eq!(planted_read_prefix(123), "ev0123_hap_");
    }
}
