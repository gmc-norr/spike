import re, sys, collections
txt = open(sys.argv[1]).read()
samples = re.split(r"^=== SAMPLE \d+ t=[0-9.]+\n", txt, flags=re.M)[1:]
idle_pat = re.compile(r"syscall \(\)|futex|rayon_core::sleep|wait_until|__GI___clock_nanosleep|pthread_cond")
def busy_stack(s):
    threads = re.split(r"^Thread \d+ ", s, flags=re.M)[1:]
    stacks = []
    for t in threads:
        frames = re.findall(r"^#\d+\s+(?:0x[0-9a-f]+ in )?(.*)$", t, flags=re.M)
        stacks.append(frames)
    for f in stacks:
        if f and not idle_pat.search(f[0]):
            return f
    return None
phases = [
 ("dupscan", r"removed_duplicates|duplicates_of"),
 ("loh_count", r"count_alleles"),
 ("loh_collect", r"collect_snp_alleles"),
 ("phasing", r"pick_target_alleles|classify_from_collected|copies_from_snps|call_snps"),
 ("class_mix", r"class_mix"),
 ("census", r"census::"),
 ("extract_pass1", r"pass1_bam"),
 ("extract_other", r"extract::"),
 ("fastq_write", r"write_paired_fastq|spike::fastq"),
 ("quality_startup", r"spike::quality|QualityProfile"),
 ("tiling", r"simulate::|synth::|haplotype::"),
 ("post", r"combine_event_outputs|consumed_original|fastq_removed_names"),
 ("gather/origin", r"spike::origin"),
]
leafcats = [
 ("inflate", r"inflate|libdeflate|miniz|flate2|zlib|crc32|adler"),
 ("recordbuf_convert", r"RecordBuf|record_buf|try_from_alignment_record|try_into_alignment_record"),
 ("origin_record_build", r"origin_record|parse_xa|five_prime|placements|reference_length"),
 ("alloc/free/copy", r"malloc|free|realloc|memcpy|memmove|cfree|_int_|alloc::"),
 ("btree", r"BTree|btree"),
 ("hash", r"hashbrown|HashMap|HashSet|sip|Hasher|hash::"),
 ("walk_cigar", r"walk_cigar"),
 ("bam_record_read", r"bam::io::reader|bgzf|read_record|Query|query"),
 ("seq_decode", r"Sequence|sequence"),
]
phase_counts = collections.Counter()
leaf = collections.defaultdict(collections.Counter)
for s in samples:
    f = busy_stack(s)
    if f is None:
        phase_counts["idle"] += 1; continue
    joined = "\n".join(f)
    ph = "other"
    for name, pat in phases:
        if re.search(pat, joined):
            ph = name; break
    phase_counts[ph] += 1
    if ph in ("dupscan", "loh_count", "loh_collect", "extract_pass1", "census", "class_mix", "post"):
        cat = "other:" + f[0][:60]
        for depth in range(min(6, len(f))):
            hit = None
            for name, pat in leafcats:
                if re.search(pat, f[depth]):
                    hit = name; break
            if hit:
                cat = hit; break
        leaf[ph][cat] += 1
tot = sum(phase_counts.values())
print("samples", tot)
for k, v in phase_counts.most_common():
    print("%-16s %3d  %5.1f%%" % (k, v, 100 * v / tot))
for ph in leaf:
    print("\n--", ph)
    for k, v in leaf[ph].most_common(12):
        print("   %-40s %3d" % (k, v))
