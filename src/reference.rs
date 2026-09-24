use anyhow::{Context, Result};
use noodles::fasta;
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::Read;
use std::path::Path;

/// True if `path`'s first two bytes are the gzip magic number, regardless of
/// what its extension says. Any error opening/reading it is treated as "no"
/// so the normal open path below still runs and reports its own failure.
fn looks_gzip_compressed(path: &Path) -> bool {
    let mut magic = [0u8; 2];
    File::open(path)
        .and_then(|mut f| f.read_exact(&mut magic))
        .is_ok_and(|_| magic == [0x1f, 0x8b])
}

/// Read the .fai index of a FASTA file.
fn read_fai(path: &Path) -> Result<fasta::fai::Index> {
    let index_path = path.with_extension("fa.fai");
    let index_path = if index_path.exists() {
        index_path
    } else {
        let alt = format!("{}.fai", path.display());
        Path::new(&alt).to_path_buf()
    };
    fasta::fai::read(&index_path)
        .with_context(|| format!("failed to read FASTA index: {}", index_path.display()))
}

/// Contig names and lengths of a FASTA file, in .fai order.
pub fn fasta_contigs(fasta_path: &str) -> Result<Vec<(String, u64)>> {
    let index = read_fai(Path::new(fasta_path))?;
    Ok(index
        .as_ref()
        .iter()
        .map(|r| {
            let name: &[u8] = r.name();
            (String::from_utf8_lossy(name).into_owned(), r.length())
        })
        .collect())
}

/// Indexed reference FASTA reader with region caching.
struct ReferenceReader {
    index: fasta::fai::Index,
    reader: fasta::IndexedReader<fasta::io::BufReader<File>>,
    cache: HashMap<(String, u64, u64), Vec<u8>>,
    cache_order: VecDeque<(String, u64, u64)>,
    max_cache_entries: usize,
}

impl ReferenceReader {
    /// Open a reference FASTA with its .fai index.
    ///
    /// Detects a bgzipped FASTA (`.gz`/`.bgz`) by extension and decompresses
    /// it via its `.gzi` index; a plain FASTA is read as raw bytes as
    /// before. Reading a bgzipped FASTA as raw bytes (the previous bug here)
    /// silently returns garbage/truncated sequence instead of failing, which
    /// then surfaces later as a misleading "beyond chromosome length" error.
    ///
    /// noodles itself only looks at the extension to decide bgzf vs. raw, so
    /// a bgzip-compressed FASTA under any other name (`.fa`, `.fna`, ...)
    /// would still hit that same bug. Peek the gzip magic bytes ourselves
    /// and fail naming the real cause -- we don't attempt to *handle* such a
    /// file, since it needs a `.gzi` index that, under a name samtools/bgzip
    /// never produced, almost certainly isn't sitting next to it.
    fn open<P: AsRef<Path>>(fasta_path: P) -> Result<Self> {
        let path = fasta_path.as_ref();
        let index = read_fai(path)?;

        let is_bgzip_named = matches!(
            path.extension().and_then(|e| e.to_str()),
            Some("gz") | Some("bgz")
        );
        if !is_bgzip_named && looks_gzip_compressed(path) {
            anyhow::bail!(
                "{} is gzip-compressed (starts with the gzip magic bytes 1f 8b) but is not \
                 named .gz/.bgz, so it would be read as raw uncompressed sequence; rename it \
                 to end in .gz or .bgz with a matching .gzi index, or decompress it first",
                path.display()
            );
        }

        let reader = fasta::io::indexed_reader::Builder::default()
            .set_index(index.clone())
            .build_from_path(path)
            .with_context(|| {
                if is_bgzip_named {
                    format!(
                        "failed to open FASTA reader for {} (if this is a bgzipped FASTA, \
                         it needs a matching .gzi index alongside it)",
                        path.display()
                    )
                } else {
                    format!("failed to open FASTA reader for {}", path.display())
                }
            })?;

        Ok(Self {
            index,
            reader,
            cache: HashMap::new(),
            cache_order: VecDeque::new(),
            max_cache_entries: 256,
        })
    }

    /// Fetch a reference sequence region.
    /// `start` and `end` are 0-based, half-open coordinates [start, end).
    /// Returns the bases as uppercase ASCII bytes (A, C, G, T, N).
    fn fetch_sequence(&mut self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        let key = (chrom.to_string(), start, end);
        if let Some(seq) = self.cache.get(&key) {
            return Ok(seq.clone());
        }

        // Find the record in the index
        let _idx_record = self
            .index
            .as_ref()
            .iter()
            .find(|r| {
                let name_bytes: &[u8] = r.name();
                name_bytes == chrom.as_bytes()
            })
            .with_context(|| format!("chromosome '{}' not found in FASTA index", chrom))?;

        // noodles uses 1-based, closed coordinates for the region
        let start_usize = usize::try_from(start).context("start position exceeds platform usize")?;
        let end_usize = usize::try_from(end).context("end position exceeds platform usize")?;
        let noodles_start = noodles::core::Position::try_from(start_usize + 1)
            .context("invalid start position")?;
        let noodles_end =
            noodles::core::Position::try_from(end_usize).context("invalid end position")?;

        let region = noodles::core::Region::new(chrom, noodles_start..=noodles_end);

        let record = self
            .reader
            .query(&region)
            .with_context(|| format!("failed to query FASTA region {}:{}-{}", chrom, start, end))?;

        let seq: Vec<u8> = record
            .sequence()
            .as_ref()
            .iter()
            .map(|&b| b.to_ascii_uppercase())
            .collect();

        // Cache management (simple FIFO eviction)
        if self.cache.len() >= self.max_cache_entries {
            if let Some(oldest) = self.cache_order.pop_front() {
                self.cache.remove(&oldest);
            }
        }
        self.cache_order.push_back(key.clone());
        self.cache.insert(key, seq.clone());

        Ok(seq)
    }

    /// Get the length of a chromosome from the FASTA index.
    fn chromosome_length(&self, chrom: &str) -> Result<u64> {
        let record = self
            .index
            .as_ref()
            .iter()
            .find(|r| {
                let name_bytes: &[u8] = r.name();
                name_bytes == chrom.as_bytes()
            })
            .with_context(|| format!("chromosome '{}' not found in FASTA index", chrom))?;

        Ok(record.length())
    }
}

/// Thread-safe, read-only reference sequence store.
///
/// Pre-loads entire chromosome sequences into memory so that all threads
/// can share a single copy via `&SharedReference` (no Arc/Mutex needed
/// since it's immutable after construction).
///
/// Memory: ~64 MB per chromosome (e.g. chr20), ~3 GB for full human genome.
pub struct SharedReference {
    sequences: HashMap<String, Vec<u8>>,
}

impl SharedReference {
    /// Load specified chromosomes from a FASTA file into memory.
    pub fn load(fasta_path: &str, chromosomes: &[&str]) -> Result<Self> {
        let mut reader = ReferenceReader::open(fasta_path)?;
        let mut sequences = HashMap::new();

        for chrom in chromosomes {
            let len = reader.chromosome_length(chrom)?;
            let seq = reader.fetch_sequence(chrom, 0, len)?;
            sequences.insert(chrom.to_string(), seq);
        }

        let total_mb: usize = sequences.values().map(|s| s.len()).sum::<usize>() / (1024 * 1024);
        log::info!(
            "Loaded {} chromosome(s) into shared reference ({} MB)",
            sequences.len(),
            total_mb
        );

        Ok(Self { sequences })
    }

    /// Get the length of a loaded chromosome, or None if not loaded.
    pub fn chromosome_length(&self, chrom: &str) -> Option<u64> {
        self.sequences.get(chrom).map(|s| s.len() as u64)
    }

    /// Create a SharedReference from pre-built in-memory sequences (for testing).
    #[cfg(test)]
    pub fn from_sequences(sequences: HashMap<String, Vec<u8>>) -> Self {
        Self { sequences }
    }

    /// Fetch a reference sequence region.
    /// `start` and `end` are 0-based, half-open coordinates [start, end).
    pub fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        let seq = self
            .sequences
            .get(chrom)
            .with_context(|| format!("chromosome '{}' not in shared reference", chrom))?;
        let start_usize = usize::try_from(start)
            .with_context(|| format!("start position {} exceeds platform usize for {}", start, chrom))?;
        let end_usize = usize::try_from(end)
            .with_context(|| format!("end position {} exceeds platform usize for {}", end, chrom))?;
        let start_clamped = start_usize.min(seq.len());
        let end_clamped = end_usize.min(seq.len());
        if start_clamped >= end_clamped {
            return Ok(Vec::new());
        }
        Ok(seq[start_clamped..end_clamped].to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Write a one-contig bgzipped FASTA plus its .fai and .gzi companions to
    /// `dir`, using noodles' own bgzf writer/gzi types (no external tools).
    /// `file_name` is the FASTA's own name within `dir` (its extension
    /// controls whether `ReferenceReader::open` treats it as bgzipped).
    fn write_bgzipped_fasta(
        dir: &Path,
        file_name: &str,
        contig: &str,
        sequence: &[u8],
    ) -> std::path::PathBuf {
        let fa_path = dir.join(file_name);

        let mut raw = Vec::new();
        raw.extend_from_slice(format!(">{contig}\n").as_bytes());
        let offset = raw.len() as u64;
        raw.extend_from_slice(sequence);
        raw.push(b'\n');

        let mut writer = noodles::bgzf::Writer::new(Vec::new());
        writer.write_all(&raw).expect("write bgzf data");
        let compressed = writer.finish().expect("finish bgzf stream");
        std::fs::write(&fa_path, &compressed).expect("write FASTA");

        let record = fasta::fai::Record::new(
            contig,
            sequence.len() as u64,
            offset,
            sequence.len() as u64,
            sequence.len() as u64 + 1,
        );
        fasta::fai::fs::write(
            dir.join(format!("{file_name}.fai")),
            &fasta::fai::Index::from(vec![record]),
        )
        .expect("write .fai");

        // The whole record fits in one bgzf block, so every position we'll
        // query here truly is at compressed offset 0 -- an empty gzi::Index
        // happens to be correct for THIS fixture, but only because it is
        // single-block. It is not a general fallback: for a file with more
        // than one block, an empty index gives the *wrong* answer for any
        // position past the first block (see the two-block fixture below,
        // which builds a real, non-empty .gzi and checks that case).
        noodles::bgzf::gzi::fs::write(
            dir.join(format!("{file_name}.gzi")),
            &noodles::bgzf::gzi::Index::default(),
        )
        .expect("write .gzi");

        fa_path
    }

    #[test]
    fn fetch_sequence_reads_bgzipped_fasta_correctly() {
        let dir = std::env::temp_dir().join(format!(
            "spike_test_reference_bgzip_{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let fa_path = write_bgzipped_fasta(&dir, "ref.fa.gz", "chr1", b"ACGTACGTAC");

        let mut reader = ReferenceReader::open(&fa_path).expect("open bgzipped FASTA");
        let seq = reader
            .fetch_sequence("chr1", 0, 10)
            .expect("fetch_sequence on bgzipped FASTA");
        assert_eq!(seq, b"ACGTACGTAC");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn open_bgzipped_fasta_without_gzi_names_the_real_problem() {
        let dir = std::env::temp_dir().join(format!(
            "spike_test_reference_bgzip_missing_gzi_{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let fa_path = write_bgzipped_fasta(&dir, "ref.fa.gz", "chr1", b"ACGTACGTAC");
        std::fs::remove_file(dir.join("ref.fa.gz.gzi")).unwrap();

        let message = match ReferenceReader::open(&fa_path) {
            Ok(_) => panic!("expected opening a bgzipped FASTA without a .gzi to fail"),
            Err(err) => format!("{err:#}").to_lowercase(),
        };
        assert!(
            message.contains("gzi"),
            "error should name the missing bgzip index (.gzi), got: {message}"
        );
        assert!(
            !message.contains("beyond chromosome length"),
            "error should not be the misleading downstream message, got: {message}"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    /// `ReferenceReader::open` decides bgzf-vs-raw by extension, same as
    /// noodles. A bgzip-compressed FASTA named `.fa` (no `.gz`/`.bgz`) must
    /// not be silently read as raw sequence bytes -- it must fail, naming
    /// the real cause, not the downstream "beyond chromosome length" error.
    #[test]
    fn open_bgzip_compressed_fasta_under_a_fa_extension_fails_with_clear_message() {
        let dir = std::env::temp_dir().join(format!(
            "spike_test_reference_bgzip_wrong_ext_{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        // Genuinely bgzf-compressed content, named as a plain FASTA would be.
        let fa_path = write_bgzipped_fasta(&dir, "ref.fa", "chr1", b"ACGTACGTAC");

        let message = match ReferenceReader::open(&fa_path) {
            Ok(_) => panic!(
                "expected opening a bgzip-compressed FASTA named .fa to fail \
                 instead of silently reading it as raw sequence"
            ),
            Err(err) => format!("{err:#}").to_lowercase(),
        };
        assert!(
            message.contains("gzip") || message.contains("bgzip"),
            "error should name the real cause (gzip-compressed content), got: {message}"
        );
        assert!(
            !message.contains("beyond chromosome length"),
            "error should not be the misleading downstream message, got: {message}"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    /// `write_bgzipped_fasta` above always produces a single bgzf block, so
    /// its empty `.gzi` never exercises `gzi::Index::query`'s non-trivial
    /// branch (a position past the first block). Build a genuinely
    /// two-block bgzf file plus a real (non-empty) `.gzi`, and fetch bases
    /// that live strictly inside the second block -- not at its very first
    /// byte, so a reader that just kept reading sequentially past a
    /// mis-sought block1 would not accidentally land on the right answer.
    #[test]
    fn fetch_sequence_seeks_correctly_past_a_real_block_boundary() {
        let dir = std::env::temp_dir().join(format!(
            "spike_test_reference_bgzip_two_block_{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();

        let contig = "chr1";
        let header = format!(">{contig}\n");
        let seq_offset = header.len() as u64;
        // Non-repeating bases in each block: a seek that lands even one
        // base off gives a *different* string, not a coincidentally
        // identical one (repeated bases like all-A/all-C would hide that).
        let first_block_seq = b"ACGTACGTAC"; // sequence bases 0..10, block 1
        let second_block_seq = b"TGCATGCATG"; // sequence bases 10..20, block 2
        let fa_path = dir.join("ref.two_block.fa.gz");

        let mut writer = noodles::bgzf::Writer::new(Vec::new());
        writer.write_all(header.as_bytes()).expect("write header");
        writer.write_all(first_block_seq).expect("write block 1 bases");
        // Force block 1 to close here, so block 2 starts at a real,
        // non-zero compressed offset.
        writer.flush().expect("close block 1");
        let block2_compressed_pos = writer.position();
        let block2_uncompressed_pos = seq_offset + first_block_seq.len() as u64;
        writer.write_all(second_block_seq).expect("write block 2 bases");
        writer.write_all(b"\n").expect("write trailing newline");
        let compressed = writer.finish().expect("finish bgzf stream");
        std::fs::write(&fa_path, &compressed).expect("write FASTA");

        let seq_len = (first_block_seq.len() + second_block_seq.len()) as u64;
        let record = fasta::fai::Record::new(contig, seq_len, seq_offset, seq_len, seq_len + 1);
        fasta::fai::fs::write(
            dir.join("ref.two_block.fa.gz.fai"),
            &fasta::fai::Index::from(vec![record]),
        )
        .expect("write .fai");

        // A real, non-empty .gzi: one entry recording exactly where block 2
        // begins, in both compressed and uncompressed coordinates.
        let gzi = noodles::bgzf::gzi::Index::from(vec![(
            block2_compressed_pos,
            block2_uncompressed_pos,
        )]);
        noodles::bgzf::gzi::fs::write(dir.join("ref.two_block.fa.gz.gzi"), &gzi)
            .expect("write .gzi");

        let mut reader = ReferenceReader::open(&fa_path).expect("open two-block bgzipped FASTA");
        // Sequence bases 12..20 are 2 bases into block 2, not at its start:
        // a seek that lands at block 2's start instead (the sequential
        // fallback a broken/empty gzi would produce) reads "TGCATGCA", not
        // this.
        let seq = reader
            .fetch_sequence(contig, 12, 20)
            .expect("fetch_sequence across the block boundary");
        assert_eq!(seq, b"CATGCATG");

        std::fs::remove_dir_all(&dir).ok();
    }
}
