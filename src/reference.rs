use anyhow::{Context, Result};
use noodles::fasta;
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::path::Path;

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
    fn open<P: AsRef<Path>>(fasta_path: P) -> Result<Self> {
        let path = fasta_path.as_ref();
        let index = read_fai(path)?;

        let reader = fasta::io::indexed_reader::Builder::default()
            .set_index(index.clone())
            .build_from_path(path)
            .with_context(|| {
                format!(
                    "failed to open FASTA reader for {} (if this is a bgzipped FASTA, \
                     it needs a matching .gzi index alongside it)",
                    path.display()
                )
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
    fn write_bgzipped_fasta(
        dir: &Path,
        contig: &str,
        sequence: &[u8],
    ) -> std::path::PathBuf {
        let fa_path = dir.join("ref.fa.gz");

        let mut raw = Vec::new();
        raw.extend_from_slice(format!(">{contig}\n").as_bytes());
        let offset = raw.len() as u64;
        raw.extend_from_slice(sequence);
        raw.push(b'\n');

        let mut writer = noodles::bgzf::Writer::new(Vec::new());
        writer.write_all(&raw).expect("write bgzf data");
        let compressed = writer.finish().expect("finish bgzf stream");
        std::fs::write(&fa_path, &compressed).expect("write .fa.gz");

        let record = fasta::fai::Record::new(
            contig,
            sequence.len() as u64,
            offset,
            sequence.len() as u64,
            sequence.len() as u64 + 1,
        );
        fasta::fai::fs::write(
            dir.join("ref.fa.gz.fai"),
            &fasta::fai::Index::from(vec![record]),
        )
        .expect("write .fai");

        // The whole record fits in one bgzf block, so the block boundary
        // list is empty; gzi::Index::query() falls back to offset 0 for any
        // position in that case, which is exactly right here.
        noodles::bgzf::gzi::fs::write(dir.join("ref.fa.gz.gzi"), &noodles::bgzf::gzi::Index::default())
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
        let fa_path = write_bgzipped_fasta(&dir, "chr1", b"ACGTACGTAC");

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
        let fa_path = write_bgzipped_fasta(&dir, "chr1", b"ACGTACGTAC");
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
}
