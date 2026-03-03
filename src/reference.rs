use anyhow::{Context, Result};
use noodles::fasta;
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Indexed reference FASTA reader with region caching.
struct ReferenceReader {
    index: fasta::fai::Index,
    reader: fasta::IndexedReader<BufReader<File>>,
    cache: HashMap<(String, u64, u64), Vec<u8>>,
    cache_order: VecDeque<(String, u64, u64)>,
    max_cache_entries: usize,
}

impl ReferenceReader {
    /// Open a reference FASTA with its .fai index.
    fn open<P: AsRef<Path>>(fasta_path: P) -> Result<Self> {
        let path = fasta_path.as_ref();

        let index_path = path.with_extension("fa.fai");
        let index_path = if index_path.exists() {
            index_path
        } else {
            let alt = format!("{}.fai", path.display());
            Path::new(&alt).to_path_buf()
        };

        let index = fasta::fai::read(&index_path)
            .with_context(|| format!("failed to read FASTA index: {}", index_path.display()))?;

        let file = File::open(path)
            .with_context(|| format!("failed to open FASTA: {}", path.display()))?;
        let reader = fasta::IndexedReader::new(BufReader::new(file), index.clone());

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
        let noodles_start = noodles::core::Position::try_from((start as usize) + 1)
            .context("invalid start position")?;
        let noodles_end =
            noodles::core::Position::try_from(end as usize).context("invalid end position")?;

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
        let start = (start as usize).min(seq.len());
        let end = (end as usize).min(seq.len());
        if start >= end {
            return Ok(Vec::new());
        }
        Ok(seq[start..end].to_vec())
    }
}
