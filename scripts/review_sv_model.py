#!/usr/bin/env python3
"""Reproduce the model findings in CLINICAL_SV_REVIEW.md on synthetic data.

Requires a built spike binary, Python 3.9+, and samtools on PATH. The optional
--with-alignment also requires bwa-mem2. This is a diagnostic, not a clinical
benchmark or an assertion that a caller should pass/fail. It writes fresh data
only; --output must not already exist. No patient data or downloads are used.
"""

import argparse
import gzip
import json
import random
import subprocess
import tempfile
from pathlib import Path


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def fastq(path):
    with gzip.open(path, "rt") as stream:
        while name := stream.readline().strip():
            sequence = stream.readline().strip()
            stream.readline()
            stream.readline()
            yield name[1:].removesuffix("/1").removesuffix("/2"), sequence


class Probe:
    def __init__(self, spike, output):
        self.spike = str(spike.resolve())
        self.output = output
        self.results = {}
        rng = random.Random(571)
        self.reference = "".join(rng.choices("ACGT", k=40000))
        self.fasta = output / "ref.fa"
        self.fasta.write_text(">chrT\n" + self.reference + "\n")
        subprocess.run(["samtools", "faidx", str(self.fasta)], check=True)
        self.index = {
            self.reference[i:i + 31]: i
            for i in range(len(self.reference) - 30)
        }
        assert len(self.index) == len(self.reference) - 30

    def locate(self, sequence):
        # Unique random reference; sample several exact seeds in either strand.
        # Used for interior depth windows, never for breakpoint accuracy.
        for oriented in (sequence, reverse_complement(sequence)):
            for offset in (0, 50, 100, 119):
                position = self.index.get(oriented[offset:offset + 31])
                if position is not None:
                    return position - offset
        return None

    def make_bam(self, label, variable=False, lowmap=False, homdel=False):
        sam = self.output / (label + ".sam")
        records = []
        marker = 18000
        haplotype = self.reference
        if homdel:
            haplotype = haplotype[:marker] + haplotype[marker + 2:]
        stop = 39496 if homdel else 39500
        with sam.open("w") as stream:
            stream.write(
                "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrT\tLN:40000\n"
                "@RG\tID:test\tSM:test\n"
            )
            for index, start in enumerate(range(100, stop, 4)):
                if variable and 16000 <= start < 25000 and start % 16:
                    continue
                name = f'{"m" if homdel else "p"}{index:06}'
                mapq = 0 if lowmap and index % 2 == 0 else 60
                starts = [start, start + 250]
                positions = [s + (2 if homdel and s >= marker else 0) for s in starts]
                end = start + 400 + (2 if homdel and start + 400 > marker else 0)
                records.append((name, positions[0], end, mapq))
                for mate, (flag, hap_start) in enumerate(zip((99, 147), starts)):
                    cigar = "150M"
                    if homdel and hap_start < marker < hap_start + 150:
                        cigar = f"{marker - hap_start}M2D{hap_start + 150 - marker}M"
                    length = (end - positions[0]) * (1 if mate == 0 else -1)
                    fields = [
                        name, flag, "chrT", positions[mate] + 1, mapq, cigar,
                        "=", positions[1 - mate] + 1, length,
                        haplotype[hap_start:hap_start + 150], "]" * 150, "RG:Z:test",
                    ]
                    stream.write("\t".join(map(str, fields)) + "\n")
        bam = self.output / (label + ".bam")
        subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
        subprocess.run(["samtools", "index", str(bam)], check=True)
        return bam, records

    def run(self, label, bam, events=(), extra=()):
        destination = self.output / label
        command = [
            self.spike, "--bam", str(bam), "--reference", str(self.fasta),
            "--flank", "2000", "--seed", "17", "-o", str(destination),
        ]
        for event in events:
            command += ["--event", event]
        command.extend(extra)
        with (self.output / (label + ".log")).open("w") as log:
            subprocess.run(command, check=True, stdout=log, stderr=log)
        return destination

    def depths(self, destination, records, windows):
        # Reconstruct merge.sh's retained originals plus the emitted FASTQs.
        replaced = set((destination / "replaced_reads.txt").read_text().splitlines())
        totals = [0] * len(windows)

        def add(position):
            if position is not None:
                for index, (low, high) in enumerate(windows):
                    totals[index] += max(0, min(position + 150, high) - max(position, low))

        for name, start, end, _ in records:
            if name not in replaced:
                add(start)
                add(end - 150)
        for mate in (1, 2):
            for _, sequence in fastq(destination / f"R{mate}.fq.gz"):
                add(self.locate(sequence))
        return [round(t / (hi - lo), 4) for t, (lo, hi) in zip(totals, windows)]

    def model_probes(self):
        uniform, records = self.make_bam("uniform")
        variable, variable_records = self.make_bam("variable", variable=True)
        lowmap, lowmap_records = self.make_bam("lowmap", lowmap=True)
        windows = [(10200, 10800), (12200, 12800)]
        for fraction in (0.5, 1):
            for adjacent in (False, True):
                events = [f"del:chrT:10000-11000;af={fraction}"]
                if adjacent:
                    events += [f"del:chrT:12000-13000;af={fraction}"]
                label = f'del_{fraction}_{"adjacent" if adjacent else "single"}'
                destination = self.run(label, uniform, events)
                self.results[label] = {
                    "depths": self.depths(destination, records, windows), "baseline_depth": 75,
                }
        for label, bam, donor in (
            ("dup_uniform", uniform, records), ("dup_variable", variable, variable_records),
        ):
            destination = self.run(label, bam, ["dup:chrT:10000-28000;af=0.5"])
            self.results[label] = {
                "depths": self.depths(destination, donor, [(12000, 14000), (19000, 21000)]),
                "baseline_depths": [75, 18.75 if label == "dup_variable" else 75],
            }
        destination = self.run("del_lowmap", lowmap, ["del:chrT:10000-14000;af=1"])
        self.results["del_lowmap"] = {
            "depth": self.depths(destination, lowmap_records, [(11000, 13000)])[0],
            "baseline_depth": 75,
        }
        for size in (1000, 10000, 100000, 1000000):
            label = f"ins_{size}"
            destination = self.run(label, uniform, [f"ins:chrT:20000:{size};af=0.5"])
            synthetic = [
                (first, second)
                for (name, first), (_, second) in zip(
                    fastq(destination / "R1.fq.gz"), fastq(destination / "R2.fq.gz"),
                ) if name.startswith("ev")
            ]
            anchored = sum(
                self.locate(first) is not None or self.locate(second) is not None
                for first, second in synthetic
            )
            self.results[label] = {
                "synthetic_pairs": len(synthetic), "pairs_with_reference_31mer": anchored,
            }

        bam, _ = self.make_bam("homdel", homdel=True)
        destination = self.run("dup_homdel", bam, ["dup:chrT:10000-28000;af=0.5"])
        alt = self.reference[17985:18000] + self.reference[18002:18017]
        ref = self.reference[17985:18017]
        counts = {"alt": 0, "ref": 0}
        for mate in (1, 2):
            for _, sequence in fastq(destination / f"R{mate}.fq.gz"):
                for key, allele in (("alt", alt), ("ref", ref)):
                    if allele in sequence or reverse_complement(allele) in sequence:
                        counts[key] += 1
        self.results["dup_homdel"] = {
            "marker_read_counts": counts, "original_marker_AF": 1,
            "output_marker_AF": counts["alt"] / sum(counts.values()),
        }

        vcf = self.output / "hom_input.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.3\n##contig=<ID=chrT,length=40000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tINPUT\n"
            f"chrT\t10000\thomdel\t{self.reference[9999]}\t<DEL>\t.\tPASS\t"
            "SVTYPE=DEL;END=11000\tGT\t1/1\n"
        )
        destination = self.run("vcf_hom", uniform, extra=["--vcf", str(vcf)])
        self.results["vcf_hom"] = {
            "input_GT": "1/1",
            "output_records": [line for line in (destination / "truth.vcf").read_text().splitlines()
                               if not line.startswith("#")],
        }
        replaced = set((self.output / "del_1_single/replaced_reads.txt").read_text().splitlines())
        self.results["mate_recovery"] = {
            "pair": "p003200", "read_starts_0based": [12900, 13150],
            "query_0based_half_open": [8000, 13000],
            "in_replaced_names": "p003200" in replaced,
        }

    def alignment_probes(self):
        with (self.output / "index.log").open("w") as log:
            subprocess.run(["bwa-mem2", "index", str(self.fasta)],
                           stdout=log, stderr=log, check=True)
        results = {}
        for label, windows in (
            ("del_1_single", [(10200, 10800)]),
            ("del_1_adjacent", [(10200, 10800), (12200, 12800)]),
            ("dup_variable", [(12000, 14000), (19000, 21000)]),
            ("del_lowmap", [(11000, 13000)]),
        ):
            destination = self.output / label
            with (destination / "pipeline.log").open("w") as log:
                for script in ("align.sh", "merge.sh"):
                    subprocess.run(["bash", str(destination / script)],
                                   stdout=log, stderr=log, check=True)
            depths = []
            for low, high in windows:
                lines = subprocess.check_output([
                    "samtools", "depth", "-aa", "-r", f"chrT:{low + 1}-{high}",
                    str(destination / "merged.bam"),
                ], text=True).splitlines()
                depths.append(sum(int(line.split()[2]) for line in lines) / (high - low))
            result = subprocess.run([
                self.spike, "validate", "--bam", str(destination / "merged.bam"),
                "--truth", str(destination / "truth.vcf"),
                "--reference", str(self.fasta), "--json",
            ], capture_output=True, text=True)
            # Nonzero is expected: the controlled fixed-length donor has SD=0,
            # and these convenience scripts have not marked duplicates.
            results[label] = {
                "aligned_depths": depths, "validate_exit": result.returncode,
                "validation": json.loads(result.stdout),
            }
        (self.output / "alignment-summary.json").write_text(json.dumps(results, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spike", type=Path,
                        default=Path(__file__).resolve().parents[1] / "target/debug/spike")
    parser.add_argument("--output", type=Path, help="Fresh output directory; defaults to a new /tmp directory")
    parser.add_argument("--with-alignment", action="store_true")
    args = parser.parse_args()
    if args.output:
        output = args.output.resolve()
        output.mkdir(parents=True, exist_ok=False)
    else:
        output = Path(tempfile.mkdtemp(prefix="spike-sv-review-"))
    probe = Probe(args.spike, output)
    probe.model_probes()
    (output / "summary.json").write_text(json.dumps(probe.results, indent=2) + "\n")
    if args.with_alignment:
        probe.alignment_probes()
    print(json.dumps(probe.results, indent=2))
    print(f"Artifacts: {output}")


if __name__ == "__main__":
    main()
