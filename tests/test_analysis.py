import logging
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd
import pyensembl

from scripts import ase_cleanup
from scripts.ase_cleanup import allelic_imbalance, het_test
from scripts.gene_annotation import create_gene_table, match_annotation_contig
from scripts.hap_table import (
    filter_genotypes,
    get_genotypes,
    haplotype_imbalance,
    select_sample,
    write_empty_table,
)
from scripts.validate_inputs import (
    validate_alignment_contigs,
    validate_vcf_contigs,
)


class AseCleanupTests(unittest.TestCase):
    def test_binomial_tests_use_supported_scipy_api(self):
        ase = pd.DataFrame({
            "refCount": [10, 18],
            "altCount": [10, 2],
            "totalCount": [20, 20],
            "ref_bias": [0.5, 0.5],
        })

        result = het_test(allelic_imbalance(ase), perror=0.01)

        self.assertTrue(result["AEI_pval"].notna().all())
        self.assertTrue(result["AEI_padj"].notna().all())
        self.assertTrue(result["het_padj"].notna().all())

    def test_missing_reference_bias_is_filled(self):
        ase = pd.DataFrame({
            "contig": ["1"],
            "position": [10],
            "variantID": ["."],
            "refAllele": ["A"],
            "altAllele": ["G"],
            "refCount": [2],
            "altCount": [18],
            "totalCount": [20],
            "lowMAPQDepth": [0],
            "lowBaseQDepth": [0],
            "rawDepth": [20],
            "otherBases": [0],
            "improperPairs": [0],
        })
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_path = root / "ase.tsv"
            output_path = root / "clean.tsv"
            ase.to_csv(input_path, sep="\t", index=False)
            previous_directory = os.getcwd()
            try:
                os.chdir(root)
                with mock.patch.object(sys, "argv", [
                        "ase_cleanup.py", "--ase", str(input_path),
                        "--output", str(output_path)]):
                    ase_cleanup.main()
            finally:
                os.chdir(previous_directory)
                logging.shutdown()
                for handler in logging.root.handlers[:]:
                    logging.root.removeHandler(handler)
                    handler.close()

            result = pd.read_csv(output_path, sep="\t")

        self.assertEqual(result.loc[0, "ref_bias"], 0.5)


class GeneAnnotationTests(unittest.TestCase):
    def test_contig_aliases_match_annotation_style(self):
        self.assertEqual(match_annotation_contig("1", {"chr1", "chrX"}), "chr1")
        self.assertEqual(match_annotation_contig("chr1", {"1", "X"}), "1")
        self.assertEqual(match_annotation_contig("MT", {"chrM"}), "chrM")

    def test_gene_lookup_accepts_chr_alias(self):
        rows = [
            'chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "GENE1"; gene_name "Gene1"; gene_biotype "protein_coding";',
            'chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "Gene1"; transcript_name "Tx1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";',
            'chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; exon_id "EX1"; exon_number "1";',
        ]
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            gtf = root / "annotation.gtf"
            gtf.write_text("\n".join(rows) + "\n")
            genome = pyensembl.Genome(
                reference_name="GRCh37",
                annotation_name="test",
                gtf_path_or_url=str(gtf),
                cache_directory_path=str(root))
            genome.index()

            result = create_gene_table(
                genome, pd.DataFrame({"contig": ["1"], "position": [50]}))

            self.assertEqual(result.iloc[0]["gene_id"], "GENE1")
            genome.db._connection.close()


class HaplotypeTests(unittest.TestCase):
    def test_genotype_filter_keeps_only_phased_diploid_hets(self):
        genotypes = pd.DataFrame({"GT": ["0/1", "0|1", "1|0", "1|1", None]})

        result = filter_genotypes(genotypes)

        self.assertEqual(result["GT"].tolist(), ["0|1", "1|0"])

    def test_hse_ratio_uses_gene_level_counts(self):
        genes = pd.DataFrame({
            "gene_id": ["g1", "g2"],
            "hap1": [30, 3],
            "hap2": [20, 7],
            "totalCount": [50, 10],
            "positions": [2, 1],
        })

        haplotype_imbalance(genes)

        self.assertEqual(genes["HSE_ratio"].tolist(), [0.6, 0.3])

    def test_sample_selection_is_unambiguous(self):
        self.assertEqual(select_sample(["only"]), "only")
        self.assertEqual(select_sample(["a", "b"], bam_name="b"), "b")
        self.assertEqual(select_sample(["a", "b"], requested_sample="a"), "a")
        with self.assertRaisesRegex(ValueError, "multiple samples"):
            select_sample(["a", "b"], bam_name="unknown")

    def test_vcf_reader_uses_selected_sample(self):
        vcf = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2
1\t10\t.\tA\tG\t.\tPASS\t.\tGT:DP\t0|1:20\t1|0:18
"""
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "test.vcf"
            path.write_text(vcf)

            result = get_genotypes(path, requested_sample="s2")

        self.assertEqual(result.to_dict("records"), [
            {"position": 10, "contig": "1", "GT": "1|0"}])

    def test_empty_haplotype_table_is_still_created(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "empty.tsv"

            write_empty_table(path)

            self.assertTrue(path.is_file())
            self.assertIn("HSE_ratio", path.read_text().splitlines()[0])


class InputValidationTests(unittest.TestCase):
    def test_matching_contigs_are_accepted(self):
        reference = {"chr1": 200, "chr2": 300}
        alignment = {"chr1": 200, "chr2": 300}

        validate_alignment_contigs(alignment, reference)
        validate_vcf_contigs(
            {"chr1": 200}, {"chr1"}, alignment, reference
        )

    def test_alignment_reference_length_mismatch_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "contig length mismatch"):
            validate_alignment_contigs({"chr1": 199}, {"chr1": 200})

    def test_vcf_record_contig_missing_from_alignment_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "absent from the alignment"):
            validate_vcf_contigs(
                {"chr2": 300}, {"chr2"}, {"chr1": 200},
                {"chr1": 200, "chr2": 300}
            )


if __name__ == "__main__":
    unittest.main()
