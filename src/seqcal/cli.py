import argparse

from .core import make_match_dirs, run_seqcal


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="seqcal",
        description="seqcal: primer-guided sequence extraction and copy-number calculation tool",
    )
    parser.add_argument("-f", "--fna-input", type=str, required=True, help="Directory or single file (.fna/.fasta)")
    parser.add_argument("-o", "--output-dir", type=str, required=True, help="Directory for output results")
    parser.add_argument(
        "--temp-dir",
        type=str,
        default=None,
        help="Base directory for temporary primer match files; defaults to OUTPUT_DIR/temp",
    )
    parser.add_argument(
        "--forward-primer",
        type=str,
        required=True,
        help="Forward primer sequence (e.g., GATGAAGARCGYAGYRAA)",
    )
    parser.add_argument(
        "--reverse-primer",
        type=str,
        required=True,
        help="Reverse primer sequence (e.g., TCCTCCGCTTATTGATATGC)",
    )
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use for concurrent tasks")
    parser.add_argument(
        "--sequence-output",
        type=str,
        default="seqhit.fasta",
        help="FASTA file name for extracted target sequences",
    )
    parser.add_argument(
        "--regions-output",
        type=str,
        default="target_regions.tsv",
        help="TSV file name for extracted target coordinates",
    )
    parser.add_argument("--include-primers", action="store_true", help="Extract target sequences including primer sites")
    parser.add_argument(
        "--min-length",
        type=int,
        default=None,
        help="Only extract target sequences with length greater than or equal to this value",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=None,
        help="Only extract target sequences with length less than or equal to this value",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    primerF_output_dir, primerR_output_dir = make_match_dirs(
        output_dir=args.output_dir,
        temp_dir=args.temp_dir,
    )
    run_seqcal(
        fna_input=args.fna_input,
        temp_dir=args.temp_dir,
        primerF_output_dir=primerF_output_dir,
        primerR_output_dir=primerR_output_dir,
        output_dir=args.output_dir,
        forward_primer=args.forward_primer,
        reverse_primer=args.reverse_primer,
        threads=args.threads,
        sequence_output=args.sequence_output,
        regions_output=args.regions_output,
        include_primers=args.include_primers,
        min_length=args.min_length,
        max_length=args.max_length,
    )


if __name__ == "__main__":
    main()
