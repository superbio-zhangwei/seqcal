import concurrent.futures
import os
import shutil
import subprocess
import uuid
from datetime import datetime

import pandas as pd
from tqdm import tqdm


IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)

REGION_COLUMNS = [
    "seqID",
    "start",
    "end",
    "gene_position",
    "left_primer_start",
    "left_primer_end",
    "right_primer_start",
    "right_primer_end",
]


def reverse_complement(sequence):
    return sequence.translate(IUPAC_COMPLEMENT)[::-1]


def read_fasta(fasta_path):
    records = {}
    seq_id = None
    chunks = []

    with open(fasta_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    records[seq_id] = "".join(chunks)
                seq_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if seq_id is not None:
        records[seq_id] = "".join(chunks)

    return records


def wrap_sequence(sequence, width=80):
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def read_locate_table(path):
    if os.path.getsize(path) == 0:
        return pd.DataFrame(columns=["seqID", "patternName", "pattern", "strand", "start", "end", "matched"])
    return pd.read_csv(path, sep="\t")


def pair_hits(left_hits, right_hits, left_inner_col, right_inner_col, gene_position):
    left = left_hits.copy()
    right = right_hits.copy()
    if left.empty or right.empty:
        return pd.DataFrame()

    left["hit"] = left[left_inner_col]
    left["primer_side"] = "left"
    right["hit"] = right[right_inner_col]
    right["primer_side"] = "right"

    ordered = pd.concat([left, right], ignore_index=True).sort_values(by="hit").reset_index(drop=True)
    pair_mask = (ordered["primer_side"] == "left") & (ordered["primer_side"].shift(-1) == "right")

    pairs = []
    for idx in ordered[pair_mask].index:
        left_row = ordered.loc[idx]
        right_row = ordered.loc[idx + 1]
        pairs.append(
            {
                "seqID": left_row["seqID"],
                "start": int(left_row["hit"]),
                "end": int(right_row["hit"]),
                "gene_position": gene_position,
                "left_primer_start": int(left_row["start"]),
                "left_primer_end": int(left_row["end"]),
                "right_primer_start": int(right_row["start"]),
                "right_primer_end": int(right_row["end"]),
            }
        )

    return pd.DataFrame(pairs)


def find_target_regions(df_F, df_R):
    seq_ids = sorted(set(df_F.get("seqID", pd.Series(dtype=str))).union(df_R.get("seqID", pd.Series(dtype=str))))
    region_tables = []

    for seq_id in seq_ids:
        f_seq = df_F[df_F["seqID"] == seq_id]
        r_seq = df_R[df_R["seqID"] == seq_id]

        plus_regions = pair_hits(
            f_seq[f_seq["strand"] == "+"],
            r_seq[r_seq["strand"] == "-"],
            left_inner_col="end",
            right_inner_col="start",
            gene_position="+",
        )
        minus_regions = pair_hits(
            r_seq[r_seq["strand"] == "+"],
            f_seq[f_seq["strand"] == "-"],
            left_inner_col="end",
            right_inner_col="start",
            gene_position="-",
        )

        region_tables.extend([plus_regions, minus_regions])

    region_tables = [table for table in region_tables if not table.empty]
    if not region_tables:
        return pd.DataFrame(columns=REGION_COLUMNS)

    return pd.concat(region_tables, ignore_index=True)


def extract_region_sequence(sequence, region, include_primers):
    if include_primers:
        start = int(region["left_primer_start"])
        end = int(region["right_primer_end"])
    else:
        start = int(region["left_primer_end"]) + 1
        end = int(region["right_primer_start"]) - 1

    extracted = sequence[start - 1 : end]
    if region["gene_position"] == "-":
        extracted = reverse_complement(extracted)

    return start, end, extracted


def pass_length_filter(sequence_length, min_length, max_length):
    if min_length is not None and sequence_length < min_length:
        return False
    if max_length is not None and sequence_length > max_length:
        return False
    return True


def make_match_dirs(output_dir, temp_dir):
    temp_base = temp_dir or os.path.join(output_dir, "temp")
    run_id = datetime.now().strftime("seqcal_run_%Y%m%d_%H%M%S") + f"_{uuid.uuid4().hex[:8]}"
    run_temp_dir = os.path.join(temp_base, run_id)
    return (
        os.path.join(run_temp_dir, "primerF_matches"),
        os.path.join(run_temp_dir, "primerR_matches"),
    )


def check_seqkit():
    if shutil.which("seqkit"):
        return
    raise FileNotFoundError(
        "seqkit was not found in PATH. Please activate the environment that contains seqkit before running seqcal."
    )


def configure_pandas():
    try:
        pd.set_option("future.no_silent_downcasting", True)
    except (KeyError, ValueError, pd.errors.OptionError):
        pass


def discover_fna_files(fna_input):
    if os.path.isdir(fna_input):
        fna_files = [f for f in os.listdir(fna_input) if f.endswith((".fna", ".fasta"))]
        return [os.path.join(fna_input, f) for f in sorted(fna_files)]
    if os.path.isfile(fna_input) and fna_input.endswith((".fna", ".fasta")):
        return [fna_input]
    raise ValueError(f"Invalid input: {fna_input}. Must be a directory or a .fna/.fasta file.")


def run_seqcal(
    fna_input,
    temp_dir,
    primerF_output_dir,
    primerR_output_dir,
    output_dir,
    forward_primer,
    reverse_primer,
    threads,
    sequence_output,
    regions_output,
    include_primers,
    min_length,
    max_length,
):
    configure_pandas()

    if min_length is not None and min_length < 0:
        raise ValueError("--min-length must be greater than or equal to 0.")
    if max_length is not None and max_length < 0:
        raise ValueError("--max-length must be greater than or equal to 0.")
    if min_length is not None and max_length is not None and min_length > max_length:
        raise ValueError("--min-length cannot be greater than --max-length.")

    check_seqkit()

    os.makedirs(primerF_output_dir, exist_ok=True)
    os.makedirs(primerR_output_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, "copies.seqcal")
    sequence_output_path = os.path.join(output_dir, sequence_output)
    regions_output_path = os.path.join(output_dir, regions_output)
    with open(output_path, "w", encoding="utf-8") as output_file:
        output_file.write("genome\tcopies\n")
    open(sequence_output_path, "w", encoding="utf-8").close()

    fna_files = discover_fna_files(fna_input)
    fna_by_sample = {os.path.splitext(os.path.basename(path))[0]: path for path in fna_files}

    def match_primer_forward(fna_file):
        base_name = os.path.splitext(os.path.basename(fna_file))[0]
        output_file = os.path.join(primerF_output_dir, f"primerF_{base_name}.txt")
        cmd = ["seqkit", "locate", "-d", "-i", "-p", forward_primer, fna_file]
        with open(output_file, "w", encoding="utf-8") as outfile:
            subprocess.run(cmd, stdout=outfile, check=True)

    def match_primer_reverse(fna_file):
        base_name = os.path.splitext(os.path.basename(fna_file))[0]
        output_file = os.path.join(primerR_output_dir, f"primerR_{base_name}.txt")
        cmd = ["seqkit", "locate", "-d", "-i", "-p", reverse_primer, fna_file]
        with open(output_file, "w", encoding="utf-8") as outfile:
            subprocess.run(cmd, stdout=outfile, check=True)

    print("Forward primer matching...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        list(tqdm(executor.map(match_primer_forward, fna_files), total=len(fna_files), desc="Processing primerF"))

    print("Reverse primer matching...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        list(tqdm(executor.map(match_primer_reverse, fna_files), total=len(fna_files), desc="Processing primerR"))

    print("Copy number calculating and sequence extracting...")
    all_region_rows = []

    for sample in fna_by_sample:
        primerF_file = f"primerF_{sample}.txt"
        primerR_file = f"primerR_{sample}.txt"
        primerF_path = os.path.join(primerF_output_dir, primerF_file)
        primerR_path = os.path.join(primerR_output_dir, primerR_file)
        if not os.path.exists(primerF_path):
            raise FileNotFoundError(f"Missing forward primer match file: {primerF_path}")
        if not os.path.exists(primerR_path):
            raise FileNotFoundError(f"Missing reverse primer match file: {primerR_path}")

        df_F = read_locate_table(primerF_path)
        df_R = read_locate_table(primerR_path)
        regions = find_target_regions(df_F, df_R)
        copy_number = len(regions)

        with open(output_path, "a", encoding="utf-8") as output_file:
            output_file.write(f"{sample}\t{copy_number}\n")

        if copy_number == 0:
            continue

        fasta_records = read_fasta(fna_by_sample[sample])
        with open(sequence_output_path, "a", encoding="utf-8") as fasta_out:
            kept_index = 0
            for _, region in regions.reset_index(drop=True).iterrows():
                seq_id = region["seqID"]
                if seq_id not in fasta_records:
                    raise KeyError(f"Sequence ID {seq_id} was found by seqkit but not in {fna_by_sample[sample]}")

                extract_start, extract_end, target_seq = extract_region_sequence(
                    fasta_records[seq_id],
                    region,
                    include_primers,
                )
                if not pass_length_filter(len(target_seq), min_length, max_length):
                    continue

                kept_index += 1
                target_id = f"{sample}_target_{kept_index}"
                fasta_out.write(
                    f">{target_id} sample={sample} seqID={seq_id} "
                    f"strand={region['gene_position']} region={extract_start}-{extract_end} "
                    f"length={len(target_seq)}\n{wrap_sequence(target_seq)}\n"
                )

                region_row = region.to_dict()
                region_row.update(
                    {
                        "sample": sample,
                        "target_id": target_id,
                        "extract_start": extract_start,
                        "extract_end": extract_end,
                        "extract_length": len(target_seq),
                        "include_primers": include_primers,
                    }
                )
                all_region_rows.append(region_row)

    region_columns = [
        "sample",
        "target_id",
        "seqID",
        "gene_position",
        "start",
        "end",
        "left_primer_start",
        "left_primer_end",
        "right_primer_start",
        "right_primer_end",
        "extract_start",
        "extract_end",
        "extract_length",
        "include_primers",
    ]
    pd.DataFrame(all_region_rows, columns=region_columns).to_csv(regions_output_path, sep="\t", index=False)
