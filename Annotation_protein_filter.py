import argparse
import os
import tempfile
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO

def process_chunk(records, standard_aa):
    kept, removed = [], []
    for record in records:
        seq = str(record.seq).upper()
        if set(seq).issubset(standard_aa):
            kept.append(record)
        else:
            removed.append((record.id, seq))
    return kept, removed

def length_filter(input_file, output_file, min_len, max_len, threads):
    cmd = f"seqkit seq -g -j {threads} -m {min_len} -M {max_len} {input_file} > {output_file}"
    subprocess.run(cmd, shell=True, check=True)

def cd_hit_cluster(input_file, output_file, identity, threads):
    cmd = f"cd-hit -i {input_file} -o {output_file} -c {identity} -n 5 -M 0 -T {threads}"
    subprocess.run(cmd, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description="Comprehensive protein sequence processing pipeline")
    parser.add_argument("-i", "--input", default="input.faa", help="Input FASTA file")
    parser.add_argument("-o", "--output", default="processed.faa", help="Final output file")
    parser.add_argument("-r", "--report", default="removed.txt", help="File to record removed sequence IDs")
    parser.add_argument("-t", "--threads", type=int, default=os.cpu_count(), help="Number of threads")
    parser.add_argument("-c", "--chunk", type=int, default=50000, help="Number of records per chunk")
    parser.add_argument("--keep-uo", action="store_true", help="Keep U and O amino acids")
    parser.add_argument("--min-len", type=int, default=100, help="Minimum sequence length")
    parser.add_argument("--max-len", type=int, default=1500, help="Maximum sequence length")
    parser.add_argument("--identity", type=float, default=0.90, help="CD-HIT identity threshold")
    parser.add_argument("--skip-length", action="store_true", help="Skip length filtering")
    parser.add_argument("--skip-cdhit", action="store_true", help="Skip CD-HIT clustering")
    args = parser.parse_args()

    standard_aa = set("ACDEFGHIKLMNPQRSTVWY")
    if args.keep_uo:
        standard_aa.update("UO")

    total_kept = total_removed = 0
    unique_ids = duplicate_ids = 0
    temp_files = []
    id_count = {}
    
    step_counts = {}

    with open(args.report, 'w'):
        pass

    print("Step 1: Removing non-standard amino acids...")
    with tempfile.TemporaryDirectory() as tmpdir:
        step1_output = os.path.join(tmpdir, "step1_cleaned.fasta")
        
        input_count = 0
        
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            futures = []
            chunk = []
            for record in SeqIO.parse(args.input, "fasta"):
                input_count += 1
                chunk.append(record)
                if len(chunk) >= args.chunk:
                    futures.append(executor.submit(process_chunk, chunk, standard_aa))
                    chunk = []
            if chunk:
                futures.append(executor.submit(process_chunk, chunk, standard_aa))

            for future in as_completed(futures):
                kept, removed = future.result()
                total_kept += len(kept)
                total_removed += len(removed)

                if kept:
                    tmp_kept = tempfile.NamedTemporaryFile(mode='w', delete=False, dir=tmpdir, suffix='.fasta')
                    SeqIO.write(kept, tmp_kept, 'fasta')
                    tmp_kept.close()
                    temp_files.append(tmp_kept.name)

                if removed:
                    with open(args.report, 'a') as report_handle:
                        for rid, seq in removed:
                            report_handle.write(f"{rid}\t{seq}\n")

        with open(step1_output, 'w') as out_handle:
            for tmp_file in temp_files:
                with open(tmp_file) as f:
                    out_handle.write(f.read())
                os.unlink(tmp_file)

        step_counts['input_sequences'] = input_count
        step_counts['after_aa_filter'] = total_kept
        step_counts['removed_nonstandard_aa'] = total_removed
        current_file = step1_output
        
        if not args.skip_length:
            print("Step 2: Length filtering...")
            step2_output = os.path.join(tmpdir, "step2_length_filtered.fasta")
            length_filter(current_file, step2_output, args.min_len, args.max_len, args.threads)
            
            step_counts['after_length_filter'] = sum(1 for _ in SeqIO.parse(step2_output, "fasta"))
            step_counts['removed_by_length'] = step_counts['after_aa_filter'] - step_counts['after_length_filter']
            current_file = step2_output
        else:
            step_counts['after_length_filter'] = step_counts['after_aa_filter']
            step_counts['removed_by_length'] = 0

        if not args.skip_cdhit:
            print("Step 3: CD-HIT clustering...")
            step3_output = os.path.join(tmpdir, "step3_clustered.fasta")
            cd_hit_cluster(current_file, step3_output, args.identity, args.threads)
            
            cluster_file = f"{step3_output}.clstr"
            output_cluster_file = f"{args.output}.clstr"
            if os.path.exists(cluster_file):
                with open(cluster_file) as src, open(output_cluster_file, 'w') as dst:
                    dst.write(src.read())
            
            step_counts['after_cdhit'] = sum(1 for _ in SeqIO.parse(step3_output, "fasta"))
            step_counts['removed_by_cdhit'] = step_counts['after_length_filter'] - step_counts['after_cdhit']
            current_file = step3_output
        else:
            step_counts['after_cdhit'] = step_counts['after_length_filter']
            step_counts['removed_by_cdhit'] = 0

        print("Step 4: Renaming duplicate IDs...")
        step4_output = os.path.join(tmpdir, "step4_renamed.fasta")
        id_count = {}
        with open(step4_output, 'w') as out_handle:
            for record in SeqIO.parse(current_file, "fasta"):
                main_id = record.id
                if main_id in id_count:
                    id_count[main_id] += 1
                    record.id = f"{main_id}_{id_count[main_id]}"
                    record.description = record.id
                    duplicate_ids += 1
                else:
                    id_count[main_id] = 0
                    unique_ids += 1
                SeqIO.write(record, out_handle, 'fasta')
        current_file = step4_output

        with open(current_file) as src, open(args.output, 'w') as dst:
            dst.write(src.read())

    print("\n" + "="*60)
    print("PROCESSING SUMMARY")
    print("="*60)
    print(f"Input sequences: {step_counts['input_sequences']}")
    print()
    print("Step-by-step filtering:")
    print(f"├─ Removed non-standard AA: {step_counts['removed_nonstandard_aa']}")
    print(f"├─ After AA filtering: {step_counts['after_aa_filter']}")
    if not args.skip_length:
        print(f"├─ Removed by length filter: {step_counts['removed_by_length']}")
        print(f"├─ After length filtering: {step_counts['after_length_filter']}")
    if not args.skip_cdhit:
        print(f"├─ Removed by CD-HIT clustering: {step_counts['removed_by_cdhit']}")
        print(f"├─ After CD-HIT clustering: {step_counts['after_cdhit']}")
    print(f"├─ Duplicate IDs renamed: {duplicate_ids}")
    print(f"└─ Final sequences: {unique_ids + duplicate_ids}")
    print()
    print("Output files:")
    print(f"├─ Final output: {args.output}")
    if not args.skip_cdhit:
        print(f"├─ CD-HIT cluster file: {args.output}.clstr")
    print(f"└─ Removed records: {args.report}")
    print("="*60)

if __name__ == "__main__":
    main()
