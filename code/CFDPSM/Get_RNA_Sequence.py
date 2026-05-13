from Bio import SeqIO
import os

# ============================================================================
# CONFIGURATION
# ============================================================================
task_dir = ""      # Your task directory
genome = "hg38"
input_path = ""    # Your input VCF file
context = 70

print(f'[INFO] New task starts at {task_dir}')
os.makedirs(task_dir, exist_ok=True)


# ============================================================================
# FUNCTIONS
# ============================================================================
def fasta_dna_to_rna(task_dir):
    """
    Convert DNA FASTA sequences to RNA (replace T with U)
    """
    ref_input = os.path.join(task_dir, "ref_sequences.fasta")
    alt_input = os.path.join(task_dir, "alt_sequences.fasta")
    ref_output = os.path.join(task_dir, "ref_sequences_rna.fasta")
    alt_output = os.path.join(task_dir, "alt_sequences_rna.fasta")
    
    # Convert reference sequences
    with open(ref_output, "w") as out_handle:
        for record in SeqIO.parse(ref_input, "fasta"):
            rna_seq = str(record.seq).replace("T", "U")
            out_handle.write(f">{record.id}\n{rna_seq}\n")
    
    # Convert alternative sequences
    with open(alt_output, "w") as out_handle:
        for record in SeqIO.parse(alt_input, "fasta"):
            rna_seq = str(record.seq).replace("T", "U")
            out_handle.write(f">{record.id}\n{rna_seq}\n")
    
    print("[INFO] RNA FASTA files generated")


# ============================================================================
# MAIN EXECUTION
# ============================================================================
vcf_to_bed(input_path, context, task_dir)
bed_to_fasta(f"{task_dir}/input.bed", task_dir)
retrieve_sequence(input_path, task_dir, context)
fasta_dna_to_rna(task_dir)