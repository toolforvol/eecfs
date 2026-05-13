Your_app_path="" # Where you install the ifeatureOmega
Your_seq_path="" # The sequence input directory
Your_task_path="" # The feature output directory

##### DNA #####
python get_iFeatureOmega_EECFS.py \
--ref_fasta_file ${Your_seq_path}/DNA_seq/ref_sequences.fasta \
--alt_fasta_file ${Your_seq_path}}/DNA_seq/ref_sequences.fasta \
--output_file ${Your_task_path}/ifo.fea.DNA.csv \
--molecular DNA \
--target_features 'SCPseDNC' 'NCP' 'EIIP' 'RCKmer type 1' \
--selected_features 'GAA_diff' 'ACT_diff' 'AGC_diff' 'SCPseDNC_lamada_1_diff' 'NCP_164_alt' 'EIIP_64_alt'

##### RNA #####
python get_iFeatureOmega_EECFS.py \
--ref_fasta_file ${Your_seq_path}/RNA_seq/ref_sequences.fasta \
--alt_fasta_file ${Your_seq_path}/RNA_seq/ref_sequences.fasta \
--output_file ${Your_task_path}/RNA_seq/ifo.fea.RNA.csv \
--molecular RNA \
--target_features 'NCP' \
--selected_features 'NCP_226_alt'

##### protein #####
python get_iFeatureOmega_EECFS.py \
--ref_fasta_file ${Your_seq_path}/protein_seq/ref_sequences.fasta \
--output_file ${Your_task_path}/protein_seq/ifo.fea.protein.csv \
--molecular protein \
--target_features 'AAIndex' 'BLOSUM62' 'ZScale' \
--selected_features 'AAindex_p.9.ANDN920101' 'AAindex_p.11.BEGF750103' 'blosum62_172' 'blosum62_195' 'blosum62_198' 'blosum62_219' 'ZScale_p9.z2'