#加载ftrCOOL包
library(ftrCOOL)
library(openxlsx)


outDirPath <- YOUR_OUTPUT_DIR # Modify this 
ref_filePath <- paste0(outDirPath, '/', 'ref_sequences_rna.fasta')
alt_filePath <- paste0(outDirPath, '/', 'alt_sequences_rna.fasta')

ref_outFileMat <- paste0(outDirPath, '/', 'ExpectedValKmerNUC_RNA_ref.xlsx')
alt_outFileMat <- paste0(outDirPath, '/', 'ExpectedValKmerNUC_RNA_alt.xlsx')

ref_mat <- ExpectedValKmerNUC_RNA(seqs=ref_filePath, k=4, ORF=FALSE, reverseORF=FALSE)
alt_mat <- ExpectedValKmerNUC_RNA(seqs=alt_filePath, k=4, ORF=FALSE, reverseORF=FALSE)

write.table(ref_mat, ref_outFileMat, row.names = FALSE, col.names = TRUE)
write.table(alt_mat, alt_outFileMat, row.names = FALSE, col.names = TRUE)