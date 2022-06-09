# Load the ribowaltz package.
library(riboWaltz)

# Load in the fly GTF annotation file.
df = create_annotation(gtfpath = "../Datasets/reference_files/annotation.gtf")

# ribowaltz requires a file that has been aligned to the transcriptome rather than 
# a file that has been aligned to the genome. 

# Load up bam_list from transcriptome.
bam_list = bamtolist("../Datasets/testing_transcriptome_alignments",
                     annotation = df)

# Calculate P_site offset.
offsets = psite(data = bam_list, start = FALSE, extremity = "3end")

# Filter the offsets dataframe to only include the information needed by plastid. 
samples = unique(offsets$sample)

for (i in 1:length(samples)) {
  sam_offs = subset(offsets, sample == samples[i], 
                    select = c("length", "corrected_offset_from_3"))
  colnames(sam_offs) = c("length", "p_offset")
  write.table(sam_offs, paste("../Datasets/testing_Psite_offsets/",
                              samples[i],"_p-site-offsets", sep=""),sep = "\t" , row.names = FALSE)
  
}
