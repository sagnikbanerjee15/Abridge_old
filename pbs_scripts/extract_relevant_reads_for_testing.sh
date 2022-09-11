# Align the raw reads to the reference genome
# Select specific reads accoring to the commands below
# Create new fastq file with new sequence

FILE_NAME_BASE=/project/maizegdb/sagnik/ABRIDGE/testing/alignments/SRR13711353_PE_all_tags_sorted

# Extract headers
cat ${FILE_NAME_BASE}.sam | grep ^@ > ${FILE_NAME_BASE}_headers.sam

# Extract non-spliced alignments that are a perfect match to the reference
cat ${FILE_NAME_BASE}.sam | awk '$6 !~ /N/' | grep MD:Z:151 | shuf -n 10000 > ${FILE_NAME_BASE}_non_spliced_perfect.sam

# Extract spliced alignments that are a perfect match to the reference
cat ${FILE_NAME_BASE}.sam | awk '$6 ~ /N/' | grep MD:Z:151 | shuf -n 10000 > ${FILE_NAME_BASE}_SE_spliced_perfect.sam