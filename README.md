To do:

1. Docker and singularity support
2. Test it out with singularity on AWS
3. 



# ABRIDGE

Advancement in technology has enabled sequencing machines to produce vast amounts of genetic data, causing an increase in storage demands. Most genomic software utilizes read alignments for several purposes including transcriptome assembly and gene count estimation. Herein we present, ABRIDGE, a state-of-the-art compressor for SAM alignment files offering users both lossless and lossy compression options. This reference-based  file compressor achieves the best compression ratio among all compression software ensuring lower space demand and faster file transmission. Central to the software is a novel algorithm that retains non-redundant information. This new approach has allowed ABRIDGE to achieve a compression 16% higher than the second-best compressor for RNA-Seq reads and over 35% for DNA-Seq reads.  ABRIDGE also offers users the option to randomly access location without having to decompress the entire file. ABRIDGE is distributed under MIT license and can be obtained from GitHub and docker hub.  We anticipate that the user community will adopt ABRIDGE within their existing pipeline encouraging further research in this domain.

`abridge` is a software program to compress short-read alignments.  

# Installation

`abridge` is released through GitHub. Download the latest release from GutHub

```bash
wget https://github.com/sagnikbanerjee15/Abridge/archive/refs/tags/ABRIDGE_v1.0.0.tar.gz
tar -xvzf ABRIDGE_v1.0.0.tar.gz
cd ABRIDGE_v1.0.0
echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc
source ~/.bashrc
```




```bash
https://github.com/sagnikbanerjee15/Abridge.git
```

On High performance clusters please use singularity.

# Generating alignments

Alignments can be generated using a number of different software. 

## STAR



# Compress

## Compress a single RNA-Seq file

### Lossless compression



### Lossy compression





## Generic compression tools



## Diagnostic mode



# Decompress





# Compare compressed and decompressed files



# Frequently asked questions (FAQ)

1. I have several RNA-Seq samples aligned to a reference which is no longer available. What should I do?
2. I aligned my RNA-Seq samples to a reference without the MD tag. How do I generate MD tag for `abrigde` compression?
3. How do I know that the compression is lossless?
4. I need to view only the references and not alignments. What command should I execute?
5. What is the correct command for viewing alignments without reference headers?
6. How do I retreive alignments to a particular chromosome?
7. I have only bamfiles but `abridge` needs samfiles to compress. What should I do?
8. How do I retrieve multiple locations using random access?
9. Can I generate overlapping read coverage for random locations?



# Future upgrades

Here is a list of future upgrades to `abridge`

1. Retrieve alignments from multiple locations from compressed file

# License

License information can be accessed from 

# Contact

Please list all issues at https://github.com/sagnikbanerjee15/Abridge/issues

For all other queries please email sagnikbanerjee15@gmail.com 