# abridge

`abridge` is a software program to compress short-read alignments. 

# Installation

`abridge` is released as a docker image and can be obtained from DockerHub.

```bash
docker pull abridge:1.0.0
```

On High performance clusters please use singularity 


# Generating alignments



## STAR

# Compress

## Compress a single RNA-Seq file

## Nitty gritties of compression

Explain the process and refer to paper and supplemental information

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

For all other queries please email sagnik@iastate.edu