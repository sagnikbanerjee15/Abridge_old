# CHANGELOG

## Version 1.1.0

- Quality scores are compressed using FCLQC. Quality scores of unmapped files are now being processed similarly.
- Brotli and 7z are removed. All compressions are carried out by ZPAQ.
- ABRIDGE can now accept BAM and SAM files. Input files may be unsorted.
- Having alignment tags are no longer mandatory. ABRIDGE will generate tags if need be. 
- Requirement for XS tag has also been removed. Presence of XS tag will be ignored.
- ABRIDGE will perform shortening of read names unless otherwise requested by the user.

## Version 1.0.0
- First release