#ifndef _ABRIDGE_SRC_DATA_STRUCTURE_DEFINITIONS_H_
#define _ABRIDGE_SRC_DATA_STRUCTURE_DEFINITIONS_H_

/*
 * Defining macros
 */
# define FILENAME_LENGTH 10000
# define ROWS 10000
# define COLS 10000
# define MAX_SEQ_LEN 2000
# define MAX_CIGAR_ITEMS 5000
# define MAX_CIGAR_LENGTH 10000
# define MAX_TAG_LENGTH 25
# define MAX_POOL_SIZE 1000000
# define MIN_POOL_SIZE 500
# define MAX_UNIQUE_CIGARS 10000
# define MAX_LINE_TO_BE_WRITTEN_TO_FILE 1000000
# define MAX_LENGTH_RLE 1000
# define MAX_CIGAR_FREQ_SIZE 1000000
# define MAX_ICIGAR_LENGTH 10000
# define MAX_ICIGAR_LENGTH_PASS1_COL2 1000000
# define MAX_REFERENCE_SEQUENCES 10000
# define MAX_SYMBOLS_FOR_PASS3_COMPRESSION 44
# define MAX_SYMBOLIC_ICIGAR_LENGTH 3
# define MAX_REFERENCE_SEQ_LENGTH 1000
# define MAX_GENERAL_LEN 1000000
# define MAX_BUFFER_SIZE_FOR_READING_PASS2_FILE 1073741824

char *sam_tags[] =
{ "NH", "HI", "AS", "nM", "NM", "MD", "jM", "jI", "XS", "RG", "XT", "XM", "AM", "X0", "X1", "XO", "XG" };
// ATGCN - ASCII Codes 33-37
char insert_characters[] =
{ '!', '"', '#', '$', '%', '\0' };
// ATGCN - ASCII Codes 38-42
char mismatch_characters[] =
{ '&', '\'', '(', ')', '*', '\0' };

struct Reference_Sequence_Info
{
		char *line;
};

struct Sam_Tags
{
		/*
		 * Template to store tag information for each alignment
		 */
		char *name; // Name of the sam tag. Could be MD, NH, etc.
		char *type; // Datatype of the tag
		char *val; // Value of the tag. Must be a string even ig tag has a single character
};

struct Soft_Clippings
{
		/*
		 * Template to store soft clipped portion of reads
		 */
		char *left;
		char *right;
};

struct Cigar_Frequency
{
		/*
		 * Template to store the CIGARs and the number of times they appear
		 */
		char *cigar;
		long int freq;
};

struct Cigar_Items
{
		/*
		 * Template to store cigar items
		 */
		char def;
		int len;
};

struct Compressed_DS
{
		/*
		 * Template to store the compressed icigars
		 */
		char *icigar;
		long int num_reads;
		long long int position;
};

struct Pass1_Compressed_DS
{
		/*
		 * Template to store data from 1st pass
		 */
		char *col1;
		char *col2;
		char *col3;
};

struct Pass2_Compressed_DS
{
		/*
		 * Template to store data from 1st pass
		 */
		char *col1;
		char *col2;
};

struct Pass3_Compression_Symbol_icigar_Mapping
{
		char *symbolic_icigar;
		char *icigar;
};

struct Sam_Alignment
{
		/*
		 * Template to store information related to each alignment in sam file
		 */

		//Fields are listed in the order of appearance in SAM file
		/*
		 * Mandatory fields
		 */
		char *read_name; // Name of the read
		int samflag; // Flag produced by aligner. For more details see https://www.samformat.info/sam-format-flag
		char *reference_name; //Chromosome name
		long long int start_position; // starting position of the read
		int mapping_quality_score;
		char *cigar;
		char *reference_name_next_mate; //Reference name where the mate/next read is mapped
		long long int start_position_next; //Position of the mate/next read
		int template_length; //observed Template length
		char *seq; //The read sequence
		char *qual; //The read quality

		/*
		 * Optional fields of variable length
		 */
		struct Sam_Tags tags[100];

		/*
		 * Other fields created by abridge
		 */
		struct Soft_Clippings soft_clippings; // Stores the left and the right soft clippings
		struct Cigar_Items *cigar_items;
		int number_of_cigar_items;
		int number_of_tag_items;
		int read_seq_len; // basically strlen(seq). Field is kept to prevent extra calls to strlen

		char *temp; //For temporary operations
		char *cigar_extended;
		char *md_extended;
		char *icigar; // Stores the integrated representation comprising of all relevant information about the alignment
		char **splices;
		char *soft_clips_removed_seq;
		int soft_clips_removed_seq_len;

};

struct Abridge_Index
{
		char **chromosome;
		long long int *start;
		long long int *end;
		long long int *start_byte;
		long long int *end_byte;
		//long long int *line_number_start;
		//long long int *line_number_end;
		long long int number_of_items;
};

struct Whole_Genome_Sequence
{
		char **reference_sequence_name;
		char **nucleotides;
		unsigned long long int *reference_sequence_length;
		int number_of_reference_sequences;
};

#endif /* _ABRIDGE_SRC_DATA_STRUCTURE_DEFINITIONS_H_ */
