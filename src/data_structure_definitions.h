#ifndef _ABRIDGE_SRC_DATA_STRUCTURE_DEFINITIONS_H_
#define _ABRIDGE_SRC_DATA_STRUCTURE_DEFINITIONS_H_

/*
 * Defining macros
 */
# define FILENAME_LENGTH 10000
# define ROWS 10000
# define COLS 10000
# define MAX_SEQ_LEN 1000
# define MAX_CIGAR_ITEMS 5000
# define MAX_CIGAR_LENGTH 10000
# define MAX_TAG_LENGTH 1000
# define MAX_POOL_SIZE 1000000
# define MIN_POOL_SIZE 10000
# define MAX_UNIQUE_CIGARS 10000
# define MAX_LINE_TO_BE_WRITTEN_TO_FILE 10000000
# define MAX_LENGTH_RLE 1000
# define MAX_CIGAR_FREQ_SIZE 1000000
# define MAX_ICIGAR_LENGTH 100000
# define MAX_ICIGAR_LENGTH_PASS1_COL2 10000000
# define MAX_REFERENCE_SEQUENCES 300000
# define MAX_SYMBOLS_FOR_PASS3_COMPRESSION 44
# define MAX_SYMBOLIC_ICIGAR_LENGTH 3
# define MAX_REFERENCE_SEQ_LENGTH 1000
# define MAX_GENERAL_LEN 1000000
# define MAX_BUFFER_SIZE_FOR_READING_PASS2_FILE 1073741824
# define MAX_REFERENCE_SEQ_LEN 1000000000
# define MAX_FILES_FOR_MERGING 1000
# define MAX_READ_ID_LENGTH 1000

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
		char *reference_name;
};

struct Sam_Tags
{
		/*
		 * Template to store tag information for each alignment
		 */
		char *name; // Name of the sam tag. Could be MD, NH, etc.
		char *type; // Datatype of the tag
		char *val; // Value of the tag. Must be a string even if tag has a single character
};

struct Soft_Clippings
{
		/*
		 * Template to store soft clipped portion of reads
		 */
		char *left;
		char *right;
		char *left_qual;
		char *right_qual;
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
		char *cigar;
		char *icigar;
		char **pointers_to_qual_scores;
		long int num_reads;
		unsigned long long int position;
		char **pointers_to_read_names;
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
		long long int position;
};

struct Pass3_Compression_Symbol_icigar_Mapping
{
		char *symbolic_icigar;
		char *icigar;
};

struct Chromosome_Info
{
		char **name;
		long long int *length;
		int number_of_chromosomes;
};

struct Chromosome_Starting_Byte
{
		char **name;
		int number_of_chromosomes;
		long int *start_byte_in_pass2_file;
};

struct Merged_Compressed_DS
{
		char *col1;
		char *col2;
		char **icigars;
		int **number_of_reads;
		long long int position;
		int number_of_unique_cigars;
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
		long long int template_length; //observed Template length
		char *seq; //The read sequence
		char *qual; //The read quality

		/*
		 * Optional fields of variable length
		 */
		//struct Sam_Tags tags[100];
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
		char *soft_clips_removed_qual;
		char *selected_qual;
		int soft_clips_removed_seq_len;

		char *NH;
		char *AS;
		char *MD;

};

struct RLE_Quality_Scores
{
		long long int frequency;
		char quality_score;
};

struct Abridge_Index
{
		char **chromosome;
		long long int *start;
		long long int *end;
		long long int *start_byte;
		long long int *end_byte;
		long long int number_of_items;
		long long int *start_byte_qual;
		long long int *end_byte_qual;
};

struct Whole_Genome_Sequence
{
		char **reference_sequence_name;
		char **nucleotides;
		unsigned long long int *reference_sequence_length;
		int number_of_reference_sequences;
};

struct Quality_Score_RLE
{
		char score_character;
		long long int frequency;
};

struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list
{
		struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *prev;
		char *old_read_id;
		char *new_read_id;
		int number_of_multi_maps;
		int valid;
		struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *next;
};

struct All_Relevant_Info_PE_per_Alignment
{
		char *new_read_id;
		char *old_read_id;
		char *icigar;
		int NH_value;
		unsigned long long int position;
		short int new_read_id_assigned;
};

struct Paired_Ended_Flag_to_Single_Character
{
		char *character;
		int *samflags;
};

struct Read_Ids_to_NH
{
		char *read_id;
		unsigned int NH;
};

int isCharacterInString( char*, char );
void splitCigar( char*, int*, struct Cigar_Items* );
long long int extractNHfromicigar( char*, int );
void generateReadSequenceAndMDString(
		struct Sam_Alignment*,
		struct Whole_Genome_Sequence* );
int findSamFormatFlag( char*, int, char* );

#endif /* _ABRIDGE_SRC_DATA_STRUCTURE_DEFINITIONS_H_ */
