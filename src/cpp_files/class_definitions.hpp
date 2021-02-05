/*
 * class_definitions.hpp
 *
 *  Created on: Jan 30, 2021
 *      Author: sagnik
 */

#ifndef ABRIDGE_CLASS_DEFINITIONS_HPP_
#define ABRIDGE_CLASS_DEFINITIONS_HPP_
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
# define MAX_NUM_TAG 500
# define MAX_POOL_SIZE 100000
# define MAX_UNIQUE_CIGARS 10000
# define MAX_LINE_TO_BE_WRITTEN_TO_FILE 100000
# define MAX_LENGTH_RLE 1000
# define MAX_CIGAR_FREQ_SIZE 1000000
# define MAX_ICIGAR_LENGTH 100000
# define MAX_ICIGAR_LENGTH_PASS1_COL2 1000000
# define MAX_REFERENCE_SEQUENCES 10000
# define MAX_SYMBOLS_FOR_PASS3_COMPRESSION 44
# define MAX_SYMBOLIC_ICIGAR_LENGTH 3
# define MAX_REFERENCE_SEQ_LENGTH 1000
# define MAX_GENERAL_LEN 1000000

# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <iostream>

char *sam_tags[] =
{ "NH", "HI", "AS", "nM", "NM", "MD", "jM", "jI", "XS", "RG", "XT", "XM", "AM", "X0", "X1", "XO", "XG" };
// ATGCN - ASCII Codes 33-37
char insert_characters[] =
{ '!', '"', '#', '$', '%' };
// ATGCN - ASCII Codes 38-42
char mismatch_characters[] =
{ '&', '\'', '(', ')', '*' };

/*********************************************************************
 * Function definitions
 *********************************************************************/

int splitByDelimiter(char*, char, char**);
void extractSubString(char*, char*, int, int);
int isCharacterInString(char*, char);
void insertCharacterInString(char*, int*, char, int, int);
/********************************************************************/

class Sam_Tags
{
		/*
		 * Class to store tag information for each alignment
		 */
	public:
		char **name; // Name of the sam tag. Could be MD, NH, etc.
		char **type; // Datatype of the tag
		char **val; // Value of the tag. Must be a string even ig tag has a single character
		int number_of_tags;

		Sam_Tags()
		{
			name = new char*[MAX_NUM_TAG];
			for (int i = 0; i < MAX_NUM_TAG; i++)
				name[i] = new char[MAX_TAG_LENGTH];
			type = new char*[MAX_NUM_TAG];
			for (int i = 0; i < MAX_NUM_TAG; i++)
				type[i] = new char[MAX_TAG_LENGTH];
			val = new char*[MAX_NUM_TAG];
			for (int i = 0; i < MAX_NUM_TAG; i++)
				val[i] = new char[MAX_TAG_LENGTH];
			number_of_tags = 0;
		}
		~Sam_Tags()
		{

		}

		int locateSamTags(char *tag)
		{
			/*
			 * Returns the location of the requested samtag
			 * Returns -1 if that tag is absent
			 */
			for (int i = 0; i < number_of_tags; i++)
			{
				if (strcmp(name[i], tag) == 0) return i;
			}
			return -1;

		}

};

class Soft_Clippings
{
	public:
		char *left;
		char *right;

		Soft_Clippings()
		{
			left = new char[MAX_SEQ_LEN];
			right = new char[MAX_SEQ_LEN];
		}
		~Soft_Clippings()
		{

		}
};

class Cigar_Items
{
	public:
		char *def;
		long long int *len;
		int number_of_cigar_items;

		Cigar_Items()
		{
			def = new char[MAX_CIGAR_ITEMS];
			len = new long long int[MAX_CIGAR_ITEMS];
			number_of_cigar_items = 0;
		}
		~Cigar_Items()
		{

		}
};

class Reference_Sequence_Info
{
	public:
		char **line;
		int number_of_reference_sequences;
		Reference_Sequence_Info()
		{
			line = new char*[MAX_REFERENCE_SEQUENCES];
			for (int i = 0; i < MAX_REFERENCE_SEQUENCES; i++)
				line[i] = new char[MAX_GENERAL_LEN];
			number_of_reference_sequences = 0;
		}
		~Reference_Sequence_Info()
		{

		}

};

class Compressed_DS
{
	public:
		char *icigar;
		long int num_reads;
		long long int position;

		Compressed_DS()
		{
			icigar = new char[MAX_ICIGAR_LENGTH];
			num_reads = 0;
			position = 0;
		}
		~Compressed_DS()
		{

		}
};

class Pass1_Compressed_DS
{
	public:
		char *col1;
		char *col2;
		char *col3;

		Pass1_Compressed_DS()
		{
			col1 = new char[MAX_SEQ_LEN * 2];
			col2 = new char[MAX_ICIGAR_LENGTH_PASS1_COL2];
			col3 = new char[MAX_SEQ_LEN * 2];
		}
		~Pass1_Compressed_DS()
		{

		}
		int isCommaInLine()
		{
			/*
			 * Return 1 if there is a comma in the second column of the record
			 */
			int i;

			for (i = 0; col2[i] != '\0'; i++)
				if (col2[i] == ',') return 1;

			return 0;
		}
};

class Pass2_Compressed_DS
{
	public:
		char *col1;
		char *col2;
		Pass2_Compressed_DS()
		{
			col1 = new char[MAX_SEQ_LEN * 2];
			col2 = new char[MAX_ICIGAR_LENGTH_PASS1_COL2];
		}
		~Pass2_Compressed_DS()
		{

		}
};

class Sam_Alignments
{
		/*
		 * Class to store information related to each alignment in sam file
		 */
	public:
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
		Sam_Tags tags;

		/*
		 * Other fields created by abridge
		 */
		Soft_Clippings soft_clippings; // Stores the left and the right soft clippings
		Cigar_Items cigar_items;
		int read_seq_len; // basically strlen(seq). Field is kept to prevent extra calls to strlen

		char *temp; //For temporary operations
		char *cigar_extended;
		char *md_extended;
		char *icigar; // Stores the integrated representation comprising of all relevant information about the alignment
		//char **splices;
		char *soft_clips_removed_seq;
		int soft_clips_removed_seq_len;
		int md_extended_length;
		char **splices;

		Sam_Alignments();

		void populateSamAlignmentInstance(char**, int, char**);

		int isSequenceSoftClipped();

		friend std::ostream& operator<<(std::ostream &output, const Sam_Alignments &S)
		{
			output << "\n";
			output << "^^^^^^^^^^^^^^^^^SAM ALIGNMENT^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^";
			output << "\n";
			output << "Read name: " << S.read_name;
			output << "\n";
			output << "SamFlag: " << S.samflag;
			output << "\n";
			output << "Reference name: " << S.reference_name;
			output << "\n";
			output << "Start position: " << S.start_position;
			output << "\n";
			output << "Mapping quality score: " << S.mapping_quality_score;
			output << "\n";
			output << "Cigar: " << S.cigar;
			output << "\n";
			output << "Reference name next mate: " << S.reference_name_next_mate;
			output << "\n";
			output << "Start position next: " << S.start_position_next;
			output << "\n";
			output << "Template length: " << S.template_length;
			output << "\n";
			output << "Sequence: " << S.seq << " " << S.read_seq_len;
			output << "\n";
			output << "Quality scores: " << S.qual;
			output << "\n";
			for (int i = 0; i < S.tags.number_of_tags; i++)
			{
				output << "Tag item %s Tag value %s" << S.tags.name[i] << S.tags.val[i];
				output << "\n";
			}
			output << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
			return output;
		}

		void splitCigar();

		int isAlignmentPerfect(int, int, int);

		void replaceNucleotidesInSeqWithInsertSymbols();

		void convertRegularCIGARToStringRepresentation();

		void designIntegratedCIGAR(int);

		void generateIntegratedCigar(short int, short int, short int);
};

/*
 * Defining member functions of Sam_Alignment class
 */

Sam_Alignments::Sam_Alignments()
{
	int i;
	/*
	 * Mandatory fields
	 */
	read_name = new char[MAX_SEQ_LEN];
	samflag = 0;
	reference_name = new char[MAX_SEQ_LEN];
	cigar = new char[MAX_SEQ_LEN];
	reference_name = new char[MAX_SEQ_LEN];
	start_position = 0;
	mapping_quality_score = 0;
	cigar = new char[MAX_CIGAR_LENGTH];
	reference_name_next_mate = new char[MAX_SEQ_LEN];
	start_position_next = 0;
	template_length = 0;
	seq = new char[MAX_SEQ_LEN];
	qual = new char[MAX_SEQ_LEN];

	/*
	 * tags will be created by automatic call to constructor of Sam_Tags class
	 */

	/*
	 * soft_clippings will be created by automatic call to constructor of Soft_Clippings class
	 */

	/*
	 * cigar_items will be created by automatic call to constructor of Cigar_Items class
	 */
	read_seq_len = 0;
	temp = new char[MAX_GENERAL_LEN];
	cigar_extended = new char[MAX_CIGAR_LENGTH];
	md_extended = new char[MAX_CIGAR_LENGTH];
	icigar = new char[MAX_CIGAR_LENGTH];
	soft_clips_removed_seq = new char[MAX_SEQ_LEN];
	soft_clips_removed_seq_len = 0;
	md_extended_length = 0;
	splices = new char*[ROWS];
	for (i = 0; i < ROWS; i++)
		splices[i] = new char[COLS];
}

void Sam_Alignments::populateSamAlignmentInstance(char **src, int number_of_fields, char **split_tags)
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	char *temp;
	int i;

	/********************************************************************/

	// Assign the first 11 mandatory sam fields
	strcpy(read_name, src[0]);
	samflag = strtol(src[1], &temp, 10);
	strcpy(reference_name, src[2]);
	start_position = strtol(src[3], &temp, 10);
	mapping_quality_score = strtol(src[4], &temp, 10);
	strcpy(cigar, src[5]);
	strcpy(reference_name_next_mate, src[6]);
	start_position_next = strtol(src[7], &temp, 10);
	template_length = strtol(src[8], &temp, 10);
	strcpy(seq, src[9]);
	strcpy(qual, src[10]);

	// Assign SAM tags
	for (i = 11; i < number_of_fields; i++)
	{
		//printf("%s\n",src[i]);
		//fflush ( stdout );
		splitByDelimiter(src[i], ':', split_tags);
		if (i == number_of_fields - 1) split_tags[2][strlen(split_tags[2])] = '\0';
		strcpy(tags.name[i - 11], split_tags[0]);
		strcpy(tags.type[i - 11], split_tags[1]);
		strcpy(tags.val[i - 11], split_tags[2]);
		//printf("\n Tags %s Parts of the tag %s %s %s ", src[i], split_tags[0], split_tags[1], split_tags[2]);
	}
	read_seq_len = strlen(seq);
	tags.number_of_tags = number_of_fields - 11;
}

int Sam_Alignments::isSequenceSoftClipped()
{
	int i;
	for (i = 0; cigar[i] != '\0'; i++)
		if (cigar[i] == 'S') return 1;
	return 0;
}

void Sam_Alignments::splitCigar()
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/

	int i, j = 0;
	int current_length = 0;
	char *temp; //Useless
	/********************************************************************/

	for (i = 0; cigar[i] != '\0'; i++)
	{
		if (isdigit(cigar[i]) != 0) //cigar[i] is a digit
		{
			current_length = current_length * 10 + cigar[i] - 48;
		}
		else if (isalpha(cigar[i]) != 0) //cigar[i] is a character
		{
			cigar_items.def[j] = cigar[i];
			cigar_items.len[j] = current_length;
			current_length = 0;
			j += 1;
		}
		else
		{
			//Should not come here
		}
	}
	cigar_items.number_of_cigar_items = j;
}

int Sam_Alignments::isAlignmentPerfect(int MD_tag_index, int NM_tag_index, int nM_tag_index)
{
	int i;
	if (NM_tag_index != -1 && nM_tag_index != -1)
	{
		if (tags.val[NM_tag_index][0] == '0' && tags.val[nM_tag_index][0] == '0') return 1;
		else return 0;
	}

	/*
	 * Checking for insertions
	 */
	for (i = 0; cigar[i] != '\0'; i++)
		if (cigar[i] == 'I') return 0;

	/*
	 * Checking for deletions
	 */
	for (i = 0; tags.val[MD_tag_index][i] != '\0'; i++)
		if (tags.val[MD_tag_index][i] == '^' || (tags.val[MD_tag_index][i] >= 65 && tags.val[MD_tag_index][i] <= 90) || (tags.val[MD_tag_index][i] >= 97 && tags.val[MD_tag_index][i] <= 122)) return 0;

	return 1;
}

void Sam_Alignments::replaceNucleotidesInSeqWithInsertSymbols()
{
	int i;
	int j;
	int sequence_index_to_be_replaced = 0;
	for (i = 0; i <= cigar_items.number_of_cigar_items; i++)
	{
		switch (cigar_items.def[i])
		{
			case 'M':
				sequence_index_to_be_replaced += cigar_items.len[i];
				break;
			case 'I':
				for (j = 0; j < cigar_items.len[i]; j++)
				{
					switch (seq[sequence_index_to_be_replaced])
					{
						case 'A':
							seq[sequence_index_to_be_replaced] = insert_characters[0];
							break;
						case 'T':
							seq[sequence_index_to_be_replaced] = insert_characters[1];
							break;
						case 'G':
							seq[sequence_index_to_be_replaced] = insert_characters[2];
							break;
						case 'C':
							seq[sequence_index_to_be_replaced] = insert_characters[3];
							break;
						case 'N':
							seq[sequence_index_to_be_replaced] = insert_characters[4];
							break;
					}
					sequence_index_to_be_replaced++;
				}
				sequence_index_to_be_replaced += cigar_items.len[i];
				break;
		}
	}
}

void Sam_Alignments::convertRegularCIGARToStringRepresentation()
{
	int i;
	int j;
	int cigar_extended_index = 0;
	for (i = 0; i < cigar_items.number_of_cigar_items; i++)
	{
		switch (cigar_items.def[i])
		{
			case 'M':
				for (j = 0; j < cigar_items.len[i]; j++)
					cigar_extended[cigar_extended_index++] = 'M';
				break;
			case 'I':
				for (j = 0; j < cigar_items.len[i]; j++)
					cigar_extended[cigar_extended_index++] = 'I';
				break;
		}
	}
}

void Sam_Alignments::designIntegratedCIGAR(int MD_tag_index)
{
	int i;
	int num;
	int j;
	int md_extended_index;
	int splice_number_index;
	int md_extended_length;
	int splices_index;
	char temp_convert_int_to_string[100];
	int intron_splice_index;
	int print_outputs = 0;

	md_extended_index = 0;
	i = 0;
	num = 0;
	if (print_outputs) printf("\n%s", seq);
	replaceNucleotidesInSeqWithInsertSymbols();
	if (print_outputs) printf("\n%s", seq);
	convertRegularCIGARToStringRepresentation();

	/*
	 * Expand the md string
	 */
	i = 0;
	while (tags.val[MD_tag_index][i] != '\0')
	{
		if (isdigit(tags.val[MD_tag_index][i]) != 0) // tags.val[MD_tag_index][i] is a digit
		num = num * 10 + tags.val[MD_tag_index][i] - 48;
		else if (tags.val[MD_tag_index][i] == '^')
		{
			for (j = 0; j < num; j++)
				md_extended[md_extended_index++] = '=';
			num = 0;
			i += 1;
			while (isdigit(tags.val[MD_tag_index][i]) == 0) // Iterate till you pick up all the deletions and ignore them for now
				i += 1;
			i -= 1;
		}
		else
		{
			for (j = 0; j < num; j++)
				md_extended[md_extended_index++] = '=';
			md_extended[md_extended_index++] = 'X';
			num = 0;
		}
		i += 1;
	}
	for (j = 0; j < num; j++)
		md_extended[md_extended_index++] = '=';
	md_extended[md_extended_index++] = '\0';
	md_extended_length = md_extended_index - 1;

	if (print_outputs) printf("\n%s %d", md_extended, md_extended_length);

	/*
	 * Add insert symbols to md_extended
	 */
	for (i = 0; i < soft_clips_removed_seq_len; i++)
		if (seq[i] >= 33 && seq[i] <= 37) insertCharacterInString(md_extended, &md_extended_length, seq[i], i, 1);
	/*printf("\n%s %d", md_extended, md_extended_length);*/

	/*
	 * Add mismatch characters to md_extended
	 */

	for (i = 0; i < md_extended_length; i++)
	{
		if (md_extended[i] == 'X')
		{
			switch (seq[i])
			{
				case 'A':
					md_extended[i] = mismatch_characters[0];
					break;
				case 'T':
					md_extended[i] = mismatch_characters[1];
					break;
				case 'G':
					md_extended[i] = mismatch_characters[2];
					break;
				case 'C':
					md_extended[i] = mismatch_characters[3];
					break;
				case 'N':
					md_extended[i] = mismatch_characters[4];
					break;
			}
		}
	}

	/*
	 * Encode deletions and splice junction info in md_extended
	 */
	md_extended_index = 0;
	splice_number_index = 0;
	for (i = 0; i < cigar_items.number_of_cigar_items; i++)
	{
		switch (cigar_items.def[i])
		{
			case 'M':
				md_extended_index += cigar_items.len[i];
				break;
			case 'N':
				md_extended[md_extended_index] = splice_number_index++ + 48;
				break;
			case 'I':
				md_extended_index += cigar_items.len[i];
				break;
			case 'D':
				insertCharacterInString(md_extended, &md_extended_length, '-', md_extended_index, cigar_items.len[i]);
				break;
		}
	}
	if (print_outputs) printf("\n%s %d", md_extended, md_extended_length);

	/*
	 * Generate the informative CIGAR
	 */

	/*
	 * Combine the splice items into an array
	 */
	splices_index = 0;
	for (i = 0; i < cigar_items.number_of_cigar_items; i++)
	{
		if (cigar_items.def[i] == 'N')
		{
			sprintf(temp_convert_int_to_string, "%d", cigar_items.len[i]);
			strcpy(splices[splices_index], temp_convert_int_to_string);
			strcat(splices[splices_index], "N");
			splices_index++;
		}
	}
	num = 0;
	icigar[0] = '\0';
	intron_splice_index = 0;
	for (i = 0; i < md_extended_length; i++)
	{
		if (md_extended[i] == '=') num++;
		else
		{
			sprintf(temp_convert_int_to_string, "%d", num);
			if (num > 0)
			{
				strcat(icigar, temp_convert_int_to_string);
				strcat(icigar, "M");
			}
			num = 0;
			if ((md_extended[i] >= 33 && md_extended[i] <= 37) || (md_extended[i] >= 38 && md_extended[i] <= 42)) // Either a mismatch or an insertion
			{
				temp_convert_int_to_string[1] = '\0';
				temp_convert_int_to_string[0] = md_extended[i];
				strcat(icigar, temp_convert_int_to_string);
			}
			else if (md_extended[i] >= '0' && md_extended[i] <= '9')
			{
				num = 0;
				intron_splice_index = 0;
				while (md_extended[i] >= '0' && md_extended[i] <= '9')
				{
					intron_splice_index = intron_splice_index * 10 + md_extended[i] - 48;
					i++;
					num++;
				}
				strcat(icigar, splices[intron_splice_index]);
				i--;
			}
			else if (md_extended[i] == '-') // Deletion
			{
				while (md_extended[i] == '-')
				{
					num++;
					i++;
				}
				i--;
				sprintf(temp_convert_int_to_string, "%d", num);
				strcat(icigar, temp_convert_int_to_string);
				strcat(icigar, "D");
				num = 0;
			}
		}
	}
	sprintf(temp_convert_int_to_string, "%d", num);
	strcat(icigar, temp_convert_int_to_string);
	strcat(icigar, "M");

	/*times[1] = clock();*/
	if (print_outputs) printf("\n%s", icigar);
}

void Sam_Alignments::generateIntegratedCigar(short int flag_ignore_soft_clippings, short int flag_ignore_mismatches, short int flag_ignore_unmapped_sequences)
{
	/*
	 * Creates the integrated cigar
	 */
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	int left_soft_clip_point = 0;
	int right_soft_clip_point = 0;
	int i;
	int number_of_mismatches_in_each_pair = -1;
	int MD_tag_index;
	int nM_tag_index;
	int XS_tag_index;
	int NH_tag_index;
	int NM_tag_index;
	int print_outputs = 0;
	int perfect_alignment_indicator = 0;
	int spliced_alignment_indicator = 0;
	char str[1000];
	char temp_str[5];
	char M_replacement_character;
	/********************************************************************/
	splitCigar();

	for (i = 0; i < cigar_items.number_of_cigar_items; i++)
	{
		if (cigar_items.def[i] == 'N')
		{
			spliced_alignment_indicator = 1;
			break;
		}
	}
	/*
	 * Process each alignment to extract soft clipped portion of reads
	 */
	if (isSequenceSoftClipped() == 1)
	{
		if (cigar_items.def[0] == 'S')
		{
			left_soft_clip_point = cigar_items.len[0];
			extractSubString(seq, soft_clippings.left, 0, left_soft_clip_point - 1);
			if (print_outputs > 1) printf("\nLSC: %d %s", cigar_items.len[0], soft_clippings.left);
			soft_clips_removed_seq_len = read_seq_len - cigar_items.len[0];
			for (i = 0; i <= left_soft_clip_point - 1; i++)
				soft_clippings.left[i] = (soft_clippings.left[i] >= 65 && soft_clippings.left[i] <= 90) ? soft_clippings.left[i] + 32 : soft_clippings.left[i];
		}
		if (cigar_items.def[cigar_items.number_of_cigar_items - 1] == 'S')
		{
			right_soft_clip_point = cigar_items.len[cigar_items.number_of_cigar_items - 1];
			extractSubString(seq, soft_clippings.right, read_seq_len - right_soft_clip_point, read_seq_len);
			if (print_outputs > 1) printf("\nRSC: %d %s", cigar_items.len[cigar_items.number_of_cigar_items - 1], soft_clippings.right);
			soft_clips_removed_seq_len = read_seq_len - cigar_items.len[cigar_items.number_of_cigar_items - 1];

			for (i = 0; i <= right_soft_clip_point; i++)
				soft_clippings.right[i] = (soft_clippings.right[i] >= 65 && soft_clippings.right[i] <= 90) ? soft_clippings.right[i] + 32 : soft_clippings.right[i];
		}

		/*
		 * Remove the soft-clips and construct new sequence
		 */
		if (left_soft_clip_point != 0 || right_soft_clip_point != 0)
		{
			for (i = left_soft_clip_point; i <= read_seq_len - right_soft_clip_point; i++)
				soft_clips_removed_seq[i] = seq[i];
			if (print_outputs > 1) printf("\nSEQ: %s newSEQ: %s", seq, soft_clips_removed_seq);
		}
		else
		{
			soft_clips_removed_seq = seq;
			if (print_outputs > 1) printf("\nSEQ: %s newSEQ: %s", seq, soft_clips_removed_seq);
		}
	}
	else
	{
		soft_clips_removed_seq = seq;
		if (print_outputs > 1) printf("\nSEQ: %s newSEQ: %s", seq, soft_clips_removed_seq);
	}

	/*
	 * Find the different tags
	 */
	XS_tag_index = -1;
	NH_tag_index = -1;
	NM_tag_index = -1;
	nM_tag_index = -1;
	MD_tag_index = -1;
	for (i = 0; i < tags.number_of_tags; i++)
	{
		if (strcmp(tags.name[i], "MD") == 0) MD_tag_index = i;
		if (strcmp(tags.name[i], "XS") == 0) XS_tag_index = i;
		if (strcmp(tags.name[i], "NH") == 0) NH_tag_index = i;
		if (strcmp(tags.name[i], "NM") == 0) NM_tag_index = i;
		if (strcmp(tags.name[i], "nM") == 0) nM_tag_index = i;
	}

	perfect_alignment_indicator = isAlignmentPerfect(MD_tag_index, NM_tag_index, nM_tag_index);

	if (perfect_alignment_indicator == 0) designIntegratedCIGAR(MD_tag_index);
	else
	{
		if (cigar_items.def[0] == 'S' && cigar_items.def[cigar_items.number_of_cigar_items - 1] == 'S')
		{
			strcpy(icigar, "");
			for (i = 1; i < cigar_items.number_of_cigar_items - 1; i++)
			{
				sprintf(str, "%d", cigar_items.len[i]);
				strcat(icigar, str);
				temp_str[0] = cigar_items.def[i];
				temp_str[1] = '\0';
				strcat(icigar, temp_str);
			}
		}
		else if (cigar_items.def[0] != 'S' && cigar_items.def[cigar_items.number_of_cigar_items - 1] == 'S')
		{
			strcpy(icigar, "");
			for (i = 0; i < cigar_items.number_of_cigar_items - 1; i++)
			{
				sprintf(str, "%d", cigar_items.len[i]);
				strcat(icigar, str);
				temp_str[0] = cigar_items.def[i];
				temp_str[1] = '\0';
				strcat(icigar, temp_str);
			}
		}
		else if (cigar_items.def[0] == 'S' && cigar_items.def[cigar_items.number_of_cigar_items - 1] != 'S')
		{
			strcpy(icigar, "");
			for (i = 1; i < cigar_items.number_of_cigar_items; i++)
			{
				sprintf(str, "%d", cigar_items.len[i]);
				strcat(icigar, str);
				temp_str[0] = cigar_items.def[i];
				temp_str[1] = '\0';
				strcat(icigar, temp_str);
			}
		}
		else strcpy(icigar, cigar);
	}
	if (print_outputs > 1) printf("\nPrevious iCIGAR: %s", icigar);

	if (left_soft_clip_point != 0 && right_soft_clip_point == 0)
	{
		strcpy(temp, soft_clippings.left);
		strcat(temp, icigar);
	}
	else if (left_soft_clip_point == 0 && right_soft_clip_point != 0)
	{
		strcpy(temp, icigar);
		strcat(temp, soft_clippings.right);
	}
	else if (left_soft_clip_point != 0 && right_soft_clip_point != 0)
	{
		strcpy(temp, soft_clippings.left);
		strcat(temp, icigar);
		strcat(temp, soft_clippings.right);
	}
	else
	{
		strcpy(temp, icigar);
	}

	strcpy(icigar, temp);

	/*
	 * Add NH tag
	 */
	/*if (XS_tag_index != -1) strcat(icigar, tags[XS_tag_index].val);*/
	if (NH_tag_index != -1) strcat(icigar, tags.val[NH_tag_index]);

	/*
	 * Change the iCIGAR representation to reflect the samformatflag and XS tag
	 */
	/*
	 switch (samflag)
	 {
	 case 0:
	 M_replacement_character = 'B';
	 break;
	 case 16:
	 M_replacement_character = 'C';
	 break;
	 case 256:
	 M_replacement_character = 'E';
	 break;
	 case 272:
	 M_replacement_character = 'F';
	 break;
	 }
	 for (i = 0; i < strlen(icigar); i++)
	 if (icigar[i] == 'M') icigar[i] = M_replacement_character;
	 */
	if (spliced_alignment_indicator == 0)
	{
		switch (samflag)
		{
			case 0:
				M_replacement_character = 'B';
				break;
			case 16:
				M_replacement_character = 'C';
				break;
			case 256:
				M_replacement_character = 'E';
				break;
			case 272:
				M_replacement_character = 'F';
				break;
				break;
		}
		for (i = 0; i < strlen(icigar); i++)
			if (icigar[i] == 'M') icigar[i] = M_replacement_character;
	}
	else
	{
		if (XS_tag_index == -1)
		{
			switch (samflag)
			{
				case 0:
					M_replacement_character = 'B';
					break;
				case 16:
					M_replacement_character = 'C';
					break;
				case 256:
					M_replacement_character = 'E';
					break;
				case 272:
					M_replacement_character = 'F';
					break;
			}
		}
		else
		{
			if (strcmp(tags.val[XS_tag_index], "+") == 0)
			{
				switch (samflag)
				{
					case 0:
						M_replacement_character = 'H';
						break;
					case 16:
						M_replacement_character = 'J';
						break;
					case 256:
						M_replacement_character = 'K';
						break;
					case 272:
						M_replacement_character = 'L';
						break;
				}
			}
			else
			{
				switch (samflag)
				{
					case 0:
						M_replacement_character = 'O';
						break;
					case 16:
						M_replacement_character = 'P';
						break;
					case 256:
						M_replacement_character = 'Q';
						break;
					case 272:
						M_replacement_character = 'R';
						break;
				}
			}

		}
		for (i = 0; i < strlen(icigar); i++)
			if (icigar[i] == 'M') icigar[i] = M_replacement_character;
	}
	/*
	 switch (samflag)
	 {
	 case 0:
	 strcat(icigar, samformat_characters_SE[0]);
	 break;
	 case 16:
	 strcat(icigar, samformat_characters_SE[1]);
	 break;
	 case 256:
	 strcat(icigar, samformat_characters_SE[2]);
	 break;
	 case 272:
	 strcat(icigar, samformat_characters_SE[3]);
	 break;
	 }*/

	/*
	 if (samflag == 256)
	 {
	 printSamAlignmentInstance(curr_alignment);
	 printf("\nIntegrated CIGAR: %s CIGAR: %s MD: %s perfect_alignment_indicator: %d SEQ: %s READ_NAME: %s", icigar, cigar, tags[MD_tag_index].val, perfect_alignment_indicator, seq, read_name);
	 }*/

	if (print_outputs > 0)
		printf("\nIntegrated CIGAR: %s CIGAR: %s MD: %s perfect_alignment_indicator: %d SEQ: %s READ_NAME: %s", icigar, cigar, tags.val[MD_tag_index], perfect_alignment_indicator, seq, read_name);
	if (print_outputs > 1) printf("\n\n\n");
}

/************************************************************************************************************************************************************************************************************************/

/*
 * Definitions of non-member functions
 */

int splitByDelimiter(char *line, char delimiter, char **new_string)
{
	/*********************************************************************
	 * Splits a string on a requested character and
	 * returns the number of splits
	 *********************************************************************/

	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	int i = 0, j = 0, ctr = 0;
	/********************************************************************/

	for (i = 0; line[i] != '\0'; i++)
	{
		// if space or NULL found, assign NULL into new_string[ctr]
		if (line[i] == delimiter)
		{
			new_string[ctr][j] = '\0';
			ctr++; //for next word
			j = 0; //for next word, init index to 0
		}
		else if (line[i] == '\n') continue;
		else
		{
			new_string[ctr][j] = line[i];
			j++;
		}
	}
	new_string[ctr][j] = '\0';
	return ctr + 1;
}

void extractSubString(char *str, char *substr, int start_index, int end_index)
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	int length = strlen(str);
	int i;
	int j = 0;

	/********************************************************************/

	if (length < (end_index - start_index + 1))
	{
		return;
	}
	for (i = start_index; i <= end_index; i++)
		substr[j++] = str[i];
	substr[j] = '\0';

}

int isCharacterInString(char *str, char character)
{
	int i;
	for (i = 0; str[i] != '\0'; i++)
		if (character == str[i]) return 1;
	return 0;
}

int isNumber(char *str)
{
	int i;
	for (i = 0; str[i] != '\0'; i++)
		if (str[i] < 48 || str[i] > 57) return 0;
	return 1;
}

void insertCharacterInString(char *str, int *str_length, char ins, int loc, int number_of_insertions_to_be_made)
{
	/*
	 * Inserts a character. Assumes that a large
	 */
	int i = (*str_length) + number_of_insertions_to_be_made;
	for (; i > loc; i--)
		str[i] = str[i - number_of_insertions_to_be_made];
	(*str_length) += number_of_insertions_to_be_made;
	while (number_of_insertions_to_be_made--)
		str[i++] = ins;
}

#endif /* ABRIDGE_CLASS_DEFINITIONS_HPP_ */
