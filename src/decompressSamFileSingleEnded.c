# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void readAbridgeIndex(struct Abridge_Index *abridge_index, char *abridge_index_filename, char **split_line)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;

	size_t len = 0;
	ssize_t line_len;

	char *temp; //Useless
	char *line = NULL; // for reading each line

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen(abridge_index_filename, "r");
	if (fhr == NULL)
	{
		printf("Error! File not found");
		exit(1);
	}
	/********************************************************************/
	while ((line_len = getline(&line, &len, fhr)) != -1)
	{
		splitByDelimiter(line, '\t', split_line);
		//printf("\n%d %s", abridge_index->number_of_items, split_line[0]);
		//fflush(stdout);
		strcpy(abridge_index->chromosome[abridge_index->number_of_items], split_line[0]);
		abridge_index->start[abridge_index->number_of_items] = strtol(split_line[1], &temp, 10);
		abridge_index->end[abridge_index->number_of_items] = strtol(split_line[2], &temp, 10);
		abridge_index->start_byte[abridge_index->number_of_items] = strtol(split_line[3], &temp, 10);
		abridge_index->end_byte[abridge_index->number_of_items] = strtol(split_line[4], &temp, 10);
		abridge_index->number_of_items++;
		//printf("\n %lld", MAX_POOL_SIZE - abridge_index->number_of_items);
	}
	fclose(fhr);
}

int maxClusterSize(struct Abridge_Index *abridge_index)
{
	long long int i;
	long long int max_cluster_size = -1;
	for (i = 0; i < abridge_index->number_of_items; i++)
		if (max_cluster_size < (abridge_index->end[i] - abridge_index->start[i])) max_cluster_size = (abridge_index->end[i] - abridge_index->start[i]);
	return max_cluster_size;

}

void loadSequencesFromFile(char **sequences, char *filename)
{
	FILE *fhr;
	long long int num_sequence = 0;
	size_t len = 0;
	ssize_t line_len;
	char *line = NULL; // for reading each line

	fhr = fopen(filename, "rb");
	if (fhr == NULL)
	{
		printf("Error! File not found");
		exit(1);
	}
	while ((line_len = getline(&line, &len, fhr)) != -1)
		if (line[0] != '>') strcpy(sequences[num_sequence++], line);

}

long long int reverseInteger(long long int val)
{
	long long int temp = 0;
	while (val)
	{
		temp = temp * 10 + (val % 10);
		val /= 10;
	}
	return temp;
}

long long int extractNHfromicigar(char *icigar, int icigar_length)
{
	char nh_val_string[15];
	char nh_val_string_rev[15];
	long long int nh_val = 0;
	int i;
	int j;
	int k;
	char *t;

	j = 0;
	for (i = icigar_length - 1; i >= 0; i--)
	{
		if (isalpha(icigar[i]) != 0) break;
		//nh_val = nh_val * 10 + icigar[i] - 48;
		nh_val_string[j++] = icigar[i];
	}
	icigar[i + 1] = '\0';
	nh_val_string[j] = '\0';

	k = j - 1;
	for (i = 0; i < j; i++)
		nh_val_string_rev[i] = nh_val_string[k--];
	nh_val_string_rev[j] = '\0';
	nh_val = strtol(nh_val_string_rev, &t, 10);
	return nh_val;
}

void replaceCharacterInString(char *str, char ch_to_be_replaced, char replace_with, int str_length)
{
	int i;
	for (i = 0; i < str_length; i++)
		if (str[i] == ch_to_be_replaced) str[i] = replace_with;
}

int findSamFormatFlag(char *icigar, int icigar_length)
{
	int i;
	int samformatflag;

	for (i = 0; i < icigar_length; i++)
	{
		switch (icigar[i])
		{
			case 'B':
				//replaceCharacterInString(icigar, 'B', 'M', icigar_length);
				samformatflag = 0;
				break;
			case 'E':
				//replaceCharacterInString(icigar, 'E', 'M', icigar_length);
				samformatflag = 16;
				break;
			case 'F':
				//replaceCharacterInString(icigar, 'F', 'M', icigar_length);
				samformatflag = 256;
				break;
			case 'H':
				//replaceCharacterInString(icigar, 'H', 'M', icigar_length);
				samformatflag = 272;
				break;
			case 'J':
				//replaceCharacterInString(icigar, 'J', 'M', icigar_length);
				samformatflag = 0;
				break;
			case 'K':
				//replaceCharacterInString(icigar, 'K', 'M', icigar_length);
				samformatflag = 16;
				break;
			case 'L':
				//replaceCharacterInString(icigar, 'L', 'M', icigar_length);
				samformatflag = 256;
				break;
			case 'O':
				//replaceCharacterInString(icigar, 'O', 'M', icigar_length);
				samformatflag = 272;
				break;
			case 'P':
				//replaceCharacterInString(icigar, 'P', 'M', icigar_length);
				samformatflag = 0;
				break;
			case 'Q':
				//replaceCharacterInString(icigar, 'Q', 'M', icigar_length);
				samformatflag = 16;
				break;
			case 'R':
				//replaceCharacterInString(icigar, 'R', 'M', icigar_length);
				samformatflag = 256;
				break;
			case 'U':
				//replaceCharacterInString(icigar, 'U', 'M', icigar_length);
				samformatflag = 272;
				break;
		}
	}
	return samformatflag;
}

void convertIcigarToCigarandMD(struct Sam_Alignment *sam_alignment_instance)
{
	/*
	 * Coverts the iCIGAR into CIGAR and MD
	 * Returns the samformatflag
	 */

	int samformatflag = -1;
	int NH_value;
	int icigar_length;
	int i;
	int j;
	int processing_left_soft_clip = 1;
	int left_soft_clip_index = 0;
	int right_soft_clip_index = 0;
	int cigar_index = 0;
	int num = 0;
	int cigar_temp_index = 0;
	int cigar_items_instance_index;
	int flag = 0;

	struct Cigar_Items cigar_items_instance[MAX_SEQ_LEN];

	char temp[100];
	char cigar_temp[MAX_SEQ_LEN];
	//printf("\niCIGAR at the beginning %s ", sam_alignment_instance->icigar);
	icigar_length = strlen(sam_alignment_instance->icigar);
	NH_value = extractNHfromicigar(sam_alignment_instance->icigar, icigar_length);

	samformatflag = findSamFormatFlag(sam_alignment_instance->icigar, icigar_length);
	/*
	 * Construct the Cigar string, MD string, Soft Clips, etc.
	 */
	splitCigar(sam_alignment_instance->icigar, &cigar_items_instance_index, cigar_items_instance);
	//printf("\nICIGAR %s Number of items: %d", sam_alignment_instance->icigar, cigar_items_instance_index);

	/*
	 for (i = 0; i < cigar_items_instance_index; i++)
	 printf("\n%d %c ", cigar_items_instance[i].len, cigar_items_instance[i].def);
	 fflush(stdout);
	 */
	sam_alignment_instance->cigar[0] = '\0';
	for (i = 0; i < cigar_items_instance_index; i++)
	{
		if (processing_left_soft_clip == 1 && isCharacterInString("atgcn", cigar_items_instance[i].def))
		{
			sam_alignment_instance->soft_clippings.left[left_soft_clip_index] = cigar_items_instance[i].def - 32;
			left_soft_clip_index++;
		}
		else if (processing_left_soft_clip == 0 && isCharacterInString("atgcn", cigar_items_instance[i].def))
		{
			sam_alignment_instance->soft_clippings.left[right_soft_clip_index] = cigar_items_instance[i].def - 32;
			right_soft_clip_index++;
		}
		else // Capital letters & special characters - cannot be soft clips
		{
			flag = 0;
			processing_left_soft_clip = 0;
			if (left_soft_clip_index > 0) //There were some left soft clips
			{
				sprintf(temp, "%d", left_soft_clip_index);
				strcat(sam_alignment_instance->cigar, temp);
				strcat(sam_alignment_instance->cigar, "S");
				num = 0;
				left_soft_clip_index = 0;
			}
			if (cigar_items_instance[i].def == 'N')
			{
				sprintf(temp, "%d", cigar_items_instance[i].len);
				strcat(sam_alignment_instance->cigar, temp);
				strcat(sam_alignment_instance->cigar, "N");
				num = 0;
			}
			else if (cigar_items_instance[i].def == 'D')
			{
				sprintf(temp, "%d", cigar_items_instance[i].len);
				strcat(sam_alignment_instance->cigar, temp);
				strcat(sam_alignment_instance->cigar, "D");
				num = 0;
			}
			num = 0;
			// Merging consecutive matches, mismatches, insertions
			while (i < cigar_items_instance_index && (isCharacterInString("BEFHJKLOPQRU", cigar_items_instance[i].def) || isCharacterInString(mismatch_characters, cigar_items_instance[i].def) || isCharacterInString(insert_characters, cigar_items_instance[i].def)))
			{
				if (isCharacterInString("BEFHJKLOPQRU", cigar_items_instance[i].def)) num += cigar_items_instance[i].len;
				else if (isCharacterInString(mismatch_characters, cigar_items_instance[i].def) || isCharacterInString(insert_characters, cigar_items_instance[i].def)) num++;
				i++;
				flag = 1;
			}
			if (flag)
			{
				sprintf(temp, "%d", num);
				strcat(sam_alignment_instance->cigar, temp);
				strcat(sam_alignment_instance->cigar, "M");
				num = 0;
			}
			if (flag && i < cigar_items_instance_index) i--;
		}
	}
	if (right_soft_clip_index > 0) //There were some left soft clips
	{
		sprintf(temp, "%d", right_soft_clip_index);
		strcat(sam_alignment_instance->cigar, temp);
		strcat(sam_alignment_instance->cigar, "S");
		num = 0;
	}
	/*
	 printf("\nICIGAR %s CIGAR %s NH:%d samformat_flag: %d", sam_alignment_instance->icigar, sam_alignment_instance->cigar, NH_value, samformatflag);
	 printf("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
	 printf("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
	 */
}

void convertToAlignment(struct Sam_Alignment **sam_alignment_pool, int sam_alignment_pool_index, struct Whole_Genome_Sequence *whole_genome, char **split_on_newline, struct Sam_Alignment *sam_alignment, int cluster_index, struct Abridge_Index *abridge_index, int number_of_entries_in_cluster, char **split_on_tab, char **split_on_dash, char **split_on_comma, FILE *fhw)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	long long int start_position;

	int number_of_columns;
	int number_of_distinct_cigars_in_a_line;
	int number_of_reads;
	int samformatflag;
	int i;
	int j;

	char *temp; //Useless
	char *distinct_icigars_in_a_line;
	char *icigar;

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	distinct_icigars_in_a_line = (char*) malloc(sizeof(char) * MAX_GENERAL_LEN);

	/********************************************************************/
	for (i = 0; i < number_of_entries_in_cluster; i++)
	{
		//printf("\nProcessing this line %s", split_on_newline[i]);
		number_of_columns = splitByDelimiter(split_on_newline[i], '\t', split_on_tab);
		//printf("\nNumber of columns %d", number_of_columns);

		if (number_of_columns == 1)
		{
			strcpy(distinct_icigars_in_a_line, split_on_tab[0]);
			start_position++;
		}
		else if (number_of_columns == 2)
		{
			if (i == 0) start_position = abridge_index->start[cluster_index];
			else start_position += strtol(split_on_tab[0], &temp, 10);
			strcpy(distinct_icigars_in_a_line, split_on_tab[1]);
		}
		else
		{
			// Should never enter here
			printf("\nTrouble");
		}
		number_of_distinct_cigars_in_a_line = splitByDelimiter(distinct_icigars_in_a_line, ',', split_on_comma);
		sam_alignment_pool_index = 0;
		for (j = 0; j < number_of_distinct_cigars_in_a_line; j++)
		{
			splitByDelimiter(split_on_comma[j], '-', split_on_dash);
			number_of_reads = strtol(split_on_dash[1], &temp, 10);
			if (split_on_dash[0][1] == '-' && isalpha(split_on_dash[0][0]) != 0)
			{
				// Use the same cigar
			}
			else
			{
				//printf("\nSplit_on_dash %s %d sam_alignment_pool_index %d ", split_on_dash[0], strlen(split_on_dash[0]), sam_alignment_pool_index);
				strcpy(sam_alignment_pool[sam_alignment_pool_index]->icigar, split_on_dash[0]);
				//printf("\nsam_alignment_pool[i]->icigar %s %d Original_icigar %s %d", sam_alignment_pool[sam_alignment_pool_index]->icigar, strlen(sam_alignment_pool[sam_alignment_pool_index]->icigar), split_on_dash[0], strlen(split_on_dash[0]));
				convertIcigarToCigarandMD(sam_alignment_pool[sam_alignment_pool_index]);
				sam_alignment_pool_index++;
			}
		}
	}
}

void readInTheEntireGenome(char *genome_filename, struct Whole_Genome_Sequence *whole_genome)
{
	FILE *fhr;
	char *buffer = NULL;
	size_t len = 0;
	ssize_t line_len;
	int i;

//buffer = (char*) malloc(sizeof(char) * pow(2, 32));
	fhr = fopen(genome_filename, "rb");
	if (fhr == NULL)
	{
		printf("Error! File not found");
		exit(1);
	}

	whole_genome->nucleotides = (char**) malloc(sizeof(char*) * MAX_REFERENCE_SEQUENCES);
	whole_genome->reference_sequence_name = (char**) malloc(sizeof(char*) * MAX_REFERENCE_SEQUENCES);
	whole_genome->reference_sequence_length = (unsigned long long int*) malloc(sizeof(unsigned long long int) * MAX_REFERENCE_SEQUENCES);
	whole_genome->number_of_reference_sequences = 0;

	while ((line_len = getline(&buffer, &len, fhr)) != -1)
	{
		if (buffer[0] == '>')
		{
			whole_genome->reference_sequence_name[whole_genome->number_of_reference_sequences] = (char*) malloc(sizeof(char) * (len + 1));
			for (i = 0; buffer[i] != '\0'; i++)
				if (buffer[i] == 32)
				{
					buffer[i] = '\0';
					break;
				}
			for (i = 1; buffer[i] != '\0'; i++)
				buffer[i - 1] = buffer[i];
			buffer[i - 1] = buffer[i];
			strcpy(whole_genome->reference_sequence_name[whole_genome->number_of_reference_sequences], buffer);
		}
		else
		{
			whole_genome->nucleotides[whole_genome->number_of_reference_sequences] = (char*) malloc(sizeof(char) * (len + 1));
			strcpy(whole_genome->nucleotides[whole_genome->number_of_reference_sequences], buffer);
			whole_genome->reference_sequence_length[whole_genome->number_of_reference_sequences] = strlen(buffer);
			whole_genome->number_of_reference_sequences++;
		}
		free(buffer);
		buffer = NULL;
	}
	/*
	 for (i = 0; i < whole_genome->number_of_reference_sequences; i++)
	 printf("\n %d %s %lld %lld", i + 1, whole_genome->reference_sequence_name[i], strlen(whole_genome->nucleotides[i]), whole_genome->reference_sequence_length[i]);
	 */
	fclose(fhr);
}

void decompressFile(char *abridge_index_filename, char *genome_filename, char *output_sam_filename, char *pass2_filename, char *genome_prefix)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;

	struct Abridge_Index *abridge_index;
	struct Sam_Alignment *sam_alignment;

	int i;
	int sam_alignment_pool_index;
	int fread_ret_val;
	int fseek_ret_val;

	long long max_cluster_size;
	int number_of_entries_in_cluster;
	int number_of_elements_after_split_on_delimiter;

	char **split_on_newline;
	char **split_on_tab;
	char **split_on_dash;
	char **split_on_comma;
	char **split_on_delimiter;
	char *buffer;
	char **sequence_portions_from_reference;
	char *fasta_file_with_expressed_portions;
	char *cigar;
	char *md;
	char *output_prefix_without_path;

	struct Sam_Alignment **sam_alignment_pool;
	struct Whole_Genome_Sequence *whole_genome;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen(pass2_filename, "rb");
	if (fhr == NULL)
	{
		printf("Error! File not found");
		exit(1);
	}
	fhw = fopen(output_sam_filename, "w");
	if (fhw == NULL)
	{
		printf("Error! File cannot be opened for writing");
		exit(1);
	}

	split_on_newline = (char**) malloc(sizeof(char*) * ROWS * 10);
	for (i = 0; i < ROWS * 10; i++)
		split_on_newline[i] = (char*) malloc(sizeof(char) * COLS * 10);

	split_on_tab = (char**) malloc(sizeof(char*) * ROWS);
	for (i = 0; i < ROWS; i++)
		split_on_tab[i] = (char*) malloc(sizeof(char) * COLS);

	split_on_dash = (char**) malloc(sizeof(char*) * ROWS);
	for (i = 0; i < ROWS; i++)
		split_on_dash[i] = (char*) malloc(sizeof(char) * COLS);

	split_on_comma = (char**) malloc(sizeof(char*) * ROWS * 10);
	for (i = 0; i < ROWS * 10; i++)
		split_on_comma[i] = (char*) malloc(sizeof(char) * COLS * 10);

	split_on_delimiter = (char**) malloc(sizeof(char*) * ROWS);
	for (i = 0; i < ROWS; i++)
		split_on_delimiter[i] = (char*) malloc(sizeof(char) * COLS);

	output_prefix_without_path = (char*) malloc(sizeof(char) * MAX_SEQ_LEN);
	sequence_portions_from_reference = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	fasta_file_with_expressed_portions = (char*) malloc(sizeof(char) * FILENAME_LENGTH);

	buffer = (char*) malloc(sizeof(char) * MAX_BUFFER_SIZE_FOR_READING_PASS2_FILE);
	abridge_index = allocateMemoryAbridge_Index();
	sam_alignment = allocateMemorySam_Alignment();
	whole_genome = (struct Whole_Genome_Sequence*) malloc(sizeof(struct Whole_Genome_Sequence));
	sam_alignment_pool = (struct Sam_Alignment**) malloc(sizeof(struct Sam_Alignment*) * 10000);
	for (i = 0; i < 10000; i++)
		sam_alignment_pool[i] = allocateMemorySam_Alignment();

	/********************************************************************/

	readAbridgeIndex(abridge_index, abridge_index_filename, split_on_newline);
//readInTheEntireGenome(genome_filename, whole_genome);

	printf("\nProcessing each line from index file");
	fflush(stdout);
	for (i = 0; i < abridge_index->number_of_items; i++)
	{
		if (i % 1000 == 0)
		{
			printf("\n%d record processed", i);
			fflush(stdout);
		}
		//printf("\nFile pointer at %ld", ftell(fhr));
		//printf("\nEntry in index file %s %d %d", abridge_index->chromosome[i], abridge_index->start_byte[i], abridge_index->end_byte[i]);
		fseek_ret_val = fseek(fhr, abridge_index->start_byte[i], SEEK_SET);
		//printf("\nFile pointer moved to %ld", ftell(fhr));
		//continue;
		buffer[0] = '\0';
		fread_ret_val = fread(buffer, 1, abridge_index->end_byte[i] - abridge_index->start_byte[i], fhr);
		//printf("\n fread_ret_val %d fseek_ret_val %d diff %lld", fread_ret_val, fseek_ret_val, MAX_BUFFER_SIZE_FOR_READING_PASS2_FILE - fread_ret_val);
		//fflush(stdout);
		buffer[fread_ret_val] = '\0';
		/*
		 printf("\nfread_ret_val %d", fread_ret_val);
		 printf("\nNew Record:%d\n", i);
		 printf("%s", buffer);
		 printf("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
		 printf("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
		 fflush(stdout);
		 */
		number_of_entries_in_cluster = splitByDelimiter(buffer, '\n', split_on_newline);
		number_of_entries_in_cluster--; //Last line is always empty
		convertToAlignment(sam_alignment_pool, sam_alignment_pool_index, whole_genome, split_on_newline, sam_alignment, i, abridge_index, number_of_entries_in_cluster, split_on_tab, split_on_dash, split_on_comma, fhw);
		if (i == 10) break;
	}
}

int main(int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char abridge_index_filename[FILENAME_LENGTH];
	char genome_filename[FILENAME_LENGTH];
	char pass2_filename[FILENAME_LENGTH];
	char output_sam_filename[FILENAME_LENGTH];
	char genome_prefix[FILENAME_LENGTH];
	char *temp; //Required for strtoi

	short int flag_ignore_soft_clippings;
	short int flag_ignore_mismatches;
	short int flag_ignore_sequence_information;
	short int flag_ignore_unmapped_sequences;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy(abridge_index_filename, argv[1]);
	strcpy(genome_filename, argv[2]);
	strcpy(output_sam_filename, argv[3]);
	strcpy(pass2_filename, argv[4]);
	strcpy(genome_prefix, argv[5]);
	/********************************************************************/

	decompressFile(abridge_index_filename, genome_filename, output_sam_filename, pass2_filename, genome_prefix);
	return 0;
}

