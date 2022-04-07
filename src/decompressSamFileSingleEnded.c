# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

long long int total_mapped_reads = 0;

void decompressFile(
		char *name_of_file_with_quality_scores,
		char *genome_filename,
		char *output_sam_filename,
		char *pass1_filename,
		char *unmapped_filename,
		char *default_quality_value,
		short int flag_ignore_sequence_information,
		int max_reads_in_each_line)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;
	FILE *fhr_qual;

	struct Abridge_Index *abridge_index;
	struct Sam_Alignment *sam_alignment;

	int i, j;
	int sam_alignment_pool_index;
	int fread_ret_val;
	int fseek_ret_val;
	int line_number;

	size_t len = 0;
	ssize_t line_len;

	short int flag_ignore_mismatches;
	short int flag_ignore_soft_clippings;
	short int flag_ignore_unmapped_sequences;
	short int flag_ignore_quality_scores_for_matched_bases;
	short int flag_ignore_quality_scores_for_mismatched_bases_and_soft_clips;
	short int flag_ignore_alignment_scores;
	short int number_of_columns;
	short int read_names_stored;

	unsigned long long int max_cluster_size;
	unsigned long long int line_num = 0;
	unsigned long long int read_number = 1;
	unsigned long long int from = -1;
	unsigned long long int to = -1;
	unsigned long long int number_of_newlines = 0;
	unsigned long long int number_of_commas_in_each_line = 0;
	unsigned long long int max_number_of_commas = 0;
	unsigned long long int max_length_of_newline = 0;
	unsigned long long int length_of_newline = 0;
	unsigned long long int curr_position = 0;

	int number_of_entries_in_cluster;
	int number_of_elements_after_split_on_delimiter;
	int BUFFER_SIZE = 8 * 100 * 1024 * 1024; // 100 MB
	//int ROWS_split_on_newline = ROWS * 10; //10,000
	//int COLS_split_on_newline = COLS * 1000; //1,000,000
	int ROWS_split_on_tab = 5; //10
	int COLS_split_on_tab = COLS * 10; //100,000
	int ROWS_split_on_dash = 5; //5
	int COLS_split_on_dash = MAX_SEQ_LEN * 3; //3,000
	int ROWS_split_on_comma = ROWS * 10; //10,000
	int COLS_split_on_comma = MAX_SEQ_LEN * 3; //3,000

	//char **split_on_newline;
	char **split_on_tab;
	char **read_names;
	char **split_on_dash;
	char **split_on_comma;
	char **split_on_tilde;
	char *buffer = NULL;
	char **sequence_portions_from_reference;
	char *fasta_file_with_expressed_portions;
	char *cigar;
	char *md;
	char *output_prefix_without_path;
	char *current_chromosome;
	char *convert_to_int_temp;
	char line_to_be_written_to_file[MAX_GENERAL_LEN];
	char temp[100];
	char read_prefix[10];

	struct Sam_Alignment **sam_alignment_pool;
	struct Sam_Alignment *sam_alignment_instance;
	struct Whole_Genome_Sequence *whole_genome;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen(pass1_filename, "r");
	if (fhr == NULL)
	{
		printf("Error! File %s not found", pass1_filename);
		exit(1);
	}
	fhr_qual = fopen(name_of_file_with_quality_scores, "r");
	if (fhr_qual == NULL)
	{
		printf("Error! File %s not found", name_of_file_with_quality_scores);
		exit(1);
	}
	fhw = fopen(output_sam_filename, "w");
	if (fhw == NULL)
	{
		printf(
				"Error! File %s cannot be opened for writing",
				output_sam_filename);
		exit(1);
	}

	/*
	 split_on_newline = ( char** ) malloc (sizeof(char*) * ROWS_split_on_newline);
	 for ( i = 0 ; i < ROWS_split_on_newline ; i++ )
	 split_on_newline[i] = ( char* ) malloc (sizeof(char) * COLS_split_on_newline);
	 */

	/*
	 * Find the largest line in the compressed file. Apparently, free and realloc isn't working.
	 */
	int max_line_len = 0;
	int max_commas = 0;
	while ((line_len = getline(&buffer, &len, fhr)) != -1)
	{
		if (buffer[0] != '@')
		{
			if (max_line_len < line_len)
				max_line_len = line_len;

			int num_commas = 0;
			for (int i = 0; buffer[i] != '\0'; i++)
				if (buffer[i] == ',')
					num_commas += 1;
			if (max_commas < num_commas)
				max_commas = num_commas;
		}
	}
	COLS_split_on_tab = line_len / 2;
	ROWS_split_on_comma = max_commas;
	fclose(fhr);
	fhr = fopen(pass1_filename, "r");
	if (fhr == NULL)
	{
		printf("Error! File %s not found", pass1_filename);
		exit(1);
	}
	printf("\nAm here max_line_len %d max_commas %d", max_line_len, max_commas);
	fflush(stdout);

	split_on_tab = (char**) malloc(sizeof(char*) * ROWS_split_on_tab);
	for (i = 0; i < ROWS_split_on_tab; i++)
		split_on_tab[i] = (char*) malloc(sizeof(char) * COLS_split_on_tab);

	split_on_dash = (char**) malloc(sizeof(char*) * ROWS_split_on_dash);
	for (i = 0; i < ROWS_split_on_dash; i++)
		split_on_dash[i] = (char*) malloc(sizeof(char) * COLS_split_on_dash);

	split_on_comma = (char**) malloc(sizeof(char*) * ROWS_split_on_comma);
	for (i = 0; i < ROWS_split_on_comma; i++)
		split_on_comma[i] = (char*) malloc(sizeof(char) * COLS_split_on_comma);

	split_on_tilde = (char**) malloc(sizeof(char*) * ROWS_split_on_comma);
	for (i = 0; i < ROWS_split_on_comma; i++)
		split_on_tilde[i] = (char*) malloc(sizeof(char) * COLS_split_on_comma);

	read_names = (char**) malloc(sizeof(char*) * ROWS_split_on_comma);
	for (i = 0; i < ROWS_split_on_comma; i++)
		read_names[i] = (char*) malloc(sizeof(char) * 100);

	output_prefix_without_path = (char*) malloc(sizeof(char) * MAX_SEQ_LEN);
	sequence_portions_from_reference = (char**) malloc(
			sizeof(char*) * MAX_POOL_SIZE);
	fasta_file_with_expressed_portions = (char*) malloc(
			sizeof(char) * FILENAME_LENGTH);
	current_chromosome = (char*) malloc(sizeof(char) * 100);

	//buffer = ( char* ) malloc (sizeof(char) * BUFFER_SIZE);
	abridge_index = allocateMemoryAbridge_Index();
	sam_alignment = allocateMemorySam_Alignment();
	whole_genome = (struct Whole_Genome_Sequence*) malloc(
			sizeof(struct Whole_Genome_Sequence));

	whole_genome->number_of_reference_sequences = 0;
	whole_genome->nucleotides = (char**) malloc(sizeof(char*) * 1);
	whole_genome->reference_sequence_name = (char**) malloc(sizeof(char*) * 1);
	whole_genome->reference_sequence_length = (unsigned long long int*) malloc(
			sizeof(unsigned long long int) * 1);

	sam_alignment_instance = allocateMemorySam_Alignment();
	read_prefix[0] = '\0'; // Empty string

	printf("\nMemory has been allocated");
	fflush(stdout);
	/********************************************************************/

	writeSequenceHeaders(fhw, genome_filename, 1);

	line_num = 0;
	line_len = getline(&buffer, &len, fhr);
	splitByDelimiter(buffer, '\t', split_on_tab);
	printf("\nThis works");
	fflush(stdout);

	splitByDelimiter(split_on_tab[0], ':', split_on_tilde);
	flag_ignore_mismatches = strtol(
			split_on_tilde[1],
			&convert_to_int_temp,
			10);

	splitByDelimiter(split_on_tab[1], ':', split_on_tilde);
	flag_ignore_soft_clippings = strtol(
			split_on_tilde[1],
			&convert_to_int_temp,
			10);

	splitByDelimiter(split_on_tab[2], ':', split_on_tilde);
	flag_ignore_unmapped_sequences = strtol(
			split_on_tilde[1],
			&convert_to_int_temp,
			10);

	splitByDelimiter(split_on_tab[3], ':', split_on_tilde);
	flag_ignore_quality_scores_for_mismatched_bases_and_soft_clips = strtol(
			split_on_tilde[1],
			&convert_to_int_temp,
			10);

	splitByDelimiter(split_on_tab[4], ':', split_on_tilde);
	flag_ignore_quality_scores_for_matched_bases = strtol(
			split_on_tilde[1],
			&convert_to_int_temp,
			10);

	splitByDelimiter(split_on_tab[6], ':', split_on_tilde);
	flag_ignore_alignment_scores = strtol(
			split_on_tilde[1],
			&convert_to_int_temp,
			10);

	printf("\nflag_ignore_mismatches %d", flag_ignore_mismatches);
	printf("\nflag_ignore_soft_clippings %d", flag_ignore_soft_clippings);
	printf(
			"\nflag_ignore_unmapped_sequences %d",
			flag_ignore_unmapped_sequences);
	printf(
			"\nflag_ignore_quality_scores_for_mismatched_bases_and_soft_clips %d",
			flag_ignore_quality_scores_for_mismatched_bases_and_soft_clips);
	printf(
			"\nflag_ignore_quality_scores_for_matched_bases %d",
			flag_ignore_quality_scores_for_matched_bases);
	printf("\nflag_save_exact_quality_scores %d", flag_ignore_alignment_scores);
	fflush(stdout);

	line_num = 0;
	while ((line_len = getline(&buffer, &len, fhr)) != -1)
	{
		line_num++;

		//printf("\nline_num = %d %s", line_num, buffer);
		//if ( line_num == 10 ) break;
		number_of_columns = splitByDelimiter(buffer, '\t', split_on_tab);

		switch (number_of_columns)
		{
			case 1:
				read_names_stored = 0;
				break;
			case 2:

				if (isNumber(split_on_tab[0]) == 1)
				{
					//its an integer
					read_names_stored = 0;
				}
				else
					read_names_stored = 1;
				break;
			case 3:
				read_names_stored = 1;
				break;
		}

		//printf ("\nRead present? %d %s" , strstr (buffer , "abridge_") , buffer);
		if (buffer[0] == '@')
		{
			splitByDelimiter(split_on_tab[1], ':', split_on_dash); // Using split_on_dash so as to save memory and not create a new data structure
			strcpy(current_chromosome, split_on_dash[1]);
			readInEachChromosome(
					genome_filename,
					whole_genome,
					current_chromosome);
			do
			{
				line_len = getline(&buffer, &len, fhr);
			}
			while (buffer[0] == '@' && line_len != -1);
			curr_position = 0;
			fseek(fhr, -line_len, SEEK_CUR);
			//printf ("\nline_num = %d" , line_num);
			//fflush (stdout);
			continue;
		}

		number_of_commas_in_each_line = 0;
		for (i = 0; buffer[i] != '\0'; i++)
			if (buffer[i] == ',')
				number_of_commas_in_each_line++;
		if (max_number_of_commas < number_of_commas_in_each_line)
			max_number_of_commas = number_of_commas_in_each_line;

		printf("\nCheckpoint 1 line_num = %d", line_num);
		printf("\n2. line_len %d len %d", line_len, len);
		fflush(stdout);
		/*
		 if (line_len > COLS_split_on_tab)
		 {
		 printf(
		 "\nB--> line_len %d COLS_split_on_tab %d",
		 line_len,
		 COLS_split_on_tab);
		 fflush(stdout);
		 COLS_split_on_tab = line_len + 100;
		 for (i = 0; i < ROWS_split_on_tab; i++)
		 {

		 split_on_tab[i] = (char*) realloc(
		 split_on_tab[i],
		 sizeof(char) * COLS_split_on_tab);
		 }
		 printf(
		 "\nA--> line_len %d COLS_split_on_tab %d",
		 line_len,
		 COLS_split_on_tab);
		 fflush(stdout);

		 }
		 printf("\nCheckpoint 2 line_num = %d", line_num);
		 fflush(stdout);
		 if (max_number_of_commas > ROWS_split_on_comma)
		 {
		 //printf ("\nB--> max_number_of_commas %d ROWS_split_on_comma %d" , max_number_of_commas , ROWS_split_on_comma);
		 //fflush (stdout);
		 for (i = 0; i < ROWS_split_on_comma; i++)
		 free(split_on_comma[i]);
		 free(split_on_comma);
		 ROWS_split_on_comma = line_len / max_number_of_commas + 10;
		 split_on_comma = (char**) malloc(
		 sizeof(char*) * ROWS_split_on_comma);
		 for (i = 0; i < ROWS_split_on_comma; i++)
		 split_on_comma[i] = (char*) malloc(
		 sizeof(char) * COLS_split_on_comma);
		 printf(
		 "\nA--> max_number_of_commas %d ROWS_split_on_comma %d",
		 max_number_of_commas,
		 ROWS_split_on_comma);
		 fflush(stdout);
		 }
		 */

		if (read_names_stored == 0)
		{
			if (number_of_columns == 1)
				curr_position++;
			else
				curr_position += strtol(
						split_on_tab[0],
						&convert_to_int_temp,
						10);
		}
		else if (read_names_stored == 1)
		{
			if (number_of_columns == 2)
				curr_position++;
			else
				curr_position += strtol(
						split_on_tab[0],
						&convert_to_int_temp,
						10);
		}
		//printf("\nCheckpoint 3 line_num = %d", line_num);
		/*
		 if (strstr(buffer, "3Zg,2hWh,dpt,DCZz") == NULL)
		 continue;
		 */
		printf("\nThe compressed data %s", buffer);
		convertToAlignmentSingleEnded(
				sam_alignment_instance,
				whole_genome,
				split_on_tab,
				split_on_dash,
				split_on_comma,
				split_on_tilde,
				default_quality_value,
				flag_ignore_alignment_scores,
				flag_ignore_mismatches,
				flag_ignore_soft_clippings,
				flag_ignore_unmapped_sequences,
				flag_ignore_quality_scores_for_mismatched_bases_and_soft_clips,
				flag_ignore_sequence_information,
				&read_number,
				&total_mapped_reads,
				read_prefix,
				fhw,
				fhr_qual,
				flag_ignore_quality_scores_for_matched_bases,
				number_of_columns,
				curr_position,
				current_chromosome,
				read_names,
				read_names_stored);
		//printf ("\nCheckpoint 4 line_num = %d" , line_num);
		//fflush (stdout);
	}

	//printf ("\nTotal mapped reads: %d" , total_mapped_reads);
	//fflush (stdout);
	/*
	 * Write all unmapped reads to samfile
	 */
	fhr = fopen(unmapped_filename, "r");
	if (fhr == NULL)
	{
		printf("Error! File not found");
		exit(1);
	}
	line_number = 1;
	free(buffer);
	buffer = NULL;
	while ((line_len = getline(&buffer, &len, fhr)) != -1)
	{
		strcpy(sam_alignment->seq, buffer);
		line_len = getline(&buffer, &len, fhr);
		strcpy(sam_alignment->qual, buffer);

		line_to_be_written_to_file[0] = '\0';
		strcat(line_to_be_written_to_file, "unmapped_");
		sprintf(temp, "%d", line_number);
		strcat(line_to_be_written_to_file, temp);
		strcat(line_to_be_written_to_file, "\t");

		strcat(line_to_be_written_to_file, "4");
		strcat(line_to_be_written_to_file, "\t");

		strcat(line_to_be_written_to_file, "*");
		strcat(line_to_be_written_to_file, "\t");

		strcat(line_to_be_written_to_file, "0");
		strcat(line_to_be_written_to_file, "\t");

		strcat(line_to_be_written_to_file, "0");

		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, "*");

		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, "*");

		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, "0");

		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, "0");

		strcat(line_to_be_written_to_file, "\t");
		sam_alignment->seq[strlen(sam_alignment->seq) - 1] = '\0';
		strcat(line_to_be_written_to_file, sam_alignment->seq);

		strcat(line_to_be_written_to_file, "\t");
		sam_alignment->qual[strlen(sam_alignment->qual) - 1] = '\0';
		strcat(line_to_be_written_to_file, sam_alignment->qual);

		strcat(line_to_be_written_to_file, "\tNH:i:0\tHI:i:0\tnM:i:1\tuT:A:1");

		strcat(line_to_be_written_to_file, "\n");
		fprintf(fhw, "%s", line_to_be_written_to_file);

		line_number++;
	}

	fclose(fhw);
	fclose(fhr);
}

int main(int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char genome_filename[FILENAME_LENGTH];
	char pass1_filename[FILENAME_LENGTH];
	char output_sam_filename[FILENAME_LENGTH];
	char genome_prefix[FILENAME_LENGTH];
	char name_of_file_with_quality_scores[FILENAME_LENGTH];
	char default_quality_value[10];
	char unmapped_filename[FILENAME_LENGTH];
	char *temp; //Required for strtoi
	int max_reads_in_each_line;

	short int flag_ignore_sequence_information;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy(genome_filename, argv[1]);
	strcpy(output_sam_filename, argv[2]);
	strcpy(pass1_filename, argv[3]);
	strcpy(default_quality_value, argv[4]);
	flag_ignore_sequence_information = strtol(argv[5], &temp, 10);
	strcpy(unmapped_filename, argv[6]);
	strcpy(name_of_file_with_quality_scores, argv[7]);
	max_reads_in_each_line = strtol(argv[8], &temp, 10);

	/********************************************************************/
	decompressFile(
			name_of_file_with_quality_scores,
			genome_filename,
			output_sam_filename,
			pass1_filename,
			unmapped_filename,
			default_quality_value,
			flag_ignore_sequence_information,
			max_reads_in_each_line);
	return 0;
}
