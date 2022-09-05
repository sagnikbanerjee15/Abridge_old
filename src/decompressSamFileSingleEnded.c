# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"
# include <argp.h>

long long unsigned int total_mapped_reads = 0;

// Set up the argument parser
const char *argp_program_version = "abridge decompressSamFileSingleEnded 1.2.0";
const char *argp_program_bug_address = "sagnikbanerjee15@gmail.com";
static char doc[] = "decompressSamFileSingleEnded will accept an compressed file and associated quality scores and decompress those to SAM alignments";
static char args_doc[] = "";  // No standard arguments
							  // (i.e. arguments without "names")

/*
 * Options.  Field 1 in ARGP.
 * Order of fields: {NAME, KEY, ARG, FLAGS, DOC, GROUP}.
 */

static struct argp_option options[] =
{
{ "reference" , 'r' , "referencefilename" , 0 , "Name of the reference file in fasta format" , 0 } ,
{ "outputfilename" , 'o' , "outputfilename" , 0 , "Name of the samfile where the output will be stored" , 0 } ,
{ "compressedfile" , 'c' , "compressedfilename" , 0 , "Name of the compressed file generated by the compression algorithm" , 0 } ,
{ "mockquality" , 'm' , "mockqualityscore" , 0 , "Value of the mock quality score" , 0 } ,
{ "ignoresequence" , 's' , 0 , 0 , "Flag to ignore sequence generation during decompression" , 0 } ,
{ "unmappedreadsfilename" , 'u' , "unmappedreadsfilename" , 0 , "Name of the file with unmapped reads" , 0 } ,
{ "qualityscoresfilename" , 'q' , "qualityscoresfilename" , 0 , "Name of the file with quality scores" , 0 } ,
{ "maxreadsineachline" , 'x' , "maxreadsineachline" , 0 , "Maximum number of reads mapped to a single reference nucleotide position" , 0 } ,
{ 0 , 0 , 0 , 0 , 0 , 0 } // Last entry should be all zeros in all fields
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
	// char *args[0];   // No standard arguments (without flags)
	char *reference;	  // Argument for --reference / -r
	char *outputfilename;  // Argument for --outputfilename / -o
	char *compressedfile;
	char *mockquality;
	int ignoresequence;
	char *unmappedreadsfilename;
	char *qualityscoresfilename;
	char *dictionaryfilename;
	unsigned long long maxreadsineachline;
};

/*
 * Parser. Field 2 in ARGP.
 * Order of parameters: KEY, ARG, STATE.
 * Parse a single option.
 */

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
	/* Get the input argument from argp_parse, which we
	 know is a pointer to our arguments structure. */
	struct arguments *arguments = state->input;
	char *temp;

	// Figure out which option we are parsing, and decide how to store it
	switch ( key )
	{
		case 'r':
			arguments->reference = arg;
			break;
		case 'o':
			arguments->outputfilename = arg;
			break;
		case 'c':
			arguments->compressedfile = arg;
			break;
		case 'm':
			arguments->mockquality = arg;
			break;
		case 's':
			arguments->ignoresequence = 1;
			break;
		case 'u':
			arguments->unmappedreadsfilename = arg;
			break;
		case 'q':
			arguments->qualityscoresfilename = arg;
			break;
		case 'x':
			arguments->maxreadsineachline = strtoull (arg , &temp , 10);
			break;

		case ARGP_KEY_END:
			// Reached the last key.
			// Check if our inputsamfilename and outputfilename REQUIRED "options" have been set to non-default values
			if ( strcmp (arguments->reference , "") == 0 || strcmp (arguments->outputfilename ,
					"") == 0 || strcmp (arguments->compressedfile , "") == 0 || strcmp (arguments->unmappedreadsfilename ,
					"") == 0 || strcmp (arguments->qualityscoresfilename , "") == 0 || arguments->maxreadsineachline == 0 )
			{
				argp_usage (state);
			}
			break;

		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

static struct argp argp =
{ options , parse_opt , args_doc , doc , 0 , 0 , 0 };

void writeAlignmentToFileSingleEnded (
		struct Sam_Alignment *sam_alignment,
		struct Cigar_Items *cigar_items_instance,
		short int flag_ignore_sequence_information,
		int number_of_repititions_of_the_same_reads,
		char *read_prefix,
		FILE *fhw,
		FILE *fhr_qual,
		short int flag_ignore_quality_scores_for_matched_bases,
		char **read_names,
		short int flag_ignore_alignment_scores)
{
	int i, j;
	int read_length_calculated_from_cigar_string = 0;
	int num_of_types;

	char line_to_be_written_to_file[MAX_GENERAL_LEN];
	char temp[100];
	char *buffer;
	/*
	 struct Cigar_Items *cigar_items_instance;
	 cigar_items_instance = (struct Cigar_Items*) malloc(
	 sizeof(struct Cigar_Items) * 100);
	 */
	size_t len = 0;
	ssize_t line_len;

	int read_name_index = 0;
	for ( i = 0 ; i < number_of_repititions_of_the_same_reads ; i++ )
	{
		line_to_be_written_to_file[0] = '\0';
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 0");
		 fflush (stdout);
		 */
		/*
		 * Set up the read name
		 */
		if ( sam_alignment->tags[0].val[0] == '1' && sam_alignment->tags[0].val[1] == '\0' )
		{
			strcat(line_to_be_written_to_file , read_prefix);
			strcat(line_to_be_written_to_file , sam_alignment->read_name);
			sprintf(temp , "%d" , i + 1);
			strcat(line_to_be_written_to_file , "_");
			strcat(line_to_be_written_to_file , temp);
		}
		else
		{
			strcat(line_to_be_written_to_file , read_names[read_name_index]);
			//printf ("\nRead name %s " , read_names[read_name_index]);
			read_name_index++;
		}
		strcat(line_to_be_written_to_file , "\t");
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 01");
		 fflush(stdout);
		 */

		sprintf(temp , "%d" , sam_alignment->samflag);
		strcat(line_to_be_written_to_file , temp);
		strcat(line_to_be_written_to_file , "\t");
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 02");
		 fflush(stdout);
		 */
		strcat(line_to_be_written_to_file , sam_alignment->reference_name);
		strcat(line_to_be_written_to_file , "\t");
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 03");
		 fflush(stdout);
		 */
		sprintf(temp , "%d" , sam_alignment->start_position);
		strcat(line_to_be_written_to_file , temp);
		strcat(line_to_be_written_to_file , "\t");
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 04");
		 fflush(stdout);
		 */
		if ( flag_ignore_alignment_scores == 1 )
			strcat(line_to_be_written_to_file , "255");
		else
		{
			sprintf(temp , "%d" , sam_alignment->mapping_quality_score);
			strcat(line_to_be_written_to_file , temp);
		}
		strcat(line_to_be_written_to_file , "\t");
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 05");
		 fflush(stdout);
		 */
		splitCigar (sam_alignment->cigar ,
				&num_of_types ,
				cigar_items_instance);
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 06");
		 fflush(stdout);

		 printf(
		 "\nwriteAlignmentToFileSingleEnded Checkpoint 2 number_of_repititions_of_the_same_reads %d Value of i %d",
		 number_of_repititions_of_the_same_reads,
		 i);
		 fflush(stdout);
		 */
		strcat(line_to_be_written_to_file , sam_alignment->cigar);
		strcat(line_to_be_written_to_file , "\t");

		strcat(line_to_be_written_to_file , "*");
		strcat(line_to_be_written_to_file , "\t");

		strcat(line_to_be_written_to_file , "0");
		strcat(line_to_be_written_to_file , "\t");

		strcat(line_to_be_written_to_file , "0");
		strcat(line_to_be_written_to_file , "\t");

		strcat(line_to_be_written_to_file , sam_alignment->seq);
		strcat(line_to_be_written_to_file , "\t");
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 3");
		 fflush(stdout);
		 */
		if ( flag_ignore_quality_scores_for_matched_bases == 0 && ( line_len = getline ( &buffer ,
				&len ,
				fhr_qual) ) != -1 )
		{
			if ( buffer[strlen (buffer) - 1] == '\n' )
				buffer[strlen (buffer) - 1] = '\0';
			//strcat(line_to_be_written_to_file, buffer);
			strcpy(sam_alignment->qual , buffer);
			free (buffer);
			buffer = NULL;

		}
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 3.5");
		 fflush(stdout);
		 */

		read_length_calculated_from_cigar_string = 0;
		for ( j = 0 ; j < num_of_types ; j++ )
		{
			if ( cigar_items_instance[j].def != 'N' && cigar_items_instance[j].def != 'D' )
				read_length_calculated_from_cigar_string += cigar_items_instance[j].len;
		}
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 3.75");
		 printf(
		 "\n%d  %d",
		 strlen(sam_alignment->qual),
		 read_length_calculated_from_cigar_string);
		 fflush (stdout);
		 */
		if ( strlen (sam_alignment->qual) != read_length_calculated_from_cigar_string )
		{
			/*
			 printf(
			 "\nread_length_calculated_from_cigar_string %d %s",
			 read_length_calculated_from_cigar_string,
			 sam_alignment->cigar);
			 fflush (stdout);

			 printf(
			 "\nDifference detected %d",
			 strlen(sam_alignment->qual)
			 - read_length_calculated_from_cigar_string);
			 fflush(stdout);
			 */
			sam_alignment->qual[read_length_calculated_from_cigar_string] = '\0';
		}
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 3.9");
		 fflush (stdout);
		 */
		strcat(line_to_be_written_to_file , sam_alignment->qual);
		strcat(line_to_be_written_to_file , "\t");
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 4");
		 fflush(stdout);
		 */
		//Tags
		strcat(line_to_be_written_to_file , "NH:i:");
		strcat(line_to_be_written_to_file , sam_alignment->tags[0].val);
		strcat(line_to_be_written_to_file , "\t");

		if ( strcmp (sam_alignment->tags[1].val , ".") != 0 && strchr (sam_alignment->cigar ,
				'N') != NULL )
		{
			strcat(line_to_be_written_to_file , "XS:A:");
			strcat(line_to_be_written_to_file , sam_alignment->tags[1].val);
			strcat(line_to_be_written_to_file , "\t");
		}
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 5");
		 fflush(stdout);
		 */
		if ( flag_ignore_sequence_information == 0 )
		{
			strcat(line_to_be_written_to_file , "MD:Z:");
			strcat(line_to_be_written_to_file , sam_alignment->tags[2].val);
			strcat(line_to_be_written_to_file , "\t");
		}
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 6");
		 fflush(stdout);
		 */
		if ( flag_ignore_alignment_scores == 0 && strcmp (sam_alignment->tags[3].val ,
				"X") != 0 )
		{
			strcat(line_to_be_written_to_file , "AS:i:");
			strcat(line_to_be_written_to_file , sam_alignment->tags[3].val);
			strcat(line_to_be_written_to_file , "\t");
		}
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 7");
		 fflush(stdout);
		 */
		strcat(line_to_be_written_to_file , "\n");
		fprintf (fhw , "%s" , line_to_be_written_to_file);
		/*
		 printf("\nwriteAlignmentToFileSingleEnded Checkpoint 8");
		 fflush(stdout);
		 */
	}
	/*
	 printf("\nWriting is completed");
	 fflush (stdout);
	 */
}

void convertToAlignmentSingleEnded (
		struct Sam_Alignment *sam_alignment_instance,
		struct Cigar_Items *cigar_items_instance_for_writing_to_file,
		struct Whole_Genome_Sequence *whole_genome,
		char **split_on_tab,
		char **split_on_dash,
		char **split_on_comma,
		char **split_on_tilde,
		char *default_quality_value,
		short int flag_ignore_alignment_scores,
		short int flag_ignore_mismatches,
		short int flag_ignore_soft_clippings,
		short int flag_ignore_unmapped_sequences,
		short int flag_ignore_all_quality_scores,
		short int flag_ignore_sequence_information,
		unsigned long long int *read_number,
		unsigned long long int *total_mapped_reads,
		char *read_prefix,
		FILE *fhw,
		FILE *fhr_qual,
		short int flag_ignore_quality_scores_for_matched_bases,
		short int number_of_columns,
		unsigned long long int curr_position,
		char *chromosome,
		char **read_names,
		short int read_names_stored)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/

	int number_of_distinct_cigars_in_a_line;
	int number_of_repititions_of_the_same_reads;
	int samformatflag;

	int i, j;

	char *temp; //Useless
	char *distinct_icigars_in_a_line;
	char *icigar;
	char str_sprintf[50];

	int read_names_index = 0;

	short int number_of_items_separated_by_underscore;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	/********************************************************************/

	if ( read_names_stored == 1 )
	{
		if ( number_of_columns == 2 )
			read_names_index = splitByDelimiter (split_on_tab[1] ,
					',' ,
					read_names);
		else if ( number_of_columns == 3 )
			read_names_index = splitByDelimiter (split_on_tab[2] ,
					',' ,
					read_names);
		/*
		 printf("\n Reads read");
		 for (i = 0; i < read_names_index; i++)
		 printf("\n%s", read_names[i]);
		 */

	}

	if ( read_names_stored == 0 )
	{
		int max_number_of_commas = 0, number_of_commas = 0;
		switch ( number_of_columns )
		{
			case 1:
			{
				for ( i = 0 ; split_on_tab[0][i] != '\0' ; i++ )
					if ( split_on_tab[0][i] == ',' ) number_of_commas++;

				number_of_distinct_cigars_in_a_line = splitByDelimiter (split_on_tab[0] ,
						',' ,
						split_on_comma);
				/*
				 printf(
				 "\nInside here 01 number_of_distinct_cigars_in_a_line %d",
				 number_of_distinct_cigars_in_a_line);
				 */
			}
				break;
			case 2:
			{
				for ( i = 0 ; split_on_tab[1][i] != '\0' ; i++ )
					if ( split_on_tab[1][i] == ',' ) number_of_commas++;

				number_of_distinct_cigars_in_a_line = splitByDelimiter (split_on_tab[1] ,
						',' ,
						split_on_comma);
				/*
				 printf(
				 "\nInside here 02 number_of_distinct_cigars_in_a_line %d",
				 number_of_distinct_cigars_in_a_line);
				 */
			}
				break;
		}
	}
	else if ( read_names_stored == 1 )
	{
		int max_number_of_commas = 0, number_of_commas = 0;
		switch ( number_of_columns )
		{
			case 2:
			{
				for ( i = 0 ; split_on_tab[0][i] != '\0' ; i++ )
					if ( split_on_tab[0][i] == ',' ) number_of_commas++;

				number_of_distinct_cigars_in_a_line = splitByDelimiter (split_on_tab[0] ,
						',' ,
						split_on_comma);
				/*
				 printf(
				 "\nInside here 12 number_of_distinct_cigars_in_a_line %d",
				 number_of_distinct_cigars_in_a_line);
				 */
			}
				break;
			case 3:
			{
				for ( i = 0 ; split_on_tab[1][i] != '\0' ; i++ )
					if ( split_on_tab[1][i] == ',' ) number_of_commas++;

				number_of_distinct_cigars_in_a_line = splitByDelimiter (split_on_tab[1] ,
						',' ,
						split_on_comma);
				/*
				 printf(
				 "\nInside here 13 number_of_distinct_cigars_in_a_line %d",
				 number_of_distinct_cigars_in_a_line);
				 */
			}
				break;
		}

	}
	/*
	 printf(
	 "\nnumber_of_columns %d read_names_stored %d number_of_distinct_cigars_in_a_line %d",
	 number_of_columns,
	 read_names_stored,
	 number_of_distinct_cigars_in_a_line);
	 */
	//printf("\nconvertToAlignmentSingleEnded Checkpoint 1");
	//fflush (stdout);
	for ( j = 0 ; j < number_of_distinct_cigars_in_a_line ; j++ )
	{
		splitByDelimiter (split_on_comma[j] , '-' , split_on_dash);
		if ( flag_ignore_alignment_scores == 0 )
			splitByDelimiter (split_on_dash[0] , '~' , split_on_tilde);
		else
		strcpy(split_on_tilde[0] , split_on_dash[0]);
		printf ("\nI am here");
		stdout (fflush);
		if ( split_on_comma[j][strlen (split_on_comma[j]) - 1] != '-' )
		{
			number_of_repititions_of_the_same_reads = strtol (split_on_dash[1] ,
					&temp ,
					10);
		}
		else number_of_repititions_of_the_same_reads = 1;

		/*
		 printf(
		 "\nnumber_of_repititions_of_the_same_reads %d split_on_comma %s",
		 number_of_repititions_of_the_same_reads,
		 split_on_comma[j]);

		 printf("\nconvertToAlignmentSingleEnded Checkpoint 2");
		 fflush (stdout);
		 */
		if ( ! ( split_on_comma[j][1] == '-' && isalpha (split_on_dash[0][0]) != 0 ) )
		{
			if ( flag_ignore_alignment_scores == 0 )
			{
				sam_alignment_instance->mapping_quality_score = strtol (split_on_tilde[1] ,
						&temp ,
						10);
				strcpy(sam_alignment_instance->tags[3].val , split_on_tilde[2]);
			}
			else
			{
				sam_alignment_instance->mapping_quality_score = 255;
				strcpy(sam_alignment_instance->tags[3].val , "X");
			}
		}
		/*
		 printf("\nconvertToAlignmentSingleEnded Checkpoint 3");
		 fflush(stdout);
		 */
		//printf ("\n%s %d" , split_on_comma[j] , number_of_repititions_of_the_same_reads);
		sam_alignment_instance->start_position = curr_position;

		if ( split_on_comma[j][1] == '-' && isalpha (split_on_dash[0][0]) != 0 )
		{
			// Use the same cigar
			sprintf(temp , "%d" , *read_number);
			( *read_number )++;
			strcpy(sam_alignment_instance->read_name , temp);
			/*if ( sam_alignment_instance->start_position == 27381 && strcmp (sam_alignment_instance->reference_name , "1") == 0 )
			 {
			 printf ("\nSame iCIGAR");
			 printf ("\nsplit_on_comma[j] = %s" , split_on_comma[j]);
			 printf ("\nWeird Location cigar %s number_of_repititions_of_the_same_reads %d" , sam_alignment_instance->icigar , number_of_repititions_of_the_same_reads);
			 printf ("\nMD String: %s" , sam_alignment_instance->tags[2].val);
			 printf ("\n==============================================================================================================================");
			 fflush (stdout);
			 }*/
			/*
			 printf("\nconvertToAlignmentSingleEnded Checkpoint 4");
			 fflush(stdout);
			 */
		}
		else
		{
			strcpy(sam_alignment_instance->icigar , split_on_tilde[0]);
			//printf ("\nj=%d number_of_distinct_cigars_in_a_line=%d Inside ICIGAR %s" , j , number_of_distinct_cigars_in_a_line , sam_alignment_instance->icigar);
			//fflush (stdout);
			//printf ("\nConvertion started");
			convertIcigarToCigarandMDSingleEnded (whole_genome ,
					sam_alignment_instance ,
					chromosome ,
					flag_ignore_mismatches ,
					flag_ignore_soft_clippings ,
					flag_ignore_unmapped_sequences ,
					flag_ignore_all_quality_scores ,
					flag_ignore_quality_scores_for_matched_bases ,
					flag_ignore_sequence_information ,
					default_quality_value);
			/*
			 printf ("\nConvertion completed");
			 fflush (stdout);
			 */
			/*if ( sam_alignment_instance->start_position == 27381 && strcmp (sam_alignment_instance->reference_name , "1") == 0 )
			 {
			 printf ("\nsplit_on_comma[j] = %s" , split_on_comma[j]);
			 printf ("\nWeird Location cigar %s number_of_repititions_of_the_same_reads %d" , sam_alignment_instance->icigar , number_of_repititions_of_the_same_reads);
			 printf ("\nMD String: %s" , sam_alignment_instance->tags[2].val);
			 printf ("\n split_on_dash[0] %s split_on_dash[1] %s" , split_on_dash[0] , split_on_dash[1]);
			 printf ("\nsplit_on_dash[0][1] == '-' %d isalpha (split_on_dash[0][0])  %d" , split_on_dash[0][1] == '-' , isalpha (split_on_dash[0][0]));
			 printf ("\n==============================================================================================================================");
			 fflush (stdout);
			 }*/
			/*printf ("\nread_number = %d" , *read_number);
			 fflush (stdout);*/

			sprintf(str_sprintf , "%llu" , *read_number);

			/*printf ("\ntemp = %s" , str_sprintf);
			 fflush (stdout);*/
			( *read_number )++;
			/*printf ("\ncopying into read name");
			 fflush (stdout);*/
			strcpy(sam_alignment_instance->read_name , str_sprintf);

			/*printf ("\nconvertToAlignmentSingleEnded Checkpoint 5");
			 fflush (stdout);*/

		}
		fflush (stdout);
		//printSamAlignmentInstance (sam_alignment_instance , 0);
		//continue;
		/*
		 printf(
		 "\nconvertToAlignmentSingleEnded Checkpoint 6 Value of j %d number_of_distinct_cigars_in_a_line %d",
		 j,
		 number_of_distinct_cigars_in_a_line);
		 fflush(stdout);
		 */

		writeAlignmentToFileSingleEnded (sam_alignment_instance ,
				cigar_items_instance_for_writing_to_file ,
				flag_ignore_sequence_information ,
				number_of_repititions_of_the_same_reads ,
				read_prefix ,
				fhw ,
				fhr_qual ,
				flag_ignore_quality_scores_for_matched_bases ,
				read_names ,
				flag_ignore_alignment_scores);
		fflush (fhw);
		//(*total_mapped_reads) += number_of_repititions_of_the_same_reads;

		printf ("\nconvertToAlignmentSingleEnded Checkpoint 7");
		fflush (stdout);
	}
	printf ("\nLeaving function convertToAlignmentSingleEnded");
	fflush (stdout);
}

void decompressFile (
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
	short int flag_ignore_all_quality_scores;
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
	int ROWS_split_on_tab = 10; //10
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
	struct Cigar_Items *cigar_items_instance_for_writing;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen (pass1_filename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File %s not found" , pass1_filename);
		exit (1);
	}
	fhr_qual = fopen (name_of_file_with_quality_scores , "r");
	if ( fhr_qual == NULL )
	{
		printf ("Error! File %s not found" , name_of_file_with_quality_scores);
		exit (1);
	}
	fhw = fopen (output_sam_filename , "w");
	if ( fhw == NULL )
	{
		printf ("Error! File %s cannot be opened for writing" ,
				output_sam_filename);
		exit (1);
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
	while ( ( line_len = getline ( &buffer , &len , fhr) ) != -1 )
	{
		if ( buffer[0] != '@' )
		{
			if ( max_line_len < strlen (buffer) )
				max_line_len = strlen (buffer);

			int num_commas = 0;
			for ( int i = 0 ; buffer[i] != '\0' ; i++ )
				if ( buffer[i] == ',' ) num_commas += 1;
			if ( max_commas < num_commas ) max_commas = num_commas;
		}
	}

	COLS_split_on_tab = max_line_len;
	ROWS_split_on_comma = max_commas;
	fclose (fhr);
	fhr = fopen (pass1_filename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File %s not found" , pass1_filename);
		exit (1);
	}
	/*
	 printf(
	 "\nAm here line_len %d COLS_split_on_tab %d ROWS_split_on_comma %d max_commas*1000 %d",
	 max_line_len,
	 COLS_split_on_tab,
	 ROWS_split_on_comma,
	 max_commas * 1000);
	 fflush(stdout);
	 */

	split_on_tab = ( char** ) malloc (sizeof(char*) * ROWS_split_on_tab);
	for ( i = 0 ; i < ROWS_split_on_tab ; i++ )
		split_on_tab[i] = ( char* ) malloc (sizeof(char) * COLS_split_on_tab);

	split_on_dash = ( char** ) malloc (sizeof(char*) * ROWS_split_on_dash);
	for ( i = 0 ; i < ROWS_split_on_dash ; i++ )
		split_on_dash[i] = ( char* ) malloc (sizeof(char) * COLS_split_on_dash);

	split_on_comma = ( char** ) malloc (sizeof(char*) * ROWS_split_on_comma);
	for ( i = 0 ; i < ROWS_split_on_comma ; i++ )
		split_on_comma[i] = ( char* ) malloc (sizeof(char) * COLS_split_on_comma);

	split_on_tilde = ( char** ) malloc (sizeof(char*) * ROWS_split_on_comma);
	for ( i = 0 ; i < ROWS_split_on_comma ; i++ )
		split_on_tilde[i] = ( char* ) malloc (sizeof(char) * COLS_split_on_comma);

	read_names = ( char** ) malloc (sizeof(char*) * ROWS_split_on_comma);
	for ( i = 0 ; i < ROWS_split_on_comma ; i++ )
		read_names[i] = ( char* ) malloc (sizeof(char) * 100);

	output_prefix_without_path = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	sequence_portions_from_reference = ( char** ) malloc (sizeof(char*) * MAX_POOL_SIZE);
	fasta_file_with_expressed_portions = ( char* ) malloc (sizeof(char) * FILENAME_LENGTH);
	current_chromosome = ( char* ) malloc (sizeof(char) * 100);

	cigar_items_instance_for_writing = ( struct Cigar_Items* ) malloc (sizeof(struct Cigar_Items) * 100);

	//buffer = ( char* ) malloc (sizeof(char) * BUFFER_SIZE);
	abridge_index = allocateMemoryAbridge_Index ();
	sam_alignment = allocateMemorySam_Alignment ();
	whole_genome = ( struct Whole_Genome_Sequence* ) malloc (sizeof(struct Whole_Genome_Sequence));

	whole_genome->number_of_reference_sequences = 0;
	whole_genome->nucleotides = ( char** ) malloc (sizeof(char*) * 1);
	whole_genome->reference_sequence_name = ( char** ) malloc (sizeof(char*) * 1);
	whole_genome->reference_sequence_length = ( unsigned long long int* ) malloc (sizeof(unsigned long long int) * 1);

	sam_alignment_instance = allocateMemorySam_Alignment ();
	read_prefix[0] = '\0'; // Empty string
	/*
	 printf("\nMemory has been allocated");
	 fflush(stdout);
	 */
	/********************************************************************/

	writeSequenceHeaders (fhw , genome_filename , 1);

	line_num = 0;
	line_len = getline ( &buffer , &len , fhr);
	/*
	 printf("\nThis works");
	 fflush(stdout);
	 */
	splitByDelimiter (buffer , '\t' , split_on_tab);

	splitByDelimiter (split_on_tab[0] , ':' , split_on_tilde);
	flag_ignore_mismatches = strtol (split_on_tilde[1] ,
			&convert_to_int_temp ,
			10);

	splitByDelimiter (split_on_tab[1] , ':' , split_on_tilde);
	flag_ignore_soft_clippings = strtol (split_on_tilde[1] ,
			&convert_to_int_temp ,
			10);

	splitByDelimiter (split_on_tab[2] , ':' , split_on_tilde);
	flag_ignore_unmapped_sequences = strtol (split_on_tilde[1] ,
			&convert_to_int_temp ,
			10);

	splitByDelimiter (split_on_tab[3] , ':' , split_on_tilde);
	flag_ignore_all_quality_scores = strtol (split_on_tilde[1] ,
			&convert_to_int_temp ,
			10);

	splitByDelimiter (split_on_tab[4] , ':' , split_on_tilde);
	flag_ignore_quality_scores_for_matched_bases = strtol (split_on_tilde[1] ,
			&convert_to_int_temp ,
			10);

	splitByDelimiter (split_on_tab[5] , ':' , split_on_tilde);
	flag_ignore_alignment_scores = strtol (split_on_tilde[1] ,
			&convert_to_int_temp ,
			10);

	/*
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
	 */
	line_num = 0;
	while ( ( line_len = getline ( &buffer , &len , fhr) ) != -1 )
	{
		line_num++;
		//if ( line_num == 5 ) return;
		/*
		 printf("\nStarting new line");
		 printf("\nline_num = %d %d %s", line_num, strlen(buffer), buffer);
		 fflush(stdout);
		 */
		//if ( line_num == 10 ) break;
		number_of_columns = splitByDelimiter (buffer , '\t' , split_on_tab);

		switch ( number_of_columns )
		{
			case 1:
				read_names_stored = 0;
				break;
			case 2:

				if ( isNumber (split_on_tab[0]) == 1 ) // Change this - check for - in the column 1
				{
					//its an integer
					read_names_stored = 0;
				}
				else read_names_stored = 1;
				break;
			case 3:
				read_names_stored = 1;
				break;
		}

		if ( buffer[0] == '@' )
		{
			splitByDelimiter (split_on_tab[1] , ':' , split_on_dash); // Using split_on_dash so as to save memory and not create a new data structure
			strcpy(current_chromosome , split_on_dash[1]);
			readInEachChromosome (genome_filename ,
					whole_genome ,
					current_chromosome);
			do
			{
				line_len = getline ( &buffer , &len , fhr);
			} while ( buffer[0] == '@' && line_len != -1 );
			curr_position = 0;
			fseek (fhr , -line_len , SEEK_CUR);
			//printf ("\nline_num = %d" , line_num);
			//fflush (stdout);
			continue;
		}
		/*
		 number_of_commas_in_each_line = 0;
		 for (i = 0; buffer[i] != '\0'; i++)
		 if (buffer[i] == ',')
		 number_of_commas_in_each_line++;
		 if (max_number_of_commas < number_of_commas_in_each_line)
		 max_number_of_commas = number_of_commas_in_each_line;
		 */
		/*
		 printf("\nCheckpoint 1 line_num = %d", line_num);
		 printf("\n2. line_len %d len %d", line_len, len);
		 fflush(stdout);
		 */
		if ( read_names_stored == 0 )
		{
			if ( number_of_columns == 1 )
				curr_position++;
			else curr_position += strtol (split_on_tab[0] ,
					&convert_to_int_temp ,
					10);
		}
		else if ( read_names_stored == 1 )
		{
			if ( number_of_columns == 2 )
				curr_position++;
			else curr_position += strtol (split_on_tab[0] ,
					&convert_to_int_temp ,
					10);
		}
		/*
		 printf("\nCheckpoint 3 line_num = %d", line_num);
		 fflush(stdout);
		 */
		/*
		 if (strstr(buffer, "3Zg,2hWh,dpt,DCZz") == NULL)
		 continue;
		 */

		convertToAlignmentSingleEnded (sam_alignment_instance ,
				cigar_items_instance_for_writing ,
				whole_genome ,
				split_on_tab ,
				split_on_dash ,
				split_on_comma ,
				split_on_tilde ,
				default_quality_value ,
				flag_ignore_alignment_scores ,
				flag_ignore_mismatches ,
				flag_ignore_soft_clippings ,
				flag_ignore_unmapped_sequences ,
				flag_ignore_all_quality_scores ,
				flag_ignore_sequence_information ,
				&read_number ,
				&total_mapped_reads ,
				read_prefix ,
				fhw ,
				fhr_qual ,
				flag_ignore_quality_scores_for_matched_bases ,
				number_of_columns ,
				curr_position ,
				current_chromosome ,
				read_names ,
				read_names_stored);
		/*
		 printf("\nCheckpoint 4 line_num = %d", line_num);
		 fflush(stdout);
		 */
		free (buffer);
		buffer = NULL;
	}

	//printf ("\nTotal mapped reads: %d" , total_mapped_reads);
	//fflush (stdout);
	/*
	 * Write all unmapped reads to samfile
	 */
	fhr = fopen (unmapped_filename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File not found");
		exit (1);
	}
	line_number = 1;
	free (buffer);
	buffer = NULL;
	while ( ( line_len = getline ( &buffer , &len , fhr) ) != -1 )
	{
		strcpy(sam_alignment->seq , buffer);
		line_len = getline ( &buffer , &len , fhr);
		strcpy(sam_alignment->qual , buffer);

		line_to_be_written_to_file[0] = '\0';
		strcat(line_to_be_written_to_file , "unmapped_");
		sprintf(temp , "%d" , line_number);
		strcat(line_to_be_written_to_file , temp);
		strcat(line_to_be_written_to_file , "\t");

		strcat(line_to_be_written_to_file , "4");
		strcat(line_to_be_written_to_file , "\t");

		strcat(line_to_be_written_to_file , "*");
		strcat(line_to_be_written_to_file , "\t");

		strcat(line_to_be_written_to_file , "0");
		strcat(line_to_be_written_to_file , "\t");

		strcat(line_to_be_written_to_file , "0");

		strcat(line_to_be_written_to_file , "\t");
		strcat(line_to_be_written_to_file , "*");

		strcat(line_to_be_written_to_file , "\t");
		strcat(line_to_be_written_to_file , "*");

		strcat(line_to_be_written_to_file , "\t");
		strcat(line_to_be_written_to_file , "0");

		strcat(line_to_be_written_to_file , "\t");
		strcat(line_to_be_written_to_file , "0");

		strcat(line_to_be_written_to_file , "\t");
		sam_alignment->seq[strlen (sam_alignment->seq) - 1] = '\0';
		strcat(line_to_be_written_to_file , sam_alignment->seq);

		strcat(line_to_be_written_to_file , "\t");
		sam_alignment->qual[strlen (sam_alignment->qual) - 1] = '\0';
		strcat(line_to_be_written_to_file , sam_alignment->qual);

		strcat(line_to_be_written_to_file , "\tNH:i:0\tHI:i:0\tnM:i:1\tuT:A:1");

		strcat(line_to_be_written_to_file , "\n");
		fprintf (fhw , "%s" , line_to_be_written_to_file);

		line_number++;
	}

	fclose (fhw);
	fclose (fhr);
}

int main (int argc, char *argv[])
{
	struct arguments arguments;

	// Parse our arguments; every option seen by parse_opt will be reflected in arguments.
	// Default values.
	arguments.reference = ""; // Empty string - only contains null character
	arguments.outputfilename = "";
	arguments.compressedfile = "";
	arguments.ignoresequence = 0;
	arguments.mockquality = "I";
	arguments.qualityscoresfilename = "";
	arguments.unmappedreadsfilename = "";
	arguments.maxreadsineachline = 0;

	argp_parse ( &argp , argc , argv , 0 , 0 , &arguments);

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
	strcpy(genome_filename , arguments.reference);
	strcpy(output_sam_filename , arguments.outputfilename);
	strcpy(pass1_filename , arguments.compressedfile);
	strcpy(default_quality_value , arguments.mockquality);
	flag_ignore_sequence_information = arguments.ignoresequence;
	strcpy(unmapped_filename , arguments.unmappedreadsfilename);
	strcpy(name_of_file_with_quality_scores , arguments.qualityscoresfilename);
	max_reads_in_each_line = arguments.maxreadsineachline;

	/********************************************************************/
	decompressFile (name_of_file_with_quality_scores ,
			genome_filename ,
			output_sam_filename ,
			pass1_filename ,
			unmapped_filename ,
			default_quality_value ,
			flag_ignore_sequence_information ,
			max_reads_in_each_line);
	return 0;
}
