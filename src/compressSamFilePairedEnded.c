# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include <stdbool.h>
# include <argp.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

// Set up the argument parser
const char *argp_program_version = "abridge compressSamFileSingleEnded 1.2.0";
const char *argp_program_bug_address = "sagnikbanerjee15@gmail.com";
static char doc[] = "compressSamFileSingleEnded will accept an alignment file in SAM format and remove all redundant information. It will preserve only the information that has been requested by the user.";
static char args_doc[] = "";  // No standard arguments
							  // (i.e. arguments without "names")

/*
 * Options.  Field 1 in ARGP.
 * Order of fields: {NAME, KEY, ARG, FLAGS, DOC, GROUP}.
 */

static struct argp_option options[] =
{
{ "genome_filename", 'g', "GENOME_FILENAME", 0, "Enter the name of the genome file in fasta format", 0 },
{ "input_sam_filename", 'i', "SAM_FILENAME", 0, "Enter the name of the SAM file to be compressed", 0 },
{ "output_abridge_filename", 'o', "TEXT_FILENAME", 0, "Enter the name of the compressed file (please note that this is not the final compressed file)", 0 },
{ "unmapped_filename", 'u', "UNMAPPED_READS_FILENAME", 0, "Enter the name of the file where the unmapped reads will be stored", 0 },
{ "name_of_file_with_max_commas", 'c', "MAX_COMMAS_FILENAME", 0, "Enter the name of the file that contains the value of maximum number of commas", 0 },
{ "name_of_file_with_quality_scores", 'q', "QUALITY_SCORES_FILENAME", 0, "Enter the name of the file where the quality scores will be stored. This file will be compressed later", 0 },
{ "name_of_file_with_read_names_to_short_read_names_and_NH", 'r', "SHORT_NAMES_NH_FILENAME", 0, "Enter the name of the file that contains the mapping between the long name to the short name and the NH values", 0 },

{ "frequency_of_flags_filename", 'j', "FREQUENCY_OF_FLAGS_FILENAME", 0, "Enter the name of the file which has the frequency of flags", 0 },
{ "name_of_file_dictionary", 't', "DICTIONARY_FILENAME", 0, "Enter the name of the dictionary", 0 },

{ "flag_ignore_soft_clippings", 's', 0, 0, "Set this flag to ignore soft clippings", 0 },
{ "flag_ignore_mismatches", 'm', 0, 0, "Set this flag to ignore mismatches", 0 },
{ "flag_ignore_all_quality_scores", 'p', 0, 0, "Set this flag to ignore quality scores for mismatched bases and soft clips", 0 },
{ "flag_ignore_unmapped_sequences", 'e', 0, 0, "Set this flag to ignore unmapped sequences along with their quality scores", 0 },
{ "flag_ignore_quality_scores_for_matched_bases", 'b', 0, 0, "Set this flag to ignore quality scores for nucleotide bases that match to the provided reference", 0 },
{ "flag_ignore_alignment_scores", 'a', 0, 0, "Set this flag to ignore the alignment scores (Column 5 of SAM file)", 0 },
{ "skip_shortening_read_names", 'f', 0, 0, "Set this flag to skip shortening read names", 0 },

{ "max_input_reads_in_a_single_nucl_loc", 'n', "MAX_READS_IN_ONE_NUCL", 0, "Enter the value of the maximum number of input reads mapped to a single nucleotide", 0 },
{ "run_diagnostics", 'd', 0, 0, "Set this flag to run diagnostics and print out a verbose report", 0 },

{ "max_number_of_alignments", 'k', "MAX_NUM_OF_ALIGNMENTS", 0, "Enter the value of the maximum number of alignments", 0 },
{ "max_read_length", 'l', "MAX_READS_LENGTH", 0, "Enter the value of the maximum length of read", 0 },

{ 0, 0, 0, 0, 0, 0 } // Last entry should be all zeros in all fields
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
	// char *args[0];   // No standard arguments (without flags)
	char *input_sam_filename; // Empty string - only contains null character
	char *output_abridge_filename;
	char *genome_filename;
	char *unmapped_filename;
	char *name_of_file_with_max_commas;
	char *name_of_file_with_quality_scores;
	char *name_of_file_with_read_names_to_short_read_names_and_NH;
	int flag_ignore_soft_clippings;
	int flag_ignore_mismatches;
	int flag_ignore_all_quality_scores;
	int flag_ignore_unmapped_sequences;
	int skip_shortening_read_names;
	int run_diagnostics;
	int flag_ignore_quality_scores_for_matched_bases;
	int flag_ignore_alignment_scores;
	unsigned long long int max_input_reads_in_a_single_nucl_loc;

	char *frequency_of_flags_filename;
	char *name_of_file_dictionary;

	unsigned long long int max_number_of_alignments;
	unsigned long long int max_read_length;
};

/*
 * Parser. Field 2 in ARGP.
 * Order of parameters: KEY, ARG, STATE.
 * Parse a single option.
 */

static error_t parse_opt( int key, char *arg, struct argp_state *state )
{
	/* Get the input argument from argp_parse, which we
	 know is a pointer to our arguments structure. */
	struct arguments *arguments = state->input;
	char *eptr;

	// Figure out which option we are parsing, and decide how to store it
	switch ( key )
	{
		case 'a':
			arguments->flag_ignore_alignment_scores = 1;
			break;
		case 'b':
			arguments->flag_ignore_quality_scores_for_matched_bases = 1;
			break;
		case 'c':
			arguments->name_of_file_with_max_commas = arg;
			break;
		case 'd':
			arguments->run_diagnostics = 1;
			break;
		case 'e':
			arguments->flag_ignore_unmapped_sequences = 1;
			break;
		case 'f':
			arguments->skip_shortening_read_names = 1;
			break;
		case 'g':
			arguments->genome_filename = arg;
			break;
		case 'i':
			arguments->input_sam_filename = arg;
			break;
		case 'j':
			arguments->frequency_of_flags_filename = arg;
			break;
		case 'k':
			arguments->max_number_of_alignments = strtoull( arg, &eptr, 10 );
			break;
		case 'l':
			arguments->max_read_length = strtoull( arg, &eptr, 10 );
			break;
		case 'm':
			arguments->flag_ignore_mismatches = 1;
			break;
		case 'n':
			arguments->max_input_reads_in_a_single_nucl_loc = strtoull(
					arg,
					&eptr,
					10 );
			break;
		case 'o':
			arguments->output_abridge_filename = arg;
			break;
		case 'p':
			arguments->flag_ignore_all_quality_scores = 1;
			break;
		case 'q':
			arguments->name_of_file_with_quality_scores = arg;
			break;
		case 'r':
			arguments->name_of_file_with_read_names_to_short_read_names_and_NH = arg;
			break;
		case 's':
			arguments->flag_ignore_soft_clippings = 1;
			break;
		case 't':
			arguments->name_of_file_dictionary = arg;
			break;
		case 'u':
			arguments->unmapped_filename = arg;
			break;

		case ARGP_KEY_END:
			// Reached the last key.
			// Check if our inputsamfilename and outputfilename REQUIRED "options" have been set to non-default values
			if ( strcmp( arguments->input_sam_filename, "" ) == 0
					|| strcmp( arguments->output_abridge_filename, "" ) == 0
					|| strcmp( arguments->genome_filename, "" ) == 0
					|| strcmp( arguments->unmapped_filename, "" ) == 0
					|| strcmp( arguments->name_of_file_with_max_commas, "" )
							== 0
					|| strcmp( arguments->name_of_file_with_quality_scores, "" )
							== 0
					|| strcmp(
							arguments->name_of_file_with_read_names_to_short_read_names_and_NH,
							"" ) == 0
					|| arguments->max_input_reads_in_a_single_nucl_loc == 0 )
			{
				argp_usage( state );
			}
			break;

		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

// Our argp parser.
static struct argp argp =
{ options, parse_opt, args_doc, doc, 0, 0, 0 };

char findMatchCharacterIcigar(
		char *icigar,
		char samformatflag_replacer_characters[] )
{
	int i;

	for ( i = 0; icigar[i] != '\0'; i++ )
	{
		if ( strchr( samformatflag_replacer_characters, icigar[i] ) != NULL )
			return icigar[i];
	}
	return ' ';
}

void generateNextReadID( char *alphabets, int *read_id, int *read_length )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	/********************************************************************/

	if ( *read_length == 0 )
	{
		read_id[0] = 0;
		(*read_length)++;
	}
	else
	{
		/*
		 * Check if the read is the last element of the maximum read_length
		 */
		for ( i = 0; i < *read_length; i++ )
			if ( read_id[i] != strlen( alphabets ) - 1 )
				break;
		if ( i == *read_length )
		{
			(*read_length)++;
			for ( i = 0; i < *read_length; i++ )
				read_id[i] = 0;
		}
		else
		{
			/*
			 * Increment the read_id
			 */
			for ( i = *read_length - 1; i >= 0; i-- )
			{
				if ( read_id[i] == (strlen( alphabets ) - 1) )
					read_id[i] = 0;
				else
				{
					read_id[i]++;
					break;
				}
			}
		}
	}
}

void convertReadIdToString(
		int *read_id,
		char *read_id_string,
		int read_length,
		char *alphabets )
{
	int i;
	for ( i = 0; i < read_length; i++ )
		read_id_string[i] = alphabets[read_id[i]];
	read_id_string[i] = '\0';
}

void writeToFile(
		short int flag_ignore_quality_scores_for_matched_bases,
		short int flag_ignore_all_quality_scores,
		FILE *fhw_qual,
		FILE *fhw_pass1,
		struct Compressed_DS **compressed_ds_pool,
		int compressed_ds_pool_total,
		char *write_to_file_col1,
		char *write_to_file_col2,
		char *write_to_file_col3,
		char *encoded_string,
		long long int *count,
		char **qual_Scores,
		int quality_score_index,
		char samformatflag_replacer_characters[],
		int number_of_unique_samformatflags,
		struct Paired_Ended_Flag_to_Single_Character *samflag_dictionary,
		short int flag_ignore_soft_clippings,
		struct Cigar_Items *cigar_items_instance,
		char *line_to_be_written_to_file,
		char *list_of_read_names,
		char *list_of_qual_scores,
		char *qual )
{
	int i, j, k, l;
	char str[1000];
	//char qual[MAX_SEQ_LEN];
	//char line_to_be_written_to_file[MAX_LINE_TO_BE_WRITTEN_TO_FILE];
	//char list_of_read_names[MAX_LINE_TO_BE_WRITTEN_TO_FILE];
	int num_of_types;
	//char investigate_qual[1000] = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFFFFFFFFFFF,FFFFFFFFF:FFFFFFFF:FFFFFFFFFF";

	line_to_be_written_to_file[0] = '\0';
	list_of_read_names[0] = '\0';
	for ( i = 0; i < compressed_ds_pool_total; i++ )
	{
		if ( i == 0 )
		{
			if ( compressed_ds_pool[i]->position != 1 )
			{
				sprintf( str, "%lld", compressed_ds_pool[i]->position );
				strcat( line_to_be_written_to_file, str );
				strcat( line_to_be_written_to_file, "\t" );
			}
			else
				str[0] = '\0'; // empty string
		}
		strcat( line_to_be_written_to_file, compressed_ds_pool[i]->icigar );
		strcat( line_to_be_written_to_file, "-" );
		if ( compressed_ds_pool[i]->num_reads > 1 )
		{
			sprintf( str, "%ld", compressed_ds_pool[i]->num_reads );
			strcat( line_to_be_written_to_file, str );
		}
		if ( i != compressed_ds_pool_total - 1 )
			strcat( line_to_be_written_to_file, "," );

		for ( j = 0; j < compressed_ds_pool[i]->num_reads; j++ )
		{
			strcat(
					list_of_read_names,
					compressed_ds_pool[i]->pointers_to_read_names[j] );
			if ( i != compressed_ds_pool_total - 1
					|| j != compressed_ds_pool[i]->num_reads - 1 )
				strcat( list_of_read_names, "," );
		}

		if ( flag_ignore_all_quality_scores == 0
				&& flag_ignore_quality_scores_for_matched_bases == 0 )
		{
			for ( j = 0; j < compressed_ds_pool[i]->num_reads; j++ )
			{
				/*
				 // Reversing quality scores if the read was mapped to the reverse strand
				 if ( returnDirection (findMatchCharacterIcigar (compressed_ds_pool[i]->icigar ,
				 samformatflag_replacer_characters) ,
				 samflag_dictionary ,
				 number_of_unique_samformatflags) == '+' )
				 {
				 for ( k = 0 ;
				 compressed_ds_pool[i]->pointers_to_qual_scores[j][k] != '\0' ;
				 k++ )
				 qual[k] = compressed_ds_pool[i]->pointers_to_qual_scores[j][k] - 90;
				 qual[k] = '\0';
				 }
				 else if ( returnDirection (findMatchCharacterIcigar (compressed_ds_pool[i]->icigar ,
				 samformatflag_replacer_characters) ,
				 samflag_dictionary ,
				 number_of_unique_samformatflags) == '-' )
				 {
				 for ( k = strlen (compressed_ds_pool[i]->pointers_to_qual_scores[j]) - 1 ;
				 k >= 0 ; k-- )
				 qual[strlen (compressed_ds_pool[i]->pointers_to_qual_scores[j]) - 1 - k] = compressed_ds_pool[i]->pointers_to_qual_scores[j][k] - 90;
				 qual[strlen (compressed_ds_pool[i]->pointers_to_qual_scores[j])] = '\0';
				 }*/

				for ( k = 0;
						compressed_ds_pool[i]->pointers_to_qual_scores[j][k]
								!= '\0'; k++ )
					qual[k] = compressed_ds_pool[i]->pointers_to_qual_scores[j][k]
							- 90;
				qual[k] = '\0';
				fprintf( fhw_qual, "%s", "\n" );
				fprintf( fhw_qual, "%s", "\n" );
				fprintf( fhw_qual, "%s", "\n" );
				fprintf( fhw_qual, "%s", qual );
				fprintf( fhw_qual, "%s", "\n" );
			}
		}
	}
	strcat( line_to_be_written_to_file, "\t" );
	strcat( line_to_be_written_to_file, list_of_read_names );
	strcat( line_to_be_written_to_file, "\n" );
	fprintf( fhw_pass1, "%s", line_to_be_written_to_file );
	*count = compressed_ds_pool_total;
}

void prepareIcigarForComparison(
		char *icigar1,
		char *icigar,
		char samformatflag_replacer_characters[] )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i, j, k;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	for ( i = 0; icigar[i] != '\0'; i++ )
	{
		if ( strchr( samformatflag_replacer_characters, icigar[i] ) == NULL )
			icigar1[i] = icigar[i];
		else
			icigar1[i] = 'M';
	}
	icigar1[i] = '\0';
}

void reModeliCIGARSPairedEnded(
		struct Compressed_DS **compressed_ds_pool,
		struct Compressed_DS **compressed_ds_pool_rearranged,
		short *already_processed,
		int compressed_ds_pool_index,
		char **modified_icigars,
		char samformatflag_replacer_characters[] )
{
	//printf ("\nInside reModeliCIGARSPairedEnded");
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i, j, k;
	int compressed_ds_pool_rearranged_index = 0;

	//char icigar1[MAX_SEQ_LEN];
	//char icigar2[MAX_SEQ_LEN];
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	for ( i = 0; i < compressed_ds_pool_index; i++ )
		already_processed[i] = 0;

	/********************************************************************/
	for ( i = 0; i < compressed_ds_pool_index; i++ )
		prepareIcigarForComparison(
				modified_icigars[i],
				compressed_ds_pool[i]->icigar,
				samformatflag_replacer_characters );

	for ( i = 0; i < compressed_ds_pool_index; i++ )
	{
		if ( already_processed[i] == 1 )
			continue;
		already_processed[i] = 1;
		//prepareIcigarForComparison ( icigar1 , compressed_ds_pool[i]->icigar,samformatflag_replacer_characters );
		/*
		 * Copy the icigar entry into the rearranged pool
		 */
		//icigar1 = modified_icigars[i];
		strcpy(
				compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->icigar,
				compressed_ds_pool[i]->icigar );
		strcpy(
				compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->cigar,
				compressed_ds_pool[i]->cigar );
		compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->num_reads = compressed_ds_pool[i]->num_reads;
		compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->position = compressed_ds_pool[i]->position;
		for ( k = 0; k < compressed_ds_pool[i]->num_reads; k++ )
		{
			compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->pointers_to_qual_scores[k] = compressed_ds_pool[i]->pointers_to_qual_scores[k];
			compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->pointers_to_read_names[k] = compressed_ds_pool[i]->pointers_to_read_names[k];
		}
		compressed_ds_pool_rearranged_index++;

		for ( j = i + 1; j < compressed_ds_pool_index; j++ )
		{
			if ( already_processed[j] == 1 )
				continue;
			//prepareIcigarForComparison ( icigar2 , compressed_ds_pool[j]->icigar );
			//icigar2 = modified_icigars[j];
			if ( strcmp( modified_icigars[i], modified_icigars[j] ) == 0 )
			{
				//printf ( "\n%s %s %d" , modified_icigars[i] , modified_icigars[j] , strcmp ( modified_icigars[i] , modified_icigars[j] ) );
				//printf ( "\n%s %s %d" , compressed_ds_pool[i]->icigar , compressed_ds_pool[j]->icigar , strcmp ( compressed_ds_pool[i]->icigar , compressed_ds_pool[j]->icigar ) );
				//strcpy( compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->icigar , compressed_ds_pool[i]->icigar );
				compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->icigar[0] = findMatchCharacterIcigar(
						compressed_ds_pool[j]->icigar,
						samformatflag_replacer_characters );
				compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->icigar[1] = '\0';
				strcpy(
						compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->cigar,
						compressed_ds_pool[j]->cigar );
				compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->num_reads = compressed_ds_pool[j]->num_reads;
				compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->position = compressed_ds_pool[j]->position;
				for ( k = 0; k < compressed_ds_pool[j]->num_reads; k++ )
				{
					compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->pointers_to_qual_scores[k] = compressed_ds_pool[j]->pointers_to_qual_scores[k];
					compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->pointers_to_read_names[k] = compressed_ds_pool[j]->pointers_to_read_names[k];
				}
				compressed_ds_pool_rearranged_index++;
				already_processed[j] = 1;
			}
		}
	}
	/*
	 printf ( "\n %d %d" , compressed_ds_pool_rearranged_index , compressed_ds_pool_index );
	 if ( compressed_ds_pool_index > 25000 )
	 {
	 for ( i = 0 ; i < compressed_ds_pool_index ; i++ )
	 {
	 printf ( "\n%s %d %s %d" , compressed_ds_pool[i]->icigar , compressed_ds_pool[i]->num_reads , compressed_ds_pool_rearranged[i]->icigar , compressed_ds_pool_rearranged[i]->num_reads );
	 }
	 printf ( "\n==============================================================================================================================================================================================" );
	 }
	 */
}

void compressPairedEndedAlignments(
		char *frequency_of_flags_filename,
		char *name_of_file_with_quality_scores,
		char *name_of_file_with_max_commas,
		char *input_samfilename,
		char *output_abridgefilename,
		char *unmapped_filename,
		char *genome_filename,
		char *name_of_file_with_read_names_to_short_read_names_and_NH,
		short int flag_ignore_soft_clippings,
		short int flag_ignore_mismatches,
		short int flag_ignore_unmapped_sequences,
		short int flag_ignore_all_quality_scores,
		short int run_diagnostics,
		long long int max_input_reads_in_a_single_nucl_loc,
		short int flag_ignore_quality_scores_for_matched_bases,
		//short int flag_save_exact_quality_scores,
		long long int max_number_of_alignments,
		int max_read_length,
		char *dictionary_filename,
		short int flag_ignore_alignment_scores,
		short int skip_shortening_read_names )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw_pass1;
	FILE *fhw_unmapped;
	FILE *fhw_name_of_file_with_max_commas;
	FILE *fhw_qual;
	FILE *fhw_dictionary;
	FILE *fhr_freq_samflags;
	FILE *fhr_name_of_file_with_read_names_to_short_read_names_and_NH;

	char **qual_scores;
	char **read_names;
	char **split_line; // List of strings to store each element of a single alignment
	char **split_tags; // List of strings to store tag information
	char **split_reference_info;
	char alphabets[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-='{}[]|?<>,.";
	char samformatflag_replacer_characters[] = "BEFHIJKLMOPQRSUVWXYZbdefhijklmopqrsuvwxyz";
	char *temp; //Useless
	char *line = NULL; // for reading each line
	char *line_name_of_file_with_read_names_to_short_read_names_and_NH = NULL;
	char *entry_in_output_file; //entry in output file
	char *prev_reference_name;
	char *curr_reference_name;
	char *reference_id_quick_read;
	char *samflag_quick_read;
	char *write_to_file_col1;
	char *write_to_file_col2;
	char *write_to_file_col3;
	char *encoded_string;
	char **modified_icigars;
	char *list_of_read_names;
	char *list_of_qual_scores;
	char *qual_for_writeToFile;
	char str[100];
	char read_id_string[100];
	char *line_to_be_written_to_file;
	char *convert_to_int_temp;

	size_t len = 0;
	ssize_t line_len;

	short *already_processed;
	short int NH_tag_index;

	int flag;
	int i, j, k; // Required in loops
	int number_of_tags;
	int number_of_fields; // Number of fields in each sam alignment entry
	int sam_tag_index;
	int tab_number;
	int num_items_in_alignment_pool = 0; // Items in pool
	int samflag_quick_read_index = 0;
	int compressed_ds_pool_index = 0;
	int quality_score_index = 0;
	int number_of_reference_sequences = 0;
	int reference_sequence_index = 0;
	int number_of_repetitions = 0;
	int read_length;
	int read_id[100];
	int number_of_unique_samformatflags;
	int current_length_of_line_to_be_written_to_file;
	int NH_val;

	long long int relative_position_to_previous_read_cluster;
	long long int previous_position = -1;
	long long int current_position;
	long long int number_of_records_written = 0;
	long long int number_of_records_read = 0;
	long long int num_pools_written = 0;
	long long int max_commas = 0;
	long long int curr_commas = 0;
	unsigned long long int mega_array_index;
	unsigned long long int line_number;

	struct Sam_Alignment *prev_alignment;
	struct Sam_Alignment *curr_alignment;
	struct Sam_Alignment *sam_alignment_instance_diagnostics;
	struct Sam_Alignment *temp_alignment;
	struct Sam_Alignment **alignment_pool_same_position;
	struct Compressed_DS **compressed_ds_pool;
	struct Compressed_DS **compressed_ds_pool_rearranged;
	struct Reference_Sequence_Info **reference_info;
	struct Whole_Genome_Sequence *whole_genome;
	struct Paired_Ended_Flag_to_Single_Character *samflag_dictionary;
	struct Cigar_Items *cigar_items_instance;
	cigar_items_instance = ( struct Cigar_Items* ) malloc(
			sizeof(struct Cigar_Items) * 50 );
	//printf ("\nmax_number_of_alignments %d" , max_number_of_alignments);

	//printf ("\nSize of mega array %d" , sizeof ( mega_array ));

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen( input_samfilename, "r" );
	if ( fhr == NULL )
	{
		printf( "1. Error! File %s not found", input_samfilename );
		exit( 1 );
	}
	fhr_freq_samflags = fopen( frequency_of_flags_filename, "r" );
	if ( fhr_freq_samflags == NULL )
	{
		printf( "2. Error! File %s not found", frequency_of_flags_filename );
		exit( 1 );
	}
	fhw_pass1 = fopen( output_abridgefilename, "w" );
	if ( fhw_pass1 == NULL )
	{
		printf( "3. %s File cannot be created", output_abridgefilename );
		exit( 1 );
	}
	fhw_unmapped = fopen( unmapped_filename, "w" );
	if ( fhw_unmapped == NULL )
	{
		printf( "4. %s File cannot be created", unmapped_filename );
		exit( 1 );
	}
	fhw_name_of_file_with_max_commas = fopen(
			name_of_file_with_max_commas,
			"w" );
	if ( fhw_name_of_file_with_max_commas == NULL )
	{
		printf( "5. %s File cannot be created", name_of_file_with_max_commas );
		exit( 1 );
	}
	fhw_qual = fopen( name_of_file_with_quality_scores, "w" );
	if ( fhw_qual == NULL )
	{
		printf(
				"6. %s File cannot be created",
				name_of_file_with_quality_scores );
		exit( 1 );
	}
	fhw_dictionary = fopen( dictionary_filename, "w" );
	if ( fhw_dictionary == NULL )
	{
		printf( "7. %s File cannot be created", dictionary_filename );
		exit( 1 );
	}
	fhr_name_of_file_with_read_names_to_short_read_names_and_NH = fopen(
			name_of_file_with_read_names_to_short_read_names_and_NH,
			"r" );
	if ( fhr_name_of_file_with_read_names_to_short_read_names_and_NH == NULL )
	{
		printf(
				"Error! File %s not found",
				name_of_file_with_read_names_to_short_read_names_and_NH );
		exit( 1 );
	}

	split_line = ( char** ) malloc( sizeof(char*) * ROWS );
	for ( i = 0; i < ROWS; i++ )
		split_line[i] = ( char* ) malloc( sizeof(char) * COLS );

	split_tags = ( char** ) malloc( sizeof(char*) * ROWS );
	for ( i = 0; i < ROWS; i++ )
		split_tags[i] = ( char* ) malloc( sizeof(char) * COLS );

	split_reference_info = ( char** ) malloc( sizeof(char*) * ROWS );
	for ( i = 0; i < ROWS; i++ )
		split_reference_info[i] = ( char* ) malloc( sizeof(char) * COLS );

	already_processed = ( short* ) malloc(
			sizeof(short) * max_input_reads_in_a_single_nucl_loc );
	max_input_reads_in_a_single_nucl_loc += 5;
	compressed_ds_pool = ( struct Compressed_DS** ) malloc(
			sizeof(struct Compressed_DS*)
					* max_input_reads_in_a_single_nucl_loc );
	for ( i = 0; i < max_input_reads_in_a_single_nucl_loc; i++ )
		compressed_ds_pool[i] = allocateMemoryCompressed_DS(
				max_input_reads_in_a_single_nucl_loc );

	compressed_ds_pool_rearranged = ( struct Compressed_DS** ) malloc(
			sizeof(struct Compressed_DS*)
					* max_input_reads_in_a_single_nucl_loc );
	for ( i = 0; i < max_input_reads_in_a_single_nucl_loc; i++ )
		compressed_ds_pool_rearranged[i] = allocateMemoryCompressed_DS(
				max_input_reads_in_a_single_nucl_loc );

	write_to_file_col1 = ( char* ) malloc(
			sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE );
	write_to_file_col2 = ( char* ) malloc(
			sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE );
	write_to_file_col3 = ( char* ) malloc(
			sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE );
	list_of_read_names = ( char* ) malloc(
			sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE );
	list_of_qual_scores = ( char* ) malloc(
			sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE );
	line_to_be_written_to_file = ( char* ) malloc(
			sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE );
	qual_for_writeToFile = ( char* ) malloc( sizeof(char) * MAX_SEQ_LEN );
	encoded_string = ( char* ) malloc(
			sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE );
	write_to_file_col1[0] = '\0';
	write_to_file_col2[0] = '\0';
	write_to_file_col3[0] = '\0';

	reference_id_quick_read = ( char* ) malloc( sizeof(char) * 1000 );
	samflag_quick_read = ( char* ) malloc( sizeof(char) * 1000 );
	prev_reference_name = ( char* ) malloc( sizeof(char) * 1000 );
	prev_reference_name[0] = '\0';
	curr_reference_name = ( char* ) malloc( sizeof(char) * 1000 );
	curr_reference_name[0] = '\0';

	curr_alignment = allocateMemorySam_Alignment();
	prev_alignment = allocateMemorySam_Alignment();
	temp_alignment = allocateMemorySam_Alignment();
	sam_alignment_instance_diagnostics = allocateMemorySam_Alignment();
	reference_info = ( struct Reference_Sequence_Info** ) malloc(
			sizeof(struct Reference_Sequence_Info*) * MAX_REFERENCE_SEQUENCES );
	for ( i = 0; i < MAX_REFERENCE_SEQUENCES; i++ )
		reference_info[i] = allocateMemoryReference_Sequence_Info();

	temp = ( char* ) malloc( sizeof(char) * MAX_GENERAL_LEN );
	whole_genome = ( struct Whole_Genome_Sequence* ) malloc(
			sizeof(struct Whole_Genome_Sequence) );
	qual_scores = ( char** ) malloc(
			sizeof(char*) * max_input_reads_in_a_single_nucl_loc );
	for ( i = 0; i < max_input_reads_in_a_single_nucl_loc; i++ )
		qual_scores[i] = ( char* ) malloc( sizeof(char) * MAX_SEQ_LEN );
	read_names = ( char** ) malloc(
			sizeof(char*) * max_input_reads_in_a_single_nucl_loc );
	for ( i = 0; i < max_input_reads_in_a_single_nucl_loc; i++ )
		read_names[i] = ( char* ) malloc( sizeof(char) * MAX_SEQ_LEN );
	modified_icigars = ( char** ) malloc(
			sizeof(char*) * max_input_reads_in_a_single_nucl_loc );
	for ( i = 0; i < max_input_reads_in_a_single_nucl_loc; i++ )
		modified_icigars[i] = ( char* ) malloc( sizeof(char) * MAX_SEQ_LEN );
	/********************************************************************/

	/*
	 * Write the first line in output file
	 */
	/*
	 * Write the first line in output file
	 */
	temp[0] = '\0';
	strcat( temp, "flag_ignore_mismatches:" );
	sprintf( str, "%lld", flag_ignore_mismatches );
	strcat( temp, str );
	strcat( temp, "\t" );
	strcat( temp, "flag_ignore_soft_clippings:" );
	sprintf( str, "%lld", flag_ignore_soft_clippings );
	strcat( temp, str );
	strcat( temp, "\t" );
	strcat( temp, "flag_ignore_unmapped_sequences:" );
	sprintf( str, "%lld", flag_ignore_unmapped_sequences );
	strcat( temp, str );
	strcat( temp, "\t" );
	strcat( temp, "flag_ignore_all_quality_scores:" );
	sprintf( str, "%lld", flag_ignore_all_quality_scores );
	strcat( temp, str );
	strcat( temp, "\t" );
	strcat( temp, "flag_ignore_quality_scores_for_matched_bases:" );
	sprintf( str, "%lld", flag_ignore_quality_scores_for_matched_bases );
	strcat( temp, str );
	strcat( temp, "\t" );
	strcat( temp, "flag_ignore_alignment_scores:" );
	sprintf( str, "%lld", flag_ignore_alignment_scores );
	strcat( temp, str );
	strcat( temp, "\t" );
	strcat( temp, "flag_skip_shortening_read_names:" );
	sprintf( str, "%lld", skip_shortening_read_names );
	strcat( temp, str );
	strcat( temp, "\n" );
	fprintf( fhw_pass1, "%s", temp );

	/*
	 * For diagnostics
	 */
	if ( run_diagnostics == 1 )
		readInTheEntireGenome( genome_filename, whole_genome );
	/*
	 *	Construct the mapping between samformatflag and the Single Character
	 */
	number_of_unique_samformatflags = 0;
	while ( (line_len = getline( &line, &len, fhr_freq_samflags )) != -1 )
		number_of_unique_samformatflags++;
	samflag_dictionary = allocateMemoryPaired_Ended_Flag_to_Single_Character(
			number_of_unique_samformatflags );
	rewind( fhr_freq_samflags );
	i = 0;
	j = 0;
	while ( (line_len = getline( &line, &len, fhr_freq_samflags )) != -1 )
	{
		if ( strtol( line, &temp, 10 ) == 77
				|| strtol( line, &temp, 10 ) == 141 )
		{
			number_of_unique_samformatflags--;
			continue;
		}
		samflag_dictionary->character[i + 1] = samformatflag_replacer_characters[j++];
		samflag_dictionary->samflags[i + 1] = strtol( line, &temp, 10 );
		i += 1;
	}

	for ( i = 0; i < number_of_unique_samformatflags; i++ )
	{
		//printf ("\n%c %d %c" , samflag_dictionary->direction[i] , samflag_dictionary->samflags[i] , samflag_dictionary->character[i]);
		line_to_be_written_to_file[0] = '\0';
		sprintf( temp, "%lld", samflag_dictionary->samflags[i] );
		strcat( line_to_be_written_to_file, temp );
		strcat( line_to_be_written_to_file, "\t" );
		current_length_of_line_to_be_written_to_file = strlen(
				line_to_be_written_to_file );
		line_to_be_written_to_file[current_length_of_line_to_be_written_to_file
				+ 2] = samflag_dictionary->character[i];
		line_to_be_written_to_file[current_length_of_line_to_be_written_to_file
				+ 3] = '\n';
		line_to_be_written_to_file[current_length_of_line_to_be_written_to_file
				+ 4] = '\0';
		fprintf( fhw_dictionary, "%s", line_to_be_written_to_file );
	}
	/*
	 * Read in the reference sequence information
	 */
	line_number = 0;
	while ( (line_len = getline( &line, &len, fhr )) != -1 )
	{
		line_number++;
		if ( line[0] == '@' )
		{
			if ( line[1] == 'S' && line[2] == 'Q' )
			{
				//printf("\n Reference: %s %d", line, strlen(line));
				//fflush(stdout);
				splitByDelimiter( line, '\t', split_line );
				splitByDelimiter( split_line[1], ':', split_tags );
				//printf ("\nLoading chromosome %s" , split_tags[1]);
				strcpy(
						reference_info[number_of_reference_sequences]->reference_name,
						split_tags[1] );
				strcpy(
						reference_info[number_of_reference_sequences]->line,
						line );
				number_of_reference_sequences++;

			}
		}
		else
			break;
	}
	do
	{
		//if ( line_number > 55398570 ) continue;
		/*
		 if ( line_number > 0 )
		 {
		 printf ("\nline number %llu" , line_number);
		 fflush (stdout);
		 }*/

		//printf ("\n%s" , line);
		number_of_fields = splitByDelimiter( line, '\t', split_line );
		populateSamAlignmentInstance(
				curr_alignment,
				split_line,
				number_of_fields,
				split_tags );
		strcpy( curr_reference_name, curr_alignment->reference_name );

		/***************************************************************************************
		 * Read a line from the short read names file
		 ****************************************************************************************/
		/*
		 if ( skip_shortening_read_names == 0 )
		 {
		 */
		getline(
				&line_name_of_file_with_read_names_to_short_read_names_and_NH,
				&len,
				fhr_name_of_file_with_read_names_to_short_read_names_and_NH );

		splitByDelimiter(
				line_name_of_file_with_read_names_to_short_read_names_and_NH,
				'\t',
				split_line );
		strcpy( curr_alignment->read_name, split_line[3] );

		//printf("\nRead name: %s", split_line[3]);
		NH_val = strtol( split_line[2], &convert_to_int_temp, 10 );
		NH_val = NH_val / 2;
		sprintf( str, "%lld", NH_val );
		strcpy( curr_alignment->NH, str );
		/*
		 }
		 */

		/****************************************************************************************/

		if ( curr_alignment->samflag == 77 || curr_alignment->samflag == 141 )
		{
			if ( flag_ignore_unmapped_sequences == 0 )
			{
				//Write the unmapped reads into file
				fprintf( fhw_unmapped, "%s", curr_alignment->seq );
				fprintf( fhw_unmapped, "%s", "\n" );
				for ( i = 0; curr_alignment->qual[i] != '\0'; i++ )
					curr_alignment->qual[i] -= 100;
				fprintf( fhw_unmapped, "%s", curr_alignment->qual );
				fprintf( fhw_unmapped, "%s", "\n" );
			}
			continue;
		}
		current_position = curr_alignment->start_position;
		/*
		 * Change this function
		 */

		generateIntegratedCigarPairedEnded(
				curr_alignment,
				flag_ignore_alignment_scores,
				flag_ignore_soft_clippings,
				flag_ignore_mismatches,
				flag_ignore_unmapped_sequences,
				flag_ignore_all_quality_scores,
				flag_ignore_quality_scores_for_matched_bases,
				whole_genome,
				sam_alignment_instance_diagnostics,
				number_of_records_read,
				run_diagnostics,
				samflag_dictionary,
				number_of_unique_samformatflags,
				samformatflag_replacer_characters );

		if ( strlen( prev_reference_name ) == 0 ) // 1st chromosome - initialize stuffs
		{
			//continue;
			//printf ("\n1. compressed_ds_pool_index %d" , compressed_ds_pool_index);
			//fflush (stdout);
			previous_position = current_position;
			strcpy( prev_reference_name, curr_reference_name );
			strcpy( qual_scores[quality_score_index], curr_alignment->qual );
			strcpy(
					read_names[quality_score_index],
					curr_alignment->read_name );
			strcpy(
					compressed_ds_pool[compressed_ds_pool_index]->icigar,
					curr_alignment->icigar );
			strcpy(
					compressed_ds_pool[compressed_ds_pool_index]->cigar,
					curr_alignment->cigar );
			compressed_ds_pool[compressed_ds_pool_index]->num_reads = 1;
			compressed_ds_pool[compressed_ds_pool_index]->pointers_to_qual_scores[compressed_ds_pool[compressed_ds_pool_index]->num_reads
					- 1] = qual_scores[quality_score_index];
			compressed_ds_pool[compressed_ds_pool_index]->pointers_to_read_names[compressed_ds_pool[compressed_ds_pool_index]->num_reads
					- 1] = read_names[quality_score_index];
			compressed_ds_pool[compressed_ds_pool_index]->position = curr_alignment->start_position;
//printf ("\n1. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d" , compressed_ds_pool[compressed_ds_pool_index]->num_reads , curr_alignment->reference_name , curr_alignment->start_position , compressed_ds_pool_index);
			quality_score_index++;
			compressed_ds_pool_index++;
//printf ("\n Writing Reference to file %s %d" , reference_info[reference_sequence_index]->line , reference_sequence_index);
//fflush (stdout);
			reference_sequence_index = findChromosomeIndex(
					reference_info,
					prev_reference_name,
					number_of_reference_sequences );
			fprintf(
					fhw_pass1,
					"%s",
					reference_info[reference_sequence_index]->line );
			reference_sequence_index++;
		}
		else if ( strcmp( prev_reference_name, curr_reference_name ) != 0 ) // New chromosome
		{
//continue;
//printf ("\2. ncompressed_ds_pool_index %d" , compressed_ds_pool_index);
//fflush (stdout);
			reModeliCIGARSPairedEnded(
					compressed_ds_pool,
					compressed_ds_pool_rearranged,
					already_processed,
					compressed_ds_pool_index,
					modified_icigars,
					samformatflag_replacer_characters );
			writeToFile(
					flag_ignore_quality_scores_for_matched_bases,
					flag_ignore_all_quality_scores,
					fhw_qual,
					fhw_pass1,
					compressed_ds_pool_rearranged,
					compressed_ds_pool_index,
					write_to_file_col1,
					write_to_file_col2,
					write_to_file_col3,
					encoded_string,
					&curr_commas,
					qual_scores,
					quality_score_index,
					samformatflag_replacer_characters,
					number_of_unique_samformatflags,
					samflag_dictionary,
					flag_ignore_soft_clippings,
					cigar_items_instance,
					line_to_be_written_to_file,
					list_of_read_names,
					list_of_qual_scores,
					qual_for_writeToFile );
			if ( max_commas < curr_commas )
				max_commas = curr_commas;
//printf ( "\n%lld %lld" , curr_commas , max_commas );
			compressed_ds_pool_index = 0;
			previous_position = current_position;
			strcpy( prev_reference_name, curr_reference_name );
			strcpy(
					compressed_ds_pool[compressed_ds_pool_index]->icigar,
					curr_alignment->icigar );
			strcpy(
					compressed_ds_pool[compressed_ds_pool_index]->cigar,
					curr_alignment->cigar );
			compressed_ds_pool[compressed_ds_pool_index]->num_reads = 1;
			compressed_ds_pool[compressed_ds_pool_index]->position = curr_alignment->start_position;
			strcpy( qual_scores[quality_score_index], curr_alignment->qual );
			strcpy(
					read_names[quality_score_index],
					curr_alignment->read_name );
			quality_score_index++;
//printf ("\n2. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d" , compressed_ds_pool[compressed_ds_pool_index]->num_reads , curr_alignment->reference_name , curr_alignment->start_position , compressed_ds_pool_index);
			compressed_ds_pool_index++;
//printf ("\n Writing Reference to file %s %d" , reference_info[reference_sequence_index]->line , reference_sequence_index);
//fflush (stdout);
			reference_sequence_index = findChromosomeIndex(
					reference_info,
					prev_reference_name,
					number_of_reference_sequences );
			fprintf(
					fhw_pass1,
					"%s",
					reference_info[reference_sequence_index]->line );
			reference_sequence_index++;

		}
		else // Same chromosome
		{
			if ( previous_position == current_position )
			{
				//printf ("\n3. compressed_ds_pool_index %d" , compressed_ds_pool_index);
				//fflush (stdout);
				for ( i = 0; i < compressed_ds_pool_index; i++ )
				{
					if ( strcmp(
							compressed_ds_pool[i]->icigar,
							curr_alignment->icigar ) == 0 )
					{
						compressed_ds_pool[i]->num_reads++;
						strcpy(
								qual_scores[quality_score_index],
								curr_alignment->qual );
						strcpy(
								read_names[quality_score_index],
								curr_alignment->read_name );
						compressed_ds_pool[i]->pointers_to_qual_scores[compressed_ds_pool[i]->num_reads
								- 1] = qual_scores[quality_score_index];
						compressed_ds_pool[i]->pointers_to_read_names[compressed_ds_pool[i]->num_reads
								- 1] = read_names[quality_score_index];
						quality_score_index++;
						//printf ("\n3. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d" , compressed_ds_pool[i]->num_reads , curr_alignment->reference_name , curr_alignment->start_position , i);
						break;
					}
				}
				if ( i == compressed_ds_pool_index ) // New icigar encountered
				{
					strcpy(
							compressed_ds_pool[compressed_ds_pool_index]->icigar,
							curr_alignment->icigar );
					strcpy(
							compressed_ds_pool[compressed_ds_pool_index]->cigar,
							curr_alignment->cigar );
					compressed_ds_pool[compressed_ds_pool_index]->num_reads = 1;
					compressed_ds_pool[compressed_ds_pool_index]->position = compressed_ds_pool[0]->position;
					strcpy(
							qual_scores[quality_score_index],
							curr_alignment->qual );
					strcpy(
							read_names[quality_score_index],
							curr_alignment->read_name );
					compressed_ds_pool[compressed_ds_pool_index]->pointers_to_qual_scores[compressed_ds_pool[compressed_ds_pool_index]->num_reads
							- 1] = qual_scores[quality_score_index];
					compressed_ds_pool[compressed_ds_pool_index]->pointers_to_read_names[compressed_ds_pool[compressed_ds_pool_index]->num_reads
							- 1] = read_names[quality_score_index];
					quality_score_index++;
					//printf ("\n4. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d" , compressed_ds_pool[compressed_ds_pool_index]->num_reads , curr_alignment->reference_name , curr_alignment->start_position , compressed_ds_pool_index);
					compressed_ds_pool_index++;
				}
			}
			else
			{

				//printf ("\n4. compressed_ds_pool_index %d" ,
				//		compressed_ds_pool_index);
				//fflush (stdout);
				reModeliCIGARSPairedEnded(
						compressed_ds_pool,
						compressed_ds_pool_rearranged,
						already_processed,
						compressed_ds_pool_index,
						modified_icigars,
						samformatflag_replacer_characters );

				//printf ("\nReturned from reModeliCIGARSPairedEnded");
				//fflush (stdout);

				writeToFile(
						flag_ignore_quality_scores_for_matched_bases,
						flag_ignore_all_quality_scores,
						fhw_qual,
						fhw_pass1,
						compressed_ds_pool_rearranged,
						compressed_ds_pool_index,
						write_to_file_col1,
						write_to_file_col2,
						write_to_file_col3,
						encoded_string,
						&curr_commas,
						qual_scores,
						quality_score_index,
						samformatflag_replacer_characters,
						number_of_unique_samformatflags,
						samflag_dictionary,
						flag_ignore_soft_clippings,
						cigar_items_instance,
						line_to_be_written_to_file,
						list_of_read_names,
						list_of_qual_scores,
						qual_for_writeToFile );
				//printf ("\nRead name: %s" , curr_alignment->read_name);
				//printf ("\nReturned from writeToFile");
				//fflush (stdout);

				//printf ("\n%s" , line);

				//fflush (stdout);
				//printSamAlignmentInstance (curr_alignment , 1);
				//fflush (stdout);
				//printf ("\nReturned from printSamAlignmentInstance");
				//fflush (stdout);
				//continue;
				if ( max_commas < curr_commas )
					max_commas = curr_commas;
				//printf ( "\n%lld %lld" , curr_commas , max_commas );
				compressed_ds_pool_index = 0;
				quality_score_index = 0;
				strcpy(
						qual_scores[quality_score_index],
						curr_alignment->qual );
				strcpy(
						read_names[quality_score_index],
						curr_alignment->read_name );
				strcpy(
						compressed_ds_pool[compressed_ds_pool_index]->icigar,
						curr_alignment->icigar );
				strcpy(
						compressed_ds_pool[compressed_ds_pool_index]->cigar,
						curr_alignment->cigar );
				compressed_ds_pool[compressed_ds_pool_index]->num_reads = 1;
				compressed_ds_pool[compressed_ds_pool_index]->pointers_to_qual_scores[compressed_ds_pool[compressed_ds_pool_index]->num_reads
						- 1] = qual_scores[quality_score_index];
				compressed_ds_pool[compressed_ds_pool_index]->pointers_to_read_names[compressed_ds_pool[compressed_ds_pool_index]->num_reads
						- 1] = read_names[quality_score_index];
				compressed_ds_pool[compressed_ds_pool_index]->position = curr_alignment->start_position
						- previous_position;
				quality_score_index++;
				//printf ("\n5. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d" ,
				//		compressed_ds_pool[compressed_ds_pool_index]->num_reads ,
				//		curr_alignment->reference_name ,
				//		curr_alignment->start_position ,
				//		compressed_ds_pool_index);
				compressed_ds_pool_index++;
			}
			previous_position = current_position;
		}
//printf ("\nMax_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d" , compressed_ds_pool[compressed_ds_pool_index]->num_reads , curr_alignment->reference_name , curr_alignment->start_position , compressed_ds_pool_index);
		reInitializeSamAlignmentInstance( curr_alignment );
		line_number++;
	} while ( (line_len = getline( &line, &len, fhr )) != -1 );

	reModeliCIGARSPairedEnded(
			compressed_ds_pool,
			compressed_ds_pool_rearranged,
			already_processed,
			compressed_ds_pool_index,
			modified_icigars,
			samformatflag_replacer_characters );
	writeToFile(
			flag_ignore_quality_scores_for_matched_bases,
			flag_ignore_all_quality_scores,
			fhw_qual,
			fhw_pass1,
			compressed_ds_pool_rearranged,
			compressed_ds_pool_index,
			write_to_file_col1,
			write_to_file_col2,
			write_to_file_col3,
			encoded_string,
			&curr_commas,
			qual_scores,
			quality_score_index,
			samformatflag_replacer_characters,
			number_of_unique_samformatflags,
			samflag_dictionary,
			flag_ignore_soft_clippings,
			cigar_items_instance,
			line_to_be_written_to_file,
			list_of_read_names,
			list_of_qual_scores,
			qual_for_writeToFile );
	if ( max_commas < curr_commas )
		max_commas = curr_commas;
	sprintf( temp, "%lld", max_commas );
	strcat( temp, "\n" );
	fprintf( fhw_name_of_file_with_max_commas, "%s", temp );

	fclose( fhr );
	fclose( fhw_pass1 );
	fclose( fhw_unmapped );
	fclose( fhw_name_of_file_with_max_commas );
	fclose( fhw_dictionary );
}

int main( int argc, char *argv[] )
{
	/********************************************************************
	 * Named CLI
	 ********************************************************************/
	struct arguments arguments;

	// Parse our arguments; every option seen by parse_opt will be reflected in arguments.
	// Default values.
	arguments.input_sam_filename = ""; // Empty string - only contains null character
	arguments.output_abridge_filename = "";
	arguments.genome_filename = "";
	arguments.unmapped_filename = "";
	arguments.name_of_file_with_max_commas = "";
	arguments.name_of_file_with_quality_scores = "";
	arguments.name_of_file_with_read_names_to_short_read_names_and_NH = "";
	arguments.frequency_of_flags_filename = "";
	arguments.name_of_file_dictionary = "";

	arguments.flag_ignore_soft_clippings = 0;
	arguments.flag_ignore_mismatches = 0;
	arguments.flag_ignore_all_quality_scores = 0;
	arguments.flag_ignore_unmapped_sequences = 0;
	arguments.flag_ignore_quality_scores_for_matched_bases = 0;
	arguments.flag_ignore_alignment_scores = 0;
	arguments.run_diagnostics = 0;
	arguments.max_input_reads_in_a_single_nucl_loc = 0;
	arguments.max_number_of_alignments = 0;
	arguments.max_read_length = 0;
	arguments.skip_shortening_read_names = 0;

	argp_parse( &argp, argc, argv, 0, 0, &arguments );

	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char input_samfilename[FILENAME_LENGTH];
	char output_abridgefilename[FILENAME_LENGTH];
	char genome_filename[FILENAME_LENGTH];
	char unmapped_filename[FILENAME_LENGTH];
	char name_of_file_with_max_commas[FILENAME_LENGTH];
	char name_of_file_with_quality_scores[FILENAME_LENGTH];
	char frequency_of_flags_filename[FILENAME_LENGTH];
	char dictionary_filename[FILENAME_LENGTH];
	char name_of_file_with_read_names_to_short_read_names_and_NH[FILENAME_LENGTH];
	char *temp; //Required for strtoi

	short int flag_ignore_soft_clippings;
	short int flag_ignore_mismatches;
	short int flag_ignore_all_quality_scores;
	short int flag_ignore_unmapped_sequences;
	short int flag_ignore_quality_scores_for_matched_bases;
	short int run_diagnostics;
	short int flag_ignore_alignment_scores;
	short int skip_shortening_read_names;

	long long int max_input_reads_in_a_single_nucl_loc;
	long long int max_number_of_alignments;

	int i;
	int max_read_length;

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy( genome_filename, arguments.genome_filename );
	strcpy( input_samfilename, arguments.input_sam_filename );
	strcpy( output_abridgefilename, arguments.output_abridge_filename );
	strcpy( unmapped_filename, arguments.unmapped_filename );
	strcpy(
			name_of_file_with_max_commas,
			arguments.name_of_file_with_max_commas );
	strcpy(
			name_of_file_with_quality_scores,
			arguments.name_of_file_with_quality_scores );
	strcpy(
			frequency_of_flags_filename,
			arguments.frequency_of_flags_filename );
	strcpy( dictionary_filename, arguments.name_of_file_dictionary );
	strcpy(
			name_of_file_with_read_names_to_short_read_names_and_NH,
			arguments.name_of_file_with_read_names_to_short_read_names_and_NH );

	flag_ignore_soft_clippings = arguments.flag_ignore_soft_clippings;
	flag_ignore_mismatches = arguments.flag_ignore_mismatches;
	flag_ignore_all_quality_scores = arguments.flag_ignore_all_quality_scores;
	flag_ignore_unmapped_sequences = arguments.flag_ignore_unmapped_sequences;
	flag_ignore_quality_scores_for_matched_bases = arguments.flag_ignore_quality_scores_for_matched_bases;
	flag_ignore_alignment_scores = arguments.flag_ignore_alignment_scores;

	run_diagnostics = arguments.run_diagnostics;
	max_input_reads_in_a_single_nucl_loc = arguments.max_input_reads_in_a_single_nucl_loc;
	max_number_of_alignments = arguments.max_number_of_alignments;
	max_read_length = arguments.max_read_length;
	skip_shortening_read_names = arguments.skip_shortening_read_names;

	/********************************************************************/

	compressPairedEndedAlignments(
			frequency_of_flags_filename,
			name_of_file_with_quality_scores,
			name_of_file_with_max_commas,
			input_samfilename,
			output_abridgefilename,
			unmapped_filename,
			genome_filename,
			name_of_file_with_read_names_to_short_read_names_and_NH,
			flag_ignore_soft_clippings,
			flag_ignore_mismatches,
			flag_ignore_unmapped_sequences,
			flag_ignore_all_quality_scores,
			run_diagnostics,
			max_input_reads_in_a_single_nucl_loc,
			flag_ignore_quality_scores_for_matched_bases,
			max_number_of_alignments,
			max_read_length,
			dictionary_filename,
			flag_ignore_alignment_scores,
			skip_shortening_read_names );
	return 0;
}

