# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
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
{ "input_sam_filename", 'i', "SAM_FILENAME", 0, "Enter the name of the SAM file to be compressed", 0 },
{ "output_abridge_filename", 'o', "TEXT_FILENAME", 0, "Enter the name of the compressed file (please note that this is not the final compressed file)", 0 },
{ "genome_filename", 'g', "GENOME_FILENAME", 0, "Enter the name of the genome file in fasta format", 0 },
{ "unmapped_filename", 'u', "UNMAPPED_READS_FILENAME", 0, "Enter the name of the file where the unmapped reads will be stored", 0 },
{ "name_of_file_with_max_commas", 'c', "MAX_COMMAS_FILENAME", 0, "Enter the name of the file that contains the value of maximum number of commas", 0 },
{ "name_of_file_with_quality_scores", 'q', "QUALITY_SCORES_FILENAME", 0, "Enter the name of the file where the quality scores will be stored. This file will be compressed later", 0 },
{ "name_of_file_with_read_names_to_short_read_names_and_NH", 'r', "SHORT_NAMES_NH_FILENAME", 0, "Enter the name of the file that contains the mapping between the long name to the short name and the NH values", 0 },

{ "flag_ignore_soft_clippings", 's', 0, 0, "Set this flag to ignore soft clippings", 0 },
{ "flag_ignore_mismatches", 'm', 0, 0, "Set this flag to ignore mismatches", 0 },
{ "flag_ignore_all_quality_scores", 'p', 0, 0, "Set this flag to ignore quality scores for mismatched bases and soft clips", 0 },
{ "flag_ignore_unmapped_sequences", 'e', 0, 0, "Set this flag to ignore unmapped sequences along with their quality scores", 0 },
{ "flag_ignore_quality_scores_for_matched_bases", 'b', 0, 0, "Set this flag to ignore quality scores for nucleotide bases that match to the provided reference", 0 },
{ "flag_ignore_alignment_scores", 'a', 0, 0, "Set this flag to ignore the alignment scores (Column 5 of SAM file)", 0 },
{ "skip_shortening_read_names", 'f', 0, 0, "Set this flag to skip shortening read names", 0 },
{ "run_diagnostics", 'd', 0, 0, "Set this flag to run diagnostics and print out a verbose report", 0 },

{ "max_input_reads_in_a_single_nucl_loc", 'n', "MAX_READS_IN_ONE_NUCL", 0, "Enter the value of the maximum number of input reads mapped to a single nucleotide", 0 },
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
	int run_diagnostics;
	int flag_ignore_quality_scores_for_matched_bases;
	int flag_ignore_alignment_scores;
	int skip_shortening_read_names;
	unsigned long long int max_input_reads_in_a_single_nucl_loc;
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

int total_mapped_reads;

char findMatchCharacterIcigar( char *icigar )
{
	int i;

	for ( i = 0; icigar[i] != '\0'; i++ )
	{
		switch ( icigar[i] )
		{
			case 'B':
				return 'B';
			case 'E':
				return 'E';
			case 'F':
				return 'F';
			case 'H':
				return 'H';
			case 'J':
				return 'J';
			case 'K':
				return 'K';
			case 'L':
				return 'L';
			case 'O':
				return 'O';
			case 'P':
				return 'P';
			case 'Q':
				return 'Q';
			case 'R':
				return 'R';
			case 'U':
				return 'U';
		}
	}
	return ' ';
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
		short int flag_ignore_soft_clippings,
		struct Cigar_Items *cigar_items_instance,
		char *line_to_be_written_to_file,
		char *list_of_read_names,
		char *list_of_qual_scores,
		char *qual )
{
	//printf ("\nInside writeToFile\n");
	//fflush (stdout);

	int i, j, k, l, m;
	char str[1000];
	//char qual[MAX_SEQ_LEN];
	//char line_to_be_written_to_file[MAX_LINE_TO_BE_WRITTEN_TO_FILE];
	//char list_of_read_names[MAX_LINE_TO_BE_WRITTEN_TO_FILE];
	//char list_of_qual_scores[MAX_LINE_TO_BE_WRITTEN_TO_FILE];
	int num_of_types;
	int qual_score_length, cigar_length;

	line_to_be_written_to_file[0] = '\0';
	list_of_read_names[0] = '\0';
	list_of_qual_scores[0] = '\0';

	//printf ("\ncompressed_ds_pool_total %d" , compressed_ds_pool_total);
	//fflush (stdout);

	for ( i = 0; i < compressed_ds_pool_total; i++ )
	{
		if ( i == 0 )
		{
			if ( compressed_ds_pool[i]->position != 1 )
			{
				sprintf( str, "%lld", compressed_ds_pool[i]->position );
				//printf ("\nstr=%s" , str);
				//fflush (stdout);
				strcat( line_to_be_written_to_file, str );
				strcat( line_to_be_written_to_file, "\t" );
			}
			else
				str[0] = '\0'; // empty string
		}

		strcat( line_to_be_written_to_file, compressed_ds_pool[i]->icigar );
		strcat( line_to_be_written_to_file, "-" );
		// Write the number of reads only if it is greater than one - will save some space
		if ( compressed_ds_pool[i]->num_reads > 1 )
		{
			sprintf( str, "%ld", compressed_ds_pool[i]->num_reads );
			strcat( line_to_be_written_to_file, str );
		}
		if ( i != compressed_ds_pool_total - 1 )
			strcat( line_to_be_written_to_file, "," );

		for ( j = 0; j < compressed_ds_pool[i]->num_reads; j++ )
		{
			if ( compressed_ds_pool[i]->pointers_to_read_names[j][0] != ' '
					&& compressed_ds_pool[i]->pointers_to_read_names[j][1]
							!= '\0' )
			{
				printf( "\nInside here" );
				//strcat(list_of_read_names, "brdg_");
				strcat(
						list_of_read_names,
						compressed_ds_pool[i]->pointers_to_read_names[j] );
				if ( i != compressed_ds_pool_total - 1
						|| j != compressed_ds_pool[i]->num_reads - 1 )
					strcat( list_of_read_names, "," );
			}
		}
		if ( flag_ignore_all_quality_scores == 0
				&& flag_ignore_quality_scores_for_matched_bases == 0 ) //Write out the entire quality score
		{
			for ( j = 0; j < compressed_ds_pool[i]->num_reads; j++ )
			{
				qual_score_length = 0;
				cigar_length = 0;
				switch ( findMatchCharacterIcigar(
						compressed_ds_pool[i]->icigar ) )
				{
					case 'B':
					case 'F':
					case 'J':
					case 'L':
					case 'P':
					case 'R':
						for ( k = 0;
								compressed_ds_pool[i]->pointers_to_qual_scores[j][k]
										!= '\0'; k++ )
							qual[k] = compressed_ds_pool[i]->pointers_to_qual_scores[j][k]
									- 90;
						qual[k] = '\0';
						break;
					case 'E':
					case 'H':
					case 'K':
					case 'O':
					case 'Q':
					case 'U':
						/*
						 //Reverse the quality scores if alignment was to the other strand
						 for ( k = strlen (compressed_ds_pool[i]->pointers_to_qual_scores[j]) - 1 ;
						 k >= 0 ; k-- )
						 {
						 qual[strlen (compressed_ds_pool[i]->pointers_to_qual_scores[j]) - 1 - k] = compressed_ds_pool[i]->pointers_to_qual_scores[j][k] - 90;
						 //printf ("\n Read name check %s" , compressed_ds_pool[i]->pointers_to_read_names[j]);
						 }
						 qual[strlen (compressed_ds_pool[i]->pointers_to_qual_scores[j])] = '\0';
						 */

						for ( k = 0;
								compressed_ds_pool[i]->pointers_to_qual_scores[j][k]
										!= '\0'; k++ )
							qual[k] = compressed_ds_pool[i]->pointers_to_qual_scores[j][k]
									- 90;
						qual[k] = '\0';
						break;
				}
				// Writing dummy lines for fclqc
				fprintf( fhw_qual, "%s", "\n" );
				fprintf( fhw_qual, "%s", "\n" );
				fprintf( fhw_qual, "%s", "\n" );

				fprintf( fhw_qual, "%s", qual );
				fprintf( fhw_qual, "%s", "\n" );
				qual_score_length = strlen( qual );

				splitCigar(
						compressed_ds_pool[i]->cigar,
						&num_of_types,
						cigar_items_instance );
				for ( m = 0; m < num_of_types; m++ )
					if ( cigar_items_instance[m].def != 'D'
							&& cigar_items_instance[m].def != 'N' )
						cigar_length += cigar_items_instance[m].len;

				if ( qual_score_length != cigar_length )
				{
					printf( "\nUNEQUAL LENGTHS" );
					printf( "\ncigar %s", compressed_ds_pool[i]->cigar );
					printf( "\nqual score %s", qual );
					printf(
							"\ncigar length %d qual score length %d",
							cigar_length,
							qual_score_length );
				}
				/*
				 if (flag_save_exact_quality_scores == 0)
				 {
				 fprintf(fhw_qual, "%s", "\t");
				 if (flag_ignore_soft_clippings == 1)
				 {
				 splitCigar(
				 compressed_ds_pool[i]->cigar,
				 &num_of_types,
				 cigar_items_instance);
				 if (cigar_items_instance[0].def == 'S'
				 && compressed_ds_pool[i]->icigar[1] != '\0') // Left soft clip exists
				 {
				 sprintf(str, "%ld", cigar_items_instance[0].len);
				 fprintf(fhw_qual, "%s", str);
				 fprintf(fhw_qual, "%s", "S");
				 }
				 }
				 if (compressed_ds_pool[i]->icigar[1] != '\0')
				 {
				 for (k = 0;
				 compressed_ds_pool[i]->icigar[k + 1] != '~'
				 && compressed_ds_pool[i]->icigar[k + 1]
				 != '\0'; k++)
				 fputc(compressed_ds_pool[i]->icigar[k], fhw_qual);
				 }
				 else
				 {

				 //Find the icigar before the ith one which is not compressed

				 l = i - 1;
				 while (l >= 0)
				 {
				 if (compressed_ds_pool[l]->icigar[1] != '\0')
				 break;
				 else
				 l--;
				 }

				 for (k = 0;
				 compressed_ds_pool[l]->icigar[k + 1] != '~'
				 && compressed_ds_pool[l]->icigar[k + 1]
				 != '\0'; k++)
				 fputc(compressed_ds_pool[l]->icigar[k], fhw_qual);
				 }
				 if (flag_ignore_soft_clippings == 1)
				 {
				 //splitCigar (compressed_ds_pool[i]->cigar , &num_of_types , cigar_items_instance);
				 if (cigar_items_instance[num_of_types - 1].def == 'S'
				 && compressed_ds_pool[i]->icigar[1] != '\0') // Right soft clip exists
				 {
				 sprintf(
				 str,
				 "%ld",
				 cigar_items_instance[num_of_types - 1].len);
				 fprintf(fhw_qual, "%s", str);
				 fprintf(fhw_qual, "%s", "S");
				 }
				 }
				 fprintf(fhw_qual, "%s", "\t");
				 // Print whether the read was mapped in forward or the reverse direction
				 switch (findMatchCharacterIcigar(
				 compressed_ds_pool[i]->icigar))
				 {
				 case 'B':
				 case 'F':
				 case 'J':
				 case 'L':
				 case 'P':
				 case 'R':
				 fprintf(fhw_qual, "%s", "1");
				 break;
				 case 'E':
				 case 'H':
				 case 'K':
				 case 'O':
				 case 'Q':
				 case 'U':
				 fprintf(fhw_qual, "%s", "2");
				 break;
				 }
				 }
				 */
			}
		}
	}

	//  Check if the last element of list_of_read_names is a comma
	if ( list_of_read_names[strlen( list_of_read_names ) - 1] == ',' )
		list_of_read_names[strlen( list_of_read_names ) - 1] = '\0';
	if ( list_of_read_names[0] != '\0' )
		strcat( line_to_be_written_to_file, "\t" );
	strcat( line_to_be_written_to_file, list_of_read_names );
	strcat( line_to_be_written_to_file, "\n" );
	fprintf( fhw_pass1, "%s", line_to_be_written_to_file );
	*count = compressed_ds_pool_total;

//strcat(line_to_be_written_to_file , "\n");
//fprintf (fhw_pass1 , "%s" , line_to_be_written_to_file);
//*count = compressed_ds_pool_total;
}

void prepareIcigarForComparison( char *icigar1, char *icigar )
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
		switch ( icigar[i] )
		{
			case 'B':
			case 'E':
			case 'F':
			case 'H':
			case 'J':
			case 'K':
			case 'L':
			case 'O':
			case 'P':
			case 'Q':
			case 'R':
			case 'U':
				icigar1[i] = 'M';
				break;
			default:
				icigar1[i] = icigar[i];
				break;
		}
	}
	icigar1[i] = '\0';
}

void reModeliCIGARSSingleEnded(
		struct Compressed_DS **compressed_ds_pool,
		struct Compressed_DS **compressed_ds_pool_rearranged,
		short *already_processed,
		int compressed_ds_pool_index,
		char **modified_icigars,
		struct Cigar_Items *cigar_items_instance )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i, j, k, m;
	int compressed_ds_pool_rearranged_index = 0;
	int change_flag = 0;
	int num_of_types;

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
				compressed_ds_pool[i]->icigar );

	for ( i = 0; i < compressed_ds_pool_index; i++ )
	{
		if ( already_processed[i] == 1 )
			continue;
		already_processed[i] = 1;
		//prepareIcigarForComparison ( icigar1 , compressed_ds_pool[i]->icigar );
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
						compressed_ds_pool[j]->icigar );
				compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->icigar[1] = '\0';
				strcpy(
						compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->cigar,
						compressed_ds_pool[j]->cigar );
				compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->num_reads = compressed_ds_pool[j]->num_reads;
				compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->position = compressed_ds_pool[j]->position;
				for ( k = 0; k < compressed_ds_pool[j]->num_reads; k++ )
				{
					//printf ("\nEntering here");
					compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->pointers_to_qual_scores[k] = compressed_ds_pool[j]->pointers_to_qual_scores[k];
					compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->pointers_to_read_names[k] = compressed_ds_pool[j]->pointers_to_read_names[k];

					//printf ("\n%s \n%s" , compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->pointers_to_qual_scores[k] , compressed_ds_pool[j]->pointers_to_qual_scores[k]);
				}
				compressed_ds_pool_rearranged_index++;
				already_processed[j] = 1;
			}
		}
	}

	for ( i = 0; i < compressed_ds_pool_index; i++ )
	{
		int cigar_length = 0;
		splitCigar(
				compressed_ds_pool_rearranged[i]->cigar,
				&num_of_types,
				cigar_items_instance );
		for ( m = 0; m < num_of_types; m++ )
			cigar_length += cigar_items_instance[m].len;
		for ( j = 0; j < compressed_ds_pool_rearranged[i]->num_reads; j++ )
		{
			if ( strlen(
					compressed_ds_pool_rearranged[i]->pointers_to_qual_scores[j] )
					!= cigar_length )
			{
				change_flag = 1;
				/*
				 printf ("\nUnequal Lengths inside reModeliCIGARSSingleEnded");
				 printf ("\ncigar %s" , compressed_ds_pool_rearranged[i]->cigar);
				 printf ("\nicigar %s" , compressed_ds_pool_rearranged[i]->icigar);
				 printf ("\nQual length %d" , strlen (compressed_ds_pool_rearranged[i]->pointers_to_qual_scores[j]));
				 printf ("\n");
				 for ( k = 0 ;
				 compressed_ds_pool_rearranged[i]->pointers_to_qual_scores[j][k] != '\0' ;
				 k++ )
				 printf ("%c" , compressed_ds_pool_rearranged[i]->pointers_to_qual_scores[j][k] - 90);
				 printf ("\n========================================================================================================================================================================");
				 printf ("\n========================================================================================================================================================================");
				 */
			}
		}
	}
	/*
	 if ( change_flag == 1 )
	 {
	 printf ("\nPrinting everything");
	 for ( i = 0 ; i < compressed_ds_pool_index ; i++ )
	 {
	 printf ("\ni=%d" , i);
	 printf ("\ncigar = %s icigar = %s" , compressed_ds_pool[i]->cigar , compressed_ds_pool[i]->icigar);
	 for ( j = 0 ; j < compressed_ds_pool[i]->num_reads ; j++ )
	 {
	 printf ("\n%d Qual length = %d" , j , strlen (compressed_ds_pool[i]->pointers_to_qual_scores[j]));
	 printf ("\n");
	 for ( k = 0 ;
	 compressed_ds_pool[i]->pointers_to_qual_scores[j][k] != '\0' ;
	 k++ )

	 printf ("%c" , compressed_ds_pool[i]->pointers_to_qual_scores[j][k] - 90);
	 printf ("\n");
	 }
	 }
	 printf ("\n********************************************************************************************************************************************************************************************");
	 for ( i = 0 ; i < compressed_ds_pool_index ; i++ )
	 {
	 printf ("\ni=%d" , i);
	 printf ("\ncigar = %s icigar = %s" , compressed_ds_pool_rearranged[i]->cigar , compressed_ds_pool_rearranged[i]->icigar);
	 for ( j = 0 ; j < compressed_ds_pool_rearranged[i]->num_reads ;
	 j++ )
	 {
	 printf ("\n%d Qual length = %d" , j , strlen (compressed_ds_pool_rearranged[i]->pointers_to_qual_scores[j]));
	 printf ("\n");
	 for ( k = 0 ;
	 compressed_ds_pool_rearranged[i]->pointers_to_qual_scores[j][k] != '\0' ;
	 k++ )

	 printf ("%c" , compressed_ds_pool_rearranged[i]->pointers_to_qual_scores[j][k] - 90);
	 printf ("\n");
	 }
	 }

	 printf ("\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
	 printf ("\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
	 }

	 */
	//printf ("\nReturning from reModeliCIGARSSingleEnded");
	//fflush (stdout);
}

void readAlignmentsAndCompress(
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
	FILE *fhr_name_of_file_with_read_names_to_short_read_names_and_NH;

	char **qual_scores;
	char **read_names;
	char **split_line; // List of strings to store each element of a single alignment
	char **split_tags; // List of strings to store tag information
	char **split_reference_info;
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
	char *line_to_be_written_to_file;
	char *list_of_read_names;
	char *list_of_qual_scores;
	char *qual_for_writeToFile;
	char str[100];

	size_t len = 0;
	ssize_t line_len;

	short *already_processed;

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

	long long int relative_position_to_previous_read_cluster;
	long long int previous_position = -1;
	long long int current_position;
	long long int number_of_records_written = 0;
	long long int number_of_records_read = 0;
	long long int num_pools_written = 0;
	long long int max_commas = 0;
	long long int curr_commas = 0;

	struct Sam_Alignment *prev_alignment;
	struct Sam_Alignment *curr_alignment;
	struct Sam_Alignment *sam_alignment_instance_diagnostics;
	struct Sam_Alignment *temp_alignment;
	struct Sam_Alignment **alignment_pool_same_position;
	struct Compressed_DS **compressed_ds_pool;
	struct Compressed_DS **compressed_ds_pool_rearranged;
	struct Reference_Sequence_Info **reference_info;
	struct Whole_Genome_Sequence *whole_genome;
	struct Cigar_Items *cigar_items_instance;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen( input_samfilename, "r" );
	if ( fhr == NULL )
	{
		printf( "Error! File %s not found", input_samfilename );
		exit( 1 );
	}
	fhw_pass1 = fopen( output_abridgefilename, "w" );
	if ( fhw_pass1 == NULL )
	{
		printf( "%s File cannot be created", output_abridgefilename );
		exit( 1 );
	}
	fhw_unmapped = fopen( unmapped_filename, "w" );
	if ( fhw_unmapped == NULL )
	{
		printf( "%s File cannot be created", unmapped_filename );
		exit( 1 );
	}
	fhw_name_of_file_with_max_commas = fopen(
			name_of_file_with_max_commas,
			"w" );
	if ( fhw_name_of_file_with_max_commas == NULL )
	{
		printf( "%s File cannot be created", name_of_file_with_max_commas );
		exit( 1 );
	}
	fhw_qual = fopen( name_of_file_with_quality_scores, "w" );
	if ( fhw_qual == NULL )
	{
		printf( "%s File cannot be created", name_of_file_with_quality_scores );
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
	line_to_be_written_to_file = ( char* ) malloc(
			sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE );
	list_of_read_names = ( char* ) malloc(
			sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE );
	list_of_qual_scores = ( char* ) malloc(
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
	cigar_items_instance = ( struct Cigar_Items* ) malloc(
			sizeof(struct Cigar_Items) * 50 );
	sam_alignment_instance_diagnostics = allocateMemorySam_Alignment();
	reference_info = ( struct Reference_Sequence_Info** ) malloc(
			sizeof(struct Reference_Sequence_Info*) * MAX_REFERENCE_SEQUENCES );
	for ( i = 0; i < MAX_REFERENCE_SEQUENCES; i++ )
		reference_info[i] = allocateMemoryReference_Sequence_Info();

	read_names = ( char** ) malloc(
			sizeof(char*) * max_input_reads_in_a_single_nucl_loc );
	for ( i = 0; i < max_input_reads_in_a_single_nucl_loc; i++ )
		read_names[i] = ( char* ) malloc( sizeof(char) * MAX_SEQ_LEN );
	temp = ( char* ) malloc( sizeof(char) * MAX_GENERAL_LEN );
	whole_genome = ( struct Whole_Genome_Sequence* ) malloc(
			sizeof(struct Whole_Genome_Sequence) );
	qual_scores = ( char** ) malloc(
			sizeof(char*) * max_input_reads_in_a_single_nucl_loc );
	for ( i = 0; i < max_input_reads_in_a_single_nucl_loc; i++ )
		qual_scores[i] = ( char* ) malloc( sizeof(char) * MAX_SEQ_LEN );
	modified_icigars = ( char** ) malloc(
			sizeof(char*) * max_input_reads_in_a_single_nucl_loc );
	for ( i = 0; i < max_input_reads_in_a_single_nucl_loc; i++ )
		modified_icigars[i] = ( char* ) malloc( sizeof(char) * MAX_SEQ_LEN );
	/********************************************************************/

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
	 * Read in the reference sequence information
	 */
	while ( (line_len = getline( &line, &len, fhr )) != -1 )
	{
		if ( line[0] == '@' )
		{
			if ( line[1] == 'S' && line[2] == 'Q' )
			{
				//printf("\n Reference: %s %d", line, strlen(line));
				//fflush(stdout);
				strcpy(
						reference_info[number_of_reference_sequences]->line,
						line );
				splitByDelimiter( line, '\t', split_line );
				splitByDelimiter( split_line[1], ':', split_tags );
				//printf ("\nLoading chromosome %s" , split_tags[1]);
				strcpy(
						reference_info[number_of_reference_sequences]->reference_name,
						split_tags[1] );
				number_of_reference_sequences++;
			}
		}
		else
			break;
	}
	//return;
	do
	{
		number_of_records_read += 1;
		/*if (number_of_records_read % 1000000 == 0)
		 {
		 printf("\nNumber of records read: %lld Million %s", number_of_records_read / 1000000, output_abridgefilename);
		 fflush(stdout);
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
		if ( skip_shortening_read_names == 0 )
		{
			getline(
					&line_name_of_file_with_read_names_to_short_read_names_and_NH,
					&len,
					fhr_name_of_file_with_read_names_to_short_read_names_and_NH );

			splitByDelimiter(
					line_name_of_file_with_read_names_to_short_read_names_and_NH,
					'\t',
					split_line );
			strcpy( curr_alignment->read_name, split_line[3] );
		}
		//printf( "\nRead Name: %s", curr_alignment->read_name );

		for ( i = 11; i < number_of_fields; i++ )
		{
			if ( strcmp( curr_alignment->tags[i - 11].name, "NH" ) == 0 )
			{
				strcpy( curr_alignment->tags[i - 11].type, "i" );
				strcpy( curr_alignment->tags[i - 11].val, split_line[2] );
				break;
			}
		}
		if ( i == number_of_fields )
		{
			strcpy(
					curr_alignment->tags[curr_alignment->number_of_tag_items].name,
					"NH" );
			strcpy(
					curr_alignment->tags[curr_alignment->number_of_tag_items].type,
					"i" );
			strcpy(
					curr_alignment->tags[curr_alignment->number_of_tag_items].val,
					split_line[2] );
			curr_alignment->number_of_tag_items++;

		}

		/****************************************************************************************/

		if ( curr_alignment->samflag == 4 )
		{
			if ( flag_ignore_unmapped_sequences == 0 )
			{
				//Write the unmapped reads into file
				fprintf( fhw_unmapped, "%s", curr_alignment->seq );
				fprintf( fhw_unmapped, "%s", "\n" );
				for ( i = 0; curr_alignment->qual[i] != '\0'; i++ )
					curr_alignment->qual[i] -= 90;
				fprintf( fhw_unmapped, "%s", curr_alignment->qual );
				fprintf( fhw_unmapped, "%s", "\n" );
			}
			continue;
		}

		current_position = curr_alignment->start_position;
		//printSamAlignmentInstance(curr_alignment,0);
		generateIntegratedCigarSingleEnded(
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
				run_diagnostics );
		//printf ("\n Position:%lld iCIGAR: %s" , curr_alignment->start_position , curr_alignment->icigar);
		if ( strlen( prev_reference_name ) == 0 ) // 1st chromosome - initialize stuffs
		{
			//printf("\n1. compressed_ds_pool_index %d", compressed_ds_pool_index);
			//fflush(stdout);
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
			//printf("\n1. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d", compressed_ds_pool[compressed_ds_pool_index]->num_reads, curr_alignment->reference_name, curr_alignment->start_position, compressed_ds_pool_index);
			quality_score_index++;
			compressed_ds_pool_index++;
			//printf("\n Writing Reference to file %s %d", reference_info[reference_sequence_index]->line, reference_sequence_index);
			//fflush(stdout);
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
			//printf("\2. ncompressed_ds_pool_index %d", compressed_ds_pool_index);
			//fflush(stdout);
			reModeliCIGARSSingleEnded(
					compressed_ds_pool,
					compressed_ds_pool_rearranged,
					already_processed,
					compressed_ds_pool_index,
					modified_icigars,
					cigar_items_instance );
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
			quality_score_index = 0;
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
			//printf("\n2. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d", compressed_ds_pool[compressed_ds_pool_index]->num_reads, curr_alignment->reference_name, curr_alignment->start_position, compressed_ds_pool_index);
			compressed_ds_pool_index++;
			//printf("\n Writing Reference to file %s %d", reference_info[reference_sequence_index]->line, reference_sequence_index);
			//fflush(stdout);
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
				//printf("\n3. compressed_ds_pool_index %d", compressed_ds_pool_index);
				//fflush(stdout);

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
						//printf("\n3. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d", compressed_ds_pool[i]->num_reads, curr_alignment->reference_name, curr_alignment->start_position, i);
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
					//printf("\n4. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d", compressed_ds_pool[compressed_ds_pool_index]->num_reads, curr_alignment->reference_name, curr_alignment->start_position, compressed_ds_pool_index);
					compressed_ds_pool_index++;
				}
			}
			else
			{

				//printf ("\n4. compressed_ds_pool_index %d" , compressed_ds_pool_index);
				//fflush (stdout);
				reModeliCIGARSSingleEnded(
						compressed_ds_pool,
						compressed_ds_pool_rearranged,
						already_processed,
						compressed_ds_pool_index,
						modified_icigars,
						cigar_items_instance );
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
						flag_ignore_soft_clippings,
						cigar_items_instance,
						line_to_be_written_to_file,
						list_of_read_names,
						list_of_qual_scores,
						qual_for_writeToFile );

				if ( max_commas < curr_commas )
					max_commas = curr_commas;
				//printf ("\n%lld %lld" , curr_commas , max_commas);
				//fflush (stdout);
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
				//printf ("\n5. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d" , compressed_ds_pool[compressed_ds_pool_index]->num_reads , curr_alignment->reference_name , curr_alignment->start_position , compressed_ds_pool_index);
				//fflush (stdout);
				compressed_ds_pool_index++;
			}
			previous_position = current_position;
		}
		//printf("\nMax_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d", compressed_ds_pool[compressed_ds_pool_index]->num_reads, curr_alignment->reference_name, curr_alignment->start_position, compressed_ds_pool_index);
		reInitializeSamAlignmentInstance( curr_alignment );
	} while ( (line_len = getline( &line, &len, fhr )) != -1 );

	/*
	 *Write final data to file
	 */
	reModeliCIGARSSingleEnded(
			compressed_ds_pool,
			compressed_ds_pool_rearranged,
			already_processed,
			compressed_ds_pool_index,
			modified_icigars,
			cigar_items_instance );
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
	arguments.flag_ignore_soft_clippings = 0;
	arguments.flag_ignore_mismatches = 0;
	arguments.flag_ignore_all_quality_scores = 0;
	arguments.flag_ignore_unmapped_sequences = 0;
	arguments.flag_ignore_quality_scores_for_matched_bases = 0;
	arguments.flag_ignore_alignment_scores = 0;
	arguments.run_diagnostics = 0;
	arguments.max_input_reads_in_a_single_nucl_loc = 0;
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
	char name_of_file_with_read_names_to_short_read_names_and_NH[FILENAME_LENGTH];
	char *temp; //Required for strtoi

	short int flag_ignore_soft_clippings;
	short int flag_ignore_mismatches;
	short int flag_ignore_all_quality_scores;
	short int flag_ignore_unmapped_sequences;
	short int run_diagnostics;
	short int flag_ignore_quality_scores_for_matched_bases;
	short int flag_ignore_alignment_scores;
	short int skip_shortening_read_names;

	long long int max_input_reads_in_a_single_nucl_loc;
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
			name_of_file_with_read_names_to_short_read_names_and_NH,
			arguments.name_of_file_with_read_names_to_short_read_names_and_NH );

	flag_ignore_alignment_scores = arguments.flag_ignore_alignment_scores; // Ignore the column 5 of SAM alignment file which is often set to 255 and also the AS tag if one is provided
	flag_ignore_soft_clippings = arguments.flag_ignore_soft_clippings;
	flag_ignore_mismatches = arguments.flag_ignore_mismatches;
	flag_ignore_all_quality_scores = arguments.flag_ignore_all_quality_scores;
	flag_ignore_unmapped_sequences = arguments.flag_ignore_unmapped_sequences;
	flag_ignore_quality_scores_for_matched_bases = arguments.flag_ignore_quality_scores_for_matched_bases;
	run_diagnostics = arguments.run_diagnostics;
	max_input_reads_in_a_single_nucl_loc = arguments.max_input_reads_in_a_single_nucl_loc;
	skip_shortening_read_names = arguments.skip_shortening_read_names;
	printf( "\nskip_shortening_read_names = %d", skip_shortening_read_names );

	/********************************************************************/

	/*
	 * If user requests no sequence information then everything else is also ignored
	 */
	readAlignmentsAndCompress(
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
			flag_ignore_alignment_scores,
			skip_shortening_read_names );
	return 0;
}
