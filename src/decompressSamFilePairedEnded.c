# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include <argp.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

long long int total_mapped_reads = 0;

// Set up the argument parser
const char *argp_program_version = "abridge decompressSamFilePairedEnded 1.2.0";
const char *argp_program_bug_address = "sagnikbanerjee15@gmail.com";
static char doc[] = "decompressSamFilePairedEnded will accept an compressed file and associated quality scores and decompress those to SAM alignments";
static char args_doc[] = "";  // No standard arguments
							  // (i.e. arguments without "names")

static struct argp_option options[] =
{
{ "reference", 'r', "referencefilename", 0, "Name of the reference file in fasta format", 0 },
{ "outputfilename", 'o', "outputfilename", 0, "Name of the samfile where the output will be stored", 0 },
{ "compressedfile", 'c', "compressedfilename", 0, "Name of the compressed file generated by the compression algorithm", 0 },
{ "mockquality", 'm', "mockqualityscore", 0, "Value of the mock quality score", 0 },
{ "ignoresequence", 's', 0, 0, "Flag to ignore sequence generation during decompression", 0 },
{ "dictionaryfilename", 'd', "dictionaryfilename", 0, "Name of the dictionary", 0 },
{ "unmappedreadsfilename", 'u', "unmappedreadsfilename", 0, "Name of the file with unmapped reads", 0 },
{ "qualityscoresfilename", 'q', "qualityscoresfilename", 0, "Name of the file with quality scores", 0 },
{ "maxreadsineachline", 'x', "maxreadsineachline", 0, "Maximum number of reads mapped to a single reference nucleotide position", 0 },
{ 0, 0, 0, 0, 0, 0 } // Last entry should be all zeros in all fields
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

static error_t parse_opt( int key, char *arg, struct argp_state *state )
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
		case 'd':
			arguments->dictionaryfilename = arg;
			break;
		case 'x':
			arguments->maxreadsineachline = strtoull( arg, &temp, 10 );
			break;

		case ARGP_KEY_END:
			// Reached the last key.
			// Check if our inputsamfilename and outputfilename REQUIRED "options" have been set to non-default values
			if ( strcmp( arguments->reference, "" ) == 0
					|| strcmp( arguments->outputfilename, "" ) == 0
					|| strcmp( arguments->compressedfile, "" ) == 0
					|| strcmp( arguments->unmappedreadsfilename, "" ) == 0
					|| strcmp( arguments->qualityscoresfilename, "" ) == 0
					|| strcmp( arguments->dictionaryfilename, "" ) == 0
					|| arguments->maxreadsineachline == 0 )
			{
				argp_usage( state );
			}
			break;

		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

static struct argp argp =
{ options, parse_opt, args_doc, doc, 0, 0, 0 };

void writeAlignmentToFilePairedEnded(
		struct Sam_Alignment *sam_alignment,
		short int flag_ignore_sequence_information,
		int number_of_repititions_of_the_same_reads,
		FILE *fhw,
		FILE *fhr_qual,
		short int flag_ignore_quality_scores_for_matched_bases,
		char **read_names,
		int *read_names_index,
		short int flag_ignore_alignment_scores,
		int number_of_reads )
{
	int i;

	char line_to_be_written_to_file[MAX_GENERAL_LEN];
	char temp[100];
	char *buffer;

	size_t len = 0;
	ssize_t line_len;

	for ( i = 0; i < number_of_repititions_of_the_same_reads; i++ )
	{
		line_to_be_written_to_file[0] = '\0';
		strcat( line_to_be_written_to_file, read_names[*read_names_index] );
		(*read_names_index)++;
		if ( *read_names_index > number_of_reads )
		{
			printf( "\nRead index exceeded" );
			exit( 1 );
		}
		else
		{
			//printf ("\nIndex=%d Total=%d" , *read_names_index , number_of_reads);
		}
		//sprintf (temp , "%d" , i + 1);
		//strcat (line_to_be_written_to_file , "_");
		//strcat (line_to_be_written_to_file , temp);
		strcat( line_to_be_written_to_file, "\t" );

		sprintf( temp, "%d", sam_alignment->samflag );
		strcat( line_to_be_written_to_file, temp );

		strcat( line_to_be_written_to_file, "\t" );
		strcat( line_to_be_written_to_file, sam_alignment->reference_name );

		strcat( line_to_be_written_to_file, "\t" );
		sprintf( temp, "%lld", sam_alignment->start_position );
		strcat( line_to_be_written_to_file, temp );

		strcat( line_to_be_written_to_file, "\t" );
		if ( flag_ignore_alignment_scores == 1 )
			strcat( line_to_be_written_to_file, "255" );
		else
		{
			sprintf( temp, "%d", sam_alignment->mapping_quality_score );
			strcat( line_to_be_written_to_file, temp );
		}

		strcat( line_to_be_written_to_file, "\t" );
		strcat( line_to_be_written_to_file, sam_alignment->cigar );

		strcat( line_to_be_written_to_file, "\t" );
		strcat( line_to_be_written_to_file, "*" );

		strcat( line_to_be_written_to_file, "\t" );
		strcat( line_to_be_written_to_file, "0" );

		strcat( line_to_be_written_to_file, "\t" );
		strcat( line_to_be_written_to_file, "0" );

		strcat( line_to_be_written_to_file, "\t" );
		strcat( line_to_be_written_to_file, sam_alignment->seq );

		strcat( line_to_be_written_to_file, "\t" );
		if ( flag_ignore_quality_scores_for_matched_bases == 0
				&& (line_len = getline( &buffer, &len, fhr_qual )) != -1 )
		{
			buffer[strlen( buffer ) - 1] = '\0';
			strcat( line_to_be_written_to_file, buffer );
		}
		else
			strcat( line_to_be_written_to_file, sam_alignment->qual );

		strcat( line_to_be_written_to_file, "\t" );
		//Tags
		strcat( line_to_be_written_to_file, "NH:i:" );
		strcat( line_to_be_written_to_file, sam_alignment->tags[0].val );
		strcat( line_to_be_written_to_file, "\t" );

		if ( strcmp( sam_alignment->tags[1].val, "." ) != 0
				&& strchr( sam_alignment->cigar, 'N' ) != NULL )
		{
			strcat( line_to_be_written_to_file, "XS:A:" );
			strcat( line_to_be_written_to_file, sam_alignment->tags[1].val );
			strcat( line_to_be_written_to_file, "\t" );
		}

		if ( flag_ignore_sequence_information == 0 )
		{
			strcat( line_to_be_written_to_file, "MD:Z:" );
			strcat( line_to_be_written_to_file, sam_alignment->tags[2].val );
			strcat( line_to_be_written_to_file, "\t" );
		}
		if ( flag_ignore_alignment_scores == 0
				&& strcmp( sam_alignment->tags[3].val, "X" ) != 0 )
		{
			strcat( line_to_be_written_to_file, "AS:i:" );
			strcat( line_to_be_written_to_file, sam_alignment->tags[3].val );
			strcat( line_to_be_written_to_file, "\t" );
		}

		if ( line_to_be_written_to_file[strlen( line_to_be_written_to_file ) - 1]
				== '\t' )
			line_to_be_written_to_file[strlen( line_to_be_written_to_file ) - 1] = '\0';
		strcat( line_to_be_written_to_file, "\n" );
		fprintf( fhw, "%s", line_to_be_written_to_file );
	}
}

void convertToAlignmentPairedEnded(
		struct Sam_Alignment *sam_alignment_instance,
		struct Whole_Genome_Sequence *whole_genome,
		char **split_on_tab,
		char **split_on_dash,
		char **split_on_comma,
		char **split_on_tilde,
		char **read_names,
		char *default_quality_value,

		short int flag_ignore_alignment_scores,
		short int flag_ignore_mismatches,
		short int flag_ignore_soft_clippings,
		short int flag_ignore_unmapped_sequences,
		short int flag_ignore_all_quality_scores,
		short int flag_ignore_sequence_information,
		short int flag_ignore_quality_scores_for_matched_bases,

		unsigned long long int *read_number,
		unsigned long long int *total_mapped_reads,
		FILE *fhw,
		FILE *fhr_qual,
		short int number_of_columns,
		unsigned long long int curr_position,
		char *chromosome,
		struct Paired_Ended_Flag_to_Single_Character *samflag_dictionary,
		int number_of_unique_samformatflags,
		char samformatflag_replacer_characters[] )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/

	int number_of_distinct_cigars_in_a_line;
	int number_of_repititions_of_the_same_reads;
	int samformatflag;
	int number_of_reads;

	int i, j;
	int read_names_index = 0;

	char *temp; //Useless
	char *distinct_icigars_in_a_line;
	char *icigar;
	char str_sprintf[50];

	short int number_of_items_separated_by_underscore;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	/********************************************************************/

	if ( number_of_columns == 2 )
	{
		//return;
		int max_number_of_commas = 0, number_of_commas = 0;
		for ( i = 0; split_on_tab[0][i] != '\0'; i++ )
			if ( split_on_tab[0][i] == ',' )
				number_of_commas++;
		number_of_distinct_cigars_in_a_line = splitByDelimiter(
				split_on_tab[0],
				',',
				split_on_comma );
		number_of_reads = splitByDelimiter( split_on_tab[1], ',', read_names );
	}
	else if ( number_of_columns == 3 )
	{
		int max_number_of_commas = 0, number_of_commas = 0;
		for ( i = 0; split_on_tab[1][i] != '\0'; i++ )
			if ( split_on_tab[1][i] == ',' )
				number_of_commas++;
		/*if ( strcmp (chromosome , "Pt") == 0 )
		 {
		 printf ("\nsplit_on_tab %s" , split_on_tab[1]);
		 fflush (stdout);
		 }*/
		number_of_distinct_cigars_in_a_line = splitByDelimiter(
				split_on_tab[1],
				',',
				split_on_comma );
		number_of_reads = splitByDelimiter( split_on_tab[2], ',', read_names );
	}

	//return;
	for ( j = 0; j < number_of_distinct_cigars_in_a_line; j++ )
	{
		splitByDelimiter( split_on_comma[j], '-', split_on_dash );
		if ( flag_ignore_alignment_scores == 0 )
			splitByDelimiter( split_on_dash[0], '~', split_on_tilde );
		else
			strcpy( split_on_tilde[0], split_on_dash[0] );
		if ( split_on_comma[j][strlen( split_on_comma[j] ) - 1] != '-' )
		{
			number_of_repititions_of_the_same_reads = strtol(
					split_on_dash[1],
					&temp,
					10 );
		}
		else
			number_of_repititions_of_the_same_reads = 1;

		if ( !(split_on_comma[j][1] == '-'
				&& isalpha( split_on_dash[0][0] ) != 0) )
		{
			if ( flag_ignore_alignment_scores == 0 )
			{
				sam_alignment_instance->mapping_quality_score = strtol(
						split_on_tilde[1],
						&temp,
						10 );
				strcpy(
						sam_alignment_instance->tags[3].val,
						split_on_tilde[2] );
			}
			else
			{
				sam_alignment_instance->mapping_quality_score = 255;
				strcpy( sam_alignment_instance->tags[3].val, "X" );
			}
		}
		sam_alignment_instance->start_position = curr_position;

		if ( split_on_comma[j][1] == '-'
				&& isalpha( split_on_dash[0][0] ) != 0 )
		{
			// Use the same cigar - basically do nothing

			/*
			 sprintf (temp , "%d" , *read_number);
			 ( *read_number )++;
			 strcpy (sam_alignment_instance->read_name , temp);
			 */
			/*if ( sam_alignment_instance->start_position == 27381 && strcmp (sam_alignment_instance->reference_name , "1") == 0 )
			 {
			 printf ("\nSame iCIGAR");
			 printf ("\nsplit_on_comma[j] = %s" , split_on_comma[j]);
			 printf ("\nWeird Location cigar %s number_of_repititions_of_the_same_reads %d" , sam_alignment_instance->icigar , number_of_repititions_of_the_same_reads);
			 printf ("\nMD String: %s" , sam_alignment_instance->tags[2].val);
			 printf ("\n==============================================================================================================================");
			 fflush (stdout);
			 }*/
		}
		else
		{
			strcpy( sam_alignment_instance->icigar, split_on_tilde[0] );
			//printf ("\nj=%d number_of_distinct_cigars_in_a_line=%d Inside ICIGAR %s" , j , number_of_distinct_cigars_in_a_line , sam_alignment_instance->icigar);
			//fflush (stdout);

			convertIcigarToCigarandMDPairedEnded(
					whole_genome,
					sam_alignment_instance,
					chromosome,
					flag_ignore_mismatches,
					flag_ignore_soft_clippings,
					flag_ignore_unmapped_sequences,
					flag_ignore_all_quality_scores,
					flag_ignore_quality_scores_for_matched_bases,
					flag_ignore_sequence_information,
					default_quality_value,
					samflag_dictionary,
					number_of_unique_samformatflags,
					samformatflag_replacer_characters );

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

			/*
			 sprintf (temp , "%d" , *read_number);
			 ( *read_number )++;
			 strcpy (sam_alignment_instance->read_name , temp);
			 */
		}
		writeAlignmentToFilePairedEnded(
				sam_alignment_instance,
				flag_ignore_sequence_information,
				number_of_repititions_of_the_same_reads,
				fhw,
				fhr_qual,
				flag_ignore_quality_scores_for_matched_bases,
				read_names,
				&read_names_index,
				flag_ignore_alignment_scores,
				number_of_reads );
		(*total_mapped_reads) += number_of_repititions_of_the_same_reads;
	}
}

void decompressFile(
		char *name_of_file_with_quality_scores,
		char *genome_filename,
		char *output_sam_filename,
		char *pass1_filename,
		char *unmapped_filename,
		char *default_quality_value,
		short int flag_ignore_sequence_information,
		char *dictionary_name,
		int max_reads_in_each_line )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhr_dictionary;
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
	short int flag_ignore_quality_score;
	short int flag_ignore_quality_scores_for_matched_bases;
	short int flag_ignore_all_quality_scores;
	short int flag_ignore_alignment_scores;
	short int flag_skip_shortening_read_names;
	short int number_of_columns;

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
	int ROWS_split_on_tilde = 5;
	int COLS_split_on_tilde = MAX_SEQ_LEN * 3; //3000
	int number_of_unique_samformatflags;

	//char **split_on_newline;
	char **split_on_tab;
	char **split_on_dash;
	char **split_on_comma;
	char **split_on_tilde;
	char **read_names;
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
	char samformatflag_replacer_characters[] = "BEFHIJKLMOPQRSUVWXYZbdefhijklmopqrsuvwxyz";

	struct Sam_Alignment **sam_alignment_pool;
	struct Sam_Alignment *sam_alignment_instance;
	struct Whole_Genome_Sequence *whole_genome;
	struct Paired_Ended_Flag_to_Single_Character *samflag_dictionary;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen( pass1_filename, "r" );
	if ( fhr == NULL )
	{
		printf( "1. Error! File %s not found", pass1_filename );
		exit( 1 );
	}
	fhr_dictionary = fopen( dictionary_name, "r" );
	if ( fhr_dictionary == NULL )
	{
		printf( "2. Error! File %s not found", dictionary_name );
		exit( 1 );
	}
	fhr_qual = fopen( name_of_file_with_quality_scores, "r" );
	if ( fhr_qual == NULL )
	{
		printf(
				"3. Error! File not found %s",
				name_of_file_with_quality_scores );
		exit( 1 );
	}
	fhw = fopen( output_sam_filename, "w" );
	if ( fhw == NULL )
	{
		printf(
				"4. Error! File %s cannot be opened for writing",
				output_sam_filename );
		exit( 1 );
	}

	/*
	 split_on_newline = ( char** ) malloc (sizeof(char*) * ROWS_split_on_newline);
	 for ( i = 0 ; i < ROWS_split_on_newline ; i++ )
	 split_on_newline[i] = ( char* ) malloc (sizeof(char) * COLS_split_on_newline);
	 */
	split_on_tab = ( char** ) malloc( sizeof(char*) * ROWS_split_on_tab );
	for ( i = 0; i < ROWS_split_on_tab; i++ )
		split_on_tab[i] = ( char* ) malloc( sizeof(char) * COLS_split_on_tab );

	max_reads_in_each_line += 10;
	read_names = ( char** ) malloc( sizeof(char*) * max_reads_in_each_line );
	for ( i = 0; i < max_reads_in_each_line; i++ )
		read_names[i] = ( char* ) malloc( sizeof(char) * 100 );

	split_on_dash = ( char** ) malloc( sizeof(char*) * ROWS_split_on_dash );
	for ( i = 0; i < ROWS_split_on_dash; i++ )
		split_on_dash[i] = ( char* ) malloc(
				sizeof(char) * COLS_split_on_dash );

	split_on_comma = ( char** ) malloc(
			sizeof(char*) * max_reads_in_each_line );
	for ( i = 0; i < max_reads_in_each_line; i++ )
		split_on_comma[i] = ( char* ) malloc(
				sizeof(char) * COLS_split_on_comma );

	split_on_tilde = ( char** ) malloc( sizeof(char*) * ROWS_split_on_tilde );
	for ( i = 0; i < ROWS_split_on_tilde; i++ )
		split_on_tilde[i] = ( char* ) malloc(
				sizeof(char) * COLS_split_on_tilde );

	output_prefix_without_path = ( char* ) malloc( sizeof(char) * MAX_SEQ_LEN );
	sequence_portions_from_reference = ( char** ) malloc(
			sizeof(char*) * MAX_POOL_SIZE );
	fasta_file_with_expressed_portions = ( char* ) malloc(
			sizeof(char) * FILENAME_LENGTH );
	current_chromosome = ( char* ) malloc( sizeof(char) * 100 );

	//buffer = ( char* ) malloc (sizeof(char) * BUFFER_SIZE);
	abridge_index = allocateMemoryAbridge_Index();
	sam_alignment = allocateMemorySam_Alignment();
	whole_genome = ( struct Whole_Genome_Sequence* ) malloc(
			sizeof(struct Whole_Genome_Sequence) );

	whole_genome->number_of_reference_sequences = 0;
	whole_genome->nucleotides = ( char** ) malloc( sizeof(char*) * 1 );
	whole_genome->reference_sequence_name = ( char** ) malloc(
			sizeof(char*) * 1 );
	whole_genome->reference_sequence_length = ( unsigned long long int* ) malloc(
			sizeof(unsigned long long int) * 1 );

	sam_alignment_instance = allocateMemorySam_Alignment();
	read_prefix[0] = '\0'; // Empty string

	samflag_dictionary = allocateMemoryPaired_Ended_Flag_to_Single_Character(
			10000 );
	number_of_unique_samformatflags = fillUpDictionary(
			samflag_dictionary,
			fhr_dictionary,
			10000 );
	number_of_unique_samformatflags /= 2;
	/********************************************************************/

	writeSequenceHeaders( fhw, genome_filename, 1 );
	line_num = 0;
	line_len = getline( &buffer, &len, fhr );
	splitByDelimiter( buffer, '\t', split_on_tab );

	splitByDelimiter( split_on_tab[0], ':', split_on_tilde );
	flag_ignore_mismatches = strtol(
			split_on_tilde[1],
			&convert_to_int_temp,
			10 );

	splitByDelimiter( split_on_tab[1], ':', split_on_tilde );
	flag_ignore_soft_clippings = strtol(
			split_on_tilde[1],
			&convert_to_int_temp,
			10 );

	splitByDelimiter( split_on_tab[2], ':', split_on_tilde );
	flag_ignore_unmapped_sequences = strtol(
			split_on_tilde[1],
			&convert_to_int_temp,
			10 );

	splitByDelimiter( split_on_tab[3], ':', split_on_tilde );
	flag_ignore_all_quality_scores = strtol(
			split_on_tilde[1],
			&convert_to_int_temp,
			10 );

	splitByDelimiter( split_on_tab[4], ':', split_on_tilde );
	flag_ignore_quality_scores_for_matched_bases = strtol(
			split_on_tilde[1],
			&convert_to_int_temp,
			10 );

	splitByDelimiter( split_on_tab[6], ':', split_on_tilde );
	flag_ignore_alignment_scores = strtol(
			split_on_tilde[1],
			&convert_to_int_temp,
			10 );

	splitByDelimiter( split_on_tab[6], ':', split_on_tilde );
	flag_skip_shortening_read_names = strtol(
			split_on_tilde[1],
			&convert_to_int_temp,
			10 );
	/*
	 printf ("\nflag_ignore_mismatches %d" , flag_ignore_mismatches);
	 printf ("\nflag_ignore_soft_clippings %d" , flag_ignore_soft_clippings);
	 printf ("\nflag_ignore_unmapped_sequences %d" ,
	 flag_ignore_unmapped_sequences);
	 printf ("\nflag_ignore_quality_score %d" ,
	 flag_ignore_all_quality_scores);
	 printf ("\nflag_save_all_quality_scores %d" ,
	 flag_ignore_quality_scores_for_matched_bases);
	 printf ("\nflag_save_exact_quality_scores %d" ,
	 flag_ignore_alignment_scores);
	 fflush (stdout);
	 */
	line_num = 0;
	while ( (line_len = getline( &buffer, &len, fhr )) != -1 )
	{
		line_num++;
		if ( buffer[0] == '@' )
		{
			splitByDelimiter( buffer, '\t', split_on_tab );
			splitByDelimiter( split_on_tab[1], ':', split_on_dash ); // Using split_on_dash so as to save memory and not create a new data structure
			strcpy( current_chromosome, split_on_dash[1] );
			//printf ("\nProcessing chromosome %s" , current_chromosome);
			readInEachChromosome(
					genome_filename,
					whole_genome,
					current_chromosome );
			do
			{
				line_len = getline( &buffer, &len, fhr );
			} while ( buffer[0] == '@' );
			curr_position = 0;
		}
		number_of_commas_in_each_line = 0;
		for ( i = 0; buffer[i] != '\0'; i++ )
			if ( buffer[i] == ',' )
				number_of_commas_in_each_line++;
		if ( max_number_of_commas < number_of_commas_in_each_line )
			max_number_of_commas = number_of_commas_in_each_line;

		//printf ("\nline_len %d len %d" , line_len , len);
		if ( line_len > COLS_split_on_tab )
		{
			//printf ("\nB--> line_len %d COLS_split_on_tab %d" , line_len , COLS_split_on_tab);
			//fflush (stdout);
			for ( i = 0; i < ROWS_split_on_tab; i++ )
				free( split_on_tab[i] );
			COLS_split_on_tab = line_len + 100;
			for ( i = 0; i < ROWS_split_on_tab; i++ )
				split_on_tab[i] = ( char* ) malloc(
						sizeof(char) * COLS_split_on_tab );
			//printf ("\nA--> line_len %d COLS_split_on_tab %d" , line_len , COLS_split_on_tab);
			//fflush (stdout);

		}
		/*
		 if ( max_number_of_commas > ROWS_split_on_comma )
		 {
		 //printf ("\nB--> max_number_of_commas %d ROWS_split_on_comma %d" , max_number_of_commas , ROWS_split_on_comma);
		 //fflush (stdout);
		 for ( i = 0 ; i < ROWS_split_on_comma ; i++ )
		 free (split_on_comma[i]);
		 free (split_on_comma);
		 ROWS_split_on_comma = line_len / max_number_of_commas + 10;
		 split_on_comma = ( char** ) malloc (sizeof(char*) * ROWS_split_on_comma);
		 for ( i = 0 ; i < ROWS_split_on_comma ; i++ )
		 split_on_comma[i] = ( char* ) malloc (sizeof(char) * COLS_split_on_comma);
		 //printf ("\nA--> max_number_of_commas %d ROWS_split_on_comma %d" , max_number_of_commas , ROWS_split_on_comma);
		 //fflush (stdout);
		 }*/
		number_of_columns = splitByDelimiter( buffer, '\t', split_on_tab );
		if ( number_of_columns == 2 )
			curr_position++;
		else if ( number_of_columns == 3 )
			curr_position += strtol(
					split_on_tab[0],
					&convert_to_int_temp,
					10 );
		/*if ( strcmp ("Pt" , current_chromosome) == 0 )
		 {
		 printf ("\n%s" , buffer);
		 fflush (stdout);
		 }*/
		convertToAlignmentPairedEnded(
				sam_alignment_instance,
				whole_genome,
				split_on_tab,
				split_on_dash,
				split_on_comma,
				split_on_tilde,
				read_names,
				default_quality_value,

				flag_ignore_alignment_scores,
				flag_ignore_mismatches,
				flag_ignore_soft_clippings,
				flag_ignore_unmapped_sequences,
				flag_ignore_all_quality_scores,
				flag_ignore_sequence_information,
				flag_ignore_quality_scores_for_matched_bases,

				&read_number,
				&total_mapped_reads,
				fhw,
				fhr_qual,
				number_of_columns,
				curr_position,
				current_chromosome,
				samflag_dictionary,
				number_of_unique_samformatflags,
				samformatflag_replacer_characters );
	}

	/*
	 * Write all unmapped reads to samfile
	 */
	fhr = fopen( unmapped_filename, "r" );
	if ( fhr == NULL )
	{
		printf( "Error! File not found" );
		exit( 1 );
	}
	line_number = 1;
	free( buffer );
	buffer = NULL;
	while ( (line_len = getline( &buffer, &len, fhr )) != -1 )
	{
		strcpy( sam_alignment->seq, buffer );
		line_len = getline( &buffer, &len, fhr );
		strcpy( sam_alignment->qual, buffer );

		line_to_be_written_to_file[0] = '\0';
		strcat( line_to_be_written_to_file, "unmapped_" );
		sprintf( temp, "%d", line_number );
		strcat( line_to_be_written_to_file, temp );
		strcat( line_to_be_written_to_file, "\t" );

		strcat( line_to_be_written_to_file, "4" );
		strcat( line_to_be_written_to_file, "\t" );

		strcat( line_to_be_written_to_file, "*" );
		strcat( line_to_be_written_to_file, "\t" );

		strcat( line_to_be_written_to_file, "0" );
		strcat( line_to_be_written_to_file, "\t" );

		strcat( line_to_be_written_to_file, "0" );

		strcat( line_to_be_written_to_file, "\t" );
		strcat( line_to_be_written_to_file, "*" );

		strcat( line_to_be_written_to_file, "\t" );
		strcat( line_to_be_written_to_file, "*" );

		strcat( line_to_be_written_to_file, "\t" );
		strcat( line_to_be_written_to_file, "0" );

		strcat( line_to_be_written_to_file, "\t" );
		strcat( line_to_be_written_to_file, "0" );

		strcat( line_to_be_written_to_file, "\t" );
		sam_alignment->seq[strlen( sam_alignment->seq ) - 1] = '\0';
		strcat( line_to_be_written_to_file, sam_alignment->seq );

		strcat( line_to_be_written_to_file, "\t" );
		sam_alignment->qual[strlen( sam_alignment->qual ) - 1] = '\0';
		strcat( line_to_be_written_to_file, sam_alignment->qual );

		strcat(
				line_to_be_written_to_file,
				"\tNH:i:0\tHI:i:0\tnM:i:1\tuT:A:1" );

		strcat( line_to_be_written_to_file, "\n" );
		fprintf( fhw, "%s", line_to_be_written_to_file );

		line_number++;
	}

	fclose( fhw );
	fclose( fhr );
}

int main( int argc, char *argv[] )
{

	struct arguments arguments;

// Parse our arguments; every option seen by parse_opt will be reflected in arguments.
// Default values.
	arguments.reference = ""; // Empty string - only contains null character
	arguments.outputfilename = "";
	arguments.compressedfile = "";
	arguments.dictionaryfilename = "";
	arguments.ignoresequence = 0;
	arguments.mockquality = "I";
	arguments.qualityscoresfilename = "";
	arguments.unmappedreadsfilename = "";
	arguments.maxreadsineachline = 0;

	argp_parse( &argp, argc, argv, 0, 0, &arguments );

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
	char dictionary_name[FILENAME_LENGTH];
	char *temp; //Required for strtoi

	short int flag_ignore_sequence_information;
	int max_reads_in_each_line;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	strcpy( genome_filename, arguments.reference );
	strcpy( output_sam_filename, arguments.outputfilename );
	strcpy( pass1_filename, arguments.compressedfile );
	strcpy( default_quality_value, arguments.mockquality );
	flag_ignore_sequence_information = arguments.ignoresequence;
	strcpy( dictionary_name, arguments.dictionaryfilename );
	strcpy( unmapped_filename, arguments.unmappedreadsfilename );
	strcpy( name_of_file_with_quality_scores, arguments.qualityscoresfilename );
	max_reads_in_each_line = arguments.maxreadsineachline;

	/********************************************************************/

	decompressFile(
			name_of_file_with_quality_scores,
			genome_filename,
			output_sam_filename,
			pass1_filename,
			unmapped_filename,
			default_quality_value,
			flag_ignore_sequence_information,
			dictionary_name,
			max_reads_in_each_line );

	return 0;
}
