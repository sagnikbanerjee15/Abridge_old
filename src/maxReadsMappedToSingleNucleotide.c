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
{ "output_filename", 'o', "TEXT_FILENAME", 0, "Enter the name of the output file that will contain the value of maximum number of reads mapped to a single nucleotide", 0 },
{ "name_of_total_number_of_alignments_filename", 'a', "TEXT_FILENAME", 0, "Enter the name of the total number of alignments file", 0 },
{ "name_of_file_max_read_length", 'm', "TEXT_FILENAME", 0, "Enter the name of the file where the maximum read length will be recorded", 0 },
{ 0, 0, 0, 0, 0, 0 } // Last entry should be all zeros in all fields
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
	// char *args[0];   // No standard arguments (without flags)
	char *input_sam_filename; // Empty string - only contains null character
	char *output_filename;
	char *name_of_total_number_of_alignments_filename;
	char *name_of_file_max_read_length;
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

	// Figure out which option we are parsing, and decide how to store it
	switch ( key )
	{
		case 'i':
			arguments->input_sam_filename = arg;
			break;
		case 'o':
			arguments->output_filename = arg;
			break;
		case 'a':
			arguments->name_of_total_number_of_alignments_filename = arg;
			break;
		case 'm':
			arguments->name_of_file_max_read_length = arg;
			break;

		case ARGP_KEY_END:
			// Reached the last key.
			// Check if our inputsamfilename and outputfilename REQUIRED "options" have been set to non-default values
			if ( strcmp( arguments->input_sam_filename, "" ) == 0
					|| strcmp( arguments->output_filename, "" ) == 0
					|| strcmp(
							arguments->name_of_total_number_of_alignments_filename,
							"" ) == 0
					|| strcmp( arguments->name_of_file_max_read_length, "" )
							== 0 )
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

void findMaximumNumberOfReadsMappedToOneNucleotide(
		char *input_samfilename,
		char *output_filename,
		char *name_of_total_number_of_alignments_filename,
		char *name_of_file_max_read_length )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i, j, k;
	int number_of_fields;
	int max_read_length;

	long long int max_position, max_value;
	long long int prev_position, prev_value;
	long long int curr_position, curr_value;
	long long int total_number_of_alignments;

	FILE *fhr;
	FILE *fhw;
	FILE *fhw_tot_alignments;
	FILE *fhw_max_read_length;

	size_t len = 0;
	ssize_t line_len;

	char str[100];
	char *temp; //Useless
	char *line = NULL; // for reading each line
	char **split_line; // List of strings to store each element of a single alignment
	char **split_tags;

	//struct Sam_Alignment *curr_alignment;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	//curr_alignment = allocateMemorySam_Alignment ();
	fhr = fopen( input_samfilename, "r" );
	if ( fhr == NULL )
	{
		printf( "Error! File %s not found", input_samfilename );
		exit( 1 );
	}
	fhw = fopen( output_filename, "w" );
	if ( fhw == NULL )
	{
		printf( "%s File cannot be created", output_filename );
		exit( 1 );
	}
	fhw_max_read_length = fopen( name_of_file_max_read_length, "w" );
	if ( fhw_max_read_length == NULL )
	{
		printf( "%s File cannot be created", name_of_file_max_read_length );
		exit( 1 );
	}
	fhw_tot_alignments = fopen(
			name_of_total_number_of_alignments_filename,
			"w" );
	if ( fhw_tot_alignments == NULL )
	{
		printf(
				"%s File cannot be created",
				name_of_total_number_of_alignments_filename );
		exit( 1 );
	}

	max_position = 0;
	max_value = 0;
	curr_position = 0;
	curr_value = 0;
	prev_position = 0;
	prev_value = 0;
	max_read_length = 0;

	split_line = ( char** ) malloc( sizeof(char*) * ROWS );
	for ( i = 0; i < ROWS; i++ )
		split_line[i] = ( char* ) malloc( sizeof(char) * COLS );

	split_tags = ( char** ) malloc( sizeof(char*) * ROWS );
	for ( i = 0; i < ROWS; i++ )
		split_tags[i] = ( char* ) malloc( sizeof(char) * COLS );
	/********************************************************************/

	while ( (line_len = getline( &line, &len, fhr )) != -1 )
		if ( line[0] != '@' )
			break;

	total_number_of_alignments = 0;
	do
	{
		total_number_of_alignments += 1;
		number_of_fields = splitByDelimiter( line, '\t', split_line );
		if ( max_read_length < strlen( split_line[9] ) )
			max_read_length = strlen( split_line[9] );
		//populateSamAlignmentInstance ( curr_alignment , split_line , number_of_fields , split_tags );

		curr_position = strtol( split_line[3], &temp, 10 );
		if ( curr_position == 0 )
			continue;
		if ( max_position == 0 )
		{
			max_position = curr_position;
			prev_position = curr_position;

			max_value = 1;
			curr_value = 1;
			prev_value = 1;
		}
		else
		{
			if ( prev_position == curr_position )
				prev_value += 1;
			else
			{
				if ( prev_value > max_value )
				{
					max_value = prev_value;
					max_position = prev_position;
				}
				prev_position = curr_position;
				prev_value = 1;
			}
		}

	} while ( (line_len = getline( &line, &len, fhr )) != -1 );
	if ( prev_value > max_value )
	{
		max_value = prev_value;
		max_position = prev_position;
		prev_position = curr_position;
		prev_value = 1;
	}

	sprintf( str, "%lld", max_value );
	strcat( str, "\n" );
	fprintf( fhw, "%s", str );

	sprintf( str, "%lld", total_number_of_alignments );
	strcat( str, "\n" );
	fprintf( fhw_tot_alignments, "%s", str );

	sprintf( str, "%lld", max_read_length );
	strcat( str, "\n" );
	fprintf( fhw_max_read_length, "%s", str );

	fclose( fhw );
	fclose( fhw_tot_alignments );
	fclose( fhw_max_read_length );
	fclose( fhr );
}

int main( int argc, char *argv[] )
{
	/********************************************************************
	 * Named CLI
	 ********************************************************************/
	struct arguments arguments;

	// Parse our arguments; every option seen by parse_opt will be reflected in arguments.
	// Default values.
	arguments.inputsamfilename = ""; // Empty string - only contains null character
	arguments.outputfilename = "";

	argp_parse( &argp, argc, argv, 0, 0, &arguments );
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char input_samfilename[FILENAME_LENGTH];
	char output_filename[FILENAME_LENGTH];
	char name_of_total_number_of_alignments_filename[FILENAME_LENGTH];
	char name_of_file_max_read_length[FILENAME_LENGTH];

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	strcpy( input_samfilename, arguments.input_sam_filename );
	strcpy( output_filename, arguments.output_filename );
	strcpy(
			name_of_total_number_of_alignments_filename,
			arguments.name_of_total_number_of_alignments_filename );
	strcpy(
			name_of_file_max_read_length,
			arguments.name_of_file_max_read_length );

	/********************************************************************/
	findMaximumNumberOfReadsMappedToOneNucleotide(
			input_samfilename,
			output_filename,
			name_of_total_number_of_alignments_filename,
			name_of_file_max_read_length );

}
