# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <argp.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

// Set up the argument parser
const char *argp_program_version = "abridge findAllTags 1.2.0";
const char *argp_program_bug_address = "sagnikbanerjee15@gmail.com";
static char doc[] = "findAllTags will scan a sam file and output the unique list of alignment tags";
static char args_doc[] = "";  // No standard arguments
							  // (i.e. arguments without "names")

/*
 * Options.  Field 1 in ARGP.
 * Order of fields: {NAME, KEY, ARG, FLAGS, DOC, GROUP}.
 */

static struct argp_option options[] =
{
{ "inputsamfilename", 'i', "SAMFILENAME", 0, "Enter the name of the SAM filename", 0 },
{ "outputfilename", 'o', "TEXTFILENAME", 0, "Enter the name of the output filename", 0 },
{ 0, 0, 0, 0, 0, 0 } // Last entry should be all zeros in all fields
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
	// char *args[0];   // No standard arguments (without flags)
	char *inputsamfilename;	  // Argument for --inputsamfilename / -i
	char *outputfilename;  // Argument for --outputfilename / -o
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
			arguments->inputsamfilename = arg;
			break;
		case 'o':
			arguments->outputfilename = arg;
			break;

		case ARGP_KEY_END:
			// Reached the last key.
			// Check if our inputsamfilename and outputfilename REQUIRED "options" have been set to non-default values
			if ( strcmp( arguments->inputsamfilename, "" ) == 0
					|| strcmp( arguments->outputfilename, "" ) == 0 )
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

void extractTagsFromSAMFile( char *inputfilename, char *outputfilename )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;

	size_t len = 0;
	ssize_t line_len;

	char *line = NULL; // for reading each line
	char str[100];
	char temp[100];

	char **split_line; // List of strings to store each element of a single alignment
	char **split_tags;

	int number_of_fields; // Number of fields in each sam alignment entry

	int i, j, k; // Required in loops
	int NH_present, MD_present, XS_present; // Flags to indicate whether these tags are present

	struct Sam_Alignment *curr_alignment;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen( inputfilename, "r" );
	if ( fhr == NULL )
	{
		printf( "Error! File %s not found", inputfilename );
		exit( 1 );
	}
	fhw = fopen( outputfilename, "w" );
	if ( fhw == NULL )
	{
		printf( "%s File cannot be created", outputfilename );
		exit( 1 );
	}

	curr_alignment = allocateMemorySam_Alignment();

	split_line = ( char** ) malloc( sizeof(char*) * ROWS );
	for ( i = 0; i < ROWS; i++ )
		split_line[i] = ( char* ) malloc( sizeof(char) * COLS );

	split_tags = ( char** ) malloc( sizeof(char*) * ROWS );
	for ( i = 0; i < ROWS; i++ )
		split_tags[i] = ( char* ) malloc( sizeof(char) * COLS );
	/********************************************************************/

	// Assume that no tags are present
	NH_present = 0;
	MD_present = 0;
	XS_present = 0;

	while ( (line_len = getline( &line, &len, fhr )) != -1 )
		if ( line[0] != '@' )
			break;

	do
	{
		number_of_fields = splitByDelimiter( line, '\t', split_line );
		populateSamAlignmentInstance(
				curr_alignment,
				split_line,
				number_of_fields,
				split_tags );
		for ( i = 11; i < number_of_fields; i++ )
		{
			if ( strcmp( curr_alignment->tags[i - 11].name, "NH" ) == 0 )
				NH_present = 1;
			else if ( strcmp( curr_alignment->tags[i - 11].name, "MD" ) == 0 )
				MD_present = 1;
			else if ( strcmp( curr_alignment->tags[i - 11].name, "XS" ) == 0 )
				XS_present = 1;
		}
		if ( NH_present * MD_present * XS_present == 1 )
			break;

	} while ( (line_len = getline( &line, &len, fhr )) != -1 );

	fprintf( fhw, "%s", "NH:" );
	sprintf( temp, "%d", NH_present );
	fprintf( fhw, "%s", temp );
	fprintf( fhw, "%s", "\n" );

	fprintf( fhw, "%s", "MD:" );
	sprintf( temp, "%d", MD_present );
	fprintf( fhw, "%s", temp );
	fprintf( fhw, "%s", "\n" );

	fprintf( fhw, "%s", "XS:" );
	sprintf( temp, "%d", XS_present );
	fprintf( fhw, "%s", temp );
	fprintf( fhw, "%s", "\n" );

	fclose( fhr );
	fclose( fhw );
}

int main( int argc, char **argv )
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
	/********************************************************************/
	extractTagsFromSAMFile(
			arguments.inputsamfilename,
			arguments.outputfilename );
	return 0;
}
