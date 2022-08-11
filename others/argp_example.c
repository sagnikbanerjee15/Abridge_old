# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <argp.h>

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
{ "inputsamfilename", 'i', "inputsamfilename", 0, "inputsamfilename of object to download", 0 },
{ "outputfilename", 'o', "outputfilename", 0, "Directory to save downloaded object to", 0 },
{ 0, 0, 0, 0, 0, 0 } // Last entry should be all zeros in all fields
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
	// char *args[0];   // No standard arguments (without flags)
	char *inputsamfilename;	  // Argument for --inputsamfilename / -d
	char *outputfilename;  // Argument for --outputfilename / -d
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

int main( int argc, char **argv )
{
	struct arguments arguments;

	// Parse our arguments; every option seen by parse_opt will be reflected in arguments.
	// Default values.
	arguments.inputsamfilename = ""; // Empty string - only contains null character
	arguments.outputfilename = "";

	argp_parse( &argp, argc, argv, 0, 0, &arguments );

	printf( "User arguments:\n"
			" inputsamfilename=%s\n"
			" outputfilename=%s\n", arguments.inputsamfilename, // This is a pointer to the start of the inputsamfilename char array
			arguments.outputfilename ); // This is a pointer to the start of the destir char array

	return 0;
}
