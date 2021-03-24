# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void convertRLEtoQualValues ( char *input_qualityscore_filename, char *output_quality_score_filename )
{

}

int main ( int argc, char *argv[] )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char *temp;
	char input_qualityscore_filename[FILENAME_LENGTH];
	char output_quality_score_filename[FILENAME_LENGTH];

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy( input_qualityscore_filename , argv[1] );
	strcpy( output_quality_score_filename , argv[2] );
	/********************************************************************/

	convertRLEtoQualValues ( input_qualityscore_filename , output_quality_score_filename );
}
