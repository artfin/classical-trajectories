#pragma once
#include "gridparameters.h"

#include <iostream>
#include <string>
#include <fstream>

// throwing exception
#include <stdexcept>

using std::cout;
using std::endl;
using std::string;
using std::ifstream;

#define MAXLINE 64

class FileReader
{
	public:
		FileReader( string filename, GridParameters& grid );
		~FileReader();

		void parse_file( ifstream& infile, GridParameters& grid );
		
		void parse_string( string curr_str, 
						   string& keyword,
						   string& variable,
						   string& value,
						   bool& is_doolar,
						   bool& is_empty,
						   int& line 
						 );

		void analyse_grid_group_line( string& variable,
						 			  string& value,
									  int& line,
					   				  GridParameters& grid
									);
};

