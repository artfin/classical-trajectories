#pragma once
#include "parameters.h"

#include <iostream>
#include <string>
#include <fstream>

// throwing exception
#include <stdexcept>

using std::cout;
using std::endl;
using std::string;
using std::ifstream;

#define MAXLINE 100 

class FileReader
{
	Parameters* parameters;

	public:
		FileReader( string filename, Parameters* parameters );
		~FileReader();

		void parse_file( ifstream& infile );
	
		double string_to_double( string& value, int& line );
		int string_to_int( string& value, int& line );
		bool string_to_bool( string& value, int& line );

		void parse_string( string curr_str, 
						   string& keyword,
						   string& variable,
						   string& value,
						   bool& is_doolar,
						   bool& is_empty,
						   bool& is_assignment,
						   int& line 
						 );

		void analyse_grid_group_line( string& variable,
						 			  string& value,
									  int& line
					   				 );

		void analyse_mcparameters_group_line( string& variable,
											  string& value,
											  int& line
											 );
		void analyse_files_group_line( string& variable,
									   string& value,
								   	   int& line
									 );
	   	void analyse_trajectory_group_line( string& variable, 
											string& value,
											int& line 
										   );
		void analyse_conditions_group_line( string& variable,
											string& value,
											int& line 
										  );
		void analyse_desymmetrization_group_line( string& variable,
												  string& value,
											  	  int& line
												);	  
};

