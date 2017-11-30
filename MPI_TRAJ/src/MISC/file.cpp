#include "file.h"

FileReader::FileReader( string filename, Parameters* parameters )
{
	this->parameters = parameters;

	ifstream infile( filename );
	if ( !infile )
	{
		throw std::invalid_argument( "Can't open the file!" );
	}	
	else
	{
		parse_file( infile );
		infile.close();
	}
}

void FileReader::parse_file( ifstream& infile )
{
	string current_string;
	string keyword;
	string variable = "";
	string value;
    
	char buf[ MAXLINE ];
    int line = 1;
	
	bool is_assignment;
	bool is_dollar;
	bool is_empty;	

	string current_group = "";

	while( infile.getline( buf, MAXLINE ) ) 
	{
		current_string = buf;
		parse_string( current_string, keyword, variable, value, is_dollar, is_empty, is_assignment, line );

		//cout << "current_string: " << current_string << endl;
		//cout << "current_group: " << current_group << endl;
		//cout << "keyword: " << keyword << endl;
		//cout << "variable: " << variable << endl;
		//cout << "value: " << value << endl;
		//cout << "is_dollar: " << is_dollar << endl;
		//cout << "is_empty: " << is_empty << endl;
		//cout << "is_assignment: " << is_assignment << endl;
		//cout << "line: " << line << endl << endl;

		// if line is empty then proceeding to the next one
		if ( is_empty )
		{
			continue;
		}

		// all assignments should be enclosed in some group
		if ( is_dollar && current_group == "" ) 
		{
			current_group = keyword;
			keyword = "";
		} 
	
		// if assignment occured outside of any group block then throw an exception
		if ( is_assignment && current_group == "" ) 
		{
			throw std::invalid_argument( "Expecting '$group' keyword!" );
		}

		// if keyword on the current line is end
		// then current group should be closed (if only it wasn't already closed, in which case we should throw an error) 
		if ( keyword == "$end" )
		{
			if ( current_group == "" )
			{
				string line_number_string = std::to_string( line );
				throw std::invalid_argument( "Invalid syntax: Occured end-group without begin-group! Line: " + line_number_string );
			}

			current_group = "";
			continue;
		}
		
		if ( current_group == "$gridparameters" )
		{
			if ( is_assignment ) 
			{
				analyse_grid_group_line( variable, value, line );

				variable.erase();
				value.erase();	
			}
			else continue;
		}
		else if ( current_group == "$mcparameters" )
		{
			if ( is_assignment )
			{
				analyse_mcparameters_group_line( variable, value, line );

				variable.erase();
				value.erase();
			}
		}
		else if ( current_group == "$files" )
		{
			if ( is_assignment )
			{
				analyse_files_group_line( variable, value, line );

				variable.erase();
				value.erase();
			}
		}
		else if ( current_group == "$trajectory" )
		{
			if ( is_assignment )
			{
				analyse_trajectory_group_line( variable, value, line );

				variable.erase();
				value.erase();
			}
		}
		else if ( current_group == "$conditions" )
		{
			if ( is_assignment )
			{
				analyse_conditions_group_line( variable, value, line );

				variable.erase();
				value.erase();
			}
		}
		else 
		{
			string line_number_string = std::to_string( line );
			throw std::invalid_argument( "Unexpected group name: " + current_group + " on line " + line_number_string );
		}
		
		line++;
    }
}
 
double FileReader::string_to_double( string& value, int& line ) 
{
	double result;
	if ( sscanf( value.c_str(), "%lg", &result ) != 1 ) 
	{
		string line_number_string = std::to_string( line );
		throw std::invalid_argument("Can't transform string to double in line " + line_number_string );
	}

	return result;
}

int FileReader::string_to_int( string& value, int& line )
{
	int result;
	if ( sscanf( value.c_str(), "%d", &result ) != 1 )
	{
		string line_number_string = std::to_string( line );
		throw std::invalid_argument( "Can't transform string to int in line " + line_number_string ); 
	}

	return result;
}

void FileReader::analyse_grid_group_line( string& variable, string& value, int& line )
{
	if ( variable == "V0_MIN" ) this->parameters->V0_MIN = string_to_double( value, line );
	else if ( variable == "V0_MAX" ) this->parameters->V0_MAX = string_to_double( value, line );
	else if ( variable == "V0_PARTS" ) this->parameters->V0_PARTS = string_to_int( value, line );
	else if ( variable == "B_MIN" ) this->parameters->B_MIN = string_to_double( value, line );
	else if ( variable == "B_MAX" ) this->parameters->B_MAX = string_to_double( value, line );
	else if ( variable == "B_PARTS" ) this->parameters->B_PARTS = string_to_int( value, line );
	else
	{
		string line_number_string = std::to_string( line );
		throw std::invalid_argument( "Invalid variable name " + variable + " in $grid group! Line number: " + line_number_string );
	}
}

void FileReader::analyse_mcparameters_group_line( string& variable, string& value, int& line )
{
	if ( variable == "CYCLES" ) this->parameters->CYCLES = string_to_double( value, line );
	else if ( variable == "CYCLE_POINTS" ) this->parameters->CYCLE_POINTS = string_to_double( value, line );
	else if ( variable == "STDDEV_MAX" ) this->parameters->STDDEV_MAX = string_to_double( value, line );
	else
	{
		string line_number_string = std::to_string( line );
		throw std::invalid_argument( "Invalid variable name " + variable + " in $mcparameters group! Line number: " + line_number_string );
	}
}

void FileReader::analyse_files_group_line( string& variable, string& value, int& line )
{
	if ( variable == "specfunc_filename" ) this->parameters->specfunc_filename.assign( value );
	else if ( variable == "spectrum_filename" ) this->parameters->spectrum_filename.assign( value ); 
	else if ( variable == "m2_filename" ) this->parameters->m2_filename.assign( value );
	else
	{
		string line_number_string = std::to_string( line );
		throw std::invalid_argument( "Invalid variable name " + variable + " in $files group! Line number: " + line_number_string );
	}
}

void FileReader::analyse_trajectory_group_line( string& variable, string& value, int& line )
{
	if ( variable == "RDIST" ) this->parameters->RDIST = string_to_double( value, line );
	else if ( variable == "sampling_time" ) this->parameters->sampling_time = string_to_double( value, line );
	else if ( variable == "MaxTrajectoryLength" ) this->parameters->MaxTrajectoryLength = string_to_int( value, line );
	else if( variable == "FREQ_MAX" ) this->parameters->FREQ_MAX = string_to_double( value, line );
	else
	{
		string line_number_string = std::to_string( line );
		throw std::invalid_argument( "Invalid variable name " + variable + " in $trajectory group! Line number: " + line_number_string );
	}
}

void FileReader::analyse_conditions_group_line( string& variable, string& value, int& line )
{
	if ( variable == "Temperature" ) this->parameters->Temperature = string_to_double( value, line );
	else
	{
		string line_number_string = std::to_string( line );
		throw std::invalid_argument( "Invalid variable name " + variable + " in $conditions group! Line number: " + line_number_string );
	}
}

void FileReader::parse_string( string curr_str, string& keyword, string& variable, string& value, bool& is_dollar, bool& is_empty, bool& is_assignment, int& line ) 
{
	// remove comments from line
	size_t pos = curr_str.find("%");
	if( pos != string::npos ) 
	{
    	curr_str.erase( pos, curr_str.size() - pos );
	};

    string space_symbols = "\n \t";

	// goes through string searching for everything except for space symbols
	pos = curr_str.find_first_not_of( space_symbols );

	// if it stopped in some place (not in the end) then string is not empty
    if( pos == string::npos ) 
	{
    	is_empty = true;
    	is_dollar = false;
		is_assignment = false;
		
		keyword.erase();
		variable.erase();
		value.erase();

		return;
    }
	else
	{
		// if we progressed to this line, then the string is not empty
    	is_empty = false;
	}

	// searching for dollar sign
    pos = curr_str.find("$");

    size_t start_keyword, end_keyword;

    // check if "$" occured
    if( pos != string::npos ) 
	{
        // found $
        is_dollar = true;
		is_assignment = false;

		// isolating keyword 
		start_keyword = curr_str.find_first_not_of( space_symbols );
        end_keyword = curr_str.find_last_not_of( space_symbols );

        keyword = curr_str.substr( pos, end_keyword - start_keyword + 1 ); 
        value = "";

        return;
    }
   	// if we are in 'else', then the dollar is not found	
	// because string is not empty, it should contain some assignment
	// if there is not assignment symbol(=) then throw an exception 
	else 
	{
        is_dollar = false;

		// looking for equality sign
        pos = curr_str.find("=");

		// suppose we found assignment
        if ( pos != string::npos ) 
		{
			// setting an assignment indicator to true
			is_assignment = true;

			// left-hand side must be a variable
			// right-hand side must be it's value
            string lhs, rhs;

			// putting left-hand side of equation into lhs
            lhs = curr_str.substr( 0, pos );
			
			size_t variable_start, variable_end;
			
			// here starts variable
            variable_start = lhs.find_first_not_of( space_symbols );
			// here it ends
            variable_end = lhs.find_last_not_of( space_symbols );
           
		    // so between variable_start and variable_end lies the variable 	
			variable = lhs.substr( variable_start, variable_end - variable_start + 1 );
            
			
			// analyzing right-hand side
            rhs = curr_str.substr( pos + 1,  curr_str.length() - pos - 1 );
            
			size_t value_start, value_end;
			
			value_start = rhs.find_first_not_of( space_symbols );
            value_end = rhs.find_last_not_of( space_symbols );
            value = rhs.substr( value_start, value_end - value_start + 1 );

            return;
        } 
		else 
		{
			string line_number_string = std::to_string( line );
            throw std::invalid_argument( "There is no '$' or '=' symbol in line " + line_number_string );
        }
    }
}

FileReader::~FileReader()
{
	//cout << "File destructor" << endl;
}
