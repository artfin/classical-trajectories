#include "file.h"

FileReader::FileReader( string filename, GridParameters& grid )
{
	ifstream infile( filename );
	if ( !infile )
	{
		throw std::invalid_argument( "Can't open the file!" );
	}	
	else
	{
		parse_file( infile, grid );
		infile.close();
	}
}

void FileReader::parse_file( ifstream& infile, GridParameters& grid )
{
	cout << "Inside parse_file" << endl;
    
	string current_string;
	string keyword;
	string variable = "";
	string value;
    
	char buf[ MAXLINE ];
    int line = 1;
	
	bool is_dollar;
	bool is_empty;	

	string current_group = "";

	while( infile.getline( buf, MAXLINE ) ) 
	{
	   	is_dollar = false;
		is_empty = false;
       
        current_string = buf;
        parse_string( current_string, keyword, variable, value, is_dollar, is_empty, line );

		cout << "current_string: " << current_string << endl;
		cout << "keyword: " << keyword << endl;
		cout << "variable: " << variable << endl;
		cout << "value: " << value << endl;
		cout << "is_dollar: " << is_dollar << endl;
		cout << "is_empty: " << is_empty << endl << endl;
		
		// all assignments should be enclosed in some group
		if ( is_dollar && current_group == "" ) 
		{
			current_group = keyword;
			keyword = "";
		} 
		
		if ( &variable && current_group == "" ) 
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

		if ( current_group == "$grid" )
		{
			analyse_grid_group_line( variable, value, line, grid ); 
		}
		else 
		{
			string line_number_string = std::to_string( line );
			throw std::invalid_argument( "Unexpected group name: " + current_group + " on line " + line_number_string );
		}
		
		line++;
    }

}

void FileReader::analyse_grid_group_line( string& variable, string& value, int& line, GridParameters &grid )
{
	cout << "Inside parser of grid group" << endl;
}

void FileReader::parse_string( string curr_str, string& keyword, string& variable, string& value, bool& is_dollar, bool& is_empty, int& line ) 
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

		// isolating keyword 
		start_keyword = curr_str.find_first_not_of( space_symbols );
        //if (start_col != col)
	   	//{
            //throw std::invalid_argument("Non space symbols before $ occured.", line);
        //}
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
	cout << "File destructor" << endl;
}
