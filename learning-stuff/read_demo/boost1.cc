#include <boost/spirit/include/qi.hpp>
#include <string>
#include <iostream>

using namespace boost::spirit;
using namespace std;

int main() 
{
	string s;
	getline(cin, s);
	// passing an iterator of a string which is read from std::cin
	// also it is not passed directly, but through interim variable 'it', because 'parse' may modify the iterator

	auto it = s.begin(); // auto will become the needed type (based on some template magic)
	bool match = qi::parse(it, s.end(), ascii::digit); // expects two iterators -- it and s.end -- and a parser: ascii::digit
	// this parser tests whether a character is a digit between 0 and 9
	// tests only the first digit. so constructions like '1sdforwngm' will pass, but words like 'digit' will not.
	
	if ( it != s.end() )
	{
		cout << string{it, s.end()} << endl;
	}

	return 0;
}
