#include <iostream>

#include <vector>
#include <boost/spirit/include/qi.hpp>

using namespace std;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

struct parser : qi::grammar<string::const_iterator, vector<double>(), ascii::space_type>
{
	parser() : parser::base_type( v )
	{
		v %= +(qi::double_);
	}

	qi::rule<string::const_iterator, vector<double>(), ascii::space_type> v;
};

int main()
{
	string const x( "1 2 3 4" );
	string::const_iterator b = x.begin();
	string::const_iterator e = x.end();

	parser p;
	bool const r = qi::phrase_parse( b, e, p, ascii::space );
	cout << ( (b == e && r) ? "PASSED" : "FAILED" ) << endl;
}
