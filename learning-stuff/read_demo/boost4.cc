#include <boost/spirit/include/qi.hpp>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>

using namespace boost::spirit;
using namespace std;

int main()
{
	string s;
	getline(cin, s);

	auto it = s.begin();
	vector<int> v;

	if ( qi::phrase_parse(it, s.end(), qi::int_ % ',', ascii::space, v) )
	{
		ostream_iterator<int> out{cout, ";"};
		copy(v.begin(), v.end(), out);
	}
	
	return 0;
}

