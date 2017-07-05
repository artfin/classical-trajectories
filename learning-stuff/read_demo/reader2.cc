#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>

#include <iostream>
#include <string>
#include <vector>

namespace client
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;

    ///////////////////////////////////////////////////////////////////////////
    //  Our number list parser
    ///////////////////////////////////////////////////////////////////////////
    //[tutorial_numlist1
    template <typename Iterator>
    bool parse_numbers(Iterator first, Iterator last)
    {
        using qi::double_;
        using qi::phrase_parse;
        using ascii::space;

        bool r = phrase_parse(
            first,                          /*< start iterator >*/
            last,                           /*< end iterator >*/
            double_ >> *(',' >> double_),   /*< the parser >*/
            space                           /*< the skip-parser >*/
        );

        if (first != last) // fail if we did not get a full match
            return false;
        return r;
    }
    //]
}

////////////////////////////////////////////////////////////////////////////
//  Main program
////////////////////////////////////////////////////////////////////////////
int main()
{

	std::ifstream ifs("data.txt");
	ifs >> std::noskipws;
	
	std::string str = "12.0";
	client::parse_numbers(str.begin(), str.end());
	std::cout << "------------------------\n";
	std::cout << "Parsing succeeded\n";
	std::cout << str << " Parses OK: " << std::endl;

	boost::spirit::istream_iterator f(ifs), l;

	bool ok = phrase_parse(
		f, 
		l,
		double_ >> *(',' >> double_),
		space
	);

    // std::string str;
    // while (getline(std::cin, str))
    // {
    //    if (str.empty() || str[0] == 'q' || str[0] == 'Q')
    //        break;

    //    if (client::parse_numbers(str.begin(), str.end()))
    //    {
    //        std::cout << "-------------------------\n";
    //        std::cout << "Parsing succeeded\n";
    //        std::cout << str << " Parses OK: " << std::endl;
    //    }
    //    else
    //    {
    //        std::cout << "-------------------------\n";
    //        std::cout << "Parsing failed\n";
    //        std::cout << "-------------------------\n";
    //    }
    //}

    std::cout << "Bye... :-) \n\n";
    return 0;
}

