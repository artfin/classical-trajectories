#include <iostream>
#include <functional>

using namespace std;

double f( double x )
{
	return x * x;
}

class Tester
{
public:
	Tester( function<double(double)> f ) : f(f) { } 

	function<double(double)> f;
};

int main()
{
	Tester tester( f );

	cout << "tester.f(2): " << tester.f( 2 ) << endl;

	return 0;
}
