#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

int pow10( int n )
{
	if ( n > 0 )
		return 10 * pow10(n - 1);
	else
		return 1;
}

bool read_int_unlocked(int &out)
{	
	int c = getchar_unlocked();

	int x = 0;
	int n = 1;

	int neg = 0;
	bool isdot = false;
	
	cout << "c: " << c << endl;

	for ( ; !('0' <= c && c <= '9') && c != '-' && c != '.'; c = getchar_unlocked() )
	{
		if ( c == EOF )
		{
			return false;
		}
	}

	if ( c == '.' )
	{
		isdot = true;
		cout << "Found dot!" << endl;
		c = getchar_unlocked();
	}

	if ( c == '-' )
	{
		neg = 1;
		c = getchar_unlocked();
	}

	if ( c == EOF )
	{
		return false;
	}

	for ( ; '0' <= c && c <= '9'; c = getchar_unlocked() ) 
	{
		cout << "inside second cycle" << endl;
		if ( isdot == false ) 
		{
			x = x * 10 + c - '0';
			cout << "isdot is false; x = " << x << endl;
		} else {
			x = x + c / pow10(n);
			cout << "isdot is true; x = " << x << endl;
			n++;
		}
	}

	out = neg ? - x : x;
	return true;
}

int max_getchar_unlocked()
{
	int x;
	int _max = -1;
	
	while ( read_int_unlocked(x) )
	{
		cout << "returned x: " << x << endl;
		_max = max(x, _max);
	}

	return _max;
}

int main()
{
	cout << max_getchar_unlocked() << endl;
	return 0;
}
