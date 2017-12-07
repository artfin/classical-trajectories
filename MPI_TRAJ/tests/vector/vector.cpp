#include <iostream>
#include <vector>
#include <stdlib.h>

using namespace std;

int main()
{
	vector<int> v1;
	v1.resize( 100 );

	for ( size_t i = 0; i < 10; i++ )
	{
		v1[i] = rand();
	}

	cout << "v1.size(): " << v1.size() << endl;
	cout << "v1.capacity(): " << v1.capacity() << endl;

	return 0;
}
