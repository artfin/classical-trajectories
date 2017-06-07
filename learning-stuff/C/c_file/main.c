#include <stdio.h>

int main() {
	FILE *ofp;
	
	ofp = fopen("output.dat", "w");

	if (ofp == NULL) {
		fprintf(stderr, "Can't open output file!\n");
		return 1;
	}
	else {
		fprintf(ofp, "HEY!");
	}

	return 0;
}
