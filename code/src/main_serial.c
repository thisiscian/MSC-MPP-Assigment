#include <stdlib.h>
#include <stdio.h>
#include <math.h>
int main(int argc, char* argv[])
{
	if(argc != 5)
	{
		printf("Usage is %s limit input_folder input_file output_folder\n", argv[0]);
		return 1;
	}
	float limit = atof(argv[1]);
	char *input_folder = argv[2];
	char *input_file = argv[3];
	char *output_folder = argv[4];

	// limit >= 1 means that you limiting the work done by the number of iterations
	// limit < 1 means that you are limiting the work done by the delta change	
	reconstruct(input_folder, input_file, output_folder,  &limit);
	return 0;
}
