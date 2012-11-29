#include <stdlib.h>
#include <stdio.h>
#include <math.h>
int main(int argc, char* argv[])
{
	if(argc != 6)
	{
		printf("Usage is %s limit_type limit input_folder input_file output_folder\n", argv[0]);
		return 1;
	}
	int limit_type = atoi(argv[1]);
	float limit = atof(argv[2]);
	char *input_folder = argv[3];
	char *input_file = argv[4];
	char *output_folder = argv[5];

	// limit >= 1 means that you limiting the work done by the number of iterations
	// limit < 1 means that you are limiting the work done by the delta change	
	reconstruct(input_folder, input_file, output_folder, limit_type, &limit);
	return 0;
}
