#include <stdlib.h>
#include <stdio.h>
#include <math.h>
int main()
{
	float iterations = 100;
	float accuracy = 0.1;
	//char infilename[256] = "../input_files/edge192x128.pgm";
	//char infilename[256] = "../input_files/edge256x192.pgm";
	//char infilename[256] = "../input_files/edge512x384.pgm";
	//char infilename[256] = "../input_files/edge768x768.pgm";
	char infilename[256] = "../input_files/edge1024x1408.pgm";
	
	// type = 0 means that you limiting the work done by the number of iterations
	// type = 1 means that you are limiting the work done by the delta change	
	p_reconstruct(infilename, 0, &iterations);
	//p_reconstruct(infilename, 1, &accuracy);
	//reconstruct(infilename, iterations);
	return 0;
}
