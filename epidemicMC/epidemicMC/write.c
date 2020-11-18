#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void WriteResults(char *outfile, double **data, int nmax, int mmax)
{
	FILE *output;

	output = fopen(outfile, "w");

	for (int i = 0; i < nmax; i++)
	{
		for (int j = 0; j < mmax-1; j++)
		{
			fprintf(output, "%lf, ", data[i][j]);
		}
		fprintf(output, "%lf\n", data[i][mmax - 1]);
	}

	fclose(output);

	return;
}