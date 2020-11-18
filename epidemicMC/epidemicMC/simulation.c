#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector_mtx.h"
#include "write.h"

#define RS_SCALE (1.0 / (1.0 + RAND_MAX))

double drand(void) 
{
	double d;
	do {
		d = (((rand() * RS_SCALE) + rand()) * RS_SCALE + rand()) * RS_SCALE;
	} while (d >= 1);
	return d;
}

double zero_to_one(void)
{
	return ((double)rand()) / RAND_MAX;
}

int is_already_infected(int x, int y, int* Yx, int* Yy, int nmax)
{
	int n_inds = 0;
	int retval = 0;
	int* inds = int_vector_malloc(nmax);

	for (int i = 0; i < nmax; i++)
	{
		if (Yx[i] == x)
		{
			inds[n_inds] = i;
			n_inds++;
		}
	}

	for (int i = 0; i < n_inds; i++)
	{
		if (Yy[inds[i]] == y)
		{
			retval = 1;
			break;
		}
	}

	free(inds);
	return retval;
}

double **MC(double p, int N, int tmax, int Lx, int Ly)
{
	double **Z_avg = mtx_malloc(tmax, Lx);
	double **Z = mtx_malloc(Lx, Ly);
	int *Yx = int_vector_malloc((int)(Lx*Ly)); // Co-ords of infected sites
	int *Yy = int_vector_malloc((int)(Lx*Ly)); // There are at most Lx*Ly infected sites
	int *Yx_next = int_vector_malloc((int)(Lx*Ly)); // Co-ords of infected sites at t+1
	int *Yy_next = int_vector_malloc((int)(Lx*Ly));
	int n_infected, next_n_infected;
	int do_infect;
	double roll;
	int site_x, site_y;
	int skip;
	int time_left;

	// initialize the array we save to
	for (int i = 0; i < tmax; i++)
	{
		for (int j = 0; j < Lx; j++)
		{
			Z_avg[i][j] = 0;
		}
	}

	for (int n = 0; n < N; n++)
	{
		printf("Realization: %d\n", n);

		// initialize

		for (int i = 0; i < Ly; i++)
		{
			Yx[i] = 0;
			Yy[i] = i;
		}
		n_infected = Ly;

		for (int i = 0; i < Lx; i++)
		{
			for (int j = 0; j < Ly; j++)
			{
				Z[i][j] = 0;
			}
		}

		time_left = 0;

		// iterate
		for (int t = 0; t < tmax; t++)
		{
			next_n_infected = 0;
			for (int site = 0; site < n_infected; site++)
			{
				
				// infect to the right
				site_y = Yy[site] + 1;
				// enforce bdy conditions
				if (site_y >= Ly)
					site_y = 0;
				// if not already infected
				if (is_already_infected(Yx[site], site_y, Yx, Yy, n_infected) + is_already_infected(Yx[site], site_y, Yx_next, Yy_next, next_n_infected) == 0)
				{
					// if not immune and we roll under p
					if ((Z[Yx[site]][site_y] == 0) & (drand() <= p))
					{
						Yx_next[next_n_infected] = Yx[site];
						Yy_next[next_n_infected] = site_y;
						// tick up the number of newly infected sites
						next_n_infected++;
					}
				}

				// infect to the left
				site_y = Yy[site] - 1;
				// enforce bdy conditions
				if (site_y < 0)
					site_y = Ly-1;
				// if not already infected
				if (is_already_infected(Yx[site], site_y, Yx, Yy, n_infected) + is_already_infected(Yx[site], site_y, Yx_next, Yy_next, next_n_infected) == 0)
				{
					// if not immune and we roll under p
					if ((Z[Yx[site]][site_y] == 0) & (drand() <= p))
					{
						Yx_next[next_n_infected] = Yx[site];
						Yy_next[next_n_infected] = site_y;
						// tick up the number of newly infected sites
						next_n_infected++;
					}
				}
				
				// infect to the bottom
				skip = 0; // for if we're looking to infect out of bounds
				site_x = Yx[site] + 1;
				// enforce bdy conditions
				if (site_x >= Lx)
					skip = 1;
				// if not already infected and in bounds
				if ((is_already_infected(site_x, Yy[site], Yx, Yy, n_infected) + is_already_infected(site_x, Yy[site], Yx_next, Yy_next, next_n_infected) == 0) & (skip == 0))
				{
					// if not immune and we roll under p
					if ((Z[site_x][Yy[site]] == 0) & (drand() <= p))
					{
						Yx_next[next_n_infected] = site_x;
						Yy_next[next_n_infected] = Yy[site];
						// tick up the number of newly infected sites
						next_n_infected++;
					}
				}

				// infect to the top
				skip = 0; // for if we're looking to infect out of bounds
				site_x = Yx[site] - 1;
				// enforce bdy conditions
				if (site_x < 0)
					skip = 1;
				// if not already infected and in bounds
				if ((is_already_infected(site_x, Yy[site], Yx, Yy, n_infected) + is_already_infected(site_x, Yy[site], Yx_next, Yy_next, next_n_infected) == 0) & (skip == 0))
				{
					// if not immune and we roll under p
					if ((Z[site_x][Yy[site]] == 0) & (drand() <= p))
					{
						Yx_next[next_n_infected] = site_x;
						Yy_next[next_n_infected] = Yy[site];
						// tick up the number of newly infected sites
						next_n_infected++;
					}
				}

			} // end site loop

			// add sum of every column to Zavg at time t
			for (int i = 0; i < Lx; i++)
			{
				for (int j = 0; j < Ly; j++)
				{
					Z_avg[t][i] += Z[i][j];
				}
			}


			// infected sites become immune
			for (int i = 0; i < n_infected; i++)
			{
				Z[Yx[i]][Yy[i]] = 1;
			}

			// update Y for next time (get rid of duplicates)
			n_infected = next_n_infected;
			for (int i = 0; i < n_infected; i++)
			{
				Yx[i] = Yx_next[i];
				Yy[i] = Yy_next[i];
			}

			// if no more infected, we know Z won't change
			// so we can fill out the rest of Zavg to save time
			if (n_infected == 0)
			{
				time_left = tmax - t;
				for (int tl = 1; tl < time_left; tl++)
				{
					for (int i = 0; i < Lx; i++)
					{
						for (int j = 0; j < Ly; j++)
						{
							Z_avg[t + tl][i] += Z[i][j];
						}
					} // end set Z_avg loop
				} // end time_left loop
				break;
			}

		} // end time loop

	} // end realisations loop

	// Divide by N to get average
	for (int i = 0; i < tmax; i++)
	{
		for (int j = 0; j < Lx; j++)
		{
			Z_avg[i][j] /= N * Ly;
		}
	}

	return Z_avg;
}