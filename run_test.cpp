/* Copyright 2009, 2010 Brendan Tauras */

/* run_test.cpp is part of FlashSim. */

/* FlashSim is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version. */

/* FlashSim is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License
 * along with FlashSim.  If not, see <http://www.gnu.org/licenses/>. */

/****************************************************************************/

/* Basic test driver
 * Brendan Tauras 2009-11-02
 *
 * driver to create and run a very basic test of writes then reads */

#include "ssd.h"
#include <stdlib.h>
#include <time.h>
#define SIZE 1500

using namespace ssd;
int read_request=0;
int write_request = 0;
int main()
{
	load_config();
	print_config(NULL);
	srand(2);
	Ssd *ssd = new Ssd();

	double result;
	/*for (int i = 0; i < 3500; i++)
	{
		result = ssd -> event_arrive(WRITE, i, 1, 0);
		printf("Write time: %.20lf\n", result);
	}*/
	for (int i = 0; i < 30000; i++)
	{
		printf("\n------HOST WRITE REQUEST------\n", write_request);
		/* event_arrive(event_type, logical_address, size, start_time) */
		int j = i%12000;
		result = ssd -> event_arrive(WRITE, j, 1, 0);
		write_request++;
		//printf("Write time: %.20lf\n", result);
	}
	for (int i = 0; i<15000; ++i)
	{
		/*printf("\n------HOST WRITE REQUEST------\n", write_request);
		int lpn = rand()%9996;
		for(int j = lpn; j<lpn+4; ++j)
		{
			result = ssd -> event_arrive(WRITE, j, 1, 0);
			write_request++;
		}

		lpn = rand()%9996;
		for(int k = lpn; k<lpn+4; ++k)
		{
			printf("\n------%d HOST READ REQUEST------\n", read_request);
			result = ssd -> event_arrive(READ, k, 1, 0);
			read_request++;
		}*/
		printf("\n------HOST WRITE REQUEST------\n", write_request);
		int lpn = rand()%12000;
		result = ssd -> event_arrive(WRITE, lpn, 1, 0);
		write_request++;
		lpn = rand()%12000;
		printf("\n------%d HOST READ REQUEST------\n", read_request);
		result = ssd -> event_arrive(READ, lpn, 1, 0);
		read_request++;
	}
	delete ssd;
	ssd = NULL;
	return 0;
}