/* Copyright 2009, 2010 Brendan Tauras */

/* ssd_page.cpp is part of FlashSim. */

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

/* Page class
 * Brendan Tauras 2009-04-06
 *
 * The page is the lowest level data storage unit that is the size unit of
 * requests (events).  Pages maintain their state as events modify them. */

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdexcept>

#include "ssd.h"

namespace ssd {
	/*
	 * Buffer used for accessing data pages.
	 */
	void *global_buffer;
	int data_read_count = 0;
	int ecc_read_count = 0;
	int gc_read_count = 0;
	int count4_20 = 0;
	int count8_24 = 0;
	int count12_28 = 0;
	int count_decodingError = 0;
	int count_readgenError = 0;
	int count_4_20_decoding_error = 0;
}

using namespace ssd;

Page::Page(const Block &parent, double read_delay, double write_delay):
	state(EMPTY),
	parent(parent),
	read_delay(read_delay),
	write_delay(write_delay)
{
	if(read_delay < 0.0){
		fprintf(stderr, "Page warning: %s: constructor received negative read delay value\n\tsetting read delay to 0.0\n", __func__);
		this -> read_delay = 0.0;
	}

	if(write_delay < 0.0){
		fprintf(stderr, "Page warning: %s: constructor received negative write delay value\n\tsetting write delay to 0.0\n", __func__);
		this -> write_delay = 0.0;
	}
	eccState = new enum page_state[4];
	for (int i = 0; i < 4; ++i)
	{
		eccState[i] = EMPTY;
	}
	spare = new Spare();
	spare->set_spare_lpn(-1);
	spare->set_spare_ppn(-1);
	return;
}

Page::~Page(void)
{
	delete eccState;
	delete spare;
	return;
}

enum status Page::_read(Event &event)
{
	printf("4 20 = %d, 8 24 = %d, 12 28 = %d, 4_20fail %d, fail = %d\n",count4_20,count8_24,count12_28,count_4_20_decoding_error, count_decodingError);
	assert(read_delay >= 0.0);
	int returnvalue = 0;
	if(event.get_event_type() == ECCREAD_STEP1)
	{

		int eccnum = event.get_eccnum();
		/*global_buffer = (char*)page_data + (event.get_address().get_linear_address() * PAGE_SIZE);
		char *Buff;
		Buff = (char *)calloc(5120,sizeof(char));

		char *parity;
		parity = (char *)malloc(sizeof(char)*NUMBER_OF_PARITY);
		char *message;
		message = (char *)malloc(sizeof(char)*NUMBER_OF_MESSAGE);

		FILE *readfile;
		readfile = fopen("read.rec","r");
		FILE *readfile_2;
		readfile_2 = fopen("read_2.rec","w");

		for(int i = 0;i<PAGE_SIZE/NUMBER_OF_MESSAGE; ++i)
		{

			//spare ecc
			fgets(parity, NUMBER_OF_PARITY+1, readfile);
			strcat(Buff,parity);
			for(int j=0;j<NUMBER_OF_PARITY;++j)
			{
				fprintf(readfile_2,"%c",Buff[((i*(NUMBER_OF_PARITY+NUMBER_OF_MESSAGE))+j)]);
			}

			//exten ecc
			for(int k=0;k<NUMBER_OF_PARITY;++k)
			{
				fprintf(readfile_2,"%c",((char *)global_buffer+(PAGE_SIZE/2*eccnum))[(i*NUMBER_OF_PARITY)+k]);
			}

			//data
			fgets(message,NUMBER_OF_MESSAGE+1,readfile);
			strcat(Buff,message);
			for(int l=0;l<NUMBER_OF_MESSAGE;++l)
			{
				fprintf(readfile_2,"%c",Buff[((i*(NUMBER_OF_PARITY+NUMBER_OF_MESSAGE))+l+NUMBER_OF_PARITY)] );
			}
		}
		free(Buff);
		free(parity);
		free(message);
		fclose(readfile);
		fclose(readfile_2);*/
		
		int value1;
		value1 = (rand()%1000000);
		if(value1 >= 0 && value1 <= ERROR_8_24)
		{
			returnvalue = 1;
		}
		else
		{
			returnvalue = 0;
		}
		
		if(returnvalue == 1)
		{
			printf("EXTEND DECODING SUCCESS\n");
			count8_24 += 1;
		
			return SUCCESS;
		}
		else if(returnvalue == 0)
		{
			count_decodingError += 1;
			return SUCCESS;
			printf("EXTEND DECODING FAIL\n");
			return FAILURE;
		}
		else
		{
			printf("READ GEN ERROR\n");
			return READFAIL;
		}
	}
	else if(event.get_event_type() == ECCREAD_STEP2)
	{
		int eccnum = event.get_eccnum();
		/*global_buffer = (char*)page_data + (event.get_address().get_linear_address() * PAGE_SIZE);
		char *Buff;
		Buff = (char *)calloc(6144,sizeof(char));

		char *parity;
		parity = (char *)malloc(sizeof(char)*NUMBER_OF_PARITY);
		char *parity_2;
		parity_2 = (char *)malloc(sizeof(char)*NUMBER_OF_PARITY);
		char *message;
		message = (char *)malloc(sizeof(char)*NUMBER_OF_MESSAGE);

		FILE *readfile_2;
		readfile_2 = fopen("read_2.rec","r");
		FILE *readfile_3;
		readfile_3 = fopen("read_3.rec","w");

		for(int i = 0;i<PAGE_SIZE/NUMBER_OF_MESSAGE; ++i)
		{

			//spare ecc
			fgets(parity, NUMBER_OF_PARITY+1, readfile_2);
			strcat(Buff,parity);
			for(int j=0;j<NUMBER_OF_PARITY;++j)
			{
				fprintf(readfile_3,"%c",Buff[((i*(NUMBER_OF_PARITY+NUMBER_OF_PARITY+NUMBER_OF_MESSAGE))+j)]);
			}
			//exten ecc
			fgets(parity_2, NUMBER_OF_PARITY+1, readfile_2);
			strcat(Buff,parity_2);
			for(int j=0;j<NUMBER_OF_PARITY;++j)
			{
				fprintf(readfile_3,"%c",Buff[((i*(NUMBER_OF_PARITY+NUMBER_OF_PARITY+NUMBER_OF_MESSAGE))+j+NUMBER_OF_PARITY)]);
			}

			//exten2 ecc
			for(int k=0;k<NUMBER_OF_PARITY;++k)
			{
				fprintf(readfile_3,"%c",((char *)global_buffer+((PAGE_SIZE/2)*eccnum)+(PAGE_SIZE/4))[(i*NUMBER_OF_PARITY)+k]);
			}

			//data
			fgets(message,NUMBER_OF_MESSAGE+1,readfile_2);
			strcat(Buff,message);
			for(int l=0;l<NUMBER_OF_MESSAGE;++l)
			{
				fprintf(readfile_3,"%c",Buff[((i*(NUMBER_OF_PARITY+NUMBER_OF_PARITY+NUMBER_OF_MESSAGE))+l+NUMBER_OF_PARITY+NUMBER_OF_PARITY)] );
			}
		}
		free(Buff);
		free(parity);
		free(parity_2);
		free(message);
		fclose(readfile_2);
		fclose(readfile_3);*/
		int value1;
		/*value1 = (rand()%1000000);
		if(value1 >= 0 && value1 <= 931396) // 40
		{
			returnvalue = 1;
		}
		else
		{
			returnvalue = 0;
		}*/	
		returnvalue = 0;
		if(returnvalue == 1)
		{
			printf("EXTEND_2 DECODING SUCCESS\n");
			count12_28 += 1;
			printf("4 20 = %d, 8 24 = %d, 12 28 = %d, 4_20fail %d, fail = %d\n",count4_20,count8_24,count12_28,count_4_20_decoding_error, count_decodingError);			return SUCCESS;
		}
		else if(returnvalue == 0)
		{
			printf("EXTEND_2 DECODING FAIL\n");
			count_decodingError += 1;
			printf("4 20 = %d, 8 24 = %d, 12 28 = %d, 4_20fail %d, fail = %d\n",count4_20,count8_24,count12_28,count_4_20_decoding_error, count_decodingError);			return FAILURE;
		}

	}
	else if(event.get_event_type() == GCREAD)
	{
		gc_read_count += 1;
		event.incr_time_taken(read_delay);
		if (!event.get_noop() && PAGE_ENABLE_DATA)
		{
			//global_buffer = (char*)page_data + event.get_address().get_linear_address() * PAGE_SIZE;
		}
			return SUCCESS;
	}
	else
	{
		/*FILE *readfile;
		readfile = fopen("read.rec","w");
		char *Buff;
		Buff = (char *)calloc(1024,sizeof(char));
		Buff= spare->get_spare_ecc();
		
		if (!event.get_noop() && PAGE_ENABLE_DATA)
		{
			global_buffer = (char*)page_data + event.get_address().get_linear_address() * PAGE_SIZE;
		}
		for(int i = 0 ; i < PAGE_SIZE/NUMBER_OF_MESSAGE; ++i)
		{
			for(int j = 0;j<NUMBER_OF_PARITY;++j)
			{
				fprintf(readfile,"%c",Buff[(i*NUMBER_OF_PARITY)+j]);
			}
			for(int k = 0;k<NUMBER_OF_MESSAGE;++k)
			{
				fprintf(readfile,"%c",((char *)global_buffer)[(i*NUMBER_OF_MESSAGE)+k]);
			}
		}
		fclose(readfile);*/

		event.incr_time_taken(read_delay);
		int value1;
		value1 = (rand()%1000000);
		if(value1 >= 0 && value1 <= ERROR_4_20)
		{
			returnvalue = 1;
		}
		else
		{
			returnvalue = 0;
		}
		if(returnvalue == 1)
		{
			count4_20 += 1;
			return SUCCESS;
		}
		else if(returnvalue == 0)
		{
			count_4_20_decoding_error += 1;
			printf("SPARE ECC DECODING ERROR\n");
			
			return FAILURE;
		}
		
	}
}

enum status Page::_write(Event &event)
{
	assert(write_delay >= 0.0);
	event.incr_time_taken(write_delay);
	if(event.get_event_type() == ECCWRITE)
	{
		//void *data = (char*)page_data + event.get_address().get_linear_address() * PAGE_SIZE;
		//memcpy (data, event.get_payload(), PAGE_SIZE);
		this->set_state(VALID);
		return SUCCESS;
	}
	else
	{
		if (PAGE_ENABLE_DATA && event.get_payload() != NULL && event.get_noop() == false)
		{
			void *data = (char*)page_data + event.get_address().get_linear_address() * PAGE_SIZE;
			//memcpy (data, event.get_payload(), PAGE_SIZE);
			//memcpy (spare->get_spare_ecc(), event.get_payload_parity(), PAGE_SIZE/4);
			this->set_state(VALID);
		}
	}
	if (event.get_noop() == false)
	{
		//assert(state == EMPTY);
		state = VALID;
	}
	return SUCCESS;
}
Spare* Page::get_spare()
{
	return spare;
}
const Block &Page::get_parent(void) const
{
	return parent;
}

enum page_state Page::get_state(void) const
{
	return state;
}
enum page_state Page::get_ecc_state(int eccnum)
{
	return eccState[eccnum];
}
void Page::set_state(enum page_state state)
{
	this -> state = state;
}
void Page::set_ecc_state(int eccnum, enum page_state state)
{
	this->eccState[eccnum] = state;
}
void Page::add_peCycle(int num)
{
	this->peCycle += num;
}
int Page::get_peCycle()
{
	return this->peCycle;
}
