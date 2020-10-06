/* Copyright 2009, 2010 Brendan Tauras */

/* ssd_ram.cpp is part of FlashSim. */

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

/* Ram class
 *
 * Brendan Tauras 2009-06-03
 *
 * This is a basic implementation that only provides delay updates to events
 * based on a delay value multiplied by the size (number of pages) needed to
 * be read or written.
 */

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ssd.h"

using namespace ssd;

 
Ram::Ram(double read_delay, double write_delay, Ssd &parent):
	read_delay(read_delay),
	write_delay(write_delay),
	ssd(parent)
{
	if(read_delay <= 0)
	{
		fprintf(stderr, "RAM: %s: constructor received negative read delay value\n\tsetting read delay to 0.0\n", __func__);
		read_delay = 0.0;
	}
	if(write_delay <= 0)
	{
		fprintf(stderr, "RAM: %s: constructor received negative write delay value\n\tsetting write delay to 0.0\n", __func__);
		write_delay = 0.0;
	}
	usedBlock = 0;
	current_ecc_page = 0;
	blockTable = (Ram_block *)malloc(sizeof(Ram_block) * (NUM_OF_CACHE_DATA_BLOCK + NUM_OF_CACHE_ECC_BLOCK));
	for(int i = 0; i < (NUM_OF_CACHE_DATA_BLOCK + NUM_OF_CACHE_ECC_BLOCK); i++)//DATA CACHE
	{
		(void) new (&blockTable[i]) Ram_block();
		blockTable[i].set_block_type(CACHE_DATA);
	}
	dram = (void *)malloc( PAGE_SIZE*BLOCK_SIZE*(NUM_OF_CACHE_DATA_BLOCK + NUM_OF_CACHE_ECC_BLOCK) );//block_size*page_size*
	memset(dram,NULL,PAGE_SIZE*BLOCK_SIZE*(NUM_OF_CACHE_DATA_BLOCK + NUM_OF_CACHE_ECC_BLOCK));
	return;
}
Ram::~Ram(void)
{
	free(dram);
	for(int i = 0; i < NUM_OF_CACHE_DATA_BLOCK + NUM_OF_CACHE_ECC_BLOCK; i++)
		blockTable[i].~Ram_block();
	return;
}
int Ram::choose_cache_block(enum ram_type type, int cache_block_index, int cache_block_tag, Event &event, FtlParent& ftl)
{
	bool ecc_exist = false;
	if (type == RAM_WRITE)
	{
		for (int i = 0; i < NUM_OF_CACHE_DATA_BLOCK; ++i)//cache 내부에 요청된 event에 대한 block이 존재하는 경우.
		{
			if (cache_block_index == blockTable[i].get_block_lpn_pointer())
			{
				blockTable[i].pageEvent[cache_block_tag] = event;
				return i;
			}
		}
		if (usedBlock == NUM_OF_CACHE_DATA_BLOCK)//캐시 내부의 모든 데이터블락이 사용중인경우
		{
			int victim_cache_block = cache_flush(ftl);
			blockTable[victim_cache_block].pageEvent[cache_block_tag] = event;
			return victim_cache_block;
		}
		else
		{
			for (int i = 0; i < NUM_OF_CACHE_DATA_BLOCK; ++i)
			{
				if (blockTable[i].get_block_valid() == false)
				{
					blockTable[i].pageEvent[cache_block_tag] = event;
					blockTable[i].set_block_valid(true);
					blockTable[i].set_block_lpn_pointer(cache_block_index);
					usedBlock++;
					return i;
				}
			}
		}
	}
	else if (type == RAM_READ)
	{
		for(int i = 0; i < NUM_OF_CACHE_DATA_BLOCK; ++i)
		{
			if (cache_block_index == blockTable[i].get_block_lpn_pointer())
			{
				if(blockTable[i].get_valid(cache_block_tag) == true)
				{
					return i;
				}
			}
		}
		bool isExist = false;
		int eccNum[2];
		for(int i = NUM_OF_CACHE_DATA_BLOCK; i<NUM_OF_CACHE_DATA_BLOCK + NUM_OF_CACHE_ECC_BLOCK;++i)
		{
			for(int j=0;j<RAM_BLOCK_SIZE;++j)
			{
				for(int k = 0; k<2;++k)
				{
					if((int)event.get_logical_address() == blockTable[i].pageEvent[j].get_eccLpn(k))
					{
						isExist == true;
						eccNum[0] = j;//cache ecc page num
						eccNum[1] = k;//ecc num
						EccHit++;
						ecc_exist = true;
					}
				}
			}
		}
		if(ssd.dltodp[event.get_logical_address()]==-1)
		{
			return 999;
		}
		else
		{
			if(ecc_exist == true)
			{
				return -9999;
			}
			else{
				return -999;
			}
		}
		
	}
}
int Ram::cache_flush(FtlParent& ftl)
{
	int victim = 0;
	int max_valid = 0;
	printf("cache flush\n");
	for (int i = 0; i < NUM_OF_CACHE_DATA_BLOCK; ++i)
	{	
		if (blockTable[i].get_num_valid() > max_valid)
		{
			victim = i;
			max_valid = blockTable[i].get_num_valid();
		}
	}	
	for (int i = 0; i < RAM_BLOCK_SIZE; ++i)
	{
		if(blockTable[victim].get_valid(i) == true)
		{
			//printf("%d\n", blockTable[victim].pageEvent[i].get_logical_address());
			ftl.write(blockTable[victim].pageEvent[i], OPDATA);
			blockTable[victim].set_valid(i, false);
			blockTable[victim].add_num_valid(-1);
		}
	}
	print();
	usedBlock--;
	blockTable[victim].clear_block();
	return victim;
}
enum status Ram::read(Event &event, FtlParent &ftl)
{
	assert(read_delay >= 0.0);
	if(event.get_event_type()!=READ)
	{
		return SUCCESS;
	}
	(void) event.incr_time_taken(read_delay * event.get_size());
	printf("Ram read\n");
	print();
	int event_lpn = (int)event.get_logical_address();
	int cache_block_index = event_lpn / RAM_BLOCK_SIZE;
	int cache_block_tag = event_lpn % RAM_BLOCK_SIZE;

	int current_cache_block = choose_cache_block(RAM_READ, cache_block_index, cache_block_tag, event, ftl);//현재 요청된 이벤트에 대한 캐쉬 내부 블락 받아옴.
	if (current_cache_block == 9999)//칩 내부에 요청된 이벤트에대한 데이터가 존재하지 않을 경우(SSD내부에 없는 데이터를 읽으려고 한 경우)
	{
		printf("miss\n");
		Miss++;
		return SUCCESS;
	}
	else if (current_cache_block == -9999)
	{
		printf("miss\n");
		Miss++;
		printf("cache write\n");
		event.set_payload(global_buffer);
		ftl.read(event,ECCEXIST);
		write(event, ftl); 
		return SUCCESS;
	}
	
	else if (current_cache_block == -999)
	{
		printf("miss\n");
		Miss++;
		printf("cache write\n");
		event.set_payload(global_buffer);
		ftl.read(event,OPDATA);
		write(event, ftl); 
		return SUCCESS;
	}
	else if(current_cache_block == 999)
	{
		Miss++;
		return SUCCESS;
	}
	else
	{
		printf("hit\n");
		Hit++;
		global_buffer = (char *)dram + ((PAGE_SIZE * RAM_BLOCK_SIZE * cache_block_index) + (PAGE_SIZE * cache_block_tag));
		return SUCCESS;
	}
}
enum status Ram::write(Event &event, FtlParent &ftl)
{

	assert(write_delay >= 0.0);
	if(event.get_event_type()!=WRITE)
	{
		return SUCCESS;
	}
	printf("ram write\n");
	(void) event.incr_time_taken(write_delay * event.get_size());
	
	int event_lpn = (int)event.get_logical_address();
	int cache_block_index = event_lpn / RAM_BLOCK_SIZE;
	int cache_block_tag = event_lpn % RAM_BLOCK_SIZE;
	
	int current_cache_block = choose_cache_block(RAM_WRITE, cache_block_index, cache_block_tag, event, ftl);//현재 요청된 이벤트가 쓰여질 캐쉬 내부 블락 결정.
	
	if(blockTable[current_cache_block].get_valid(cache_block_tag) == true)//현재 요청된 이벤트가 업데이트 요청일 경우.
	{
			//memcpy((char *)dram + ( (PAGE_SIZE * RAM_BLOCK_SIZE * current_cache_block) + (PAGE_SIZE * cache_block_tag) ), event.get_payload(), PAGE_SIZE);
			blockTable[current_cache_block].set_valid(cache_block_tag, true);
			blockTable[current_cache_block].add_num_valid(1);
	}
	else//업데이트 요청과 동일하게 동작.
	{
			//memcpy((char *)dram + ( (PAGE_SIZE * RAM_BLOCK_SIZE * current_cache_block) + (PAGE_SIZE * cache_block_tag) ), event.get_payload(), PAGE_SIZE);
			blockTable[current_cache_block].set_valid(cache_block_tag, true);
			blockTable[current_cache_block].add_num_valid(1);
	}
	return SUCCESS;
}
enum status Ram::ecc_write(Event &event, FtlParent &ftl)
{
	for(int i=0;i<4;++i)
	{
		printf("<RAM ecc_write>\nEcc no.%d lpn = %d\n", i, event.get_eccLpn(i));
	}
	int cache_block_index = 0;
	int block_num = current_ecc_page / RAM_BLOCK_SIZE;
	cache_block_index = block_num+NUM_OF_CACHE_DATA_BLOCK;
	int cache_tag[4];
	for(int i=0;i<4;++i)
	{
		blockTable[cache_block_index].set_block_type(CACHE_ECC);
		int ecc_page_num = current_ecc_page / RAM_BLOCK_SIZE;
		cache_tag[i] = event.get_eccLpn(i);
		blockTable[cache_block_index].pageEvent[current_ecc_page%RAM_BLOCK_SIZE].set_eccLpn(cache_tag[i], i);

	}
	printf("current_ecc_page = %d\n", current_ecc_page);
	blockTable[cache_block_index].set_valid(current_ecc_page%RAM_BLOCK_SIZE,true);
	//memcpy((char *)dram + ( (PAGE_SIZE * RAM_BLOCK_SIZE * 64) + (PAGE_SIZE * current_ecc_page) ), event.get_payload(), PAGE_SIZE);
	++current_ecc_page;
	if(current_ecc_page == BLOCK_SIZE * NUM_OF_CACHE_ECC_BLOCK)
	{
		current_ecc_page = 0;
	}
}
void Ram::print()
{
	printf("HIT = %d, EccHit = %d, Miss = %d\n", Hit, EccHit, Miss);
}
