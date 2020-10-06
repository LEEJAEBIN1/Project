#include <new>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <queue>
#include "ssd.h"

using namespace ssd;

Ftl_block_manager::Ftl_block_manager(FtlImpl_Fast* ftl):
	ftl(ftl)
{
	num_valid = new int[SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE];
	num_invalid = new int[SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE];
	num_free = new int[SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE];
	current_page = new int[SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE];
	block_state = new enum block_state[SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE];
	block_type = new enum block_type[SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE];
	used_num = new int[SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE];
	for (int i = 0; i < SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE; ++i)
	{
		num_valid[i] = 0;
		num_invalid[i] = 0;
		num_free[i] = 0;
		current_page[i] = 0;
		block_state[i] = FREE;
		block_type[i] = NOTHING;
		used_num[i] = 0;
	}
	free_block_num = SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE;
	return;
}
Ftl_block_manager::~Ftl_block_manager()
{
	delete num_valid;
	delete num_invalid;
	delete num_free;
	delete current_page;
	delete block_state;
	delete block_type;
	delete used_num;
	num_valid = NULL;
	num_invalid = NULL;
	num_free = NULL;
	current_page = NULL;
	block_state = NULL;
	block_type = NULL;
	used_num = NULL;
	return;
}

void ssd::Ftl_block_manager::add_valid_num(int block_num, int add_num)
{
	//printf("%d block add valid num + %d\n",block_num, add_num);
	num_valid[block_num] += add_num;
}

int ssd::Ftl_block_manager::get_valid_num(int block_num)
{
	return num_valid[block_num];
}

void ssd::Ftl_block_manager::add_invalid_num(int block_num, int add_num)
{
	num_invalid[block_num] += add_num;
}

int ssd::Ftl_block_manager::get_invalid_num(int block_num)
{
	return num_invalid[block_num];
}

void ssd::Ftl_block_manager::add_free_num(int block_num, int add_num)
{
	num_free[block_num] += add_num;
}

int ssd::Ftl_block_manager::get_free_num(int block_num)
{
	return num_free[block_num];
}

void ssd::Ftl_block_manager::set_current_page(int block_num, int page_num)
{
	current_page[block_num] = page_num;
}

int ssd::Ftl_block_manager::get_current_page(int block_num)
{
	return current_page[block_num];
}

enum block_state ssd::Ftl_block_manager::get_block_state(int block_num)
{
	return block_state[block_num];
}

void ssd::Ftl_block_manager::set_block_state(int block_num, enum block_state bstate)
{
	block_state[block_num] = bstate;
}
enum block_type ssd::Ftl_block_manager::get_block_type(int block_num)
{
	return block_type[block_num];
}

void ssd::Ftl_block_manager::set_block_type(int block_num, enum block_type btype)
{
	block_type[block_num] = btype;
}
void ssd::Ftl_block_manager::clean_block(int block_num)
{
	init_block(block_num);
}
void ssd::Ftl_block_manager::init_block(int block_num)
{
	num_valid[block_num] = 0;
	num_invalid[block_num] = 0;
	num_free[block_num] = 0;
	current_page[block_num] = 0;
	used_num[block_num] = 0;
	block_state[block_num] = FREE;
}

int ssd::Ftl_block_manager::get_free_block_num()
{
	return free_block_num;
}

void ssd::Ftl_block_manager::set_free_block_num(int num)
{
	free_block_num += num;
}

int ssd::Ftl_block_manager::get_valid_page_num()
{
	int num = 0;
	for (int i = 0; i < SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE; ++i)
	{
		num += get_valid_num(i);
	}
	return num;
}
void ssd::Ftl_block_manager::set_used_page_num(int num)
{
	used_page_num += num;
}

int ssd::Ftl_block_manager::get_used_page_num()
{
	int num = 0;
	for(int i = 0;i<SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE; ++i)
	{
		num += get_used_num(i);
	}
	return num;
}

void ssd::Ftl_block_manager::add_used_num(int block_num,int add_num)
{
	used_num[block_num] += add_num;
}

int ssd::Ftl_block_manager::get_used_num(int block_num)
{
		return used_num[block_num];
}
