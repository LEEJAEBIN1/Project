#include <new>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <queue>
#include "ssd.h"

using namespace ssd;


Ram_block::Ram_block(void)
{
	pageEvent = (Event *)malloc(sizeof(Event) * 64);
	block_lpn_pointer = -1;
	used_block_number = 0;

	block_valid = false;	
	start_address = 0;
	num_valid = 0;
	num_dirty = 0;	
	for(int i=0;i<(RAM_BLOCK_SIZE);++i)
	{
		valid[i] = dirty[i] = false;
	}
	for(int i = 0; i < (RAM_BLOCK_SIZE); i++)
	{
		(void) new (&pageEvent[i]) Event();
	}

	PPN = 0;
}
Ram_block::~Ram_block(void)
{
	free(pageEvent);
	for(int i = 0; i < (RAM_BLOCK_SIZE); i++)
		pageEvent[i].~Event();
	return;
}
void Ram_block::clear_block()
{
	this->block_valid = false;
	this->start_address = 0;
	this->num_valid = 0;
	this->num_dirty = 0;
	for(int i=0;i<RAM_BLOCK_SIZE;++i)
	{
		this->valid[i] = false;
		this->dirty[i] = false;
	}
	this->PPN = 0;
}
bool Ram_block::get_block_valid()
{
	return this->block_valid;
}
int Ram_block::get_start_address()
{	
	return start_address;
}
char Ram_block::get_num_valid()
{
	return num_valid;
}
char Ram_block::get_num_dirty()
{
	return num_dirty;
}
char Ram_block::get_ppn()
{
	return PPN;
}
bool Ram_block::get_valid(char index)
{
	return valid[index];
}
bool Ram_block::get_dirty(char index)
{
	return dirty[index];
}
void Ram_block::set_block_valid(bool valid)
{
	this->block_valid = valid;
}
void Ram_block::set_start_address(int address)
{
	this->start_address = address;
}
void Ram_block::set_ppn(char ppn)
{
	this->PPN = ppn;
}
void Ram_block::set_valid(char index, bool valid)
{
	this->valid[index] = valid;
}
void Ram_block::set_dirty(char index, bool dirty)
{
	this->dirty[index] = dirty;
}
void Ram_block::add_num_valid(char num)
{
	this->num_valid = this->num_valid + num;
}
void Ram_block::add_num_dirty(char num)
{
	this->num_dirty = this->num_dirty + num;
}
void Ram_block::set_tag(int index, int tag)
{
	this->tag[index] = tag;
}
int Ram_block::get_tag(int index)
{
	return this->tag[index];
}
enum ram_type Ram_block::get_block_type()
{
	return this->state;
}
void Ram_block::set_block_type(enum ram_type state)
{
	this->state = state;
}
void Ram_block::set_block_lpn_pointer(int pointer)
{
	this->block_lpn_pointer = pointer;
}
int Ram_block::get_block_lpn_pointer()
{
	return this->block_lpn_pointer;
}
