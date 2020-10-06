#include <new>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <queue>
#include "ssd.h"
using namespace ssd;

Ecc_block_page::Ecc_block_page()
{
}
Ecc_block_page::~Ecc_block_page()
{
}
int Ecc_block_page::get_ecc_count()
{
	return ecc_count;
}
void Ecc_block_page::set_ecc_count(int count)
{
	ecc_count=0;
}
void Ecc_block_page::add_ecc_count()
{
	ecc_count += 1;
}
int Ecc_block_page::get_ecc_ppn(int index)
{
	return ecc_ppn[index];
}
void Ecc_block_page::set_ecc_ppn(int index, int ppn)
{
	ecc_ppn[index] = ppn;
}