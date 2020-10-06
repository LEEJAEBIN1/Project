#include <new>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <queue>
#include "ssd.h"

using namespace ssd;

EccBlock::EccBlock(int size):
	size(size),
	data((EccPage *) malloc(size * sizeof(EccPage)))
{
	for(int i = 0; i < size; i++)
		(void) new (&data[i]) EccPage(PAGE_READ_DELAY, PAGE_WRITE_DELAY);

}
EccBlock::~EccBlock(void)
{
	for(int i = 0; i < size; ++i)
	{
		data[i].~EccPage();
	}
	free(data);
}
enum status EccBlock::read(int EccPageAddress)
{
	return SUCCESS;
}
enum status EccBlock::write(int EccPageAddress)
{
	return SUCCESS;
}
