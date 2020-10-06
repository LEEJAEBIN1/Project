#include <new>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <queue>
#include "ssd.h"

using namespace ssd;

EccPage::EccPage(double read_delay, double write_delay):
	read_delay(read_delay),
	write_delay(write_delay)
{
}
EccPage::~EccPage(void)
{	
}
enum status EccPage::read()
{
	return SUCCESS;
}
enum status EccPage::write()
{
	return SUCCESS;
}
