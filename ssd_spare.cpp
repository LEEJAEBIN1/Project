#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdexcept>
#include <time.h>
#include <stdlib.h>
#include "ssd.h"

using namespace ssd;
Spare::Spare()
{
	ecc = (char *)calloc(PAGE_SIZE/4,sizeof(char));
	return;
}
Spare::~Spare()
{
	free(ecc);
	return;
}
void Spare::set_spare_lpn(uint data)
{
	this->lpn = data;
}
void Spare::set_spare_ppn(uint data)
{
	this->ppn = data;
}
void Spare::set_spare_ecc(char *data)
{

}
uint Spare::get_spare_lpn(void)
{
	return this->lpn;
}
uint Spare::get_spare_ppn(void)
{
	return this->ppn;
}
char* Spare::get_spare_ecc(void)
{
	return this->ecc;
}
void Spare::print_ecc(void)
{
	printf("parity = ");
	printf("%s\n", ecc);
}
