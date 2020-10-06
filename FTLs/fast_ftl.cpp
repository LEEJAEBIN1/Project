#include <new>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <queue>
#include <iostream>
#include <signal.h>
#include "../ssd.h"

using namespace ssd;

FtlImpl_Fast::FtlImpl_Fast(Controller* controller, Ssd& ssd, Ram& ram) :
	controller(controller),
	ssd(ssd),
	FtlParent(controller, ssd),
	ram(ram)
{
	gcRatio = 0;

	numCells = SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE * BLOCK_SIZE;
	//************mapping table

	//**************************
	//*************block manager
	bm = new Ftl_block_manager(this);
	eccBuf = (void*)malloc(PAGE_SIZE);
	eccBuffer = (void*)malloc(PAGE_SIZE);
	memset(eccBuf, NULL, PAGE_SIZE);
	memset(eccBuffer, NULL, PAGE_SIZE);

	return;
}

FtlImpl_Fast::~FtlImpl_Fast()
{


	free(eccBuf);
	return;
}

enum status FtlImpl_Fast::read(Event& event, optype type)
{
	int pecycle = 0;
	print_statistic();
	if(type == ECCEXIST)
	{
		if(bm->get_block_type(ssd.dltodp[event.get_logical_address()]/BLOCK_SIZE) == ECC)
		{
			printf("CAN'T READ ECC PAGE FIRST\n");
			return FAILURE;
		}
		if (ssd.dltodp[event.get_logical_address()] == -1)
		{
			return FAILURE;
		}
		else
		{
			pecycle = ssd.get_tail((resolve_logical_address(ssd.dltodp[event.get_logical_address()])))->get_peCycle();
			event.set_address(resolve_logical_address(ssd.dltodp[event.get_logical_address()]));
			printf("read lpn = %d, ppn = %d\n",event.get_logical_address(), ssd.dltodp[event.get_logical_address()]);
			enum status read_state = ssd.read(event);
			read_count += 1;
			if(read_state == FAILURE)
			{
				enum status read_state_extend;
				printf("extend 1 step decoding\n");
				event.set_event_type(ECCREAD_STEP1);
				int ecc_page_ppn = ssd.dltoep[event.get_logical_address()] >> 2;
				if(ecc_page_ppn == -1)
				{
					return SUCCESS;
				}
				int ecc_num = ssd.dltoep[event.get_logical_address()] - (ecc_page_ppn << 2);
				event.set_eccnum(ecc_num);
				event.set_address(resolve_logical_address(ecc_page_ppn));
				read_state_extend = ssd.read(event);
				if(read_state_extend == FAILURE)
				{
					printf("extend 2 step decoding\n");
					event.set_event_type(ECCREAD_STEP2);
					return ssd.read(event);
				}
			}
			else if(read_state == READFAIL)
			{
				return SUCCESS;
			}
			else
			{
				return read_state;
			}
		}
	}
	if (type == OPDATA)
	{
		if(bm->get_block_type(ssd.dltodp[event.get_logical_address()]/BLOCK_SIZE) == ECC)
		{
			printf("CAN'T READ ECC PAGE FIRST\n");
			return FAILURE;
		}
		if (ssd.dltodp[event.get_logical_address()] == -1)
		{
			return FAILURE;
		}
		else
		{
			pecycle = ssd.get_tail((resolve_logical_address(ssd.dltodp[event.get_logical_address()])))->get_peCycle();
			event.set_address(resolve_logical_address(ssd.dltodp[event.get_logical_address()]));
			printf("read lpn = %d, ppn = %d\n",event.get_logical_address(), ssd.dltodp[event.get_logical_address()]);
			enum status read_state = ssd.read(event);
			read_count += 1;
			if(read_state == FAILURE)
			{
				enum status read_state_extend;
				printf("extend 1 step decoding\n");
				event.set_event_type(ECCREAD_STEP1);
				int ecc_page_ppn = ssd.dltoep[event.get_logical_address()] >> 2;
				if(ecc_page_ppn == -1)
				{
					return SUCCESS;
				}
				int ecc_num = ssd.dltoep[event.get_logical_address()] - (ecc_page_ppn << 2);
				event.set_eccnum(ecc_num);
				event.set_address(resolve_logical_address(ecc_page_ppn));
				read_state_extend = ssd.read(event);
				ecc_read_count += 1;
				if(read_state_extend == FAILURE)
				{
					printf("extend 2 step decoding\n");
					event.set_event_type(ECCREAD_STEP2);
					return ssd.read(event);
				}
			}
			else if(read_state == READFAIL)
			{
				return SUCCESS;
			}
			else
			{
				return read_state;
			}
		}
	}
	return SUCCESS;
}

enum status FtlImpl_Fast::write(Event& event, optype type)
{
	print_statistic();
	if (type == OPDATA)
	{
		int returnvalue;

		//returnvalue = Ldpc_encode(2);
		if(returnvalue == 3)
		{
			printf("read_pchk error\n");
			return SUCCESS;
		}
		double error_prob = 0;
		int pecycle = 0;
		printf("----------------------------------------------------------------------------\n");
		if (ssd.dltodp[event.get_logical_address()] != -1)
		{
			bm->add_invalid_num(ssd.dltodp[event.get_logical_address()] / BLOCK_SIZE, 1);
			ssd.get_tail(resolve_logical_address(ssd.dltodp[event.get_logical_address()]))->set_state(INVALID);
		}
		int update_ppn = insert_event(event, OPDATA);
		Address address = resolve_logical_address(update_ppn);
		event.set_address(address);
		ssd.dltodp[event.get_logical_address()] = update_ppn;
		//pecycle = ssd.get_tail((resolve_logical_address(ssd.dltodp[event.get_logical_address()])))->get_peCycle();

	
		error_prob = 0.001255;

		/*ldpc_transmit(error_prob);

		parity = (char *)calloc(PAGE_SIZE/(NUMBER_OF_MESSAGE)*NUMBER_OF_PARITY,sizeof(char));
		parity_2 = (char *)calloc(PAGE_SIZE/(NUMBER_OF_MESSAGE)*NUMBER_OF_PARITY,sizeof(char));
		parity_3 = (char *)calloc(PAGE_SIZE/(NUMBER_OF_MESSAGE)*NUMBER_OF_PARITY,sizeof(char));
		message = (char *)calloc(PAGE_SIZE/(NUMBER_OF_MESSAGE)*NUMBER_OF_MESSAGE,sizeof(char));

		parityBuff = (char *)malloc(sizeof(char)*NUMBER_OF_PARITY);
		parity_2Buff = (char *)malloc(sizeof(char)*NUMBER_OF_PARITY);
		parity_3Buff = (char *)malloc(sizeof(char)*NUMBER_OF_PARITY);
		messageBuff = (char *)malloc(sizeof(char)*NUMBER_OF_MESSAGE);
		FILE *f = fopen("./write.rec", "r");
		for(int i=0;i<PAGE_SIZE/NUMBER_OF_MESSAGE;++i)
		{
			fgets(parityBuff, NUMBER_OF_PARITY+1, f);
			fgets(parity_2Buff, NUMBER_OF_PARITY+1, f);
			fgets(parity_3Buff, NUMBER_OF_PARITY+1, f);
			fgets(messageBuff, NUMBER_OF_MESSAGE+1, f);
			strcat(parity,parityBuff);
			strcat(parity_2,parity_2Buff);
			strcat(parity_3,parity_3Buff);
			strcat(message,messageBuff);
		}
		*/
		char *a = "a";
		event.set_payload(a);
		event.set_payload_parity(a);
		ssd.write(event);
		write_count += 1;
		//fclose(f);
		//only spare write.
		if(ssd.get_tail((resolve_logical_address(ssd.dltodp[event.get_logical_address()])))->get_peCycle() > -1)
		{
	
			if(eccBufferCount == 4)
			{
				event.set_payload(eccBuffer);
				event.set_event_type(ECCWRITE);
				printf("%d\n", eccBufferCount);
				write(event, OPECC);

				eccBufferCount = 0;
				ecc_write_count += 1;
				//memcpy((char*)eccBuffer + (eccBufferCount * PAGE_SIZE / 4), parity_2, PAGE_SIZE / 4);
				eccBufferPpnIndex[eccBufferCount] = ssd.dltodp[event.get_logical_address()];
				eccBufferLpnIndex[eccBufferCount] = event.get_logical_address();
			}
			else
			{

				//memcpy((char*)eccBuffer + (eccBufferCount * PAGE_SIZE / 4), parity_2, PAGE_SIZE / 4);
				eccBufferPpnIndex[eccBufferCount] = ssd.dltodp[event.get_logical_address()];
				eccBufferLpnIndex[eccBufferCount] = event.get_logical_address();
			}
			eccBufferCount++;
		}
		/*free(parity);
		free(parity_2);
		free(message);
		free(parityBuff);
		free(parity_2Buff);
		free(messageBuff);*/
		print_statistic();
		return SUCCESS;
	}
	else if (type == OPECC)
	{
		printf("EXTEND WRITE\n");
		for(int i=0; i<4; ++i)
		{
			if(ssd.dltoep[eccBufferLpnIndex[i]] != -1)
			{
				int eccInvalidCount = 0;
				bm->add_invalid_num( (ssd.dltoep[eccBufferLpnIndex[i]]>>2)/BLOCK_SIZE, 1);
				int ppn = ssd.dltoep[eccBufferLpnIndex[i]]>>2;
				int ecc_num = ssd.dltoep[eccBufferLpnIndex[i]] - (ppn<<2);
				ssd.get_tail(resolve_logical_address(ssd.dltoep[eccBufferLpnIndex[i]]>>2))->set_ecc_state(ecc_num, INVALID);
				for(int j = 0; j < 4;++j)
				{
					if(ssd.get_tail(resolve_logical_address(ssd.dltoep[eccBufferLpnIndex[i]]>>2))->get_ecc_state(j) == INVALID)
					{
						eccInvalidCount+=1;
					}
					if(eccInvalidCount == 4)
					{
						ssd.get_tail(resolve_logical_address(ssd.dltoep[eccBufferLpnIndex[i]]>>2))->set_state(INVALID);
					}
				}
			}
		}
		int update_ecc_ppn = insert_event(event, OPECC);

		for(int i=0;i<4;++i)
		{
			ssd.dltoep[eccBufferLpnIndex[i]] = update_ecc_ppn;
			ssd.dltoep[eccBufferLpnIndex[i]] = ssd.dltoep[eccBufferLpnIndex[i]]<<2;
			ssd.dltoep[eccBufferLpnIndex[i]] += i;
			ssd.get_tail(resolve_logical_address(ssd.dltoep[eccBufferLpnIndex[i]]>>2))->set_ecc_state(i, VALID);
			printf("no.%d lpn = %d, ppn = %d, eccnum = %d\n", i,eccBufferLpnIndex[i], ssd.dltoep[eccBufferLpnIndex[i]]>>2, ssd.dltoep[eccBufferLpnIndex[i]]-((ssd.dltoep[eccBufferLpnIndex[i]]>>2)<<2));
		}
		//***************//
		Address address = resolve_logical_address(update_ecc_ppn);
		event.set_address(address);
		ssd.write(event);
		for(int i=0;i<4;++i)
		{
			event.set_eccLpn(eccBufferLpnIndex[i], i);
		}
		ram.ecc_write(event, *this);
		return SUCCESS;
	}
}
int FtlImpl_Fast::insert_event(Event& event, optype type)
{
	if (type == OPDATA)
	{
		//******first require
		if (bm->get_free_block_num() == SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE)
		{
			printf("First SSD Write\n");
			currentBlock = 0;
			bm->set_current_page(0, ++currentPage);
			bm->set_block_state(0, ACTIVE);
			bm->set_block_type(0, DATA);
			bm->set_free_block_num(-1);
			bm->add_free_num(0, -1);
			bm->add_valid_num(0, 1);
			return 0;
		}
		if (bm->get_free_block_num()<7 && currentPage == 63)
		{
			printf("NEED TO DATA GC\n");
			currentBlock = garbage_collect(type);
		}
		if (currentPage == 63)
		{
			for (int i = 0; i < SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE; ++i)
			{
				if(bm->get_block_state(i) == ACTIVE || bm->get_block_state(i) == INACTIVE)
					continue;
				currentBlock = i;
				currentPage = -1;
				break;
			}
			bm->set_current_page(currentBlock, ++currentPage);
			bm->set_block_type(currentBlock, DATA);
			bm->set_block_state(currentBlock, ACTIVE);
			bm->set_free_block_num(-1);
			bm->add_free_num(currentBlock, -1);
			bm->add_valid_num(currentBlock, 1);
		}
		else
		{
			bm->set_current_page(currentBlock, ++currentPage);
			bm->add_free_num(currentBlock, -1);
			bm->add_valid_num(currentBlock, 1);
		}
		return (currentBlock * BLOCK_SIZE) + currentPage;
	}
	else if (type == OPECC)
	{
		if (currentEccBlock == -1)
		{
			for (int i = 0; i < SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE; ++i)
			{
				if(bm->get_block_state(i) == ACTIVE || bm->get_block_state(i) == INACTIVE)
					continue;
				currentEccBlock = i;
				currentEccPage = -1;
				break;
			}
			bm->set_current_page(currentEccBlock, ++currentEccPage);
			bm->set_block_type(currentEccBlock, ECC);
			bm->set_block_state(currentEccBlock, ACTIVE);
			bm->set_free_block_num(-1);
			bm->add_free_num(currentEccBlock, -1);
			bm->add_valid_num(currentEccBlock, 1);
			eccBufferPieceIndex = 0;
			event.set_eccLpn(event.get_logical_address(), 0);
			return currentEccBlock * BLOCK_SIZE;
		}
		if (bm->get_free_block_num()<7  && currentEccPage == 63)
		{
			printf("NEED TO ECC GC\n");
			currentEccBlock = garbage_collect(type);
		}
		if (currentEccPage == BLOCK_SIZE - 1)
		{
			for (int i = 0; i < SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE; ++i)
			{
				if(bm->get_block_state(i) == ACTIVE || bm->get_block_state(i) == INACTIVE)
					continue;
				currentEccBlock = i;
				currentEccPage = -1;
				break;
			}
			bm->set_current_page(currentEccBlock, ++currentEccPage);
			bm->set_block_type(currentEccBlock, ECC);
			bm->set_block_state(currentEccBlock, ACTIVE);
			bm->set_free_block_num(-1);
			bm->add_free_num(currentEccBlock, -1);
			bm->add_valid_num(currentEccBlock, 1);
		}
		else
		{
			bm->set_current_page(currentEccBlock, ++currentEccPage);
			bm->add_free_num(currentEccBlock, -1);
			bm->add_valid_num(currentEccBlock, 1);
		}
		return (currentEccBlock * BLOCK_SIZE) + currentEccPage;
	}
	else if (type == OPERASE)
	{
	}
}
Address FtlImpl_Fast::resolve_logical_address(int ppn)
{
	Address phyAddress;
	phyAddress.set_linear_address(ppn);
	phyAddress.valid = PAGE;

	//fprintf(stderr, "numCells: %i package: %i die: %i plane: %i block: %i page: %i\n", numCells, phyAddress.package, phyAddress.die, phyAddress.plane, phyAddress.block, phyAddress.page);

	return phyAddress;
}

enum status FtlImpl_Fast::erase(Event& event)
{
	return SUCCESS;
}

enum status FtlImpl_Fast::merge(Event& event)
{
	return SUCCESS;
}

int FtlImpl_Fast::garbage_collect(optype type)
{
	int victim_block = 0;
	int free_block = -1;
	int invalid_num = 0;
	int page_offset = 0;
	int ecc_page_offset = 0;
	int old_ppn = 0;
	print_statistic();
	if(gc_free_block == SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE-1)
	{
		gc_free_block = 0;
	}

	if (type == OPDATA)
	{
		data_gc_count += 1;
		//printf("*****DATA GC START******\n\n");
			for (int i = gc_free_block; i < SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE; ++i)
			{
				if(bm->get_block_state(i) == ACTIVE || bm->get_block_state(i) == INACTIVE)
					continue;
				free_block = i;
				break;
			}
		if(free_block == -1)
		{
			gc_free_block = 0;
			for (int i = gc_free_block; i < SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE; ++i)
			{
				if(bm->get_block_state(i) == ACTIVE || bm->get_block_state(i) == INACTIVE)
					continue;
				free_block = i;
				break;
			}
		}
		if(free_block == -1)
		{
			free_block = rand()%256;
		}
		
		gc_free_block = free_block;
		for (int i = 0; i < SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE; ++i)//choose victim block number
		{
			if (bm->get_block_state(i) != ACTIVE || bm->get_block_type(i) == ECC)
			{
				continue;
			}
			else
			{
				if (bm->get_invalid_num(i) > invalid_num)
				{
					invalid_num = bm->get_invalid_num(i);
					victim_block = i;
				}
			}
		}
		if(free_block == -1)
		{
			printf("free block -1\n");
			exit(1);
		}
		if(victim_block == -1)
		{
			printf("victim_block -1\n");
			exit(1);
		}
		printf("VICTIM DATA BLOCK = %d\n", victim_block);
		old_ppn = victim_block * BLOCK_SIZE;

		bm->set_block_state(free_block, ACTIVE);
		bm->set_block_type(free_block, DATA);
		bm->set_free_block_num(-1);
		for (int i = 0; i < BLOCK_SIZE; ++i)
		{
			if (ssd.get_tail(resolve_logical_address((victim_block * BLOCK_SIZE) + i))->get_state() == VALID)
			{

				Event* event;
				Event* gc_read_event;
				int lpn = 0;
				for (int j = 0; j < numCells; ++j)
				{
					if (ssd.dltodp[j] == (victim_block * BLOCK_SIZE) + i)//old ppn
					{
						lpn = j;
						break;
					}
				}
				int oldppn = ssd.dltodp[lpn];
				ssd.dltodp[lpn] = (free_block * BLOCK_SIZE) + page_offset;
				event = new Event(WRITE, lpn, 1, 0);
				gc_read_event = new Event(GCREAD, oldppn, 1, 0);
				Address address = resolve_logical_address((free_block * BLOCK_SIZE) + page_offset);//new ppn
				Address gc_address = resolve_logical_address((oldppn));
				bm->set_current_page(free_block, page_offset);
				bm->add_free_num(free_block, -1);
				bm->add_valid_num(free_block, 1);
				bm->add_used_num(free_block, 1);
				event->set_address(address);
				gc_read_event->set_address(gc_address);
				ssd.read(*gc_read_event);
				event->set_payload(global_buffer);
				event->set_payload_parity(ssd.get_tail(resolve_logical_address(oldppn))->get_spare()->get_spare_ecc());
				controller->issue(*event);
				data_gc_read_count += 1;
				data_gc_write_count += 1;
				page_offset++;
				//printf("old ppn = %d, new ppn = %d\n", victim_block * BLOCK_SIZE + i, (free_block * BLOCK_SIZE) + page_offset);
				free(event);
				free(gc_read_event);
			}
		}
		bm->clean_block(victim_block);
		Event* erase_event;
		erase_event = new Event(ERASE, 0, 1, 0);
		Address erase_address = new Address(old_ppn, PAGE);
		erase_event->set_address(erase_address);
		enum status erase_status = ssd.erase(*erase_event);
		bm->set_free_block_num(1);
		bm->set_block_state(victim_block, FREE);
		bm->set_block_type(victim_block, NOTHING);
		currentPage = page_offset-1;
		free(erase_event);
		//printf("*************\n\n");
		return free_block;
	}
	if (type == OPECC)
	{
		ecc_gc_count += 1;
			int count = 0;
	//printf("*****ECC GC START******\n\n");
			for (int i = 0; i < SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE; ++i)
			{
				if(bm->get_block_state(i) == ACTIVE || bm->get_block_state(i) == INACTIVE)
					continue;
				free_block = i;
				break;
			}
			printf("FREE BLOCK FOR ECC = %d\n", free_block);
		for (int i = 0; i < SSD_SIZE * PACKAGE_SIZE * DIE_SIZE * PLANE_SIZE; ++i)//choose victim block number
		{
			if (bm->get_block_state(i) != ACTIVE || bm->get_block_type(i) == DATA)
			{
				continue;
			}
			else
			{
				if (bm->get_invalid_num(i) > invalid_num)
				{
					invalid_num = bm->get_invalid_num(i);
					victim_block = i;
				}
			}
		}
		if(free_block == -1)
		{
			free_block = rand()%256;
		}
		printf("VICTIM ECC BLOCK = %d\n", victim_block);
		old_ppn = victim_block * BLOCK_SIZE;

		bm->set_block_state(free_block, ACTIVE);
		bm->set_block_type(free_block, ECC);
		bm->set_free_block_num(-1);




		bool ecc_state[4];
		int lpn_buff[4];
		for (int i = 0; i < BLOCK_SIZE; ++i)
		{
			for (int eccCount = 0; eccCount < 4; ++eccCount)
			{
				if (ssd.get_tail(resolve_logical_address((victim_block * BLOCK_SIZE) + i))->get_ecc_state(eccCount) == VALID)
				{
					printf("%d ECC PAGE %d PIECE = VALID\n",i, eccCount);
					ecc_state[eccCount] = true;
				}
				else
				{
					printf("%d ECC PAGE %d PIECE = INVALID\n",i, eccCount);
					ecc_state[eccCount] = false;
				}
			}
			Event* ecc_event;
			ecc_event = new Event(GCREAD, 0, 1, 0);
			Address address = resolve_logical_address((victim_block * BLOCK_SIZE) + i);
			ecc_event->set_address(address);
			ssd.read(*ecc_event);
			for (int j = 0; j < 4; j++)
			{
				if (ecc_state[j] == true)
				{
					for (int k = 0; k < numCells; k++)
					{
						if (ssd.dltoep[k]>>2 == (victim_block * BLOCK_SIZE) + i && (ssd.dltoep[k] - ((ssd.dltoep[k]>>2)<<2)) == j)
						{
							lpn_buff[count] = k;
							//memcpy((char*)eccBuf + (count * PAGE_SIZE / 2), (char*)global_buffer + (j * PAGE_SIZE / 2), PAGE_SIZE / 2);
							ecc_gc_read_count += 1;
							count++;
							break;
						}
					}
				}
				if (count == 4 || i == 63)
				{
					if(i != 63)
					{
						printf("ECC PAGE FLUSH\n");
						for (int l = 0; l < count; ++l)
						{
							ssd.get_tail(resolve_logical_address((free_block * BLOCK_SIZE) + ecc_page_offset))->set_ecc_state(l,VALID);
							ssd.dltoep[lpn_buff[l]] = (free_block * BLOCK_SIZE) + ecc_page_offset;
							ssd.dltoep[lpn_buff[l]] = ssd.dltoep[lpn_buff[l]] << 2;
							ssd.dltoep[lpn_buff[l]] += l;
							lpn_buff[l] = -1;
						}
						Event* event;
						event = new Event(ECCWRITE, 0, 1, 0);
						Address address = resolve_logical_address((free_block * BLOCK_SIZE) + ecc_page_offset);//new ppn
						bm->set_current_page(free_block, ecc_page_offset);
						bm->add_free_num(free_block, -1);
						bm->add_valid_num(free_block, 1);
						bm->add_used_num(free_block, 1);
						event->set_address(address);
						event->set_payload(eccBuf);
				  		ssd.write(*event);
						ecc_gc_write_count += 1;
						printf("old ppn = %d, new ppn = %d\n", victim_block * BLOCK_SIZE + i, (free_block * BLOCK_SIZE) + ecc_page_offset);
						ecc_page_offset++;
						free(event);
						count = 0;
					}
					else
					{
						printf("ECC PAGE FLUSH\n");
						for (int l = 0; l < count; ++l)
						{
							if(ecc_state[l] == true)
							{
								ssd.get_tail(resolve_logical_address((free_block * BLOCK_SIZE) + ecc_page_offset))->set_ecc_state(l,VALID);
								ssd.dltoep[lpn_buff[l]] = (free_block * BLOCK_SIZE) + ecc_page_offset;
								ssd.dltoep[lpn_buff[l]] = ssd.dltoep[lpn_buff[l]] << 2;
								ssd.dltoep[lpn_buff[l]] += l;
							}
						}
						Event* event;
						event = new Event(ECCWRITE, 0, 1, 0);
						ecc_gc_write_count += 1;
						Address address = resolve_logical_address((free_block * BLOCK_SIZE) + ecc_page_offset);//new ppn
						bm->set_current_page(free_block, ecc_page_offset);
						bm->add_free_num(free_block, -1);
						bm->add_valid_num(free_block, 1);
						bm->add_used_num(free_block, 1);
						event->set_address(address);
						event->set_payload(eccBuf);
						ssd.write(*event);
						ecc_gc_write_count += 1;
						printf("old ppn = %d, new ppn = %d\n", victim_block * BLOCK_SIZE + i, (free_block * BLOCK_SIZE) + ecc_page_offset);
						ecc_page_offset++;
						free(event);
						break;
					}
				}
			}
			free(ecc_event);
		}
		bm->clean_block(victim_block);
		Event* erase_event;
		erase_event = new Event(ERASE, 0, 1, 0);
		Address erase_address = new Address(old_ppn, PAGE);
		erase_event->set_address(erase_address);
		enum status erase_status = ssd.erase(*erase_event);
		bm->set_free_block_num(1);
		bm->set_block_state(victim_block, FREE);
		bm->set_block_type(victim_block, NOTHING);
		currentEccPage = ecc_page_offset-1;
		free(erase_event);

		return free_block;
	}
}
void FtlImpl_Fast::print_statistic()
{
	printf("READ DATA : %d\nREAD ECC : %d\nWRITE DATA : %d\nWRITE ECC : %d\nGC DATA READ : %d\nGC ECC READ : %d\nGC DATA WRITE : %d\nGC ECC WRITE : %d\nDATA GC COUNT : %d\nECC GC COUNT : %d\n",
	read_count, ecc_read_count, write_count, ecc_write_count, data_gc_read_count, ecc_gc_read_count, data_gc_write_count, ecc_gc_write_count, data_gc_count, ecc_gc_count);
}
