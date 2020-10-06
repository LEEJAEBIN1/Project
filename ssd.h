/* Copyright 2009, 2010 Brendan Tauras */

/* ssd.h is part of FlashSim. */

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

/* ssd.h
 * Brendan Tauras 2010-07-16
 * Main SSD header file
 * 	Lists definitions of all classes, structures,
 * 		typedefs, and constants used in ssd namespace
 *		Controls options, such as debug asserts and test code insertions
 */

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <queue>
#include <map>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/global_fun.hpp>
#include <boost/multi_index/random_access_index.hpp>
 
#ifndef _SSD_H
#define _SSD_H

extern "C"
{
extern int Ldpc_encode(int extend_no);
extern int Ldpc_decode(int exten_no);
extern void ldpc_transmit(double error);
}


namespace ssd {
#define NUMBER_OF_PARITY 4
#define	NUMBER_OF_MESSAGE 16

/* define exit codes for errors */
#define MEM_ERR -1
#define FILE_ERR -2

/* Uncomment to disable asserts for production */
#define NDEBUG


/* some obvious typedefs for laziness */
typedef unsigned int uint;
typedef unsigned long ulong;

//extern int c = 0;
/* Simulator configuration from ssd_config.cpp */

/* Configuration file parsing for extern config variables defined below */
void load_entry(char *name, double value, uint line_number);
void load_config(void);
void print_config(FILE *stream);

/* Ram class:
 * 	delay to read from and write to the RAM for 1 page of data */
extern const double RAM_READ_DELAY;
extern const double RAM_WRITE_DELAY;

/* Bus class:
 * 	delay to communicate over bus
 * 	max number of connected devices allowed
 * 	flag value to detect free table entry (keep this negative)
 * 	number of time entries bus has to keep track of future schedule usage
 * 	number of simultaneous communication channels - defined by SSD_SIZE */
extern const double BUS_CTRL_DELAY;
extern const double BUS_DATA_DELAY;
extern const uint BUS_MAX_CONNECT;
extern const double BUS_CHANNEL_FREE_FLAG;
extern const uint BUS_TABLE_SIZE;
/* extern const uint BUS_CHANNELS = 4; same as # of Packages, defined by SSD_SIZE */

/* Ssd class:
 * 	number of Packages per Ssd (size) */
extern const uint SSD_SIZE;
extern const uint NUM_OF_CACHE_DATA_BLOCK;
extern const uint NUM_OF_CACHE_ECC_BLOCK;

/* Package class:
 * 	number of Dies per Package (size) */
extern const uint PACKAGE_SIZE;

/* Die class:
 * 	number of Planes per Die (size) */
extern const uint DIE_SIZE;

/* Plane class:
 * 	number of Blocks per Plane (size)
 * 	delay for reading from plane register
 * 	delay for writing to plane register
 * 	delay for merging is based on read, write, reg_read, reg_write 
 * 		and does not need to be explicitly defined */
extern const uint PLANE_SIZE;
extern const double PLANE_REG_READ_DELAY;
extern const double PLANE_REG_WRITE_DELAY;

/* Block class:
 * 	number of Pages per Block (size)
 * 	number of erases in lifetime of block
 * 	delay for erasing block */
extern const uint BLOCK_SIZE;
extern const uint BLOCK_ERASES;
extern const double BLOCK_ERASE_DELAY;

/* Page class:
 * 	delay for Page reads
 * 	delay for Page writes */
extern const double PAGE_READ_DELAY;
extern const double PAGE_WRITE_DELAY;
extern const uint PAGE_SIZE;
extern const bool PAGE_ENABLE_DATA;
extern const uint ERROR_4_20;
extern const uint ERROR_8_24;
extern const uint ERROR_12_28;
/*
 * Mapping directory
 */
extern const uint MAP_DIRECTORY_SIZE;

/*
 * FTL Implementation
 */
extern const uint FTL_IMPLEMENTATION;

/*
 * LOG page limit for BAST.
 */
extern const uint BAST_LOG_PAGE_LIMIT;

/*
 * LOG page limit for FAST.
 */
extern const uint FAST_LOG_PAGE_LIMIT;

/*
 * Number of blocks allowed to be in DFTL Cached Mapping Table.
 */
extern const uint CACHE_DFTL_LIMIT;

/*
 * Parallelism mode
 */
extern const uint PARALLELISM_MODE;

/* Virtual block size (as a multiple of the physical block size) */
extern const uint VIRTUAL_BLOCK_SIZE;

/* Virtual page size (as a multiple of the physical page size) */
extern const uint VIRTUAL_PAGE_SIZE;

extern const uint NUMBER_OF_ADDRESSABLE_BLOCKS;

/* RAISSDs: Number of physical SSDs */
extern const uint RAID_NUMBER_OF_PHYSICAL_SSDS;

/* number of blokcs in RAM */
extern const uint RAM_BLOCK_SIZE;
/*
 * Memory area to support pages with data.
 */
extern void *page_data;
extern void *ecc_page_data;
extern void *global_buffer;

/* Enumerations to clarify status integers in simulation
 * Do not use typedefs on enums for reader clarity */

/* Page states
 * 	empty   - page ready for writing (and contains no valid data)
 * 	valid   - page has been written to and contains valid data
 * 	invalid - page has been written to and does not contain valid data */
enum page_state{EMPTY, VALID, INVALID};

/* Block states
 * 	free     - all pages in block are empty
 * 	active   - some pages in block are valid, others are empty or invalid
 * 	inactive - all pages in block are invalid */
enum block_state{FREE = 0, ACTIVE = 1, INACTIVE = 2};

/* I/O request event types
 * 	read  - read data from address
 * 	write - write data to address (page state set to valid)
 * 	erase - erase block at address (all pages in block are erased - 
 * 	                                page states set to empty)
 * 	merge - move valid pages from block at address (page state set to invalid)
 * 	           to free pages in block at merge_address */
enum event_type{READ, WRITE, ERASE, MERGE, TRIM, ECCWRITE, ECCREAD_STEP1, ECCREAD_STEP2, GCREAD};

/* General return status
 * return status for simulator operations that only need to provide general
 * failure notifications */
enum status{FAILURE, SUCCESS, READFAIL};


/* Address valid status
 * used for the valid field in the address class
 * example: if valid == BLOCK, then
 * 	the package, die, plane, and block fields are valid
 * 	the page field is not valid */
enum address_valid{NONE, PACKAGE, DIE, PLANE, BLOCK, PAGE};

/*
 * Block type status
 * used for the garbage collector specify what pool
 * it should work with.
 * the block types are log, data and map (Directory map usually)
 */
enum block_type {LOG, DATA, LOG_SEQ, ECC, NOTHING};

/*
 * Enumeration of the different FTL implementations.
 */
enum ftl_implementation {IMPL_PAGE, IMPL_BAST, IMPL_FAST, IMPL_DFTL, IMPL_BIMODAL};

enum ram_type {CACHE_DATA, CACHE_ECC, RAM_ECC, RAM_WRITE, RAM_READ};
enum optype {OPDATA, OPECC, OPERASE, ECCEXIST};
#define BOOST_MULTI_INDEX_ENABLE_SAFE_MODE 1

/* List classes up front for classes that have references to their "parent"
 * (e.g. a Package's parent is a Ssd).
 *
 * The order of definition below follows the order of this list to support
 * cases of agregation where the agregate class should be defined first.
 * Defining the agregate class first enables use of its non-default
 * constructors that accept args
 * (e.g. a Ssd contains a Controller, Ram, Bus, and Packages). */
class Address;
class Stats;
class Event;
class Channel;
class Bus;
class Ecc;
class Spare; //***
class Page;
class Block;
class EccBlock;
class EccPage;
class Plane;
class Die;
class Package;
class Garbage_Collector;
class Wear_Leveler;
class Block_manager;
class FtlParent;
class FtlImpl_Page;
class FtlImpl_Bast;
class FtlImpl_Fast;
class FtlImpl_DftlParent;
class FtlImpl_Dftl;
class FtlImpl_BDftl;
class Ftl_block_manager;
class Ram;
class Ram_block;
class Controller;
class Ssd;
class Ecc_block_page;



/* Class to manage physical addresses for the SSD.  It was designed to have
 * public members like a struct for quick access but also have checking,
 * printing, and assignment functionality.  An instance is created for each
 * physical address in the Event class. */
class Address
{
public:
	uint package;
	uint die;
	uint plane;
	uint block;
	uint page;
	ulong real_address;
	enum address_valid valid;
	Address(void);
	Address(const Address &address);
	Address(const Address *address);
	Address(uint package, uint die, uint plane, uint block, uint page, enum address_valid valid);
	Address(uint address, enum address_valid valid);
	~Address();
	enum address_valid check_valid(uint ssd_size = SSD_SIZE, uint package_size = PACKAGE_SIZE, uint die_size = DIE_SIZE, uint plane_size = PLANE_SIZE, uint block_size = BLOCK_SIZE);
	enum address_valid compare(const Address &address) const;
	void print(FILE *stream = stdout);

	void operator+(int);
	void operator+(uint);
	Address &operator+=(const uint rhs);
	Address &operator=(const Address &rhs);

	void set_linear_address(ulong address, enum address_valid valid);
	void set_linear_address(ulong address);
	ulong get_linear_address() const;
};

class Stats
{
public:
	// Flash Translation Layer
	long numFTLRead;
	long numFTLWrite;
	long numFTLErase;
	long numFTLTrim;

	// Garbage Collection
	long numGCRead;
	long numGCWrite;
	long numGCErase;

	// Wear-leveling
	long numWLRead;
	long numWLWrite;
	long numWLErase;

	// Log based FTL's
	long numLogMergeSwitch;
	long numLogMergePartial;
	long numLogMergeFull;

	// Page based FTL's
	long numPageBlockToPageConversion;

	// Cache based FTL's
	long numCacheHits;
	long numCacheFaults;

	// Memory consumptions (Bytes)
	long numMemoryTranslation;
	long numMemoryCache;

	long numMemoryRead;
	long numMemoryWrite;

	// Advance statictics
	double translation_overhead() const;
	double variance_of_io() const;
	double cache_hit_ratio() const;

	// Constructors, maintainance, output, etc.
	Stats(void);

	void print_statistics();
	void reset_statistics();
	void write_statistics(FILE *stream);
	void write_header(FILE *stream);
private:
	void reset();
};

/* Class to emulate a log block with page-level mapping. */
class LogPageBlock
{
public:
	LogPageBlock(void);
	~LogPageBlock(void);

	int *pages;
	long *aPages;
	Address address;
	int numPages;

	LogPageBlock *next;

	bool operator() (const ssd::LogPageBlock& lhs, const ssd::LogPageBlock& rhs) const;
	bool operator() (const ssd::LogPageBlock*& lhs, const ssd::LogPageBlock*& rhs) const;
};


/* Class to manage I/O requests as events for the SSD.  It was designed to keep
 * track of an I/O request by storing its type, addressing, and timing.  The
 * SSD class creates an instance for each I/O request it receives. */
class Event
{
public:
	Event(enum event_type type, ulong logical_address, uint size, double start_time);
	Event();
	~Event(void);
	void consolidate_metaevent(Event &list);
	ulong get_logical_address(void) const;
	const Address &get_address(void) const;
	const Address &get_merge_address(void) const;
	const Address &get_log_address(void) const;
	const Address &get_replace_address(void) const;
	uint get_size(void) const;
	enum event_type get_event_type(void) const;
	double get_start_time(void) const;
	double get_time_taken(void) const;
	double get_bus_wait_time(void) const;
	bool get_noop(void) const;
	Event *get_next(void) const;
	void set_address(const Address &address);
	void set_merge_address(const Address &address);
	void set_log_address(const Address &address);
	void set_replace_address(const Address &address);
	void set_next(Event &next);
	void set_payload(void *payload);
	void set_payload_parity(void *payload_parity);
	void set_event_type(const enum event_type &type);
	void set_noop(bool value);
	void *get_payload(void) const;
	void *get_payload_parity(void) const;
	void set_eccLpn(int lpn, int eccnum);
	int get_eccLpn(int eccnum);
	double incr_bus_wait_time(double time);
	double incr_time_taken(double time_incr);
	void set_eccnum(int eccnum);
	int get_eccnum();
	void print(FILE *stream = stdout);
private:
	double start_time;
	double time_taken;
	double bus_wait_time;
	enum event_type type;

	int eccLpn[4];
	ulong logical_address;
	Address address;
	Address merge_address;
	Address log_address;
	Address replace_address;
	uint size;
	void *payload;
	void *payload_parity;
	Event *next;
	bool noop;
	int ecc_num;
};

/* Single bus channel
 * Simulate multiple devices on 1 bus channel with variable bus transmission
 * durations for data and control delays with the Channel class.  Provide the 
 * delay times to send a control signal or 1 page of data across the bus
 * channel, the bus table size for the maximum number channel transmissions that
 * can be queued, and the maximum number of devices that can connect to the bus.
 * To elaborate, the table size is the size of the channel scheduling table that
 * holds start and finish times of events that have not yet completed in order
 * to determine where the next event can be scheduled for bus utilization. */
class Channel
{
public:
	Channel(double ctrl_delay = BUS_CTRL_DELAY, double data_delay = BUS_DATA_DELAY, uint table_size = BUS_TABLE_SIZE, uint max_connections = BUS_MAX_CONNECT);
	~Channel(void);
	enum status lock(double start_time, double duration, Event &event);
	enum status connect(void);
	enum status disconnect(void);
	double ready_time(void);
private:
	void unlock(double current_time);

	struct lock_times {
		double lock_time;
		double unlock_time;
	};

	static bool timings_sorter(lock_times const& lhs, lock_times const& rhs);
	std::vector<lock_times> timings;

	uint table_entries;
	uint selected_entry;
	uint num_connected;
	uint max_connections;
	double ctrl_delay;
	double data_delay;

	// Stores the highest unlock_time in the vector timings list.
	double ready_at;
};

/* Multi-channel bus comprised of Channel class objects
 * Simulates control and data delays by allowing variable channel lock
 * durations.  The sender (controller class) should specify the delay (control,
 * data, or both) for events (i.e. read = ctrl, ctrl+data; write = ctrl+data;
 * erase or merge = ctrl).  The hardware enable signals are implicitly
 * simulated by the sender locking the appropriate bus channel through the lock
 * method, then sending to multiple devices by calling the appropriate method
 * in the Package class. */
class Bus
{
public:
	Bus(uint num_channels = SSD_SIZE, double ctrl_delay = BUS_CTRL_DELAY, double data_delay = BUS_DATA_DELAY, uint table_size = BUS_TABLE_SIZE, uint max_connections = BUS_MAX_CONNECT);
	~Bus(void);
	enum status lock(uint channel, double start_time, double duration, Event &event);
	enum status connect(uint channel);
	enum status disconnect(uint channel);
	Channel &get_channel(uint channel);
	double ready_time(uint channel);
private:
	uint num_channels;
	Channel * const channels;
};


/* The page is the lowest level data storage unit that is the size unit of
 * requests (events).  Pages maintain their state as events modify them. */
class Page
{
public:
	Page(const Block &parent, double read_delay = PAGE_READ_DELAY, double write_delay = PAGE_WRITE_DELAY);
	~Page(void);
	enum status _read(Event &event);
	enum status _write(Event &event);
	const Block &get_parent(void) const;
	enum page_state get_state(void) const;
	void set_state(enum page_state state);
	enum page_state get_ecc_state(int eccnum);
	void set_ecc_state(int eccnum, enum page_state state);
	void add_peCycle(int num);
	int get_peCycle();
	Spare* get_spare();
private:
	enum page_state* eccState;
	enum page_state state;
	const Block &parent;
	double read_delay;
	double write_delay;
	int peCycle = 0;
	Spare * spare;
};

class Spare
{
private:
	uint lpn;
	uint ppn;
	char* ecc;
public:
	Spare(void);
	~Spare(void);
	void set_spare_lpn(uint data);
	void set_spare_ppn(uint data);
	void set_spare_ecc(char *data);
	uint get_spare_lpn();
	uint get_spare_ppn();
	char* get_spare_ecc();
	void print_ecc();
};

/* The block is the data storage hardware unit where erases are implemented.
 * Blocks maintain wear statistics for the FTL. */
class Block 
{
public:
	long physical_address;
	uint pages_invalid;
	Block(const Plane &parent, uint size = BLOCK_SIZE, ulong erases_remaining = BLOCK_ERASES, double erase_delay = BLOCK_ERASE_DELAY, long physical_address = 0);
	~Block(void);
	enum status read(Event &event);
	enum status write(Event &event);
	enum status replace(Event &event);
	enum status _erase(Event &event);
	const Plane &get_parent(void) const;
	uint get_pages_valid(void) const;
	uint get_pages_invalid(void) const;
	enum block_state get_state(void) const;
	enum page_state get_state(uint page) const;
	enum page_state get_state(const Address &address) const;
	double get_last_erase_time(void) const;
	double get_modification_time(void) const;
	ulong get_erases_remaining(void) const;
	uint get_size(void) const;
	enum status get_next_page(Address &address) const;
	void invalidate_page(uint page);
	long get_physical_address(void) const;
	Block *get_pointer(void);
	block_type get_block_type(void) const;
	void set_block_type(block_type value);
	Page *get_tail(Address address);

private:
	uint size;
	Page * const data;
	const Plane &parent;
	uint pages_valid;
	enum block_state state;
	ulong erases_remaining;
	double last_erase_time;
	double erase_delay;
	double modification_time;

	block_type btype;
};

class EccBlock
{
public:
	EccBlock(int size);
	~EccBlock(void);
	enum status read(int EccPageAddress);
	enum status write(int EccPageAddress);
private:	
	EccPage * data;
	int size;
};

class EccPage
{
public:
	EccPage(double read_delay = PAGE_READ_DELAY, double write_delay = PAGE_WRITE_DELAY);
	~EccPage();
	enum status write();
	enum status read();
private:
	double read_delay;
	double write_delay;
};
/* The plane is the data storage hardware unit that contains blocks.
 * Plane-level merges are implemented in the plane.  Planes maintain wear
 * statistics for the FTL. */
class Plane 
{
public:
	Plane(const Die &parent, uint plane_size = PLANE_SIZE, double reg_read_delay = PLANE_REG_READ_DELAY, double reg_write_delay = PLANE_REG_WRITE_DELAY, long physical_address = 0);
	~Plane(void);
	enum status read(Event &event);
	enum status write(Event &event);
	enum status erase(Event &event);
	enum status replace(Event &event);
	enum status _merge(Event &event);
	const Die &get_parent(void) const;
	double get_last_erase_time(const Address &address) const;
	ulong get_erases_remaining(const Address &address) const;
	void get_least_worn(Address &address) const;
	uint get_size(void) const;
	enum page_state get_state(const Address &address) const;
	enum block_state get_block_state(const Address &address) const;
	void get_free_page(Address &address) const;
	ssd::uint get_num_free(const Address &address) const;
	ssd::uint get_num_valid(const Address &address) const;
	ssd::uint get_num_invalid(const Address &address) const;
	Block *get_block_pointer(const Address & address);
	Page *get_tail(Address address);
private:
	void update_wear_stats(void);
	enum status get_next_page(void);
	uint size;
	Block * const data;
	const Die &parent;
	uint least_worn;
	ulong erases_remaining;
	double last_erase_time;
	double reg_read_delay;
	double reg_write_delay;
	Address next_page;
	uint free_blocks;
};

/* The die is the data storage hardware unit that contains planes and is a flash
 * chip.  Dies maintain wear statistics for the FTL. */
class Die 
{
public:
	Die(const Package &parent, Channel &channel, uint die_size = DIE_SIZE, long physical_address = 0);
	~Die(void);
	enum status read(Event &event);
	enum status write(Event &event);
	enum status erase(Event &event);
	enum status replace(Event &event);
	enum status merge(Event &event);
	enum status _merge(Event &event);
	const Package &get_parent(void) const;
	double get_last_erase_time(const Address &address) const;
	ulong get_erases_remaining(const Address &address) const;
	void get_least_worn(Address &address) const;
	enum page_state get_state(const Address &address) const;
	enum block_state get_block_state(const Address &address) const;
	void get_free_page(Address &address) const;
	ssd::uint get_num_free(const Address &address) const;
	ssd::uint get_num_valid(const Address &address) const;
	ssd::uint get_num_invalid(const Address &address) const;
	Block *get_block_pointer(const Address & address);
	Page *get_tail(Address address);
private:
	void update_wear_stats(const Address &address);
	uint size;
	Plane * const data;
	const Package &parent;
	Channel &channel;
	uint least_worn;
	ulong erases_remaining;
	double last_erase_time;
};

/* The package is the highest level data storage hardware unit.  While the
 * package is a virtual component, events are passed through the package for
 * organizational reasons, including helping to simplify maintaining wear
 * statistics for the FTL. */
class Package 
{
public:
	Package (const Ssd &parent, Channel &channel, uint package_size = PACKAGE_SIZE, long physical_address = 0);
	~Package ();
	enum status read(Event &event);
	enum status write(Event &event);
	enum status erase(Event &event);
	enum status replace(Event &event);
	enum status merge(Event &event);
	const Ssd &get_parent(void) const;
	double get_last_erase_time (const Address &address) const;
	ulong get_erases_remaining (const Address &address) const;
	void get_least_worn (Address &address) const;
	enum page_state get_state(const Address &address) const;
	enum block_state get_block_state(const Address &address) const;
	void get_free_page(Address &address) const;
	ssd::uint get_num_free(const Address &address) const;
	ssd::uint get_num_valid(const Address &address) const;
	ssd::uint get_num_invalid(const Address &address) const;
	Block *get_block_pointer(const Address & address);
	Page *get_tail(Address address);
private:
	void update_wear_stats (const Address &address);
	uint size;
	Die * const data;
	const Ssd &parent;
	uint least_worn;
	ulong erases_remaining;
	double last_erase_time;
};

/* place-holder definitions for GC, WL, FTL, RAM, Controller 
 * please make sure to keep this order when you replace with your definitions */
class Garbage_collector 
{
public:
	Garbage_collector(FtlParent &ftl);
	~Garbage_collector(void);
private:
	void clean(Address &address);
};

class Wear_leveler 
{
public:
	Wear_leveler(FtlParent &FTL);
	~Wear_leveler(void);
	enum status insert(const Address &address);
};
class Ftl_block_manager
{
public:
	Ftl_block_manager(FtlImpl_Fast* ftl);
	~Ftl_block_manager();
	void set_block_type(int block_num, enum block_type btype);
	enum block_type get_block_type(int block_num);
	void set_block_state(int block_num, enum block_state bstate);
	enum block_state get_block_state(int block_num);
	void add_valid_num(int block_num, int add_num);
	int get_valid_num(int block_num);
	void add_invalid_num(int block_num, int add_num);
	int get_invalid_num(int block_num);
	void add_free_num(int block_num, int add_num);
	int get_free_num(int block_num);
	void set_current_page(int block_num, int page_num);
	int get_current_page(int block_num);
	void add_used_num(int block_num, int add_num);
	int get_used_num(int block_num);
	void clean_block(int block_num);
	void init_block(int block_num);
	int get_free_block_num();
	void set_free_block_num(int num);
	int get_valid_page_num();
	void set_used_page_num(int num);
	int get_used_page_num();

private:
	int* num_valid;
	int* num_invalid;
	int* num_free;
	int* current_page;
	int* used_num;
	FtlImpl_Fast* ftl;
	enum block_state* block_state;
	enum block_type* block_type;
	int free_block_num;
	int used_page_num=0;
};

class Block_manager
{
public:
	Block_manager(FtlParent *ftl);
	~Block_manager(void);




	// Usual suspects
	Address get_free_block(Event &event);
	Address get_free_block(block_type btype, Event &event);
	void invalidate(Address address, block_type btype);
	void print_statistics();
	void insert_events(Event &event);
	void promote_block(block_type to_type);
	bool is_log_full();
	void erase_and_invalidate(Event &event, Address &address, block_type btype);
	int get_num_free_blocks();

	// Used to update GC on used pages in blocks.
	void update_block(Block * b);

	// Singleton
	static Block_manager *instance();
	static void instance_initialize(FtlParent *ftl);
	static Block_manager *inst;

	void cost_insert(Block *b);

	void print_cost_status();



private:
	void get_page_block(Address &address, Event &event);
	static bool block_comparitor_simple (Block const *x,Block const *y);

	FtlParent *ftl;

	ulong data_active;
	ulong log_active;
	ulong logseq_active;

	ulong max_log_blocks;
	ulong max_blocks;

	ulong max_map_pages;
	ulong map_space_capacity;

	// Cost/Benefit priority queue.
	typedef boost::multi_index_container<
			Block*,
			boost::multi_index::indexed_by<
				boost::multi_index::random_access<>,
				boost::multi_index::ordered_non_unique<BOOST_MULTI_INDEX_MEMBER(Block,uint,pages_invalid) >
		  >
		> active_set;

	typedef active_set::nth_index<0>::type ActiveBySeq;
	typedef active_set::nth_index<1>::type ActiveByCost;

	active_set active_cost;

	// Usual block lists
	std::vector<Block*> active_list;
	std::vector<Block*> free_list;
	std::vector<Block*> invalid_list;

	// Counter for returning the next free page.
	ulong directoryCurrentPage;
	// Address on the current cached page in SRAM.
	ulong directoryCachedPage;

	ulong simpleCurrentFree;

	// Counter for handling periodic sort of active_list
	uint num_insert_events;

	uint current_writing_block;

	bool inited;

	bool out_of_blocks;





	int numInvalid;
	int address;


};


class FtlParent
{
public:
	FtlParent(Controller *controller, Ssd &ssd);
	virtual ~FtlParent () {};
	virtual enum status read(Event& event, enum optype type)=0;
	virtual enum status write(Event& event, enum optype type)=0;
	friend class Block_manager;


	ulong get_erases_remaining(const Address &address) const;
	void get_least_worn(Address &address) const;
	enum page_state get_state(const Address &address) const;
	enum block_state get_block_state(const Address &address) const;
	Block *get_block_pointer(const Address & address);
	Address resolve_logical_address(unsigned int logicalAddress);
	void print_ftl_statistics();
protected:
	Controller *controller;
	Ssd &ssd;
};

class FtlImpl_Fast : public FtlParent
{
public:
	FtlImpl_Fast(Controller *controller, Ssd &ssd, Ram &ram);
	~FtlImpl_Fast();
	enum status read(Event& event, enum optype type);
	enum status write(Event& event, enum optype type);
	Address resolve_logical_address(int ppn);
	enum status erase(Event& event);
	enum status merge(Event& event);
	int garbage_collect(optype type);
	ssd::ulong get_erases_remaining(const Address& address) const;
	void get_least_worn(Address& address) const;
	enum page_state get_state(const Address& address) const;
	int insert_event(Event& event, optype type);
	void print_statistic();
private:
	int numCells;
	int gcRatio;

	int* blockInvalid;
	char *parity;
	char *parity_2;
	char *parity_3;
	char *message;
	char *parityBuff;
	char *parity_2Buff;
	char *parity_3Buff;
	char *messageBuff;
	int validNum;
	int currentPage = -1;
	int currentEccPage = -1;
	int currentBlock = -1;
	int currentEccBlock = -1;
	int ecc_gc_write_count = 0;
	int ecc_gc_read_count = 0;
	int data_gc_write_count = 0; 
	int data_gc_read_count = 0;
	int write_count = 0;
	int read_count = 0;
	int ecc_write_count = 0;
	int ecc_read_count = 0;
	int data_gc_count = 0;
	int ecc_gc_count = 0;

	void* eccBuf;
	void* eccBuffer;
	int eccBufferLpnIndex[4];
	int eccBufferPpnIndex[4];
	int eccBufferCount = 0;
	int eccBufCount = 0;
	int eccBufferPieceIndex = 0;
  int gc_free_block = 0;
	Controller *controller;
	Ssd &ssd;
	Ram &ram;
	Ftl_block_manager *bm;
};

class Ram 
{
public:
	Ram(double read_delay, double write_delay, Ssd &parent);
	~Ram(void);
	enum status read(Event &event, FtlParent &ftl);
	enum status write(Event &event, FtlParent &ftl);
	enum status ecc_write(Event &event, FtlParent &ftl);
	//enum ram_return block_replace(Event &event);
	int choose_block();
	int choose_cache_block(enum ram_type type, int cache_block_index, int cache_block_tag, Event &event, FtlParent &ftl);
	int cache_flush(FtlParent& ftl);
	void print();
private:
	short current_ecc_page;
	char currentBlock=0;
	char usedBlock=0;
	int blockPointer=0;

	int dirtyMiss=0;
	int Hit=0;
	int EccHit=0;
	int nondirtyMiss=0;
	int Miss=0;
	double read_delay;
	double write_delay;
	void* dram;
	char* dirty;
	char* valid;
	short* tag;
	Ssd &ssd;
	Ram_block *blockTable;

};
class Ram_block
{
public:
	Ram_block(void);
	Ram_block(int i);
	~Ram_block(void);
	
	Event *pageEvent;
	Ecc_block_page *ecc_block_page;
	void clear_block();
	bool get_block_valid();
	int get_start_address();
	char get_num_valid();
	char get_num_dirty();
	char get_ppn();
	bool get_valid(char index);
	bool get_dirty(char index);

	void set_block_type(enum ram_type state);
	enum ram_type get_block_type();
	void set_block_valid(bool valid);
	void set_start_address(int address);
	void set_ppn(char ppn);
	void set_valid(char index, bool valid);
	void set_dirty(char index, bool dirty);
	void add_num_valid(char num);
	void add_num_dirty(char num);
	void set_tag(int index, int tag);
	int get_tag(int index);
	void set_block_lpn_pointer(int pointer);
	int get_block_lpn_pointer();
private:
	int block_lpn_pointer;
	int used_block_number;
	enum ram_type state;
	bool dirty[64];
	bool valid[64];
	int tag[64];
	bool block_valid;	
	int start_address;
	char num_valid;
	char num_dirty;
	char PPN;
};
/* The controller accepts read/write requests through its event_arrive method
 * and consults the FTL regarding what to do by calling the FTL's read/write
 * methods.  The FTL returns an event list for the controller through its issue
 * method that the controller buffers in RAM and sends across the bus.  The
 * controller's issue method passes the events from the FTL to the SSD.
 *
 * The controller also provides an interface for the FTL to collect wear
 * information to perform wear-leveling.  */
class Controller 
{
public:
	Controller(Ssd &parent);
	~Controller(void);
	enum status event_arrive(Event &event);
	friend class FtlParent;
	friend class FtlImpl_Page;
	friend class FtlImpl_Bast;
	friend class FtlImpl_Fast;
	friend class FtlImpl_DftlParent;
	friend class FtlImpl_Dftl;
	friend class FtlImpl_BDftl;
	friend class Block_manager;

	Stats stats;
	void print_ftl_statistics();
	const FtlParent &get_ftl(void) const;
private:
	enum status issue(Event &event_list);	
	void ecc_issue(void *ecc, int ecc_block_ppn, int ecc_page_ppn, enum event_type type);
	void translate_address(Address &address);
	ssd::ulong get_erases_remaining(const Address &address) const;
	void get_least_worn(Address &address) const;
	double get_last_erase_time(const Address &address) const;
	enum page_state get_state(const Address &address) const;
	enum block_state get_block_state(const Address &address) const;
	void get_free_page(Address &address) const;
	ssd::uint get_num_free(const Address &address) const;
	ssd::uint get_num_valid(const Address &address) const;
	ssd::uint get_num_invalid(const Address &saddress) const;
	Block *get_block_pointer(const Address & address);
	Ssd &ssd;
	Ecc *ecc;
	FtlParent *ftl;
	char *parity;
	char *message;
	char *parityBuff;
	char *messageBuff;
};

/* The SSD is the single main object that will be created to simulate a real
 * SSD.  Creating a SSD causes all other objects in the SSD to be created.  The
 * event_arrive method is where events will arrive from DiskSim. */
class Ssd
{
public:
	Ssd (uint ssd_size = SSD_SIZE);
	~Ssd(void);
	enum status read(Event &event);
	enum status write(Event &event);
	enum status erase(Event &event);
	double event_arrive(enum event_type type, ulong logical_address, uint size, double start_time);
	double event_arrive(enum event_type type, ulong logical_address, uint size, double start_time, void *buffer);
	void *get_result_buffer();
	friend class Controller;
	void print_statistics();
	void reset_statistics();
	void write_statistics(FILE *stream);
	void write_header(FILE *stream);
	const Controller &get_controller(void) const;
	int* dltodp;
	int* dltoep;
	void print_ftl_statistics();
	double ready_at(void);
	Page *get_tail(Address address);
private:

	enum status merge(Event &event);
	enum status replace(Event &event);
	enum status merge_replacement_block(Event &event);
	ulong get_erases_remaining(const Address &address) const;
	void update_wear_stats(const Address &address);
	void get_least_worn(Address &address) const;
	double get_last_erase_time(const Address &address) const;
	Package &get_data(void);
	enum page_state get_state(const Address &address) const;
	enum block_state get_block_state(const Address &address) const;
	void get_free_page(Address &address) const;
	ssd::uint get_num_free(const Address &address) const;
	ssd::uint get_num_valid(const Address &address) const;
	ssd::uint get_num_invalid(const Address &address) const;
	Block *get_block_pointer(const Address & address);

	uint size;
	Controller controller;
	Ram ram;
	Bus bus;
	Package * const data;
	ulong erases_remaining;
	ulong least_worn;
	double last_erase_time;
};

class RaidSsd
{
public:
	RaidSsd (uint ssd_size = SSD_SIZE);
	~RaidSsd(void);
	double event_arrive(enum event_type type, ulong logical_address, uint size, double start_time);
	double event_arrive(enum event_type type, ulong logical_address, uint size, double start_time, void *buffer);
	void *get_result_buffer();
	friend class Controller;
	void print_statistics();
	void reset_statistics();
	void write_statistics(FILE *stream);
	void write_header(FILE *stream);
	const Controller &get_controller(void) const;

	void print_ftl_statistics();
private:
	uint size;

	Ssd *Ssds;

};
class Ecc_block_page
{
public:
	Ecc_block_page();
	~Ecc_block_page();
	void add_ecc_count();
	int get_ecc_count();
	void set_ecc_ppn(int index, int ppn);
	void set_ecc_count(int count);
	int get_ecc_ppn(int index);
private:
	int ecc_ppn[4];
	int ecc_count=0;
	
};
} /* end namespace ssd */

#endif
