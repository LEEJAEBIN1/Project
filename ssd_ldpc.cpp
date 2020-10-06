#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <limits.h>
#include <sys/types.h>
#include "ssd.h"
#include <time.h>
using namespace ssd;

extern "C"
{
double error_prob;
typedef struct mod2entry /* Structure representing a non-zero entry, or
			      the header for a row or column               */
{
  int row, col;		  /* Row and column indexes of this entry, starting
                             at 0, and with -1 for a row or column header  */

  struct mod2entry *left, *right,  /* Pointers to entries adjacent in row  */
                   *up, *down;     /*   and column, or to headers.  Free   */
                                   /*   entries are linked by 'left'.      */

  double pr, lr;	  /* Probability and likelihood ratios - not used  */
			  /*   by the mod2sparse module itself             */
} mod2entry;

#define Mod2sparse_block 10  /* Number of entries to block together for
                                memory allocation */

typedef struct mod2block /* Block of entries allocated all at once */
{
  struct mod2block *next;  /* Next block that has been allocated */

  mod2entry entry[Mod2sparse_block]; /* Entries in this block */

} mod2block;

typedef struct		/* Representation of a sparse matrix */
{
  int n_rows;		  /* Number of rows in the matrix */
  int n_cols;		  /* Number of columns in the matrix */

  mod2entry *rows;	  /* Pointer to array of row headers */
  mod2entry *cols;	  /* Pointer to array of column headers */

  mod2block *blocks;	  /* Blocks that have been allocated */
  mod2entry *next_free;	  /* Next free entry */

} mod2sparse;


/* MACROS TO GET AT ELEMENTS OF A SPARSE MATRIX.  The 'first', 'last', 'next',
   and 'prev' macros traverse the elements in a row or column.  Moving past
   the first/last element gets one to a header element, which can be identified
   using the 'at_end' macro.  Macros also exist for finding out the row
   and column of an entry, and for finding out the dimensions of a matrix. */

#define mod2sparse_first_in_row(m,i) ((m)->rows[i].right) /* Find the first   */
#define mod2sparse_first_in_col(m,j) ((m)->cols[j].down)  /* or last entry in */
#define mod2sparse_last_in_row(m,i) ((m)->rows[i].left)   /* a row or column  */
#define mod2sparse_last_in_col(m,j) ((m)->cols[j].up)

#define mod2sparse_next_in_row(e) ((e)->right)  /* Move from one entry to     */
#define mod2sparse_next_in_col(e) ((e)->down)   /* another in any of the four */
#define mod2sparse_prev_in_row(e) ((e)->left)   /* possible directions        */
#define mod2sparse_prev_in_col(e) ((e)->up)

#define mod2sparse_at_end(e) ((e)->row<0) /* See if we've reached the end     */

#define mod2sparse_row(e) ((e)->row)      /* Find out the row or column index */
#define mod2sparse_col(e) ((e)->col)      /* of an entry (indexes start at 0) */

#define mod2sparse_rows(m) ((m)->n_rows)  /* Get the number of rows or columns*/
#define mod2sparse_cols(m) ((m)->n_cols)  /* in a matrix                      */


/* POSSIBLE LU DECOMPOSITION STRATEGIES.  For use with mod2sparse_decomp. */

typedef enum
{ Mod2sparse_first,
  Mod2sparse_mincol,
  Mod2sparse_minprod
} mod2sparse_strategy;

void mod2sparse_mulvec    (mod2sparse *m, char *u, char *v);
/*
PROCEDURES TO MANIPULATE SPARSE MATRICES.
//void mod2sparse_clear(mod2sparse *r);
//mod2sparse *mod2sparse_allocate(int n_rows,int n_cols);
//void mod2sparse_free(mod2sparse *);



//void mod2sparse_copy     (mod2sparse *, mod2sparse *);
//void mod2sparse_copyrows (mod2sparse *, mod2sparse *, int *);
//void mod2sparse_copycols (mod2sparse *, mod2sparse *, int *);

void mod2sparse_print       (FILE *, mod2sparse *);
int  mod2sparse_write       (FILE *, mod2sparse *);
//mod2sparse *mod2sparse_read (FILE *);

mod2entry *mod2sparse_find   (mod2sparse *, int, int);
//mod2entry *mod2sparse_insert (mod2sparse *, int, int);
//void mod2sparse_delete       (mod2sparse *, mod2entry *);

void mod2sparse_transpose (mod2sparse *, mod2sparse *);
void mod2sparse_add       (mod2sparse *, mod2sparse *, mod2sparse *);
void mod2sparse_multiply  (mod2sparse *, mod2sparse *, mod2sparse *);


int mod2sparse_equal (mod2sparse *, mod2sparse *);

//int mod2sparse_count_row (mod2sparse *, int);
//int mod2sparse_count_col (mod2sparse *, int);

//void mod2sparse_add_row (mod2sparse *, int, mod2sparse *, int);
void mod2sparse_add_col (mod2sparse *, int, mod2sparse *, int);

int mod2sparse_decomp (mod2sparse *, int, mod2sparse *, mod2sparse *,
                       int *, int *, mod2sparse_strategy, int, int);

int mod2sparse_forward_sub  (mod2sparse *, int *, char *, char *);
int mod2sparse_backward_sub (mod2sparse *, int *, char *, char *);
*/

//mod2dense

/* PACKING OF BITS INTO WORDS.  Bits are packed into 32-bit words, with
   the low-order bit coming first. */

typedef uint32_t mod2word;	/* Data type that holds packed bits. If uint32_t
		 		   doesn't exist, change it to unsigned long */

#define mod2_wordsize 32	/* Number of bits that fit in a mod2word. Can't
				   be increased without changing intio module */

#define mod2_wordsize_shift 5	/* Amount to shift by to divide by wordsize */
#define mod2_wordsize_mask 0x1f /* What to AND with to produce mod wordsize */

/* Extract the i'th bit of a mod2word. */

#define mod2_getbit(w,i) (((w)>>(i))&1)

/* Make a word like w, but with the i'th bit set to 1 (if it wasn't already). */

#define mod2_setbit1(w,i) ((w)|(1<<(i)))

/* Make a word like w, but with the i'th bit set to 0 (if it wasn't already). */

#define mod2_setbit0(w,i) ((w)&(~(1<<(i))))


/* STRUCTURE REPRESENTING A DENSE MATRIX.  These structures are dynamically
   allocated using mod2dense_allocate (or by other procedures that call
   mod2dense_allocate).  They should be freed with mod2dense_free when no
   longer required.

   Direct access to this structure should be avoided except in low-level
   routines.  Use the macros and procedures defined below instead. */

typedef struct
{
  int n_rows;		/* Number of rows in the matrix */
  int n_cols;		/* Number of columns in the matrix */

  int n_words;		/* Number of words used to store a column of bits */

  mod2word **col;	/* Pointer to array of pointers to columns */

  mod2word *bits;	/* Pointer to storage block for bits in this matrix
                           (pieces of this block are pointed to from col) */
} mod2dense;


/* MACROS. */

#define mod2dense_rows(m) ((m)->n_rows)  /* Get the number of rows or columns */
#define mod2dense_cols(m) ((m)->n_cols)  /* in a matrix                       */


/* PROCEDURES. */

/*mod2dense *mod2dense_allocate (int, int);
void mod2dense_free           (mod2dense *);

void mod2dense_clear    (mod2dense *);
void mod2dense_copy     (mod2dense *, mod2dense *);
void mod2dense_copyrows (mod2dense*, mod2dense *, int *);
void mod2dense_copycols (mod2dense*, mod2dense *, int *);

void mod2dense_print      (FILE *, mod2dense *);
int  mod2dense_write      (FILE *, mod2dense *);
mod2dense *mod2dense_read (FILE *);

int  mod2dense_get (mod2dense *, int, int);
void mod2dense_set (mod2dense *, int, int, int);
int  mod2dense_flip(mod2dense *, int, int);

void mod2dense_transpose (mod2dense *, mod2dense *);
void mod2dense_add       (mod2dense *, mod2dense *, mod2dense *);
void mod2dense_multiply  (mod2dense *, mod2dense *, mod2dense *);

int mod2dense_equal (mod2dense *, mod2dense *);

int mod2dense_invert          (mod2dense *, mod2dense *);
int mod2dense_forcibly_invert (mod2dense *, mod2dense *, int *, int *);
int mod2dense_invert_selected (mod2dense *, mod2dense *, int *, int *);
*/

//rcode

/* PROCEDURES FOR READING DATA. */

int read_pchk (int);
int read_gen  (int, int, int);

//rand
#define N_tables 5		/* Number of tables of real random numbers */

typedef struct
{ int seed;			/* Seed state derives from */
  int ptr[N_tables];		/* Pointers for tables of real random numbers */
  unsigned short state48[3];	/* State of 'rand48' pseudo-random generator */
} rand_state;


/* BASIC PSEUDO-RANDOM GENERATION PROCEDURES. */

void rand_seed (int);		/* Initialize current state structure by seed */

void rand_use_state (rand_state *); /* Start using given state structure */
rand_state *rand_get_state (void);  /* Return pointer to current state */

int rand_word (void);		/* Generate random 31-bit positive integer */


/* GENERATORS FOR VARIOUS DISTRIBUTIONS. */

double rand_uniform (void);	/* Uniform from [0,1) */
double rand_uniopen (void);	/* Uniform from (0,1) */

int rand_int (int);		/* Uniform from 0, 1, ... (n-1) */
int rand_pickd (double *, int);	/* From 0 ... (n-1), with given distribution */
int rand_pickf (float *, int);	/* Same as above, but with floats */
void rand_permutation (int *, int); /* Random permutation */

int rand_poisson (double);	/* Poisson with given mean */
double rand_gaussian (void);	/* Gaussian with mean zero and unit variance */
double rand_logistic (void);	/* Logistic centred at zero with unit width */
double rand_cauchy (void);	/* Cauchy centred at zero with unit width */
double rand_gamma (double);	/* Gamma with given shape parameter */
double rand_exp (void);		/* Exponential with mean one */
double rand_beta (double, double); /* Beta with given parameters */


//open
FILE *open_file_std (char *, char *);




//mode2convert
void mod2sparse_to_dense (mod2sparse *, mod2dense *);
void mod2dense_to_sparse (mod2dense *, mod2sparse *);


//intio
int  intio_read  (FILE *);	/* Read an integer */
void intio_write (FILE *, int);	/* Write an integer */


//enc
void sparse_encode (char *, char *);
void dense_encode  (char *, char *, mod2dense *, mod2dense *);
void mixed_encode  (char *, char *, mod2dense *, mod2dense *);


//distrib
typedef struct distrib_entry
{ int num;			/* A positive number */
  double prop;			/* Proportion for this number */
} distrib_entry;

typedef struct distrib
{ struct distrib_entry *list;	/* The list of numbers and proportions */
  int size;			/* Number of entries in the list */
} distrib;


/* MACROS TO ACCESS ELEMENTS OF A DISTRIBUTION LIST.  Note that indexes for
   entries start at 0. */

#define distrib_num(d,i) \
  ((d)->list[i].num)		/* The number for the i'th entry */

#define distrib_prop(d,i) \
  ((d)->list[i].prop)		/* The i'th entry's proportion [probability] */

#define distrib_size(d) \
  ((d)->size)			/* The length of the list (integer) */


/* PROCEDURES FOR DISTRIBUTION LISTS. */

distrib *distrib_create	(char *);
void distrib_free (distrib *);

int distrib_max(distrib *);


mod2sparse *H;		/* Parity check matrix */

int M;			/* Number of rows in parity check matrix */
int N;			/* Number of columns in parity check matrix */

char type;		/* Type of generator matrix representation (s/d/m) */
int *cols;		/* Ordering of columns in generator matrix */

mod2sparse *L, *U;	/* Sparse LU decomposition, if type=='s' */
int *rows;		/* Ordering of rows in generator matrix (type 's') */

mod2dense *G;		/* Dense or mixed representation of generator matrix,
			   if type=='d' or type=='m' */
         static long int this_nrand48 (unsigned short int xsubi[3]);
         #ifndef M_PI
         #define M_PI 3.14159265358979323846
         #endif

         #define Table_size 5000			/* Number of words in each table */

         static int rn[N_tables][Table_size];	/* Random number tables */

         static int initialized = 0;		/* Has module been initialized? */

         static rand_state state0;		/* Default state structure */

         static rand_state *state;		/* Pointer to current state */
//dec
typedef enum
{ Enum_block, Enum_bit, Prprp
} decoding_method;

extern decoding_method dec_method; /* Decoding method to use */

extern int table;	/* Trace option, 2 for a table of decoding details */
extern int block_no;	/* Number of current block, from zero */

extern int max_iter;	/* Maximum number of iteratons of decoding to do */
extern char *gen_file;	/* Generator file for Enum_block and Enum_bit */


/* PROCEDURES RELATING TO DECODING METHODS. */

int enum_decode_setup (int exten_no);
unsigned enum_decode (double *, char *, double *, int);

void prprp_decode_setup (void);
unsigned prprp_decode
(mod2sparse *, double *, char *, char *, double *);

void initprp (mod2sparse *, double *, char *, double *);
void iterprp (mod2sparse *, double *, char *, double *);


//check
int check (mod2sparse *, char *, char *);


double expected_parity_errors (mod2sparse *, double *);

double loglikelihood (double *, char *, int);
double expected_loglikelihood (double *, double *, int);

double entropy (double *, int);


//channel
typedef enum { BSC, AWGN, AWLN } channel_type;

extern channel_type channel;	/* Type of channel */

extern double std_dev;		/* Noise standard deviation for AWGN */
extern double lwidth;		/* Width of noise distributoin for AWLN */


/* PROCEDURES TO DO WITH CHANNELS. */

int  channel_parse (char **, int);
void channel_usage (void);


//blockio
extern int blockio_flush;

int  blockio_read  (FILE *, char *, int);
void blockio_write (FILE *, char *, int);


//alloc
void *chk_alloc (unsigned, unsigned);

extern mod2sparse *H;

extern int M;
extern int N;

extern char type;
extern int *cols;

extern mod2sparse *L, *U;
extern int *rows;
extern mod2dense *G;

void mod2dense_free
( mod2dense *m		/* Matrix to free */
)
{ free(m->bits);
  free(m->col);
  free(m);
}

mod2dense *mod2dense_allocate
( int n_rows, 		/* Number of rows in matrix */
  int n_cols		/* Number of columns in matrix */
)
{
  mod2dense *m;
  int j;

  if (n_rows<=0 || n_cols<=0)
  { fprintf(stderr,"mod2dense_allocate: Invalid number of rows or columns\n");
    exit(1);
  }

  m = (mod2dense *)chk_alloc (1, sizeof *m);

  m->n_rows = n_rows;
  m->n_cols = n_cols;

  m->n_words = (n_rows+mod2_wordsize-1) >> mod2_wordsize_shift;

  m->col = (mod2word **)chk_alloc (m->n_cols, sizeof *m->col);

  m->bits = (mod2word *)chk_alloc(m->n_words*m->n_cols, sizeof *m->bits);

  for (j = 0; j<m->n_cols; j++)
  { m->col[j] = m->bits + j*m->n_words;
  }

  return m;
}
mod2dense *mod2dense_read
( FILE *f
)
{
  int n_rows, n_cols;
  mod2dense *m;
  int j, k;

  n_rows = intio_read(f);
  if (feof(f) || ferror(f) || n_rows<=0) return 0;

  n_cols = intio_read(f);
  if (feof(f) || ferror(f) || n_cols<=0) return 0;

  m = mod2dense_allocate(n_rows,n_cols);

  for (j = 0; j<mod2dense_cols(m); j++)
  {
    for (k = 0; k<m->n_words; k++)
    { m->col[j][k] = intio_read(f);
      if (feof(f) || ferror(f))
      { mod2dense_free(m);
        return 0;
      }
    }
  }

  return m;
}
int check
( mod2sparse *H,	/* Parity check matrix */
  char *dblk,		/* Guess for codeword */
  char *pchk		/* Place to store parity checks */
)
{
  int M, i, c;

  M = mod2sparse_rows(H);

  mod2sparse_mulvec (H, dblk, pchk);

  c = 0;
  for (i = 0; i<M; i++)
  { c += pchk[i];
  }

  return c;
}
double changed( double *lratio,	/* Likelihood ratios for bits */char *dblk,		/* Candidate decoding */int N			/* Number of bits */)
{
  double changed;
  int j;
  changed = 0;
  for (j = 0; j<N; j++)
  {
	if(lratio[j] == 1)
	{
		changed += 0.5;
	}
	else
	{
		changed += (dblk[j] != (lratio[j]>1));
	}
  }
  return changed;
}
decoding_method dec_method;	/* Decoding method to use */
int table;	/* Trace option, 2 for a table of decoding details */
int block_no;	/* Number of current block, from zero */

int max_iter;	/* Maximum number of iteratons of decoding to do */
char *gen_file;	/* Generator file for Enum_block and Enum_bit */

int enum_decode_setup(int exten_no)
{
  int returnvalue = 0;
  if(exten_no == 0)
  {
    returnvalue = read_gen(0,0,0);
  }
  else if(exten_no == 1)
  {
    returnvalue = read_gen(0,0,1);
  }
  else if(exten_no == 2)
  {
    returnvalue = read_gen(0,0,2);
  }


  if (N-M>31)
  { fprintf(stderr,"Trying to decode messages with %d bits by exhaustive enumeration is absurd!\n", N-M);
    //exit(1);
  }

  if (table==2)
  { printf("  block   decoding  likelihood\n");
  }
  return returnvalue;
}

unsigned enum_decode
( double *lratio,	/* Likelihood ratios for bits */
  char *dblk, 		/* Place to stored decoded message */
  double *bitpr,	/* Place to store marginal bit probabilities */
  int max_block		/* Maximize probability of whole block being correct? */
)
{
  mod2dense *u, *v;
  double lk, maxlk, tpr;
  double *bpr, *lk0, *lk1;
  char sblk[33];
  char *cblk;
  unsigned d;
  int i, j;

  if (N-M>33) abort();

  /* Allocate needed space. */

  bpr = bitpr;
  if (bpr==0 && max_block==0)
  { bpr = (double *)chk_alloc (N, sizeof *bpr);
  }

  cblk = (char *)chk_alloc (N, sizeof *cblk);

  if (type=='d')
  { u = mod2dense_allocate(N-M,1);
    v = mod2dense_allocate(M,1);
  }

  if (type=='m')
  { u = mod2dense_allocate(M,1);
    v = mod2dense_allocate(M,1);
  }

  lk0 = (double *)chk_alloc (N, sizeof *lk0);
  lk1 = (double *)chk_alloc (N, sizeof *lk1);

  /* Pre-compute likelihoods for bits. */

  for (j = 0; j<N; j++)
  { lk0[j] = 1/(1+lratio[j]);
    lk1[j] = 1 - lk0[j];
  }

  /* Initialize marginal bit probabilities. */

  if (bpr)
  { for (j = 0; j<N; j++) bpr[j] = 0.0;
  }

  /* Exhaustively try all possible decoded messages. */

  tpr = 0.0;

  for (d = 0; d<=(1<<(N-M))-1; d++)
  {
    /* Unpack message into source block. */

    for (i = N-M-1; i>=0; i--)
    { sblk[i] = (d>>i)&1;
    }

    /* Find full codeword for this message. */

    switch (type)
    { case 's':
      { sparse_encode (sblk, cblk);
        break;
      }
      case 'd':
      { dense_encode (sblk, cblk, u, v);
        break;
      }
      case 'm':
      { mixed_encode (sblk, cblk, u, v);
        break;
      }
    }

    /* Compute likelihood for this decoding. */

    lk = 1;
    for (j = 0; j<N; j++)
    { lk *= cblk[j]==0 ? lk0[j] : lk1[j];
    }

    /* Update maximum likelihood decoding. */

    if (max_block)
    { if (d==0 || lk>maxlk)
      { for (j = 0; j<N; j++)
        { dblk[j] = cblk[j];
        }
        maxlk = lk;
      }
    }

    /* Update bit probabilities. */

    if (bpr)
    { for (j = 0; j<N; j++)
      { if (cblk[j]==1)
        { bpr[j] += lk;
        }
      }
      tpr += lk;
    }

    /* Output data to trace file. */

    if (table==2)
    { printf("%7d %10x  %10.4e\n",block_no,d,lk);
    }
  }

  /* Normalize bit probabilities. */

  if (bpr)
  { for (j = 0; j<N; j++) bpr[j] /= tpr;
  }

  /* Decoding to maximize bit-by-bit success, if that's what's wanted.
     In case of a tie, decode to a 1. */

  if (!max_block)
  { for (j = 0; j<N; j++)
    { dblk[j] = bpr[j]>=0.5;
    }
  }

  /* Free space. */

  if (bpr!=0 && bpr!=bitpr) free(bpr);
  free(cblk);
  free(lk0);
  free(lk1);
  return 1<<(N-M);
}
int Ldpc_decode(int exten_no)
{
  char *pchk_file, *rfile, *dfile, *pfile;
  char **meth;
  FILE *rf, *df, *pf;

  char *dblk, *pchk;
  double *lratio;
  double *bitpr;

  double *awn_data;		/* Places to store channel data */
  int *bsc_data;

  unsigned iters;		/* Unsigned because can be huge for enum */
  double tot_iter;		/* Double because can be huge for enum */
  double chngd, tot_changed;	/* Double because can be fraction if lratio==1*/

  int tot_valid;
  char junk;
  int valid;
  int return_pchk;
  int i, j, k;

  /* Look at initial flag arguments. */

  table = 0;
  blockio_flush = 0;

  /* Look at arguments up to the decoding method specification. */

  pfile = 0;
  dec_method = Enum_bit;

	printf("read_pchk\n");
  if(exten_no == 0)
  {
    return_pchk = read_pchk(0);
  }
  else if(exten_no == 1)
  {
    return_pchk = read_pchk(1);
  }
  else if(exten_no == 2)
  {
    return_pchk = read_pchk(2);
  }
  else if(return_pchk == 3)
  {
    return 3;
  }
  if (N<=M)
  { fprintf(stderr,
     "Number of bits (%d) should be greater than number of checks (%d)\n",N,M);
    exit(1);
  }
	//printf("load receive file\n");
  /* Open file of received data. */
  if(exten_no==0)
  {
    rf = open_file_std("read.rec","rb");
    if (rf==NULL)
    { fprintf(stderr,"Can't open file of received data: %s\n",rfile);
      exit(1);
    }
  }
  else if(exten_no==1)
  {
    rf = open_file_std("read_2.rec","rb");
    if (rf==NULL)
    { fprintf(stderr,"Can't open file of received data: %s\n",rfile);
      exit(1);
    }
  }
  else if(exten_no==2)
  {
    rf = open_file_std("read_3.rec","rb");
    if (rf==NULL)
    { fprintf(stderr,"Can't open file of received data: %s\n",rfile);
      exit(1);
    }
  }
	printf("make decode file\n");
  /* Create file for decoded data. */
  df = open_file_std("test.dec","w");

  if (df==NULL)
  { fprintf(stderr,"Can't create file for decoded data: %s\n",dfile);
    exit(1);
  }
  /* Allocate space for data from channel. */
  bsc_data = (int *)chk_alloc (N, sizeof *bsc_data);

  dblk   = (char *)chk_alloc (N, sizeof *dblk);
  lratio = (double *)chk_alloc (N, sizeof *lratio);
  pchk   = (char *)chk_alloc (M, sizeof *pchk);
  bitpr  = (double *)chk_alloc (N, sizeof *bitpr);

  /* Print header for summary table. */
  int returnvalue = 0;
  if(exten_no == 0)
  {
    returnvalue = enum_decode_setup(0);
  }
  else if(exten_no == 1)
  {
    returnvalue = enum_decode_setup(1);
  }
  else if(exten_no == 2)
  {
    returnvalue = enum_decode_setup(2);
  }
  if(returnvalue == 3)
  {
    return 3;
  }
  /* Read received blocks, decode, and write decoded blocks. */

  tot_iter = 0;
  tot_valid = 0;
  tot_changed = 0;
  for (block_no = 0; ; block_no++)
  {
    /* Read block from received file, exit if end-of-file encountered. */

    for (i = 0; i<N; i++)
    {
      char c;
      c = fscanf(rf,"%1d",&bsc_data[i]);
      if (c==EOF)
      {
        if (i>0)
        {
          fprintf(stderr, "Warning: Short block (%d long) at end of received file ignored\n",i);
        }
        goto done;
      }
      if (c<1 || bsc_data[i]!=0 && bsc_data[i]!=1)
      { fprintf(stderr,"File of received data is garbled\n");
        exit(1);
      }
    }

    /* Find likelihood ratio for each bit. */
    for (i = 0; i<N; i++)
    {
      lratio[i] = bsc_data[i]==1 ? (1-error_prob) / error_prob : error_prob / (1-error_prob);
    }
    /* Try to decode using the specified method. */
    iters = enum_decode (lratio, dblk, bitpr, dec_method==Enum_block);


    /* See if it worked, and how many bits were changed. */

    valid = check(H,dblk,pchk)==0;

    chngd = changed(lratio,dblk,N);
    tot_iter += iters;
    tot_valid += valid;
    tot_changed += chngd;

    /* Print summary table entry. */

    /* Write decoded block. */
    blockio_write(df,dblk,N);

    /* Write bit probabilities, if asked to. */
}
done:
  fprintf(stderr,
  "Decoded %d blocks, %d valid.  Average %.1f iterations, %.03f%% bit changes\n",
   block_no, tot_valid, (double)tot_iter/block_no,
   100.0*(double)tot_changed/(N*block_no));
   if(block_no == tot_valid)
   {
     return 1;
   }
   else
   {
     return 0;
   }
}
int Ldpc_encode(int extend_no)
{
  mod2dense *u, *v;
  //int extend_number = extend_no;
  FILE *srcf, *encf;
  char *sblk, *cblk, *chks;
  int i, n;
  blockio_flush = 0;
  int return_pchk;
  int return_gen;
  /* Read parity check file */
  return_pchk = read_pchk(extend_no);
  if(return_pchk == 3)
  {
    return 3;
  }

  if (N<=M)
  { fprintf(stderr,
 "Can't encode if number of bits (%d) not greater than number of checks (%d)\n",
      N,M);
    exit(1);
  }
  /* Read generator matrix file. */
  
  return_gen = read_gen(0,0,extend_no);
  if(return_gen == 3)
  {
    return 3;
  }

  /* Allocate needed space. */

  if (type=='d')
  { u = mod2dense_allocate(N-M,1);
    v = mod2dense_allocate(M,1);
  }

  if (type=='m')
  { u = mod2dense_allocate(M,1);
    v = mod2dense_allocate(M,1);
  }

  /* Open source file. */
  int num = rand()%3;
  if(num == 0)
  {
  	srcf = fopen("./src_1.src","r");
  }
  else if (num == 1)
  {
	srcf = fopen("./src_2.src","r");
  }
  else if(num==2)
  {
	srcf = fopen("./src_3.src","r");
  }
  /* Create encoded output file. */
  encf = fopen("./12_28.enc","w");

  sblk = (char *)chk_alloc (N-M, sizeof *sblk);
  cblk = (char *)chk_alloc (N, sizeof *cblk);
  chks = (char *)chk_alloc (M, sizeof *chks);

  /* Encode successive blocks. */
  for (n = 0; ; n++)
  {
    /* Read block from source file. */

    if (blockio_read(srcf,sblk,N-M)==EOF)
    {

    break;
    }

    /* Compute encoded block. */
    switch (type)
    { case 's':
      { sparse_encode (sblk, cblk);
        break;
      }
      case 'd':
      { dense_encode (sblk, cblk, u, v);
        break;
      }
      case 'm':
      { mixed_encode (sblk, cblk, u, v);
        break;
      }
    }
    /* Check that encoded block is a code word. */

    mod2sparse_mulvec (H, cblk, chks);

    for (i = 0; i<M; i++)
    { if (chks[i]==1)
      { fprintf(stderr,"Output block %d is not a code word!  (Fails check %d)\n",n,i);
        abort();
      }
    }

    /* Write encoded block to encoded output file. */

    blockio_write(encf,cblk,N);
    if (ferror(encf))
    { break;
    }
  }
  fclose(srcf);
  fclose(encf);
  fprintf(stderr,
    "Encoded %d blocks, source block size %d, encoded block size %d\n",n,N-M,N);
    return 1;


}
static mod2entry *alloc_entry( mod2sparse *m)
{
  mod2block *b;
  mod2entry *e;
  int k;

  if (m->next_free==0)
  {
    b = (mod2block *)chk_alloc (1, sizeof *b);

    b->next = m->blocks;
    m->blocks = b;

    for (k = 0; k<Mod2sparse_block; k++)
    { b->entry[k].left = m->next_free;
      m->next_free = &b->entry[k];
    }
  }

  e = m->next_free;
  m->next_free = e->left;

  e->pr = 0;
  e->lr = 0;

  return e;
}
void ldpc_transmit(double error)
{
  char *tfile, *rfile;
FILE *tf, *rf;
int block_size, n_bits;
char junk;
int seed;
int cnt;
int n, b;

/* Look at arguments.  The arguments specifying the channel are looked
   at by channel_parse in channel.c */

tfile = "12_28.enc";
rfile = "write.rec";
seed = 2;
error_prob = error;
junk = error_prob;
/* See if the source is all zeros or a file. */

if (sscanf(tfile,"%d%c",&n_bits,&junk)==1 && n_bits>0)
{ block_size = 1;
  tf = NULL;
}
else if (sscanf(tfile,"%dx%d%c",&block_size,&n_bits,&junk)==2
          && block_size>0 && n_bits>0)
{ n_bits *= block_size;
  tf = NULL;
}
else
{ tf = open_file_std(tfile,"r");
  if (tf==NULL)
  { fprintf(stderr,"Can't open encoded file to transmit: %s\n",tfile);
    exit(1);
  }
}

/* Open output file. */

rf = open_file_std(rfile,"w");
if (rf==NULL)
{ fprintf(stderr,"Can't create file for received data: %s\n",rfile);
  exit(1);
}

/* Set random seed to avoid duplications with other programs. */
  /* Transmit bits. */
  srand(time(NULL));
  int randnum;
  for (cnt = 0; ; cnt++)
  {

    randnum = rand()%1000000;
    /* Get next bit to transmit. */

    if (tf) /* Data comes from a file */
    {
      for (;;)
      { b = getc(tf);
        if (b!=' ' && b!='\t' && b!='\n' && b!='\r')
        { break;
        }
        putc(b,rf);
      }

      if (b==EOF) break;

      if (b!='0' && b!='1')
      { fprintf(stderr,"Bad character (code %d) file being transmitted\n",b);
        exit(1);
      }
    }

    else /* Data is all zeros */
    {
      if (cnt>0 && cnt%block_size==0)
      { putc('\n',rf);
      }

      if (cnt==n_bits) break;

      b = '0';
    }

    b = b=='1';
    int bsc_noise;
           //bsc_noise = rand_uniform() < error_prob;


         //  printf("%d ", randnum);
           if(randnum < (error_prob * 1000000))
           {
             bsc_noise = 1;
           }
           else
           {
             bsc_noise = 0;
           }

           fprintf (rf, "%d", b^bsc_noise);
 }
fclose(tf);
fclose(rf);
}

/* ALLOCATE SPACE FOR A SPARSE MOD2 MATRIX.  */

mod2sparse *mod2sparse_allocate(int n_rows,int n_cols)
{
  mod2sparse *m;
  mod2entry *e;
  int i, j;

  if (n_rows<=0 || n_cols<=0)
  { fprintf(stderr,"mod2sparse_allocate: Invalid number of rows or columns\n");
    exit(1);
  }

  m = (mod2sparse *)chk_alloc (1, sizeof *m);

  m->n_rows = n_rows;
  m->n_cols = n_cols;

  m->rows = (mod2entry *)chk_alloc (n_rows, sizeof *m->rows);
  m->cols = (mod2entry *)chk_alloc (n_cols, sizeof *m->cols);

  m->blocks = 0;
  m->next_free = 0;

  for (i = 0; i<n_rows; i++)
  { e = &m->rows[i];
    e->left = e->right = e->up = e->down = e;
    e->row = e->col = -1;
  }

  for (j = 0; j<n_cols; j++)
  { e = &m->cols[j];
    e->left = e->right = e->up = e->down = e;
    e->row = e->col = -1;
  }

  return m;
}
int mod2dense_get
( mod2dense *m, 	/* Matrix to get element from */
  int row,		/* Row of element (starting with zero) */
  int col		/* Column of element (starting with zero) */
)
{
  if (row<0 || row>=mod2dense_rows(m) || col<0 || col>=mod2dense_cols(m))
  { fprintf(stderr,"mod2dense_get: row or column index out of bounds\n");
    exit(1);
  }

  return mod2_getbit (m->col[col][row>>mod2_wordsize_shift],
                      row&mod2_wordsize_mask);
}

/* FREE SPACE OCCUPIED BY A SPARSE MOD2 MATRIX. */

void mod2sparse_free
( mod2sparse *m		/* Matrix to free */
)
{
  mod2block *b;

  free(m->rows);
  free(m->cols);

  while (m->blocks!=0)
  { b = m->blocks;
    m->blocks = b->next;
    free(b);
  }
}
mod2entry *mod2sparse_insert
( mod2sparse *m,
  int row,
  int col
)
{
  mod2entry *re, *ce, *ne;

  if (row<0 || row>=mod2sparse_rows(m) || col<0 || col>=mod2sparse_cols(m))
  { fprintf(stderr,"mod2sparse_insert: row or column index out of bounds\n");
    exit(1);
  }

  /* Find old entry and return it, or allocate new entry and insert into row. */

  re = mod2sparse_last_in_row(m,row);

  if (!mod2sparse_at_end(re) && mod2sparse_col(re)==col)
  { return re;
  }

  if (mod2sparse_at_end(re) || mod2sparse_col(re)<col)
  { re = re->right;
  }
  else
  {
    re = mod2sparse_first_in_row(m,row);

    for (;;)
    {
      if (!mod2sparse_at_end(re) && mod2sparse_col(re)==col)
      { return re;
      }

      if (mod2sparse_at_end(re) || mod2sparse_col(re)>col)
      { break;
      }

      re = mod2sparse_next_in_row(re);
    }
  }

  ne = alloc_entry(m);

  ne->row = row;
  ne->col = col;

  ne->left = re->left;
  ne->right = re;
  ne->left->right = ne;
  ne->right->left = ne;

  /* Insert new entry into column.  If we find an existing entry here,
     the matrix must be garbled, since we didn't find it in the row. */

  ce = mod2sparse_last_in_col(m,col);

  if (!mod2sparse_at_end(ce) && mod2sparse_row(ce)==row)
  { fprintf(stderr,"mod2sparse_insert: Garbled matrix\n");
    exit(1);
  }

  if (mod2sparse_at_end(ce) || mod2sparse_row(ce)<row)
  { ce = ce->down;
  }
  else
  {
    ce = mod2sparse_first_in_col(m,col);

    for (;;)
    {
      if (!mod2sparse_at_end(ce) && mod2sparse_row(ce)==row)
      { fprintf(stderr,"mod2sparse_insert: Garbled matrix\n");
        exit(1);
      }

      if (mod2sparse_at_end(ce) || mod2sparse_row(ce)>row)
      { break;
      }

      ce = mod2sparse_next_in_col(ce);
    }
  }

  ne->up = ce->up;
  ne->down = ce;
  ne->up->down = ne;
  ne->down->up = ne;

  /* Return the new entry. */

  return ne;
}
void mod2dense_set
( mod2dense *m, 	/* Matrix to modify element of */
  int row,		/* Row of element (starting with zero) */
  int col,		/* Column of element (starting with zero) */
  int value		/* New value of element (0 or 1) */
)
{
  mod2word *w;

  if (row<0 || row>=mod2dense_rows(m) || col<0 || col>=mod2dense_cols(m))
  { fprintf(stderr,"mod2dense_set: row or column index out of bounds\n");
    exit(1);
  }

  w = &m->col[col][row>>mod2_wordsize_shift];

  *w = value ? mod2_setbit1(*w,row&mod2_wordsize_mask)
             : mod2_setbit0(*w,row&mod2_wordsize_mask);
}
/* CLEAR A SPARSE MATRIX TO ALL ZEROS. */

void mod2sparse_clear(mod2sparse *r)
{
  mod2block *b;
  mod2entry *e;
  int i, j;

  for (i = 0; i<mod2sparse_rows(r); i++)
  { e = &r->rows[i];
    e->left = e->right = e->up = e->down = e;
  }

  for (j = 0; j<mod2sparse_cols(r); j++)
  { e = &r->cols[j];
    e->left = e->right = e->up = e->down = e;
  }

  while (r->blocks!=0)
  { b = r->blocks;
    r->blocks = b->next;
    free(b);
  }
}


/* COPY A SPARSE MATRIX. */

void mod2sparse_copy
( mod2sparse *m,	/* Matrix to copy */
  mod2sparse *r		/* Place to store copy of matrix */
)
{
  mod2entry *e, *f;
  int i;

  if (mod2sparse_rows(m)>mod2sparse_rows(r)
   || mod2sparse_cols(m)>mod2sparse_cols(r))
  { fprintf(stderr,"mod2sparse_copy: Destination matrix is too small\n");
    exit(1);
  }

  mod2sparse_clear(r);

  for (i = 0; i<mod2sparse_rows(m); i++)
  {
    e = mod2sparse_first_in_row(m,i);

    while (!mod2sparse_at_end(e))
    { f = mod2sparse_insert(r,e->row,e->col);
      f->lr = e->lr;
      f->pr = e->pr;
      e = mod2sparse_next_in_row(e);
    }
  }
}


/* COPY ROWS OF A SPARSE MOD2 MATRIX. */

void mod2sparse_copyrows
( mod2sparse *m,	/* Matrix to copy */
  mod2sparse *r,	/* Place to store copy of matrix */
  int *rows		/* Indexes of rows to copy, from 0 */
)
{
  mod2entry *e;
  int i;

  if (mod2sparse_cols(m)>mod2sparse_cols(r))
  { fprintf(stderr,
     "mod2sparse_copyrows: Destination matrix has fewer columns than source\n");
    exit(1);
  }

  mod2sparse_clear(r);

  for (i = 0; i<mod2sparse_rows(r); i++)
  { if (rows[i]<0 || rows[i]>=mod2sparse_rows(m))
    { fprintf(stderr,"mod2sparse_copyrows: Row index out of range\n");
      exit(1);
    }
    e = mod2sparse_first_in_row(m,rows[i]);
    while (!mod2sparse_at_end(e))
    { mod2sparse_insert(r,i,e->col);
      e = mod2sparse_next_in_row(e);
    }
  }
}


/* COPY COLUMNS OF A SPARSE MOD2 MATRIX. */

void mod2sparse_copycols
( mod2sparse *m,	/* Matrix to copy */
  mod2sparse *r,	/* Place to store copy of matrix */
  int *cols		/* Indexes of columns to copy, from 0 */
)
{
  mod2entry *e;
  int j;

  if (mod2sparse_rows(m)>mod2sparse_rows(r))
  { fprintf(stderr,
      "mod2sparse_copycols: Destination matrix has fewer rows than source\n");
    exit(1);
  }

  mod2sparse_clear(r);

  for (j = 0; j<mod2sparse_cols(r); j++)
  { if (cols[j]<0 || cols[j]>=mod2sparse_cols(m))
    { fprintf(stderr,"mod2sparse_copycols: Column index out of range\n");
      exit(1);
    }
    e = mod2sparse_first_in_col(m,cols[j]);
    while (!mod2sparse_at_end(e))
    { mod2sparse_insert(r,e->row,j);
      e = mod2sparse_next_in_col(e);
    }
  }
}


/* PRINT A SPARSE MOD2 MATRIX IN HUMAN-READABLE FORM. */

void mod2sparse_print
( FILE *f,
  mod2sparse *m
)
{
  int rdigits, cdigits;
  mod2entry *e;
  int i;

  rdigits = mod2sparse_rows(m)<=10 ? 1
          : mod2sparse_rows(m)<=100 ? 2
          : mod2sparse_rows(m)<=1000 ? 3
          : mod2sparse_rows(m)<=10000 ? 4
          : mod2sparse_rows(m)<=100000 ? 5
          : 6;

  cdigits = mod2sparse_cols(m)<=10 ? 1
          : mod2sparse_cols(m)<=100 ? 2
          : mod2sparse_cols(m)<=1000 ? 3
          : mod2sparse_cols(m)<=10000 ? 4
          : mod2sparse_cols(m)<=100000 ? 5
          : 6;

  for (i = 0; i<mod2sparse_rows(m); i++)
  {
    fprintf(f,"%*d:",rdigits,i);

    e = mod2sparse_first_in_row(m,i);
    while (!mod2sparse_at_end(e))
    { fprintf(f," %*d",cdigits,mod2sparse_col(e));
      e = mod2sparse_next_in_row(e);
    }

    fprintf(f,"\n");
  }
}


/* WRITE A SPARSE MOD2 MATRIX TO A FILE IN MACHINE-READABLE FORM. */

int mod2sparse_write
( FILE *f,
  mod2sparse *m
)
{
  mod2entry *e;
  int i;

  intio_write(f,m->n_rows);
  if (ferror(f)) return 0;

  intio_write(f,m->n_cols);
  if (ferror(f)) return 0;

  for (i = 0; i<mod2sparse_rows(m); i++)
  {
    e = mod2sparse_first_in_row(m,i);

    if (!mod2sparse_at_end(e))
    {
      intio_write (f, -(i+1));
      if (ferror(f)) return 0;

      while (!mod2sparse_at_end(e))
      {
        intio_write (f, mod2sparse_col(e)+1);
        if (ferror(f)) return 0;

        e = mod2sparse_next_in_row(e);
      }
    }
  }

  intio_write(f,0);
  if (ferror(f)) return 0;

  return 1;
}


/* READ A SPARSE MOD2 MATRIX STORED IN MACHINE-READABLE FORM FROM A FILE. */

mod2sparse *mod2sparse_read( FILE *f)
{

  int n_rows, n_cols;
  mod2sparse *m;
  int v, row, col;
  n_rows = intio_read(f);

  if (feof(f) || ferror(f) || n_rows<=0) return 0;
  n_cols = intio_read(f);

  if (feof(f) || ferror(f) || n_cols<=0) return 0;
  m = mod2sparse_allocate(n_rows,n_cols);

  row = -1;
  for (;;)
  {
    v = intio_read(f);

    if (feof(f) || ferror(f)) break;

    if (v==0)
    {

	return m;
    }
    else if (v<0)
    {
	row = -v-1;
      if (row>=n_rows) break;
    }
    else
    {
	 col = v-1;
      if (col>=n_cols) break;
      if (row==-1) break;
      mod2sparse_insert(m,row,col);
    }

  }

  /* Error if we get here. */

  mod2sparse_free(m);
  return 0;
}


/* LOOK FOR AN ENTRY WITH GIVEN ROW AND COLUMN. */

mod2entry *mod2sparse_find
( mod2sparse *m,
  int row,
  int col
)
{
  mod2entry *re, *ce;

  if (row<0 || row>=mod2sparse_rows(m) || col<0 || col>=mod2sparse_cols(m))
  { fprintf(stderr,"mod2sparse_find: row or column index out of bounds\n");
    exit(1);
  }

  /* Check last entries in row and column. */

  re = mod2sparse_last_in_row(m,row);
  if (mod2sparse_at_end(re) || mod2sparse_col(re)<col)
  { return 0;
  }
  if (mod2sparse_col(re)==col)
  { return re;
  }

  ce = mod2sparse_last_in_col(m,col);
  if (mod2sparse_at_end(ce) || mod2sparse_row(ce)<row)
  { return 0;
  }
  if (mod2sparse_row(ce)==row)
  { return ce;
  }

  /* Search row and column in parallel, from the front. */

  re = mod2sparse_first_in_row(m,row);
  ce = mod2sparse_first_in_col(m,col);

  for (;;)
  {
    if (mod2sparse_at_end(re) || mod2sparse_col(re)>col)
    { return 0;
    }
    if (mod2sparse_col(re)==col)
    { return re;
    }

    if (mod2sparse_at_end(ce) || mod2sparse_row(ce)>row)
    { return 0;
    }
    if (mod2sparse_row(ce)==row)
    { return ce;
    }

    re = mod2sparse_next_in_row(re);
    ce = mod2sparse_next_in_col(ce);
  }
}


/* INSERT AN ENTRY WITH GIVEN ROW AND COLUMN. */




/* DELETE AN ENTRY FROM A SPARSE MATRIX. */

void mod2sparse_delete
( mod2sparse *m,
  mod2entry *e
)
{
  if (e==0)
  { fprintf(stderr,"mod2sparse_delete: Trying to delete a null entry\n");
    exit(1);
  }

  if (e->row<0 || e->col<0)
  { fprintf(stderr,"mod2sparse_delete: Trying to delete a header entry\n");
    exit(1);
  }

  e->left->right = e->right;
  e->right->left = e->left;

  e->up->down = e->down;
  e->down->up = e->up;

  e->left = m->next_free;
  m->next_free = e;
}


/* TEST WHETHER TWO SPARSE MATRICES ARE EQUAL. */

int mod2sparse_equal
( mod2sparse *m1,
  mod2sparse *m2
)
{
  mod2entry *e1, *e2;
  int i;

  if (mod2sparse_rows(m1)!=mod2sparse_rows(m2)
   || mod2sparse_cols(m1)!=mod2sparse_cols(m2))
  { fprintf(stderr,"mod2sparse_equal: Matrices have different dimensions\n");
    exit(1);
  }

  for (i = 0; i<mod2sparse_rows(m1); i++)
  {
    e1 = mod2sparse_first_in_row(m1,i);
    e2 = mod2sparse_first_in_row(m2,i);

    while (!mod2sparse_at_end(e1) && !mod2sparse_at_end(e2))
    {
      if (mod2sparse_col(e1)!=mod2sparse_col(e2))
      { return 0;
      }

      e1 = mod2sparse_next_in_row(e1);
      e2 = mod2sparse_next_in_row(e2);
    }

    if (!mod2sparse_at_end(e1) || !mod2sparse_at_end(e2))
    { return 0;
    }
  }

  return 1;
}


/* COMPUTE THE TRANSPOSE OF A SPARSE MOD2 MATRIX. */

void mod2sparse_transpose
( mod2sparse *m,	/* Matrix to compute transpose of (left unchanged) */
  mod2sparse *r		/* Result of transpose operation */
)
{
  mod2entry *e;
  int i;

  if (mod2sparse_rows(m)!=mod2sparse_cols(r)
   || mod2sparse_cols(m)!=mod2sparse_rows(r))
  { fprintf(stderr,
     "mod2sparse_transpose: Matrices have incompatible dimensions\n");
    exit(1);
  }

  if (r==m)
  { fprintf(stderr,
     "mod2sparse_transpose: Result matrix is the same as the operand\n");
    exit(1);
  }

  mod2sparse_clear(r);

  for (i = 0; i<mod2sparse_rows(m); i++)
  {
    e = mod2sparse_first_in_row(m,i);

    while (!mod2sparse_at_end(e))
    { mod2sparse_insert(r,mod2sparse_col(e),i);
      e = mod2sparse_next_in_row(e);
    }
  }
}


/* ADD TWO SPARSE MOD2 MATRICES. */

void mod2sparse_add
( mod2sparse *m1,	/* Left operand of add */
  mod2sparse *m2,	/* Right operand of add */
  mod2sparse *r		/* Place to store result of add */
)
{
  mod2entry *e1, *e2;
  int i;

  if (mod2sparse_rows(m1)!=mod2sparse_rows(r)
   || mod2sparse_cols(m1)!=mod2sparse_cols(r)
   || mod2sparse_rows(m2)!=mod2sparse_rows(r)
   || mod2sparse_cols(m2)!=mod2sparse_cols(r))
  { fprintf(stderr,"mod2sparse_add: Matrices have different dimensions\n");
    exit(1);
  }

  if (r==m1 || r==m2)
  { fprintf(stderr,
     "mod2sparse_add: Result matrix is the same as one of the operands\n");
    exit(1);
  }

  mod2sparse_clear(r);

  for (i = 0; i<mod2sparse_rows(r); i++)
  {
    e1 = mod2sparse_first_in_row(m1,i);
    e2 = mod2sparse_first_in_row(m2,i);

    while (!mod2sparse_at_end(e1) && !mod2sparse_at_end(e2))
    {
      if (mod2sparse_col(e1)==mod2sparse_col(e2))
      { e1 = mod2sparse_next_in_row(e1);
        e2 = mod2sparse_next_in_row(e2);
      }

      else if (mod2sparse_col(e1)<mod2sparse_col(e2))
      {
	mod2sparse_insert(r,i,mod2sparse_col(e1));
        e1 = mod2sparse_next_in_row(e1);
      }

      else
      { mod2sparse_insert(r,i,mod2sparse_col(e2));
        e2 = mod2sparse_next_in_row(e2);
      }
    }

    while (!mod2sparse_at_end(e1))
    { mod2sparse_insert(r,i,mod2sparse_col(e1));
      e1 = mod2sparse_next_in_row(e1);
    }

    while (!mod2sparse_at_end(e2))
    { mod2sparse_insert(r,i,mod2sparse_col(e2));
      e2 = mod2sparse_next_in_row(e2);
    }
  }
}


/* MULTIPLY TWO SPARSE MOD2 MATRICES. */

void mod2sparse_multiply
( mod2sparse *m1, 	/* Left operand of multiply */
  mod2sparse *m2,	/* Right operand of multiply */
  mod2sparse *r		/* Place to store result of multiply */
)
{
  mod2entry *e1, *e2;
  int i, j, b;

  if (mod2sparse_cols(m1)!=mod2sparse_rows(m2)
   || mod2sparse_rows(m1)!=mod2sparse_rows(r)
   || mod2sparse_cols(m2)!=mod2sparse_cols(r))
  { fprintf (stderr,
      "mod2sparse_multiply: Matrices have incompatible dimensions\n");
    exit(1);
  }

  if (r==m1 || r==m2)
  { fprintf(stderr,
     "mod2sparse_multiply: Result matrix is the same as one of the operands\n");
    exit(1);
  }

  mod2sparse_clear(r);

  for (i = 0; i<mod2sparse_rows(m1); i++)
  {
    if (mod2sparse_at_end(mod2sparse_first_in_row(m1,i)))
    { continue;
    }

    for (j = 0; j<mod2sparse_cols(m2); j++)
    {
      b = 0;

      e1 = mod2sparse_first_in_row(m1,i);
      e2 = mod2sparse_first_in_col(m2,j);

      while (!mod2sparse_at_end(e1) && !mod2sparse_at_end(e2))
      {
        if (mod2sparse_col(e1)==mod2sparse_row(e2))
        { b ^= 1;
          e1 = mod2sparse_next_in_row(e1);
          e2 = mod2sparse_next_in_col(e2);
        }

        else if (mod2sparse_col(e1)<mod2sparse_row(e2))
        { e1 = mod2sparse_next_in_row(e1);
        }

        else
        { e2 = mod2sparse_next_in_col(e2);
        }
      }

      if (b)
      { mod2sparse_insert(r,i,j);
      }
    }
  }
}


/* MULTIPLY VECTOR BY SPARSE MATRIX. */

void mod2sparse_mulvec
( mod2sparse *m,	/* The sparse matrix, with M rows and N columns */
  char *u,		/* The input vector, N long */
  char *v		/* Place to store the result, M long */
)
{
  mod2entry *e;
  int M, N;
  int i, j;

  M = mod2sparse_rows(m);
  N = mod2sparse_cols(m);

  for (i = 0; i<M; i++) v[i] = 0;

  for (j = 0; j<N; j++)
  { if (u[j])
    { for (e = mod2sparse_first_in_col(m,j);
           !mod2sparse_at_end(e);
           e = mod2sparse_next_in_col(e))
      { v[mod2sparse_row(e)] ^= 1;
      }
    }
  }
}


/* COUNT ENTRIES IN A ROW. */

int mod2sparse_count_row
( mod2sparse *m,
  int row
)
{
  mod2entry *e;
  int count;

  if (row<0 || row>=mod2sparse_rows(m))
  { fprintf(stderr,"mod2sparse_count_row: row index out of bounds\n");
    exit(1);
  }

  count = 0;

  for (e = mod2sparse_first_in_row(m,row);
       !mod2sparse_at_end(e);
       e = mod2sparse_next_in_row(e))
  { count += 1;
  }

  return count;
}


/* COUNT ENTRIES IN A COLUMN. */

int mod2sparse_count_col
( mod2sparse *m,
  int col
)
{
  mod2entry *e;
  int count;

  if (col<0 || col>=mod2sparse_cols(m))
  { fprintf(stderr,"mod2sparse_count_col: column index out of bounds\n");
    exit(1);
  }

  count = 0;

  for (e = mod2sparse_first_in_col(m,col);
       !mod2sparse_at_end(e);
       e = mod2sparse_next_in_col(e))
  { count += 1;
  }

  return count;
}


/* ADD TO A ROW. */

void mod2sparse_add_row
( mod2sparse *m1,	/* Matrix containing row to add to */
  int row1,		/* Index in this matrix of row to add to */
  mod2sparse *m2,	/* Matrix containing row to add from */
  int row2		/* Index in this matrix of row to add from */
)
{
  mod2entry *f1, *f2, *ft;

  if (mod2sparse_cols(m1)<mod2sparse_cols(m2))
  { fprintf (stderr,
     "mod2sparse_add_row: row added to is shorter than row added from\n");
    exit(1);
  }

  if (row1<0 || row1>=mod2sparse_rows(m1)
   || row2<0 || row2>=mod2sparse_rows(m2))
  { fprintf (stderr,"mod2sparse_add_row: row index out of range\n");
    exit(1);
  }

  f1 = mod2sparse_first_in_row(m1,row1);
  f2 = mod2sparse_first_in_row(m2,row2);

  while (!mod2sparse_at_end(f1) && !mod2sparse_at_end(f2))
  { if (mod2sparse_col(f1)>mod2sparse_col(f2))
    { mod2sparse_insert(m1,row1,mod2sparse_col(f2));
      f2 = mod2sparse_next_in_row(f2);
    }
    else
    { ft = mod2sparse_next_in_row(f1);
      if (mod2sparse_col(f1)==mod2sparse_col(f2))
      { mod2sparse_delete(m1,f1);
        f2 = mod2sparse_next_in_row(f2);
      }
      f1 = ft;
    }
  }

  while (!mod2sparse_at_end(f2))
  { mod2sparse_insert(m1,row1,mod2sparse_col(f2));
    f2 = mod2sparse_next_in_row(f2);
  }
}


/* ADD TO A COLUMN. */

void mod2sparse_add_col
( mod2sparse *m1,	/* Matrix containing column to add to */
  int col1,		/* Index in this matrix of column to add to */
  mod2sparse *m2,	/* Matrix containing column to add from */
  int col2		/* Index in this matrix of column to add from */
)
{
  mod2entry *f1, *f2, *ft;

  if (mod2sparse_rows(m1)<mod2sparse_rows(m2))
  { fprintf (stderr,
     "mod2sparse_add_col: Column added to is shorter than column added from\n");
    exit(1);
  }

  if (col1<0 || col1>=mod2sparse_cols(m1)
   || col2<0 || col2>=mod2sparse_cols(m2))
  { fprintf (stderr,"mod2sparse_add_col: Column index out of range\n");
    exit(1);
  }

  f1 = mod2sparse_first_in_col(m1,col1);
  f2 = mod2sparse_first_in_col(m2,col2);

  while (!mod2sparse_at_end(f1) && !mod2sparse_at_end(f2))
  { if (mod2sparse_row(f1)>mod2sparse_row(f2))
    { mod2sparse_insert(m1,mod2sparse_row(f2),col1);
      f2 = mod2sparse_next_in_col(f2);
    }
    else
    { ft = mod2sparse_next_in_col(f1);
      if (mod2sparse_row(f1)==mod2sparse_row(f2))
      { mod2sparse_delete(m1,f1);
        f2 = mod2sparse_next_in_col(f2);
      }
      f1 = ft;
    }
  }

  while (!mod2sparse_at_end(f2))
  { mod2sparse_insert(m1,mod2sparse_row(f2),col1);
    f2 = mod2sparse_next_in_col(f2);
  }
}


/* FIND AN LU DECOMPOSITION OF A SPARSE MATRIX. */



/* SOLVE A LOWER-TRIANGULAR SYSTEM BY FORWARD SUBSTITUTION. */

int mod2sparse_forward_sub
( mod2sparse *L,	/* Matrix that is lower triangular after reordering */
  int *rows,		/* Array of indexes (from 0) of rows for new order */
  char *x,		/* Vector on right of equation, also reordered */
  char *y		/* Place to store solution */
)
{
  int K, i, j, ii, b, d;
  mod2entry *e;

  K = mod2sparse_cols(L);

  /* Make sure that L is lower-triangular, after row re-ordering. */

  for (i = 0; i<K; i++)
  { ii = rows ? rows[i] : i;
    e = mod2sparse_last_in_row(L,ii);
    if (!mod2sparse_at_end(e) && mod2sparse_col(e)>i)
    { fprintf(stderr,
        "mod2sparse_forward_sub: Matrix is not lower-triangular\n");
      exit(1);
    }
  }

  /* Solve system by forward substitution. */

  for (i = 0; i<K; i++)
  {
    ii = rows ? rows[i] : i;

    /* Look at bits in this row, forming inner product with partial
       solution, and seeing if the diagonal is 1. */

    d = 0;
    b = 0;

    for (e = mod2sparse_first_in_row(L,ii);
         !mod2sparse_at_end(e);
         e = mod2sparse_next_in_row(e))
    {
      j = mod2sparse_col(e);

      if (j==i)
      { d = 1;
      }
      else
      { b ^= y[j];
      }
    }

    /* Check for no solution if the diagonal isn't 1. */

    if (!d && b!=x[ii])
    { return 0;
    }

    /* Set bit of solution, zero if arbitrary. */

    y[i] = b^x[ii];
  }

  return 1;
}


/* SOLVE AN UPPER-TRIANGULAR SYSTEM BY BACKWARD SUBSTITUTION. */

int mod2sparse_backward_sub
( mod2sparse *U,	/* Matrix that is upper triangular after reordering */
  int *cols,		/* Array of indexes (from 0) of columns for new order */
  char *y,		/* Vector on right of equation */
  char *z		/* Place to store solution, also reordered */
)
{
  int K, i, j, ii, b, d;
  mod2entry *e;

  K = mod2sparse_rows(U);

  /* Make sure that U is upper-triangular, after column re-ordering. */

  for (i = 0; i<K; i++)
  { ii = cols ? cols[i] : i;
    e = mod2sparse_last_in_col(U,ii);
    if (!mod2sparse_at_end(e) && mod2sparse_row(e)>i)
    { fprintf(stderr,
        "mod2sparse_backward_sub: Matrix is not upper-triangular\n");
      exit(1);
    }
  }

  /* Solve system by backward substitution. */

  for (i = K-1; i>=0; i--)
  {
    ii = cols ? cols[i] : i;

    /* Look at bits in this row, forming inner product with partial
       solution, and seeing if the diagonal is 1. */

    d = 0;
    b = 0;

    for (e = mod2sparse_first_in_row(U,i);
         !mod2sparse_at_end(e);
         e = mod2sparse_next_in_row(e))
    {
      j = mod2sparse_col(e);

      if (j==ii)
      { d = 1;
      }
      else
      { b ^= z[j];
      }
    }

    /* Check for no solution if the diagonal isn't 1. */

    if (!d && b!=y[i])
    { return 0;
    }

    /* Set bit of solution, zero if arbitrary. */

    z[ii] = b^y[i];
  }

  return 1;
}


/* READ PARITY CHECK MATRIX.  Sets the H, M, and N global variables.  If an
   error is encountered, a message is displayed on standard error, and the
   program is terminated. */

int read_pchk(int extend_no)
{
  FILE *fr;
  if(extend_no == 0)
  {
  	fr = open_file_std("4_20.pchk","rb");
  }
  else if(extend_no == 1)
  {
  	fr = open_file_std("8_24.pchk","rb");
  }
  else if(extend_no == 2)
  {
    fr = open_file_std("12_28.pchk", "rb");
  }

  if (fr==NULL)
  { fprintf(stderr,"Can't open parity check file\n");
    return 3;
  }

  if (intio_read(fr)!=('P'<<8)+0x80)
  { fprintf(stderr,"File doesn't contain a parity check matrix\n");
    exit(1);
  }

  H = mod2sparse_read(fr);

  if (H==0)
  { fprintf(stderr,"Error reading parity check matrix from\n");
    exit(1);
  }

  M = mod2sparse_rows(H);
  N = mod2sparse_cols(H);

  fclose(fr);
}



/* READ GENERATOR MATRIX.  The parity check matrix must have already been
   read, unless the last argument is set to 1.  The generator matrix must be
   compatible with the parity check matrix, if it has been read.  If the
   second argument is 1, only the column ordering (the last N-M of which are
   the indexes of the message bits) is read, into the 'cols' global variable.
   Otherwise, everything is read, into the global variables appropriate
   to the representation.  The 'type' global variable is set to a letter
   indicating which represention is used.

   If an error is encountered, a message is displayed on standard error,
   and the program is terminated. */

int read_gen
(
  int cols_only,	/* Read only column ordering? */
  int no_pchk_file,	/* No parity check file used? */
  int extend_no
)
{
  int M2, N2;
  FILE *f;
  int i;
  if(extend_no == 0)
  {
	  f = open_file_std("4_20.gen","rb");
  }
  else if(extend_no == 1)
  {
		f = open_file_std("8_24.gen","rb");
    printf("open 8_24.gen\n");
  }
  else if(extend_no == 2)
  {
    f = open_file_std("12_28.gen", "rb");
    printf("open 12_28.gen\n");
  }
  

  if (f==NULL)
  { 
    fprintf(stderr,"Can't open generator matrix file: %s\n","16_32.gen");
    return 3;
  }

  if (intio_read(f)!=('G'<<8)+0x80)
  { fprintf(stderr,"File %s doesn't contain a generator matrix\n","16_32.gen");
    exit(1);
  }

  if (fread (&type, 1, 1, f) != 1) goto error;

  M2 = intio_read(f);
  N2 = intio_read(f);
  if (feof(f) || ferror(f)) goto error;

  if (no_pchk_file)
  { M = M2;
    N = N2;
  }
  else
  { if (M2!=M || N2!=N)
    { fprintf(stderr,
              "Generator matrix and parity-check matrix are incompatible\n");
      exit(1);
    }
  }

  cols = (int *)chk_alloc (N, sizeof *cols);
  rows = (int *)chk_alloc (M, sizeof *rows);

  for (i = 0; i<N; i++)
  { cols[i] = intio_read(f);
    if (feof(f) || ferror(f)) goto error;
  }

  if (!cols_only)
  {
    switch (type)
    {
      case 's':
      {
        for (i = 0; i<M; i++)
        { rows[i] = intio_read(f);
          if (feof(f) || ferror(f)) goto error;
        }

        if ((L = mod2sparse_read(f)) == 0) goto error;
        if ((U = mod2sparse_read(f)) == 0) goto error;

        if (mod2sparse_rows(L)!=M || mod2sparse_cols(L)!=M) goto garbled;
        if (mod2sparse_rows(U)!=M || mod2sparse_cols(U)<M) goto garbled;

        break;
      }

      case 'd':
      {
        if ((G = mod2dense_read(f)) == 0) goto error;

        if (mod2dense_rows(G)!=M || mod2dense_cols(G)!=N-M) goto garbled;

        break;
      }

      case 'm':
      {
        if ((G = mod2dense_read(f)) == 0) goto error;

        if (mod2dense_rows(G)!=M || mod2dense_cols(G)!=M) goto garbled;

        break;
      }

      default:
      { fprintf(stderr,
         "Unknown type of generator matrix in file \n");
        exit(1);
      }
    }
  }
  fclose(f);
  return 1;

error:
  fprintf(stderr,"Error reading generator matrix from file \n");
  exit(1);

garbled:
  fprintf(stderr,"Garbled generator matrix in file \n");
  exit(1);
}

void mod2dense_clear
( mod2dense *r
)
{
  int k, j;

  for (j = 0; j<mod2dense_cols(r); j++)
  { for (k = 0; k<r->n_words; k++)
    { r->col[j][k] = 0;
    }
  }

}


/* COPY A DENSE MOD2 MATRIX. */

void mod2dense_copy
( mod2dense *m,		/* Matrix to copy */
  mod2dense *r		/* Place to store copy of matrix */
)
{
  int k, j;

  if (mod2dense_rows(m)>mod2dense_rows(r)
   || mod2dense_cols(m)>mod2dense_cols(r))
  { fprintf(stderr,"mod2dense_copy: Destination matrix is too small\n");
    exit(1);
  }

  for (j = 0; j<mod2dense_cols(m); j++)
  { for (k = 0; k<m->n_words; k++)
    { r->col[j][k] = m->col[j][k];
    }
    for ( ; k<r->n_words; k++)
    { r->col[j][k] = 0;
    }
  }

  for ( ; j<mod2dense_cols(r); j++)
  { for (k = 0; k<r->n_words; k++)
    { r->col[j][k] = 0;
    }
  }
}


/* COPY ROWS OF A DENSE MOD2 MATRIX. */

void mod2dense_copyrows
( mod2dense *m,		/* Matrix to copy */
  mod2dense *r,		/* Place to store copy of matrix */
  int *rows		/* Indexes of rows to copy, from 0 */
)
{
  int i, j;

  if (mod2dense_cols(m)>mod2dense_cols(r))
  { fprintf(stderr,
      "mod2dense_copyrows: Destination matrix has fewer columns than source\n");
    exit(1);
  }

  mod2dense_clear(r);

  for (i = 0; i<mod2dense_rows(r); i++)
  { if (rows[i]<0 || rows[i]>=mod2dense_rows(m))
    { fprintf(stderr,"mod2dense_copyrows: Row index out of range\n");
      exit(1);
    }
    for (j = 0; j<mod2dense_cols(m); j++)
    { mod2dense_set(r,i,j,mod2dense_get(m,rows[i],j));
    }
  }
}


/* COPY COLUMNS OF A DENSE MOD2 MATRIX. */

void mod2dense_copycols
( mod2dense *m,		/* Matrix to copy */
  mod2dense *r,		/* Place to store copy of matrix */
  int *cols		/* Indexes of columns to copy, from 0 */
)
{
  int k, j;

  if (mod2dense_rows(m)>mod2dense_rows(r))
  { fprintf(stderr,
      "mod2dense_copycols: Destination matrix has fewer rows than source\n");
    exit(1);
  }

  for (j = 0; j<mod2dense_cols(r); j++)
  { if (cols[j]<0 || cols[j]>=mod2dense_cols(m))
    { fprintf(stderr,"mod2dense_copycols: Column index out of range\n");
      exit(1);
    }
    for (k = 0; k<m->n_words; k++)
    { r->col[j][k] = m->col[cols[j]][k];
    }
    for ( ; k<r->n_words; k++)
    { r->col[j][k] = 0;
    }
  }
}


/* PRINT A DENSE MOD2 MATRIX IN HUMAN-READABLE FORM. */

void mod2dense_print
( FILE *f,
  mod2dense *m
)
{
  int i, j;

  for (i = 0; i<mod2dense_rows(m); i++)
  { for (j = 0; j<mod2dense_cols(m); j++)
    { fprintf(f," %d",mod2dense_get(m,i,j));
    }
    fprintf(f,"\n");
  }
}


/* WRITE A DENSE MOD2 MATRIX TO A FILE IN MACHINE-READABLE FORM.

   Data is written using intio_write, so that it will be readable on a machine
   with a different byte-ordering.  At present, this assumes that the words
   used to pack bits into are no longer than 32 bits. */

int mod2dense_write
( FILE *f,
  mod2dense *m
)
{
  int j, k;

  intio_write(f,m->n_rows);
  if (ferror(f)) return 0;

  intio_write(f,m->n_cols);
  if (ferror(f)) return 0;

  for (j = 0; j<mod2dense_cols(m); j++)
  {
    for (k = 0; k<m->n_words; k++)
    { intio_write(f,m->col[j][k]);
      if (ferror(f)) return 0;
    }
  }

  return 1;
}


/* READ A DENSE MOD2 MATRIX STORED IN MACHINE-READABLE FORM FROM A FILE. */



/* GET AN ELEMENT FROM A DENSE MOD2 MATRIX. */




/* SET AN ELEMENT IN A DENSE MOD2 MATRIX. */




/* FLIP AN ELEMENT OF A DENSE MOD2 MATRIX. */

int mod2dense_flip
( mod2dense *m, 	/* Matrix to flip element in */
  int row,		/* Row of element (starting with zero) */
  int col		/* Column of element (starting with zero) */
)
{
  mod2word *w;
  int b;

  if (row<0 || row>=mod2dense_rows(m) || col<0 || col>=mod2dense_cols(m))
  { fprintf(stderr,"mod2dense_flip: row or column index out of bounds\n");
    exit(1);
  }

  b = 1 ^ mod2_getbit (m->col[col][row>>mod2_wordsize_shift],
                       row&mod2_wordsize_mask);

  w = &m->col[col][row>>mod2_wordsize_shift];

  *w = b ? mod2_setbit1(*w,row&mod2_wordsize_mask)
         : mod2_setbit0(*w,row&mod2_wordsize_mask);

  return b;
}


/* COMPUTE THE TRANSPOSE OF A DENSE MOD2 MATRIX. */

void mod2dense_transpose
( mod2dense *m,		/* Matrix to compute transpose of (left unchanged) */
  mod2dense *r		/* Result of transpose operation */
)
{
  mod2word w, v, *p;
  int k1, j1, i2, j2;

  if (mod2dense_rows(m)!=mod2dense_cols(r)
   || mod2dense_cols(m)!=mod2dense_rows(r))
  { fprintf(stderr,
     "mod2dense_transpose: Matrices have incompatible dimensions\n");
    exit(1);
  }

  if (r==m)
  { fprintf(stderr,
     "mod2dense_transpose: Result matrix is the same as the operand\n");
    exit(1);
  }

  mod2dense_clear(r);

  for (j1 = 0; j1<mod2dense_cols(m); j1++)
  {
    i2 = j1 >> mod2_wordsize_shift;
    v = 1 << (j1 & mod2_wordsize_mask);

    p = m->col[j1];
    k1 = 0;

    for (j2 = 0; j2<mod2dense_cols(r); j2++)
    { if (k1==0)
      { w = *p++;
        k1 = mod2_wordsize;
      }
      if (w&1)
      { r->col[j2][i2] |= v;
      }
      w >>= 1;
      k1 -= 1;
    }
  }
}


/* ADD TWO DENSE MOD2 MATRICES. */

void mod2dense_add
( mod2dense *m1,	/* Left operand of add */
  mod2dense *m2,	/* Right operand of add */
  mod2dense *r		/* Place to store result of add */
)
{
  int j, k;

  if (mod2dense_rows(m1)!=mod2dense_rows(r)
   || mod2dense_cols(m1)!=mod2dense_cols(r)
   || mod2dense_rows(m2)!=mod2dense_rows(r)
   || mod2dense_cols(m2)!=mod2dense_cols(r))
  { fprintf(stderr,"mod2dense_add: Matrices have different dimensions\n");
    exit(1);
  }

  for (j = 0; j<mod2dense_cols(r); j++)
  { for (k = 0; k<r->n_words; k++)
    { r->col[j][k] = m1->col[j][k] ^ m2->col[j][k];
    }
  }
}


/* MULTIPLY TWO DENSE MOD2 MATRICES.

   The algorithm used runs faster if the second matrix (right operand of the
   multiply) is sparse, but it is also appropriate for dense matrices.  This
   procedure could be speeded up a bit by replacing the call of mod2dense_get
   with in-line code that avoids division, but this doesn't seem worthwhile
   at the moment.
*/

void mod2dense_multiply
( mod2dense *m1, 	/* Left operand of multiply */
  mod2dense *m2,	/* Right operand of multiply */
  mod2dense *r		/* Place to store result of multiply */
)
{
  int i, j, k;

  if (mod2dense_cols(m1)!=mod2dense_rows(m2)
   || mod2dense_rows(m1)!=mod2dense_rows(r)
   || mod2dense_cols(m2)!=mod2dense_cols(r))
  { fprintf(stderr,
     "mod2dense_multiply: Matrices have incompatible dimensions\n");
    exit(1);
  }

  if (r==m1 || r==m2)
  { fprintf(stderr,
      "mod2dense_multiply: Result matrix is the same as one of the operands\n");
    exit(1);
  }

  mod2dense_clear(r);

  for (j = 0; j<mod2dense_cols(r); j++)
  { for (i = 0; i<mod2dense_rows(m2); i++)
    { if (mod2dense_get(m2,i,j))
      { for (k = 0; k<r->n_words; k++)
        { r->col[j][k] ^= m1->col[i][k];
        }
      }
    }
  }
}


/* SEE WHETHER TWO DENSE MOD2 MATRICES ARE EQUAL. */

int mod2dense_equal
( mod2dense *m1,
  mod2dense *m2
)
{
  int k, j, w;
  mod2word m;

  if (mod2dense_rows(m1)!=mod2dense_rows(m2)
   || mod2dense_cols(m1)!=mod2dense_cols(m2))
  { fprintf(stderr,"mod2dense_equal: Matrices have different dimensions\n");
    exit(1);
  }

  w = m1->n_words;

  /* Form a mask that has 1s in the lower bit positions corresponding to
     bits that contain information in the last word of a matrix column. */

  m = (1 << (mod2_wordsize - (w*mod2_wordsize-m1->n_rows))) - 1;

  for (j = 0; j<mod2dense_cols(m1); j++)
  {
    for (k = 0; k<w-1; k++)
    { if (m1->col[j][k] != m2->col[j][k]) return 0;
    }

    if ((m1->col[j][k]&m) != (m2->col[j][k]&m)) return 0;
  }

  return 1;
}


/* INVERT A DENSE MOD2 MATRIX. */

int mod2dense_invert
( mod2dense *m,		/* The matrix to find the inverse of (destroyed) */
  mod2dense *r		/* Place to store the inverse */
)
{
  mod2word *s, *t;
  int i, j, k, n, w, k0, b0;

  if (mod2dense_rows(m)!=mod2dense_cols(m))
  { fprintf(stderr,"mod2dense_invert: Matrix to invert is not square\n");
    exit(1);
  }

  if (r==m)
  { fprintf(stderr,
      "mod2dense_invert: Result matrix is the same as the operand\n");
    exit(1);
  }

  n = mod2dense_rows(m);
  w = m->n_words;

  if (mod2dense_rows(r)!=n || mod2dense_cols(r)!=n)
  { fprintf(stderr,
     "mod2dense_invert: Matrix to receive inverse has wrong dimensions\n");
    exit(1);
  }

  mod2dense_clear(r);
  for (i = 0; i<n; i++)
  { mod2dense_set(r,i,i,1);
  }

  for (i = 0; i<n; i++)
  {
    k0 = i >> mod2_wordsize_shift;
    b0 = i & mod2_wordsize_mask;

    for (j = i; j<n; j++)
    { if (mod2_getbit(m->col[j][k0],b0)) break;
    }

    if (j==n) return 0;

    if (j!=i)
    {
      t = m->col[i];
      m->col[i] = m->col[j];
      m->col[j] = t;

      t = r->col[i];
      r->col[i] = r->col[j];
      r->col[j] = t;
    }

    for (j = 0; j<n; j++)
    { if (j!=i && mod2_getbit(m->col[j][k0],b0))
      { s = m->col[j];
        t = m->col[i];
        for (k = k0; k<w; k++) s[k] ^= t[k];
        s = r->col[j];
        t = r->col[i];
        for (k = 0; k<w; k++) s[k] ^= t[k];
      }
    }
  }

  return 1;
}


/* INVERT A DENSE MOD2 MATRIX WITH ROWS & COLUMNS SELECTED FROM BIGGER MATRIX.*/

int mod2dense_invert_selected
( mod2dense *m,		/* Matrix from which to pick a submatrix to invert */
  mod2dense *r,		/* Place to store the inverse */
  int *rows,		/* Set to indexes of rows used and not used */
  int *cols		/* Set to indexes of columns used and not used */
)
{
  mod2word *s, *t;
  int i, j, k, n, n2, w, k0, b0, c, R;

  if (r==m)
  { fprintf(stderr,
      "mod2dense_invert_selected2: Result matrix is the same as the operand\n");
    exit(1);
  }

  n = mod2dense_rows(m);
  w = m->n_words;

  n2 = mod2dense_cols(m);

  if (mod2dense_rows(r)!=n || mod2dense_cols(r)!=n2)
  { fprintf(stderr,
"mod2dense_invert_selected2: Matrix to receive inverse has wrong dimensions\n");
    exit(1);
  }

  mod2dense_clear(r);

  for (i = 0; i<n; i++)
  { rows[i] = i;
  }

  for (j = 0; j<n2; j++)
  { cols[j] = j;
  }

  R = 0;
  i = 0;

  for (;;)
  {
    while (i<n-R)
    {
      k0 = rows[i] >> mod2_wordsize_shift;
      b0 = rows[i] & mod2_wordsize_mask;

      for (j = i; j<n2; j++)
      { if (mod2_getbit(m->col[cols[j]][k0],b0)) break;
      }

      if (j<n2) break;

      R += 1;
      c = rows[i];
      rows[i] = rows[n-R];
      rows[n-R] = c;

    }

    if (i==n-R) break;

    c = cols[j];
    cols[j] = cols[i];
    cols[i] = c;

    mod2dense_set(r,rows[i],c,1);

    for (j = 0; j<n2; j++)
    { if (j!=c && mod2_getbit(m->col[j][k0],b0))
      { s = m->col[j];
        t = m->col[c];
        for (k = 0; k<w; k++) s[k] ^= t[k];
        s = r->col[j];
        t = r->col[c];
        for (k = 0; k<w; k++) s[k] ^= t[k];
      }
    }

    i += 1;
  }

  for (j = n-R; j<n; j++)
  { s = r->col[cols[j]];
    for (k = 0; k<w; k++) s[k] = 0;
  }

  return R;
}


/* FORCIBLY INVERT A DENSE MOD2 MATRIX. */

int mod2dense_forcibly_invert
( mod2dense *m, 	/* The matrix to find the inverse of (destroyed) */
  mod2dense *r,		/* Place to store the inverse */
  int *a_row,		/* Place to store row indexes of altered elements */
  int *a_col		/* Place to store column indexes of altered elements */
)
{
  mod2word *s, *t;
  int i, j, k, n, w, k0, b0;
  int u, c;

  if (mod2dense_rows(m)!=mod2dense_cols(m))
  { fprintf(stderr,
      "mod2dense_forcibly_invert: Matrix to invert is not square\n");
    exit(1);
  }

  if (r==m)
  { fprintf(stderr,
      "mod2dense_forcibly_invert: Result matrix is the same as the operand\n");
    exit(1);
  }

  n = mod2dense_rows(m);
  w = m->n_words;

  if (mod2dense_rows(r)!=n || mod2dense_cols(r)!=n)
  { fprintf(stderr,
 "mod2dense_forcibly_invert: Matrix to receive inverse has wrong dimensions\n");
    exit(1);
  }

  mod2dense_clear(r);
  for (i = 0; i<n; i++)
  { mod2dense_set(r,i,i,1);
  }

  for (i = 0; i<n; i++)
  { a_row[i] = -1;
    a_col[i] = i;
  }

  for (i = 0; i<n; i++)
  {
    k0 = i >> mod2_wordsize_shift;
    b0 = i & mod2_wordsize_mask;

    for (j = i; j<n; j++)
    { if (mod2_getbit(m->col[j][k0],b0)) break;
    }

    if (j==n)
    { j = i;
      mod2dense_set(m,i,j,1);
      a_row[i] = i;
    }

    if (j!=i)
    {
      t = m->col[i];
      m->col[i] = m->col[j];
      m->col[j] = t;

      t = r->col[i];
      r->col[i] = r->col[j];
      r->col[j] = t;

      u = a_col[i];
      a_col[i] = a_col[j];
      a_col[j] = u;
    }

    for (j = 0; j<n; j++)
    { if (j!=i && mod2_getbit(m->col[j][k0],b0))
      { s = m->col[j];
        t = m->col[i];
        for (k = k0; k<w; k++) s[k] ^= t[k];
        s = r->col[j];
        t = r->col[i];
        for (k = 0; k<w; k++) s[k] ^= t[k];
      }
    }
  }

  c = 0;
  for (i = 0; i<n; i++)
  { if (a_row[i]!=-1)
    { a_row[c] = a_row[i];
      a_col[c] = a_col[i];
      c += 1;
    }
  }

  return c;
}

void sparse_encode
( char *sblk,
  char *cblk
)
{
  int i, j;

  mod2entry *e;
  char *x, *y;

  x = (char *)chk_alloc (M, sizeof *x);
  y = (char *)chk_alloc (M, sizeof *y);

  /* Multiply the vector of source bits by the systematic columns of the
     parity check matrix, giving x.  Also copy these bits to the coded block. */

  for (i = 0; i<M; i++) x[i] = 0;

  for (j = M; j<N; j++)
  {
    cblk[cols[j]] = sblk[j-M];

    if (sblk[j-M]==1)
    { for (e = mod2sparse_first_in_col(H,cols[j]);
           !mod2sparse_at_end(e);
           e = mod2sparse_next_in_col(e))
      { x[mod2sparse_row(e)] ^= 1;
      }
    }
  }

  /* Solve Ly=x for y by forward substitution, then U(cblk)=y by backward
     substitution. */

  if (!mod2sparse_forward_sub(L,rows,x,y)
   || !mod2sparse_backward_sub(U,cols,y,cblk))
  {
    abort(); /* Shouldn't occur, even if the parity check matrix has
                redundant rows */
  }

  free(x);
  free(y);
}


/* ENCODE A BLOCK USING DENSE REPRESENTATION OF GENERATOR MATRIX. */

void dense_encode
( char *sblk,
  char *cblk,
  mod2dense *u,
  mod2dense *v
)
{
  int j;

  /* Copy source bits to the systematic part of the coded block. */

  for (j = M; j<N; j++)
  { cblk[cols[j]] = sblk[j-M];
  }

  /* Multiply by Inv(A) X B to produce check bits. */

  for (j = M; j<N; j++)
  { mod2dense_set(u,j-M,0,sblk[j-M]);
  }

  mod2dense_multiply(G,u,v);

  /* Copy check bits to the right places in the coded block. */

  for (j = 0; j<M; j++)
  { cblk[cols[j]] = mod2dense_get(v,j,0);
  }
}


/* ENCODE A BLOCK USING MIXED REPRESENTATION OF GENERATOR MATRIX. */

void mixed_encode
( char *sblk,
  char *cblk,
  mod2dense *u,
  mod2dense *v
)
{
  mod2entry *e;
  int j;

  /* Multiply the vector of source bits by the message bit columns of the
     parity check matrix.  Also copy these bits to the coded block.  Take
     account of how columns have been reordered. */

  mod2dense_clear(u);

  for (j = M; j<N; j++)
  {
    cblk[cols[j]] = sblk[j-M];

    if (sblk[j-M]==1)
    { for (e = mod2sparse_first_in_col(H,cols[j]);
           !mod2sparse_at_end(e);
           e = mod2sparse_next_in_col(e))
      { (void) mod2dense_flip(u,mod2sparse_row(e),0);
      }
    }
  }

  /* Multiply by Inv(A) to produce check bits. */

  mod2dense_multiply(G,u,v);

  /* Copy check bits to the right places in the coded block. */

  for (j = 0; j<M; j++)
  { cblk[cols[j]] = mod2dense_get(v,j,0);
  }
}
int intio_read
( FILE *f   /* File to read from */
)
{

  unsigned char b[4];
  int top;
  int i;

  for (i = 0; i<4; i++)
  { if (fread(&b[i],1,1,f) != 1) return 0;
  }

  top = b[3]>127 ? (int)b[3] - 256 : b[3];
  return (top<<24) + (b[2]<<16) + (b[1]<<8) + b[0];
}


/* WRITE AN INTEGER ONE BYTE AT A TIME.  Four bytes are written, ordered from
   low to high order.  These are considered to represent a signed integer,
   in two's complement form.  This should work as long as the integer passed
   can be represented in four bytes, even if a C "int" is longer than this.

   The file written to should have been opened as "binary".
*/

void intio_write
( FILE *f,  /* File to write to */
  int v     /* Value to write to file */
)
{
  unsigned char b;
  int i;

  for (i = 0; i<3; i++)
  { b = v&0xff;
    fwrite(&b,1,1,f);
    v >>= 8;
  }

  b = v>0 ? v : v+256;
  fwrite(&b,1,1,f);
}
void *chk_alloc
( unsigned n,		/* Number of elements */
  unsigned size		/* Size of each element */
)
{
  void *p;

  p = calloc(n,size);

  if (p==0)
  { fprintf(stderr,"Ran out of memory (while trying to allocate %d bytes)\n",
      n*size);
    exit(1);
  }

  return p;
}
FILE *open_file_std
( char *fname,	/* Name of file to open, or "-" for stdin/stdout */
  char *mode	/* Mode for opening: eg, "r" or "w" */
)
{
  if (strcmp(fname,"-")==0)
  { switch (mode[0])
    { case 'r':
      { return stdin;
      }
      case 'w':
      { return stdout;
      }
      default:
      { fprintf(stderr,"Bad mode passed to open_file_std: %s\n",mode);
        exit(1);
      }
    }
  }
  else
  {
    return fopen(fname,mode);
  }
}
int blockio_flush = 0;	/* Should blocks written be immediately flushed? */

int blockio_read
( FILE *f,    /* File to read from */
  char *b,    /* Place to store bits read */
  int l       /* Length of block */
)
{
  int i, c;
  for (i = 0; i<l; i++)
  {

    do
    {
    	c = getc(f);
      if (c==EOF)
      { if (i>0)
        { fprintf(stderr,
           "Warning: Short block (%d long) at end of input file ignored\n",i);
        }
        return EOF;
      }
    } while (c==' ' || c=='\t' || c=='\n' || c=='\r');

    if (c!='0' && c!='1')
    { fprintf(stderr,"Bad character in binary file (not '0' or '1')\n");
      exit(1);
    }

    b[i] = c=='1';
  }

  return 0;
}


/* WRITE A BLOCK OF BITS.  Bits are written as '0' and '1' characters, with
   no spaces between them, followed by a newline. */

void blockio_write( FILE *f,     /* File to write to */  char *b,     /* Block of bits to write */  int l        /* Length of block */)
{
  unsigned int i;
  for (i = 0; i<l; i++)
  {
	if (b[i]!=0 && b[i]!=1) abort();
    	putc("01"[b[i]],f);
  }
}
}
static void initialize (void)
{
  int i, j, k, w;
  char b;
  FILE *f;

  if (!initialized)
  {
    f = fopen("./randfile","rb");


    for (i = 0; i<N_tables; i++)
    { for (j = 0; j<Table_size; j++)
      { w = 0;
        for (k = 0; k<4; k++)
        { if (fread(&b,1,1,f)!=1)
          {
            exit(1);
          }
          w = (w<<8) | (b&0xff);
        }
        rn[i][j] = w;
      }
    }

    state = &state0;

    initialized = 1;

    rand_seed(1);
  }
}


/* SET CURRENT STATE ACCORDING TO SEED. */

void rand_seed
( int seed
)
{
  int j;

  if (!initialized) initialize();

  state->seed = seed;

  state->state48[0] = seed>>16;
  state->state48[1] = seed&0xffff;
  state->state48[2] = rn[0][(seed&0x7fffffff)%Table_size];

  for (j = 0; j<N_tables; j++)
  { state->ptr[j] = seed%Table_size;
    seed /= Table_size;
  }
}
double rand_uniform (void)
{
  return (double)rand_word() / (1.0+(double)0x7fffffff);
}
void rand_use_state
( rand_state *st
)
{
  if (!initialized) initialize();

  state = st;
}


/* RETURN POINTER TO CURRENT STATE. */
double rand_uniopen (void)
{
  return (0.5+(double)rand_word()) / (1.0+(double)0x7fffffff);
}


/* GENERATE RANDOM INTEGER FROM 0, 1, ..., (n-1). */

int rand_int
( int n
)
{
  return (int) (n * rand_uniform());
}
rand_state *rand_get_state (void)
{
  if (!initialized) initialize();

  return state;
}
int rand_word(void)
{
  int v;
  int j;

  if (!initialized) initialize();

  v = this_nrand48(state->state48);

  for (j = 0; j<N_tables; j++)
  { v ^= rn[j][state->ptr[j]];
  }

  for (j = 0; j<N_tables && state->ptr[j]==Table_size-1; j++)
  { state->ptr[j] = 0;
  }

  if (j<N_tables)
  { state->ptr[j] += 1;
  }

  return v & 0x7fffffff;
}

int rand_pickd
( double *p,
  int n
)
{
  double t, r;
  int i;

  t = 0;
  for (i = 0; i<n; i++)
  { if (p[i]<0) abort();
    t += p[i];
  }

  if (t<=0) abort();

  r = t * rand_uniform();

  for (i = 0; i<n; i++)
  { r -= p[i];
    if (r<0) return i;
  }

  /* Return value with non-zero probability if we get here due to roundoff. */

  for (i = 0; i<n; i++)
  { if (p[i]>0) return i;
  }

  abort();
}


/* SAME PROCEDURE AS ABOVE, BUT WITH FLOAT ARGUMENT. */

int rand_pickf
( float *p,
  int n
)
{
  double t, r;
  int i;

  t = 0;
  for (i = 0; i<n; i++)
  { if (p[i]<=0) abort();
    t += p[i];
  }

  if (t<=0) abort();

  r = t * rand_uniform();

  for (i = 0; i<n; i++)
  { r -= p[i];
    if (r<0) return i;
  }

  /* Return value with non-zero probability if we get here due to roundoff. */

  for (i = 0; i<n; i++)
  { if (p[i]>0) return i;
  }

  abort();
}


/* GENERATE RANDOM PERMUTATION OF INTEGERS FROM 1 TO N. */

void rand_permutation
( int *perm,		/* Place to store permutation */
  int n			/* Number of integers to permute */
)
{
  int i, j, t;

  for (i = 0; i<n; i++)
  { perm[i] = i+1;
  }

  for (i = 0; i<n; i++)
  { t = perm[i];
    j = i + rand_int(n-i);
    perm[i] = perm[j];
    perm[j] = t;
  }
}


/* POISSON GENERATOR.  The method used is simple, but not very fast.  See
   Devroye, p. 503.  Very large means are done using Gaussian approximation. */

int rand_poisson
( double lambda
)
{ int v;
  if (lambda>10000)
  { v = (int) (lambda + rand_gaussian()*sqrt(lambda) + 0.5);
  }
  else
  { v = 0;
    for (;;)
    { lambda -= rand_exp();
      if (lambda<=0) break;
      v += 1;
    }
  }
  return v;
}


/* GAUSSIAN GENERATOR.  Done by using the Box-Muller method, but only one
   of the variates is retained (using both would require saving more state).
   See Devroye, p. 235.

   As written, should never deliver exactly zero, which may sometimes be
   helpful. */

double rand_gaussian (void)
{
  double a, b;

  a = rand_uniform();
  b = rand_uniopen();

  return cos(2.0*M_PI*a) * sqrt(-2.0*log(b));
}


/* EXPONENTIAL GENERATOR.  See Devroye, p. 29.  Written so as to never
   return exactly zero. */

double rand_exp (void)
{
  return -log(rand_uniopen());
}


/* LOGISTIC GENERATOR.  Just inverts the CDF. */

double rand_logistic (void)
{ double u;
  u = rand_uniopen();
  return log(u/(1-u));
}


/* CAUCHY GENERATOR.  See Devroye, p. 29. */

double rand_cauchy (void)
{
  return tan (M_PI * (rand_uniopen()-0.5));
}


/* GAMMA GENERATOR.  Generates a positive real number, r, with density
   proportional to r^(a-1) * exp(-r).  See Devroye, p. 410 and p. 420.
   Things are fiddled to avoid ever returning a value that is very near
   zero. */

double rand_gamma
( double a
)
{
  double b, c, X, Y, Z, U, V, W;

  if (a<0.00001)
  { X = a;
  }

  else if (a<=1)
  {
    U = rand_uniopen();
    X = rand_gamma(1+a) * pow(U,1/a);
  }

  else if (a<1.00001)
  { X = rand_exp();
  }

  else
  {
    b = a-1;
    c = 3*a - 0.75;

    for (;;)
    {
      U = rand_uniopen();
      V = rand_uniopen();

      W = U*(1-U);
      Y = sqrt(c/W) * (U-0.5);
      X = b+Y;

      if (X>=0)
      {
        Z = 64*W*W*W*V*V;

        if (Z <= 1 - 2*Y*Y/X || log(Z) <= 2 * (b*log(X/b) - Y)) break;
      }
    }
  }

  return X<1e-30 && X<a ? (a<1e-30 ? a : 1e-30) : X;
}


/* BETA GENERATOR. Generates a real number, r, in (0,1), with density
   proportional to r^(a-1) * (1-r)^(b-1).  Things are fiddled to avoid
   the end-points, and to make the procedure symmetric between a and b. */

double rand_beta
( double a,
  double b
)
{
  double x, y, r;

  do
  { x = rand_gamma(a);
    y = rand_gamma(b);
    r = 1.0 + x/(x+y);
    r = r - 1.0;
  } while (r<=0.0 || r>=1.0);

  return r;
}


/* ROUTINES FROM THE GNU C LIBRARY.  These were modified to extract
   only the routines used here, and to allow them to be included in
   this module without any possible name conflict with other modules.
   Inclusion here ensures that these routines are always available, and
   operate in exactly the same way on all systems.  The routines as copied
   below are still easily useable by other programs by simply inserting
   this source code into an appropriate source file.

   The following is the copyright notice for these routines:

     Copyright (C) 1995, 1996, 1997, 2002 Free Software Foundation, Inc.
     This file is part of the GNU C Library.
     Contributed by Ulrich Drepper <drepper@gnu.ai.mit.edu>, August 1995.

     The GNU C Library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public
     License as published by the Free Software Foundation; either
     version 2.1 of the License, or (at your option) any later version.

     The GNU C Library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Lesser General Public License for more details.

     You should have received a copy of the GNU Lesser General Public
     License along with the GNU C Library; if not, write to the Free
     Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
     02111-1307 USA.

   The GNU Lesser General Public License is included with these source
   files in the file LGPL. */

struct this_drand48_data
  {
    unsigned short int x[3];	/* Current state.  */
    unsigned short int old_x[3]; /* Old state.  */
    unsigned short int c;	/* Additive const. in congruential formula.  */
    unsigned short int init;	/* Flag for initializing.  */
    unsigned long long int a;	/* Factor in congruential formula.  */
  };
  struct this_drand48_data libc_this_drand48_data;


  static int this_drand48_iterate (unsigned short int xsubi[3], struct this_drand48_data *buffer)
  {
    uint64_t X;
    uint64_t result;

    /* Initialize buffer, if not yet done.  */
    if (!buffer->init)
      {
        buffer->a = 0x5deece66dull;
        buffer->c = 0xb;
        buffer->init = 1;
      }

    /* Do the real work.  We choose a data type which contains at least
       48 bits.  Because we compute the modulus it does not care how
       many bits really are computed.  */

    X = (uint64_t) xsubi[2] << 32 | (uint32_t) xsubi[1] << 16 | xsubi[0];

    result = X * buffer->a + buffer->c;

    xsubi[0] = result & 0xffff;
    xsubi[1] = (result >> 16) & 0xffff;
    xsubi[2] = (result >> 32) & 0xffff;

    return 0;
  }


  static int this_nrand48_r (unsigned short int xsubi[3],struct this_drand48_data *buffer, long int *result)
  {
    /* Compute next state.  */
    if (this_drand48_iterate (xsubi, buffer) < 0)
      return -1;

    /* Store the result.  */
    if (sizeof (unsigned short int) == 2)
      *result = xsubi[2] << 15 | xsubi[1] >> 1;
    else
      *result = xsubi[2] >> 1;

    return 0;
  }





static long int this_nrand48 (unsigned short int xsubi[3])
  {
    long int result;

    (void) this_nrand48_r (xsubi, &libc_this_drand48_data, &result);

    return result;
  }
