/******************************************************************************
*     FILE: pgapack.h: This file contains all constant and structure
*                      definitions definitions for PGAPack as welll as all
*                      function declarations.
*     Authors: David M. Levine, Philip L. Hallstrom, David M. Noelle
******************************************************************************/
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <string.h> 
#include <ctype.h>

#ifdef OPTIMIZE                      /* remove debug prints for */
#define PGADebugPrint(a,b,c,x,y,z)   /* optimized version       */
#endif

/*****************************************
*           BINARY   MACROS              *
*****************************************/
#if defined(NeXT) || defined(sun4) || defined(AIX) || defined(IRIX)\
|| defined(freebsd) || defined(hp) || defined(sun)
#define WL 32
#else
#error Your architecture has an unknown word length
#endif
#define BIT(x,y)    (y&(1<<((WL-1)-(x))))       /* true if bit is 1,         */
#define SET(x,y)    (y|=(1<<((WL-1)-(x))))      /* set a bit to 1            */
#define UNSET(x,y)  (y&=(~(1<<((WL-1)-(x)))))   /* set a bit to 0, clear     */
#define TOGGLE(x,y) (y^=(1<<((WL-1)-(x))))      /* complement a bits value   */
#define INDEX(ix,bx,bit,WL) ix=bit/WL;bx=bit%WL /* map global column (bit)   */
                                                /* to word (ix) and bit (bx) */

/*****************************************
*       ABSTRACT DATA TYPES              *
*****************************************/
#define PGA_DATATYPE_BINARY      1    /* Array of unsigned ints          */
                                      /* parsed into bits    : binary.c  */
#define PGA_DATATYPE_INTEGER     2    /* Array of ints       : integer.c */
#define PGA_DATATYPE_REAL        3    /* Array of doubles    : real.c    */
#define PGA_DATATYPE_USER        4    /*  --user defined--               */
    
#define PGABinary                unsigned int
#define PGAInteger               long int
#define PGAReal                  double

#define PGA_INT                   1
#define PGA_DOUBLE                2
#define PGA_CHAR                  3
#define PGA_VOID                  4

    
/*****************************************
*       BOOLEANS &  FLAGS                *
*****************************************/
#define PGATRUE                   1
#define PGAFALSE                  0
#define PGAERROR                 -1

#define PGA_FATAL                 1
#define PGA_WARNING               2

#define PGA_UNINITIALIZED_INT    -3827
#define PGA_UNINITIALIZED_DOUBLE -968.3827

/*****************************************
*    TEMP & POP REFERENT CONSTANTS       *
*****************************************/
#define PGA_TEMP1                -1138
#define PGA_TEMP2                -4239

#define PGA_OLDPOP               -6728
#define PGA_NEWPOP               -8376


/*****************************************
*        DEBUG LEVELS                    *
*****************************************/
#define PGADEBUG_ENTERED        11
#define PGADEBUG_EXIT           12
#define PGADEBUG_MALLOC         14
#define PGADEBUG_PRINTVAR       20
#define PGADEBUG_ENTEREDHIGH    31
#define PGADEBUG_EXITHIGH       32
#define PGADEBUG_NUMFLAGS	1000
#define PGADEBUG_MAXPROCS	300


/*****************************************
*           DIRECTION                    *
*****************************************/
#define PGA_MAXIMIZE            1    /* specify direction for fitness calc  */
#define PGA_MINIMIZE            2    /* specify direction for fitness calc  */
    
/*****************************************
*         STOPPING CRITERIA              *
*****************************************/
#define PGA_STOP_MAXITER        1    /* Stop: for maximum iterations      */
#define PGA_STOP_NOCHANGE       2    /* Stop: no change in best string    */
#define PGA_STOP_TOOSIMILAR     4    /* Stop: homogeneous population      */

/*****************************************
*            CROSSOVER                   *
*****************************************/
#define PGA_CROSSOVER_ONEPT     1    /* One point crossover                */
#define PGA_CROSSOVER_TWOPT     2    /* Two point crossover                */
#define PGA_CROSSOVER_UNIFORM   3    /* Uniform   crossover                */

/*****************************************
*            SELECTION                   *
*****************************************/
#define PGA_SELECT_PROPORTIONAL 1    /* proportional selection              */
#define PGA_SELECT_SUS          2    /* stochastic universal selection      */
#define PGA_SELECT_TOURNAMENT   3    /* tournament selection                */
#define PGA_SELECT_PTOURNAMENT  4    /* probabilistic tournament selection  */

/*****************************************
*            FITNESS                     *
*****************************************/
#define PGA_FITNESS_RAW         1    /* use raw fitness (evaluation)        */
#define PGA_FITNESS_NORMAL      2    /* linear normalization fitness        */
#define PGA_FITNESS_RANKING     3    /* linear ranking fitness              */

/*****************************************
*            FITNESS (MINIMIZATION)      *
*****************************************/
#define PGA_FITNESSMIN_RECIPROCAL  1 /* reciprocal fitness                  */
#define PGA_FITNESSMIN_CMAX        2 /* cmax fitness                        */

/*****************************************
*               MUTATION                 *
*****************************************/
#define PGA_MUTATION_FIXED      1    /* User specified mutation rate        */
#define PGA_MUTATION_STRINGLEN  2    /* 1/L is mutation rate                */
#define PGA_MUTATION_HAMMING    3    /* mutation rate = f(Hamming distance) */
#define PGA_MUTATION_UNIFORM    4    /* +- uniform random no. (Real only)  */
#define PGA_MUTATION_GAUSSIAN   5    /* +- Gaussian random no. (Real only) */
    
/*****************************************
*        POPULATION REPLACEMENT          *
*****************************************/
#define PGA_POPREPL_BEST         1   /* Select best   string                */
#define PGA_POPREPL_RANDOM_NOREP 2   /* Select random string w/o replacement*/
#define PGA_POPREPL_RANDOM_REP   3   /* Select random string w/  replacement*/

/****************************************
 *       REPORT OPTIONS                 *
 ****************************************/
#define PGA_ONLINE              1    /* Print the online analysis           */
#define PGA_OFFLINE             2    /* Print the offline analysis          */
#define PGA_HAMMING             4    /* Print the Hamming distance          */
#define PGA_STRING              8    /* Print the string                    */
#define PGA_WORST               16   /* Print the worst individual          */

/*****************************************
*            RANDOMIZER                  *
*****************************************/
#define PGA_PERMUTE             1  /* unqique integer values in 1, PopSize  */
#define PGA_RANGE               2  /* pick (nonunique) with prob=1/rangelen */

/*****************************************
*         SET USER FUNCTION              *
*****************************************/
#define PGA_USERFUNCTION_CREATESTRING 0
#define PGA_USERFUNCTION_MUTATION     1
#define PGA_USERFUNCTION_CROSSOVER    2
#define PGA_USERFUNCTION_WRITESTRING  3
#define PGA_USERFUNCTION_COPYSTRING   4
#define PGA_USERFUNCTION_DUPLICATE    5
#define PGA_USERFUNCTION_RANDOMIZE    6
      
/*****************************************
*       INDIVIDUAL STRUTURE              *
*****************************************/

typedef struct {                    /* primary population data structure   */
  double evalfunc;                  /* evaluation function value           */
  double fitness;                   /* fitness    function value           */
  int    evaluptodate;              /* flag whether evalfunc is current    */
  void   *chrom;                    /* pointer to the GA string            */
} PGAIndividual;


/*****************************************
*          GA ALGORITHM STRUCTURE        *
*****************************************/
typedef struct {
    int datatype;            /* data type: binary, integer, or real       */
    int optdir;              /* direction of optimization                 */
    int tw;                  /* total number of words, full + partial     */
    int fw;                  /* number of full (WL length) words          */
    int eb;                  /* number of extra bits in last NOT full word*/
    int PopSize;             /* Number of strings to use                  */
    int StringLen;           /* string lengths                            */
    int StoppingRule;        /* Termination Criteria                      */
    int MaxIter;             /* Maximum number of iterations to run       */
    int MaxNoChange;         /* # of iters with no change before stopping */
    int MaxSimilarity;       /* % of pop the same before stopping         */
    int NumReplace;          /* Number of string to replace each gen      */
    int CrossoverType;       /* Type of crossover for genetic algorithm   */
    int SelectType;          /* Type of selection for genetic algorithm   */
    int FitnessType;         /* Type of fitness transformation used       */
    int FitnessMinType;      /* Transformation for minimization problems  */
    int MutationType;        /* Type of mutation used                     */
    int MutateOnlyNoCross;   /* Mutate only strings not from crossover    */
    double MutateRealVal;    /* Multiplier to mutate Real strings with    */
    int MutateIntegerVal;    /* Multiplier to mutate Integer strings with */
    int NoDuplicates;        /* Don't allow duplicate strings             */
    double MutationProb;     /* Starting mutation probability             */
    double CrossoverProb;    /* Crossover probability                     */
    double UniformCrossProb; /* Prob of bit select in uniform crossover   */
    double PTournamentProb;  /* Prob of selection in Prob. Tournament     */
    double FitnessRankMax;   /* MAX value for use in ranking              */
    int PopReplace;          /* Method of choosing ind.s to copy to newpop*/
    PGAIndividual *oldpop;   /* pointer to population (old)               */
    PGAIndividual *newpop;   /* pointer to population (new)               */
    int iter;                /* iteration (generation) counter            */
    int ItersOfSame;         /* # iters with no change in best            */
    int PercentSame;         /* % of pop that is homogeneous              */
    int *selected;           /* array of indices for selection            */
    int *sorted;             /* array of sorted individual indices        */
    int SelectIndex;         /* index of Select for next two individuals  */
} PGAAlgorithm;


struct str_PGAContext;

/*****************************************
*        OPERATIONS STRUCTURE            *
*****************************************/
typedef struct {
    void (*CreateString)(struct str_PGAContext *, int, int, int);
    int  (*Mutation)(struct str_PGAContext *, int, int, double);
    void (*Crossover)(struct str_PGAContext *, int, int, int, int, int, int);
    void (*WriteString)(struct str_PGAContext *, FILE *, int, int);
    void (*CopyString)(struct str_PGAContext *, int, int, int, int);
    int  (*Duplicate)(struct str_PGAContext *, int, int, int, int);
    void (*Randomize)(struct str_PGAContext *, int, int);
} PGAOperations;


/*****************************************
*          REPORT STRUCTURE              *
*****************************************/
typedef struct {
  int    PrintFreq;                 /* How often to print statistics reports*/
  int    PrintOptions;
  double Offline;
  double Online;
  double Best;
} PGAReport;


/*****************************************
*          SYSTEM STRUCTURE              *
*****************************************/
typedef struct {
    int    fortran;                 /* user routines in Fortran or C?        */
    int    SetUpCalled;             /* has PGASetUp been called?             */
    int    PGAMaxInt;               /* largest  int     of machine           */
    int    PGAMinInt;               /* smallest int     of machine           */
    double PGAMaxDouble;            /* largest  double  of machine           */
    double PGAMinDouble;            /* smallest double  of machine           */
} PGASystem;


/*****************************************
*          DEBUG STRUCTURE               *
*****************************************/
typedef struct {
    int PGADebugFlags[PGADEBUG_NUMFLAGS + 1];
} PGADebug;

/*****************************************
*      INITIALIZATION STRUCTURE          *
*****************************************/
typedef struct {
    int    RandomInit;             /* flag whether to randomize strings    */
    double BinaryProbability;      /* probability that a Bit will be 1     */
    int    IntegerType;            /* type of integer randomization        */
    int    *IntegerMin;            /* minimum of range of integers         */
    int    *IntegerMax;            /* maximum of range of integers         */
    double *RealMin;               /* minimum of range of reals            */
    double *RealMax;               /* maximum of range of reals            */
    int  RandomSeed;             /* integer to seed random numbers with  */
} PGAInitialize;

/*****************************************
*      SCRATCH DATA STRUCTURES           *
*****************************************/
typedef struct {
    int *intscratch;                /* integer-scratch space                 */
    double *dblscratch;             /* double- scratch space                 */
} PGAScratch;


/*****************************************
*          CONTEXT STRUCTURE             *
*****************************************/
typedef struct str_PGAContext {
    PGAAlgorithm  ga;
    PGAOperations ops;
    PGAReport     rep;
    PGASystem     sys;
    PGADebug      debug;
    PGAInitialize init;
    PGAScratch    scratch;
} PGAContext;

/*****************************************
*          binary.c
*****************************************/

void PGABinaryCreateString(PGAContext *ctx, int p, int pop, int initflag) ;
void PGABinaryRandomize(PGAContext *ctx, int p, int pop);
int PGAGetBinaryAllele ( PGAContext *ctx, int p, int pop, int i );
void PGASetBinaryAllele ( PGAContext *ctx, int p, int pop, int i, int val );
int PGABinaryMutation( PGAContext *ctx, int p, int pop, double mr );
void PGABinaryOneptCrossover(PGAContext *ctx, int p1, int p2, int pop1, int c1,
			     int c2, int pop2);
void PGABinaryTwoptCrossover(PGAContext *ctx, int p1, int p2, int pop1, int c1,
			     int c2, int pop2);
void PGABinaryUniformCrossover(PGAContext *ctx, int p1, int p2, int pop1,
			       int c1, int c2, int pop2);
int PGABinaryHammingDistance ( PGAContext *ctx, PGABinary *s1, PGABinary *s2 );
void PGABinaryWriteString( PGAContext *ctx, FILE *fp, int p, int pop );
void PGABinaryWrite( PGAContext *ctx, FILE *fp, PGABinary *chrom, int nb );
void PGABinaryCopyString (PGAContext *ctx, int p1, int pop1, int p2, int pop2);
int PGABinaryDuplicate( PGAContext *ctx, int p1, int pop1, int p2, int pop2);
double PGAGetRealFromBinary(PGAContext *ctx, int p, int pop, int start,
			    int end, double lower, double upper);
double PGAGetRealFromGrayCode(PGAContext *ctx, int p, int pop, int start,
				  int end, double lower, double upper);
void PGAEncodeRealAsBinary(PGAContext *ctx, int p, int pop, int start,
			       int end, double low, double high, double val);
void PGAEncodeRealAsGrayCode(PGAContext *ctx, int p, int pop, int start,
			      int end, double low, double high, double val);
double PGAMapIntegerToReal (PGAContext *ctx, int v, int a, int b, double l,
			    double u);
int PGAMapRealToInteger(PGAContext *ctx, double r, double l, double u, int a,
			int b);
void PGAEncodeIntegerAsBinary(PGAContext *ctx, int p, int pop, int start,
			      int end, int val);
void PGAEncodeIntegerAsGrayCode(PGAContext *ctx, int p, int pop, int start,
				int end, int val);
int PGAGetIntegerFromBinary(PGAContext *ctx, int p, int pop, int start,
				 int end);
int PGAGetIntegerFromGrayCode(PGAContext *ctx, int p, int pop, int start,
				   int end);

/*****************************************
*          crossover.c
*****************************************/

void PGACrossover ( PGAContext *ctx, int p1, int p2, int pop1,
                    int c1, int c2, int pop2 );

/*****************************************
*          debug.c
*****************************************/

#ifndef OPTIMIZE
void PGADebugPrint( PGAContext *ctx, int level, char *funcname,
		   char *msg, int datatype, void *data );
#endif
int PGAGetDebugFlag(PGAContext *ctx, char *funcname);
void PGASetSupportingDebugFlags(PGAContext *ctx);

/*****************************************
*          duplicate.c
*****************************************/

int PGADuplicate(PGAContext *ctx, int p, int pop1, int pop2, int n);
void PGAChange( PGAContext *ctx, int p, int pop );

/*****************************************
*          fitness.c
*****************************************/

void PGAEvaluate ( PGAContext *ctx, int pop,
                   double (*f)(PGAContext *c, int p, int pop)  );
void PGASetEvaluate ( PGAContext *ctx, int p, int pop, double val );
void PGASetEvalUpToDate ( PGAContext *ctx, int p, int pop, int status );
void PGAFitness ( PGAContext *ctx, int popindex );
void PGAFitnessLinNor ( PGAContext *ctx, PGAIndividual *pop );
void PGAFitnessLinRank ( PGAContext *ctx, PGAIndividual *pop );
void PGAFitnessMinRecprl ( PGAContext *ctx, PGAIndividual *pop );
void PGAFitnessMinCmax ( PGAContext *ctx, PGAIndividual *pop );
double PGAMean ( PGAContext *ctx, double *a, int n);
double PGAStddev ( PGAContext *ctx, double *a, int n, double mean);
int PGARank( PGAContext *ctx, int p, int *order, int n );

/*****************************************
*          get.c
*****************************************/

int PGAGetDataType (PGAContext *ctx);
int PGAGetOptDir (PGAContext *ctx);
int PGAGetPopSize (PGAContext *ctx);
int PGAGetStringLen (PGAContext *ctx);
int PGAGetStoppingRule (PGAContext *ctx);
int PGAGetMaxIter (PGAContext *ctx);
int PGAGetNumReplace (PGAContext *ctx);
int PGAGetCrossoverType (PGAContext *ctx);
int PGAGetSelectType (PGAContext *ctx);
int PGAGetFitnessType (PGAContext *ctx);
int PGAGetFitnessMinType (PGAContext *ctx);
int PGAGetMutationType (PGAContext *ctx);
int PGAGetMutateOnlyNoCross (PGAContext *ctx);
double PGAGetMutateRealVal (PGAContext *ctx);
int PGAGetMutateIntegerVal (PGAContext *ctx);
int PGAGetNoDuplicates (PGAContext *ctx);
double PGAGetMutationProb (PGAContext *ctx);
double PGAGetCrossoverProb (PGAContext *ctx);
double PGAGetUniformCrossProb (PGAContext *ctx);
double PGAGetFitnessRankMax (PGAContext *ctx);
int PGAGetPopReplace (PGAContext *ctx);
int PGAGetIter (PGAContext *ctx);
int PGAGetMaxInt (PGAContext *ctx);
int PGAGetMinInt (PGAContext *ctx);
double PGAGetMaxDouble (PGAContext *ctx);
double PGAGetMinDouble (PGAContext *ctx);
int PGAGetPrintFreq (PGAContext *ctx);
int PGAGetRandomInit (PGAContext *ctx);
double PGAGetBinaryInitProb (PGAContext *ctx);
int PGAGetIntegerType (PGAContext *ctx);
int PGAGetInitIntegerMin (PGAContext *ctx, int i);
int PGAGetInitIntegerMax (PGAContext *ctx, int i);
double PGAGetInitRealMin (PGAContext *ctx, int i);
double PGAGetInitRealMax (PGAContext *ctx, int i);
int PGAGetRandomSeed(PGAContext *ctx);
double PGAGetEvaluate ( PGAContext *ctx, int p, int pop );
double PGAGetFitness ( PGAContext *ctx, int p, int pop );
int PGAGetEvalUpToDate ( PGAContext *ctx, int p, int pop );
int PGAGetWorst(PGAContext *ctx, int pop);
int PGAGetBest(PGAContext *ctx, int pop);
double PGAGetPTournamentProb(PGAContext *ctx);

/*****************************************
*          hamming.c
*****************************************/

double PGAHammingDistance( PGAContext *ctx, int popindex);

/*****************************************
*          heap.c
*****************************************/

void PGASortPop ( PGAContext *ctx, int pop );
int PGAGetSortPop ( PGAContext *ctx, int n );
void PGADblHeapSort ( PGAContext *ctx, double *a, int *idx, int n );
void PGADblHeapify( PGAContext *ctx, double *a, int *idx, int n );
void PGADblAdjustHeap ( PGAContext *ctx, double *a, int *idx, int i, int n );
void PGAIntHeapSort ( PGAContext *ctx, int *a, int *idx, int n );
void PGAIntHeapify( PGAContext *ctx, int *a, int *idx, int n );
void PGAIntAdjustHeap ( PGAContext *ctx, int *a, int *idx, int i, int n );

/*****************************************
*          integer.c
*****************************************/

void PGAIntegerCreateString (PGAContext *ctx, int p, int pop, int InitFlag);
void PGASetIntegerInitPermute ( PGAContext *ctx, int min, int max);
void PGASetIntegerInitLU (PGAContext *ctx, int *min, int *max);
void PGAIntegerRandomize(PGAContext *ctx, int p, int pop);
int PGAGetIntegerAllele (PGAContext *ctx, int p, int pop, int i);
void PGASetIntegerAllele (PGAContext *ctx, int p, int pop, int i, int value);
void PGAIntegerWriteString ( PGAContext *ctx, FILE *fp, int p, int pop);
void PGAIntegerCopyString (PGAContext *ctx, int p1, int pop1, int p2, int pop2);
int PGAIntegerMutation( PGAContext *ctx, int p, int pop, double mr );
void PGAIntegerOneptCrossover(PGAContext *ctx, int p1, int p2, int pop1,
			      int c1, int c2, int pop2);
void PGAIntegerTwoptCrossover( PGAContext *ctx, int p1, int p2, int pop1,
			      int c1, int c2, int pop2);
void PGAIntegerUniformCrossover(PGAContext *ctx, int p1, int p2, int pop1,
				int c1, int c2, int pop2);
int PGAIntegerDuplicate( PGAContext *ctx, int p1, int pop1, int p2, int pop2);

/*****************************************
*          mutation.c
*****************************************/

int PGAMutate(PGAContext *ctx, int p, int pop);

/*****************************************
*          pga.c
*****************************************/

void PGARun(PGAContext *ctx,
	    double (*evaluate)(PGAContext *c, int p, int pop));
void PGARunMutateAndCross (PGAContext *ctx, int oldpop, int newpop);
void PGARunMutateOrCross ( PGAContext *ctx, int oldpop, int newpop );
PGAIndividual *PGAGetIndividual ( PGAContext *ctx, int p, int pop);
void PGACopyIndividual( PGAContext *ctx, int p1, int pop1, int p2, int pop2);
void PGAUpdateGeneration ( PGAContext *ctx );
int PGADone( PGAContext *ctx);
PGAIndividual *PGAPopAddress( PGAContext *ctx, int popix );
void PGADestroy (PGAContext *ctx);
PGAContext *PGACreate ( int *argc, char **argv,
                        int datatype, int len, int maxormin );
void PGASetUp ( PGAContext *ctx );
void PGACreatePop (PGAContext *ctx, int pop);
void PGACreateIndividual (PGAContext *ctx, int p, int pop, int initflag);
void PGAUpdateOnline(PGAContext *ctx, int pop);
void PGAUpdateOffline(PGAContext *ctx, int pop);
int PGAComputeSimilarity(PGAContext *ctx, PGAIndividual *pop);

/*****************************************
*          real.c
*****************************************/

void PGARealCreateString (PGAContext *ctx, int p, int pop, int initflag);
void PGASetRealInitPercent ( PGAContext *ctx, double *median, double *percent);
void PGASetRealInitLU (PGAContext *ctx, double *min, double *max);
void PGARealRandomize ( PGAContext *ctx, int p, int pop);
double PGAGetRealAllele (PGAContext *ctx, int p, int pop, int i);
void PGASetRealAllele (PGAContext *ctx, int p, int pop, int i, double value);
void PGARealWriteString (PGAContext *ctx, FILE *fp, int p, int pop);
void PGARealCopyString ( PGAContext *ctx, int p1, int pop1, int p2, int pop2);
int PGARealMutation( PGAContext *ctx, int p, int pop, double mr );
void PGARealOneptCrossover( PGAContext *ctx, int p1, int p2, int pop1,
			   int c1, int c2, int pop2);
void PGARealTwoptCrossover( PGAContext *ctx, int p1, int p2, int pop1,
			   int c1, int c2, int pop2);
void PGARealUniformCrossover( PGAContext *ctx, int p1, int p2, int pop1,
			     int c1, int c2, int pop2);
int PGARealDuplicate( PGAContext *ctx, int p1, int pop1, int p2, int pop2);

/*****************************************
*          report.c
*****************************************/

void PGAPrintPopulation ( PGAContext *ctx, FILE *fp, int pop );
void PGAPrintIndividual ( PGAContext *ctx, FILE *fp, int p, int pop );
void PGAPrintReport(PGAContext *ctx, FILE *fp, int pop);
void PGAPrintContext ( PGAContext *ctx, FILE *fp );

/*****************************************
*          select.c
*****************************************/

void PGASelect( PGAContext *ctx, int popix );
int PGASelectProportional(PGAContext *ctx, PGAIndividual *pop);
void PGASelectSUS( PGAContext *ctx, PGAIndividual *pop );
int PGASelectTournament( PGAContext *ctx, PGAIndividual *pop );
int PGASelectPTournament( PGAContext *ctx, PGAIndividual *pop );
int PGASelectNext ( PGAContext *ctx );

/*****************************************
*          set.c
*****************************************/

void PGASetPopSize (PGAContext *ctx, int popsize);
void PGASetStoppingRule (PGAContext *ctx, int stoprule);
void PGASetMaxIter(PGAContext *ctx, int maxiter);
void PGASetMaxNoChange(PGAContext *ctx, int max_no_change);
void PGASetMaxSimilarity(PGAContext *ctx, int max_similarity);
void PGASetNumReplace( PGAContext *ctx, int pop_replace);
void PGASetCrossoverType (PGAContext *ctx, int crossover_type);
void PGASetCrossoverProb( PGAContext *ctx, double crossover_prob);
void PGASetUniformCrossProb( PGAContext *ctx, double uniform_cross_prob);
void PGASetSelectType( PGAContext *ctx, int select_type);
void PGASetPTournamentProb(PGAContext *ctx, double ptournament_prob);
void PGASetFitnessType( PGAContext *ctx, int fitness_type);
void PGASetFitnessMinType( PGAContext *ctx, int fitness_type);
void PGASetFitnessRankMax( PGAContext *ctx, double fitness_rank_max);
void PGASetMutationType( PGAContext *ctx, int mutation_type);
void PGASetMutateOnlyNoCross( PGAContext *ctx, int mutate_only_no_cross);
void PGASetMutationRealVal( PGAContext *ctx, double val);
void PGASetMutationIntegerVal( PGAContext *ctx, int val);
void PGASetNoDuplicates( PGAContext *ctx, int no_dup);
void PGASetMutationProb(PGAContext *ctx, double mutation_prob);
void PGASetPopReplacement( PGAContext *ctx, int pop_replace);
void PGASetPrintFreq( PGAContext *ctx, int print_freq);
void PGASetRandomInit(PGAContext *ctx, int RandomBoolean);
void PGASetDebugLevel(PGAContext *ctx, int level);
void PGASetPrintOptions (PGAContext *ctx, int option);
void PGASetBinaryInitProb ( PGAContext *ctx, double probability );
void PGASetRandomSeed(PGAContext *ctx, int seed);
void PGASetUserFunction(PGAContext *ctx, int constant, void *f);

/*****************************************
*          system.c
*****************************************/

void PGAError( PGAContext *ctx, char *msg,
               int level, int datatype, void *data );
void  PGAExit( PGAContext *ctx, int exittype );
void PGAWriteString ( PGAContext *ctx, FILE *fp, int p, int pop );
int PGARandomFlip ( PGAContext *ctx, double p );
int PGARandomInterval( PGAContext *ctx, int start, int end);
double PGARandom01( PGAContext *ctx, int newseed );
double PGARandomUniform( PGAContext *ctx, double start, double end);
double PGARandomGaussian( PGAContext *ctx, double mean, double sigma);
int PGARound(PGAContext *ctx, double x);
void PGAReadCmdLine( PGAContext *ctx, int *argc, char **argv );
void PGAParseDebugArg(PGAContext *ctx, char * st);
void PGAStripArgs(char **argv, int *argc, int *c, int num);
void PGAUsage( PGAContext *ctx );
void PGAPrintVersion( PGAContext *ctx );
void PGAPrintDebugOptions(PGAContext *ctx);
