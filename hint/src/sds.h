#include <math.h>

#define TRUE 1
#define FALSE 0
#define DESC_COL 20
#define CUNDEF 127
#define CCLASH 101
#define CLOG2  0.69314718056

#define M_NONE    0
#define M_INFORM  1 
#define M_GINI    2
#define M_GR      3
#define M_RELIEFF 4
#define M_MDL 5

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#ifdef sun
#define LOG2(x) log2((double)(x))
#else
#define LOG2(x) log((double)(x))/CLOG2
#endif

#define MMAX(x,y) ((x)<(y)?(y):(x))
#define MMIN(x,y) ((x)<(y)?(x):(y))
#define XOR(a,b) ((a)&&(!(b)) || (!(a))&&(b))
#define ABS(x) ((x) < 0 ? -(x) : x)
#define FREE(x) {if (x==NULL) free(x); x=NULL; }

#ifdef sun
#define NLOGN(x) (((x)==0)?0:((x) * log2((double)x)))
#else
#define NLOGN(x) (((x)==0)?0:((x) * log((double)x) / CLOG2))
#endif

#define FIND_VAR(x,s) {if (!((x=find_var(variables, s))!=NULL)) { \
		       printf("error: variable %s not found\n", s); \
		       return; \
		      }}
#define FIND_TABLE(t,s) {if ((t=find_table(s)) == NULL) { \
                        printf("error: table %s node defined\n", s); \
                        return; \
		       }}

/****************************************************************************
TYPE DEFINITIONS
****************************************************************************/

typedef char Str255[255];
typedef char Str15[15];

typedef struct list_of_vars list_of_vars;
typedef struct fams fams;
typedef struct list_of_fams list_of_fams;
typedef struct list_of_str list_of_str;
typedef struct list_of_num list_of_num;
typedef struct desc_t desc_t;
typedef struct list_of_q list_of_q;
typedef struct list_of_opt list_of_opt;
typedef struct var_type var_type;
typedef struct list_of_itables list_of_itables;
typedef struct rule_list rule_list;
typedef struct dec_tree dec_tree;

typedef enum {ct_nominal, ct_contin} ctype_type;
typedef enum {regular, left, right, none} trapz_type;
typedef enum {ts_sigm, ts_lin} ts_mem_type;
typedef enum {col_optimal,col_heuristic,col_ga,col_set,col_soft} coloring_type;

struct var_type {		/* DEFINITION OF A VARIABLE */
  Str255 name;			/* name */
  ctype_type ctype;		/* type of the variable */
  int ndesc;			/* number of descriptors */
  desc_t *desc;			/* array of qualitative descriptions */

  double *opt;			/* array of len=ndesc of degrees of mem */
  char class;			/* class of this variable */

  double val;			/* numerical val of current option */
  char valdef;			/* is num value defined? */
  double expect;		/* expected num value */
  char expect_def;		/* is expected defined? */
  fams *famsout;		/* fam that is used to compute variable */
  list_of_fams *famsin;		/* fams that var is an input variable to */
  char mark;			/* a mark needed for derive */
  double infor;			/* informativity */
  double gain;			/* the gain of the attribute */
  double vmin, vmax;		/* boundaries in double */
  int des_min, des_max;		/* descriptors cannot cross this boundary */
				/* used withing GA to learn descriptors */
  char c_default;		/* default class */
  double *c_apriory;		/* apriory probabilities of classes */
  char decomposable;		/* is famsout of this var decomposable? */

  struct {			/* TS model */
    var_type **var;		/* independent variables */
    int d;			/* their number */
    double *p;			/* parameters for mem, d x k x 4 */
    ts_mem_type type;		/* membership function type */
    int k;			/* number of rules/clusters */
    double *a;			/* a coeficients, k * d  */
    double *b;			/* b coeficients, k */
  } ts;
};

struct list_of_vars {		/* LIST OF VARIABLES */
  var_type *var;		/* a pointer to a variable */
  struct list_of_vars *prev, *next; /* a pointer to next and prev el of list */
};

struct list_of_opt {		/* LIST OF OPTIONS */
  Str255 name;			/* option's name */
  double *val;			/* numerical values */
  char *valdef;			/* is num value defined? */
  double **opt;			/* degrees of membership's */
  double *expect;		/* expected num value */
  char *expect_def;		/* is expected defined? */
  char mark;
  struct list_of_opt *prev, *next;
};

struct list_of_itables {		/* LIST OF INSTANCE TABLES */
  Str255 name;
  int n_in;			/* number of in variables */
  int n_inst;			/* number of instances */

  var_type **in;		/* input variables */
  var_type *out;		/* output variable */

  double **val;			/* continous value, [n_inst][n_in+1] */
  char **qval;			/* qualitative value as defined from input, 
				   CUNDEF if not defined, [n_inst][n_in+1] */

  char **tval;			/* a mapping from cont to qval using 
				   intervals, while deriving rules */

  char *estdef;			/* estimated output value defined? */
  double *est;			/* estimated output */
  char *estq;			/* estimated qualitative output */

  char *mark;			/* is instance already used? */
  struct list_of_itables *next;
};

struct desc_t {			/* VARIABLE'S DESCRIPTION ENTRY */
  Str255 name;			/* qualitative mark */
  double apriory;		/* apriory probability */
				/* INTERVAL */
  double start, delta;		/* its start and width */
  char int_def;			/* is interval defined? */
				/* FUZZY */
  trapz_type tp;		/* trapezoidal type */
  double a, b, c, d;		/* vertices of trapezoid */
  double cent;			/* centeroid (of mass) */
  double m;			/* mass */
  double l_slope, r_slope;	/* two slopes */
  double *cap;			/* apriory probability for given class*/
};


struct list_of_str {		/* LIST OF STRINGS */
  Str15 str;			/* string */
  struct list_of_str *prev, *next; /* prev and next element in a list */
};

struct list_of_num {		/* LIST OF NUMBERS */
  double num;			/* number */
  struct list_of_num *prev, *next; /* prev and next element in a list */
};

typedef enum {undef, man, autom, manfix, autofix} rule_type;
				/* undef - not defined yet */
				/* man - manually defined */
				/* auto - automatically induced */
				/* *fix - fixed not to alter when GA */
typedef enum {er_table, er_list, er_tree} er_type;

typedef struct hash_entry hash_entry;
struct hash_entry {
  rule_list *rl;
  hash_entry *next;
};

struct fams {			/* FAMS */
  Str255 name;			/* FAM's name */
  var_type **in;		/* array of pointers to input variables */
  var_type *out;		/* pointer to output variable */
  int n_in;			/* number of input variables */

  int n_rules;			/* number of q rules */
  int n_unfixed;		/* number of rules not fixed */

  er_type er;			/* how the rules are encoded? */

				/* if rules are defined as table */
  char *rule;			/* qualitative mapping of input to output */
  double *imp;			/* importance of rule, [0,1] */
  int *usage;			/* # rule was used in dec process */
  rule_type *rtp;		/* type of rule in FAM */

  rule_list *lrule;		/* rules are defined as list */
  dec_tree *dt;			/* decision tree for these rules */

  double *degrees;		/* used for evaluation, degrees of out */
				/* variable or this fam */
  fams *prev, *next;		/* pointers to prev and next el of list */

  hash_entry **hash_table;	/* hash table */  
};

struct list_of_fams {		/* LIST OF FAMS */
  fams *fam;			/* pointer to fam */
  list_of_fams *next, *prev;	/* next and prev el of list */
};

struct rule_list {		/* LIST OF RULES */
  char *att;			/* attributes */
  char class;			/* corresponding class */
  double *dist;			/* class probabilities P(c_i) */
  rule_type rtp;		/* rule type */
  double imp;			/* importance of rule */
  int usage;			/* its usage */
  int n;			/* index of this rule */
  rule_list *next, *prev;	/* next rule in the list */

  char unused;			/* flag for decomposition */
  char train, test;		/* to which set does rule belong? nonexcl */
  char mark;			/* used when sorting, indicates last in row */
  int cr;			/* which colon/row belongs to */
};

struct dec_tree {		/* DECISION TREE */
  char is_leaf;			/* is node a leaf? */
  char class;			/* class (if leaf) */
  var_type *v;			/* variable in internal node */
  fams *fam;			/* fam with rules to decompose */
  int nrules;			/* number of defined rules used */
  dec_tree **next;		/* successors */
};

struct list_of_q {		/* LIST OF QVALS WITH DEGREES */
  Str255 desc;			/* name of descriptor */
  double degree;		/* membership degree */
  list_of_q *next;
};

/* what to list when DAG of variables is displayed */
typedef enum {desc, curropt, opt, infor, empty} ltype;

GL enum {latex, graph} plot_type;	/* type of plotting to be used */

/****************************************************************************
GLOBAL VARIABLES
****************************************************************************/

typedef enum {famin, famax, famult} fa_type ;
typedef enum {gae_sqr, gae_perc, gae_sqrperc, gae_norm,
	      gae_maxperc, gae_abs} ga_error_type;
typedef enum {e_crisp, e_fuzzy, e_numerical, e_interval} eval_method_type;
typedef enum {c_cm, c_r, c_cr, c_sic, c_dfc, c_error, c_m} dec_crit_type;
typedef enum {n_laplace, n_cluster, n_entropy, n_mprob} noise_handling_type;
typedef enum {m_cm, m_error, m_cm_error} m_dec_type;
typedef enum {t_cv, t_sample} test_type_type;
typedef enum {dc_dont_care, dc_dont_know, dc_apriory} dont_care_type;

GL fa_type fa_method;		/* how to determine */
				/* the rule strength from the premises */

GL noise_handling_type noise_handling; /* noise handling type */
GL m_dec_type m_dec;		/* when using m, what is part sel measure? */
GL test_type_type test_type;	/* how to test the class error for m-param */
GL dont_care_type dont_care;	/* don't care handling */

GL list_of_vars *variables;	/* root of the list of variables */
GL list_of_fams *tables;	/* root of the list of tables */
GL list_of_itables *itables;	/* root of list of instance tables */
GL list_of_opt *options;	/* options */

GL int n_dis_lo, n_dis_hi;	/* decomposition, size of bound set */
GL int n_ndis_lo, n_ndis_hi;	/* # non-disjoint vars of free and bound set */

GL dec_crit_type dec_crit;	/* what is the criteria for partition select */
GL char dec_global;		/* global decomposition? */

GL char debug_l;		/* lexical debug, display what lex reads in */
GL char debug_e;		/* debug evaluation process */
GL char debug_g;		/* GA debugging */
GL char debug_d;		/* check descriptors consistency */
GL char print_short;		/* long or short descriptions for output */
				/* also tells how many chars to print BBB */
GL char save_im;		/* for coloring seminar - do we save im */
GL Str255 im_fname;		/* file name for saving im */

GL char interact_mode;		/* is reading from file or tt */
GL int n_errors;		/* # errors detected while parsing */
GL fams *cfam;			/* current fam, to be set when add rules */
GL char pref_qual;		/* qual/quant preference for entries */
GL eval_method_type eval_method; /* evaluation method */
GL int repeat_tests;		/* number of repetition for all tests */
GL char save_tests;		/* save the tests to c4.5 like file */
GL er_type encode_rules;	/* how are rules encoded? */
GL char use_dm;			/* use decomposition matrix */

GL int g_argc;			/* copy of argc */
GL char **g_argv;
GL FILE *lfile;			/* log file */
GL coloring_type coloring;	/* what algorithm to use to color graph */
GL double mcriteria;		/* m-criteria, Cestnik */

				/* GENETIC ALGORITHM PARAMETERS */
GL ga_error_type ga_error_method;	/* an arror measure used in ga */
GL int ga_maxiter;		/* max number of iterations */
GL int ga_population_size;	/* population size */
GL int ga_desc_pts;		/* # num landmarks for membership function */
GL int ga_print_freq;		/* debug print frequency */
GL int ga_weight_pts;		/* # points the weights can be at */
GL double ga_width;		/* description enlargment */
GL double ga_max_error;		/* max error to be displayed on plots */

GL Str255 prefix_var, prefix_opt; /* generic prefix for vars and opts */
GL int iprefix_var, iprefix_opt; /* their numeric extensions */

typedef enum {gap_and, gap_or, gap_custom} ga_policy_type;
GL char ga_policy;		/* mutate and cross policy */

				/* mutations and crossover probabilities */
GL double ga_mut_fams, ga_mut_desc, ga_mut_w;
GL double ga_cross_fams, ga_cross_desc, ga_cross_w;
GL char ga_learn_fams, ga_learn_desc, ga_learn_w;

extern FILE *yyin;		/* input file for yacc */
extern int yylineno;

#define MY_MAX_FILES 50		/* max files that can be opened */
GL struct {			/* save line number and f pointer */
  FILE *f;			/* of the files being loaded */
  int lineno;
} myfiles[MY_MAX_FILES];
GL int myno_files;		/* number of files loading */

				/* FUZZY IDENTIFICATION */
GL int cl_K;			/* number of clusters, initial */
GL char deb_cl_cp;		/* debug cluster prototypes */
GL char deb_cl_pm;		/* debug partition matrix */
GL char deb_cl_f;		/* debug F */
GL double cl_gamma;		/* merging treshold */
GL double cl_m;			/* fuzziness */
GL double cl_e;			/* convergence limit */
GL int cl_redo;			/* loops to redo the algorithm */
GL char cl_sigm;		/* sigmoid memberships? */

				/* FUNCTION DECOMPOSITION */
GL int deb_dec;			/* debug, ie show decomposition matrix */
GL int krelieff;		/* k in NN algorithm for relieff */
GL char heur_comp, heur_split;	/* compute heuristics, use it to split */
GL char log_inform, log_gini, log_gr, log_relieff, log_mdl;
GL char deb_inform, deb_gini, deb_gr, deb_relieff, deb_mdl;
GL int nrow, ncol;		/* number of rows/cols decomposition matrix */
GL int nvrow, nvcol;		/* number of vars in decomposition matrix */
GL char *drow, *dcol;		/* index(color) for different rows/columns */
GL double m_param;		/* m-parameter of Cestnik & Bratko */
GL char use_distribution;	/* use distribution for classes */

/****************************************************************************
PROTOTYPES
****************************************************************************/

list_of_vars *add_var(list_of_vars **root, Str255 name, char b, int ndesc);
var_type *find_var(list_of_vars *root, Str255 name);
var_type *find_var_num(list_of_vars *root, Str255 name, int *i);
list_of_itables *find_itable(Str255 name);
void free_list_of_vars(list_of_vars *lv, char dest);

void add_fdesc_var(Str255 dname, Str255 vname, trapz_type tp,
		   double a, double b, double c, double d);
fams *find_fam(Str255 vname, Str255 fname);
fams *find_fam_vartable(char *idname);
fams *find_table(Str255 name);
list_of_opt *find_opt(Str255 oname, char d);
void mark_var(var_type *v, char m);
int find_var_pos(list_of_vars *v, Str255 vname);
void select_opt(Str255 oname);
void derive_var(var_type *v);
void set_opt(Str255 oname);
void add_fdesc(desc_t *des, trapz_type tp,
	       double a, double b, double c, double d);
void real_to_fuzzy(var_type *v);
void list_struct(list_of_vars *root, ltype lt);
double get_opt_error(double ex, double der);
char *get_q_val(var_type *v);
double rnd1e();
double rnd1();
char *strn(char *s, int n);
char *strnc(char *s, int n);
double get_opt_error(double ex, double der);
double stat_for_vnum(int vnum);

void indx2att(fams *f, int indx, char *att);
int att2indx(fams *f, char *att);

list_of_vars *find_leaves(var_type *v);

/* memory allocation */

int *i_vector(int n);
char *c_vector(int n);
double *d_vector(int n);
int **i_matrix(int nr, int nc);
char **c_matrix(int nr, int nc);
double **d_matrix(int nr, int nc);

int *i_vector_ini(int n, int ini);
double *d_vector_ini(int n, double ini);
int **i_matrix_ini(int nr, int nc, int ini);
char *c_vector_ini(int n, char ini);
char **c_matrix_ini(int nr, int nc, char ini);
double d_vector_max(double *v, int len, int *pos);
double d_vector_max_max(double *v, double *v1, int len, int *pos);
int i_vector_max(int *v, int len, int *pos);
double d_vector_sum(double *v, int len);

double bin_compare_pos_heur_real(int **cr, double **h, int n);
double bin_compare_splits_heur_real(int **cr, double **h, int n);
char copy_table(Str255 tname, Str255 vname, double perc);
char *gen_opt_name();
char *gen_var_name();
double sstat_for_itable(list_of_itables *it);
double stat_for_itable(Str255 iname);
char get_nrule_e(fams *f, int n, rule_type *rtp);
char get_nrule(fams *f, int n);
int count_rules(fams *f);

				/* misc.c */
void indx2att(fams *f, int indx, char *att);
void indx2satt(fams *f, int indx, char *att, char *sel);
int att2indx(fams *f, char *att);
int satt2indx(fams *f, char *att, char *sel);
void free_c_matrix(char **m, int row);
void free_i_matrix(int **m, int row);

void test_hash(Str255 name);
int get_hrule(fams *f, char *att);
int hash(fams *f, char *att);
rule_list *find_hrule(fams *f, char *att);

char get_next_rule(fams *f, char *class, rule_type *rtp);
char get_next_rulea(fams *f, char *class, rule_type *rtp, char **att);

double get_code_measure(Str255 vname);
double logfact(int i);

double get_class_error();
int rule_cmp_rows(rule_list *r1, rule_list *r2);
double compare_struct_dist(char *id1, char *id2);
