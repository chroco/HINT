#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#define GL extern
#include "sds.h"

void usage(char *progname, char *str)
{
  if (str != NULL)
    printf("%s", str);
  printf("Usage: %s [options]\n", progname);
  printf(" -l    lexical debugging\n");
  printf(" -g    ga debugging\n");
  printf(" -e    evaluation debugging\n");
  printf(" -d    check descriptor consistency\n");
  printf(" -c f  read cmds from command file f\n");
  exit(0);
}

void init_init(int argc, char *argv[])
{
  lfile = fopen("log","w");
  if (lfile==NULL) printf("warning: cannot open lfile file\n");
  
  g_argc = argc;
  g_argv = argv;
  variables = NULL;
  itables = NULL;
  tables = NULL;
  options = NULL;
  n_errors = 0;
  cfam = NULL;
  myno_files = -1;
  plot_type = graph;
  print_short = 1;
  pref_qual = TRUE;
  eval_method = e_crisp;
  drow=NULL, dcol=NULL;

  fa_method = famin;
  repeat_tests = 10;
  save_tests = FALSE;

  /* parameters, that can be set from the script */

  use_distribution = FALSE;
  m_param = 0.;
  encode_rules = er_list;
  strcpy(prefix_var,"c");
  strcpy(prefix_opt,"o");
  iprefix_var = iprefix_opt = 1;
  use_dm = FALSE;
  dec_crit = c_cm;
  dec_global = FALSE;

  debug_l = debug_e = debug_g = debug_d = FALSE;

  n_dis_lo = n_dis_hi = 2;
  n_ndis_lo = n_ndis_hi = 0;

  noise_handling = n_laplace;
  m_dec = m_cm;
  test_type = t_sample;
  dont_care = dc_dont_know;

  ga_maxiter = 100;
  ga_population_size = 100;
  ga_desc_pts = 50;
  ga_max_error = 500;
  ga_print_freq = 1;
  ga_weight_pts = 10;
  ga_error_method = gae_abs;
  ga_policy = gap_or;
  ga_width = 1.5;

  ga_mut_fams = ga_mut_desc = ga_mut_w = 1;
  ga_cross_fams = ga_cross_desc = ga_cross_w = 0.85;
  ga_learn_fams = ga_learn_desc = TRUE; ga_learn_w = FALSE;

  /* fuzzy identification */
  cl_K = 5;
  deb_cl_cp = FALSE;
  deb_cl_pm = FALSE;
  deb_cl_f = FALSE;
  cl_m = 2.5;			/* fuzziness */
  cl_gamma = 0.1;		/* merging treshold */
  cl_e = 0.05;			/* clustering loop */
  cl_redo = 1;			/* loops to redo the algorithm */
  cl_sigm = FALSE;

  /* decomposition */
  coloring = col_heuristic;
  mcriteria = 0.;
  deb_dec = 0;
  krelieff = 10;
  heur_comp = heur_split = FALSE;
  log_inform = log_gini = log_gr = log_mdl = TRUE;
  log_relieff = TRUE;
  deb_inform = deb_gini = deb_gr = deb_relieff = deb_mdl = FALSE;
}


char read_options(int argc, char *argv[])
{
  int i, j;
  char b;

  if (argc > 1)
    for (i=1; i<argc; i++) {
      if (argv[i][0] != '-') {
	/* instead of data file here should be the storage file BBB */
/*	if (data_file != NULL)
	  usage(argv[0], "error: input file specified more than once\n");
	data_file = fopen(argv[i], "r");
	if (data_file == NULL) {
	  printf("error: can't open %s\n", argv[i]);
	  exit(0);
	} */
      }
      
      for (j=1, b=TRUE; j<strlen(argv[i]) && b; j++) {
	switch (argv[i][j]) {
	case 'h':
	  usage(argv[0], NULL);
	  break;
	case 'l':
	  debug_l = TRUE;
	  break;
	case 'e':
	  debug_e = TRUE;
	  break;
	case 'g':
	  debug_g = TRUE;
	  break;
	case 'd':
	  debug_d = TRUE;
	  break;
	case 'c':
	  yyin = fopen(argv[++i], "r");
	  if (yyin == NULL)
	    usage(argv[0], "error: could not open cmd file\n");
	  else
	    fprintf(lfile,"FILE %s\n", argv[i]);
	  b = FALSE;
	break;
	}
      }
    }
}

/* load_cmds: loads the commands, so that it takes care for the files
   that is currently being loaded and those from which the call load
   was made. A list of uncompletely loaded files is maintained by
   myfiles. */

void load_cmds(char mode)
{
  interact_mode = mode;
  yyparse();

  while (myno_files != -1) {
    fclose(yyin);
    yyin = myfiles[myno_files].f;
    ungetc('\n',yyin);
    yylineno = myfiles[myno_files].lineno;
    myno_files--;
    if (mode && myno_files == -1) {
      interact_mode = TRUE;
      ungetc('\n',yyin);
    }
    yyparse();
  }
  fclose(yyin);      
}

/****************************************************************************
Usefull rutines (random, strn, strnc)
****************************************************************************/

/* rnd1: returns a random number [0..1] */

double rnd1()
{
  double d;
#ifdef hp
  d = (double) rand()/32767.0;
#else
/*  return (double) rand()/32767.0; */
  d = (double) rand()/2147483647.0;
#endif
  if (d<0 || d>1.0) {
    printf("Error in function rnd1 in misc.c\n");
    printf("Change the division constant\n");
    exit(0);
  }
  return d;
}

/* rnd1e: returns a random number [0..1) BBBB */

double rnd1e()
{
  double d;
#ifdef hp
  d = rand()/32768.0;
#else
/*  return (double) rand()/32768.0; */
  d = rand()/2147483648.0;
#endif
  if (d<0 || d>=1.0) {
    printf("Error in function rnd1 in misc.c\n");
    printf("Change the division constant\n");
    exit(0);
  }
  return d;
}

Str255 rs;

char *strn(char *s, int n)
{
  int i;

  if (strlen(s)<n && FALSE) {
    strcpy(rs, s);
  }
  else {
    strncpy(rs, s, n);
    if (strlen(rs)<n) {
      for(i=strlen(rs); i<n; i++) rs[i]=' ';
      rs[i]='\0';
    }
    else rs[n]='\0';
    /*  printf("xxx %s\n", rs); */
  }
  return rs;
}

char *strnc(char *s, int n)
{
  int i, k;

  k = (n-MMIN(strlen(s),n))/2;
  for (i=0; i<k; i++) rs[i]=' ';
  strncpy(&rs[i], s, n);
  for (i=k+MMIN(strlen(s),n); i<n; i++) rs[i]=' ';
  rs[i]='\0';
  return rs;
}

void pspace(int n)
{
  int i;
  for (i=0; i<n; i++) printf(" ");
}

/****************************************************************************
Generation of names of variables and options
****************************************************************************/

Str255 xname;

char *gen_opt_name()
{
  char b=TRUE;
  while (b) {
    sprintf(xname, "%s%d", prefix_opt, iprefix_opt++);
    b = find_opt(xname, FALSE) != NULL;
  }
  return xname;
}

char *gen_var_name()
{
  char b=TRUE;
  while (b) {
    sprintf(xname, "%s%d", prefix_var, iprefix_var++);
    b = find_var(variables, xname) != NULL;
  }
  return xname;
}

/****************************************************************************
Memory allocation rutines
****************************************************************************/

void myfree(void *p)
{
  if (p!=NULL) { free(p); p=NULL; }
}

void free_str_list(list_of_str *s)
{
  list_of_str *s1;
  while (s!=NULL) {
    s1 = s->next;
    free(s);
    s = s1;
  }
}

int *i_vector(int n)
{
  int *v;
 
  v=(int *)malloc((unsigned) n*sizeof(int));
  if (!v) {
    printf("allocation failure in i_vector()\n");
    exit(0);
  }
  return v;
}
 
int *i_vector_ini(int n, int ini)
{
  int *v, i;
 
  v=(int *)malloc((unsigned) n*sizeof(int));
  if (!v) printf("allocation failure in i_vector()\n");
  else for (i=0; i<n; i++) v[i] = ini;
  return v;
}
 
char *c_vector(int n)
{
  char *v;
 
  v=(char *)malloc((unsigned) n*sizeof(char));
  if (!v) printf("allocation failure in i_vector()\n");
  return v;
}
 
char *c_vector_ini(int n, char ini)
{
  char *v;
  int i;
 
  v=(char *)malloc((unsigned) n*sizeof(char));
  if (!v) printf("allocation failure in i_vector()\n");
  else for (i=0; i<n; i++) v[i] = ini;
  return v;
}
 
double *d_vector(int n)
{
  double *v;
 
  v=(double *)malloc((unsigned) n*sizeof(double));
  if (!v) printf("allocation failure in d_vector()\n");
  return v;
}

double *d_vector_ini(int n, double ini)
{
  double *v;
  int i;
 
  v=(double *)malloc((unsigned) n*sizeof(double));
  if (!v) printf("allocation failure in d_vector()\n");
  for (i=0; i<n; i++) v[i] = ini;
  return v;
}

int **i_matrix(int nr, int nc)
{
  int i;
  int **m;
 
  m=(int **) malloc((unsigned) nr * sizeof(int*));
  if (!m) printf("allocation failure 1 in i_matrix() - %d\n", nr);
  for(i=0; i<nr; i++) {
    m[i]=(int *) malloc((unsigned) nc * sizeof(int));
    if (!m[i]) printf("allocation failure 2 in i_matrix()\n");
  }
  return m;
}

int **i_matrix_ini(int nr, int nc, int ini)
{
  int i, j;
  int **m;
 
  m=(int **) malloc((unsigned) nr * sizeof(int*));
  if (!m) printf("allocation failure 1 in i_matrix()\n");
  for(i=0; i<nr; i++) {
    m[i]=(int *) malloc((unsigned) nc * sizeof(int));
    if (!m[i]) printf("allocation failure 2 in i_matrix()\n");
  }
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
      m[i][j] = ini;
  return m;
}

char **c_matrix(int nr, int nc)
{
  int i;
  char **m;
 
  m=(char **) malloc((unsigned) nr * sizeof(char*));
  if (!m) printf("allocation failure 1 in c_matrix()\n");
  for(i=0; i<nr; i++) {
    m[i]=(char *) malloc((unsigned) nc * sizeof(char));
    if (!m[i]) printf("allocation failure 2 in c_matrix()\n");
  }
  return m;
}

char **c_matrix_ini(int nr, int nc, char ini)
{
  int i, j;
  char **m;
 
  m=(char **) malloc((unsigned) nr * sizeof(char*));
  if (!m) printf("allocation failure 1 in c_matrix()\n");
  for(i=0; i<nr; i++) {
    m[i]=(char *) malloc((unsigned) nc * sizeof(char));
    if (!m[i]) printf("allocation failure 2 in c_matrix() - %d x %d\n",nr,nc);
  }
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
      m[i][j] = ini;
  return m;
}

double **d_matrix(int nr, int nc)
{
  int i;
  double **m;
 
  m=(double **) malloc((unsigned) nr * sizeof(double*));
  if (!m) printf("allocation failure 1 in d_matrix()\n");
  for(i=0; i<nr; i++) {
    m[i]=(double *) malloc((unsigned) nc * sizeof(double));
    if (!m[i]) printf("allocation failure 2 in d_matrix()\n");
  }
  return m;
}

double **d_matrix_ini(int nr, int nc, double ini)
{
  int i, j;
  double **m;
 
  m=(double **) malloc((unsigned) nr * sizeof(double*));
  if (!m) printf("allocation failure 1 in d_matrix()\n");
  for(i=0; i<nr; i++) {
    m[i]=(double *) malloc((unsigned) nc * sizeof(double));
    if (!m[i]) printf("allocation failure 2 in d_matrix()\n");
  }
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
      m[i][j] = ini;  
  return m;
}

void free_c_matrix(char **m, int row)
{
  int i;

  if (m!=NULL) {
    for (i=0; i<row; i++)
      free(m[i]);
    free(m);
    m = NULL;
  }
}

void free_i_matrix(int **m, int row)
{
  int i;

  if (m!=NULL) {
    for (i=0; i<row; i++)
      free(m[i]);
    free(m);
    m = NULL;
  }
}

void c_matrix_copy(char **from, char ***to, int a, int b)
{
  int i, j, k;
  char **ttt;
  ttt = c_matrix(a,b);
  for (i=0; i<a; i++)
    for (j=0; j<b; j++)
      ttt[i][j] = from[i][j]; 
  *to = ttt;
}

void c_vector_copy(char *from, char **to, int a)
{
  int i;
  char *ttt;
  ttt = c_vector(a);
  for (i=0; i<a; i++)
    ttt[i] = from[i];
  *to = ttt;
}


int i_vector_max(int *v, int len, int *pos)
{
  int i, max = v[0];
  *pos = 0;
  for (i=1; i<len; i++)
    if (v[i]>max) {
      max = v[i];
      *pos = i;
    }
  return max;
}

void i_vector_set(int *v, int len, int val) 
{
  int i;
  for (i=0; i<len; i++) v[i] = val;
}

void c_vector_set(char *v, int len, char val) 
{
  int i;
  for (i=0; i<len; i++) v[i] = val;
}

void d_vector_set(double *v, int len, double val) 
{
  int i;
  for (i=0; i<len; i++) v[i] = val;
}

double d_vector_sum(double *v, int len)
{
  int i;
  double d=.0;
  for (i=0; i<len; i++) d += v[i];
  return d;
}

void d_vector_plus(double *v1, double *v2, int len)
{
  int i;
  for (i=0; i<len; i++) v1[i] += v2[i];
}

void c_vector_print(char *v, int len) 
{
  int i;
  for (i=0; i<len; i++) printf("%d ", v[i]);
}

void i_vector_print(int *v, int len) 
{
  int i;
  for (i=0; i<len; i++) printf("%d ", v[i]);
}

void d_vector_print(double *v, int len) 
{
  int i;
  for (i=0; i<len; i++) printf("%6.3lf ", v[i]);
}

void d_vector_normalize(double *v, int len)
{
  int i;
  double d=0.;

  for (i=0; i<len; i++) d += v[i];
  for (i=0; i<len; i++) v[i] /= d;
}

double d_vector_max(double *v, int len, int *pos)
{
  int i;
  double max = v[0];
  *pos = 0;
  for (i=1; i<len; i++)
    if (v[i]>max) {
      max = v[i];
      *pos = i;
    }
  return max;
}

double d_vector_max_max(double *v, double *v1, int len, int *pos)
{
  int i;
  double max = v[0], max1 = v1[0];
  *pos = 0;
  for (i=1; i<len; i++)
    if (v[i]>max || (v[i]==max && v1[i]>max1)) {
      max = v[i];
      max1 = v1[i];
      *pos = i;
    }
  return max;
}

double d_vector_min_max(double *v, double *v1, int len, int *pos)
{
  int i;
  double min = v[0], max = v1[0];
  *pos = 0;
  for (i=1; i<len; i++)
    if (v[i]<min || (v[i]==min && v1[i]>max)) {
      min = v[i];
      max = v1[i];
      *pos = i;
    }
  return min;
}

/****************************************************************************
knn PROBLEM - general solution
****************************************************************************/

double *dist = NULL;
int *nn = NULL;

void gen_knn_ini(int k, int n)
{
  dist = d_vector(k);
  nn = i_vector(k);
}

void gen_knn_term()
{
  free(dist); free(nn);
}

void gen_print_knn(int kk, Str255 s)
{
  int i;
  printf("%s ", s);
  for (i=0; i<kk; i++) printf("%d ", nn[i]);
  printf("\n");
}

void gen_find_knn(int indx, int kk, int n, char cond(int, int),
		   double distance(int,int))
{
  int i, j, k, l, last;
  double d;
				/* initialization */
  for (j=0, i=0; i<kk; j++)
    if (cond(indx, j)) {
      dist[i] = distance(indx, j);
      nn[i++] = j;
    }
  last = j;

				/* sort */
  if (kk!=1)
    for (i=0; i<kk-1; i++)
      for (j=i+1; j<kk; j++)
	if (dist[i]>dist[j]) {
	  k = dist[i]; l = nn[i];
	  dist[i] = dist[j]; nn[i] = nn[j];
	  dist[j] = k; nn[j] = l;
	}

  for (i=last; i<n; i++) 
    if (cond(indx, i)) {
      d = distance(indx, i);
      if (d<dist[kk-1]) {	/* insert */
	for (j=0; dist[j] <= d; j++);
	for (l=kk-1; l>j; l--) {
	  dist[l] = dist[l-1]; nn[l] = nn[l-1];
	}
	dist[j] = d; nn[j] = i;
      }
      else if (d==dist[kk-1]) {	/* replace */
	for (j=0; j<kk && dist[kk-1-j] == d; j++);
	k = (int)(rnd1() * (double)(j+1));
	if (k<j) {
	  dist[kk-1-k] = d; nn[kk-1-k] = i;
	}
      }
    }
}

char my_cond(int i, int j) 
{
  return i!=j;
}

double my_dist(int i, int j)
{
  return ABS(i-j);
}

/****************************************************************************
finding the leaves of the structure for a variables or a single root
variable
****************************************************************************/

list_of_vars *lv_desc;		/* this stores a list of structure's leafs
				   that influence the variable to learn for */

void find_leafs(var_type *v)
{
  int i;
  fams *f;
  list_of_vars *tmp;

  if (!v->mark) {
    v->mark = TRUE;
    if (v->famsout != NULL) {
      for (f=v->famsout; f!=NULL; f=f->next)
	for (i=0; i<f->n_in; i++)
	  find_leafs(f->in[i]);
    }
    else {
      tmp = (list_of_vars *) malloc(sizeof(*tmp));
      tmp->var = v;
      tmp->prev = lv_desc;
      lv_desc = tmp;
/*      printf("fl %s\n", v->name);  */
    }
  }
}

void find_leafs_lv(list_of_vars *lv)
{
  for (; lv!= NULL; lv=lv->next)
    find_leafs(lv->var);
}

list_of_vars *find_leaves(var_type *v)
{
  list_of_vars *l1, *l2;

  mark_var(v, FALSE);
  lv_desc = NULL;
  find_leafs(v);
  l2 = NULL;			/* reverse the list */
  for (l1=lv_desc; l1!=NULL; l1=l1->prev) {
    l1->next = l2;
    l2 = l1;
  }
  lv_desc = l2;
  return lv_desc;
}

/****************************************************************************
RULE INDEX AND ATTRIBUTE VALUES MANIPULATION
****************************************************************************/

void indx2att(fams *f, int indx, char *att)
{ 
  int n, j, nin;

  nin = f->n_in;
  for (n=1, j=nin-1; j >= 0; n*=f->in[j]->ndesc, j--)
    att[j] = (indx / n) % f->in[j]->ndesc;
}

void indx2satt(fams *f, int indx, char *att, char *sel)
{ 
  int n, j, nin;

  nin =   f->n_in;
  for (n=1, j=nin-1; j >= 0; j--)
    if (sel[j]) {
      att[j] = (indx / n) % f->in[j]->ndesc;
      n*=f->in[j]->ndesc;
    }
}

int att2indx(fams *f, char *att)
{
  int j, m, nin;
  int indx = 0;

  nin =   f->n_in;
  for (j = nin-1, m = 1; j >= 0; m *= f->in[j]->ndesc, j--)
    indx += m * att[j];
  return indx;
}

int satt2indx(fams *f, char *att, char *sel)
{
  int j, m, nin;
  int indx = 0;

/*  printf("yyy "); for (j=0; j<f->n_in; j++) printf("%d ", att[j]); */
  nin = f->n_in;
  for (j = nin-1, m = 1; j >= 0; j--)
    if (sel[j]) {
      indx += m * att[j];
      m *= f->in[j]->ndesc;
    }
/*  printf("->%d\n", indx); */
  return indx;
}

/****************************************************************************
HASH TABLE FOR LIST OF RULES
used when lists of rules are used and when they are frequently accessed
(e.g. for decomposition)
****************************************************************************/

int hash_size = 2048;		/* size of hash table */
int hash_limit = 1000 * 2048; /* limit used not to get int overflow */

void free_hash_table(fams *f)
{
  int i;
  hash_entry *he, *ohe;

  if (f->hash_table==NULL) return;
  for (i=0; i<hash_size; i++)
    for (ohe=he=f->hash_table[i]; ohe!=NULL;) {
      ohe = ohe->next;
      free(he);
      he = ohe;
    }
  FREE(f->hash_table);
}

int hash(fams *f, char *att)
{
  int j, m, nin;
  int indx = 0;

  nin =   f->n_in;
/*  for (j = nin-1, m = 1; j >= 0 && indx < hash_limit; */
  for (j = nin-1, m = 1; j >= 0 && m < hash_limit;
       m *= f->in[j]->ndesc, j--)
    indx += m * att[j];
  return (indx % hash_size);
}

void build_hash_table(fams *f)
{
  int i, h;
  rule_list *rl;
  hash_entry *he;

/*  int kk=0;  */

  free_hash_table(f);
  f->hash_table = (hash_entry **) malloc(sizeof(*(f->hash_table)) * hash_size);
  if (f->hash_table==NULL) printf("memory overrun for ht\n");
  for (i=0; i<hash_size; i++) f->hash_table[i] = NULL;

  for (rl=f->lrule; rl!=NULL; rl=rl->next) {
    rl->unused = TRUE;
    h = hash(f, rl->att);
    he = (hash_entry *) malloc(sizeof(*he));
    if (he==NULL) printf("memory overrun for ht 1\n");
    he->rl = rl;
    he->next = f->hash_table[h];
    f->hash_table[h] = he;
/*    printf("%3d ", h); kk++; if (kk==18) {printf("\n"); kk=0;}  */
  }
/*  printf("\n");  */
}

extern char bbb;
rule_list *find_hrule(fams *f, char *att)
{
  int i, h;
  hash_entry *ha;
  char b;


  h = hash(f, att);
  for (ha=f->hash_table[h]; ha!=NULL; ha=ha->next) {
/*    printf("%d(%d)\n", h, f->n_rules);  */
    if (ha->rl == NULL) printf("error in ht\n");
    if (ha->rl->unused) {	/* this is because of decomposition */
      for (b=TRUE, i=0; i<f->n_in && b; i++) {
	b = att[i] == ha->rl->att[i];
      }
      if (b) return ha->rl;
    }
  }
  return NULL;
}

int get_hrule(fams *f, char *att)
{
  rule_list *rl;

  if (f->hash_table==NULL) {printf("hash error 1\n"); exit(0);}
  rl = find_hrule(f, att);
  if (rl==NULL) return CUNDEF; else return rl->class;
}

void test_hash(Str255 name)
{
  fams *f;
  char *att;
  int i;
  rule_list *rl;

  f = find_fam(name, name);
  if (f==NULL) {printf("not found\n"); return; }

  build_hash_table(f);
  att = c_vector(f->n_in);
  for (i=0; i<f->n_rules; i++) {
    indx2att(f, i, att);
    rl = find_hrule(f, att);
    if (rl!=NULL)
      printf("%2d: %3d %d\n", i, hash(f, att), rl->class);
    else
      printf("%2d: %3d ?\n", hash(f, att), i);
  }
  free_hash_table(f);
}

/****************************************************************************
COMBINATORIAL RUTINES

These rutines are useful for selection wich decomposition to perform,
i.e., which set of variables to join. From original set, a disjunct or
non-disjunct set has to be created. Sets are represented as char
vector (if TRUE, then i-th var is in the set). There are three sets: a
(set of vars to join), b (negated a or c), and c (indicating vars that
appear in a and b). For decomposition, a=is_col, b=is_row, c=is_both.

Rutines are included that give the number of all possible combination,
reset a combination vectors, and on a basis of existing valid
combination create a next one. Also, rutines are provided that create
i-th possible combination.
****************************************************************************/

int n_dcomb(int k, int n)
{
  int c1=1, c2=1, i;
  for (i=0; i<k; i++) {
    c1 *= n-i; c2 *= k-i;
  }
  return c1 / c2;
}

int n_ndcomb(int dk, int ndk, int n)
{
  return n_dcomb(dk, n) * n_dcomb(ndk, dk);
}

void neg_vector(int n, char *b, char *a)
{
  int i;
  for (i=0; i<n; i++) b[i] = ! a[i];
}

void reset_dcomb(int k, int n, char *a, char *b)
{
  int i;
  for (i=0; i<n; i++) a[i] = i<k;
  neg_vector(n, b, a); 
}

void or_vector(int n, char *b, char *a, char *c)
{
  int i;
  for (i=0; i<n; i++) b[i] = a[i] || c[i];
}

reset_ndcomb(int ndk, int n, char *a, char *c)
{
  int i, j;
  for (i=0, j=0; i<n; i++)
    if (a[i]) c[i] = j++<ndk;
    else c[i] = FALSE;
}

void reset_comb(int dk, int ndk, int n, char *a, char *b, char *c)
{
  reset_dcomb(dk, n, a, b);
  reset_ndcomb(ndk, n, a, c);
  or_vector(n, b, b, c);
}

void print_v(int n, char *b)
{
  int i;
  for (i=0; i<n; i++) printf("%c ", b[i]?'1':'0');
}

char next_dcomb(int dk, int n, char *a, char *b)
{
  int i, j;
  char bb=FALSE;

  for (i=n-2; !bb && i>=0; i--)
    if (a[i] && !a[i+1]) {
      bb = TRUE;
      a[i] = FALSE; a[i+1] = TRUE;
    }
  if (!bb) return FALSE;
  for (i+=2, j=0; i<n; i++)
    if (a[i]) {
      a[i] = FALSE; a[i-j] = TRUE;
    }
    else j++;
  neg_vector(n, b, a);
  return TRUE;
}

char check_inclusion(int n, char *a, char *b)
{
  int i;
  for (i=0; i<n; i++)
    if (b[i] && !a[i]) return FALSE;
  return TRUE;
}

char next_ndcomb(int dk, int ndk, int n, char *a, char *b, char *c)
{
  char b2, b1;

  do {
    b2 = next_dcomb(ndk, n, c, b);
    b1 = check_inclusion(n, a, c);
  } while (b2 && !b1);
  if (b2) {
    neg_vector(n, b, a); or_vector(n, b, b, c); return TRUE;
  }

  if (next_dcomb(dk, n, a, b)) {
    reset_ndcomb(ndk, n, a, c);
    or_vector(n, b, b, c);
    return TRUE;
  }
  return FALSE;
}

void ith_dcomb(int ii, int dk, int n, char *a, char *b)
{
  int i;
  reset_dcomb(dk, n, a, b);
  for (i=0; i<ii; i++) next_dcomb(dk, n, a, b);
}

void ith_ndcomb(int ii, int dk, int ndk, int n, char *a, char *b, char *c)
{
  int i;
  reset_comb(dk, ndk, n, a, b, c);
  for (i=0; i<ii; i++) next_ndcomb(dk, ndk, n, a, b, c);
}

void test_comb_rutines(int dk, int ndk, int n)
{
  int i, j, k;
  char *a, *b, *c;

  printf("This is a test for %d %d %d\n", dk, ndk, n);
  i=n_dcomb(dk, n);
  printf("n_dcomb %d\n", i);
  k = n_ndcomb(dk, ndk, n);
  printf("n_ndcomb %d\n", k);
  a = c_vector(n);
  b = c_vector(n);
  c = c_vector(n);
  reset_comb(dk, ndk, n, a, b, c);

/*  for (j=0; j<i; j++) {
    ith_dcomb(j, dk, n, a, b);
    printf("%d: ", j); print_v(n, a); printf("\n");
    next_dcomb(dk, n, a, b);
  } */

  for (j=0; j<k; j++) {
    ith_ndcomb(j, dk, ndk, n, a, b, c);
    printf("%3d: ", j); print_v(n, a); printf(" .. "); 
    print_v(n, b); printf(" .. "); 
    print_v(n, c); printf("\n");
/*    next_ndcomb(dk, ndk, n, a, b, c);   */
  } 
}

/****************************************************************************
COMBINATORIAL FUNCTIONS
****************************************************************************/

int fact(int i)
{
  int j=i-1;
  for (; j>0; --j) i*=j;
  return i;
}

double dfact(int i)
{
  int j=i-1;
  double d;

  d = (double) i;
  for (; j>0; --j) d = d * (double) j;
  return d;
}

int icomb(int n, int i)
{
  return fact(n) / (fact(i) * fact(n-i));
}

double dcomb_i(int n, int i)
{
  return dfact(n) / (dfact(i) * dfact(n-i));
}

double logfact(int i)
{
  double s=0.;

  for (; i>1; --i) 
    s += LOG2((double)i);
  return s;
}

/* dpow_i: returns a power a^b where a and b are integer and the result
   is in double precision */

double dpow_i(int a, int b)
{
  return exp((double)b * log((double)a));
}

/* dtic: returns DTIC(F), the input parameters are |c| and ||F|| */

double dtic(int c, int f)
{
  return (double)f * LOG2(c);
}

/* dticp: returns DTIC'(H), the input parameters are |c| and ||H|| */

double dcticp(int c, int h)
{
  double comp, d, d1, d2, d3;
  int i;

  if (c==1) return 1.0;
  comp = dpow_i(c, h);
  for (i=1; i<=c-1; i++) {
    d1 = dcomb_i(c,i);
    d2 = dcticp(c-i, h);
    d3 = d1 * d2;
    comp = comp - d3;
  }
  return comp;
}

double dticp(int c, int h)
{
  if (c>12) return 1e300;
  /*  printf("XXX %lf\n", dcticp(c,h)); */
  return LOG2(dcticp(c,h) / dfact(c));
}
