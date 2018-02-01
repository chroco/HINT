/****************************************************************************
decomp.c

Function decomposition as defined by Biermann.
****************************************************************************/

#include <stdio.h>
#include <math.h>
/* #include <limits.h> */
#define GL extern
#include "sds.h"

/* debugging, deb_dec:
   0 no debugging
   1 show best partition and basic info
   2 show all partitions and basic info
   3 coloring
   4 show partition matrix */

var_type *sv;			/* variable we are splitting on */
fams *sf;			/* fam that holds the data */
int nin;			/* number of input variables to that fam */
				/* change this to char */
char **dm=NULL;			/* decomposition matrix, CUNDEF is undef */
int r0, r1;			/* two variables to decompose on */
char *xis_col = NULL;		/* variable is on the column */
char *is_row=NULL, *is_col=NULL; /* variable is on the row/column */
char *is_both=NULL;		/* variable is both in row and col */
char *rv=NULL, *rvc=NULL, *rvr=NULL; /* value of i-th var for the rule */
char nondis_dec = FALSE;	/* do we have non-disjunct decomposition? */

int nboth;			/* number of vars in both sets (free,bounded)*/
int *ncomb=NULL;		/* number of combinations for bound sets */
double *bcomb=NULL;		/* best criteria for a set of combinations */

typedef struct {
  double c;			/* criteria use for split */
  char ok;			/* ok to decompose? */
  int cm;			/* column multiplicity */
  int dfc;			/* single dfc */
  int dfc_lim;			/* limit for dfc decomposable */
  double r;			/* PHI_SC, SDTIC, dfc defided by fact */
  double r_lim;			/* limit for r decomposable */
  double cr;			/* PHI_C, DTIC, functions with lower colors
				   are subtracted, uses r_lim */
  double cr_lim;
  double error;			/* estimated error when merging */
  double m;			/* best m */
  
  char *dcol;			/* coloring */
} split_info;

extern double dtic(int c, int f);
extern double dticp(int c, int h);

typedef struct {		/* INFO ABOUT BEST PARTITION */
				/* USED TO BE ABLE TO DECOMPOSE ON ITS BASIS */
  var_type *v;
  char *is_row, *is_col, *is_both;
  split_info sinfo;
} best_part_info_type;

best_part_info_type best_part_info;

split_info **spl_cr;		/* split criteria (#cols) */

char *tdcol;			/* temporary (start) index(color) for 
				   different rows/columns in case of ND */

Str255 joined_var_name;		/* global storage for new node */

char nc;			/* number of colors returned by coloring */

rule_list **rules = NULL;	/* sorted rules to learn and to test */
int nn_rules;			/* number of the rules processing */

extern double *im_one, **im_two; /* information measure for attributes */
extern double *w_one, **w_two;	/* relieff measure for attributes */
extern double *gr_one, **gr_two; /* gain ratio measure for attributes */
extern double *gini_one, **gini_two; /* gini index for attributes */
extern double *mdl_one, **mdl_two; /* MDL for attributes */

FILE *log_file;			/* log file for report */

char check_inclusion(int n, char *a, char *b);
char next_dcomb(int dk, int n, char *a, char *b);

/****************************************************************************
UTILITY FUNCTIONS
****************************************************************************/

/* copy_split_info: copy one split info to the other */

void copy_split_info(split_info f, split_info *t)
{
  memcpy(t, &f, sizeof(f));
  t->dcol = f.dcol;
}

void reset_best_part()
{
  best_part_info.v = NULL;
  best_part_info.is_row = best_part_info.is_col = best_part_info.is_both= NULL;
}

void test_best_part(split_info si)
{
  if ((best_part_info.v == NULL || 
      best_part_info.sinfo.c > si.c) &&  si.ok) {
    copy_split_info(si, &(best_part_info.sinfo));
    best_part_info.v = sf->out;
    FREE(best_part_info.is_row); FREE(best_part_info.is_col);
    FREE(best_part_info.is_both);
    c_vector_copy(is_row, &best_part_info.is_row, sf->n_in);
    c_vector_copy(is_col, &best_part_info.is_col, sf->n_in);
    c_vector_copy(is_both, &best_part_info.is_both, sf->n_in);
  }
}

/****************************************************************************
SINGLE NODE DECOMPOSITION
tries to decompose a single node on variables r0 and r1 and counts
the number of different entries in decomposition matrix
****************************************************************************/

void print_dm(int ncol, int nrow)
{
  int i, j, k, l, m;

  for (nvrow=0, i=0; i<nin; i++) if (is_row[i]) nvrow++;
  pspace(nvrow+1);
  for (i=0; i<nin; i++)
    if (is_col[i])  printf("%c", sf->in[i]->name[0]);
  printf("\n");
  
  for (m=0, k=0; k<nin; k++) if (is_col[k]) m++;
  for (k=0, i=0; i<nin; i++)
    if (is_col[i]) {
      k++;
      if (k<m) pspace(nvrow+1);
      else {
	for (l=0; l<nin; l++)
	  if (is_row[l]) printf("%c", sf->in[l]->name[0]);
	printf(" ");
      }
      for (j=0; j<ncol; j++) {
	indx2satt(sf, j, rv, is_col);
	printf("%c ", sf->in[i]->desc[rv[i]].name[0]);
      }
      printf("\n");
    }

  for (i=0; i<nrow; i++) {
    for (j=0; j<nin; j++)
      if (is_row[j]) {
	indx2satt(sf, i, rv, is_row);
	printf("%c", sf->in[j]->desc[rv[j]].name[0]);
      }
    printf(" ");

    for (j=0; j<ncol; j++) {
      indx2satt(sf, j, rv, is_col);
      indx2satt(sf, i, rv, is_row);
      printf("%c ", dm[i][j]==CUNDEF ? '-' : (dm[i][j]==CCLASH ? 
	     'X' : sv->desc[dm[i][j]].name[0])); 
/*      printf("%d ", dm[i][j]); */
    }
    printf("\n");
  }
}

/****************************************************************************
DERIVATION OF NUMBER OF DIFFERENT ROWS/COLUMNS
dm is a matrix with unknown elements (dm[i][j]==CUNDEF)
coloring is used, these are more complex rutines
****************************************************************************/

/* cound_diff_X(): cound different columns or rows, previous
   algorithms checked if columns are full or empty, but this rutine
   does not and leaves everything to coloring. This makes the rutine
   simple and easy to understand.

   There are two methods, one is complete and the other direct. First
   one generates dm table, and, in case rule tables are used, checks
   all the rules. Second one uses just the rules defined and buils im
   matrix directly. */

int ncn;			/* number of nodes to color */
char **im = NULL;		/* incompatibility matrix */
extern int color_graph(char **im, int nn);
extern char *colors;		/* colors of nodes */

int count_diff_col(int ncol, int nrow, int ncn, char **dm)
{
  int i, j, k, m;
  char b;
  static int old_ncn;

  ncn = ncol;

  free_c_matrix(im, old_ncn);
  old_ncn = ncn;
  FREE(colors);
  im = c_matrix_ini(ncn, ncn, FALSE);

  for (i=0; i<ncol; i++)
    for (j=0; j<ncol; j++)
      if (i!=j) {
	for (b=TRUE, k=0; k<nrow && b; k++)
	  b = (dm[k][i]==CUNDEF || dm[k][j]==CUNDEF) || (dm[k][i] == dm[k][j]);
	im[i][j] = im[j][i] = (!b);
      }
      else im[i][j] = FALSE;
  m = color_graph(im, ncn);
  for (i=0; i<ncol; i++) dcol[i] = colors[i];
  return m;
}

int count_diff_row(int ncol, int nrow)
{
  int i, j, k, m;
  char b;

  ncn = nrow;
  im = c_matrix_ini(ncn, ncn, FALSE);

  for (i=0; i<nrow; i++)
    for (j=0; j<nrow; j++)
      if (i!=j) {
	for (b=TRUE, k=0; k<ncol && b; k++)
	  b = (dm[i][k]==CUNDEF || dm[j][k]==CUNDEF) || (dm[i][k] == dm[j][k]);
	im[i][j] = im[j][i] = (!b);
      }
      else im[i][j] = FALSE;
  m = color_graph(im, ncn);
  for (i=0; i<nrow; i++) drow[i] = colors[i];
  free_c_matrix(im, ncn);
  FREE(colors);
  return m;
}

rule_list *find_hrule(fams *f, char *att);

int count_diff_col_direct(int ncol, int nrow)
{
  int i, j, k, m;
  char b;			/* still any columns compatible? */
  char b1;
  char b2;			/* any rules in same column found? */
  char ndb;			/* used with non-disjunct decomposition */
  rule_list *rl1, *rl2, *rl;
  char *xrow=NULL;		/* represents a row in dm matrix */
  int *n_incom=NULL;		/* how many cols is incompatible
				   with i-th, initially 0*/
  char *b_com=NULL;		/* is i-th col still compatible with any? */
  int nr_com = ncol;		/* # still compatible columns */
  char *att=NULL, *batt;

  free_c_matrix(im, ncn);
  FREE(colors);

  ncn = ncol;
  im = c_matrix_ini(ncn, ncn, FALSE);
  xrow = c_vector(ncol);
  n_incom = i_vector_ini(ncol, 0);
  b_com = c_vector_ini(ncol, TRUE);
  att = c_vector(nin);

  for (rl1=sf->lrule; rl1!=NULL; rl1=rl1->next) {
    rl1->unused=TRUE; 
    rl1->cr = satt2indx(sf, rl1->att, is_col); /* used also in split node */
  }
  for (b=TRUE, rl1=sf->lrule; rl1->next!=NULL && b; rl1=rl1->next) 
    if (rl1->unused) {
      batt = rl1->att;
				/* fill xrow vector */
      for (i=0; i<nin; i++) if (is_row[i]) att[i]=rl1->att[i];
      for (i=0; i<ncol; i++) {
	indx2satt(sf, i, att, is_col);
	ndb = TRUE;

/* 	if (nondis_dec)
	  for (j=0; j<nin && ndb; j++)
	    if (is_both[j]) ndb = att[j] == batt[j];  */
	
	if (ndb) {
	  rl = find_hrule(sf, att);
	  if (rl==NULL) xrow[i] = CUNDEF;
	  else {
	    xrow[i] = rl->class;
	    rl->unused = FALSE;
	    b2 = TRUE;
	  }
	} 
	else xrow[i] = CUNDEF;
      }
				/* check incompatibilities */
      if (b2) {
	xrow[rl1->cr] = rl1->class;
	for (i=0; i<ncol-1; i++)
	  if (b_com[i])
	    for (j=i+1; j<ncol; j++)
	      if (b_com[j] && !im[i][j]) {
		b1 = (xrow[i]!=CUNDEF) && (xrow[j]!=CUNDEF)
		  && (xrow[i]!=xrow[j]);
		if (b1) {
		  im[i][j] = im[j][i] = TRUE;
		  ++n_incom[i]; ++n_incom[j];
		  if (n_incom[i]==ncol-1) b_com[i]=FALSE;
		  if (n_incom[j]==ncol-1) b_com[j]=FALSE;
		  if (!b_com[i]) --nr_com;
		  if (!b_com[j]) --nr_com;
		  if (nr_com<=0) b = FALSE;
		}
	      }
      }
    }
  FREE(xrow); FREE(n_incom); FREE(b_com); FREE(att);
  m = color_graph(im, ncn);
  for (i=0; i<ncol; i++) dcol[i] = colors[i];
  return m;
}

/* count_diff_col_direct: Uses Demsar's idea to sort the rules
   first. See page 24 of his BSc Thesis. */

int my_rcmp(rule_list *r1, rule_list *r2)
{
  int i, j;

  for (i=0; i<sf->n_in; i++) 
    if (is_row[i]) {
      if (r1->att[i] < r2->att[i]) return -1;
      else if (r1->att[i] > r2->att[i])	return 1;
    }
  if (r1->class < r2->class) return -1;
  if (r1->class > r2->class) return 1;
  return 0;
}

int my_part(int p, int r)
{
  rule_list *x = rules[p], *y;
  int i=p-1, j=r+1;
  while (TRUE) {
    do {
      j--;
    } while (my_rcmp(rules[j],x)>0);
    do {
      i++;
    } while (my_rcmp(rules[i],x)<0);
    if (i<j) {y=rules[i]; rules[i]=rules[j]; rules[j]=y;}
    else return j;
  }  
}

void my_quick(int p, int r)
{
  int q;
  if (p<r) {
    q = my_part(p, r);
    my_quick(p, q);
    my_quick(q+1, r);
  }
}

void ppr()
{
  int i;
  printf("nn_rules %d\n", nn_rules);

  for (i=0; i<sf->n_in; i++) printf("%s", is_row[i]?"X":"_"); printf("\n");
  for (i=0; i<nn_rules; i++) {
    printf("%5d %s%s %5d", i, 
	   rules[i]->mark?"T":"F",
	   rules[i]->unused?"T":"F", rules[i]->cr);
    list_rule(sf, rules[i]->class, rules[i]->att);
  }
}

int rule_cmp_rows_class_s(rule_list **a, rule_list **b)
{
  int i, j;
  rule_list *r1 = *a, *r2 = *b;

  for (i=0; i<sf->n_in; i++) 
    if (is_row[i]) {
      if (r1->att[i] < r2->att[i]) return -1;
      else if (r1->att[i] > r2->att[i])	return 1;
    }
  if (r1->class < r2->class) return -1;
  if (r1->class > r2->class) return 1;
  return 0;
}


char use_evidence = FALSE;	/* when coloring with no noise, prefer
				   to merge the columns that have more
				   evidence for compatibility */

int count_diff_col_fast(int ncol, int nrow)
{
  int i, j, k;
  int i1, i2;			/* idices where group starts and where
				   subgroup starts */
  rule_list *rl;
  int c1, c2;			/* class of the group */
  static int old_ncol;
  char *col_used;			/* is the column used? */

  col_used = c_vector_ini(ncol, FALSE);
  for (i=0; i<nn_rules; i++) {
    rules[i]->cr = satt2indx(sf, rules[i]->att, is_col);
    col_used[rules[i]->cr] = TRUE;
  }
  
  free_c_matrix(im, old_ncol);
  old_ncol = ncol;
  FREE(colors);

  qsort(&rules[0], nn_rules, sizeof(rule_list **), rule_cmp_rows_class_s); 
  im = c_matrix_ini(ncol, ncol, FALSE);
  for (i=0; i<nn_rules-1; i++) {
    rules[i]->mark = !rule_cmp_rows(rules[i], rules[i+1]);
    rules[i]->unused = rules[i]->mark && rules[i]->class == rules[i+1]->class;
  }
  rules[nn_rules-1]->mark = rules[nn_rules-1]->unused = FALSE;

  for (i1=0, i=0; i<nn_rules; i++) {
    i2 = i;
    for (; rules[i]->unused; i++) {
    }

    if (i2>i1)
      for (j=i2; j<=i; j++) {
	c1 = rules[j]->cr;
	for (k=i1; k<i2; k++) {
	  c2 = rules[k]->cr;
	  if (c1>=ncol || c2>=ncol || c1<0 || c2<0) {
	    printf("sss %d (%d-%d) (%d %d) - %d\n", nn_rules, j, k, c1, c2, ncol);
	    ppr(); 
	    printf("Internal error in decomp.c - blame the author\n");
	    exit(0);
	  }
	  im[c1][c2] = im[c2][c1] = TRUE;
	}
      }

    if (!rules[i]->mark) i1 = i+1;
  }

  if (use_evidence) derive_evidence();

  k = color_graph(im, ncol);
  for (i=0; i<ncol; i++)
    if (col_used[i]) dcol[i] = colors[i];
    else dcol[i] = CUNDEF;
  FREE(col_used);

/*  printf("colors=%d\n", k);
  for (i=0; i<ncol; i++) printf("%d ", dcol[i]); printf("\n"); */

  return k;
}

/* count_entries: builds a decomposition matrix and counts the number
   of different entries. For now, assumes that all the combinations
   are specified (all possible rules are given).  If rules are given
   as table, than this builds first dm matrix, then im matrix, and
   then uses coloring. Else, it derives im matrix directly. */

split_info sinfo;

#define DEC_SLOW 0
#define DEC_FAST 1

char dec_speed = DEC_FAST;

void count_entries()
{
  int r;			/* rule number */
  int irow, icol;		/* indices of decomposition matrix */
  rule_type rtp;
  rule_list *rl;
  int i, ncol, m, n;
  double e, mm;

  for (ncol=1, i=0; i<sf->n_in; i++) if (is_col[i]) ncol *= sf->in[i]->ndesc;

  if (dec_crit == c_error) {
    int i;
    i = n = get_partition_error(&e, &mm);
    sinfo.cm = i;
    sinfo.m = mm;
    sinfo.error = e;
  }
  else if (dec_crit == c_m) {
    int i;
    i = n = get_partition_m(&e);
    sinfo.error = e;
    sinfo.cm = i;
    sinfo.m = m_param;
  }
  else if (use_dm) {
				/* fill dm with its entries */
				/* it is also valid for non-disjunct dec */
    if (sf->er==er_table)
      for (r=0; r < sf->n_rules; r++) {
	indx2att(sf, r, rv);
	irow = satt2indx(sf, rv, is_row);
	icol = satt2indx(sf, rv, is_col);
	m = get_nrule_e(sf, r, &rtp);
	if (rtp != undef) dm[irow][icol] = m;
      }
    else if (sf->er==er_list) {
      for (rl=sf->lrule; rl!=NULL; rl=rl->next) {
	irow = satt2indx(sf, rl->att, is_row);
	icol = satt2indx(sf, rl->att, is_col);
	dm[irow][icol] = rl->class;
      }
    }
    if (deb_dec>3) print_dm(ncol, nrow); /* print out dm */
    /* count entries in dm */
    sinfo.cm = count_diff_col(ncol, nrow, ncn, dm);
  }
  else {
    if (dec_speed == DEC_SLOW) {
      printf("SLOW\n");
      sinfo.cm = count_diff_col_direct(ncol, nrow);
    }
    else {
      sinfo.cm = count_diff_col_fast(ncol, nrow);
    }
  }
    
  sinfo.dfc = nrow * sinfo.cm + ncol;
  sinfo.dfc_lim = nrow * ncol;

  sinfo.cr =  dtic(sv->ndesc, nrow * sinfo.cm) + dticp(sinfo.cm, ncol);
  sinfo.cr_lim = dtic(sv->ndesc, nrow * ncol);

  sinfo.r = (double) nrow * sinfo.cm * LOG2(sv->ndesc) +
    (double) ncol * LOG2(sinfo.cm) - logfact(sinfo.cm);
				/* new from 20 Apr 97 */
  sinfo.r_lim = (double) nrow * ncol * LOG2(sv->ndesc);
/*  sinfo.r_lim = (double) nrow * ncol * logfact(sv->ndesc); */

  
  switch (dec_crit) {
  case c_cm:
    sinfo.c = (double) sinfo.cm;
    sinfo.ok = sinfo.r < sinfo.r_lim;
    break;
  case c_m:
    if (m_dec == m_cm) {
      sinfo.c = (double) sinfo.cm;
      sinfo.ok = sinfo.r < sinfo.r_lim;
    }
    else {
      if (m_dec == m_error) {
	sinfo.c = sinfo.error;
	sinfo.ok = (n<ncol) && (n<100);
      }
      else if (m_dec == m_cm_error) {
	sinfo.c = sinfo.cm + sinfo.error;
	sinfo.ok = sinfo.r < sinfo.r_lim;
      }

      for (ncol=1, i=0; i<nin; i++)
	if (is_col[i]) ncol *= sf->in[i]->ndesc;
    }

    break;
  case c_dfc:
    sinfo.c = (double) sinfo.dfc;
    sinfo.ok = sinfo.dfc < sinfo.dfc_lim;
    break;
  case c_r:
    sinfo.c = sinfo.r;
    sinfo.ok = sinfo.r < sinfo.r_lim;
    break;
  case c_cr:
    sinfo.c = sinfo.cr;
    sinfo.ok = sinfo.cr < sinfo.cr_lim;
    break;
  case c_error:
    sinfo.c = sinfo.error;
    sinfo.ok = sinfo.r < sinfo.r_lim;
    break;
  }
  sinfo.dcol = dcol;
  sinfo.ok = sinfo.ok && sinfo.cm < CUNDEF;
  test_best_part(sinfo);
}

/* nd_count_entries(ncolors): counts number of diff cols for non-disjunct
   decomposition, assumes corresponding disjunct decomposition was
   done already. ncolors is # colors in colors[].

   Uses dcol[ncol]
*/

char **nim = NULL;		/* incompatibility matrix for nd case*/
char **ndm = NULL;		/* decomposition matrix for nd case*/
char *ndcol;			/* colors of nodes in nim */
char nnc;			/* number of colors */

int nd_count_diff_col(int ncol, int nrow, int ncn, char **dm)
{
  int i, j, k, m;
  char b;
  static int old_ncn;

  ncn = ncol;

  free_c_matrix(nim, old_ncn);
  old_ncn = ncn;
  FREE(colors);
  nim = c_matrix_ini(ncn, ncn, FALSE);

  for (i=0; i<ncol; i++)
    for (j=0; j<ncol; j++)
      if (i!=j) {
	for (b=TRUE, k=0; k<nrow && b; k++)
	  b = (dm[k][i]==CUNDEF || dm[k][j]==CUNDEF) || (dm[k][i]==dm[k][j]);
	nim[i][j] = nim[j][i] = (!b);
      }
      else nim[i][j] = FALSE;
  m = color_graph(nim, ncn);
  printf("Pending number of colors = %d\n", m);

  FREE(ndcol); 
  ndcol = c_vector(ncn); for (i=0; i<ncn; i++) ndcol[i] = colors[i];
  return m;
}

int nd_count_entries(int ncolors)
{
  int i, j;
  static int nboth_old=0;
  char *att, *tmp;

  att = c_vector(nin);
  for (nboth=1, i=0; i<nin; i++) 
    if (is_both[i]) nboth *= sf->in[i]->ndesc;

/*  free_c_matrix(ndm, nboth_old);  */
/*    THIS IS A COMPLETE SHIT. BOBS HERE, DON'T KNOW WHY. */
  ndm = c_matrix_ini(nboth, ncolors, CUNDEF);
  nboth_old = nboth;
  for (i=0; i<ncol; i++) {
    indx2satt(sf, i, att, is_col);
    ndm[satt2indx(sf, att, is_both)][dcol[i]] = dcol[i];
  }

  if (deb_dec>3 || TRUE) {
    printf("ndm:\n");
    for (j=0; j<nboth; j++) {
      for (i=0; i<ncolors; i++)
	if (ndm[j][i]==CUNDEF) printf("- "); else printf("%d ", ndm[j][i]);
      printf("\n");
    }
  }

  i = nd_count_diff_col(ncolors, nboth, ncolors, ndm);
  free(att);
  return i;
}

prepare_partition()
{
  int i;
  static int old_nrow;

  /* BBBthis is stupid, but when removing red this would crash w/o old_* */
  if (use_dm && old_nrow) free_c_matrix(dm, nrow);
  FREE(dcol);
  for (ncol=1, nrow=1, i=0; i<nin; i++) {
    if (is_col[i]) ncol *= sf->in[i]->ndesc;
    if (is_row[i]) nrow *= sf->in[i]->ndesc;
  }
  dcol = c_vector(ncol);
  
  if (use_dm) {
    dm = c_matrix_ini(nrow, ncol, CUNDEF);
    old_nrow = nrow;
    FREE(drow); 
    drow = c_vector(nrow);
  }
}

/****************************************************************************
Node decomposition
****************************************************************************/

void print_partition(FILE *f)
{
  int i;

  for (i=0; i<nin; i++) 
    if (is_row[i]
	&& !is_both[i]) 
      fprintf(f, "%s ", strn(sf->in[i]->name,print_short));
  for (i=0; i<nin; i++) 
    if (is_both[i]) {fprintf(f, "+ "); break;}
  for (i=0; i<nin; i++) 
    if (is_both[i]) fprintf(f, "%s ", strn(sf->in[i]->name,print_short));
  fprintf(f, "// ");
  for (i=0; i<nin; i++) 
    if (is_col[i]) fprintf(f, "%s ", strn(sf->in[i]->name,print_short));
}

/* split_node: based on is_row, is_col, and dm splits a node sv by
   inserting the new node and updating the information about sv. Works
   also for non-binary trees. It removes the dependecies if remove_dep
   is set. Requires following to be set:

  sv, sf, nin
  is_col, is_row
  ncol, nrow
  rv = c_vector(nin);
  dm = c_matrix_ini(nrow, ncol, CUNDEF);
  dcol = c_vector(ncol);
  drow = c_vector(nrow);
  count_entries();
  nentries = sinfo.cm;
  
  Initially in ND (Col && Both) is non-empty and (Row && Both is
  empty). */

var_type *split_node(int nentries, Str255 newname, char remove_dep)
{
  var_type *nv;
  int i, j, k, l;
  int nvcol, nvrow, nvboth=0;	/* number of vars in columns/rows */
  Str255 s;
  list_of_fams *lf;
  list_of_vars *lv;
  fams *f, *tmpf;		/* tmpf is a copy of sf being decomposed */
  char b, *xboth, *xxboth, *xrow, *xclass;

  rule_list *rls;		/* copy of original rules if list exists */
  rule_list *rl, *rl1;
  char *att, *atto;

				/* stuff to update sv, sf */
  var_type **in;
  int n_rules, n_in;		/* don't mix with nin */

  tmpf = (fams *) malloc(sizeof(*tmpf));
  tmpf->hash_table = sf->hash_table;
  sf->hash_table = NULL;
  tmpf->n_in = sf->n_in;
  tmpf->in = (var_type **) malloc(sizeof(*tmpf->in) * tmpf->n_in);
  for (i=0; i<sf->n_in; i++) tmpf->in[i] = sf->in[i];

				/* new node/variable and its rule table */
  for (nvcol=0, i=0; i<nin; i++) if (is_col[i]) nvcol++;
  for (nvrow=0, i=0; i<nin; i++) if (is_row[i]) nvrow++;
  nboth = 1;
  if (nondis_dec) {
    for (nvboth=0, i=0; i<nin; i++) if (is_both[i]) nvboth++;
    for (nboth=1, i=0; i<nin; i++) if (is_both[i]) nboth*=sf->in[i]->ndesc;
  }

/*  printf("CRB: "); 
  for (i=0; i<nin; i++)
    printf("%s-%d%d%d ", sf->in[i]->name, is_col[i], is_row[i], is_both[i]);
  printf("\n"); */

  if (!remove_dep) {
    nv = add_var(&variables, newname, TRUE, nentries)->var;

    lf = (list_of_fams *) malloc(sizeof(*lf));
    lf->next = NULL; lf->prev = NULL; lf->fam = sf;
    nv->famsin = lf;
    
    f = (fams*) malloc(sizeof(*f));
    f->hash_table = NULL;
    strcpy(f->name, nv->name);
    f->n_in = nvcol;
    f->in = (var_type **) malloc(sizeof(*f->in)*f->n_in);
    for (i=0, j=0; i<nin; i++) if (is_col[i]) f->in[j++] = sf->in[i];
    f->out = nv;
    


				/* here we set the function H */
    allocate_rules(f, ncol);

    /* here we compute the evidence of each column */
/* printf("START (ncol=%d, class=%d)\n", ncol, f->out->ndesc);
    { int i; for (i=0; i<ncol; i++) printf("%d ", dcol[i]); printf("\n"); } */
    if (use_distribution) {
      int i,j;
      double *e = d_vector_ini(ncol,0);
      extern int *ccolor;

      for (i=0; i<nn_rules; i++) 
	for (j=0; j<sf->out->ndesc; j++) {
    /* printf("xx %d", rules[i]->cr); d_vector_print(rules[i]->dist,f->out->ndesc); printf("\n"); */
	  e[rules[i]->cr] += rules[i]->dist[j];
	}

      for (i=0; i<ncol; i++) 
	if (dcol[i]!=CUNDEF) 
	  set_nrule_ndis(f, i, dcol[i], autom, e[i]);

      free(e);
/* printf("END\n"); */
    }
    else
      for (i=0; i<ncol; i++) 
	if (dcol[i]!=CUNDEF) set_nrule(f, i, dcol[i], autom);


/*    printf("New colors\n");
    pspace(nvcol); for (i=0; i<ncol; i++) 
      if (dcol[i]!=CUNDEF) printf("%d ", dcol[i]); else printf("- "); */

    f->prev = f->next = NULL;
    f->degrees = d_vector(f->out->ndesc);
    nv->famsout = f;

    /* list of input fams should now include f instead of sf */
    for (i=0; i<f->n_in; i++) {
      for (lf=f->in[i]->famsin; lf->fam != sf; lf=lf->next) {
/*	printf("check %s %s\n", lf->fam->name, lf->fam == sf?"y":"n"); */
      }
      lf->fam = f;
    }
  }
				/* update old node sv */
				/* discovered node is the last one */

  sf->n_in = n_in = nvrow + nvboth + (remove_dep ? 0 : 1);

  xrow = c_vector_ini(n_in, FALSE);
  xclass = c_vector_ini(n_in, FALSE);
  in = (var_type **) malloc(sizeof(*sf->in)*n_in);

				/* BBBB */
  for (j=0, i=0; i<nin; i++) if (is_row[i]) {
/*    xrow[j]=TRUE;  */
    in[j++]=sf->in[i];
  }
/*  xclass[j]=TRUE; */

  if (!remove_dep) in[j++] = nv;

  if (nondis_dec) {
    xboth = c_vector_ini(n_in, FALSE);
    for (i=0; i<nin; i++) if (is_both[i]) in[j++]=sf->in[i];
  }

  sf->in = in;

  n_rules = nrow * nentries * nboth;

  rls = sf->lrule;	/* copy of original rules if list exists */
  sf->lrule = NULL;

  free_and_allocate_rules(sf, n_rules);
  sf->degrees = d_vector(sf->out->ndesc);

  if (use_dm) {
    if (nondis_dec) {

      xboth = c_vector_ini(f->n_in, FALSE);
      for (i=0, j=0; i<nin; i++)
	if (is_col[i]) {
	  if (is_both[i]) xboth[j]=TRUE;
	  j++;
	}

      att = c_vector(n_in);
      for (i=0; i<nrow; i++)
	for (j=0; j<nentries; j++) {
	  for (l=0; l<nboth; l++) {
	    for (b=FALSE, k=0; k<ncol && !b; k++) {
	      indx2att(f, k, att);
	      if (satt2indx(f, att, xboth) == l)
		b = (dcol[k]==j) && (dm[i][k]!=CUNDEF);
	    }
	    if (b) {
	      k--;
	      set_nrule(sf, i*nentries*nboth + j*nboth + l, dm[i][k], autom);
	    }
	  }
	}
    }
    else {
      for (i=0; i<nrow; i++)
	for (j=0; j<nentries; j++) {
	  for (b=FALSE, k=0; k<ncol && !b; k++)
	    b = (dcol[k]==j) && (dm[i][k]!=CUNDEF);
	  if (b) {
	    k--;
	    set_nrule(sf, i*nentries + j, dm[i][k], autom);
	  }
	}
    }
  }

  else {
    if (dec_crit == c_error || dec_crit == c_m) {
      derive_G(sf, f, tmpf, remove_dep);
      free_lrule(rls);		/* BBB here most of the time is spent */
    }
    else if (dec_speed == DEC_SLOW) {
    /* for each rule, find which row it belongs, which are the values of 
       ND vars, use it for new rule,
       and mark all rules that belong different column and have same ND vars
       as used
       (unused=FALSE) BBBBB */
      att = c_vector(n_in);

      atto = c_vector(nin);
      for (rl=rls; rl!=NULL; rl=rl->next) rl->unused=TRUE; 

      k=0;
      for (rl=rls; rl!=NULL; rl=rl->next) 
	if (rl->unused) {
	  for (j=0, i=0; i<nin; i++)
	    if (is_row[i]) {
	      att[j++] = rl->att[i];
	      atto[i] = rl->att[i];
	    }
	  if (!remove_dep) att[j++] = dcol[rl->cr];
	  
	  if (nondis_dec)
	    for (i=0; i<nin; i++)
	      if (is_both[i]) att[j++] = rl->att[i];
	
	  set_rule(sf, att, rl->class, rl->rtp);
	  
	  k++;
	  for (i=0; i<ncol; i++) {
	    
	    if (nondis_dec) {
	      b = satt2indx(tmpf, rl->att, is_both) == i;
	    }
	    else b = TRUE;
	    if (i != rl->cr && dcol[rl->cr]==dcol[i] && b) {
	      indx2satt(tmpf, i, atto, is_col);
	      rl1 = find_hrule(tmpf, atto);
	      if (rl1!=NULL) rl1->unused = FALSE;
	    }
	  }
	  rl->unused = FALSE;
	}
      free_hash_table(tmpf);
      
      if (deb_dec) printf("%d created\n", k);
      /*    list_rule_table(sf, 2); */
      FREE(att); FREE(atto);
      free_lrule(rls);		/* BBB here most of the time is spent */
    }
    else if (dec_speed == DEC_FAST) {
				/* here the rules are sorted by row */
				/* attributes and class */
      int i, j, k, l;
      double ndis;
      char *c_free = c_vector(nentries);
      att = c_vector(n_in);

/*      for (i=0; i<nn_rules; i++) {
	printf("%s ", rules[i]->mark?"T":"F");
	for (j=0; j<tmpf->n_in; j++) printf("%2d ", rules[i]->att[j]);
	printf("\n");
      } 
      exit(0); */

      for (i=0; i<nn_rules; i++) {
	c_vector_set(c_free, nentries, TRUE);
	for (; i<nn_rules; i++) {
	  if (c_free[j = dcol[rules[i]->cr]]) {
	    c_free[j] = FALSE;
	    for (l=0, k=0; k<nin; k++)
	      if (is_row[k])
		att[l++] = rules[i]->att[k];
	    if (!remove_dep) att[l++] = j;
	    if (use_distribution) {
	      for (k=i, ndis=0.; TRUE; k++) {
		/* BBB this should work since they are sorted by class */
/*		if (dcol[rules[i]->cr] != dcol[rules[k]->cr]) break; */
		if (dcol[rules[i]->cr] == dcol[rules[k]->cr])
		  for (l=0; l<sf->out->ndesc; l++)
		    ndis += rules[k]->dist[l];
		if (!rules[k]->mark) break;
	      }
	      set_rule_ndis(sf, att, rules[i]->class, rules[i]->rtp, ndis);
	    }
	    else
	      set_rule(sf, att, rules[i]->class, rules[i]->rtp);
	  }
/*	  if (dcol[rules[i]->cr]==CUNDEF) {
	    printf("error in fast decompose: undefined\n");
	    exit(0);
	  } */
	  if (!rules[i]->mark) break; /* new row */
	}
      }
      FREE(att); FREE(c_free);
      free_lrule(rls);		/* BBB here most of the time is spent */
    }
  }

  for (i=0; i<nin; i++) if (is_both[i]) {
    lf = (list_of_fams *) malloc(sizeof(*lf));
    lf->next = tmpf->in[i]->famsin;
    lf->fam = sf;
    tmpf->in[i]->famsin = lf;
  }

  if (!remove_dep) set_class_distr(nv);
  set_class_distr(sf->out);
  return nv;
}

int orig_n_dis_hi;

void init_decompose_node()
{
  int i, j, bn, k;
  int dk, ndk;

  orig_n_dis_hi = n_dis_hi;
  if (nin <= n_dis_lo) return;
  n_dis_hi = MMIN(n_dis_hi, nin-1); /* bound set cant be same size the split */
  
  if (deb_dec) {
    printf("Decompose: %s = %s(", sv->name, sv->name);
    for (i=0; i<nin; i++) printf("%s ", sf->in[i]->name); 
    for (i=0, j=1; i<sf->n_in; i++) j*=sf->in[i]->ndesc;
    printf(")\ndfc=%d sdtic=%lf ", j, (double) j * LOG2(sv->ndesc));
    printf(" dtic=%lf\n", dtic(sv->ndesc,j));
    fflush(stdout);
  }
  bn = (n_dis_hi - n_dis_lo + 1) * (n_ndis_hi - n_ndis_lo + 2);
  spl_cr = (split_info **) malloc((unsigned) bn * sizeof(*spl_cr));
  ncomb = i_vector(bn);
  bcomb = d_vector(bn);

  for (bn=0, dk=n_dis_lo; dk<=n_dis_hi; dk++) {
    k = n_dcomb(dk, nin);
    ncomb[bn] = k;
    spl_cr[bn++] = (split_info *) malloc(k * sizeof(**spl_cr));

    for (ndk=n_ndis_lo; ndk<=n_ndis_hi; ndk++) {
      k = n_ndcomb(dk, ndk, nin);
      ncomb[bn] = k;
      spl_cr[bn] = (split_info *) malloc(k * sizeof(**spl_cr));
      bn++;
    }
  }

  is_col = c_vector(nin); is_row = c_vector(nin); is_both = c_vector(nin);
  rv = c_vector(nin);
  if (!use_dm) {
    if (dec_speed == DEC_SLOW) build_hash_table(sf);
    else {
      rule_list *rl;
      nn_rules = count_rules(sf);
      rules = (rule_list **) malloc(nn_rules * sizeof(*rules));
      for (i=0, rl=sf->lrule; rl!=NULL; i++, rl=rl->next)
	rules[i] = rl;
    }
  }

  for (nvcol=0, i=0; i<nin; i++) if (is_col[i]) nvcol++;
  for (nvrow=0, i=0; i<nin; i++) if (is_row[i]) nvrow++;
}

void term_decompose_node()
{
  int i, k;

  if (!use_dm) {
/*    if (dec_speed == DEC_SLOW) free_hash_table(); */
    if (dec_speed == DEC_FAST) FREE(rules);
  }
  FREE(rv); FREE(is_col); FREE(is_row); FREE(is_both);
  FREE(ncomb); FREE(bcomb);
  k = (n_dis_hi-n_dis_lo+1) * (n_ndis_hi-n_ndis_lo+1);
/*  if (spl_cr!=NULL) {
    for (i=0; i<k; i++) FREE(spl_cr[i]);
    FREE(spl_cr);
  } */
  /* BBB above should be after the results are reported */
  n_dis_hi = orig_n_dis_hi;
}

void log_best_partitions(int bn, int dk, int ndk)
{
  int i;

  bcomb[bn] = spl_cr[bn][0].c;
  for (i=1; i<ncomb[bn]; i++) {
    if (bcomb[bn] > spl_cr[bn][i].c) bcomb[bn] = spl_cr[bn][i].c;
    next_dcomb(dk, nin, is_col, is_row);
  }
  fprintf(lfile, "DECH best for (%d)(%d) with %d diff cols\n",
	  dk, ndk, (int) bcomb[bn]);
  for (i=0; i<ncomb[bn]; i++)
    if (spl_cr[bn][i].c == bcomb[bn]) {
      ith_ndcomb(i, dk, ndk, nin, is_col, is_row, is_both);
      fprintf(lfile, "DECB %d: ", (int) spl_cr[bn][i].c);
      print_partition(lfile);
      fprintf(lfile, " (%dx%d) (%d,%d)\n", nrow, ncol, bn, i);
      fflush(lfile);
    }
}

/* select_best_partition: finds the best partition under a criteria in
   spl_cr, and returns number of values of the new variables */

double select_best_partition(char dec_test)
{
  int i, bn, dk, ndk, k;
  int best_i, best_dk, best_ndk, best_bn;
  double best_bcomb, tmp;

  /* report best partition of each group */
  for (bn=0, dk=n_dis_lo; dk<=n_dis_hi; dk++) {
    log_best_partitions(bn, dk, 0);
    bn++;
    for (ndk=n_ndis_lo; ndk<=n_ndis_hi; ndk++, bn++) 
      if (ndk<dk && ndk!=0) {
	log_best_partitions(bn, dk, ndk);
      }
  }

  /* select best partition */
  /* prefer those with smaller nondisjunct set and 
     those with bigger bound set */

  if (!dec_test) {

    best_bcomb = 1e300;	/* should be DBL_MAX really */

				/* check for disjunctive */
    best_ndk = 0; 
    for (bn=0, dk=n_dis_lo; dk<=n_dis_hi; dk++) {
      if (best_bcomb >= bcomb[bn]) {
	best_bcomb=bcomb[bn]; best_bn=bn; best_dk = dk;
      }
      bn += 2 + n_ndis_hi - n_ndis_lo;
    }
    tmp = best_bcomb;

				/* check for non-disjunctitve */
    for (bn=0, dk=n_dis_lo; dk<=n_dis_hi; dk++) {
      bn++;
      for (ndk=n_ndis_lo; ndk<=n_ndis_hi; ndk++, bn++) 
	if (ndk<dk && ndk!=0)
	  if (tmp > bcomb[bn] && best_bcomb >= bcomb[bn]) {
	    best_bcomb=bcomb[bn]; best_bn=bn; best_dk = dk; best_ndk = ndk; 
	  }
    }

				/* selecting the best one from a group */
    best_i = -1;
    for (i=0; i<ncomb[best_bn]; i++)
      if (spl_cr[best_bn][i].ok && spl_cr[best_bn][i].c == best_bcomb) {
	best_i = i;
      }
    if (best_i == -1) return -1;

    ith_ndcomb(best_i, best_dk, best_ndk, nin, is_col, is_row, is_both);
    prepare_partition();

    fprintf(lfile, "DECS %d: ", (int) spl_cr[best_bn][best_i].c);
    print_partition(lfile); fprintf(lfile, " (%dx%d)\n", nrow, ncol);
    fflush(lfile);
    if (deb_dec) {
      printf("Best: "); print_partition(stdout);
      printf(" == %9.4lf\n", spl_cr[best_bn][best_i].c);
      fflush(stdout);
    }

    count_entries();
    return sinfo.cm;
  }
}

void log_decompose_node()
{
  int bn, dk, ndk, comb;

  for (bn=0, dk=n_dis_lo; dk<=n_dis_hi; dk++) {
/*   for (dkbn=0; bn < (n_dis_hi - n_dis_lo + 2) * (n_ndis_hi - n_ndis_lo + 1);) { */
    reset_comb(dk, 0, nin, is_col, is_row, is_both);
    for (comb=0; comb<ncomb[bn]; comb++) {
      fprintf(lfile, "DECP %d: ", (int)spl_cr[bn][comb].c); 
      print_partition(lfile);
      fprintf(lfile, " (%dx%d) (%d,%d)\n", nrow, ncol, bn, comb);
      fflush(lfile);
      next_dcomb(dk, nin, is_col, is_row);
    }
    bn++;
    for (ndk=n_ndis_lo; ndk<=n_ndis_hi; ndk++, bn++) 
      if (ndk<dk && ndk!=0) {
	reset_comb(dk, ndk, nin, is_col, is_row, is_both);
	for (comb=0; comb<ncomb[bn]; comb++) {
	  fprintf(lfile, "DECP %d: ", (int)spl_cr[bn][comb].c);
	  print_partition(lfile);
	  fprintf(lfile, " (%dx%d) (%d,%d)\n", nrow, ncol, bn, comb);
	  fflush(lfile);
	  next_ndcomb(dk, ndk, nin, is_col, is_row, is_both);  
	}
      }
  }
}

char decompose_node(char dec_test)
{
  char b1, b2;
  int i, j, a, b;
  int dk, ndk;			/* size of bound set and #non-disjoint vars */
  int bn=0;			/* size of bound set and id num */
  int ndbn, ndcomb;		/* bn where ndk is 0 */
  int comb;			/* combination number and its max */
  char *tmp = c_vector(nin);
  char use_ev_tmp = use_evidence;

  use_evidence = FALSE;
  init_decompose_node();

  fprintf(lfile, "DECI one decompose %s = %s(", sv->name, sv->name);
  for (i=0; i<nin; i++) fprintf(lfile, "%s ", sf->in[i]->name);
  fprintf(lfile, ")\n");
  fprintf(lfile, "DECI disjunct (%d,%d), non-disjunct (%d,%d)\n",
	  n_dis_lo, n_dis_hi, n_ndis_lo, n_ndis_hi);


  for (bn=0, dk=n_dis_lo; dk<=n_dis_hi; dk++) {
    a = n_dcomb(dk, nin);
/*    printf("n_dcomb %d\n", a); */

    reset_comb(dk, 0, nin, is_col, is_row, is_both);
    for (comb = 0; comb < ncomb[bn]; comb++) {
      for (i=0; i<nin; i++) is_both[i]=0;
      prepare_partition();
      count_entries();
      copy_split_info(sinfo, &(spl_cr[bn][comb]));

      if (deb_dec>1) {

	printf("Partx: (%2d,%2d) ", bn, comb); print_partition(stdout);
	printf(": %dx%d= %d, ", nrow, ncol, spl_cr[bn][comb].cm);

	printf("%d, ", spl_cr[bn][comb].dfc);
	printf("%6.5g, ", spl_cr[bn][comb].r);
	printf("%9.8g, ", spl_cr[bn][comb].cr);
	printf("%8.6lf(%4.1lf), ", spl_cr[bn][comb].error, 
	       spl_cr[bn][comb].m);

	printf("%s\n", spl_cr[bn][comb].ok?"ok":"not");
	fflush(stdout);
      }

      /* here we now test for non-disjoint decompositions */
      for (ndbn=bn+1, ndk=n_ndis_lo; ndk<=n_ndis_hi; ndbn++, ndk++)
	if (ndk<dk && ndk!=0) {
	  b = n_ndcomb(dk, ndk, nin) / a;
	  for (i=0, j=0; i<nin; i++) {
	    is_both[i] = is_col[i] && j<ndk;
	    if (is_col[i]) j++;
	  }
	  for (ndcomb = comb*b; ndcomb < comb*b+b; ndcomb++) {
	    spl_cr[ndbn][ndcomb].c = nd_count_entries(spl_cr[bn][comb].cm);

	    printf("Part: (%2d,%2d) ", ndbn, ndcomb); print_partition(stdout);
	    printf(": Comb=%dx%d == %d\n", nrow, ncol, spl_cr[ndbn][ndcomb].cm);
	    fflush(stdout);

	    do {
	      b2 = next_dcomb(ndk, nin, is_both, tmp);
	      b1 = check_inclusion(nin, is_col, is_both);
	    } while (b2 && !b1);
	  }
	}
      next_dcomb(dk, nin, is_col, is_row);
    }

    bn += n_ndis_hi - n_ndis_lo + 2;
  }

  log_decompose_node();

				/* generate and insert new node */
  use_evidence = use_ev_tmp;
  if (!dec_test) {
    i = select_best_partition(dec_test);
    if (i>0) {split_node(i, gen_var_name(), i<=1);}
    else return FALSE;
  }
  return TRUE;
  /* BBB this should be here but causes a bomb when global
     decomposition because init_dec_node was not called before */

/*  term_decompose_node(); */
}

/****************************************************************************
REPORT
Report the results of split and of heuristic approximation to the selection
of attributes on which to split on.
****************************************************************************/

void report_split_results_head()
{
  int i;

  if (log_file==NULL) {printf("Failed to open dec_log\n"); exit(0);}
  fprintf(log_file, "Attribute selection measures for:\n    ");
  for (i=0; i<nin; i++) fprintf(log_file, "%s ", strn(sf->in[i]->name,7));
  fprintf(log_file, "\n");
}

void report_done(double *one, Str255 s)
{
  int i;
  fprintf(log_file, "%s\n", s);
  for (i=0; i<nin; i++) fprintf(log_file, "%7.3lf ", one[i]);
  fprintf(log_file, "\n");
}

void report_dtwo(double **two)
{
  int i, j;
  for (i=0; i<nin-1; i++) {
    for (j=0; j<nin; j++)
      if (j>i) fprintf(log_file, "%7.3lf ", two[i][j]);
      else fprintf(log_file, "        ");
    fprintf(log_file, "\n");
  }
}

void report_donetwo(double *one, double **two, Str255 s, Str255 code)
{
  int i, j;
  double d;

  fprintf(log_file, "%s\n", s);
  for (i=0; i<nin-1; i++) {
    for (j=0; j<nin; j++)
      if (j>i) {
	two[i][j] -= one[i]+one[j]; /* distortion!!! */
	fprintf(log_file,"%7.3lf ", two[i][j]);
      }
      else fprintf(log_file, "        ");
    fprintf(log_file, "\n");
  }
/*  d = bin_compare_pos_heur_real(spl_cr, two, nin);
  fprintf(log_file, "lower %6.2lf%%\n", d);
  fprintf(lfile, "HP%s %6.2lf%% better to split\n", code, d);
  d = bin_compare_splits_heur_real(spl_cr, two, nin);
  fprintf(log_file, "splits %6.2lf%%\n", d);
  fprintf(lfile, "HS%s %6.2lf%% splits to find optimal (%6.2lf%% min)\n",
	  code, d, 1./(double) (nin * (nin-1.) / 2.) * 100.); */
}

void report_split_results()
{
  int i, j;
  double d;
  
 return;				/* BBBBB */
  if (log_inform) {
    report_done(im_one, "Informativity measure");
    report_dtwo(im_two);
    report_donetwo(im_one, im_two, "IM(x,y) - IM(x) - IM(y)", "IM");
  }
  if (log_gr) {
    report_done(gr_one, "Gain ratio");
    report_dtwo(gr_two);
    report_donetwo(gr_one, gr_two, "GR(x,y) - GR(x) - GR(y)", "GR");
  }
  if (log_gini) {
    report_done(gini_one, "GINI index");
    report_dtwo(gini_two);
    report_donetwo(gini_one, gini_two, "GI(x,y) - GI(x) - GI(y)", "GI");
  }
  if (log_relieff) {
    report_done(w_one, "Relieff");
    report_dtwo(w_two);
    report_donetwo(w_one, w_two, "R(x,y) - R(x) - R(y)", "RF");
  }

  fprintf(log_file, "Decomposition criteria C(x,y):\n");
 /* for (i=0; i<nin-1; i++) {
    for (j=0; j<nin; j++) {
      printf("%d %d\n", i, j);
      if (j>i) fprintf(log_file, "  %4d  ", (int)spl_cr[i][j].c);
      else fprintf(log_file, "        ");
    }
    fprintf(log_file, "\n"); 
  } */
}

/****************************************************************************
DECOMPOSITION (PARTITION) TABLE
This part is ment for manual exploration of what we get for certain 
decomposition. Input is a node in a decision model, and a split
defined by two subsets, which are not neccessary disjunctive.
****************************************************************************/

void decomposition_table(Str255 vname, list_of_vars *lv1, list_of_vars *lv2)
{
  int i, j, k, m, indx;
  char b;			/* clash of attribute values, such dm field 
				   is then indicated with CCLASH */
  list_of_vars *lv;
  rule_type rtp;
				/* does the node exist? */

  FIND_VAR(sv, vname);
  sf = sv->famsout;
  nin = sf->n_in;
				/* preset the  */
  for (nvrow=0, lv=lv1; lv!=NULL; lv=lv->next, nvrow++);
  for (nvcol=0, lv=lv2; lv!=NULL; lv=lv->next, nvcol++);
  is_col = c_vector_ini(nin,FALSE);
  is_row = c_vector_ini(nin,FALSE);
  is_both = c_vector(nin);
  for (i=0, lv=lv1; lv!=NULL; lv=lv->next) {
    for (j=0; j<nin && sf->in[j]!=lv->var; j++);
    is_row[j]=TRUE;
  }
  for (i=0, lv=lv2; lv!=NULL; lv=lv->next) {
    for (j=0; j<nin && sf->in[j]!=lv->var; j++);
    is_col[j]=TRUE;
  }
  for (i=0; i<nin; i++) is_both[i] = is_col[i] && is_row[i];

  for (ncol=1, nrow=1, i=0; i<nin; i++) {
    if (is_col[i]) ncol *= sf->in[i]->ndesc;
    if (is_row[i]) nrow *= sf->in[i]->ndesc;
  }
  dm = c_matrix_ini(nrow, ncol, CUNDEF);
  dcol = c_vector(ncol);
  drow = c_vector(nrow);
  rv = c_vector(nin); rvr = c_vector(nin); rvc = c_vector(nin);

  printf("Decompose: %s = %s(", sv->name, sv->name);
  for (i=0; i<nin; i++) printf("%s ", sf->in[i]->name);
  printf(")\n");
  printf("Part: "); print_partition(stdout);
  printf("\n");

  for (i=0; i<nrow; i++)
    for (j=0; j<ncol; j++) {
      for (k=0; k<nin; k++) rvr[k]=CUNDEF;
      for (k=0; k<nin; k++) rvc[k]=CUNDEF;
      indx2satt(sf, i, rvr, is_row);
      indx2satt(sf, j, rvc, is_col);

      for (b=TRUE, k=0; k<nin && b; k++) if (is_both[k]) b = rvr[k] == rvc[k];
      for (k=0; k<nin; k++) if (rvr[k]==CUNDEF) rvr[k] = rvc[k];
      indx = att2indx(sf, rvr);
      m = get_nrule_e(sf, indx, &rtp);
      if (b) {
	if (rtp != undef) 
	  dm[i][j] = m;
      }
      else {
	for (k=0; k<nin; k++) if (rvc[k]==CUNDEF) rvc[k] = rvr[k];
	dm[i][j] = m;
	indx = att2indx(sf, rvc);
	m = get_nrule_e(sf, indx, &rtp);
	if (dm[i][j] != m)
	  dm[i][j] = CCLASH;
      }
    }
  
  if (sf->er==er_table) {
    print_dm(ncol, nrow);
    printf("diff col %d\n", count_diff_col(ncol, nrow, ncn, dm));
  }
  else  if (sf->er==er_list) {
    printf("here diff dcol\n");
    printf("diff dcol %d\n", count_diff_col_direct(ncol, nrow));
  }
/*  printf("diff row %d\n", count_diff_row(ncol, nrow)); */
  
				/* count diff rows and cols */
				/* implement this for empty matrix using
				   coloring */
				/* print estimates (heuristics) */

  FREE(rv); FREE(rvr); FREE(rvc); FREE(is_col); FREE(is_row);
  FREE(dcol); FREE(drow);
  /*  free_c_matrix(dm, nrow); */
}

/****************************************************************************
BOTTOM UP DECOMPOSITION -- MAIN FUNCTION
****************************************************************************/

void rule_table_statistics()
{
  int i, n=0;
  n = count_rules(sf);
  fprintf(lfile, "DRUL %d rules of %d possible (%6.2lf%%)\n",
	  n, sf->n_rules, (double) n / sf->n_rules * 100.);
}

/* bottom_up_decompose: dec_multi? if true, then decompose until
   nothing is left to decompose, i.e., until we get a binary tree. */

void bottom_up_decompose(Str255 vname, char dec_multi, char dec_test)
{
  char do_decomp = TRUE, c;

  FIND_VAR(sv, vname);
  log_file = fopen("log_dec", "w");
  sf = sv->famsout;
  while (sv->famsout->n_in > 2 && do_decomp) {
    nin = sf->n_in;
    if (nin <= n_dis_hi) break;	/* this can be commented */
    report_split_results_head();
    rule_table_statistics();
    if (log_inform || log_gr || log_gini) {
/*      printf("IM & GR & GINI one ...\n"); */
/*      informativity_one(sf); */
/*      printf("IM & GR & GINI two ...\n"); */
/*      informativity_two(sf); */
    }
    if (log_relieff) {
/*      printf("relieff ...\n");
      relieff(sf); */
    }
/*    printf("decomposition ... \n"); */

    c = decompose_node(dec_test);
    report_split_results();
    if (log_inform) {FREE(im_one); FREE(im_two);}
    if (log_gr) {FREE(gr_one); FREE(gr_two);}
    if (log_gini) {FREE(gini_one); FREE(gini_two);}
    if (log_relieff) {FREE(w_one); FREE(w_two);}
/*     free_i_matrix(spl_cr, nin); BBB ? why bombs?? */
    do_decomp = dec_multi && !dec_test && c;
  }
  fclose(log_file);
}


void global_decompose(Str255 vname, char dec_multi, char dec_test)
{
  char do_decomp = TRUE, do_loop = TRUE;
  list_of_vars *lv, *lvt, *l2;
  var_type *v;
  char use_ev_tmp = use_evidence;

  printf("global decompose ... \n");

  FIND_VAR(sv, vname);
				/* QWERT */
  lv = (list_of_vars *) malloc(sizeof(*lv));
  lv->var = sv;
  lv->next = NULL;

  for (lvt=lv; lvt!=NULL; lvt=lvt->next) lvt->var->decomposable = TRUE;

  do {
    use_evidence = FALSE;
    reset_best_part();
    printf("--------------- ");
    for (lvt=lv; lvt!=NULL; lvt=lvt->next)
      printf("%s ", lvt->var->name);
    printf("\n");

    do_loop = do_decomp = FALSE;
    for (lvt=lv; lvt!=NULL; lvt=lvt->next) 
      if (lvt->var->decomposable) {
	do_loop = TRUE;
	sv = lvt->var;
	sf = sv->famsout;
	nin = sf->n_in;

	if (nin > n_dis_lo) {
	  printf("EXAMINE %s %d\n", sv->name, sf->n_in);
	  decompose_node(TRUE);
	  do_decomp = (best_part_info.v != NULL) && (best_part_info.sinfo.ok);
	}
	if (do_decomp) break;
	else sv->decomposable = FALSE;
      }
    use_evidence = use_ev_tmp;

    if (do_decomp) {

    /* after this, best_part_info holds the data about the best
       partition to use for decomposition */

      sv = best_part_info.v;
      sf = sv->famsout;
      nin = sf->n_in;

      c_vector_copy(best_part_info.is_row, &is_row, nin);
      c_vector_copy(best_part_info.is_col, &is_col, nin);
      c_vector_copy(best_part_info.is_both, &is_both, nin);
    
      printf("Best: %s %d, %d (", 
	     best_part_info.v->name, (int) best_part_info.sinfo.cm,
	     (int) best_part_info.sinfo.c);
      print_partition(stdout);
      printf(")\n");


      prepare_partition();
      if (!use_dm) {
	if (dec_speed==DEC_SLOW) build_hash_table(sf);
	else {
	  int i;
	  rule_list *rl;
	  nn_rules = count_rules(sf);
	  rules = (rule_list **) malloc(nn_rules * sizeof(*rules));
	  for (i=0, rl=sf->lrule; rl!=NULL; i++, rl=rl->next)
	    rules[i] = rl;
	}
      }

      count_entries();

      /* note that if use_evidence, int the loop the distribution was
         not used, and outside it it was, so the number of colors may
         be different (and actually is in some cases) */

      copy_split_info(sinfo, &(best_part_info.sinfo));
/*      best_part_info.sinfo.cm = sinfo.cm; */

      split_node(best_part_info.sinfo.cm, gen_var_name(), best_part_info.sinfo.cm<=1);

				/* add the new variable to the list */
      v = sv->famsout->in[sv->famsout->n_in-1];
      v->decomposable = TRUE;
      printf("New node: %s, %d\n", v->name, v->ndesc);

      lvt = (list_of_vars *) malloc(sizeof(*lvt));
      lvt->var = v;
      lvt->next = NULL;
      for (l2=lv; l2->next!=NULL; l2=l2->next);
      l2->next = lvt;
    }
  } while (do_loop);
}

/****************************************************************************
COMPOSITION
composes a node so that this has only one layer of successors. This is
the reverse process of decomposition.
****************************************************************************/

void compose(Str255 vname)
{
  var_type *v;
  char b = TRUE;
  int i;

  FIND_VAR(v, vname);
  while (b) {
    for (b=FALSE, i=0; (i<v->famsout->n_in) && (!b); i++)
      b = v->famsout->in[i]->famsout != NULL;
    i--;
    if (b) {
/*      printf("del var %s\n", v->famsout->in[i]->name); */
      del_var(&variables, v->famsout->in[i]->name);
    }
  }
}

/****************************************************************************
Checking for and removing redundant variables or redundant values of
variables

This test is done through decomposition. Having a set of variables
y=f(a,b,c) and decomp({b,c},{a}) yields, for a, all columns that are
the same (column count == 1), then this variable is redundant.

If, for a, decomp({b,c},{a}) yields column count n, and n<|a|, then
some descriptors of a are redundant.
****************************************************************************/

void count_col_prepare()
{
  int i;

/*  for (i=0; i<nin && !nondis_dec; i++) nondis_dec = is_both[i]; */

  for (ncol=1, nrow=1, i=0; i<nin; i++) {
    if (is_col[i]) ncol *= sf->in[i]->ndesc;
    if (is_row[i]) nrow *= sf->in[i]->ndesc;
  }

  if (use_dm) {
    for (nrow=1, i=0; i<nin; i++)
      if (is_row[i]) nrow *= sf->in[i]->ndesc;
    dm = c_matrix_ini(nrow, ncol, CUNDEF);
    drow = c_vector(nrow);
  }
  dcol = c_vector(ncol);
  rv = c_vector(nin); rvr = c_vector(nin); rvc = c_vector(nin);
  if (deb_dec) {printf("Part: "); print_partition(stdout);  printf("\n");}
}

void count_col_terminate()
{
  if (use_dm) {
    free_c_matrix(dm, nrow); FREE(drow);
  }

  FREE(dcol); FREE(rv); FREE(rvr); FREE(rvc); 
}

/* remove_>dep_fn: this changes the dependecy of fam f, so that it does
   not include its n-th dependent variable any more */

void remove_dep_fn(fams *f, int n)
{
  int i, m;
  char *att;
  fams *tmpf;			/* backup fam */
  char b;
}

/* check_redundant: displays if variables are redundant or have
   redundant descriptors. If remove_red is TRUE, then redundancies are
   removed. Special care has to be taken, because sf changes, and
   number of variables changes. Note that new variables, which are not
   redundant, are added to the end of sf->in list and these are not
   needed to be checked. */

char check_redundant(Str255 vname, char remove_red, list_of_vars *lv)
{
  int n, i, nc;
  Str255 s;
  list_of_vars *ll;
  char b=TRUE;
  int nreplaced = 0;		/* how many vars which included
				   redundant descriptors were replaced */
  char found = FALSE;

				/* does the node exist? */
  FIND_VAR(sv, vname);
  sf = sv->famsout;
  if (!use_dm) {
    if (dec_speed==DEC_SLOW) build_hash_table(sf);
    else {
      int i;
      rule_list *rl;
      nn_rules = count_rules(sf);
      rules = (rule_list **) malloc(nn_rules * sizeof(*rules));
      for (i=0, rl=sf->lrule; rl!=NULL; i++, rl=rl->next)
	rules[i] = rl;
    }
  }
  nin = sf->n_in;
  is_col = c_vector_ini(nin,FALSE); is_row = c_vector_ini(nin,FALSE);
  is_both = c_vector_ini(nin,FALSE);
  
  for (n=0; n<nin-nreplaced; n++) {
    if (lv!=NULL) {
      b = FALSE;
      for (ll=lv; ll!=NULL && !b; ll=ll->next)
	b = ll->var == sv->famsout->in[n];
    }
    if (b) {
      for (i=0; i<nin; i++) {
	is_col[i] = (i==n) ? TRUE : FALSE;
	is_row[i] = (i!=n) ? TRUE : FALSE;
      }
      
      count_col_prepare();
      count_entries();
      nc = sinfo.cm;
      
      printf("%s ", sf->in[n]->name);
      if (nc == 1) printf("redundant\n");
      else if (nc == sf->in[n]->ndesc) printf("ok\n");
      else if (nc < sf->in[n]->ndesc)
	printf("has %d of %d desc redundant\n",
	       sf->in[n]->ndesc - nc, sf->in[n]->ndesc);
      else printf("coloring error, ncol %d instead of %d\n", 
		  nc, sf->in[n]->ndesc);
      
      if (remove_red && nc < sf->in[n]->ndesc) {
	found = TRUE;
	if (nc == 1) split_node(nc, gen_var_name(), TRUE);
	else {
	  sprintf(s, "%s'", sf->in[n]->name);
	  split_node(nc, s, FALSE);
	  nreplaced++;
	}	
	nin = sf->n_in;
	--n;
	if (!use_dm) {
	  if (dec_speed==DEC_SLOW) build_hash_table(sf);
	  else {
	    int i;
	    rule_list *rl;
	    nn_rules = count_rules(sf);
	    rules = (rule_list **) malloc(nn_rules * sizeof(*rules));
	    for (i=0, rl=sf->lrule; rl!=NULL; i++, rl=rl->next)
	      rules[i] = rl;
	  }
	}
      }
      count_col_terminate();
    }
  }
  FREE(is_col); FREE(is_row);
  if (!use_dm) {
/*    if (dec_speed == DEC_SLOW) free_hash_table();
    else */ FREE(rules);
  }
  return found;
}

/* del_redundant_sorted: removes the redundancy by trying first the
   least informative attributes */



int vcmp_infor(var_type **a, var_type **b)
{
  var_type *v1 = *a, *v2 = *b;

  if (v1->infor > v2->infor) return 1;
  if (v1->infor < v2->infor) return -1;
  return 0;
}

/* del_redundant_sorted: removes redundancy by a predifined order of
   variables. There are two variants: the order is related to apriory
   information measure (or gini or gain), or these are checked every
   time something is removed. */

double *measure;
void set_measure(int mtype)
{
  switch (mtype) {
  case M_INFORM:
    measure = im_one;
    break;
  case M_GINI:
    measure = gini_one;
    break;
  case M_GR:
    printf("SET\n");
    measure = gr_one;
    break;
  case M_RELIEFF:
    measure = w_one;
    break;
  }
}

void del_redundant_sorted(Str255 idname, int mtype, list_of_vars *lv)
{
  fams *f;
  int n, i;
  var_type **v;
  list_of_vars llv;
  char c;

  if ((f = find_fam_vartable(idname)) == NULL) return;
  v = (var_type **) malloc(sizeof(*v) * f->n_in);
  if (FALSE) {
    informativity_one(f);
    set_measure(mtype);
    for (i=0; i<f->n_in; i++) {
      v[i] = f->in[i];
      v[i]->infor = measure[i];
    }
    qsort(&v[0], f->n_in, sizeof(var_type **), vcmp_infor);
    llv.prev = llv.next = NULL;
    n = f->n_in;
    for (i=0; i<n; i++) {
      llv.var = v[i];
      check_redundant(idname, TRUE, &llv);
    }
  }
  else {
    double min=1;
    int imin, nn;

    llv.prev = llv.next = NULL;
    nn = f->n_in;

				/* mark all to be checked */
    for (i=0; i<f->n_in; i++) f->in[i]->mark = TRUE;

    for (n=0; n<nn; n++) {
      informativity_one(f);
      set_measure(mtype);
      for (i=0; i<f->n_in; i++) {
	v[i] = f->in[i];
	v[i]->infor = measure[i];
      }
      qsort(&v[0], f->n_in, sizeof(var_type **), vcmp_infor);
      for (i=0; i<f->n_in; i++) 
	if (v[i]->mark) {
	  llv.var = v[i];
	  c = check_redundant(idname, TRUE, &llv);
	  if (!c) v[i]->mark = FALSE;
	  else break;
	}
      if (!c) break;
    }
  }
}

/****************************************************************************
join_nodes: joins several vars into a single var as specified
****************************************************************************/

list_of_vars *g_ndlv;

/* join_type: 0 both, 1 prepare, 2 do */

int ttry=0;

map_nd_colors(char *dcol, char *ndcol)
{
  int i;
/*  for (i=0; i<nc; i++) printf("%d", ndcol[i]); printf("\n"); */
  for (i=0; i<ncol; i++) if (dcol[i]!=CUNDEF) {
    dcol[i] = ndcol[dcol[i]];
/*    printf("%d", dcol[i]); */
  }
/*  printf("\n"); */
}

void join_nodes(list_of_vars *lv, list_of_vars *ndlv, 
		Str255 vname, Str255 newname, int join_type)
{
  list_of_vars *l;
  var_type *nv;
  int i;
  char b;

  if (join_type!=2) {
    nondis_dec = ndlv!=NULL;
    g_ndlv = ndlv;
    strcpy(joined_var_name, newname);
    orig_n_dis_hi = n_dis_hi;
    FIND_VAR(sv, vname);
    sf = sv->famsout;
    if (!use_dm) {
      if (dec_speed == DEC_SLOW) build_hash_table(sf);
      else {
	int i;
	rule_list *rl;
	nn_rules = count_rules(sf);
	rules = (rule_list **) malloc(nn_rules * sizeof(*rules));
	for (i=0, rl=sf->lrule; rl!=NULL; i++, rl=rl->next)
	  rules[i] = rl;
      }
    }
    nin = sf->n_in;

    is_col = c_vector_ini(nin,FALSE); is_row = c_vector_ini(nin,TRUE);
    is_both = c_vector_ini(nin,FALSE);

    for (l=lv; l!=NULL; l=l->next) {
      for (b=FALSE, i=0; i<nin && !b; i++) b = sf->in[i] == l->var;
      if (b) {is_col[i-1] = TRUE; is_row[i-1] = FALSE;}
      else printf("warning: %s not an immediate successor of %s\n",
		  l->var->name, sv->name);
    }
  
    for (l=ndlv; l!=NULL; l=l->next) {
      for (b=FALSE, i=0; i<nin && !b; i++) b = sf->in[i] == l->var;
      if (b) {is_both[i-1] = is_col[i-1] = is_row[i-1] = TRUE;}
      else printf("warning: %s not an immediate successor of %s\n",
		  l->var->name, sv->name);
    }

    if (deb_dec) {printf("Join: "); print_partition(stdout); printf("\n");}
    for (i=0; i<nin; i++) if (is_both[i]) is_row[i]=FALSE;

    for (nvcol=0, i=0; i<nin; i++) if (is_col[i]) nvcol++;
    count_col_prepare();
    count_entries();
    nc = sinfo.cm;

    printf("xxx %d\n", nc);
    list_color("Join uses color", sinfo.dcol);

    if (nondis_dec) {
      nnc = nd_count_entries(nc);
      c_vector_copy(dcol, &tdcol, ncol);
/*      list_color("Before", dcol); */
      map_nd_colors(dcol, ndcol);
/*      list_color("After", dcol); */
    }
    color_analyze_start();
/*    if (join_type==1)  report_color(ncol, nrow, dcol);  */
  }

  if (join_type!=1) {
    if (nondis_dec) printf("colors %d -> %d\n", nc, nnc);
    else printf("colors %d\n", nc);
    list_color("Join uses color", dcol);
/*    nv = split_node(nondis_dec ? nnc : nc, joined_var_name, FALSE);  */
    nv = split_node(nondis_dec ? nnc : nc, joined_var_name, 
		    (nondis_dec ? nnc : nc) <= 1); 
    term_decompose_node();
				/* BBBB add the stuff for disjoint */
  }
}

/****************************************************************************
Fast decomposition using heuristic measures
****************************************************************************/
double get_gain(int i, int j, int measure, int m_type)
{
  if (m_type == 0)
    switch (measure) {
    case M_INFORM: return im_two[i][j] - (im_one[i]+im_one[j]);
    case M_GINI:   return gini_two[i][j] - (gini_one[i]+gini_one[j]);
    case M_GR:     return gr_two[i][j] - (gr_one[i]+gr_one[j]);
    case M_MDL:    return mdl_two[i][j] - (mdl_one[i]+mdl_one[j]);
    }
  else
    switch (measure) {
    case M_INFORM: return im_two[i][j]/(im_one[i]+im_one[j]);
    case M_GINI:   return gini_two[i][j]/(gini_one[i]+gini_one[j]);
    case M_GR:     return gr_two[i][j]/(gr_one[i]+gr_one[j]);
    case M_MDL:    return mdl_two[i][j]/(mdl_one[i]+mdl_one[j]);
    }
}

void fast_decompose(Str255 idname, int measure, int m_type)
{
  fams *f;
  int best_i, best_j, i, j;
  double m, best_m;
  list_of_vars lv1, lv2;

  lv1.next = &lv2;
  lv2.next = NULL;
  if ((f = find_fam_vartable(idname)) == NULL) return;
  while (f->n_in>2) {
    f = find_fam_vartable(idname);
    informativity_one(f);
    informativity_two(f);
    best_i=0; best_j=1; best_m = get_gain(best_i, best_j, measure, m_type);
    for (i=0; i<f->n_in-1; i++)
      for (j=i+1; j<f->n_in; j++)
	if (best_m < (m = get_gain(i, j, measure, m_type))) {
	  best_i = i; best_j = j; best_m = m;
	}
    lv1.var = f->in[best_i];
    lv2.var = f->in[best_j];
    join_nodes(&lv1, NULL, f->out->name, gen_var_name(), 0);
  }
}
