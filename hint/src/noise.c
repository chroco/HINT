/****************************************************************************
noise.c

Rutines to experiment with noise removal.
****************************************************************************/

#define CHECK(a,b,c) {if(a<=b || b<0) {printf("%s\n", c); exit(0); }}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#define GL extern
#include "sds.h"

char sel_mfreq = FALSE;		/* sele the most frequent class, or
				   the one with the smallest error */
double merror[255];		/* saves class errors */
double d_vector_min_max(double *v, double *v1, int len, int *pos);

int fold_n;
				/* this is from decomp.c */
extern fams *sf;
extern char *is_row, *is_col;
int ncol;			/* number of columns */
int ncolors;			/* number of current colors */
int nclass;			/* number of classes */
int *ccolor, *ncolor;		/* current and new set of colors */
int ncolors;			/* number of different colors */
double e_margin;		/* error of the previous table */
char *is_defined;		/* there is at least one entry in column? */
double *distr = NULL;		/* distribution by colors */
double *color_err;		/* error by color */
double *prior_class = NULL;		/* prior probability for classes */
char *changed = NULL;		/* is column's color changed? */

int n_rules =-1;		/* number of rules */
int n_rules_train;		/* number of rules used for merging, */
int n_rules_test;		/* other are for testing */

double n_items_train;		/* rule can contain several items */
double n_items_test;		/* this are computed through distributions */

rule_list **rules_l = NULL;	/* sorted rules to learn and to test */
rule_list **rules_t = NULL;

void pp()
{
  int i;
if (fold_n==2 && TRUE){
  for (i=0; i<n_rules_train; i++) {
    printf("%s%s    %s   %2d(%2d)  > ", 
	   rules_l[i]->train?"T":"F",
	   rules_l[i]->test?"T":"F",
	   rules_l[i]->mark?"T":"F",
	   ccolor[rules_l[i]->cr], rules_l[i]->cr);
    list_rule(sf, rules_l[i]->class, rules_l[i]->att);
  } 
  printf("--------\n");
  for (i=n_rules_train; i<n_rules; i++) {
    printf("%s%s    %s   %2d(%2d)  > ", 
	   rules_l[i]->train?"T":"F",
	   rules_l[i]->test?"T":"F",
	   rules_l[i]->mark?"T":"F",
	   ccolor[rules_l[i]->cr], rules_l[i]->cr);
    list_rule(sf, rules_l[i]->class, rules_l[i]->att);
  } 
  exit(0);
}
}

/* rule_cmp_used: compare two rules based on their usage info. This is
   used to sort the rules so that train come first. */

int rule_cmp_used(rule_list **a, rule_list **b)
{
  rule_list *r1 = *a, *r2 = *b;
  if (r1->train == r2->train) return 0;
  if (r1->train && !r2->train) return -1;
  return 1;
}

/* rule_cmp_rows: compare two rules (indexes), so that to take the order
   for is_row. if r1 is smaller */

int rule_cmp_rows(rule_list *r1, rule_list *r2)
{
  int i, j;

  for (i=0; i<sf->n_in; i++) 
    if (is_row[i]) {
      if (r1->att[i] < r2->att[i])
	return -1;
      else if (r1->att[i] > r2->att[i])
	return 1;
    }
  return 0;
}

int rule_cmp_rows_s(rule_list **a, rule_list **b)
{
  int i, j;
  rule_list *r1 = *a, *r2 = *b;

  for (i=0; i<sf->n_in; i++) 
    if (is_row[i]) {
      if (r1->att[i] < r2->att[i])
	return -1;
      else if (r1->att[i] > r2->att[i])
	return 1;
    }
  return 0;
}

int rule_cmp_rcs(rule_list **a, rule_list **b)
{
  int i, j;
  rule_list *r1 = *a, *r2 = *b;

  for (i=0; i<sf->n_in; i++) 
    if (is_row[i]) {
      if (r1->att[i] < r2->att[i])
	return -1;
      else if (r1->att[i] > r2->att[i])
	return 1;
    }
  if (ccolor[r1->cr] > ccolor[r2->cr]) return 1;
  if (ccolor[r1->cr] < ccolor[r2->cr]) return -1;
  return 0;
}

void set_row_markers(int from, int to, rule_list **rules)
{
  int i;

  for (i=from; i<to-1; i++)
    rules[i]->mark = rule_cmp_rows(rules[i], rules[i+1]) != 0;
  rules[to-1]->mark = TRUE;
}

void sort_rules_rows(int from, int to)
{
  qsort(&rules_l[from], to-from, sizeof(rule_list **), rule_cmp_rows_s);
  set_row_markers(from, to, rules_l);
}

int *sort_color;

int rule_cmp_colors(rule_list *r1, rule_list *r2)
{
  int i, j;

  if (sort_color[r1->cr] > sort_color[r2->cr]) return 1;
  if (sort_color[r1->cr] < sort_color[r2->cr]) return -1;
  return 0;
}

int rule_cmp_colors_s(rule_list **a, rule_list **b)
{
  rule_list *r1 = *a, *r2 = *b;

  if (sort_color[r1->cr] > sort_color[r2->cr]) return 1;
  if (sort_color[r1->cr] < sort_color[r2->cr]) return -1;
  return 0;
}

sort_rules_color(int *color, int from, int to)
{
  int i=from;
  int first = from, last;

  sort_color = color;

  while (i<to) {
    for (last=first; last<to && !rules_l[last]->mark; last++);
    i = last;
    if ((last>=to) || (first>=to)) {
      printf("s error\n"); exit(0);
    }

/*    printf("xxx %d %d (%d %d)\n", first, last, from, to); */
    qsort(&rules_l[first], last-first+1, sizeof(rules_l), rule_cmp_colors_s);
    first = ++i;
  }
				/* needs to readjust row markers */
  set_row_markers(from, to, rules_l);

/*  for (i=0; i<n; i++) {
    printf(" %2d %2d %s   ", sort_color[rules_l[i]->cr], 
	   rules_l[i]->cr, rules_l[i]->mark?"T":"F");
    list_rule(sf, rules_l[i]->class, rules_l[i]->att);
  } */
}

/* laplace_err: derive laplace error for given color matrix */

#define PLOGP(x) (((x)==0)?0:(- (x) * log(x) / CLOG2))


double laplace_color_err(int *color, int scolor)
{
  int i, j, k, pos;
  double e, max, error = 0., d;


  d_vector_set(distr, nclass, 0.);
  for (i=0; i<n_rules_train; i++) {
    if (color[rules_l[i]->cr] == scolor)
      d_vector_plus(distr, rules_l[i]->dist, nclass);

    if (rules_l[i]->mark) {
      max = d_vector_max_max(distr, prior_class, nclass, &pos);
      if (sel_mfreq) {
	if (max>0) {
	  /* probability of the cell * error */
	  d = d_vector_sum(distr, nclass);
	  switch (noise_handling) {
	  case n_laplace:
	    e = 1. - (max+1.) / (double)(d + nclass);
	    break;
	  case n_entropy:
	    for (j=0, e=0; j<nclass; j++) e += PLOGP(distr[j]/d);
	    break;
	  case n_mprob:
	    e = 1. -(max + prior_class[pos] * m_param) / (double)(d + m_param);
	    break;
	  }
	  error += d * e;
	  if (deb_dec>10) 
	    printf("  %3.1lf (%3.1lf) of %3.1lf is %6.2lf, item %lf\n", 
		   max, prior_class[pos], d, e, n_items_train);
	  if (e<=0) {printf("major error e=0\n"); exit(0);}
	}
      }
      else {
	if (max>0) {
	  d = d_vector_sum(distr, nclass);
	  for (k=0; k<nclass; k++) {
	    merror[k] = 1. - (distr[k] + prior_class[k] * m_param) 
	      / (double)(d + m_param);
	  }
	  e = d_vector_min_max(&merror[0], prior_class, nclass, &pos);
	  error += d * e;
	}
      }

      d_vector_set(distr, nclass, 0.);
    }
  }
  error = error / n_items_train;
  return error;
}

double laplace_err(int *color, int ncolors)
{
  int i;
  double e = 0., d;
  
  for (i=0; i<ncolors; i++) {
    if (deb_dec>10) printf("color %d\n", i);
    d = laplace_color_err(color, i);
    e += d;
    if (deb_dec>10) printf("color %d error %6.3lf\n", i, d);
  }
  return e;
}

/* select_pair: tries out all the couples to see which minimizes the
   error best, uses  */

void select_pair(int *best_i, int *best_j, double *new_e)
{
  int i, j, k;
  double e, d;
  double best_e = 99999.;

  if (ncolors==1) return;

				/* get laplace of original table */
  if (deb_dec>7) {
    e = laplace_err(ccolor, ncolors);
    printf("\nMargin (%d): %10.8lf (%10.8lf)\n", ncolors, e, e_margin);
    printf("       ");
    for (i=1; i<ncolors; i++) printf("  %2d   ", i); printf("\n");
  }
  for (i=0; i<ncolors-1; i++) {
    if (deb_dec>7) {printf("%-2d     ", i); pspace(7*i);}
    for (j=i+1; j<ncolors; j++) {

      for (k=0; k<ncol; k++)
	if (ccolor[k]==j) {
	  changed[k] = TRUE;
	  ccolor[k] = i;
	}
	else changed[k] = FALSE;
				/* compute the error for these two columns */
      e = laplace_color_err(ccolor, i);
				/* add the error of other columns*/
      d = e;
      for (k=0; k<ncolors; k++)
	if (k != i && k != j) e += color_err[k];

      if (e<best_e) {
	*best_i = i; *best_j = j; best_e = e; *new_e = d;
      }

      if (deb_dec>7) printf("%6.4lf ", e);
      for (k=0; k<ncol; k++) if (changed[k]) ccolor[k]=j;
    }
    if (deb_dec>7) printf("\n");
  }
  if (deb_dec>6) 
    printf("Best %6.4lf, try:  merge %d %d\n", best_e, *best_i, *best_j);
}

/* initialization: init_init_merge() should be called as soon as the
   bound and free sets and sf are set. init_merge() should be called
   after the set is split to learn and test set, i.e., 
   before every set of merges (bef. testing of each fold). */

void init_init_merge()
{
  int i, j;
  char c;
  double dd=0;
  rule_type rtp;
  rule_list *rl;

  nclass = sf->out->ndesc;
				/* count the rules and*/
				/* set prior probabilities for classes */
  FREE(prior_class); prior_class = d_vector_ini(nclass, 0);
  reset_next_rule(sf);

/*  for (i=0; get_next_rule(sf, &c, &rtp); i++)
    ++prior_class[c];  */

  for (i=0, rl=sf->lrule; rl!=NULL; i++, rl=rl->next) {
    for (j=0; j<nclass; j++) {
      prior_class[j] += rl->dist[j];
      dd += rl->dist[j];
    }
  }

  if (i!=n_rules) {
    FREE(rules_l);
    n_rules = i;
    rules_l = (rule_list **) malloc(n_rules * sizeof(*rules_l));
    rules_t = (rule_list **) malloc(n_rules * sizeof(*rules_t));
    if (rules_l == NULL) {printf("rules_l\n"); exit(0);}
  }
  if (deb_dec>6) printf("rules %d\n", n_rules);

  for (i=0, rl=sf->lrule; rl!=NULL; i++, rl=rl->next)
    rules_l[i] = rl;

/*  printf("%d %d %d\n", (int)rules_l,n_rules,rules_l[n_rules-1]->train);*/

/*  for (i=0; i<nclass; i++) prior_class[i] /= dd;  */

				/* use Laplace */
  for (i=0; i<nclass; i++) 
    prior_class[i] = (prior_class[i] + 1) / (dd + nclass);

/*  d_vector_print(prior_class, nclass); printf("\n");
  exit(0); */


  if (deb_dec>4) {
    printf("Prior: "); d_vector_print(prior_class, nclass); printf("\n");
  }

				/* find number of columns */
  for (ncol=1, i=0; i<sf->n_in; i++)
    if (is_col[i]) ncol *= sf->in[i]->ndesc;
  if (deb_dec>6) printf("columns: %d\n", ncol);

  ccolor = i_vector(ncol);
  ncolor = i_vector(ncol);

				/* some allocation stuff */
  FREE(distr); distr = d_vector_ini(nclass, 0);
  FREE(changed); changed = c_vector(ncol);
}

void init_merge()
{
  int i, j;
  double d;

				/* find n intems */
  n_items_test = n_items_train = 0.;
  for (i=0; i<n_rules; i++) {
    for (d=0., j=0; j<nclass; j++) d += rules_l[i]->dist[j];
    if (rules_l[i]->train) n_items_train += d;
    else n_items_test += d;
  }

				/* find undefined columns */
  is_defined = c_vector_ini(ncol, FALSE);
  reset_next_rule(sf);
  for (i=0; i<n_rules; i++) {
    rules_l[i]->cr = satt2indx(sf, rules_l[i]->att, is_col);
    is_defined[rules_l[i]->cr] = TRUE;
  }
  if (deb_dec>6) {
    for (i=0; i<ncol; i++) printf("%s ", is_defined[i]?"T ":"F ");
    printf("\n");
  }

				/* sort the rules by their rows */
  sort_rules_rows(0, n_rules_train);
				/* set the initial colors */
				/* each column separate color */
  for (i=0, ncolors=0; i<ncol; i++)
    if (is_defined[i]) ccolor[i] = ncolors++;
    else ccolor[i] = -1;
  if (deb_dec>6) {
    for (i=0; i<ncol; i++) if (ccolor[i]>=0)
      printf("%-2d ", ccolor[i]); else printf("-  ");
    printf("\n");
  }
/*  printf("ncol %d\n", ncol); */

				/* error for each of the colors (separate) */
  color_err = d_vector(ncolors);
  for (i=0, e_margin = 0.; i<ncolors; i++) {
    color_err[i] = laplace_color_err(ccolor, i);
    e_margin += color_err[i];
    if (deb_dec>6) printf("%6.4lf ", color_err[i]);
  }
  if (deb_dec>6) printf("\n");

/*  for (i=0; i<n_rules; i++) {
    printf("%s->%d\n", rules_l[i]->train?"T":"F", rules_l[i]->class);
  }  */
}


/* use_all_rules(): sets train so that all rules are used for
   merging. */

void use_all_rules()
{
  int i, j;

  n_rules_train = n_rules;
  n_rules_test = 0;
  for (i=0; i<n_rules; i++) {
    rules_l[i]->train = TRUE;
    rules_l[i]->test = FALSE;
  }
}

/* test: sorts the rules, and sets the inital row error */

void test(double p, double pincl)
{
  int i, j;

  init_init_merge();

				/* make a selection of rules, this
                                   will move to some other place */

				/* variant 1: exclusive */
/*  for (n_rules_train=0, n_rules_test=0, i=0; i<n_rules; i++) {
    if (rnd1e() > (1.-p)) {
      rules_l[i]->train = TRUE;
      rules_l[i]->test = FALSE;
      ++n_rules_train;
    }
    else {
      rules_l[i]->train = FALSE;
      rules_l[i]->test = TRUE;
      ++n_rules_test;
    }
  } */

				/* variant 2: inclusive by random */

/*  for (n_rules_train=0, i=0; i<n_rules; i++) {
    rules_l[i]->train = FALSE;
    rules_l[i]->test = FALSE;
  }

  for (n_rules_train=0, i=0; i<n_rules; i++)
    if (rnd1e() > (1.-p)) {
      rules_l[i]->train = TRUE;
      ++n_rules_train;
    }
  for (n_rules_test=0, i=0; i<n_rules; i++) 
    if (rnd1e() > p) {
      rules_l[i]->test = TRUE;
      ++n_rules_test;
    } */


				/* variant 3: inclusive by definition */

  for (n_rules_train=0, n_rules_test=0, i=0; i<n_rules; i++) {
    if (rnd1e() > (1.-p)) {
      rules_l[i]->train = TRUE;
      rules_l[i]->test = FALSE;
      ++n_rules_train;
    }
    else {
      rules_l[i]->train = FALSE;
      rules_l[i]->test = TRUE;
      ++n_rules_test;
    }
  }
  for (i=0; i<n_rules; i++)
    if (rules_l[i]->train && rnd1e() > (1.-pincl)) {
      rules_l[i]->test = TRUE;
      ++n_rules_test;
    }

  for (i=0, j=0; i<n_rules; i++) 
    if (rules_l[i]->test) rules_t[j++] = rules_l[i];

/*  use_all_rules(); */
  if (deb_dec>5) printf("used rules %d (%6.2lf%%)\n", 
	 n_rules_train, n_rules_train/(double)n_rules * 100);
  qsort(rules_l, n_rules, sizeof(rules_l), rule_cmp_used);

  init_merge();
}

/****************************************************************************
What about some clustering. We first have to define a measure of
likeness.
****************************************************************************/

double get_new_e_margin(int a, int b, double new_err)
{
  double d=new_err;
  int i;

  for (i=0; i<ncolors; i++)
    if (i!=a && i!=b) d+=color_err[i];
  return d;
}

void merge_colors(int a, int b, double new_err)
{
  int i, j, c, d;
  double e;

  c = MMIN(a,b); d = MMAX(a,b);
  for (i=0; i<ncol; i++) if (ccolor[i]==d) ccolor[i]=c;
  for (i=0; i<ncol; i++) if (ccolor[i]>d) --ccolor[i];
  for (i=d+1; i<ncolors; i++) color_err[i-1] = color_err[i];

  --ncolors;
  if (new_err>=0) {
    color_err[c] = new_err;
  }
  else
    color_err[c] = laplace_color_err(ccolor, c);

  for (i=0, e_margin = 0.; i<ncolors; i++)
    e_margin += color_err[i];

/*  sort_rules_color(ccolor); */

  if (deb_dec>6) {
    printf("\n");
    for (i=0; i<ncol; i++) printf("%2d ", ccolor[i]); printf("\n");
  }

  if (new_err == -1.) select_pair(&i, &j, &e);
}

/* repeat_merge: tryes to repeat the merging until some criterion is
   satisfied */

double repeat_merge(int fixed)
{
  int i, j;
  double e, e_delta;

  fprintf(lfile, "TERR %3d %6.4lf w m=%lf\n", ncolors, e_margin, m_param);
  while (ncolors>fixed && ncolors>1) {
    select_pair(&i, &j, &e);
    e_delta = e_margin - get_new_e_margin(i, j, e);

    if (e_delta>-1e-10) e_delta=0; /* can be numerical error */

    if (deb_dec>2 && e_delta<0) printf("Delta rejected: %8.5g\n", e_delta);
    if ((fixed<0) && (e_delta<0)) break; 
    merge_colors(i, j, e);
    if (deb_dec>6) printf("Delta %2d: %6.4lf\n", ncolors, e_delta);
    fprintf(lfile, "MERG %d %d\n", i, j);
    fprintf(lfile, "TERR %3d %6.4lf (d %6.4lf)\n", ncolors, e_margin,e_delta);
  }

  if (deb_dec>5) {
    for (i=0; i<ncol; i++) printf("%d ", ccolor[i]); printf("\n");
  }
  return e_margin;
}

/*
To draw the resulting graph do:

% grep 'TERR' log | awk '{print $2,$3}' | graph | xplot
% grep 'TERR' log | awk '{if (NR!=1) print $2,$3-x; x=$3}' | graph | xplot

*/

int find_rule_row_color(rule_list *rl)
{
  int lo = 0, hi = n_rules_train - 1, mid, i;
  
  if (rule_cmp_rows(rules_l[lo], rl) == 0 &&
      rule_cmp_colors(rules_l[lo], rl) == 0) return lo;

  if (rule_cmp_rows(rules_l[hi], rl) == 0 &&
      rule_cmp_colors(rules_l[hi], rl) == 0) return hi;

  while (TRUE) {
    mid = (lo + hi) / 2;
/*    printf("xx %d\n", mid); */

    if (mid == lo || mid == hi) return -1; /* not found */
    
    i = rule_cmp_rows(rules_l[mid], rl);

/*printf("xx lala\n"); */
    if (i == 0) i = rule_cmp_colors(rules_l[mid], rl);
/*printf("xx lili\n");*/
    if (i == 0) return mid;
    else if (i < 0) lo = mid;
    else hi = mid;
  }

/*  printf("find %d\n", rule_cmp_rows(rules_l[0], rl));
  printf("find %d\n", rule_cmp_rows(rules_l[n_rules_train-1], rl));*/
}

/* find_rule_row_color_start: find the rule of the list that is the first
   one that starts with the right color. It assumes that start is
   already a position that points to such rule, only it has to find
   the begining of the sequence of such rules */

int find_rule_row_color_start(rule_list *rl, int start)
{
  int last_start;

  last_start = start;
  for (--start; start >= 0; --start) {
/* printf("xx3 lolo %d %d\n", start, n_rules); */
    if (rule_cmp_colors(rules_l[start], rl)!=0 ||
	rule_cmp_rows(rules_l[start], rl)!=0)
      return last_start;
    last_start = start;
  }
  return 0;
}

/* set_distrib_row_color: for the specific row and specific color,
   this sets the distribution of classes to distr */

void set_distrib_row_color(int *start)
{
  int i = *start, j;
  d_vector_set(distr, nclass, 0.);
  for (; i<n_rules_train && 
       rule_cmp_colors(rules_l[*start], rules_l[i])== 0; i++) {
/*    distr[rules_l[i]->class]++; */
    for (j=0; j<nclass; j++)
      distr[j] += rules_l[i]->dist[j];

    if (rules_l[i]->mark) {
      i++;
      break;
    }
  }

if (fold_n==2 && FALSE){
  printf("distr to    ");
  printf("    %s   %d(%2d)  > ", rules_l[i-1]->mark?"T":"F",
	 ccolor[rules_l[i-1]->cr], rules_l[i-1]->cr);
  list_rule(sf, rules_l[i-1]->class, rules_l[i-1]->att);
}
  *start = i;
}

double get_class_error()
{
  int r, n, i, j, nerrors=0, last_r=0, max_c, default_c;
  double e;
  char found = FALSE;

  d_vector_max(prior_class, nclass, &default_c);
/*  sort_rules_rows(n_rules_train, n_rules);  */
  sort_rules_color(ccolor, 0, n_rules_train);

/*  for (i=0, j=0; i<n_rules; i++) 
    if (rules_l[i]->test) rules_t[j++] = rules_l[i]; */
  qsort(&rules_t[0], n_rules_test, sizeof(rules_t), rule_cmp_rcs);

/* printf("xxx %d\n", n_rules_test); */
  for (n = 0, r = 0; r < n_rules_test; r++)
    if (rules_t[r]->test || TRUE) {
      if (r==n_rules_train ||
	!(rule_cmp_colors(rules_t[r], rules_t[last_r])==0 &&
	  rule_cmp_rows(rules_t[r], rules_t[last_r])==0)) {
      found = FALSE;

      while (n < n_rules_train && 
	     rule_cmp_rows(rules_t[r], rules_l[n])>0)
	n++;

      if (n < n_rules_train) {
	while (n < n_rules_train && 
	       rule_cmp_rows(rules_t[r], rules_l[n])==0 &&
	       rule_cmp_colors(rules_t[r], rules_l[n])>0)
	  n++;

	found = rule_cmp_rows(rules_t[r], rules_l[n])==0 &&
	  rule_cmp_colors(rules_t[r], rules_l[n])==0;

	if (found) {
	  set_distrib_row_color(&n);
	  d_vector_max_max(distr, prior_class, nclass, &max_c);
	}
      }
    }      
    
    if (found) {
      if (rules_t[r]->class != max_c)
	++nerrors;
    }
    else {
      if (rules_t[r]->class != default_c)
	++nerrors;
    }
    last_r = r;
  }
  e = nerrors/(double)(n_rules-n_rules_train) * 100.0;
  if (deb_dec>6) 
    printf("Error: %6.2lf (%d of %d)\n", e, nerrors, n_rules-n_rules_train);

  return e;
}

/* set the fold indexes: init_folds() intializes these indices and
   stores them in fold_mark. set_fold uses this, and sets the train &
   test fields of the rule accordingly */

int n_folds = 3;		/* number of folds */
double p_sample_l = 0.7;
double p_sample_mut = 0.3; /* probab for learn and learn*test */

char *fold_mark = NULL;		/* which fold does element belong to*/

void init_folds()
{
  int *fold_elements, i, j, k, l;
  char *fold_indx;

  FREE(fold_mark); fold_mark = c_vector_ini(n_rules, 0);
  fold_elements = i_vector_ini(n_folds, n_rules / n_folds);
  for (i=n_rules-((int)n_rules/n_folds)*n_folds; i>0; i--)
    ++fold_elements[(int) (rnd1e() * n_folds)];
  
  fold_indx = c_vector(n_folds);
  for (i=0; i<n_folds; i++) fold_indx[i] = i;

  j = n_folds;
  for (i=0; i<n_rules; i++) {
    fold_mark[i] = fold_indx[k = (int)(rnd1e() * j)];

    --fold_elements[k];
    if (fold_elements[k] == 0) {
      for (l=k; l<n_folds-1; l++) {
	fold_elements[l] = fold_elements[l+1];
	fold_indx[l] = fold_indx[l+1];
      }
      --j;
    }
  }

/*  i_vector_set(fold_elements, n_folds, 0);
  for (i=0; i<n_rules; i++) fold_elements[fold_mark[i]]++;
  i_vector_print(fold_elements, n_folds);
  printf("\n"); */

  FREE(fold_elements);
}

void set_learn_test(int fold)
{
  int i, j;

  if (test_type == t_cv) {
    for (n_rules_train=0, n_rules_test=0, i=0; i<n_rules; i++) {
      if (fold_mark[i] != fold) {
	rules_l[i]->train = TRUE;
	rules_l[i]->test = FALSE;
	++n_rules_train;
      }
      else {
	rules_l[i]->train = FALSE;
	rules_l[i]->test = TRUE;
	++n_rules_test;
      }
    }
  }
  else {
    for (n_rules_train=0, n_rules_test=0, i=0; i<n_rules; i++) {
      if (rnd1e() > (1.-p_sample_l)) {
	rules_l[i]->train = TRUE;
	rules_l[i]->test = FALSE;
	++n_rules_train;
      }
      else {
	rules_l[i]->train = FALSE;
	rules_l[i]->test = TRUE;
	++n_rules_test;
      }
    }
    for (i=0; i<n_rules; i++)
      if (rules_l[i]->train && rnd1e() > (1.-p_sample_mut)) {
	rules_l[i]->test = TRUE;
	++n_rules_test;
      }
  }

  for (i=0, j=0; i<n_rules; i++) 
    if (rules_l[i]->test) rules_t[j++] = rules_l[i];
  
  if (deb_dec>6) printf("used rules %d (%6.2lf%%)\n", 
	 n_rules_train, n_rules_train/(double)n_rules * 100);
  qsort(rules_l, n_rules, sizeof(rules_l), rule_cmp_used);
}

double eval_m_param(double *c)
{
  double e, error=0;
  int fold;

  *c = 0;
  if (test_type == t_cv) init_folds();
  for (fold=0; fold<n_folds; fold++) {
    fold_n = fold;
    set_learn_test(fold);
    init_merge();
    repeat_merge(-1);
    error += e = get_class_error();
    *c += ncolors;
    if (deb_dec>4) printf("Fold %d c=%d e=%6.2lf\n", fold, ncolors, e);
  }
  error /= n_folds; *c /= n_folds;
  if (deb_dec>3) printf("M %4.2lf e=%6.2lf c=%4.1lf\n", m_param, error, *c);
  return error;
}

/* Search for the m-probability estimate that minimizes the error
   estimated by cross validation. We use the Golden Section Search
   algorithm as described in Numerical Recipes, pp 400. */

double min_e;
double min_m;
double min_c;

double gss_eval(double m)
{
  double c, e;
/*  if (m<0) return 1./(double) n_folds * 100.; */
  m_param = m;
  e = eval_m_param(&c);
  if (e<min_e) {
    min_e = e; min_m = m; min_c = c;
  }
  return e;
}

double gss_find(double ax, double bx, double *e, double *c)
{
  double cx, fa, fb, fc;
  double m, g;
/*  double golden(double ax, double bx, double cx, double (*f)(double),
		double tol, double *xmin); */

  min_e = 100.;
  if (1) {
    mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, gss_eval);
  }
  else {
    ax = 0;   fa = gss_eval(ax);
    bx = 0.5; fb = gss_eval(bx);
    cx = 3;   fc = gss_eval(cx); 
				/* sort everything! */
  }
/* printf("end bracket %lf %lf %lf\n", ax, bx, cx); */
  min_e = fb; min_m = bx;
/*  g = golden(ax, bx, cx, gss_eval, 0.3, &m); */
/*  printf("mm %lf %lf\n", m, g); */
  *e = min_e; *c = min_c; return min_m;
}

void test_classify()
{
  double c, m, e;

  init_init_merge();
/*  for (m_param=0.2; m_param<3.; m_param += 0.2)
    eval_m_param(&c); */
  m = gss_find(0.2, 2, &e, &c);
  printf("Found %5.2lf e= %6.3lf c=%4.1lf\n", m, e, c);

  use_all_rules();
  init_merge();
  repeat_merge(-1);
  if (deb_dec>3) {
    int i; for (i=0; i<ncol; i++) printf("%d ", ccolor[i]); printf("\n");
  }
}

int get_partition_error(double *e, double *m)
{
  double c;
  int i;
  
  init_init_merge();
  *m = gss_find(0.1, 1, e, &c);
  if (deb_dec>3) printf("Found %5.2lf e= %6.3lf c=%4.1lf\n", *m, *e, c);

  use_all_rules();

  init_merge();
  repeat_merge(-1);
  if (deb_dec>3) {
    for (i=0; i<ncol; i++) printf("%d ", ccolor[i]); printf("\n");
  }

  for (i=0; i<ncol; i++)
    if (ccolor[i]>=0) dcol[i] = ccolor[i]; else dcol[i]=CUNDEF;
  return ncolors;
}

int get_partition_m(double *e)
{
  double c;
  int i;

  init_init_merge();
  use_all_rules();
  init_merge();
  *e = repeat_merge(-1);
  if (deb_dec>3) {
    for (i=0; i<ncol; i++) printf("%d ", ccolor[i]); printf("\n");
  }

  for (i=0; i<ncol; i++)
    if (ccolor[i]>=0) dcol[i] = ccolor[i]; else dcol[i]=CUNDEF;
  return ncolors;
}

/****************************************************************************
test_merge: tests the merge, so that it varies m, and checks to how
many classes it derives. Expects the command 
>> test NUM
to be run before and to set the learning/test set
****************************************************************************/

void test_merge(double from, double to, double step)
{
  double e;
  for (m_param=from; m_param<=to; m_param+=step) {
    init_merge();
    repeat_merge(-1);
    e = get_class_error();
    printf("m=%5.2lf  c=%2d  e=%7.3lf\n", m_param, ncolors, e);
  }
}

/****************************************************************************
Derivation of tables y=F(A, G(B)) with distributions
****************************************************************************/

void derive_G(fams *F, fams *G, fams *oldf, char remove_dep)
{
  int i, j, k, r;
  double d;
  fams *tmp = sf;
  int n = F->n_in;
  char *att = c_vector(n);
  double *cumul;		/* cumulative distribution by color */
				/* to find out the default rule for H */

  sf = oldf;

/*  printf("Derive G %s with %d\n", G->out->name, n_rules); */
  sort_rules_color(ccolor, 0, n_rules);
  cumul = d_vector_ini(ncolors, 0.);

  i = 0;
  while (i<n_rules) {
    set_distrib_row_color(&i);
    r = i-1;
    cumul[ccolor[rules_l[r]->cr]] += d_vector_sum(distr, nclass);
/*    d_vector_normalize(distr, nclass); */
    
    for (k=0, j=0; k<oldf->n_in; k++)
      if (is_row[k]) 
	att[j++] = rules_l[r]->att[k];
    if (!remove_dep) att[j] = ccolor[rules_l[r]->cr];

    if (sel_mfreq) {
      d_vector_max(distr, nclass, &k);
    }
    else {
      d = d_vector_sum(distr, nclass);
      for (k=0; k<nclass; k++) {
	merror[k] = 1. - (distr[k] + prior_class[k] * m_param) 
	  / (double)(d + m_param);
      }
      d_vector_min_max(merror, prior_class, nclass, &k);
    }

    set_drule(F, att, k, rules_l[r]->rtp, distr);
  }

  if (!remove_dep) {
    printf("Cumul by color: "); d_vector_print(cumul, ncolors); printf("\n");

    if (sel_mfreq) d_vector_max(cumul, ncolors, &i);
    else {
      d = d_vector_sum(cumul, nclass);
      for (k=0; k<nclass; k++) {
	merror[k] = 1. - (cumul[k] + prior_class[k] * m_param) 
	  / (double)(d + m_param);
      }
      d_vector_min_max(merror, prior_class, nclass, &i);
    }


    G->out->c_default = i;
/*    d_vector_normalize(cumul, ncolors); */
    G->out->c_apriory = d_vector(ncolors);
    memcpy(G->out->c_apriory, cumul, sizeof(double) * ncolors);
  }

  sf = tmp;
  FREE(cumul); FREE(att);
}

/****************************************************************************
Get the initial error of the example set
****************************************************************************/

double estimate_set_error(var_type *v)
{
  fams *f=v->famsout;
  double d, dd=0, max, e=0;
  int pos, nclass;
  rule_list *rl;

  if (f==NULL) return 0;
  nclass = v->ndesc;
  for (rl=f->lrule; rl!=NULL; rl=rl->next) {
    max = d_vector_max_max(rl->dist, v->c_apriory, nclass, &pos);
    d = d_vector_sum(rl->dist, nclass);
    dd += d;
    e += 1. - (max + v->c_apriory[pos] * m_param) / (double)(d + m_param);
  }
  printf("yyy %lf\n", dd);
  return e / dd;
}
