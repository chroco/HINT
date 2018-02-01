/****************************************************************************
heuristics.c

Heuristic measures for decomposition or for estimation of godness of
attributes. The code includes:

  * contingency table derivation
  * derivation of information measure (from contingency table)
  * relieff estimation

****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#define GL extern
#include "sds.h"

double *im_one=NULL, **im_two=NULL;	/* information measure for attributes
				   im_one - single attribute 
				   im_two - combination of two attributes */
double *w_one=NULL, **w_two=NULL; /* relieff measure for attributes */
double *gr_one=NULL, **gr_two=NULL;	/* gain-ratio measure */
double *gini_one=NULL, **gini_two=NULL;	/* gini measure */
double *mdl_one=NULL, **mdl_two=NULL; /* MDL measure */

extern fams *sf;		/* table of rules for which to derive
				   different measures */
extern int nin;			/* number of input variables to that fam */

void relieff(fams *f);

/****************************************************************************
INFORMATIVITY & GAIN-RATIO MEASURE

The computation of the informativity of a single attribute and on a
combination of two attributes. Computation is based on condingency
tables.
See Mingles: Decision tree induction, Machine Learning 3: 319 -342
****************************************************************************/

double *clt=NULL, **con, *at=NULL;	/* components of contingency table */

/* using Kononeko's notation
   - clt[i]    = n_i.  number of instances from class Ci 
   - con[j][i] = n_ij  number of instances from class Ci with j-value of
                       given attribute
   - at[j]     = n.j   number of instances with j-th value of given att

   - mcon = H_CA
   - mclt = H_C 
   - mat  = H_A 
*/

char use_doubles = TRUE;

int nclass, naval;		/* #examples, #classes, #attribute values */
double nr;			/* #examples */
double mclt, mcon, mat;		/* components in inf measure equation */

void compute_mclt_nr()
{
  int i, j;
  char class;
  rule_type rtp;

  FREE(clt);
  clt = d_vector_ini(nclass, 0);

  if (use_distribution && use_doubles) {
    rule_list *rl;

    for (rl=sf->lrule; rl!=NULL; i++, rl=rl->next) {
      for (j=0; j<nclass; j++) {
	clt[j] += rl->dist[j];
	nr += rl->dist[j];
      }
    }
  }
  else {
    reset_next_rule(sf);
    nr = 0;
    while (get_next_rule(sf, &class, &rtp)) {
      clt[class]++;
      nr++;
    }
  }
  for (mclt=0., i=0; i<nclass; i++) mclt -= NLOGN(clt[i]);
}

void compute_mat_mcon() 
{
  int j, k;
  char class;
  rule_type rtp;

  mat = mcon = 0.;
  for (j=0; j<naval; j++) {
    mat -= NLOGN(at[j]);
    for (k=0; k<nclass; k++)
      mcon -= NLOGN(con[j][k]);
  }
}

double compute_iv()
{
  int i;
  double d = 0, d1;

  for (i=0; i<naval; i++) {
    d1 = at[i] / (double) nr;
    d += NLOGN(d1);
  }
  return -d;
}

void build_contingency_table(char *sel)
{
  int i, indx;
  char *att, class;
  rule_type rtp;
  int j;
  double **d_matrix_ini(int nr, int nc, double ini);

  FREE(at);			/* BBBBB free con */

  for (naval=1, i=0; i<nin; i++) if (sel[i]) naval *= sf->in[i]->ndesc;
  at = d_vector_ini(naval, 0);
  con = d_matrix_ini(naval, nclass, 0);

  if (use_distribution && use_doubles) {
    rule_list *rl;

    for (rl=sf->lrule; rl!=NULL; i++, rl=rl->next) {
      indx = satt2indx(sf, rl->att, sel);
      if (indx>=0 && indx<naval) {
	at[indx]++;
	for (j=0; j<nclass; j++)
	  con[indx][j] += rl->dist[j];
      }
    }
  }
  else {
    reset_next_rule(sf);
    while (get_next_rulea(sf, &class, &rtp, &att)) {
      indx = satt2indx(sf, att, sel);
      if (indx>=0 && indx<naval) {
	at[indx]++;
	(con[indx][class])++;
      }
    }
  }


/*  printf("CT:\n");
  for (i=0; i<naval; i++) {
    for (j=0; j<nclass; j++) printf("%d ", con[i][j]);
    printf("\n");
  } */
}

double compute_gini()
{
  double ic=0, ir=0, d;
  int i, j;

  for (i=0; i<naval; i++) {
    for (d=0., j=0; j<nclass; j++) d += SQR(con[i][j]);
    ic += d / (double)at[i];
  }
  for (j=0; j<nclass; j++) ir += SQR(clt[j]);
  ir = ir / (double) nr;
/* printf("%lf %lf %lf %lf", ic, ir, (double)nr, (ic - ir) / (double) nr); */
  return (ic - ir) / (double) nr;
}

void informativity_init(fams *selected_fam)
{
  sf = selected_fam;
  nin = sf->n_in;
  nclass = sf->out->ndesc;
  compute_mclt_nr();		/* determine clt and its inform */
}

double compute_informativity(char *sel)
{
  compute_mat_mcon();
/*  printf("%lf %lf %lf %lf\n", -mcon/nr, mclt/nr, mat/nr, NLOGN(nr)/nr); */
  return ((-mcon + mclt + mat + NLOGN(nr))/ (double) nr);
}


double my_comb(int n, int k)
{
  double c1=1, c2=1;
  int i;
  for (i=0; i<k; i++) {
    c1 *= n-i; c2 *= k-i;
  }
  return c1 / c2;
}

double my_comb_log(int n, int k)
{
  double c1=0, c2=0;
  int i;
  for (i=0; i<k; i++) {
    c1 += LOG2(n-i); c2 += LOG2(i+1);
  }
  return c1 - c2;
}

double compute_mdl()
{
  double d;
  double mdl;
  int i;


  mdl = (( mclt + mat - mcon + NLOGN(nr))/ (double) nr);
  mdl += my_comb_log(nr + nclass - 1, nclass - 1) / (double) nr;

  for (d=0., i=0; i<naval; i++)
    d += my_comb_log(at[i] + nclass - 1, nclass - 1);
  d = d / (double) nr;

  mdl -= d;

  return mdl;
}

void informativity_one(fams *selected_fam)
{
  int i, j;
  char *sel;			/* which attribute to derive IM for */


  informativity_init(selected_fam);
  im_one = d_vector_ini(sf->n_in,0);
  w_one = d_vector_ini(sf->n_in,0);
  gr_one = d_vector_ini(sf->n_in,0);
  gini_one = d_vector_ini(sf->n_in,0);
  mdl_one  = d_vector_ini(sf->n_in,0);

/*  if (log_inform || log_gr) im_one = d_vector(sf->n_in);
  if (log_gr) gr_one = d_vector(sf->n_in);
  if (log_gini) gini_one = d_vector(sf->n_in);
*/

  sel = c_vector(sf->n_in);
  if (deb_dec > 2) {
    pspace(print_short-1);
    printf("              IM /     Gini /    GainR /     MDL /  Relieff\n");
  }

  if (log_relieff) 
    relieff(selected_fam);

  for (i=0; i<sf->n_in; i++) {
    for (j=0; j<nin; j++) sel[j]=FALSE;
    sel[i]=TRUE;
    build_contingency_table(sel);

    if (log_inform) im_one[i] = compute_informativity(sel);
    if (log_gr) {
      if (!log_inform) im_one[i] = compute_informativity(sel);
      gr_one[i] = im_one[i] / compute_iv();
    }
    if (log_gini) gini_one[i] = compute_gini();
    if (log_mdl) mdl_one[i] = compute_mdl();

    if (deb_dec > 2) 
      printf("IM(%s) = %8.5lf / %8.5lf / %8.5lf / %8.5lf/ %8.5lf\n", 
	   strn(sf->in[i]->name, print_short), 
	   im_one[i], gini_one[i], gr_one[i], mdl_one[i], w_one[i]);
  }
  free(sel);
}

void informativity_two(fams *f)
{
  char *sel;			/* which attribute to derive IM for */
  int i, j, k;

  informativity_init(f);
  if (log_inform) im_two = d_matrix(f->n_in, f->n_in);
  if (log_gr) gr_two = d_matrix(f->n_in, f->n_in);
  if (log_gini) gini_two = d_matrix(f->n_in, f->n_in);
  if (log_mdl) mdl_two = d_matrix(f->n_in, f->n_in);
  sel = c_vector(f->n_in);
  for (i=0; i<f->n_in-1; i++) 
    for (j=i+1; j<f->n_in; j++) {
      for (k=0; k<nin; k++) sel[k]=FALSE;
      sel[i]=TRUE; sel[j]=TRUE;
      build_contingency_table(sel);

      if (log_inform) im_two[i][j] = compute_informativity(sel);      
      if (log_gr) gr_two[i][j] = im_two[i][j] / compute_iv();
      if (log_gini) gini_two[i][j] = compute_gini();
      if (log_mdl) mdl_two[i][j] = compute_mdl();

      if (deb_dec > 2)  {
	printf("IM(%s,", strn(f->in[i]->name, print_short));
	printf("%s,%3d) = ", strn(f->in[j]->name, print_short), naval);

	printf("%6.3lf / %6.3lf / %6.3lf / %6.3lf (%6.3lf / %6.3lf / %6.3lf / %6.3lf)\n", 
	       im_two[i][j]/(im_one[i]+im_one[j]), 
	       gini_two[i][j]/(gini_one[i]+gini_one[j]),
	       gr_two[i][j]/(gr_one[i]+gr_one[j]), 
	       mdl_two[i][j]/(mdl_one[i]+mdl_one[j]),
	       im_two[i][j]-(im_one[i]+im_one[j]), 
	       gini_two[i][j]-(gini_one[i]+gini_one[j]),
	       gr_two[i][j]-(gr_one[i]+gr_one[j]), 
	       mdl_two[i][j]-(mdl_one[i]+mdl_one[j]));
      }
    }
  free(sel);
  if (deb_dec > 2) 
    for (i=0; i<f->n_in-1; i++) 
      for (j=i+1; j<f->n_in; j++) {
	printf("IM(%s,", strn(f->in[i]->name, print_short));
	printf("%s) = ", strn(f->in[j]->name, print_short));
	printf("%6.3lf / %6.3lf / %6.3lf / %6.3lf \n",
	       im_two[i][j], gini_two[i][j], gr_two[i][j], mdl_two[i][j]);
      }
}

/****************************************************************************
RELIEFF
Another measure to estimate the weight of the attributes on which the split
should be made. 
See Kononenko: Estimating attributes - analysis and extensions of RELIEF, 
Proc. ECML, April 1994, Catania, Italy.
****************************************************************************/

/****************************************************************************
New relieff
****************************************************************************/

double *w, *aw;
int natt;
int nins;
rule_list **rl;

int k_relieff = 5;

typedef struct {
  int indx;
  double dist;
} dist_item;

dist_item *indx;

double diff(int at, int a, int b)
{
  if (rl[a]->att[at] == rl[b]->att[at]) return 0.;
  else return 1.;
}

double xdist(int a, int b)
{
  double d=0;
  int i;

  for (i=0; i<natt; i++) { 
    d += diff(i,a,b);
  }
  return d;
}

static int mcmp(dist_item *a, dist_item *b)
{
  if (a->dist < b->dist) return -1;
  if (a->dist > b->dist) return 1;  
}

void find_diff(int ins, int class) 
{
  int i, a, k;
  
  for (a=0; a<natt; a++) aw[a] = 0;

  for (k=0, i=1; i<nins; i++) {
    if (rl[indx[i].indx]->class == class) {
      k++;
      for (a=0; a<natt; a++) {
	aw[a] += diff(a, indx[0].indx, indx[i].indx);
/*	printf("aa %lf\n", diff(a, indx[0].indx, indx[i].indx)); */
      }
      if (k==k_relieff) break;
    }
  }
/*  printf("k=%d %lf\n", k, aw[0]); */

/*  printf("xxx %d ", k);
  for (a=0; a<natt; a++) 
    printf("%9.5lf ", aw[a]); printf("\n"); */

  if (k>0) 
    for (a=0; a<natt; a++) 
      aw[a] = aw[a]/(double) k;    
}

int n_relieff = 200;

void relieff(fams *f)
{
  rule_list *rl1;
  int i, ins, a, c, k, j, *rank, *irank, vins;
  
  set_apriory(f);
  nins = count_rules(f);
  natt = f->n_in;
/*  printf("Relieff for %s\n", f->name); */

  rl = (rule_list **) malloc(nins * sizeof(*rl));
  for (i=0, rl1=f->lrule; rl1!=NULL; i++, rl1=rl1->next)
    rl[i] = rl1;

/*  printf("%d rules\n", nins); */

  w = (double *) malloc(sizeof(double) * natt);
  aw = (double *) malloc(sizeof(double) * natt);
  for (i=0; i<natt; i++) w[i] = 0;
  
  indx = (dist_item *) malloc(sizeof(*indx) * nins);

  for (vins=0; vins<n_relieff; vins++) {
    ins = (int) (rnd1e() * nins);
/*  for (ins=0; ins<nins; ins++) { */
    for (i=0; i<nins; i++) {
      indx[i].indx = i;
      indx[i].dist = xdist(ins, i);
    }
    qsort(indx, nins, sizeof(*indx), mcmp);
/*    for (i=0; i<nins; i++) {
      printf("%d %lf\n", indx[i].indx, indx[i].dist);
    } */

				/* hits */
    find_diff(ins, rl[ins]->class);

    for (a=0; a<natt; a++)
      w[a] -= aw[a];


				/* misses */
    for (c=0; c<nclass; c++) 
      if (c!=rl[ins]->class) {
	find_diff(ins, c);
	for (a=0; a<natt; a++) 
	  w[a] += f->out->c_apriory[c] * aw[a];
      }
  }

  for (a=0; a<natt; a++) 
    w[a] = w[a] / (double) nins;

/*  printf("\nweights\n    weight   rank\n"); */

/*  rank = (int*) malloc(sizeof(int)*natt);
  irank = (int*) malloc(sizeof(int)*natt);
  for (i=0; i<natt; i++) rank[i]=i;

  for (i=0; i<natt-1; i++)
    for (j=i+1; j<natt; j++) 
      if (w[rank[j]]>w[rank[i]]) { k=rank[i]; rank[i]=rank[j]; rank[j]=k; }
  for (i=0; i<natt; i++)
    irank[rank[i]] = i;
  
  for (a=0; a<natt; a++) 
    printf("%2d %9.5lf   %2d\n", a+1, w[a], irank[a]+1); */

  for (a=0; a<natt; a++) 
    w_one[a] = w[a];
}


/****************************************************************************
EVALUATION OF GOODNESS OF HEURISTICS

Given are two vectors, cr[] and h[]. cr reflects the real criteria
(the lower the better), h the heuristics (the higher the better).

How many splits, of which h was lower, was actually better to split?

   x = h[i], where i is such that x = max(h)
   mmax = 0
   foreach y in h, y=h[i]
     if (y==x)
       m = elements lower than cr[i]
       if m>mmax then mmax = m
   return n / |h|

Sort h. How many of top h we have to evaluate (split), before we come
to the optimal split (plus one)? 
Better, and this is what we use: how many splits of all splist should 
be performed?

   sort h, so that indx[i] reflects the previous position
   min = min(cr)
   i = 0
   while (cr[indx[i]] > min) { m++; i++; }
   j = i
   while (h[j] == h[i])
     if (cr[indx[j]] > cr[indx[i]]) m++;
     i++;
   return m+1

****************************************************************************/
extern FILE *log_file;

double compare_pos_heur_real(int *cr, double *h, int n)
{
  double max = 0.;
  int i, j, m;
  int mmax = 0;			/* max number of lower elements in cr */

/*for (i=0; i<n; i++) fprintf(log_file, "%d ", cr[i]); fprintf(log_file,"\n");
 for (i=0;i<n;i++) fprintf(log_file,"%5.3lf ", h[i]);fprintf(log_file,"\n"); */

  for (max=0, i=0; i<n; i++) if (h[i]>max) max=h[i];
  for (i=0; i<n; i++)
    if (h[i] == max) {
      for (m=0, j=0; j<n; j++)
	if (cr[j] < cr[i]) m++;
      if (m>mmax) mmax = m;
    }
  fprintf(log_file, "lower %d\n", mmax);
  return (double) mmax / (double) n * 100.0;
}

int size, *ncr;
double *nh;

void matrix_to_vect(int **cr, double **h, int n)
{
  int i,j,k;

  size = n * (n-1) / 2;
  ncr = i_vector(size);
  nh = d_vector(size);

  for (k=0, i=0; i<nin-1; i++)
    for (j=i+1; j<nin; j++) {
      ncr[k] = cr[i][j];
      nh[k] = h[i][j];
      k++;
    }
}

double bin_compare_pos_heur_real(int **cr, double **h, int n)
{
  double d;

  matrix_to_vect(cr, h, n);
  d = compare_pos_heur_real(ncr, nh, size);
  free(ncr); free(nh);
  return d;
}


double compare_splits_heur_real(int *cr, double *h, int n)
{
  int *indx;
  int i, j, k, min, m;

  indx = i_vector(n);
  for (i=0; i<n; i++) indx[i] = i;

  for (i=0; i<n-1; i++)		/* the most stupid sort */
    for (j=i+1; j<n; j++)
      if (h[indx[i]] < h[indx[j]]) {
	k = indx[i]; indx[i] = indx[j]; indx[j] = k;
      }
/*  for (i=0; i<n; i++) printf("%d ", indx[i]); printf("\n"); */

  for (min=cr[0], i=1; i<n; i++) if (min>cr[i]) min = cr[i];

  for (m=0; cr[indx[m]] > min; m++);
  j = i = m;
  i++;
  for (; i<n && h[indx[i]] == h[indx[j]]; i++)
    if (cr[indx[i]] < cr[indx[j]]) m++;
  fprintf(log_file, "splits %d\n", m);
  free(indx);
  return (m+1.) / (double) n * 100.;
}


double bin_compare_splits_heur_real(int **cr, double **h, int n)
{
  double d;

  matrix_to_vect(cr, h, n);  
  d = compare_splits_heur_real(ncr, nh, size);
  free(ncr); free(nh);
  return d;
}

/****************************************************************************
DERIVATION OF COMPLEXITY OF THE HIERARCHICAL STRUCTURE
****************************************************************************/

int get_dfc_node(var_type *v)
{
  int i, dfc;
  fams *sf = v->famsout;

  if (v->mark || sf==NULL) return 0;
  for (i=0, dfc=1; i<sf->n_in; i++) dfc *= sf->in[i]->ndesc;
  for (i=0; i<sf->n_in; i++)
    dfc += get_dfc_node(sf->in[i]);
  return dfc;
}

int get_dfc_measure(Str255 vname)
{
  var_type *v;

  FIND_VAR(v, vname);
  mark_var(v, FALSE);
  return get_dfc_node(v);
}


double get_code_node(var_type *v, var_type *root)
{
  int i;
  double dfc;
  fams *sf = v->famsout;

  if (v->mark || sf==NULL) return 0;
  for (i=0, dfc=1; i<sf->n_in; i++) dfc *= sf->in[i]->ndesc;
  dfc *= LOG2(v->ndesc);
  if (v!=root) dfc = dfc - LOG2((double)fact(v->ndesc));
  for (i=0; i<sf->n_in; i++)
    dfc += get_code_node(sf->in[i], root);
  return dfc;
}

double get_code_measure(Str255 vname)
{
  var_type *v;

  FIND_VAR(v, vname);
  mark_var(v, FALSE);
  return get_code_node(v, v);
}

extern double dtic(int c, int f);
extern double dticp(int c, int h);


double get_dtic_node(var_type *v, var_type *root)
{
  int i, j;
  double dfc;
  fams *sf = v->famsout;

  if (v->mark || sf==NULL) return 0;
  for (i=0, j=1; i<sf->n_in; i++) j *= sf->in[i]->ndesc;

  if (v!=root) 
    dfc = dticp(v->ndesc, j);
  else
    dfc = dtic(v->ndesc, j);
  /*  printf("xx for %s %lf\n", v->name, dfc); */

    
  for (i=0; i<sf->n_in; i++)
    dfc += get_dtic_node(sf->in[i], root);
  return dfc;
}

double get_dtic_measure(Str255 vname)
{
  var_type *v;

  FIND_VAR(v, vname);
  mark_var(v, FALSE);
  return get_dtic_node(v, v);
}
