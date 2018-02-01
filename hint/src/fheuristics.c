/****************************************************************************
heuristics.c

Heuristic measures for decomposition or for estimation of godness of
attributes. The code includes:

  * contingency table derivation
  * derivation of information measure (from contingency table)
  * relieff estimation

OLD RUTINES WITH SPECIAL CASES FOR ONE AND TWO PARAMETERS, BUT WHEN
PERFORMANCE IS NEEDED THIS CAN BE MUCH FASTER

****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#define GL extern
#include "sds.h"

double *im_one, **im_two;	/* information measure for attributes
				   im_one - single attribute 
				   im_two - combination of two attributes */
double *w_one, **w_two;		/* relieff measure for attributes */
fams *sf;			/* table of rules for which to derive
				   different measures */
int nin;			/* number of input variables to that fam */

/****************************************************************************
INFORMATIVITY

The computation of the informativity of a single attribute and on a
combination of two attributes. Computation is based on condingency
tables.
See Mingles: Decision tree induction, Machine Learning 3: 319 -342
****************************************************************************/

int *clt, **con, *at;		/* components of contingency table */
int nclass, naval;		/* #examples, #classes, #attribute values */
double mclt, mcon, mat;		/* components in inf measure equation */

double clog2 = 0.69314718056;

#ifdef hp
#define NLOGN(x) (((x)==0)?0:((x) * log((double)x) / clog2))
#else
#define NLOGN(x) (((x)==0)?0:((x) * log2((double)x)))
#endif

void compute_mclt_nr()
{
  int i;

  clt = i_vector_ini(nclass, 0);
  for (i=0; i<sf->n_rules; i++) if (sf->rtp[i]!=undef) clt[sf->rule[i]]++;
  for (mclt=0., i=0; i<nclass; i++) mclt -= NLOGN(clt[i]);
/*  for (nr=0, j=0; j<sf->n_rules; j++) if (sf->rtp[j]!=undef) nr++; */
  free(clt);
}

void compute_mat_mcon() 
{
  int j, k;

  mat = mcon = 0.;
  for (j=0; j<naval; j++) {
    mat -= NLOGN(at[j]);
/*    printf("xxx %d %lf %lf\n", at[j], log((double)at[j]), NLOGN(at[j])); */
    for (k=0; k<nclass; k++) {
/*      printf("%d ", con[j][k]); */
      mcon += NLOGN(con[j][k]);
    }
/*    printf("\n"); */
  }
}

void informativity_one(fams *selected_fam)
{
  int i, j, nr;
  int *att;			/* values of the attributes */
  int *which;			/* which attribute to derive IM for,
				   has TRUE for that which[i] */

  sf = selected_fam;
  nin = sf->n_in;
  nclass = sf->out->ndesc;
  compute_mclt_nr();		/* determine clt and its inform */

  im_one = d_vector(sf->n_in);
  att = i_vector(sf->n_in);
  for (nr=0, j=0; j<sf->n_rules; j++) if (sf->rtp[j]!=undef) nr++;
  for (i=0; i<sf->n_in; i++) {
    naval = sf->in[i]->ndesc;
    at = i_vector_ini(naval, 0);
    con = i_matrix_ini(naval, nclass, 0);

    for (j=0; j<sf->n_rules; j++) {
      indx2att(sf, j, att);
      at[att[i]]++;
      if (sf->rtp[j]!=undef) (con[att[i]][sf->rule[j]])++;
    }
    compute_mat_mcon();

/*im_one[i] = (mcon + mclt + mat + NLOGN(sf->n_rules))/ (double) sf->n_rules;*/
    im_one[i] = (mcon + mclt + mat + NLOGN(nr))/ (double) nr;
/*    printf("IM(%s) = %6.3lf\n", sf->in[i], im_one[i]); */
    free(at); free(con);
  }
}

void informativity_two(fams *selected_fam)
{
  int i, j, k, nr;
  int *att;			/* values of the attributes */
  int indx;

  sf = selected_fam;
  nin = sf->n_in;
  nclass = sf->out->ndesc;
  compute_mclt_nr();		/* determine clt and its inform */

  for (nr=0, j=0; j<sf->n_rules; j++) if (sf->rtp[j]!=undef) nr++;
  im_two = d_matrix(sf->n_in, sf->n_in);
  att = i_vector(sf->n_in);
  for (i=0; i<sf->n_in-1; i++) 
    for (j=i+1; j<sf->n_in; j++) {
      naval = sf->in[i]->ndesc * sf->in[j]->ndesc;
      at = i_vector_ini(naval, 0);
      con = i_matrix_ini(naval, nclass, 0);

      for (k=0; k<sf->n_rules; k++) {
	indx2att(sf, k, att);
	indx = att[i] * sf->in[j]->ndesc + att[j];
/*	printf("%d %d: %d.%d->\n", att[i], att[j], indx, sf->rule[k]); */
	at[indx]++;
	if (sf->rtp[k]!=undef) (con[indx][sf->rule[k]])++;
      }
      compute_mat_mcon();
      
      im_two[i][j] = (mcon + mclt + mat + NLOGN(nr)) /
	(double) nr;
/*      printf("IM(%s,%s) = %6.3lf (mcon %6.3lf mclt %6.3lf mat %6.3lf N%d\n", 
	 sf->in[i], sf->in[j], im[i], mcon, mclt, mat, sf->n_rules); */
/*    printf("IM(%s,%s) = %6.3lf\n", sf->in[i], sf->in[j], im[i][j]); */
      free(at); free(con);
    }
}

/****************************************************************************
RELIEFF
Another measure to estimate the weight of the attributes on which the split
should be made. 
See Kononenko: Estimating attributes - analysis and extensions of RELIEF, 
Proc. ECML, April 1994, Catania, Italy.
****************************************************************************/

int *nc;			/* number of instances belonging to class */
double *pc;			/* prior probabilities of classes */
int *dist, *nn;			/* heaps for kNN, distances and inst. index */
int *att1, *att2;		/* values of attributes */
int kr;				/* adjusted number of k in NN */

int diff(int *a, int *b)
{
  int i, d=0;
  for (i=0; i<nin; i++) if (a[i]!=b[i]) d++;
  return d;
}

void print_knn(int class, char *s)
{
  int i, j;

  printf("%s\n", s);
  printf("class %d, kr %d, att ", class, kr);
  for (i=0; i<nin; i++) printf("%d ", att1[i]);
  printf("\n");
  for (i=0; i<kr; i++) {
    indx2att(sf, nn[i], att2);
    printf("dist  %d        nei ", dist[i]);
    for (j=0; j<nin; j++) printf("%d ", att2[j]);
    printf(" -> %d\n", sf->rule[nn[i]]);
  }
}

char find_knn(int indx, int class)
{
  int i, j, k, l, last, d;
  int maxdist = 9999;

  if (nc[class] < 2) return FALSE;
  if (krelieff > nc[class]/2) kr = nc[class]/2; else kr = krelieff;

				/* initialization */
  for (j=0, i=0; i<kr; j++)
    if (j != indx && sf->rule[j] == class) {
      indx2att(sf, j, att2);
      dist[i] = diff(att1, att2);
      nn[i++] = j;
    }
  last = j;
				/* sort */
  if (kr!=1)
    for (i=0; i<kr-1; i++)
      for (j=i+1; j<kr; j++)
	if (dist[i]>dist[j]) {
	  k = dist[i]; l = nn[i];
	  dist[i] = dist[j]; nn[i] = nn[j];
	  dist[j] = k; nn[j] = l;
	}

  for (i=last; i<sf->n_rules; i++) 
    if (i != indx && sf->rule[i] == class) {
      indx2att(sf, i, att2);
      d = diff(att1, att2);
      if (d<dist[kr-1]) {	/* insert */
	for (j=0; dist[j] <= d; j++);
	for (l=kr-1; l>j; l--) {
	  dist[l] = dist[l-1]; nn[l] = nn[l-1];
	}
	dist[j] = d; nn[j] = i;
      }
      else if (d=dist[kr-1]) {	/* replace */
	for (j=0; j<kr && dist[kr-1-j] == d; j++);
	k = (int)(rnd1() * (double)(j+1));
	if (k<j) {
	  dist[kr-1-k] = d; nn[kr-1-k] = i;
	}
      }
    }
/*  print_knn(1, "final"); */
  return TRUE;
}

double mean_diff_one(int at)
{
  int i;
  double d = 0.;

  for (i=0; i<kr; i++) {
    indx2att(sf, nn[i], att2);
    if (att1[at] != att2[at]) d+=1.;
  }
/*  printf("xx %lf \n", d/kr); */
  return d/kr;
}

double mean_diff_two(int at1, int at2)
{
  int i, i1, i2;
  double d = 0.;

  i1 = att1[at1] * sf->in[at2]->ndesc + att1[at2];
  for (i=0; i<kr; i++) {
    indx2att(sf, nn[i], att2);
    i2 = att2[at1] * sf->in[at2]->ndesc + att2[at2];
    if (i1 != i2) d+=1.;
  }
/*  printf("xx %lf \n", d/kr); */
  return d/kr;
}


double relieff(fams *selected_fam)
{
  int i, j, k, class;
  double d;
  double nr;

  sf = selected_fam;
  nin = sf->n_in;
  w_one = (double *) malloc(sizeof(double) * nin);  
  for (i=0; i<nin; i++) w_one[i]=i;

  nr = (double) sf->n_rules;
  w_one = d_vector_ini(nin, 0.); 
  w_two = d_matrix(nin, nin); 
  dist = i_vector(krelieff);
  nn = i_vector(krelieff);
  att1 = i_vector(nin);
  att2 = i_vector(nin);

  nc = i_vector_ini(sf->out->ndesc, 0);
  pc = d_vector(sf->out->ndesc); 
  for (i=0; i<sf->n_rules; i++) nc[sf->rule[i]]++;
  for (i=0; i<sf->out->ndesc; i++) pc[i] = (double) nc[i] / (double) sf->n_rules;

  for (i=0; i<sf->n_rules; i++) {
    indx2att(sf, i, att1);
    for (class=0; class<sf->out->ndesc; class++)
      if (find_knn(i, class)) {
	for (j=0; j<nin; j++) {
	  d = mean_diff_one(j);
	  if (class == sf->rule[i]) w_one[j] -= d;
	  else w_one[j] += pc[class] * d;
	}
	for (j=0; j<nin-1; j++) 
	  for (k=j+1; k<nin; k++) {
	    d = mean_diff_two(j, k);
	    if (class == sf->rule[i]) w_two[j][k] -= d;
	    else w_two[j][k] += pc[class] * d;
	  }
      }
  } 
  for (i=0; i<nin; i++) w_one[i] = w_one[i]/(double)nr;
  for (i=0; i<nin-1; i++) 
    for (j=i+1; j<nin; j++)
      w_two[i][j] = w_two[i][j]/(double)nr;
  free(nc); free(pc); free(att1); free(att2); 
}

