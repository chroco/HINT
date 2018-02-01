/****************************************************************************
fuzzy.c

Fuzzy identification:
 - fuzzy clustering, Gustafson & Kessel
 - compatible cluster merging, Kaymak & Babuska
 - identification of Takagi & Sugeno model
 - translation into linguistic (Mamdani) model

The rutines are provided to do an evaluation (output variable
derivation) using derived TS model.
****************************************************************************/
/* fuzzy clustering */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "nrutil.h"
#define GL extern
#include "sds.h"
#include "pgapack.h"

var_type *fv;			/* variable used to identify */

int cl_N;			/* N data pairs */
int cl_d;			/* dimensionality of data,
				   i.e. number of attributes*/
int cl_NK;			/* tmp number of clusters */

				/* CLUSTERING */
double **cl_Z;			/* data points, N x (d+1) */
double **cl_NZ;			/* normalized data points */
double **cl_U;			/* fuzzy partition matrix, K x N */
double **cl_NU;			/* new U */
double **cl_V;			/* cluster prototypes, K x (d+1) */
double **cl_F;			/* cluster covariance matrix, (d+1)x(d+1) */
double **cl_M;			/* distance matrix, (d+1)x(d+1) */
double **cl_D;			/* real distance matrix, N x K */
int *cl_indx;			/* index vector, d+1 */
double *cl_data, *cl_datb;	/* tmp data, d+1 */
int cl_step;			/* iteration number */
double *cl_min, *cl_max;	/* min,max of attribute, (d+1) */

				/* CLUSTER MERGING */
double **cl_E;			/* eigenvectors, d+1 x d+1, IN COLUMN [][i] */
double **cl_ES;			/* smallest eigenvectors, [K][n+1] */
double **cl_C1, **cl_C2;	/* criteria matrices, K * K, symmetric */
double **cl_S;			/* similarity matrix, K * K */
double **cl_dis;		/* cluster distance in premise space */
int *cl_count;			/* number of clusters to merge */
int *cl_sort;			/* clst intx sorted by number of
				   clusters to merge */
double cl_aa, cl_bb;		/* parameters of mem func for cluster sym */

extern list_of_vars *lv_desc;	/* this stores a list of structure's leafs
				   that influence the variable to learn for */


#define NFORALL(indx,limit) for((indx)=1; (indx)<=(limit); (indx)++)

void cl_allocate()
{
				/* clustering */
  cl_Z = matrix(1, cl_N, 1, cl_d+1);
  cl_NZ = matrix(1, cl_N, 1, cl_d+1);
  cl_U = matrix(1, cl_K, 1, cl_N);
  cl_NU = matrix(1, cl_K, 1, cl_N);
  cl_V = matrix(1, cl_K, 1, cl_d+1);
  cl_F = matrix(1, cl_d+1, 1, cl_d+1);
  cl_M = matrix(1, cl_d+1, 1, cl_d+1);
  cl_D = matrix(1, cl_N, 1, cl_K);
  cl_data = dvector(1, cl_d+1);
  cl_datb = dvector(1, cl_d+1);
  cl_min = dvector(1, cl_d+1);
  cl_max = dvector(1, cl_d+1);
  cl_indx = ivector(1, cl_d+1);

				/* cluster merging */
  cl_E = matrix(1, cl_d+1, 1, cl_d+1);
  cl_ES = matrix(1, cl_K, 1, cl_d+1);
  cl_C1 = matrix(1, cl_K, 1, cl_K);
  cl_C2 = matrix(1, cl_K, 1, cl_K);
  cl_S = matrix(1, cl_K, 1, cl_K);
  cl_dis = matrix(1, cl_K, 1, cl_K);
  cl_count = ivector(1, cl_K);
  cl_sort = ivector(1, cl_K);
}

void cl_deallocate()
{
				/* clustering */
  free_matrix(1, cl_N, 1, cl_d+1);
  free_matrix(1, cl_N, 1, cl_d+1);
  free_matrix(1, cl_K, 1, cl_N);
  free_matrix(1, cl_K, 1, cl_N);
  free_matrix(1, cl_K, 1, cl_d+1);
  free_matrix(1, cl_d+1, 1, cl_d+1);
  free_matrix(1, cl_d+1, 1, cl_d+1);
  free_matrix(1, cl_N, 1, cl_K);
  free_dvector(1, cl_d+1);
  free_dvector(1, cl_d+1);
  free_dvector(1, cl_d+1);
  free_dvector(1, cl_d+1);
  free_ivector(1, cl_d+1);

				/* cluster merging */
  free_matrix(1, cl_d+1, 1, cl_d+1);
  free_matrix(1, cl_K, 1, cl_d+1);
  free_matrix(1, cl_K, 1, cl_K);
  free_matrix(1, cl_K, 1, cl_K);
  free_matrix(1, cl_K, 1, cl_K);
  free_matrix(1, cl_K, 1, cl_K);
  free_ivector(1, cl_K);
  free_ivector(1, cl_K);
}

/****************************************************************************
FUZZY CLUSTERING (Gustafson & Kessel)
****************************************************************************/

/* cl_compute_norm: computes a Frobenius norm of a matrix a with
   dimensions n x m */

double cl_compute_norm(double **a, int m, int n)
{
  int i, j;
  double norm = 0.;

  NFORALL(i, m)
    NFORALL(j, n)
      norm += SQR(a[i][j]);
  return sqrt(norm);
}

void cl_set_init_partition()
{
  int i, j, l;

  NFORALL(j, cl_N) {
    l = j % cl_K + 1;
    NFORALL(i, cl_K)
      if (i==l) cl_U[i][j] = 1.;
      else cl_U[i][j] = 0.;
  }    
}

void cl_set_data()
{
  int i, j;
  list_of_opt *opt;
  list_of_vars *l;

  for (cl_N=0, opt=options; opt!=NULL; cl_N++, opt=opt->next);
  for (cl_d=0, l=lv_desc; l->next!=NULL; cl_d++, l=l->next);
  printf("%d examples of dim %d\n", cl_N, cl_d);
  cl_deallocate();
  cl_allocate();

  for (i=1, opt=options; opt!=NULL; i++, opt=opt->next) {
    select_opt(opt->name);
    for (l=lv_desc, j=1; l->next!=NULL; l=l->next, j++) {
      cl_Z[i][j] = l->var->val;
    }
    cl_Z[i][j] = l->var->expect;
  }

  NFORALL(i, cl_N) {
    printf("%2d: ", i);
    NFORALL(j, cl_d+1)
      printf("%7.4lf ", cl_Z[i][j]);
    printf("\n");
  }

  /* normalization of data */
  NFORALL(j, cl_d+1) {
    cl_min[j] = cl_max[j] = cl_Z[1][j];
    for (i=2; i<=cl_N; i++) {
      if (cl_Z[i][j] > cl_max[j]) cl_max[j] = cl_Z[i][j];
      else if (cl_Z[i][j] < cl_min[j]) cl_min[j] = cl_Z[i][j];
    }
    NFORALL(i, cl_N)
      cl_NZ[i][j] = (cl_Z[i][j] - cl_min[j])/(cl_max[j]-cl_min[j]);
  }
}

void cl_print_data()
{
  int i, j;

  printf("partition matrix, data:\n");
  NFORALL(i, cl_N) {
    NFORALL(j, cl_K)
      printf("%8.5lf ", cl_U[j][i]);
    printf("     ");
    NFORALL(j, cl_d+1)
      printf("%8.5lf ", cl_Z[i][j]);
    printf("\n");
  }
}

void cl_compute_cluster_prototypes()
{
  int i, j, k;
  double mu;

  NFORALL(i, cl_K) {
    mu = 0;
    NFORALL(j, cl_N) mu += cl_U[i][j];
    NFORALL(k, cl_d+1) cl_V[i][k] = 0.;
    NFORALL(k, cl_d+1)
      NFORALL(j, cl_N)
	cl_V[i][k] += cl_U[i][j] * cl_NZ[j][k];
    NFORALL(k, cl_d+1) cl_V[i][k] = cl_V[i][k] / mu;
  }

  if (deb_cl_cp) {
    printf("cluster prototypes:\n");
    NFORALL(i, cl_K) {
      NFORALL(k, cl_d+1)
	printf("%8.5lf ", cl_V[i][k]);
      printf("\n");
    }
  }
}

void cl_compute_unscaled_cluster_prototypes()
{
  int i, j, k;
  double mu;

  NFORALL(i, cl_K) {
    mu = 0;
    NFORALL(j, cl_N) mu += cl_U[i][j];
    NFORALL(k, cl_d+1) cl_V[i][k] = 0.;
    NFORALL(k, cl_d+1)
      NFORALL(j, cl_N)
	cl_V[i][k] += cl_U[i][j] * cl_Z[j][k];
    NFORALL(k, cl_d+1) cl_V[i][k] = cl_V[i][k] / mu;
  }
}

void cl_compute_cluster_covariance_matrix(int i)
{
  int j, k, l;
  double mu;
  
  NFORALL(j, cl_d+1)
    NFORALL(k, cl_d+1)
      cl_F[j][k] = 0.;

  mu = 0;
  NFORALL(j, cl_N) mu += cl_U[i][j];

  NFORALL(l, cl_N) {
    NFORALL(j, cl_d+1)
      cl_data[j] = cl_NZ[l][j] - cl_V[i][j];

    NFORALL(j, cl_d+1)
      NFORALL(k, cl_d+1)
	cl_F[j][k] += pow(cl_U[i][l], cl_m) * cl_data[j] * cl_data[k];
  }
  
  NFORALL(j, cl_d+1)
    NFORALL(k, cl_d+1)
      cl_F[j][k] = cl_F[j][k] / mu;

  if (deb_cl_f) {
    printf("F for %d:\n", i);
    NFORALL(j, cl_d+1) {
      NFORALL(k, cl_d+1)
	printf("%8.5lf ", cl_F[j][k]);
      printf("\n");
    }
  }
}

double cl_distance(int j, int i)
{
  int k, l;
  double d;

  NFORALL(k, cl_d+1)
    cl_data[k] = cl_NZ[j][k] - cl_V[i][k];
  
  NFORALL(k, cl_d+1) cl_datb[k] = 0;
  NFORALL(k, cl_d+1)
    NFORALL(l, cl_d+1)
      cl_datb[k] += cl_M[k][l] * cl_data[l];
  
  d = 0.;
  NFORALL(k, cl_d+1)
    d += cl_data[k] * cl_datb[k];
  return d*d;
}

void cl_compute_distance_matrix(int i)
{
  double d;
  int j, k;

  ludcmp(cl_F, cl_d+1, cl_indx, &d);
  NFORALL(j, cl_d+1)
    d *= cl_F[j][j];
  if (deb_cl_f) printf("determinant %10.6e\n", d);
  if (d<0) d=1e-20;			/* for a slight num error */
  if (d==0.) {
    printf("error: F's determinant is 0\n");
    exit(0);
  }

  NFORALL(j, cl_d+1) {
    NFORALL(k, cl_d+1) cl_data[k] = 0.;
    cl_data[j] = 1.;
    lubksb(cl_F, cl_d+1, cl_indx, cl_data);
    NFORALL(k, cl_d+1) cl_M[k][j] = cl_data[k];
  }

/* printf("a %8.5e %lf\n",d, 1./(cl_d+1));*/
  d = pow(d, 1./(cl_d+1));
/* printf("b %lf\n", d); */
  NFORALL(j, cl_d+1)
    NFORALL(k, cl_d+1)
      cl_M[j][k] *= d;

  if (deb_cl_f) {
    printf("M for %d\n", i);
    NFORALL(j, cl_d+1) {
      NFORALL(k, cl_d+1)
	printf("%8.5lf ", cl_M[j][k]);
      printf("\n");
    }
  }

  NFORALL(j, cl_N)
    cl_D[j][i] = pow(cl_distance(j, i), -1./(cl_m-1));
}

void cl_update()
{
  int i, j;
  double d;
  
  NFORALL(i, cl_K) {
    cl_compute_cluster_covariance_matrix(i);
    cl_compute_distance_matrix(i);
  }

  NFORALL(j, cl_N) {
    d = 0.;
    NFORALL(i, cl_K) d += cl_D[j][i];

    NFORALL(i, cl_K)
      cl_NU[i][j] = cl_D[j][i] / d;
  }
}

void cl_get_clusters()
{
  int i, j, k;
  double **dd;
  double e = 1.;

  for (cl_step=0; e>cl_e && cl_step<100; cl_step++) {
    if (deb_cl_pm) cl_print_data();
    cl_compute_cluster_prototypes();
    cl_update();

    NFORALL(k, cl_K)
      NFORALL(j, cl_N)
	cl_U[k][j] -= cl_NU[k][j];

    e = cl_compute_norm(cl_U, cl_K, cl_N);
    printf("e[%d] %lf\n", cl_step, e);
    
    dd = cl_NU; cl_NU = cl_U; cl_U = dd;
  }
}

/****************************************************************************
COMPATIBLE CLUSTERS MERGING
****************************************************************************/

#define MF1(x) (exp(-SQR((x-1)/(0.4*(1-cl_aa)))))
#define MF2(x) (exp(-SQR(x/(0.4*cl_bb))))

int cl_merge_compatible_clusters()
{
  int i, j, k, l;
  double min, max, d;
  int imin;
  char b;

  /* derivation of eigenvectors and eigenvalues */
  NFORALL(i, cl_K) {

    cl_compute_cluster_covariance_matrix(i);
    jacobi(cl_F, cl_d+1, cl_data, cl_E, cl_indx);

    /* one has to remember only the smallest eigenvector */

/*    printf("eigenvectors for F%d (", i);
    NFORALL(j, cl_d+1) 
      printf("%8.5lf ", cl_data[j]);
    printf(" )\n");
    NFORALL(j, cl_d+1) {
      NFORALL(k, cl_d+1)
	printf("%8.5lf ", cl_E[k][j]);
      printf("\n");
    } */

    /* find smallest eigenvalue */
    min = cl_data[1]; imin = 1;
    for (j=2; j<=cl_d+1; j++)
      if (cl_data[j]<min) {min=cl_data[j]; imin=j;}
/*    printf("index of min %d\n", imin); */

    printf("%d min vect: ", i);
    NFORALL(k, cl_d+1) printf("%8.5lf ", cl_E[k][imin]);
    printf("\n");

    /* store the eigenvector to cl_ES */
    NFORALL(j, cl_d+1) cl_ES[i][j] = cl_E[j][imin];    
  }

  /* compute C1 and C2, both symmetrical, unidiagonal */
  for (i=1; i<=cl_K-1; i++)
    for (j=i+1; j<=cl_K; j++) {
      d = 0;
      NFORALL(k, cl_d+1)
	d += cl_ES[i][k] * cl_ES[j][k];
      cl_C1[i][j] = ABS(d);

      d = 0;
      NFORALL(k, cl_d+1)
	d += SQR(cl_V[i][k] - cl_V[j][k]);
      cl_C2[i][j] = sqrt(d);
    }

  printf("C1:\n");
  NFORALL(i, cl_K) {
    NFORALL(j, cl_K)
      printf("%5.3lf ", cl_C1[i][j]);
    printf("\n");
  }

  printf("C2:\n");
  NFORALL(i, cl_K) {
    NFORALL(j, cl_K)
      printf("%5.3lf ", cl_C2[i][j]);
    printf("\n");
  }

  /* compute aa, bb */
  cl_aa = cl_bb = 0.;
  for (i=1; i<=cl_K-1; i++)
    for (j=i+1; j<=cl_K; j++) {
      cl_aa += cl_C1[i][j];
      cl_bb += cl_C2[i][j];
    }
  cl_aa = 2 * cl_aa / (cl_K * (cl_K-1));
  cl_bb = 2 * cl_bb / (cl_K * (cl_K-1));

  printf("aa = %6.4lf, bb = %6.4lf\n", cl_aa, cl_bb);

  /* compute C1' and C2' */
  for (i=1; i<=cl_K-1; i++) {
    printf("%2d: ", i);
    for (j=i+1; j<=cl_K; j++) {
      cl_C1[i][j] = MF1(cl_C1[i][j]);
      cl_C2[i][j] = MF2(cl_C2[i][j]);
      cl_S[i][j] = cl_S[j][i] = sqrt(cl_C1[i][j] * cl_C2[i][j]);
      printf("%5.3lf ", cl_S[i][j]);
    }
    cl_S[i][i] = 1;
    printf("\n");
  }

  printf("C1:\n");
  NFORALL(i, cl_K) {
    NFORALL(j, cl_K)
      printf("%5.3lf ", cl_C1[i][j]);
    printf("\n");
  }

  printf("C2:\n");
  NFORALL(i, cl_K) {
    NFORALL(j, cl_K)
      printf("%5.3lf ", cl_C2[i][j]);
    printf("\n");
  }

  /* compute cl_dis */
  printf("distance matrix\n");
  for (i=1; i<=cl_K-1; i++) {
    printf("%2d: ", i);
    for (j=i+1; j<=cl_K; j++) {
      d = 0;
      NFORALL(k, cl_d)
	d += SQR(cl_V[i][k]-cl_V[j][k]);
      cl_dis[i][j] = sqrt(d);
      printf("%5.3lf ", cl_dis[i][j]);
    }
    printf("\n");
  }
  
  /* identify groups of clusters to be merged */
  printf("to merge\n");
  for (i=1; i<=cl_K-1; i++) {
    cl_count[i]=0;
				/* are there any clusters to merge? */
    b = FALSE;
    for (k=i+1; (k<=cl_K) && (!b); k++)
      if (cl_S[i][k] > cl_gamma)
	b = TRUE;

    if (b) {
      max = 0.;
      for (k=i; k<=cl_K-1; k++)
	if (cl_S[i][k] > cl_gamma)
	  for (l=k+1; l<=cl_K; l++)
	    if (cl_S[i][l] > cl_gamma)
	      if (max < cl_dis[k][l]) max = cl_dis[k][l];

      min = 1.;
      for (k=i; k<=cl_K-1; k++)
	if (cl_S[i][k] > cl_gamma)
	  for (l=k+1; l<=cl_K; l++)
	    if (cl_S[i][l] < cl_gamma)
	      if (min > cl_dis[k][l]) min = cl_dis[k][l];
    
      printf("%2d: %s (%5.3lf %5.3lf): ", i, min>max?"yes":"no", min, max);
      if (min>max) {
	for (k=i; k<=cl_K; k++)
	  if (cl_S[i][k] > cl_gamma) {
	    printf("%d ", k);
	    cl_count[i]++;
	  }
      }
      printf("\n");
    }
  }

  /* define which are biggest sets to merge */
  NFORALL(i, cl_K)
    cl_sort[i] = i;

  for (i=1; i<=cl_K-1; i++)
    for (j=i+1; j<=cl_K; j++)
      if (cl_count[cl_sort[i]] < cl_count[cl_sort[j]]) {
	k = cl_sort[i];
	cl_sort[i] = cl_sort[j];
	cl_sort[j] = k;
      }

  printf("sorted:\n");
  NFORALL(i, cl_K) printf("%d ", cl_sort[i]);
  printf("\n");

/*start from the biggest an merge, mark the cluster that is merged, do not
merge if the cluster is already marked */

  cl_NK = cl_K;
  for (l=1; l<cl_K; l++) {
    i = cl_sort[l];
    if (cl_count[i]>0) {
      printf("merging ");
      for (k=i+1; k<=cl_K; k++)
	if (cl_S[i][k] > cl_gamma) {
	  printf("%d ", k);
	  cl_count[k] = -1;
	  NFORALL(j, cl_N)
	    cl_U[i][j] += cl_U[k][j];
	  --cl_NK;
	}
      printf("\n");
    }
  }

  /* we now have to shift the data so to delete clusters with cl_count -1 */
  j = 2;
  for (i=2; i<=cl_K; i++) {
    for (; cl_count[i] == -1; i++);
    NFORALL(k, cl_N)
      cl_U[j][k] = cl_U[i][k];
    j++;
  }

  i = cl_K - cl_NK;
  cl_K = cl_NK;
  return i;

  /* merge, compute new cl_V, cl_F, etc... */
  /* return the number of merges */
}

/****************************************************************************
MEMBERSHIP FUNCTION IDENTIFICATION
****************************************************************************/

PGAContext *ctx;		/* current ga context */
double flower[4], fupper[4];
double a, b, c, d;		/* parameters */
int der_v, der_k;		/* var and cluster to derive */

double piece_sigm(double x, double a1, double s1, double a2, double s2)
{
  if (x<a1) return exp(-(SQR((x-a1)/(2.*s1))));
  if (x>a2) return exp(-(SQR((x-a2)/(2.*s2))));
  return 1.;

/*  double d;
  if (x<a1) d= exp(-(SQR((x-a1)/(2.*s1))));
  else if (x>a2) d= exp(-(SQR((x-a2)/(2.*s2))));
  else d = 1.;
  return d; */
}

double piece_lin(double x, double a, double b, double c, double d)
{
  if (x<a) return 0.;
  if (x<b) return (x-a)/(b-a);
  if (x<c) return 1.;
  if (x<d) return (d-x)/(d-c);
  return 0.;
}

double cl_evaluate_mem(PGAContext *ctx, int p, int pop) 
{
  double err, d;
  int i, j;

  PGAIndividual *ind;
  PGAReal      *chrom;

  ind = PGAGetIndividual ( ctx, p, pop );
  chrom = (PGAReal *)ind->chrom;

  if (cl_sigm) {
    if (chrom[0] > chrom[2]) {
      d=chrom[0]; chrom[0]=chrom[2]; chrom[2]=d;
    }
  }
  else {
    for (i=0; i<3; i++)
      for (j=i+1; j<4; j++)
	if (chrom[i]>chrom[j]) {
	  d=chrom[i]; chrom[i]=chrom[j]; chrom[j]=d;
	}
  }

  a = PGAGetRealAllele(ctx, p, pop, 0);
  b = PGAGetRealAllele(ctx, p, pop, 1);
  c = PGAGetRealAllele(ctx, p, pop, 2);
  d = PGAGetRealAllele(ctx, p, pop, 3);

/*  printf("eval %lf %lf %lf %lf  ", a, b, c, d); */
  
  err = 0.;
  NFORALL(i, cl_N)
    if (cl_sigm)
      err += SQR(piece_sigm(cl_Z[i][der_v], a, b, c, d) - cl_U[der_k][i]);
    else
      err += SQR(piece_lin(cl_Z[i][der_v], a, b, c, d) - cl_U[der_k][i]);
/*  printf("err %lf\n", sqrt(err)/cl_N); */
  return sqrt(err)/cl_N;
}

/* cl_derive_mem: derives parameters for a membership function for a
   given variable v and cluster k */

#define TSP(v,ti,tj,tk) v->ts.p[(ti)*v->ts.k*4 + (tj)*4 + tk]
#define TSA(v,ti,tj) v->ts.a[(ti)*(v->ts.d) + (tj)]

void cl_derive_mem()
{
  int i;
  int best_p;			/* best offspring */
  double err;

  printf("FOR var %d cluster %d\n", der_v, der_k);

  for (i=0; i<4; i++) {
    flower[i] = cl_min[der_v];
    fupper[i] = cl_max[der_v];
  }

  ctx = PGACreate(&g_argc, g_argv, PGA_DATATYPE_REAL, 4, PGA_MINIMIZE);
  PGASetRealInitLU(ctx, flower, fupper);
  PGASetRandomSeed(ctx, 100);
  PGASetMaxIter(ctx, 30);
  PGASetPopSize(ctx, 200);
  PGASetUp(ctx);

  PGARun(ctx, cl_evaluate_mem);

  best_p = PGAGetBest(ctx, PGA_OLDPOP);

  for (i=0; i<4; i++)
    TSP(fv, der_v-1, der_k-1, i)
      = PGAGetRealAllele(ctx, best_p, PGA_OLDPOP, i);

  a = PGAGetRealAllele(ctx, best_p, PGA_OLDPOP, 0);
  b = PGAGetRealAllele(ctx, best_p, PGA_OLDPOP, 1);
  c = PGAGetRealAllele(ctx, best_p, PGA_OLDPOP, 2);
  d = PGAGetRealAllele(ctx, best_p, PGA_OLDPOP, 3);

  err = PGAGetFitness(ctx, PGAGetBest(ctx, PGA_OLDPOP), PGA_OLDPOP);
  printf("final err %8.4g for %lf,%lf,%lf,%lf\n", err, a, b, c, d);

}

void cl_derive_membership_functions()
{
  int i;

  NFORALL(der_v, cl_d)
    NFORALL(der_k, cl_K)
      cl_derive_mem();
}

void cl_derive_ante()
{
  int i, j, k;
  double d;

  NFORALL(i, cl_K) {
    NFORALL(j, cl_d) 
      TSA(fv, i-1, j-1) = - cl_ES[i][j] / cl_ES[i][cl_d+1];
    d = 0;
    NFORALL(j, cl_d+1)
      d += cl_ES[i][j] * cl_V[i][j];
    fv->ts.b[i-1] = d / cl_ES[i][cl_d+1];

/*    printf("prototype: (%8.3lf)\n", fv->ts.b[i]);
    NFORALL(j, cl_d+1)
      printf("%8.3lf ", cl_V[i][j]);
    printf("\n"); */
  }
}

void cl_allocate_ts()
{
  int i, j;
  list_of_vars *l1;

  fv->ts.d = cl_d;
  fv->ts.var = (var_type **) malloc(sizeof(*fv->ts.var) * cl_d);
  fv->ts.k = cl_K;
  fv->ts.type = cl_sigm ? ts_sigm : ts_lin;
  
  fv->ts.p = (double *) malloc(sizeof(double) * cl_d * cl_K * 4);

  for (i=0, l1=lv_desc; l1->next!=NULL; i++, l1=l1->next)
    fv->ts.var[i] = l1->var;

  fv->ts.a = (double *) malloc(sizeof(double) * cl_K * cl_d);
  fv->ts.b = (double *) malloc(sizeof(double) * cl_K);
}

void cl_print_ts(var_type *v)
{
  int i, j, k;

  for (i=0; i<v->ts.k; i++) {
    printf("%-2d If  ", i);
    for (j=0; j<v->ts.d; j++) {
      printf("%s is (", v->ts.var[j]);
      for (k=0; k<4; k++)
	printf("%8.5lf ", TSP(v,j,i,k));
      printf(")\n");
      if (j+1<v->ts.d) printf("       ");
    }
    printf("   Then %s = ", v->name);
    for (j=0; j<v->ts.d; j++)
      printf("%8.3lf %s + ", TSA(v,i,j), v->ts.var[j]);
    printf("%8.3lf\n\n", v->ts.b[i]);
  }
}

/****************************************************************************
TS EVALUATE - usage of the identified model
****************************************************************************/

/* ts_evaluate_var: inference with TS model, see Babuska: Fuzzy set
   methods (chap 8), page 5 */

double ts_evaluate_var(var_type *v)
{
  int i, j;
  double beta, b;
  double sum_beta = 0., sum_y = 0.;
  double y;
  
  for (i=0; i<v->ts.k; i++) {

/*    printf("%d: ", i); */
    if (fa_method == famax) beta=0.;
    else beta=1.;
    for (j=0; j<v->ts.d; j++) {
      if (v->ts.type == ts_sigm)
	b = piece_sigm(v->ts.var[j]->val, 
		      TSP(v,j,i,0), TSP(v,j,i,1), TSP(v,j,i,2), TSP(v,j,i,3));
      else
	b = piece_lin(v->ts.var[j]->val, 
		      TSP(v,j,i,0), TSP(v,j,i,1), TSP(v,j,i,2), TSP(v,j,i,3));
      if (fa_method == famax) beta = MMAX(beta, b);
      else if (fa_method == famin) beta = MMIN(beta, b);
      else beta *= b;
    }
    sum_beta += b;
    
    for (y=0., j=0; j<v->ts.d; j++)
      y += TSA(v,i,j) * v->ts.var[j]->val;
    y += v->ts.b[i];

/*    printf("%6.3lf*%6.3lf ", b, y); */
    sum_y += b * y;
  }
  v->val = sum_y / sum_beta;
  return v->val;
}

void ts_check_data(var_type *v)
{
  int i;
  list_of_opt *opt;
  list_of_vars *l;

  for (opt=options; opt!=NULL; opt=opt->next) {
    select_opt(opt->name);
    for (i=0; i<v->ts.d; i++) 
      printf("%8.3lf ", v->ts.var[i]->val);
    printf("    %8.3lf   %8.3lf\n", v->expect, ts_evaluate_var(v));
  }
}

/****************************************************************************
MAIN FUNCTION
****************************************************************************/
void fuzzy_identification(Str255 vname)
{
  list_of_vars lv, *l1;
  int n_merged;
  int i;

  if (!((fv=find_var(variables, vname))!=NULL)) {
    printf("error: variable %s not found\n", vname);
    return;
  }
  lv.var = fv;
  lv.next = lv.prev = NULL;

  derive_v_with_desc(&lv);

  printf("clustering for %s(", fv->name);
  for (l1=lv_desc; l1->next!=NULL; l1=l1->next)
    printf("%s ", l1->var->name);
  printf(")\n");

  cl_set_data();

  for (i=0; i<cl_redo; i++) {
    cl_set_init_partition();
    do {
      cl_get_clusters();
      cl_print_data();
      n_merged = cl_merge_compatible_clusters();
      cl_print_data();
    } while (n_merged > 0);
  }

  /* estimate the membership functions */
  cl_allocate_ts();
  cl_derive_membership_functions();

  /* estimate the local linear models with weighted points */
  cl_compute_unscaled_cluster_prototypes();
  cl_derive_ante();
  cl_print_ts(fv);

  ts_check_data(fv);
}
