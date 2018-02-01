/****************************************************************************
learn.c

Rutines that implement learning of rules. This, for now, uses genetic
algorithm optimization method.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define GL extern 
#include "sds.h"
#include "pgapack.h"

#define PUNISH_UNDEF 1e10	/* punishment for undefined */

#define PLOT_POINTS 1000
double ga_error[PLOT_POINTS];	/* error for best offspring in each
				   generation observed */
var_type **allele_var = NULL;	/* info for which var is allele representing */

PGAContext *ctx;		/* current context */

int offset_fams = 0;		/* offsets in strings */
int offset_desc;
int offset_w = 0;

int len_fams, len_desc, len_w;	/* string lengths */
int string_len;			/* length of gene */
int *lower, *upper;		/* boundaries for allele */

extern list_of_vars *lv_desc;	/* this stores a list of structure's
				   leafs that influence the variable
				   to learn for */
typedef enum {ca0, cb0, ca, cb, cc, cd, cbn, ccn} desc_chk_type;
int n_desc_v;			/* number of vars with desc to learn */
double *desc_hi, *desc_lo;	/* boundaries of numerical values of vars */
int *desc_min, *desc_max;	/* descriptors cannot cross this boundary */
				/* used withing GA to learn descriptors */
int *n_desc_pts;		/* landmark points for every variable */
desc_chk_type *chk_tp = NULL;	/* check type for points in membership func */

list_of_vars *lv_fams;		/* this stores a list of variables
				   that form the internal nodes in the
				   substructure rooting at a variable
				   for which we are interested */
list_of_vars *lv_to_learn;	/* variable for which fams to learn */
int n_vlearn;			/* number of them */
int *vnum_learn;		/* their position in opt values */

void ga_get_sites(PGAContext *ctx, int *xsite1, int *xsite2,
		  int from, int to)
{
  int temp;

  *xsite1 = PGARandomInterval(ctx, from, to);
  *xsite2 = *xsite1;
  while (*xsite2 == *xsite1)
    *xsite2 = PGARandomInterval(ctx, from, to);
  if (*xsite1 > *xsite2) {
    temp   = *xsite1;
    *xsite1 = *xsite2;
    *xsite2 = temp;
  }
}

void ga_crossover_twopoint(int xsite1, int xsite2,
		       PGAInteger *parent1, PGAInteger *parent2,
		       PGAInteger *child1, PGAInteger *child2,
		       int offset, int len)
{
  int i;

  for(i=offset; i<xsite1; i++) {
    child1[i] = parent1[i];
    child2[i] = parent2[i];
  }

  for(i=xsite1; i<xsite2; i++) {
    child1[i] = parent2[i];
    child2[i] = parent1[i];
  }

  for(i=xsite2; i<offset+len; i++) {
    child1[i] = parent1[i];
    child2[i] = parent2[i];
  }
}

void ga_crossover_fams(PGAContext *ctx,
		       PGAInteger *parent1, PGAInteger *parent2,
		       PGAInteger *child1, PGAInteger *child2)
{
  int xsite1, xsite2;

  ga_get_sites(ctx, &xsite1, &xsite2, offset_fams, offset_fams+len_fams);
  ga_crossover_twopoint(xsite1, xsite2, child1, child2, parent1, parent2,
			offset_fams, offset_fams+len_fams);
}

void ga_crossover_w(PGAContext *ctx,
		    PGAInteger *parent1, PGAInteger *parent2,
		    PGAInteger *child1, PGAInteger *child2)
{
  int xsite1, xsite2;

  ga_get_sites(ctx, &xsite1, &xsite2, offset_w, offset_w+len_w);
  ga_crossover_twopoint(xsite1, xsite2, child1, child2, parent1, parent2,
			offset_w, offset_w+len_w);
}

void ga_desc_minmax(PGAInteger *c, int i, int *min, int *max)
{
  switch (chk_tp[i]) {
  case ca:
    if (chk_tp[i-3]!=ca0) *min = MMAX(c[i-3], c[i-5]) + 1;
    else *min = c[i-3] + 1;
    *max = MMIN(c[i+1], c[i-1])-1;
    break;
  case cb:
    *min = MMAX(c[i-1], c[i-3]) + 1;
    *max = MMIN(c[i+1], c[i+3])-1;
    break;
  case cc:
    *min = MMAX(c[i-1], c[i-3]) + 1;
    *max = MMIN(c[i+1], c[i+3])-1;
    break;
  case cd:
    *min = MMAX(c[i-1], c[i+1]) + 1;
    if (chk_tp[i+3]!=ccn) *max = MMIN(c[i+3],c[i+5])-1;
    else *max = c[i+3]-1;
    break;
  case cb0:
    *min = c[i-1] + 1;
    *max = MMIN(c[i+1],c[i+3])-1;
    break;
  case cbn:
    *min = MMAX(c[i-1], c[i-3]) + 1;
    *max = c[i+1]-1;
    break;
  case ca0:			/* BBB this can be set at the begining */
    *min = 0;
    *max = MMIN(MMIN(c[i+1], c[i+3]), allele_var[i]->des_min)-1;
    /* 	printf("xx %d %d %s\n", *min, *max, allele_var[i]->name); */
    break;
  case ccn:
    *min = MMAX(MMAX(c[i-1], c[i-3]), allele_var[i]->des_max) + 1;
    *max = ga_desc_pts;
    /*	printf("yy %d %d\n", *min, *max, allele_var[i]->name); */
    break;
  }
  
  if (*min > *max) {
    printf("%d .. %d %d [%d]\n", *min, *max, chk_tp[i], i);
    exit(0);
  }
}

void ga_desc_change(PGAInteger *parent, PGAInteger *child, int i)
{
  int min, max;

  ga_desc_minmax(child, i, &min, &max);
  if (parent[i] > max)
    child[i] = max;
  else if (parent[i] < min)
    child[i] = min;
  else
    child[i] = parent[i];
/*  printf("%d .. %d -> %d\n", min, max, child[i]); */
}

void ga_crossover_desc(PGAContext *ctx,
		       PGAInteger *parent1, PGAInteger *parent2,
		       PGAInteger *child1, PGAInteger *child2)
{
  int i, temp, xsite1, xsite2;
  char b; int min, max;

  /* pick two cross sites such that xsite2 > xsite1 */
  
  if (debug_d)
    for (i=offset_desc; i<offset_desc+len_desc; i++) {
      ga_desc_minmax(parent1, i, &min, &max);
      if (parent1[i] > max || parent1[i] < min) {
	printf("parent1 error %d %d %d (%d)\n",
	       min, max, parent1[i], chk_tp[i]);
	for (i=offset_desc; i<offset_desc+len_desc; i++)
	  printf("%d ", parent1[i]);
	exit(0);
      }
    }
  
  if (debug_d)
  for (i=offset_desc; i<offset_desc+len_desc; i++) {
    ga_desc_minmax(parent2, i, &min, &max);
    if (parent2[i] > max || parent2[i] < min) {
      printf("parent2 error\n");
      printf("parent2 error %d %d %d (%d)\n", min, max, parent2[i], chk_tp[i]);
      exit(0);
    }
  }

  ga_get_sites(ctx, &xsite1, &xsite2, offset_desc, offset_desc+len_desc);

  for(i=offset_desc; i<offset_desc+len_desc; i++) {
    child1[i] = parent1[i];
    child2[i] = parent2[i];
  }

  for(i=xsite1; i<xsite2; i++) {
    ga_desc_change(parent2, child1, i);
    ga_desc_change(parent1, child2, i);
/*    child2[i] = parent2[i];
    child1[i] = parent1[i]; */
  }

  if (debug_d)
  for (i=offset_desc; i<offset_desc+len_desc; i++) {
    ga_desc_minmax(child1, i, &min, &max);
    if (child1[i] > max || child1[i] < min) {
      printf("error\n");
      exit(0);
    }
  }
}

void ga_crossover(PGAContext *ctx, int p1, int p2, int pop1,
		  int c1, int c2, int pop2)
{
  PGAInteger *parent1 = (PGAInteger *)PGAGetIndividual(ctx, p1, pop1)->chrom;
  PGAInteger *parent2 = (PGAInteger *)PGAGetIndividual(ctx, p2, pop1)->chrom;
  PGAInteger *child1  = (PGAInteger *)PGAGetIndividual(ctx, c1, pop2)->chrom;
  PGAInteger *child2  = (PGAInteger *)PGAGetIndividual(ctx, c2, pop2)->chrom;

/*  printf("x"); */
  if (ga_learn_fams) ga_crossover_fams(ctx, parent1, parent2, child1, child2);
/*  printf("in cross w\n"); */
  if (ga_learn_w) ga_crossover_w(ctx, parent1, parent2, child1, child2);
/*  printf("out cross\n"); */
  if (ga_learn_desc) ga_crossover_desc(ctx, parent1, parent2, child1, child2);

  return;
}

#define INTER(a) c[n++] = PGARandomInterval(ctx, inter[a]+1, inter[(a)+1]);

int inter[50];			/* BBB max n of descriptos */

void ga_create_desc_string(PGAContext *ctx, PGAInteger *c)
{
  list_of_vars *l;
  var_type *v;
  int i, j, k;
  int iv;			/* index of variable */
  int n;			/* index of string element */
  int min, max;

  int n_int;			/* number of intervals */
  int pts;

  n = offset_desc;
/*  if (debug_g) printf("Landmarks:\n"); */

  for (l=lv_desc, iv=0; l!=NULL; l=l->next, iv++) {
    v = l->var;

/*
    n_int = v->ndesc * 2;
    for (i=0; i<n_int; i++) inter[i]=2;
    if (n_int*2 > ga_desc_pts) {printf("create error\n"); exit(0);}

    for (i=0; i<ga_desc_pts-n_int*2; i++)
      inter[PGARandomInterval(ctx, 0, n_int-1)]++; */

/*    n_int = v->ndesc * 2;
    for (i=1; i<n_int-1; i++) inter[i]=2;
    inter[0] = v->des_min-1;
    inter[n_int-1] = ga_desc_pts - v->des_max;
    if (n_int*2 > ga_desc_pts) {printf("create error\n"); exit(0);}

    for (i=0; i<ga_desc_pts-((n_int-2)*2+v->des_min-1+
			     (ga_desc_pts-v->des_max)); i++)
      inter[PGARandomInterval(ctx, 1, n_int-2)]++; */

    n_int = v->ndesc * 2;
    for (i=1; i<n_int-1; i++) inter[i]=2;
    inter[0] = PGARandomInterval(ctx, 2, v->des_min-1);
    inter[n_int-1] = PGARandomInterval(ctx, 2, ga_desc_pts - v->des_max);
    if (n_int*2 > ga_desc_pts) {printf("create error\n"); exit(0);}

    for (pts=0, i=0; i<n_int; i++) pts+=inter[i];
    for (i=0; i<ga_desc_pts-pts; i++)
      inter[PGARandomInterval(ctx, 1, n_int-2)]++;

    for (i=1; i<n_int; i++)
      inter[i] += inter[i-1];

    if (debug_d) {
      for (i=0; i<n_int; i++)
	printf("%d ", inter[i]);
      printf("\n");
    }

    for (i=0; i<v->ndesc; i++) {
      switch (v->desc[i].tp) {
      case left:
	c[n++] = PGARandomInterval(ctx, 0, inter[0]);
	INTER(0);
	INTER(1);
	break;
      case regular:
	INTER((i-1)*2);
	INTER((i-1)*2+1);
	INTER((i-1)*2+2);
	INTER((i-1)*2+3);
	break;
      case right:
	INTER((i-1)*2);
	INTER((i-1)*2+1);
	INTER((i-1)*2+2);
	break;
      }      
    }
  }

  if (debug_d)
  for (i=offset_desc; i<offset_desc+len_desc; i++) {
    ga_desc_minmax(c, i, &min, &max);
    if (c[i] > max || c[i] < min) {
      printf("mico error %d..%d %d [%d] (%d)\n", min, max, c[i], i, chk_tp[i]);
      exit(0);
    }
  }
}

void ga_create_string(PGAContext *ctx, int p, int pop, int InitFlag)
{
  int i;
  PGAInteger *c;
  PGAIndividual *new = PGAGetIndividual(ctx, p, pop);
  
  new->chrom = (void *)malloc(ctx->ga.StringLen * sizeof(PGAInteger));
  if (new->chrom == NULL)
    PGAError(ctx, "PGAIntegerCreateString: No room to allocate "
	     "new->chrom", PGA_FATAL, PGA_VOID, NULL);
  c = (PGAInteger *)new->chrom;
  if (InitFlag) {
    for (i = 0; i < string_len; i++)
      c[i] = PGARandomInterval(ctx, ctx->init.IntegerMin[i],
			       ctx->init.IntegerMax[i]);
    ga_create_desc_string(ctx, c);
/*    (*ctx->ops.Randomize)(ctx, p, pop); */
  }
  else
    for (i=0; i<ctx->ga.StringLen; i++)
      c[i] = 0;
}

/* myPGASetUp: sets the PGA parameters */

PGAContext *myPGASetUp(int lower[], int upper[], int string_len, void *f)
{
  ctx = PGACreate(&g_argc, g_argv, PGA_DATATYPE_INTEGER, string_len,
		  PGA_MINIMIZE);
  PGASetUserFunction(ctx, PGA_USERFUNCTION_MUTATION, (void*)f);
  PGASetUserFunction(ctx, PGA_USERFUNCTION_CROSSOVER, (void*)ga_crossover);
  PGASetUserFunction(ctx, PGA_USERFUNCTION_CREATESTRING, (void*)ga_create_string);
  PGASetRandomSeed(ctx, 100);
  PGASetIntegerInitLU(ctx, lower, upper);
  PGASetMaxIter(ctx, ga_maxiter);
  PGASetCrossoverType(ctx, PGA_CROSSOVER_ONEPT);
  PGASetPopSize(ctx, ga_population_size);
  PGASetPrintFreq(ctx, 5);
  PGASetPrintOptions(ctx, PGA_STRING);
  PGASetPrintFreq(ctx, ga_print_freq);
  PGASetUp(ctx);
  return ctx;
}

/* myPGADebug: outputs the data about iteration number and
   best offspring in the generation */

int d_iter = 0;			/* iterations that where debugged */

void myPGADebug() 
{
  int best_p;
  double err;

  if (!(ctx->ga.iter % ctx->rep.PrintFreq) || ctx->ga.iter==0) {
    if (ctx->ga.iter == 0) d_iter = 0; else d_iter++;
    best_p = PGAGetBest(ctx, PGA_OLDPOP);
    err = PGAGetFitness(ctx, PGAGetBest(ctx, PGA_OLDPOP), PGA_OLDPOP);
    printf("Iter # %-5d - %8.4g\n", ctx->ga.iter, err);
    if (d_iter < PLOT_POINTS) ga_error[d_iter] = err;
  }
}

void plot_ga_error()
{
  int i;
  char fname[] = "000.plot";
  Str255 s;
  FILE *f;

  if (!d_iter) {
    printf("error: no errors to plot, use GA or set print freq\n");
    return;
  }

  f = fopen(fname, "w");
  if (f==NULL) {
    printf("error: can't open temporary %s file\n", fname);
    return;
  }

  for (i=0; i<d_iter && i<PLOT_POINTS; i++)
    fprintf(f, "%d %lf\n", i * ga_print_freq, ga_error[i]);
  fclose(f);
  
  sprintf(s, "gnugraph -L \"Fitness\" < %s | xplot -f 7x13", fname);
  system(s);
}

/****************************************************************************
  LEARNING FAMS
  learn_fams: trims the rules of the fams that are used to derive a
  given variable so that the estimation error is minimized
****************************************************************************/

void find_global_dependency(var_type *v)
{
  int i;
  fams *f;
  list_of_vars *tmp;

  if (!v->mark && v->famsout != NULL) {
/*    printf("%s\n", v->name); */

    tmp = (list_of_vars *) malloc(sizeof(*tmp));
    tmp->var = v;
    tmp->prev = lv_fams;
    lv_fams = tmp;

    v->mark = TRUE;
    for (f=v->famsout; f!=NULL; f=f->next)
      for (i=0; i<f->n_in; i++)
	find_global_dependency(f->in[i]);
  }
}

void derive_dependent_on(list_of_vars *lv)
{
  list_of_vars *l1, *l2, *tlv;

  lv_fams = NULL;
  mark_lvar(lv, FALSE);
  for (tlv=lv; tlv!=NULL; tlv=tlv->next)
    find_global_dependency(tlv->var);

  l2 = NULL;			/* reverse the list */
  for (l1=lv_fams; l1!=NULL; l1=l1->prev) {
    l1->next = l2;
    l2 = l1;
  }
  lv_fams = l2;

  if (debug_g) {
    printf("Dependency:\n");
    for (l1=lv_fams; l1!=NULL; l1=l1->next)
      printf("%s\n", l1->var->name);
  }
}

void derive_string_len_fams(int *string_len)
{
  list_of_vars *l;
  fams *f;
  int i;

  *string_len = 0;
  for (l=lv_fams; l!=NULL; l=l->next)
    for (f=l->var->famsout; f!=NULL; f=f->next)
      for (i=0; i<f->n_rules; i++)
	if ((f->rtp[i] != manfix) && (f->rtp[i] != autofix)) {
	  ++(*string_len);
	  f->rtp[i] = autom;
	}
}

void load_fams(PGAContext *ctx, int p, int pop)
{
  list_of_vars *l;
  fams *f;
  int i, n=0;

				/* load fam tables */
  for (l=lv_fams; l!=NULL; l=l->next)
    for (f=l->var->famsout; f!=NULL; f=f->next) {
/*      printf("fam %s for %s\n", f->name, l->var->name); */
      for (i=0; i<f->n_rules; i++)
	if ((f->rtp[i] != manfix) && (f->rtp[i] != autofix)) {
	  f->rule[i] = PGAGetIntegerAllele(ctx, p, pop, n);
/*	  printf("%d ", f->rule[i]); */
	  n++;
	}
    }
}

int n_eval = 0;			/* debugging purposes, eval number */

int mutate_fams(PGAContext *ctx, int p, int pop, double mr)
{
  PGAInteger *c;
  int i;
  int count = 0;

#ifndef OPTIMIZE
  PGADebugPrint(ctx, PGADEBUG_ENTERED, "PGAIntegerMutation", "Entered", 
		PGA_VOID, NULL);
#endif
  c = (PGAInteger *)PGAGetIndividual(ctx, p, pop)->chrom;
  mr = ga_mut_fams/(double)string_len; 
  for(i=offset_fams; i<(len_fams+offset_fams); i++)
    /* randomly choose an allele   */
    if (PGARandomFlip(ctx, mr)) { 
      if (PGAGetInitIntegerMin(ctx, i) == c[i])
	c[i] += ctx->ga.MutateIntegerVal;
      else  if (PGAGetInitIntegerMax(ctx, i) == c[i])
	c[i] -= ctx->ga.MutateIntegerVal;
      else {
	if (PGARandomFlip(ctx, .5)) /* add or subtract from allele */
	  c[i] += ctx->ga.MutateIntegerVal;
	else
	  c[i] -= ctx->ga.MutateIntegerVal;
      }
      count++;
    }
/*  if (count>0)   printf("f"); */
  return(count);
}

void set_lower_upper_fams(int lo[], int up[])
{
  list_of_vars *l;
  fams *f;
  int i, n=0;

  for (l=lv_fams; l!=NULL; l=l->next)

    for (f=l->var->famsout; f!=NULL; f=f->next)
      for (i=0; i<f->n_rules; i++)
	if ((f->rtp[i] != manfix) && (f->rtp[i] != autofix)) {
	  lo[n] = 0;
	  up[n] = l->var->ndesc - 1;
	  ++n;
	}
}

/****************************************************************************
  LEARNING DESCRIPTIONS

  Learns fuzzy descriptions of the output and input variable. Output
  variable is the one for which learning was entered, and input
  variables are the corresponding leafs.
****************************************************************************/

void derive_v_with_desc(list_of_vars *lv)
{
  var_type *v;
  list_of_vars *l1, *l2, *tmp;

				/* the top ones are the output variable */
  lv_desc = NULL;

  mark_lvar(lv, FALSE);
  find_leafs_lv(lv);

  l2 = lv;			/* reverse the list */
  for (l1=lv_desc; l1!=NULL; l1=l1->prev) {
    l1->next = l2;
    l2 = l1;
  }
  lv_desc = l2;
  for (n_desc_v=0,tmp=lv_desc; tmp!=NULL; tmp=tmp->next, n_desc_v++);

  if (debug_g) {
    printf("Inp/Out Nodes (%d):\n", n_desc_v);
    for (l1=lv_desc; l1!=NULL; l1=l1->next)
      printf("%s\n", l1->var->name);
  }
}

/* derive_string_len_desc: derives a string length for desc learning
   and checks if all variables from the list lv_desc have appropriate
   description to be used for GA learning. returns FALSE if at least
   one of the desc is empty. Also derives the boundaries of variable's
   numerical values */

char derive_string_len_desc(int *str_len)
{
  list_of_opt *opt;
  list_of_vars *tmp;
  var_type *v;
  int i, j, p;
  double d, dlo, dhi;		/* dummies */
  char alo, ahi;		/* automatic derivation of limits? */

				/* soft boundaries */
  desc_hi = (double *) malloc(sizeof(double) * n_desc_v);
  desc_lo = (double *) malloc(sizeof(double) * n_desc_v);
				/* hard boundaries */
  desc_min = (int *) malloc(sizeof(int) * n_desc_v);
  desc_max = (int *) malloc(sizeof(int) * n_desc_v);
  n_desc_pts = (int *) malloc(sizeof(int) * n_desc_v);

  *str_len = 0;
  for (tmp=lv_desc, j=0; tmp!=NULL; tmp=tmp->next, j++) {
    v = tmp->var;
    n_desc_pts[j] = 2*v->ndesc;
    if (!v->ndesc) {
      printf("error: descriptors for %s are not defined\n", v->name);
      return FALSE;
    }
    alo = ahi = FALSE;

    

    for (i=0; i<v->ndesc; i++)
      switch (v->desc[i].tp) {
      case right:
/*	*str_len += 2; CCC */
	*str_len += 3;
	desc_hi[j] = v->desc[i].c;
	break;
      case left:
/*	*str_len += 2; CCC */
	*str_len += 3;
	desc_lo[j] = v->desc[i].a;
	break;
      case regular:
	*str_len += 4;
	break;
      case none:
				/* this sets the q desc types, and figures
				   out the limits from the set of examples */
	if (i==0) {
	  v->desc[0].tp = left;
	  p = find_var_pos(variables, v->name);
	  opt=options;
	  desc_lo[j] = v->famsout ? opt->expect[p] : opt->val[p];
	  for (opt=opt->next; opt!=NULL; opt=opt->next) {
	    d = v->famsout ? opt->expect[p] : opt->val[p];
	    if (desc_lo[j] > d) desc_lo[j] = d;
	  }
	  v->desc[0].a = desc_lo[j]; /* this does not affect GA */
	  *str_len += 3;	/* CCC */
/*	  *str_len += 2; */
	  alo = TRUE;
	}
	else if (i==v->ndesc-1) {
	  v->desc[v->ndesc-1].tp = right;
	  p = find_var_pos(variables, v->name);
	  opt=options;
	  desc_hi[j] = v->famsout ? opt->expect[p] : opt->val[p];
	  for (opt=opt->next; opt!=NULL; opt=opt->next) {
	    d = v->famsout ? opt->expect[p] : opt->val[p];
	    if (desc_hi[j] < d) desc_hi[j] = d;
	  }
	  v->desc[v->ndesc-1].c = desc_hi[j]; /* this does not affect GA */
/* 	  *str_len += 2;*/
	  *str_len += 3;	/* CCC */
	  ahi = TRUE;
	}
	else {
	  v->desc[i].tp = regular;
	  *str_len += 4;
	}
/*	printf("error: membership function for %s of %s missing\n",
	       v->desc[i].name, v->name); */
/*	return FALSE; */
      }
    
    /* this is tricky. For the limits to be included and to even
       become possible in this sort of evaluation, we have to broaden
       to limits so that potentially mid of first and last descriptor
       would include lo and hi points in the middle */


    dlo = desc_lo[j];
    dhi = desc_hi[j];

    if (alo) {
      desc_lo[j] -= ga_width * (dhi-dlo)/(v->ndesc-1);
      v->desc[0].a = desc_lo[j];
    }
    if (ahi) {
      desc_hi[j] += ga_width * (dhi-dlo)/(v->ndesc-1);
      v->desc[v->ndesc-1].c = desc_hi[j];
    }
    v->des_min = ga_desc_pts *
      (dlo-desc_lo[j])/(desc_hi[j]-desc_lo[j]);
    v->des_max = 1 + ga_desc_pts *
      (dhi-desc_lo[j])/(desc_hi[j]-desc_lo[j]);
    if (debug_g) printf("0..%d: %d, %d\n", ga_desc_pts, v->des_min, v->des_max);

/*    printf("xx for %s, lo %10.5lf, hi %10.5lf ... +- %lf\n",
	   v->name, desc_lo[j], desc_hi[j], (dhi-dlo)/(v->ndesc-1)); */


  }
  return TRUE;
}

void set_lohi_desc(int lo[], int hi[], int *n, double dl, double dh, double d,
		   desc_chk_type tp, var_type *v)
{
  lo[*n] = dl * d + 1;		/* changed for safety */
   if (lo[3]!=0) {printf("here1 %d\n", *n); exit(0);} 
  hi[*n] = dh * d;
  allele_var[*n] = v;
  chk_tp[*n] = tp;
  if (debug_g) printf("  %d %d \t%d\n", lo[*n], hi[*n], chk_tp[*n]);
  (*n)++;
}

void set_lower_upper_desc(int lo[], int hi[])
{
  list_of_vars *l;
  var_type *v;
  int i, j, k;
  int iv;			/* index of variable */
  int n = 0;			/* index of string element */
  double d;

  n = offset_desc;
  if (debug_g) printf("Landmarks:\n");
  for (l=lv_desc, iv=0; l!=NULL; l=l->next, iv++) {
    v = l->var;
    if (debug_g) printf("var %s\n", v->name);
    d = (double) ga_desc_pts / (double)(n_desc_pts[iv]-1);

    for (i=0; i<v->ndesc; i++) {
      switch (v->desc[i].tp) {
      case left:
/*	set_lohi_desc(lo, hi, &n, 0.5, 0.7, d, cb0, v);
	set_lohi_desc(lo, hi, &n, 0.7, 2.5, d, cc, v); */
/*	set_lohi_desc(lo, hi, &n, 0, 0.1*n_desc_pts[iv], d, ca0, v);
	set_lohi_desc(lo, hi, &n, 0.01*n_desc_pts[iv], 0.15*n_desc_pts[iv], d, cb0, v);
	set_lohi_desc(lo, hi, &n, 0.15*n_desc_pts[iv], 0.20*n_desc_pts[iv], d, cc, v); */
        set_lohi_desc(lo, hi, &n, 0, 0.5, d, ca0, v); /* CCC */
        set_lohi_desc(lo, hi, &n, 0.5, 1.5, d, cb0, v);
        set_lohi_desc(lo, hi, &n, 1.5, 2.5, d, cd, v);
	break;
      case regular:
	for (j=i*2-1, k=0; j<i*2+3; j++, k++)
	  set_lohi_desc(lo, hi, &n, j-0.5, j+0.5, d, ca+k, v);
	break;
      case right:
/*	set_lohi_desc(lo,hi,&n,n_desc_pts[iv]-0.15*n_desc_pts[iv],
		      n_desc_pts[iv]-0.10*n_desc_pts[iv],d,ca,v);
	set_lohi_desc(lo,hi,&n,n_desc_pts[iv]-0.10*n_desc_pts[iv],
		      n_desc_pts[iv]-0.08*n_desc_pts[iv],d,cbn,v);
	set_lohi_desc(lo,hi,&n,n_desc_pts[iv]-0.08*n_desc_pts[iv],
		      n_desc_pts[iv]-0.01*n_desc_pts[iv],d,ccn, v); */

        set_lohi_desc(lo,hi,&n,n_desc_pts[iv]-3.5, n_desc_pts[iv]-2.5,d,ca,v);
        set_lohi_desc(lo,hi,&n,n_desc_pts[iv]-2.5, n_desc_pts[iv]-1.5,d,cbn,v);
        set_lohi_desc(lo,hi,&n,n_desc_pts[iv]-1.5, n_desc_pts[iv]-1,d,ccn, v);

	break;
      default:
	printf("something's wrong\n");
      }
    }
  }
}

int mutate_desc(PGAContext *ctx, int p, int pop, double mr)
{
  PGAInteger *c;
  int i,j;
  int count = 0;
  int min, max;

#ifndef OPTIMIZE
  PGADebugPrint(ctx, PGADEBUG_ENTERED, "PGAIntegerMutation", "Entered", 
		PGA_VOID, NULL);
#endif

  c = (PGAInteger *)PGAGetIndividual(ctx, p, pop)->chrom;

  if (debug_d)
  for (i=offset_desc; i<offset_desc+len_desc; i++) {
    ga_desc_minmax(c, i, &min, &max);
    if (c[i] > max || c[i] < min) {
      printf("\npre error %d %d %d (%d) [%d]\n", min, max, c[i], chk_tp[i], i);
      exit(0);
    }
  }


  mr = ga_mut_desc/(double)string_len; 
  for(i=offset_desc; i<(len_desc+offset_desc); i++)
    /* randomly choose an allele   */
    if (PGARandomFlip(ctx, mr)) { 
				/* derivation of min & max */
      ga_desc_minmax(c, i, &min, &max);

      if (min == max)
	c[i] = min;
      else if (min == c[i])
	c[i] += ctx->ga.MutateIntegerVal;
      else  if (max == c[i])
	c[i] -= ctx->ga.MutateIntegerVal;
      else {
	if (PGARandomFlip(ctx, .5)) /* add or subtract from allele */
	  c[i] += ctx->ga.MutateIntegerVal;
	else
	  c[i] -= ctx->ga.MutateIntegerVal;
      }
      count++;

      if (debug_d)
	for (j=MMAX(i-6,offset_desc);
	     j<(MMIN(i+6,offset_desc+len_desc-1)); j++) {
	  ga_desc_minmax(c, j, &min, &max);
	  if (c[j] > max || c[j] < min) {
	    printf("\nx error %d..%d %d (%d-%d) [%d-%d] %d\n",
		   min, max, c[j], chk_tp[j], chk_tp[i], j, i, count);
	    for (j=MMAX(i-6,offset_desc);
		 j<(MMIN(i+6,offset_desc+len_desc-1)); j++)
	      printf("%d:%d ", j, c[j]);
	  printf("\n");
	    exit(0);
	  }
	}
      
    }

  if (debug_d)
    for (i=offset_desc; i<offset_desc+len_desc; i++) {
      ga_desc_minmax(c, i, &min, &max);
      if (c[i] > max || c[i] < min) {
	printf("\nd error %d..%d %d (%d) [%d] %d\n",
	       min, max, c[i], chk_tp[i], i, count);
	exit(0);
      }
    }

/*  if (count > 0) printf("d"); */
  
  return(count);
}

void load_desc(PGAContext *ctx, int p, int pop)
{
  list_of_vars *l;
  var_type *v;
  int i;
  int iv;			/* index of variable */
  int n;			/* index of string element */
  double a, b, c, d, dd;

				/* load membership functions */
  n = offset_desc;
  for (l=lv_desc, iv=0; l!=NULL; l=l->next, iv++) {
    v = l->var;
    dd = (desc_hi[iv] - desc_lo[iv]) / ga_desc_pts;
    for (i=0; i<v->ndesc; i++) {
      switch (v->desc[i].tp) {
      case left:
/*	a = v->desc[i].a; */
				/* CCC */
	a = desc_lo[iv] + dd * PGAGetIntegerAllele(ctx, p, pop, n++);
	b = desc_lo[iv] + dd * PGAGetIntegerAllele(ctx, p, pop, n++);
	c = desc_lo[iv] + dd * PGAGetIntegerAllele(ctx, p, pop, n++);
	break;
      case regular:
	a = desc_lo[iv] + dd * PGAGetIntegerAllele(ctx, p, pop, n++);
	b = desc_lo[iv] + dd * PGAGetIntegerAllele(ctx, p, pop, n++);
	c = desc_lo[iv] + dd * PGAGetIntegerAllele(ctx, p, pop, n++);
	d = desc_lo[iv] + dd * PGAGetIntegerAllele(ctx, p, pop, n++);
	break;
      case right:
	a = desc_lo[iv] + dd * PGAGetIntegerAllele(ctx, p, pop, n++);
	b = desc_lo[iv] + dd * PGAGetIntegerAllele(ctx, p, pop, n++);
				/* CCC */
	c = desc_lo[iv] + dd * PGAGetIntegerAllele(ctx, p, pop, n++);
/*	c = v->desc[i].c; */
      }
      add_fdesc(&(v->desc[i]), v->desc[i].tp, a, b, c, d);
    }
				/* BBB could be done more efficiently */
  }
}

/****************************************************************************
  LEARNING FAMS AND DESCRIPTIONS

  Learns fams and fuzzy descriptions of the output and input
  variable. Output variable is the one for which learning was entered,
  and input variables are the corresponding leafs.

  A gene (string) in this section is thus combined of two parts, one
  representing fams and the other one representing a fuzzy
  desctiptions of a variables.
****************************************************************************/

int mutate_fams_desc(PGAContext *ctx, int p, int pop, double mr)
{
  int count;

  count = mutate_fams(ctx, p, pop, mr);
  count += mutate_desc(ctx, p, pop, mr);
  return count;
}

/****************************************************************************
  LEARNING RULE WEIGHTS
****************************************************************************/

void set_lower_upper_w(int lo[], int up[])
{
  list_of_vars *l;
  fams *f;
  int i, n=0;

  for (i=offset_w; i<(len_w+offset_w); i++) {
    lo[i] = 1;
    up[i] = ga_weight_pts;
  }
}

/* BBB this is actually the same as mutate for fams, this should be
used somehow. Note though that the offsets are different. */

int mutate_w(PGAContext *ctx, int p, int pop, double mr)
{
  PGAInteger *c;
  int i;
  int count = 0;

#ifndef OPTIMIZE
  PGADebugPrint(ctx, PGADEBUG_ENTERED, "PGAIntegerMutation", "Entered", 
		PGA_VOID, NULL);
#endif

  c = (PGAInteger *)PGAGetIndividual(ctx, p, pop)->chrom;
  mr = ga_mut_w/(double)string_len;
  for(i=offset_w; i<(len_w+offset_w); i++)
    /* randomly choose an allele   */
    if (PGARandomFlip(ctx, mr)) { 
      if (PGAGetInitIntegerMin(ctx, i) == c[i])
	c[i] += ctx->ga.MutateIntegerVal;
      else  if (PGAGetInitIntegerMax(ctx, i) == c[i])
	c[i] -= ctx->ga.MutateIntegerVal;
      else {
	if (PGARandomFlip(ctx, .5)) /* add or subtract from allele */
	  c[i] += ctx->ga.MutateIntegerVal;
	else
	  c[i] -= ctx->ga.MutateIntegerVal;
      }
      count++;
    }
  
  return(count);
}

void load_w(PGAContext *ctx, int p, int pop)
{
  list_of_vars *l;
  fams *f;
  int i, n;
				/* load weigths tables */
  n = offset_w;
  for (l=lv_fams; l!=NULL; l=l->next)
    for (f=l->var->famsout; f!=NULL; f=f->next) {
/*      printf("fam %s for %s\n", f->name, l->var->name); */
      for (i=0; i<f->n_rules; i++)
	if ((f->rtp[i] != manfix) && (f->rtp[i] != autofix)) {
	  f->imp[i] = (double) PGAGetIntegerAllele(ctx, p, pop, n) /
	    (double) ga_weight_pts;
/*	  printf("%4.2lf (%d) ", f->imp[i], PGAGetIntegerAllele(ctx, p, pop, n)); */
	  n++;
	}
/*      printf("\n"); */
    }
}

/****************************************************************************
****************************************************************************/


/* set_string_len: considers what we are learning and sets the string
   len accordingly */

void set_string_len()
{
  string_len = 0;
  if (ga_learn_fams) {
    offset_fams = string_len;
    string_len += len_fams;
  }
  if (ga_learn_desc) {
    offset_desc = string_len;
    string_len += len_desc;
  }
  if (ga_learn_w) {
    offset_w = string_len;
    len_w = len_fams;
    string_len += len_w;
  }
}
/* set_string_len: considers what we are learning and sets the
   boundaries accordingly */

void set_lower_upper()
{
  free(lower); free(upper);

  lower = (int *) malloc(sizeof(int) * string_len);
  upper = (int *) malloc(sizeof(int) * string_len);
  
  if (ga_learn_fams) set_lower_upper_fams(lower, upper);
  if (ga_learn_desc) set_lower_upper_desc(lower, upper);
  if (ga_learn_w) set_lower_upper_w(lower, upper);

/*  for (i=0; i<string_len; i++)
    printf("oo %i %d %d\n", i, lower[i], upper[i]); */
}

int ga_mutate(PGAContext *ctx, int p, int pop, double mr)
{
  int count;

  count = 0;
  if (ga_learn_fams) count += mutate_fams(ctx, p, pop, mr);
  if (ga_learn_desc) count += mutate_desc(ctx, p, pop, mr);
  if (ga_learn_w) count += mutate_w(ctx, p, pop, mr);
  return count;
}

double ga_evaluate(PGAContext *ctx, int p, int pop)
{
  double err;
  list_of_opt *opt;
  list_of_vars *l;
  int i;			/* index of variable */

  if (debug_g) printf("Evaluating .. %d\n", n_eval++);
  if (ga_learn_fams) load_fams(ctx, p, pop);
  if (ga_learn_desc) load_desc(ctx, p, pop);
  if (ga_learn_w) load_w(ctx, p, pop);

/*  plot_var_desc(variables);
  exit(0); */

  for (opt=options; opt!=NULL; opt=opt->next) {
    select_opt(opt->name);
    for (l=lv_desc; l!=NULL; l=l->next)
      if (l->var->famsout == NULL)
	real_to_fuzzy(l->var);
    mark_lvar(lv_fams, FALSE);
    derive_lvar(lv_fams);
/*
printf("Opt %s:\n", opt->name);
list_struct(variables, curropt, NULL);
*/
    set_opt(opt->name);
  }

  for (err=0, i=0; i<n_vlearn; i++)
    err += stat_for_vnum(vnum_learn[i]);
  if (debug_g) printf("-> %8.3lf\n",err);
  return err;
}

/* ga_mutate_cross: Performs crossover and mutation from one
   population to create the next.  Assumes PGASelect has been
   called. */

void ga_mutate_cross(PGAContext *ctx, int oldpop, int newpop)
{
  int i, j, n, m1, m2;
  int popsize, numreplace;
  double pc;

  PGADebugPrint(ctx, PGADEBUG_ENTERED, "PGARunMutateAndCross", "Entered", 
		PGA_VOID, NULL);
  popsize = PGAGetPopSize(ctx);
  numreplace = PGAGetNumReplace(ctx);
  /*** first, copy n best strings (sorted by fitness) to new pop ***/
  PGASortPop(ctx, oldpop);
  n = popsize - numreplace;
  for (i=0; i < n; i++) {
    j = PGAGetSortPop(ctx, i);
    PGACopyIndividual (ctx, j, oldpop, i, newpop);
  }
  pc = PGAGetCrossoverProb(ctx);
/*  pc = ga_cross_fams; */
  /*** reproduce to create the rest of the new population ***/
  while (n < popsize) {
    m1 = PGASelectNext(ctx);
    m2 = PGASelectNext(ctx);        
    if (PGARandomFlip(ctx, pc)) {
      PGACrossover (ctx, m1, m2, oldpop, PGA_TEMP1,
		    PGA_TEMP2, newpop);
             
      /*** mutate and copy first string to new population ***/
      PGAMutate (ctx, PGA_TEMP1, newpop);
      while (PGADuplicate(ctx, PGA_TEMP1, newpop, newpop, n))
	PGAChange (ctx, PGA_TEMP1, newpop);
      PGACopyIndividual (ctx, PGA_TEMP1, newpop, n, newpop);
      n++;

      if (n < popsize) {
	/*** mutate and copy second string to new population ***/
	PGAMutate (ctx, PGA_TEMP2, newpop);
	while (PGADuplicate(ctx, PGA_TEMP2, newpop, newpop, n))
	  PGAChange (ctx, PGA_TEMP2, newpop);
	PGACopyIndividual (ctx, PGA_TEMP2, newpop, n, newpop);
	n++;
      }
    }
    else {
      PGACopyIndividual (ctx, m1, oldpop, n, newpop);
      n++;
      if (n < ctx->ga.PopSize) {
	PGACopyIndividual (ctx, m2, oldpop, n, newpop);
	n++;
      }
    }
  }

  PGADebugPrint(ctx, PGADEBUG_EXIT, "PGARunMutateAndCross", "Exit", 
		PGA_VOID, NULL);
  return;
}


void ga_run()
{
  int monc, best_p;
  double err;
  
  if (debug_g) printf("Running GA ..\n");
  PGADebugPrint(ctx, PGADEBUG_ENTERED, "PGARun", "Entered", 
		PGA_VOID, NULL);
  
  PGAEvaluate (ctx, PGA_OLDPOP, ga_evaluate);
  myPGADebug();

  PGAFitness  (ctx, PGA_OLDPOP);
/*  monc = PGAGetMutateOnlyNoCross(ctx); */
  while(!PGADone(ctx)) {
    PGASelect(ctx, PGA_OLDPOP);
    switch (ga_policy) {
    case gap_or:
      PGARunMutateOrCross(ctx, PGA_OLDPOP, PGA_NEWPOP);
      break;
    case gap_and:
      PGARunMutateAndCross(ctx, PGA_OLDPOP, PGA_NEWPOP);
      break;
    case gap_custom:
      ga_mutate_cross(ctx, PGA_OLDPOP, PGA_NEWPOP);
      break;
    }
	
    PGAEvaluate(ctx, PGA_NEWPOP, ga_evaluate);
    PGAFitness(ctx, PGA_NEWPOP);
    PGAUpdateGeneration(ctx);
/*    PGAPrintReport(ctx, stdout, PGA_OLDPOP); */
    myPGADebug();
  }
  PGADebugPrint(ctx, PGADEBUG_EXIT, "PGARun", "Exit", 
		PGA_VOID, NULL);

  err = ga_evaluate(ctx, PGAGetBest(ctx, PGA_OLDPOP), PGA_OLDPOP);
  printf("Best fit %3.8lf %3.8lf\n", err, 
	 PGAGetFitness(ctx, PGAGetBest(ctx, PGA_OLDPOP), PGA_OLDPOP));

  PGADestroy              (ctx);
}

void ga_learn(list_of_vars *lv)
{
  list_of_vars *tlv;
  var_type *v;
  int i;
  double err;


  if (lv==NULL) return;
  for (n_vlearn=0, tlv=lv; tlv!=NULL; n_vlearn++, tlv=tlv->next);
  free(vnum_learn);
  vnum_learn = (int*)malloc(sizeof(int)*n_vlearn);
  for (i=0, tlv=lv; tlv!=NULL; i++, tlv=tlv->next)
    find_var_num(variables, tlv->var->name, &(vnum_learn[i]));

  lv_to_learn = lv;
  derive_v_with_desc(lv);	/* derives v's with descriptions */
  derive_dependent_on(lv);	/* derive a list of internal nodes */

  derive_string_len_fams(&len_fams); /* derives a string length for fams */
				/* derives a string length for desc */

  if (ga_learn_desc)
    if (!derive_string_len_desc(&len_desc)) return;

  set_string_len();

  printf("String len %d = %d (fams) + %d (desc) + %d (w)\n", 
	 string_len, len_fams, len_desc, len_w);

  free(chk_tp); free(allele_var);
  allele_var = (var_type **) malloc(sizeof(*allele_var) * string_len);
  chk_tp = (desc_chk_type *) malloc(sizeof(*chk_tp)*(string_len));

  /* BBB we assume desc to be sorted. better - include a sorter and
     call it here */

  set_lower_upper();

  ctx = myPGASetUp(lower, upper, string_len, (void*)ga_mutate);
/*  PGASetFitnessType(ctx, PGA_FITNESS_RANKING); */
  ga_run(ga_evaluate);
  free_list_of_vars(lv_desc, FALSE);
  free_list_of_vars(lv_fams, FALSE);
/*  if (debug_g) for (tlv=lv; tlv!=NULL; tlv=tlv->next)
    printf("to learn: %s\n", tlv->var->name); */
}

/****************************************************************************
Learning of the intervals
****************************************************************************/

double *imin, *imax;
list_of_itables *it;
var_type *v;
list_of_vars *ilv;
char cont_root;			/* is root a contninous variable? */
char add_len;

int mutate_int(PGAContext *ctx, int p, int pop, double mr)
{
  int is=0, i, k;
  list_of_vars *lv;
  PGAReal *c;
  int count = 0;
  double dd;

#ifndef OPTIMIZE
  PGADebugPrint(ctx, PGADEBUG_ENTERED, "PGAIntegerMutation", "Entered", 
		PGA_VOID, NULL);
#endif
  c = (PGAReal *)PGAGetIndividual(ctx, p, pop)->chrom;
  for (lv=ilv; lv!=NULL; k++, lv=lv->next) {
    for (i=0; i<lv->var->ndesc; i++)
      if (PGARandomFlip(ctx, mr)) { 
	c[i+is] = rnd1();
	count++;
      }
    is += lv->var->ndesc;
  }
  return(count);
}

double int_evaluate(PGAContext *ctx, int p, int pop)
{
  int is=0, i, k;
  double d;
  list_of_vars *lv;

				/* set the intervals */
  for (k=0, lv=ilv; lv!=NULL; k++, lv=lv->next) {
    for (d=0., i=0; i<lv->var->ndesc; i++) {
      d += PGAGetRealAllele(ctx, p, pop, is+i);
    }

    for (i=0; i<lv->var->ndesc; i++)
      PGASetRealAllele(ctx, p, pop, is+i,
		       PGAGetRealAllele(ctx, p, pop, is+i)/d);

    d = imin[k];
    for (i=0; i<lv->var->ndesc; i++) {
      lv->var->desc[i].start = d;
      lv->var->desc[i].delta = (imax[k]-imin[k]) * 
	PGAGetRealAllele(ctx, p, pop, is+i);
      d += lv->var->desc[i].delta;
    }
    lv->var->desc[i-1].delta *= 1.0001;

/*    for (i=0; i<lv->var->ndesc; i++)
      printf("%4.2lf..%4.2lf ", lv->var->desc[i].start,
	     lv->var->desc[i].start + lv->var->desc[i].delta);
    printf("\n"); */

    is += lv->var->ndesc;
  }
/*  for (i=0; i<string_len; i++)
    printf("%6.2g ", PGAGetRealAllele(ctx, p, pop, i));
  printf("\n"); */

  

				/* derive rules */
  sderive_rtable_from_itable(v->famsout, it);
  eval_instance_table(it);
  d=sstat_for_itable(it);
/*  printf("eval %lf\n", d); */
  return d;
}

int *s_ins;			/* index of sorted instances */
char *m_ins;			/* mark where the tval changes */
int nvals;			/* compare n values fro each instance */

int *i_class;			/* instance class distribution */
char set_rules;			/* does evaluation also need to set rules? */

/* rule_cmp_rows: compare two rules (indexes), so that to take the order
   for is_row. if r1 is smaller */

int ins_cmp(int *a, int *b)
{
  int i, j;
  int i1 = *a, i2 = *b;

  for (i=0; i<nvals; i++) 
    if (it->tval[i1][i] < it->tval[i2][i])
      return -1;
    else if (it->tval[i1][i] > it->tval[i2][i])
      return 1;
  return 0;
}

void my_list_ins()
{
  int i, j;

  for (i=0; i<it->n_inst; i++) {
/*    printf("%-3d: ", i);
    for (j=0; j<it->n_in; j++)
      print_ins_entry(it->in[j], it, s_ins[i], j);
    printf("\n"); */
    printf("%-3d:  %s ", i, m_ins[i]?"T":"F");
    for (j=0; j<it->n_in; j++) printf("%d ", it->tval[s_ins[i]][j]);
    printf(" -> %d\n", it->qval[s_ins[i]][j]);
  }
}

void set_ins_marks() 
{
  int i;
  for (i=0; i<it->n_inst-1; i++) {
    m_ins[i] = ins_cmp(&s_ins[i], &s_ins[i+1]) != 0;
  }
  m_ins[i] = TRUE;
}

double int_n_evaluate(PGAContext *ctx, int p, int pop)
{
  int is=0, i, j, k, n, e, nerror;
  double d;
  list_of_vars *lv;
  rule_type rtp = autom;

  if (set_rules) {
    del_rules(it->out->famsout);
  }

				/* set the intervals */
  for (k=0, lv=ilv; k<it->n_in; k++) 
    if (it->in[k]->ctype == ct_contin) {
    for (d=0., i=0; i<lv->var->ndesc; i++) {
      d += PGAGetRealAllele(ctx, p, pop, is+i);
    }

    for (i=0; i<lv->var->ndesc; i++)
      PGASetRealAllele(ctx, p, pop, is+i,
		       PGAGetRealAllele(ctx, p, pop, is+i)/d);

    d = imin[k];
    for (i=0; i<lv->var->ndesc; i++) {
      lv->var->desc[i].start = d;
      lv->var->desc[i].delta = (imax[k]-imin[k]) * 
	PGAGetRealAllele(ctx, p, pop, is+i);
      d += lv->var->desc[i].delta;
    }
    lv->var->desc[i-1].delta *= 1.0001;
    is += lv->var->ndesc;
    lv=lv->next;
  }

				/* derive rules */
  set_tval_ins(it);

  for (i=0; i<it->n_inst; i++) s_ins[i] = i;
  qsort(&s_ins[0], it->n_inst, sizeof(int *), ins_cmp);
  set_ins_marks();

  nerror = 0;
  for (i=0; i<it->n_inst; i++) {
    i_vector_set(i_class, it->out->ndesc, 0);
    n = 0;
    for (;;i++) {
      n++;
      ++i_class[it->qval[s_ins[i]][it->n_in]];
      if (m_ins[i]) break;
    }
    i_vector_max(i_class, it->out->ndesc, &k);
    e = n - i_class[k];
    if (set_rules) {
      set_rule(it->out->famsout, it->tval[s_ins[i]], k, rtp);
    }
    nerror += e;
  }
  return nerror;
}
  
void int_cr_str(PGAContext *ctx, int p, int pop, int InitFlag)
{
  int i;
  PGAReal *c;
  PGAIndividual *new = PGAGetIndividual(ctx, p, pop);
  
  new->chrom = (void *)malloc(ctx->ga.StringLen * sizeof(PGAReal));
  if (new->chrom == NULL)
    PGAError(ctx, "PGAIntegerCreateString: No room to allocate "
	     "new->chrom", PGA_FATAL, PGA_VOID, NULL);
  c = (PGAReal *)new->chrom;
  for (i = 0; i < string_len; i++)
    c[i] = rnd1();
}

void learn_intervals(Str255 vname, Str255 iname)
{
  list_of_vars *lv;
  int i, j;
  double err;
  
  if ((it = find_itable(iname)) == NULL) {
    printf("error: instance table %s not known\n", iname);
    return;
  }
  if ((v=find_var(variables, vname))==NULL) {
    printf("error: variable %s not found\n", vname);
    return;
  }

  cont_root = v->ctype != ct_nominal;
  add_len = cont_root?1:0;
  nvals = it->n_in+add_len;
  
  find_leafs(v);
  if (cont_root) {
    ilv = (list_of_vars *) malloc(sizeof(*ilv));
    ilv->var = v; 
    ilv->next = NULL;
  }
  else ilv = NULL;

  for (i=it->n_in-1; i>=0; i--) 
    if (it->in[i]->ctype == ct_contin) {
      lv = (list_of_vars *) malloc(sizeof(*lv));
      lv->var = it->in[i];
      lv->next = ilv; ilv = lv;
    }

  imin = d_vector(it->n_in+add_len);
  imax = d_vector(it->n_in+add_len);
  for (i=0; i<it->n_in+add_len; i++) imin[i] = imax[i] = it->val[0][i];
  for (j=1; j<it->n_inst; j++)
    for (i=0; i<it->n_in+add_len; i++) {
      if (imin[i] > it->val[j][i]) imin[i] = it->val[j][i];
      if (imax[i] < it->val[j][i]) imax[i] = it->val[j][i];
    }

  for (i=0,string_len=0, lv=ilv; i<it->n_in+add_len; i++) 
    if (it->in[i]->ctype == ct_contin) {
      printf("xx %s %lf %lf\n", lv->var->name, imin[i], imax[i]);
      string_len += lv->var->ndesc;
      lv=lv->next;
    }
  printf("string len %d\n", string_len);

  ctx = PGACreate(&g_argc,g_argv,PGA_DATATYPE_REAL,string_len,PGA_MINIMIZE);

  PGASetPopSize(ctx, ga_population_size);
  PGASetMaxIter(ctx, ga_maxiter);
  PGASetRandomSeed(ctx, 100);
/*  PGASetUserFunction(ctx, PGA_USERFUNCTION_MUTATION, (void*)mutate_int); */
/*  PGASetCrossoverType(ctx, PGA_CROSSOVER_ONEPT); */
  PGASetUserFunction(ctx, PGA_USERFUNCTION_CREATESTRING, (void*)int_cr_str);
  PGASetUp(ctx);

  if (v->ctype == ct_nominal) {
    s_ins = i_vector(it->n_inst);
    m_ins = c_vector(it->n_inst);
    i_class = c_vector(it->out->ndesc);
    set_rules = FALSE;
    PGARun(ctx, int_n_evaluate);

    set_rules = TRUE;
    err = int_n_evaluate(ctx, PGAGetBest(ctx, PGA_OLDPOP), PGA_OLDPOP);

    FREE(s_ins); FREE(m_ins); FREE(i_class);
  }
  else {
    PGARun(ctx, int_evaluate);
    err = int_evaluate(ctx, PGAGetBest(ctx, PGA_OLDPOP), PGA_OLDPOP);
  }

  PGADestroy(ctx); 
}
