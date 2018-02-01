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

#define PLOT_POINTS 1000
double ga_error[PLOT_POINTS];	/* error for best offspring in each
				   generation observed */
var_type **allele_var;		/* info for which var is allele representing */

PGAContext *ctx;		/* current context */

/* myPGASetUp: sets the PGA parameters */

PGAContext *myPGASetUp(int lower[], int upper[], int string_len, void *f)
{
  ctx = PGACreate(&g_argc, g_argv, PGA_DATATYPE_INTEGER, string_len,
		  PGA_MINIMIZE);
  PGASetUserFunction(ctx, PGA_USERFUNCTION_MUTATION, (void*)f);
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
    err = PGAGetEvaluate(ctx, best_p, PGA_OLDPOP);
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

/* myPGARun: High-level routine to execute the genetic algorithm. Copy
   of PGARun, with some modifications. */

void myPGARun(double (*evaluate)(PGAContext *c, int p, int pop))
{
  int monc, best_p;
  double err;
  

  if (debug_g) printf("Running GA ..\n");
  PGADebugPrint(ctx, PGADEBUG_ENTERED, "PGARun", "Entered", 
		PGA_VOID, NULL);
  
  PGAEvaluate (ctx, PGA_OLDPOP, (*evaluate));
  myPGADebug();

  PGAFitness  (ctx, PGA_OLDPOP);
  monc = PGAGetMutateOnlyNoCross(ctx);
  while(!PGADone(ctx)) {
    PGASelect(ctx, PGA_OLDPOP);
    if (monc == PGATRUE) {
      PGARunMutateOrCross(ctx, PGA_OLDPOP, PGA_NEWPOP);
    }
    else
      PGARunMutateAndCross(ctx, PGA_OLDPOP, PGA_NEWPOP);
	
    PGAEvaluate(ctx, PGA_NEWPOP, (*evaluate));
    PGAFitness(ctx, PGA_NEWPOP);
    PGAUpdateGeneration(ctx);
/*    PGAPrintReport(ctx, stdout, PGA_OLDPOP); */
    myPGADebug();
  }
  PGADebugPrint(ctx, PGADEBUG_EXIT, "PGARun", "Exit", 
		PGA_VOID, NULL);

  err = (*evaluate)(ctx, PGAGetBest(ctx, PGA_OLDPOP), PGA_OLDPOP);
  printf("Best fit %3.8lf\n", err);

  PGADestroy              (ctx);
}

/* stat_for_var: derives an error of the estimation err, where err is
   sqrt of sum of squares of differences btw expected and derived
   values, where they exist */

double stat_for_vnum(int vnum)
{
  list_of_opt *opt;
  int n_opt = 0;
  double sum = 0;

  for (opt=options; opt!=NULL; opt=opt->next, ++n_opt)
    if (opt->valdef[vnum] && opt->expect_def[vnum]) {
      if (ga_error_method == gae_abs)
	sum += SQR(opt->val[vnum] - opt->expect[vnum]);
      else if (ga_error_method == gae_perc) {
	if (opt->expect[vnum] != 0.0)
	  sum += 100 * 
	    ABS((opt->val[vnum] - opt->expect[vnum])/opt->expect[vnum]);
	else
	  sum += 100000;	/* BBB this is stupid */
      }
    }
  if (ga_error_method == gae_abs)
    return sqrt(sum)/n_opt;
  else if (ga_error_method == gae_perc)
    return sum/n_opt;
}

void stat_for_var(Str255 vname)
{
  int vnum;

  if (find_var_num(variables, vname, &vnum)!=NULL)
    printf("Error: %10.5lf\n", stat_for_vnum(vnum));
  else
    printf("error: variable %s not found\n", vname);
}

/****************************************************************************
  LEARNING FAMS
  learn_fams: trims the rules of the fams that are used to derive a
  given variable so that the estimation error is minimized
****************************************************************************/

list_of_vars *lv_fams;		/* this stores a list of variables
				   that form the internal nodes in the
				   substructure rooting at a variable
				   for which we are interested */
var_type *v_to_learn;		/* variable for which fams to learn */
int vnum_to_learn;		/* its ID number */

int offset_fams = 0, offset_desc;	/* offsets in strings */
int offset_w = 0;
int len_fams, len_desc, len_w;		/* string lengths */


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

void derive_dependent_on(var_type *v)
{
  list_of_vars *l1, *l2;

  lv_fams = NULL;
  mark_var(v, FALSE);
  find_global_dependency(v);

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

double evaluate_fams(PGAContext *ctx, int p, int pop)
{
  double err;
  list_of_opt *opt;

  if (debug_g) printf("Evaluating .. %d\n", n_eval++);
  load_fams(ctx, p, pop);

  for (opt=options; opt!=NULL; opt=opt->next) {
    select_opt(opt->name);
    mark_var(lv_fams->var, FALSE);
    derive_var(lv_fams->var);

/*
printf("Opt %s:\n", opt->name);
list_struct(variables, curropt, NULL);
*/

    set_opt(opt->name);
  }

  err = stat_for_vnum(vnum_to_learn);
  if (debug_g) printf("-> %8.3lf\n",err);

  return err;
}

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
  for(i=offset_fams; i<(len_fams+offset_fams); i++)
    if (PGARandomFlip(ctx, mr)) { /* randomly choose an allele   */
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

void learn_fams(Str255 vname)
{
  var_type *v;
  double err;

  int *lower, *upper;

  if ((v=find_var_num(variables, vname, &vnum_to_learn))==NULL) {
    printf("error: variable %s not found\n", vname);
    return;
  }

  v_to_learn = v;
  derive_dependent_on(v);	/* derive a list of internal nodes */
  derive_string_len_fams(&len_fams); /* derives a string length for fams */

  printf("String length: %d\n", len_fams);

  lower = (int *) malloc(sizeof(int) * len_fams);
  upper = (int *) malloc(sizeof(int) * len_fams);
  allele_var = (var_type **) malloc(sizeof(*allele_var) * len_fams);
  set_lower_upper_fams(lower, upper);

  ctx = myPGASetUp(lower, upper, len_fams, (void*)mutate_fams);
  myPGARun(evaluate_fams);
  free_list_of_vars(lv_fams, FALSE);
}

/****************************************************************************
  LEARNING DESCRIPTIONS

  Learns fuzzy descriptions of the output and input variable. Output
  variable is the one for which learning was entered, and input
  variables are the corresponding leafs.
****************************************************************************/

list_of_vars *lv_desc;		/* this stores a list of structure's leafs
				   that influence the variable to learn for */
typedef enum {ca0, cb0, ca, cb, cc, cd, cbn, ccn} desc_chk_type;
int n_desc_v;			/* number of vars with desc to learn */
double *desc_hi, *desc_lo;	/* boundaries of numerical values of vars */
int *desc_min, *desc_max;	/* descriptors cannot cross this boundary */
				/* used withing GA to learn descriptors */
int *n_desc_pts;		/* landmark points for every variable */
desc_chk_type *chk_tp;		/* check type for points in membership func */

void find_leafs(var_type *v)
{
  int i;
  fams *f;
  list_of_vars *tmp;

  if (!v->mark) {
    v->mark = TRUE;
    if (v->famsout != NULL) {
      /*    printf("%s\n", v->name); */
      for (f=v->famsout; f!=NULL; f=f->next)
	for (i=0; i<f->n_in; i++)
	  find_leafs(f->in[i]);
    }
    else {
      tmp = (list_of_vars *) malloc(sizeof(*tmp));
      tmp->var = v;
      tmp->prev = lv_desc;
      lv_desc = tmp;
    }
  }
}

void derive_v_with_desc(var_type *v)
{
  list_of_vars *l1, *l2, *tmp;

				/* the top one is the output variable */
  tmp = (list_of_vars *) malloc(sizeof(*tmp));
  tmp->var = v;
  tmp->prev = NULL;
  lv_desc = tmp;

  mark_var(v, FALSE);
  find_leafs(v);

  n_desc_v = 0;
  l2 = NULL;			/* reverse the list */
  for (l1=lv_desc; l1!=NULL; l1=l1->prev, n_desc_v++) {
    l1->next = l2;
    l2 = l1;
  }
  lv_desc = l2;

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

char derive_string_len_desc(int *string_len)
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

  *string_len = 0;
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
/*	*string_len += 2; CCC */
	*string_len += 3;
	desc_hi[j] = v->desc[i].c;
	break;
      case left:
/*	*string_len += 2; CCC */
	*string_len += 3;
	desc_lo[j] = v->desc[i].a;
	break;
      case regular:
	*string_len += 4;
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
	  *string_len += 3;	/* CCC */
/*	  *string_len += 2; */
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
/* 	  *string_len += 2;*/
	  *string_len += 3;	/* CCC */
	  ahi = TRUE;
	}
	else {
	  v->desc[i].tp = regular;
	  *string_len += 4;
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
      desc_lo[j] -= 2*(dhi-dlo)/(v->ndesc-1);
      v->desc[0].a = desc_lo[j];
    }
    if (ahi) {
      desc_hi[j] += 2*(dhi-dlo)/(v->ndesc-1);
      v->desc[v->ndesc-1].c = desc_hi[j];
    }
    v->des_min = ga_desc_pts *
      (dlo-desc_lo[j])/(desc_hi[j]-desc_lo[j]);
    v->des_max = 1 + ga_desc_pts *
      (dhi-desc_lo[j])/(desc_hi[j]-desc_lo[j]);
    printf("0..%d: %d, %d\n", ga_desc_pts, v->des_min, v->des_max);

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
    /*if (debug_g)*/ printf("var %s\n", v->name);
    d = (double) ga_desc_pts / (double)(n_desc_pts[iv]-1);

    for (i=0; i<v->ndesc; i++) {
      switch (v->desc[i].tp) {
      case left:
/*	set_lohi_desc(lo, hi, &n, 0.5, 0.7, d, cb0, v);
	set_lohi_desc(lo, hi, &n, 0.7, 2.5, d, cc, v); */
	set_lohi_desc(lo, hi, &n, 0, 0.5, d, ca0, v); /* CCC */
	set_lohi_desc(lo, hi, &n, 0.5, 1.5, d, cb0, v);
	set_lohi_desc(lo, hi, &n, 1.5, 2.5, d, cc, v);
	break;
      case regular:
	for (j=i*2-1, k=0; j<i*2+3; j++, k++)
	  set_lohi_desc(lo, hi, &n, j-0.5, j+0.5, d, ca+k, v);
	break;
      case right:
/*	set_lohi_desc(lo,hi,&n,n_desc_pts[iv]-3.5, n_desc_pts[iv]-1.7,d,ca,v);
	set_lohi_desc(lo,hi,&n,n_desc_pts[iv]-1.7, n_desc_pts[iv]-1.5,d,cbn,v);
*/
	set_lohi_desc(lo,hi,&n,n_desc_pts[iv]-3.5, n_desc_pts[iv]-2.5,d,ca,v);
	set_lohi_desc(lo,hi,&n,n_desc_pts[iv]-2.5, n_desc_pts[iv]-1.5,d,cbn,v);
	set_lohi_desc(lo,hi,&n,n_desc_pts[iv]-1.5, n_desc_pts[iv]-1,d,ccn, v);
				/* CCC */
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
  int i;
  int count = 0;
  int min, max;

#ifndef OPTIMIZE
  PGADebugPrint(ctx, PGADEBUG_ENTERED, "PGAIntegerMutation", "Entered", 
		PGA_VOID, NULL);
#endif

  c = (PGAInteger *)PGAGetIndividual(ctx, p, pop)->chrom;
  for(i=offset_desc; i<(len_desc+offset_desc); i++)
    if (PGARandomFlip(ctx, mr)) { /* randomly choose an allele   */
				/* derivation of min & max */
      switch (chk_tp[i]) {
      case ca:
/*	min = c[i-2] + 1;
	max = c[i+1]; */
	if (chk_tp[i-3]!=ca0) min = MMAX(c[i-3], c[i-5]) + 1;
	else min = c[i-3] + 1;
	max = MMIN(c[i+1], c[i-1]);
	break;
      case cb:
/*	min = MMIN(c[i-2], c[i-1]) + 1;
	max = c[i+1]; */
	min = MMAX(c[i-1], c[i-3]) + 1;
	max = MMIN(c[i+1], c[i+3]);
	break;
      case cc:
/*	min = c[i-1] + 1;
	max = MMAX(c[i+1], c[i+2]); */
	min = MMAX(c[i-1], c[i-3]) + 1;
	max = MMIN(c[i+1], c[i+3]);
	break;
      case cd:
/*	min = c[i-1] + 1;
	max = c[i+2]; */
	min = MMAX(c[i-1], c[i+1]) + 1;
	if (chk_tp[i+3]!=ccn) max = MMIN(c[i+3],c[i+5]);
	else max = c[i+3];
	break;

      case cb0:
	min = c[i-1] + 1;
	max = MMIN(c[i+1],c[i+3]);
	break;
      case cbn:
	min = MMIN(c[i-1], c[i-3]) + 1;
	max = c[i+1];
      case ca0:			/* BBB this can be set at the begining */
	min = 0;
	max = allele_var[i]->des_min;
/* 	printf("xx %d %d %s\n", min, max, allele_var[i]->name); */
	break;
      case ccn:
	min = allele_var[i]->des_max;
	max = ga_desc_pts;
/*	printf("yy %d %d\n", min, max, allele_var[i]->name); */
	break;
      }

      if (min >= c[i])
	c[i] += ctx->ga.MutateIntegerVal;
      else  if (max <= c[i])
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

double evaluate_desc(PGAContext *ctx, int p, int pop)
{
  list_of_vars *l;
  int iv;			/* index of variable */
  double err;
  list_of_opt *opt;

  if (debug_g) printf("Evaluating .. %d\n", n_eval++);

  load_desc(ctx, p, pop);
  
/*  list_var_desc(variables); */
/*  plot_var_desc(variables); */

  for (opt=options; opt!=NULL; opt=opt->next) {
    select_opt(opt->name);
				/* for all leaves we have to recompute
                                   qualitative presentation */
    for (l=lv_desc, iv=0; l!=NULL; l=l->next, iv++)
      if (l->var->famsout == NULL)
	real_to_fuzzy(l->var);

    mark_var(lv_desc->var, FALSE);
    derive_var(lv_desc->var);
/*
printf("Opt %s:\n", opt->name);
list_struct(variables, curropt, NULL);
*/
    set_opt(opt->name);
  }

  err = stat_for_vnum(vnum_to_learn);
  if (debug_g) printf("-> %8.3lf\n",err);
  
  return err;
}



void learn_desc(Str255 vname)
{
  var_type *v;
  double err;

  int *lower, *upper;

  if ((v=find_var_num(variables, vname, &vnum_to_learn))==NULL) {
    printf("error: variable %s not found\n", vname);
    return;
  }

  v_to_learn = v;
  derive_v_with_desc(v);	/* derives v's with descriptions */
  
				/* derives a string length for desc */
  if (!derive_string_len_desc(&len_desc)) 
    return;			/* return if empty desc found */
  offset_desc = 0;
  chk_tp = (desc_chk_type *) malloc(sizeof(*chk_tp)*(len_desc));

				/* BBB we assume desc to be
                                   sorted. better - include a sorter
                                   and call it here */
  printf("String length: %d\n", len_desc);

  lower = (int *) malloc(sizeof(int) * len_desc);
  upper = (int *) malloc(sizeof(int) * len_desc);
  allele_var = (var_type **) malloc(sizeof(*allele_var) * len_desc);
  set_lower_upper_desc(lower, upper);

  ctx = myPGASetUp(lower, upper, len_desc, (void*)mutate_desc);
  myPGARun(evaluate_desc);
  free_list_of_vars(lv_desc, FALSE);
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

double evaluate_fams_desc(PGAContext *ctx, int p, int pop)
{
  double err;
  list_of_opt *opt;
  list_of_vars *l;
  int iv;			/* index of variable */

  if (debug_g) printf("Evaluating .. %d\n", n_eval++);
  load_fams(ctx, p, pop);
  load_desc(ctx, p, pop);
/*  plot_var_desc(variables);
  exit(0); */

  for (opt=options; opt!=NULL; opt=opt->next) {
    select_opt(opt->name);
    for (l=lv_desc, iv=0; l!=NULL; l=l->next, iv++)
      if (l->var->famsout == NULL)
	real_to_fuzzy(l->var);
    mark_var(lv_fams->var, FALSE);
    derive_var(lv_fams->var);

/*
printf("Opt %s:\n", opt->name);
list_struct(variables, curropt, NULL);
*/

    set_opt(opt->name);
  }

  err = stat_for_vnum(vnum_to_learn);
  if (debug_g) printf("-> %8.3lf\n",err);

  return err;
}

void learn_fams_desc(Str255 vname)
{
  var_type *v;
  int i;
  double err;
  int string_len; /* length of gene */

  int *lower, *upper;

  if ((v=find_var_num(variables, vname, &vnum_to_learn))==NULL) {
    printf("error: variable %s not found\n", vname);
    return;
  }

  v_to_learn = v;
  derive_v_with_desc(v);	/* derives v's with descriptions */
  derive_dependent_on(v);	/* derive a list of internal nodes */


/*  check_empty_desc(v);	*/	/* checks for the empty desc and */
				/* builds them, using the limits */
				/* derived from given options */
  derive_string_len_fams(&len_fams); /* derives a string length for fams */
				/* derives a string length for desc */
  if (!derive_string_len_desc(&len_desc)) 
    return;			/* return if empty desc found */
  string_len = len_fams + len_desc;
  printf("String len %d = %d (fams) + %d desc\n", 
	 string_len, len_fams, len_desc);
  offset_desc = len_fams;
  chk_tp = (desc_chk_type *) malloc(sizeof(*chk_tp)*(string_len));
				/* BBB we assume desc to be
                                   sorted. better - include a sorter
                                   and call it here */

  lower = (int *) malloc(sizeof(int) * string_len);
  upper = (int *) malloc(sizeof(int) * string_len);
  allele_var = (var_type **) malloc(sizeof(*allele_var) * string_len);
  set_lower_upper_fams(lower, upper);

/*  for (i=0; i<string_len; i++)
    printf("oo %i %d %d\n", i, lower[i], upper[i]); */

  set_lower_upper_desc(lower, upper);

  ctx = myPGASetUp(lower, upper, string_len, (void*)mutate_fams_desc);
  myPGARun(evaluate_fams_desc);
  free_list_of_vars(lv_desc, FALSE);
  free_list_of_vars(lv_fams, FALSE);
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
  for(i=offset_w; i<(len_w+offset_w); i++)
    if (PGARandomFlip(ctx, mr)) { /* randomly choose an allele   */
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
  int i, n=0;
				/* load weigths tables */
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

double evaluate_w(PGAContext *ctx, int p, int pop)
{
  double err;
  list_of_opt *opt;

  if (debug_g) printf("Evaluating .. %d\n", n_eval++);
  load_w(ctx, p, pop);

  for (opt=options; opt!=NULL; opt=opt->next) {
    select_opt(opt->name);
    mark_var(lv_fams->var, FALSE);
    derive_var(lv_fams->var);

/*
printf("Opt %s:\n", opt->name);
list_struct(variables, curropt, NULL);
*/

    set_opt(opt->name);
  }

  err = stat_for_vnum(vnum_to_learn);
  if (debug_g) printf("-> %8.3lf\n",err);

  return err;
}

void learn_w(Str255 vname)
{
  var_type *v;
  double err;
  int string_len; /* length of gene */

  int *lower, *upper;

  if ((v=find_var_num(variables, vname, &vnum_to_learn))==NULL) {
    printf("error: variable %s not found\n", vname);
    return;
  }

  v_to_learn = v;
  derive_dependent_on(v);	/* derive a list of internal nodes */
  derive_string_len_fams(&len_w); /* derives a string length for weights */
				/* this is same as for fams */
  string_len = len_w;
  printf("String len %d\n", len_w);

  lower = (int *) malloc(sizeof(int) * string_len);
  upper = (int *) malloc(sizeof(int) * string_len);
  allele_var = (var_type **) malloc(sizeof(*allele_var) * string_len);
  set_lower_upper_w(lower, upper);

  ctx = myPGASetUp(lower, upper, string_len, (void*)mutate_w);
  myPGARun(evaluate_w);
  free_list_of_vars(lv_fams, FALSE);
}

void ga_learn()
{
}
