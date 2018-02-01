/****************************************************************************
eval.c

Everything about the evaluation and derivation of values of variables.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define GL extern
#include "sds.h"

/* trapz: takes a real value x as an input, and for a description trz
   derives a degree of membership */

double trapz(double x, desc_t *trz)
{
  switch (trz->tp) {
  case left:
    if (x <= trz->b)
      return 1.0;
    if (x >= trz->c)
      return 0.0;
    /* b < x < c */
    return trz->r_slope * (x - trz->c);
  case right:
    if (x <= trz->a)
      return 0.0;
    if (x >= trz->b)
      return 1.0;
    /* a < x < b */
    return trz->l_slope * (x - trz->a);
  case regular:
    if ((x <= trz->a) || (x >= trz->d))
      return 0.0;
    if ((x >= trz->b) && (x <= trz->c))
      return 1.0;
    if ((x >= trz->a) && (x <= trz->b))
      return trz->l_slope * (x - trz->a);
    if ((x >= trz->c) && (x <= trz->d))
      return  trz->r_slope * (x - trz->d);
  }
  return 0.0;  /* should not get to this point */
}

void real_to_fuzzy(var_type *v)
{
  int i;

  for (i=0; i<v->ndesc; i++) {
    v->opt[i] = trapz(v->val, &(v->desc[i]));
/*    printf("xx var %s, %s, for %lf %lf\n", v->name, v->desc[i].name, v->val, v->opt[i]); */
  }
}

void real_to_interval(var_type *v)
{
  int i;

  for (i=0; i<v->ndesc; i++) v->opt[i]=0.;
  for (i=0; i<v->ndesc; i++) 
    if (v->val >= v->desc[i].start &&
	v->val < v->desc[i].start + v->desc[i].delta) {
      v->opt[i] = 1.;
      return;
    }
}


/* var_defined: returns true if the variable is defined, i.e., at
   least one of the descriptors has nonzero membership. else, returns
   false */

char check_var_defined(var_type *v)
{
  char b;
  int i;

  b = TRUE;
  for (i=0; i<v->ndesc && b; i++)
    if (v->opt[i] > 0)
      b = FALSE;
  return (!b);	/* BBB this is stupid */
}

void fuzzy_to_real(var_type *v)
{
  double w1=0, w2=0;
  int i;

  w1=0; w2=0;
  if (!check_var_defined(v)) {
    v->valdef = FALSE;
    if (debug_e) printf("xx undefined\n");
    return;
  }
  for (i=0; i<v->ndesc; i++) {
    if (v->desc[i].tp == none) {
      v->valdef = FALSE;
      return;
    }
    w1 += v->opt[i] * v->desc[i].cent;
    if (debug_e) printf("cent %lf %lf\n", v->opt[i], v->desc[i].cent);
    w2 += v->opt[i];
  }
  v->val = w1 / w2;
  if (debug_e) printf("%lf %lf\n", w1, w2);
  v->valdef = TRUE;    
}


/****************************************************************************
Dealing with fuzzy rules
****************************************************************************/

/* evaluate: takes a variable, figures out if everything is set,
   determines its values for all the FAMs that have this variable as
   output, and the combines the values to derive a final value of the
   variable */

char default_used;

void eval_fam(fams *f)
{
  var_type *v;
  int i, j, n;
  char *att, c;
  double min, w_sum;
  char b, b1, use_hash = FALSE;
  rule_list *rl;

  if (f->n_in>100) {use_hash=TRUE; printf("use hash\n"); build_hash_table(f);}

  v = f->out;
  if (debug_e) printf("eval %s, fam %s\n", v->name, f->name);

  for (i=0; i<v->ndesc; i++) f->degrees[i] = 0.;
  w_sum = 0.;
  att = c_vector(f->n_in);

  if (eval_method == e_crisp || eval_method == e_interval) {
    if (!use_distribution) { 
      for (i=0; i<f->n_in; i++) att[i] = f->in[i]->class;
    } 
    else {
      for (b1=TRUE, i=0; i<f->n_in && b1; i++) {
	for (b=FALSE, j=0; j<f->in[i]->ndesc && !b; j++)
	  b = f->in[i]->opt[j] == 1.;
	b1 = b;
	att[i] = --j;
      }
      if (!b1) {
	printf("distribution\n");
	exit(0);
      }
    }

    for (i=0; i<f->out->ndesc; i++) f->degrees[i]=0;

    c = get_rule(f, att);
    if (c != CUNDEF) f->degrees[j=c] = 1;
    else {
      f->degrees[f->out->c_default] = 1; 
      default_used = TRUE;
    }
    
    if (debug_e) {
      for (i=0; i<f->n_in; i++) printf("%s ", f->in[i]->desc[att[i]].name);
      printf(" -> %s (%d)\n", c!=CUNDEF?f->out->desc[j].name:"?", j);
    }      
  }

  else if (eval_method == e_fuzzy) {
    for (i=0; i < f->n_rules; i++) {
      if (debug_e) printf("rule %2d (%5.2lf): ", i, f->imp[i]); 
      if (f->rtp[i] != undef) {
	if (fa_method == famin || fa_method == famult) min = 1.;
	else if (fa_method == famax) min = 0.;
	
	indx2att(f, i, att);
	
	for (n=1, j=0; j < f->n_in; n*=f->in[j]->ndesc, j++) {
	  if (fa_method == famin) min = MMIN(min, f->in[j]->opt[att[j]]);
	  else if (fa_method == famax)
	    min = MMAX(min, f->in[j]->opt[att[j]]);
	  else if (fa_method == famult) {
	    min *= (min, f->in[j]->opt[att[j]]);
	  }
	  if (debug_e) printf(" (%d, %5.2lf)", att[j], f->in[j]->opt[att[j]]); 
	}
	
	min *= f->imp[i];
	(f->degrees[f->rule[i]]) += min;
	if (min>0) f->usage[i]++;
	w_sum += min;
	if (debug_e) printf(" -> (%d, %5.2lf)", f->rule[i], min);
      }
      if (debug_e)printf("\n");
    }
    
    if (w_sum != 0.0)
      for (i=0; i<v->ndesc; i++)
	f->degrees[i] /= w_sum;
    else {
      if (debug_e) printf("error: all weights 0, no rule matched\n");
      for (i=0; i<v->ndesc; i++)
	f->degrees[i] = 0;
    }
    
    if (debug_e) {
      printf("end of evaluation: ");
      for (i=0; i<v->ndesc; i++)
	if (f->degrees[i] > 0.0)
	  printf("%s/%6.3lf\t", v->desc[i].name, f->degrees[i]);
      printf("\n");
    }
  }

  free(att);
}

void eval_var(var_type *v)
{
  fams *f;
  int i;
  double w_sum = 0., d;
  char b;			/* checks if variables are defined */

  for (i=0; i<v->ndesc; i++)
    v->opt[i] = 0;

  for (f=v->famsout; f!=NULL; f=f->next) {
				/* are all vars in defined? */
    for (b=TRUE, i=0; i<f->n_in && b; i++)
      b = check_var_defined(f->in[i]);
    if (b) {			/* if yes, evaluate */
      eval_fam(f);
      for (i=0; i<v->ndesc; i++)
	(v->opt[i]) += f->degrees[i];
    }
    else {			/* if no, undefined and return */
      v->valdef = FALSE;
      for (i=0; i<v->ndesc; i++)
	v->opt[i] = 0;
      v->class = CUNDEF;
      if (debug_e) printf("leaving %s undef\n", v->name);
      return;
    }
  }

  for (i=0; i<v->ndesc; i++) w_sum += v->opt[i];
  if (w_sum != 0)
    for (i=0; i<v->ndesc; i++)
      v->opt[i] /= w_sum;
  if (eval_method==e_fuzzy) fuzzy_to_real(v);
  
  d = d_vector_max(v->opt, v->ndesc, &i);
  if (d>0) v->class = i;
  else v->class = v->c_default;
}

void evaluate(Str255 vname)
{
  var_type *v;

  if ((v=find_var(variables, vname))!=NULL)
    eval_var(v);
  else
    printf("error: variable %s not found\n", vname);
}


/* mark_var: marks all variables that v depends on with m, v included */

void mark_var(var_type *v, char m)
{
  int i;
  fams *f;

  v->mark = m;
  for (f=v->famsout; f!=NULL; f=f->next)
    for (i=0; i<f->n_in; i++)
      mark_var(f->in[i], m);
}

/* mark a list of vars and their dependents, BBB could be done better */

void mark_lvar(list_of_vars *lv, char m)
{
  for (; lv!= NULL; lv=lv->next)
    mark_var(lv->var, m);
}


void derive_var(var_type *v)
{
  int i;
  fams *f;

  if (!v->mark && v->famsout != NULL) {
    for (f=v->famsout; f!=NULL; f=f->next)
      for (i=0; i<f->n_in; i++)
	derive_var(f->in[i]);
    eval_var(v);
    v->mark = TRUE;
  }
}

/* get the mean interval of the variables that are continous-valued */

double get_k_interval(var_type *v)
{
  list_of_vars *lv, *l;
  int i, k, n;
  double d1, d2;
  char b;

  lv = find_leaves(v);
  d2 = 0.;
  for (i=0, n=0, l=lv; l!=NULL; i++, l=l->next) {
    v = l->var;

    for (k=0; k<v->ndesc; k++) v->opt[k]=0.;
    for (b=FALSE, k=0; k<v->ndesc && !b; k++)
      b = v->val >= v->desc[k].start &&
	v->val <= v->desc[k].start + v->desc[k].delta;
    k--;
    if (b) {
      v->opt[k]=1.;
      d1 = (v->val - v->desc[k].start)/v->desc[k].delta;
      d2 += d1;
    }
    else {
      printf("warning: %s out  of bounds %lf\n", v->name, v->val);
      v->opt[0]=1;
    }
    n++;
  }
  if (n>0)
    d2 = d2 / (double) n;
  else {
    printf("warning: using interval evaluation with no cont attribute\n");
    d2 = 0.5;
  }
  return d2;
}

void set_interval_val(var_type *v, double d)
{
  int i;
  char b = FALSE;

  for (i=0; i<v->ndesc && !b; i++)
    b = v->opt[i] >= 1.;
  if (b) {
    i--;
    v->val = v->desc[i].start + d * v->desc[i].delta;
    v->valdef = TRUE;
/*    printf("%d found %lf\n", i, v->val); */
  }
  else
    v->valdef = FALSE;
/*    printf("not ok\n"); */
}

void derive_lvar(list_of_vars *lv)
{
  int i;
  fams *f;
  list_of_vars *tlv;
  double k;

  for (tlv=lv; tlv!=NULL; tlv=tlv->next)
    mark_var(tlv->var, FALSE);
  for (tlv=lv; tlv!=NULL; tlv=tlv->next) {
    if (eval_method == e_interval) {
      k = get_k_interval(tlv->var);
      mark_var(tlv->var, FALSE);	/* this causes reevaluation */
    }
    derive_var(tlv->var);
    if (eval_method == e_interval) set_interval_val(tlv->var, k);
  }
}

void derive(Str255 vname)
{
  var_type *v;

  if ((v=find_var(variables, vname))!=NULL) {
    mark_var(v, FALSE);
    derive_var(v);
  }
  else
    printf("error: variable %s not found\n", vname);
}

/* get_q_val: takes a variable, does fuzzification, and returns a
   string which describes its qualitative value */

Str255 xqname;
char *get_q_val(var_type *v)
{
  int i;

  if (!v->valdef) {
    fprintf(stderr, 
	    "error: can't get qval from nondefined variable %s\n", v->name);
    return 0;
  }
  real_to_fuzzy(v);
  sprintf(xqname, "");
  for (i=0; i<v->ndesc; i++)
    if (v->opt[i] > 0.)
      sprintf(xqname, "%s %s/%5.3lf", xqname, v->desc[i].name, v->opt[i]);
  return xqname;
}

/****************************************************************************
EVALUATION ERRORS
these rutines are used to derive estimation errors for options
****************************************************************************/

/* stat_for_var: derives an error of the estimation err, where err is
   sqrt of sum of squares of differences btw expected and derived
   values, where they exist */

double get_opt_error(double ex, double der)
{
  double a, b, c, k, s=0.;
  double N;

  switch (ga_error_method) {
  case gae_sqr:
    return SQR(der - ex);
    break;
  case gae_abs:
    return ABS(der - ex);
    break;
  case gae_norm:
    /*	a = ABS(der - ex);
	if (ex != 0.0)
	b = 100 * 
	ABS((der - ex)/ex);
	else
	b = 100000;
	k = (250.0-ctx->ga.iter)/500.0+ 0.5;
	if (k>1) k=1; if (k<0) k=0;
	return k * b + (1-k) * a;
	return ABS(der - ex); */

    /*	a = (der/ex);
	a = SQR(a-1);
	return a; */

    if (ABS(der) > ABS(ex)) {a = der; b = ex;}
    else {a = ex; b = der;}
    if (a*b<0) a = ABS(a)+ABS(a+b);
    c = ABS(ABS(a/b)-1);

    /*	if (ex > -1.)
	printf("%lf %lf -> %lf\n", ex, der, a); */
    return c;
    break;
  case gae_perc:
    if (ex != 0.0)
      return 100 * 
	ABS((der - ex)/ex);
    else
      return 100000;		/* BBB this is stupid */
    break;
  case gae_sqrperc:
    if (ex != 0.0)
      return 100 * 
	SQR(ABS((der - ex)/ex));
    else
      return 100000;
    break;
  case gae_maxperc:
    if (ex != 0.0)
      s = 100*ABS((der - ex)/ex);
    else s += 100000;		/* BBB this is stupid */
    /* return MMAX(sum, s); */	/* error here */
    break;
  }
}

double stat_for_vnum(int vnum)
{
  list_of_opt *opt;
  int n_opt = 0;
  double sum = 0, s;

  for (opt=options; opt!=NULL; opt=opt->next, ++n_opt)
    if (opt->valdef[vnum] && opt->expect_def[vnum]) {
      sum += get_opt_error(opt->expect[vnum], opt->val[vnum]);
    }
    else {			/* what to do if undefined? */
      return 100e100;
    }
  switch (ga_error_method) {
  case gae_sqr:
    return sqrt(sum)/n_opt;
  case gae_perc:
  case gae_sqrperc:
  case gae_abs:
  case gae_norm:
    return sum/n_opt;
  case gae_maxperc:
    return sum;
  }
}

void stat_for_var(Str255 vname)
{
  int vnum;

  if (find_var_num(variables, vname, &vnum)!=NULL)
    printf("Error: %10.5lf\n", stat_for_vnum(vnum));
  else
    printf("error: variable %s not found\n", vname);
}

double sstat_for_itable(list_of_itables *it)
{
  double sum=0;
  int i;

  if (it->out->ctype == ct_nominal) {
    int n=0;
    for (i=0; i<it->n_inst; i++)
      if (it->estq[i] != it->qval[i][it->n_in]) n++;
    return n / (double) it->n_inst;
  }

  for (i=0; i<it->n_inst; i++)
    if (it->estdef[i])
      sum += get_opt_error(it->val[i][it->n_in], it->est[i]);
  switch (ga_error_method) {
  case gae_sqr:
    return sum = sqrt(sum)/it->n_inst;
  case gae_perc:
  case gae_sqrperc:
  case gae_abs:
  case gae_norm:
    return sum/it->n_inst;
  case gae_maxperc:
    return sum;
  }
}

double stat_for_itable(Str255 iname)
{
  list_of_itables *it;

  if ((it = find_itable(iname)) == NULL) {
    printf("error: instance table %s not defined\n", iname);
    return;
  }
  eval_instance_table(it);
  return sstat_for_itable(it);
}


/****************************************************************************
eval.c

This is new evaluation which should be much faster than the old one
****************************************************************************/

void neval_var(var_type *v)
{
  fams *f = v->famsout;
  int i, nud, p, c, nn;
  char *att, cd = FALSE;
  char use_hash = FALSE;
  

  if (f->n_in>100) {use_hash=TRUE; printf("use hash\n"); build_hash_table(f);}

  att = c_vector(f->n_in);

  for (i=0; i<f->n_in; i++) att[i] = f->in[i]->class;

  for (nud=0, i=0; i<f->n_in; i++) 
    if (f->in[i]->class == CUNDEF) {
      nud++; att[i] = 0; p = i;
      if (debug_e) printf("undef %s\n", f->in[i]->name);
    }

  if (nud==0) {
    v->class = get_rule(f, att);
    if (v->class == CUNDEF) {
      default_used = cd = TRUE;
      v->class = v->c_default;
      /* printf("def %s\n", v->name); */
    }
    if (debug_e) printf("%s <-- %s %s\n", v->name, v->desc[v->class].name, 
			cd?"(def)":"");
  }
  else {			/* some undefined children */
    for (i=0; i<v->ndesc; i++) v->opt[i] = 0;
    nn = 0;
    while (p>=0) {
      nn++;

/*      for (i=0; i<f->n_in; i++) printf("%d ", att[i]); printf("\n"); */

      c = get_rule(f, att);
      if (c==CUNDEF) {c = v->c_default;}
      v->opt[c]++;
      
/*      printf("   %d: %d\n", nn, c); */

      while (p>=0) {
	if (att[p] < f->in[p]->ndesc-1) {
	  att[p]++;
	  for (i=p+1; i<f->n_in; i++)
	    if (f->in[i]->class==CUNDEF) { p = i; att[i]=0; }
	  break;
	}
	for (p--; p>=0; p--)
	  if (f->in[p]->class==CUNDEF)
	    break;
      }
    }
    for (i=0; i<v->ndesc; i++) v->opt[i] = v->opt[i]/(double)nn;
    v->class = CUNDEF;
    for (i=0; i<v->ndesc; i++) {
      if (v->opt[i]>=1.) {v->class=i; break;}
      if (v->opt[i]>0) break;
    }
    if (debug_e) 
      if (v->class!=CUNDEF) 
	printf("%s <--- %s\n", v->name, v->desc[v->class].name);
      else {
	printf("%s <--- ", v->name);
	for (i=0; i<v->ndesc; i++) 
	  printf("%s/%lf ", f->in[i]->name, v->opt[i]);
	printf("\n");
      }
  }
  free(att);
}

/* nrderive_var: derives the value of var and for that uses a
   recursive definition of how to derive it */

void nrderive_var(var_type *v)
{
  int i;
  fams *f = v->famsout;

  if (f != NULL) {
    for (i=0; i<f->n_in; i++)
      nrderive_var(f->in[i]);
    neval_var(v);
  }
  else {
    if (v->class == CUNDEF) {	/* in case leaf is undefined */
      for (i=0; i<v->ndesc; i++)
	v->opt[i] = 1./(double)v->ndesc;
    }
  }
}

/* nderive_var: uses nrderive_var to derive the value of v. In case
   this is fuzzy, it uses probabilities and apriory probabilities to
   still derive its unique class. The tie is resolved arbitrarily only
   if two or more classes have the same derived and apriory
   probabilities. */

void nderive_var(var_type *v)
{
  int i, c;
  double maxp=-1, maxa=-1;
  fams *f = v->famsout;

  nrderive_var(v);
  if (v->class == CUNDEF) {
    /* get the class that has the highest probability, if tie, than
       select the one which has the highest apriory */
    for (i=0; i<v->ndesc; i++) 
      if (v->opt[i]>maxp || (v->opt[i]==maxp && v->c_apriory[i] > maxa)) {
	c = i;
	maxp = v->opt[i];
	maxa = v->c_apriory[i];
      }
    v->class = c;
  }
}

