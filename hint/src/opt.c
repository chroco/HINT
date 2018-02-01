/****************************************************************************
opt.c

Rutines within this file manage the options, i.e., values of all
variables. Variable value is fuzzy, i.e., a value belongs to a fuzzy
sets (descriptors) to a certain degree.

Unlike in DEX, options store the values of all variables, not just
those in leafs of DAGs. This is useful when experimenting with a
subset of DAG. We can also set an internal node of a DAG, usually
introducing the incosistency of variable's values and rule base.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define GL extern
#include "sds.h"

/* undef_opt: makes the current option undefined */

void undefine_opt()
{
  list_of_vars *lv;
  int i;

  for (lv=variables; lv!=NULL; lv=lv->next) {
    lv->var->valdef = FALSE;
    lv->var->expect_def = FALSE;
    for (i=0; i<lv->var->ndesc; i++)
      lv->var->opt[i] = 0;
  }
}


/* set_var_val_num: set a value of variable as a number. This number
   is then immediately transfered to qualitative representation. */

void set_var_val_num(Str255 vname, double val)
{
  var_type *v;

  if ((v=find_var(variables, vname))!=NULL) {
    v->val = val;
    v->valdef = TRUE;
    if (eval_method==e_fuzzy) real_to_fuzzy(v);
    else if (eval_method==e_interval) real_to_interval(v);
				/* BBB here crisp to fuzzy */
  }
  else printf("error: variable %s not found\n", vname);
}

/* set_var_val_q: set a value of a variable to a single qualitative
   value, i.e., all other memberships are 0, and this one is 1 */

void set_var_val_q(Str255 vname, Str255 qname)
{
  int i;
  var_type *v;

  v = find_var(variables, vname);
  if (v==NULL) {
    printf("error: variable %s not found\n", vname);
    return;
  }
  
  for (i=0; i<v->ndesc && strcmp(v->desc[i].name, qname); i++);
  if (i<v->ndesc) {
    for (i=0; i<v->ndesc; i++) {
      if (strcmp(v->desc[i].name, qname))
	v->opt[i] = 0.;
      else {
	v->opt[i] = 1.;
	if (eval_method == e_fuzzy) {
	  v->val = v->desc[i].cent;
	  v->valdef = TRUE;
	}
      }
    }
  }
  else
    printf("error: %s is not legal qval for %s\n", qname, vname);
}

void set_var_val_qlist(Str255 vname, list_of_q *qroot)
{
  int i;
  var_type *v;
  list_of_q *q;
  double tot;

  v = find_var(variables, vname);
  if (v==NULL) {
    printf("error: variable %s not found\n", vname);
    return;
  }
  
  for (tot=0, q=qroot; q!=NULL; q=q->next)
    tot += q->degree;
  for (q=qroot; q!=NULL; q=q->next) {
    for (i=0; i<v->ndesc && strcmp(v->desc[i].name, q->desc); i++);
    if (i<v->ndesc)
      v->opt[i] = q->degree / tot;
    else
      printf("error: %s is not legal qval for %s\n", q->desc, vname);
  }
  fuzzy_to_real(v);
}

/* set_var_expect: set expected field for a variable vname */

void set_var_expect(Str255 vname, double expect)
{
  var_type *v;

  v = find_var(variables, vname);
  if (v==NULL)
    printf("error: variable %s not found\n", vname);
  else {
    v->expect_def = TRUE;
    v->expect = expect;
  }
}

/****************************************************************************
setting, saving, deleting, and listing the options
****************************************************************************/

int count_var(list_of_vars *v)
{
  int n;

  for (n=0; v!=NULL; n++, v=v->next);
  return n;
}

int find_var_pos(list_of_vars *v, Str255 vname)
{
  int n;

  for (n=0; v!=NULL && strcmp(v->var->name, vname); n++, v=v->next);
  return n;
}

void list_options()
{
  list_of_opt *opt;
  
  for (opt=options; opt!=NULL; opt=opt->next)
    printf("%s ", opt->name);
  printf("\n");
}

list_of_opt *find_opt(Str255 oname, char d)
{
  list_of_opt *opt;
  for (opt=options; opt!=NULL && strcmp(opt->name, oname); opt=opt->next);
  if (d && opt==NULL) 
    printf("error: option %s does not exists\n", oname);  
  return opt;
}

void del_opt(Str255 oname)
{
  list_of_opt *opt;
  list_of_vars *v;
  int i, n;
 
  if ((opt=find_opt(oname,TRUE))!=NULL) {
    if (opt == options)
      options = options->next;
    else {
      opt->prev->next = opt->next;
      if (opt->next != NULL) opt->next->prev = opt->prev;
      free(opt->val);
      n = count_var(variables);
      for (i=0; i<n; i++)
	free(opt->opt[i]);
      free(opt->opt);		/* BBB some more freeing necessary */
      free(opt);
    }
  }
}

void set_opt(Str255 oname)
{
  list_of_opt *opt;
  list_of_vars *v;
  int i, j, n;
  
  if ((opt=find_opt(oname,FALSE))==NULL) { /* have to create it */
    opt = (list_of_opt *) malloc(sizeof(*opt));
    strcpy(opt->name, oname);
    opt->next = options;
    if (options!=NULL) options->prev = opt;
    options = opt;

    n = count_var(variables);
    opt->opt = (double **) malloc(sizeof(*opt->opt)*n);
    opt->val = (double *) malloc(sizeof(*opt->val)*n);
    opt->valdef = (char *) malloc(sizeof(*opt->valdef)*n);
    opt->expect = (double *) malloc(sizeof(*opt->expect)*n);
    opt->expect_def = (char *) malloc(sizeof(*opt->expect_def)*n);

    for (i=0, v=variables; v!=NULL; i++, v=v->next)
      opt->opt[i] = (double *) malloc(sizeof(double)*v->var->ndesc);
  }
/*  else
    printf("found it\n");  */

  for (i=0, v=variables; v!=NULL; i++, v=v->next) {
    for (j=0; j<v->var->ndesc; j++)
      opt->opt[i][j] = v->var->opt[j];
    opt->val[i] = v->var->val;
    opt->valdef[i] = v->var->valdef;
    opt->expect[i] = v->var->expect;
    opt->expect_def[i] = v->var->expect_def;
/*    printf("yy %s %5.2lf\n", opt->name, opt->val[i]); */
  }
}

void select_opt(Str255 oname)
{
  list_of_opt *opt;
  list_of_vars *v;
  int i, j, n;
  
  if ((opt=find_opt(oname,TRUE))!=NULL) {
    n = count_var(variables);
    for (i=0, v=variables; v!=NULL; i++, v=v->next) {
      for (j=0; j<v->var->ndesc; j++)
	v->var->opt[j] = opt->opt[i][j];
      v->var->val = opt->val[i];
      v->var->valdef = opt->valdef[i];
      v->var->expect = opt->expect[i];
      v->var->expect_def = opt->expect_def[i];
/*      printf("xx %s %5.2lf\n", opt->name, opt->val[i]); */
    }
  }
}

void list_var_opt_all(Str255 vname)
{
  list_of_vars *lv;
  list_of_opt *opt;
  int i, nv;
  char b;

  for (nv=0, lv=variables; lv!=NULL && strcmp(lv->var->name, vname);
       nv++, lv=lv->next);
  if (lv!=NULL) {
    for (opt=options; opt!=NULL; opt=opt->next) {
      printf("%s\t", opt->name);
      if (opt->valdef[nv])
	printf("%15.5lf\t", opt->val[nv]); else printf("      \t\t\t");
      if (opt->expect_def[nv])
	printf("%15.5lf\t",opt->expect[nv]); else printf("      \t\t\t");

      if (opt->valdef[nv] && opt->expect_def[nv]) {
	printf("%8.4lf\t", get_opt_error(opt->expect[nv], opt->val[nv]));
/*	printf("%8.4lf%%\t", 100.0 * ABS((opt->val[nv]-opt->expect[nv])/
			       opt->expect[nv])); */
      }
      else printf("\t\t");

      b = TRUE;
      for (i=0; i<lv->var->ndesc; i++)
	if (opt->opt[nv][i] > 0) {
	  printf("%s/%4.2lf\t", lv->var->desc[i].name, opt->opt[nv][i]);
	  b = FALSE;
	}
      if (b) printf("?");
      printf("\n");
    }
  }
  else
    printf("error: variable %s not known\n", vname);
}

/* plot_var_opt_all: takes a variable defined with its name vname, and
   for all known options plots a graph which shows the approximation
   error in percents with the respect of the expected value */

void plot_var_opt_all(Str255 vname)
{
  list_of_vars *lv;
  list_of_opt *opt;
  int i, nv, np=0;
  char b;
  double err;
  char fname[] = "000.plot";
  Str255 s;
  FILE *f;

  for (nv=0, lv=variables; lv!=NULL && strcmp(lv->var->name, vname);
       nv++, lv=lv->next);
  if (lv==NULL) {
    printf("error: variable %s not known\n", vname);
    return;
  }

  f = fopen(fname, "w");
  if (f==NULL) {
    printf("error: can't open temporary %s file\n", fname);
    return;
  }

  for (opt=options; opt!=NULL; opt=opt->next) {
    if ((opt->valdef[nv]) && (opt->expect_def[nv])) {
      err = 100.0 * ABS((opt->val[nv]-opt->expect[nv])/opt->expect[nv]);
      if (err < ga_max_error) {
	np++;
	fprintf(f, "%lf %lf\n", opt->expect[nv], err);
      }
    }
  }
  fclose(f);
  if (np>0) {
    sprintf(s, "gnugraph -L \"Error (%s)\" -S -m \"-1\" < %s | xplot -f 7x13",
	    vname, fname);
    system(s);
  }
  else
    printf("error: no options plotted; none are defined or error too big\n");
}

void list_var_opt_str(Str255 vname, list_of_str *sroot)
{
  list_of_vars *lv;
  list_of_opt *opt;
  list_of_str *s;
  int i, nv;
  char b;

  for (nv=0, lv=variables; lv!=NULL && strcmp(lv->var->name, vname);
       nv++, lv=lv->next);
  if (lv!=NULL) {
    for (s=sroot; s!=NULL; s=s->next)
      if ((opt=find_opt(s->str, TRUE)) != NULL) {
	printf("%s\t", opt->name);
	if (opt->valdef[nv])
	  printf("%6.3lf\t", opt->val[nv]); else printf("      \t");
	if (opt->expect_def[nv])
	  printf("%6.3lf\t",opt->expect[nv]); else printf("      \t");
	for (i=0; i<lv->var->ndesc; i++)
	  if (opt->opt[nv][i] > 0) {
	    printf("%s/%4.2lf\t", lv->var->desc[i].name, opt->opt[nv][i]);
	    b = FALSE;
	  }
	if (b) printf("?");
	printf("\n");
      }
      else
	printf("error: option %s not known\n", s->str);
  }
  else
    printf("error: variable %s not known\n", vname);
}

/****************************************************************************
INSTANCE TABLES
provide more compact representation of options
****************************************************************************/

extern list_of_vars *lv_desc;	/* leafs of variable, see learn.c */

/* find_itable: returns a pointer to an instance table named name. */

list_of_itables *find_itable(Str255 name)
{
  list_of_itables *v;

  for (v=itables; v!=NULL; v=v->next)
    if (!strcmp(v->name, name))
      return v;
  return NULL;
}

/* find_qval: finds if v has a qualitative value set to s. If so, c is
   the position of this qval. */

char find_qval(var_type *v, Str255 s, char *c)
{
  int j;

  for (j=0; j<v->ndesc && strcmp(v->desc[j].name, s); j++);
    if (j<v->ndesc) {
      *c = j;
      return TRUE;
    }
    else return FALSE;
}

/* add_itable: build an instance table iname, using the leaf
   attributes of vname and list of instances ln. Instance table is
   then added to the list of instance tables. */

void add_itable(Str255 iname, Str255 vname, list_of_str *ln)
{
  list_of_itables *it;
  var_type *v;
  list_of_vars *lv;
  list_of_str *l;
  int n_in, n, i, j;
  char c;

  if ((it = find_itable(iname)) != NULL) {
    printf("error: %s defined, delete it first to redifine\n", iname);
    return;
  }
  if ((v=find_var(variables, vname))==NULL) {
    printf("error: variable %s not found\n", vname);
    return;
  }

				/* create an instance table */
  lv_desc = find_leaves(v);
  for (lv=lv_desc, n_in=0; lv!=NULL; lv=lv->next, n_in++);
  for (l=ln, n=0; l!=NULL; l=l->next, n++);
  if (n % (n_in+1) != 0) {
    printf("error: illegal specification of instances\n");
    return;
  }
  it = (list_of_itables *) malloc(sizeof(*it));
  strcpy(it->name, iname);
  it->n_in = n_in;
  it->n_inst = n / (n_in+1);
  it->out = v;
  it->mark = c_vector(it->n_inst);
  it->in = (var_type **) malloc(sizeof(*it->in)*n_in);
  for (lv=lv_desc, i=0; lv!=NULL; lv=lv->next, i++)
    it->in[i] = lv->var;  
  it->val = d_matrix(it->n_inst, it->n_in+1);
  it->qval = c_matrix(it->n_inst, it->n_in+1);
  it->tval = c_matrix(it->n_inst, it->n_in+1);
  it->estdef = c_vector_ini(it->n_inst, FALSE);
  it->estq = c_vector(it->n_inst);
  it->est = d_vector(it->n_inst);
  it->next = itables;
  itables = it;

				/* set the values of instance table */
  for (l = ln, i=0; i<it->n_inst; i++)
    for (j=0; j<it->n_in+1; j++) {
      if (j<it->n_in) v=it->in[j]; else v=it->out;
      if (find_qval(v, l->str, &c)) {
	it->qval[i][j] = c;
      }
      else if (!strcmp(l->str,"?")) {
	it->qval[i][j] = CUNDEF;
	it->val[i][j] = -1e300;
      }
      else {
	it->qval[i][j] = CUNDEF;
	it->val[i][j] = atof(l->str);
	it->in[j]->ctype = ct_contin;
      }
      l=l->next;
    }
}

void list_instance_tables()
{
  list_of_itables *it;
  for (it=itables; it!=NULL; it=it->next)
    printf("%s ", it->name);
  printf("\n");
}

void eval_instance_table(list_of_itables *it)
{
  int i, j;
  list_of_vars lv;
  void sload_instance_from_table(int inst, list_of_itables *it);

  for (i=0; i<it->n_inst; i++) {
    sload_instance_from_table(i, it);
    lv.next = NULL; lv.var = it->out;
    derive_lvar(&lv);
    if (lv.var->valdef) {
      it->estdef[i] = TRUE;
      it->est[i] = lv.var->val;
    }
    d_vector_max(lv.var->opt, lv.var->ndesc, &j);
    it->estq[i]=j;
/*    printf("-> %d\n", j); */
  }
}

print_ins_entry(var_type *v, list_of_itables *it, int i, int j)
{
  if (it->qval[i][j]==CUNDEF) printf("%8.4lf ", it->val[i][j]);
  else printf("%s ", strn(v->desc[it->qval[i][j]].name,1));
}

void list_instance_table(Str255 iname, char do_derive)
{
  list_of_itables *it;
  int i, j;
  void load_instance_from_table(int inst, Str255 iname);

  if ((it = find_itable(iname)) == NULL) {
    printf("error: instance table %s does not exist\n", iname);
    return;
  }
  printf("     ");
  for (i=0; i<it->n_in; i++)
    printf("%s ", strnc(it->in[i]->name, 8));
  printf(" -> %s  est      err\n", strnc(it->out->name, 8));

  if (do_derive) eval_instance_table(it);
  for (i=0; i<it->n_inst; i++) {
    printf("%-3d: ", i);
    for (j=0; j<it->n_in; j++)
      print_ins_entry(it->in[j], it, i, j);
    printf("    ");
    print_ins_entry(it->out, it, i, j);

    if (it->out->ctype == ct_nominal) {
      printf(" %s", it->out->desc[it->estq[i]]);
    }
    else {
      if (it->out->valdef) {
	printf(" %8.5lf %8.5lf", it->est[i], 
	       get_opt_error(it->val[i][it->n_in], it->est[i]));
      }
    }
    printf("\n");
  }
}

/* get_qval: gets the qualitative value of variable v, which has been
   stored in the instance table row i column k. If this has been
   stored as a real value, it is converted to a qualitative value. */

char get_qval(var_type *v, list_of_itables *it, int i, int k)
{
  int l;

  if (it->qval[i][k]!=CUNDEF)
    return it->qval[i][k];
  else if (it->val[i][k] == -1e300) {
    return CUNDEF;
    printf("xxxxxxx\n");
  }
  else {
    for (l=0; l<v->ndesc-1; l++)
      if (it->val[i][k] >= v->desc[l].start &&
	  it->val[i][k] < v->desc[l].start + v->desc[l].delta)
		break;
    return l;
  }
}

void sderive_rtable_from_itable(fams *rt, list_of_itables *it)
{
  char *att=NULL, *class=NULL;
  int max, imax, i, j, k, l;
  char b, c;

  /* here one should extensively check compatibility of rt and it */
  /* instead, we do this */
  if (rt->n_in == it->n_in) {
    b = TRUE;
    for (i=0; i<rt->n_in && b; i++)
      b = rt->in[i] == it->in[i];
  }
  else b = FALSE;
  if (!b) {
    printf("error: instance table and target rule table not compatible\n");
    return;
  }

  del_rules(rt);
  for (i=0; i<it->n_inst; i++) it->mark[i] = TRUE;

  att = c_vector(it->n_in);
  class = c_vector(it->out->ndesc);

  for (i=0; i<it->n_inst; i++)
    if (it->mark[i]) {
      for (j=0; j<it->out->ndesc; j++) class[j] = 0;
      it->mark[i] = FALSE;

      for (k=0; k<it->n_in; k++)
	att[k] = get_qval(it->in[k], it, i, k);
      c = get_qval(it->out, it, i, k);
      class[c]++;

/*      printf("yyy ");
      for (j=0; j<it->n_in; j++) 
	printf("%d ", att[j]); 
      printf(" -> %d\n", c); */

      for (j=i+1; j<it->n_inst; j++)
	if (it->mark[j]) {
	  for (k=0, b=TRUE; k<rt->n_in && b; k++)
	    b = att[k]==get_qval(it->in[k], it, j, k);
	  if (b) {
	    it->mark[j] = FALSE;
	    class[get_qval(it->out, it, j, k)]++;
	  }
	}
      max = class[0]; imax = 0;
      for (j=0; j<it->out->ndesc; j++)
	if (class[j] > max) {
	  max = class[j];
	  imax = j;
	}
      set_rule(rt, att, imax, autom);
    }
  FREE(att); FREE(class);
}

/* set_qval_ins: using defined intervals, sets the nominal values for
   each instance (finds the position within an interval) */

void set_tval_ins(list_of_itables *it)
{
  int i, k;

  for (i=0; i<it->n_inst; i++) {
    
    for (k=0; k<it->n_in; k++) {
      if (it->in[k]->ctype == ct_contin)
	it->tval[i][k] = get_qval(it->in[k], it, i, k);
      else
	it->tval[i][k] = it->qval[i][k];
    }
  }
}


void derive_rtable_from_itable(Str255 rname, Str255 iname, char is_table)
{
  list_of_itables *it;
  var_type *v;
  fams *rt;

  if ((it = find_itable(iname)) == NULL) {
    printf("error: instance table %s does not exist\n", iname);
    return;
  }

  if (is_table) {
    if ((rt=find_table(rname)) == NULL) {
      printf("error: rule table %s does not exist\n", rname);
      return;
    }
  }
  else {
    if ((v = find_var(variables, rname)) == NULL || v->famsout == NULL) {
      printf("error: variable %s not found or inappropriate\n", rname);
      return;
    }
    rt = v->famsout;
  }
  sderive_rtable_from_itable(rt, it);
}

void sload_instance_from_table(int inst, list_of_itables *it)
{
  int i,j;

  for (i=0; i<it->n_in; i++) {
    for (j=0; j<it->in[i]->ndesc; j++) it->in[i]->opt[j]=0;
    if (it->in[i]->ctype == ct_contin) {
      it->in[i]->val = it->val[inst][i];
      real_to_interval(it->in[i]);
    } 
    else
      it->in[i]->opt[it->qval[inst][i]]=1;
  }
/*  it->out->expect = it->qval[inst][it->n_in];
  it->out->expect_def = TRUE; */

/*  for (i=0; i<it->n_in; i++) {
    d_vector_max(it->in[i]->opt, it->in[i]->ndesc, &j);
    printf("%s ", it->in[i]->desc[j].name);
  }
  printf("\n"); */
}

void load_instance_from_table(int inst, Str255 iname)
{
  list_of_itables *it;

  if ((it = find_itable(iname)) == NULL) {
    printf("error: instance table %s does not exist\n", iname);
    return;
  }
  if (inst < 0 || inst > it->n_inst-1) {
    printf("error: no such instance in table %s\n", iname);
    return;
  }
  sload_instance_from_table(inst, it);
}

void load_instances_from_table(Str255 iname)
{
  list_of_itables *it;
  int i;
  char *name;

  if ((it = find_itable(iname)) == NULL) {
    printf("error: instance table %s does not exist\n", iname);
    return;
  }
  for (i=0; i<it->n_inst; i++) {
    name = gen_opt_name();
    printf("%s ", name); if (i%10==0 && i!=0) printf("\n");
    load_instance_from_table(i, iname);
    set_opt(name);
  }
  printf("\n");
}

void derive_int_using_itable(Str255 iname, list_of_vars *lv)
{
  list_of_vars *l;
  list_of_itables *it;
  int i, j, k;
  double min, max;
  double delta, start;
  char b;

  if ((it = find_itable(iname)) == NULL) {
    printf("error: instance table %s does not exist\n", iname);
    return;
  }
  for (l=lv; l!=NULL; l=l->next) {
    if (l->var==it->out) {
      b = TRUE; k = it->n_in;
    }
    else {
      for (b=FALSE, i=0; i<it->n_in && !b; i++) 
	b = (it->in[i]==l->var);
      k = i-1;
    }
    if (b) {
      min = max = it->val[0][k];
      for (i=1; i<it->n_inst; i++) {
	if (min > it->val[i][k]) min = it->val[i][k];
	if (max < it->val[i][k]) max = it->val[i][k];
      }
      delta = (max-min)/l->var->ndesc;
      start = min;
      for (i=0; i<l->var->ndesc; i++) {
	l->var->desc[i].start = start;
	l->var->desc[i].delta = delta;
	start += delta;
      }
/*      printf("%s %d %lf %lf\n", l->var->name, k, min, max); */
/* this is tricky, but prevents out of bounds due to numerical error */
      i--;
/*      it->in[k]->desc[i].start = max - it->in[k]->desc[i].delta; */
      l->var->desc[i].start += l->var->desc[i].delta / 100000.0;
    }
    else
      printf("error: %s not in instance table %s\n", l->var->name, iname);
  }
}

/* split_instance_table: split the instance table to two disjoint tables 
   it1 = it2 + it3
   The % of instances in it2 is given */

void sel_copy_it(list_of_itables *it, list_of_itables *itf, char *sel, char m)
{
  int i, j, k, count;

  it->n_in = itf->n_in;

  for (i=0, count=0; i<itf->n_inst; i++) if (sel[i]==m) count++;
  it->n_inst = count;

  it->out = itf->out;
  it->mark = c_vector(it->n_inst);
  it->in = (var_type **) malloc(sizeof(*it->in)*itf->n_in);

  it->val = d_matrix(it->n_inst, it->n_in+1);
  it->qval = c_matrix(it->n_inst, it->n_in+1);
  it->tval = c_matrix(it->n_inst, it->n_in+1);

  it->estdef = c_vector_ini(it->n_inst, FALSE);
  it->estq = c_vector(it->n_inst);
  it->est = d_vector(it->n_inst);

  for (i=0; i<itf->n_in; i++) it->in[i] = itf->in[i];
  for (i=0, j=0, count=0; i<itf->n_inst; i++)
    if (sel[i]==m) {
      it->estdef[j] = itf->estdef[i];
      it->estq[j] = itf->estq[i];
      it->est[j] = itf->est[i];
      for (k=0; k<it->n_in+1; k++) {
	it->val[j][k] = itf->val[i][k];
	it->qval[j][k] = itf->qval[i][k];
	it->tval[j][k] = itf->tval[i][k];
      }
      j++;
    }
}

void split_instance_table(char *id1, char *id2, char *id3, double perc)
{
  list_of_itables *it1, *it2, *it3;
  char *sel;
  int i, j, n;

printf("%s %s %s\n", id1, id2, id3);
  if ((it1 = find_itable(id1)) == NULL) {
    printf("error: instance table %s not defined\n", id1);
    return;
  }
  
  sel = c_vector_ini(it1->n_inst, FALSE);
  n = it1->n_inst * perc / 100.;
    for (i=0; i<n;) {
      j = it1->n_inst * rnd1e();
      if (!sel[j]) {
	sel[j] = TRUE;
	i++;
      }
    }

  it2 = (list_of_itables *) malloc(sizeof(*it2));
  strcpy(it2->name, id2);
  it3 = (list_of_itables *) malloc(sizeof(*it3));
  strcpy(it3->name, id3);
  
/*  memcpy(it2, it1, sizeof(*it1));
  memcpy(it3, it1, sizeof(*it1)); */

  sel_copy_it(it2, it1, sel, TRUE);
  sel_copy_it(it3, it1, sel, FALSE);


  it2->next = itables;
  itables = it2;
  it3->next = itables;
  itables = it3;  
}


/****************************************************************************
Discretization
****************************************************************************/

void discretize(Str255 id, Str255 iname)
{
  list_of_itables *it;
  var_type *v;
  fams *rt;
  char *att=NULL, b;
  int i, j, class;
  fams *f;

  if ((it = find_itable(iname)) == NULL) {
    printf("error: instance table %s does not exist\n", iname);
    return;
  }

  if ((f = find_fam_vartable(id)) == NULL) return;

  /* here one should extensively check compatibility of f and it */
  /* instead, we do this */
  if (f->n_in == it->n_in) {
    b = TRUE;
    for (i=0; i<f->n_in && b; i++)
      b = f->in[i] == it->in[i];
  }
  else b = FALSE;
  if (!b) {
    printf("error: instance table and target rule table not compatible\n");
    return;
  }

  del_rules(f);
  att = c_vector(it->n_in);

  for (i=0; i<it->n_inst; i++) {
    for (j=0; j<f->n_in; j++)
      att[j]= get_qval(it->in[j], it, i, j);
    class = get_qval(it->out, it, i, j);
    set_rule(f, att, class, autom);
  }
  FREE(att);
}
