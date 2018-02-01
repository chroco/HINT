/****************************************************************************
id3.c

Builds an id3-like tree out of a rule table.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#define GL extern
#include "sds.h"

int *maj_class;			/* vector that serves 
				   determination of majority class */

/* copy_sel_rules: from existing fam creates a fam that has all input
attributes but sindx, and all entries with attribute sindx being
sval. Can be used to copy fams completely. */

fams *copy_fam_sel(fams *f, int sindx, char sval)
{
  int i, j;
  char class, *att, *natt;
  rule_type rtp;
  fams *t;

  t = (fams *) malloc(sizeof(*f));
  t->hash_table = NULL;
  t->er = encode_rules;
  t->n_in = (sindx>=0 && sindx<f->n_in) ? f->n_in-1 : f->n_in;
  t->out = f->out;
  t->in = (var_type **) malloc(sizeof(*t->in) * t->n_in);
  for (j=0, i=0; i<f->n_in; i++) if (i!=sindx) t->in[j++] = f->in[i];
  for (j=1, i=0; i<t->n_in; i++) if (i!=sindx) j *= t->in[i]->ndesc;
  allocate_rules(t, j);

  if (sindx<0 || sindx>=f->n_in) copy_rules(f, t);
  else {
    natt = c_vector(f->n_in);
    reset_next_rule(f);
    while (get_next_rulea(f, &class, &rtp, &att)) {
/*      list_rule(f, class, att); */
      if (att[sindx] == sval) {
	for (i=0; i<sindx; i++) natt[i] = att[i];
	for (; i<t->n_in; i++) natt[i] = att[i+1];
	set_rule(t, natt, class, rtp);
      }
    }
  }
  return t;
}

extern double *im_one;		/* information measure for attribute */
extern double *w_one;		/* relieff measure for attributes */
extern double *gr_one;		/* gain-ratio measure */
extern double *gini_one;	/* gini measure */

dec_tree *tdidt(fams *f)
{
  char class, c, b;
  rule_type rtp;
  dec_tree *dt = (dec_tree *) malloc(sizeof(*dt));
  int nrules = 0, i, best;
  double *crit, bestc;			/* criteria */
  fams *nf;

  deb_dec=5;

  for (i=0; i<f->out->ndesc; i++) maj_class[i]=0;
  reset_next_rule(f);
  if (get_next_rule(f, &class, &rtp)) {
    c = class;
    maj_class[class]++;
    nrules++; b=TRUE;
    while (get_next_rule(f, &class, &rtp)) {
      maj_class[class]++;
      b = b && (c==class);
      nrules++;
    }
  }

  dt->nrules = nrules;
  dt->fam = f;

  if (nrules==0) {
    dt->is_leaf = TRUE;
    dt->class = CUNDEF;
  }
  else if (b) {
    dt->is_leaf = TRUE;
    dt->class = c;
  }
  else if (f->n_in == 0) {
    printf("Here\n");
    fprintf(lfile, "TDIN conflict, %s: ", f->out->name);
    for (i=0; i<f->out->ndesc; i++) 
      fprintf(lfile, "%s/%d ", f->out->desc[i].name, maj_class[i]);
    fprintf(lfile, "\n");
    dt->class = CUNDEF;
  }
  else {
    dt->is_leaf = FALSE;

    informativity_one(f);
    crit = log_inform ? im_one : (log_gr ? gr_one : gini_one);

    if (deb_dec>2) {
      printf("Crit: ");
      for (i=0; i<f->n_in; i++)
	printf("(%s,%8.3lf) ", f->in[i]->name, crit[i]);
      printf("\n");
    }

    for (best=0, bestc=crit[0], i=1; i<f->n_in; i++)
      if (crit[i] >= bestc) {best=i; bestc = crit[i];}
    if (deb_dec>2) printf("Best %s\n", f->in[best]->name);
    
    dt->v = f->in[best];
    dt->next = (dec_tree **) malloc(sizeof(*dt->next) * f->in[best]->ndesc);
    for (i=0; i<f->in[best]->ndesc; i++) {
      nf = copy_fam_sel(f, best, (char)i);
/*      printf("FOR %s = %s THE TABLE IS\n", f->in[best]->name, f->in[best]->desc[i].name);
      list_rule_table(f, 2); */
      dt->next[i] = tdidt(nf);
    }
  }
  return dt;
}

list_dt_rec(dec_tree *dt, int l)
{
  int i;

  if (dt->is_leaf) {
    if (dt->class!=CUNDEF) printf("%s\n", dt->fam->out->desc[dt->class]);
    else printf("?\n");
  }
  else {
    printf("%s\n", dt->v->name);
    for (i=0; i<dt->v->ndesc; i++) {
      pspace((l+1)*3); printf("%s: ", dt->v->desc[i].name);
      list_dt_rec(dt->next[i], l+1);
    }
  }
}

void list_dt(dec_tree *dt)
{
  printf("%s=?\n", dt->fam->out->name);
  list_dt_rec(dt, 0);
}

void list_dec_tree(Str255 vname)
{
  var_type *v;

  FIND_VAR(v, vname);
  list_dt(v->famsout->dt);
}

void derive_dec_tree(Str255 vname, int crit)
{
  var_type *v;
  fams *f, *nf;
  char li=log_inform, lg=log_gini, lgr=log_gr, lr=log_relieff;

  if (crit==0) log_inform=TRUE;
  else if (crit==1) log_gini=TRUE;
  else if (crit==2) log_gr=TRUE;
  else if (crit==3) log_relieff=TRUE;
    
  FIND_VAR(v, vname);
  maj_class = i_vector(v->ndesc);
  f = v->famsout;
  nf = copy_fam_sel(f, -1, 0);
  v->famsout->dt = tdidt(f);
  list_dt(v->famsout->dt); 

  log_inform=li, log_gini=lg, log_gr=lgr, log_relieff=lr;
  FREE(maj_class);
}
