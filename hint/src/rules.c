/****************************************************************************
fams.c

Rutines for assertion, modification, and deletion of FAMs. 
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#define GL extern
#include "sds.h"

void del_rules(fams *f);
void reset_next_rule(fams *f);

FILE *log_tests = NULL;		/* file to log cn2.5 data */
char do_save_tests = FALSE;	/* really do write down the rules? */
var_type **in_order = NULL;		/* order of variables as
				   saved before the test */

/****************************************************************************
set and get rules
****************************************************************************/

void set_rule(fams *f, char *att, char class, rule_type rtp)
{
  int i;
  rule_list *rl;

  if (f->er==er_table) {
    f->rule[i = att2indx(f, att)] = class;
    f->rtp[i] = rtp;
    f->usage[i] = 0; f->imp[i] = 1.;
  }
  else if (f->er==er_list) {
    if (class == CUNDEF) return;
    rl = (rule_list *) malloc(sizeof(*rl));
    if (rl==NULL) {printf("234 memory\n"); exit(0);}
    rl->att = c_vector(f->n_in);
    for (i=0; i<f->n_in; i++) rl->att[i] = att[i];
    rl->class = class; rl->rtp = rtp; rl->usage = 0; rl->imp = 1.;
    rl->n = att2indx(f, att);
    rl->next = f->lrule;
    f->lrule = rl;
    rl->dist = d_vector_ini(f->out->ndesc, 0.);
    rl->dist[class] = 1.;
  }
}

void set_rule_ndis(fams *f, char *att, char class, rule_type rtp, double ndis)
{
  int i;
  rule_list *rl;

  if (f->er==er_table) {
    f->rule[i = att2indx(f, att)] = class;
    f->rtp[i] = rtp;
    f->usage[i] = 0; f->imp[i] = 1.;
  }
  else if (f->er==er_list) {
    if (class == CUNDEF) return;
    rl = (rule_list *) malloc(sizeof(*rl));
    if (rl==NULL) {printf("234 memory\n"); exit(0);}
    rl->att = c_vector(f->n_in);
    for (i=0; i<f->n_in; i++) rl->att[i] = att[i];
    rl->class = class; rl->rtp = rtp; rl->usage = 0; rl->imp = 1.;
    rl->n = att2indx(f, att);
    rl->next = f->lrule;
    f->lrule = rl;
    rl->dist = d_vector_ini(f->out->ndesc, 0.);
    rl->dist[class] = ndis;
  }
}

void set_drule(fams *f, char *att, char class, rule_type rtp, double *dist)
{
  int i;
  rule_list *rl;

  rl = (rule_list *) malloc(sizeof(*rl));
  if (rl==NULL) {printf("234 memory\n"); exit(0);}
  rl->att = c_vector(f->n_in);
  for (i=0; i<f->n_in; i++) rl->att[i] = att[i];
  rl->class = class; rl->rtp = rtp; rl->usage = 0; rl->imp = 1.;
  rl->n = att2indx(f, att);
  rl->next = f->lrule;
  f->lrule = rl;
  rl->dist = d_vector(f->out->ndesc);
  memcpy(rl->dist, dist, sizeof(double) * f->out->ndesc);
}

void change_rule(fams *f, char *att, char class, rule_type rtp)
{
  int i;
  rule_list *rl;

  if (f->er==er_table) {
    f->rule[i = att2indx(f, att)] = class;
    f->rtp[i] = rtp;
  }
  else if (f->er==er_list) {
    i = att2indx(f, att);
    for (rl=f->lrule; rl!=NULL; rl=rl->next)
      if (i == rl->n) {
	rl->class = class; rl->rtp = rtp;
	break;
      }
  }
}

void set_nrule(fams *f, int n, char class, rule_type rtp)
{
  rule_list *rl;

  if (f->er==er_table) {
    f->rule[n] = class;
    f->rtp[n] = rtp;
    f->usage[n] = 0;
    f->imp[n] = 1.;
  }
  else if (f->er==er_list && rtp!=undef && class!=CUNDEF) {
/*    printf("%s -> %d\n", f->out->name, class); */
    rl = (rule_list *) malloc(sizeof(*rl));
    rl->att = c_vector(f->n_in);
    indx2att(f, n, rl->att);
    rl->class = class; rl->rtp = rtp; rl->usage = 0; rl->imp = 1.;
    rl->n = n;
    rl->next = f->lrule;
    f->lrule = rl;
    rl->dist = d_vector_ini(f->out->ndesc, 0.);
    if (class>=f->out->ndesc) {
      printf("c[1] error %d\n", class); exit(0);
    }
    rl->dist[class] = 1.;
  }
}

void set_nrule_ndis(fams *f, int n, char class, rule_type rtp, double ndis)
{
  rule_list *rl;

  if (f->er==er_table) {
    f->rule[n] = class;
    f->rtp[n] = rtp;
    f->usage[n] = 0;
    f->imp[n] = 1.;
  }
  else if (f->er==er_list && rtp!=undef && class!=CUNDEF) {
/*    printf("%s -> %d\n", f->out->name, class); */
    rl = (rule_list *) malloc(sizeof(*rl));
    rl->att = c_vector(f->n_in);
    indx2att(f, n, rl->att);
    rl->class = class; rl->rtp = rtp; rl->usage = 0; rl->imp = 1.;
    rl->n = n;
    rl->next = f->lrule;
    f->lrule = rl;
    rl->dist = d_vector_ini(f->out->ndesc, 0.);
    if (class>=f->out->ndesc) {
      printf("c[2] error %d\n", class); exit(0);
    }
    rl->dist[class] = ndis;
  }
}

char get_nrule(fams *f, int n)
{
  rule_list *rl;

  if (f->er==er_table) 
    if (f->rtp[n]!=undef) return f->rule[n]; else return CUNDEF;
  if (f->er==er_list) {
    for (rl=f->lrule; rl!=NULL; rl=rl->next)
      if (n == rl->n) return rl->class;
    return CUNDEF;
  }
}

char get_rule(fams *f, char *att)
{
  rule_list *rl;
  int i;
  char b, class;
  rule_type rtp;

  if (f->er==er_table) return get_nrule(f, att2indx(f, att));
  if (f->er==er_list) {
    if (f->hash_table == NULL || TRUE)
      for (rl=f->lrule; rl!=NULL; rl=rl->next) {
	for (b=TRUE, i=0; i<f->n_in && b; i++)
	  b = att[i] == rl->att[i];
	if (b) return rl->class;
      }
    else {
      return get_hrule(f, att);
    }
  }
  return CUNDEF;
}

char get_nrule_e(fams *f, int n, rule_type *rtp)
{
  rule_list *rl;

  if (f->er==er_table) {
    *rtp =  f->rtp[n];
    return f->rule[n];
  }
  if (f->er==er_list) {
    for (rl=f->lrule; rl!=NULL; rl=rl->next)
      if (n == rl->n) {
	*rtp = rl->rtp;
	return rl->class;
      }
    *rtp = undef;
    return -1;
  }
}

#define SAFE_FREE(x) {if ((x)!=NULL) {free(x); x=NULL;} }

free_lrule(rule_list *l)
{
  rule_list *rl;
  
/*  return;	*/		/* BUG BUG BUG BBBBB */
  if (l==NULL) return;
  for (rl=l, l=l->next; l!=NULL; l=l->next) {
/*    if (use_distribution) free(rl->dist); */
    free(rl->dist);
    free(rl->att);
    free(rl);
    rl = l;
  }
  l = NULL;
}

free_rules(fams *f)
{
  rule_list *rl;

  if (f->er==er_table) {
    SAFE_FREE(f->rule);
    SAFE_FREE(f->rtp);
    SAFE_FREE(f->imp); 
    SAFE_FREE(f->usage);
  }
  else if (f->er==er_list) {
    free_lrule(f->lrule);
    f->lrule = NULL;
  }
}

void allocate_rules(fams *f, int n_rules)
{
  int i;
  f->n_rules = n_rules;
  f->er = encode_rules; 
  if (f->er==er_table) {
    f->rule = c_vector(n_rules);
    f->rtp = (rule_type *) malloc(sizeof(*f->rtp)*n_rules);
    f->imp = d_vector(f->n_rules);
    f->usage = i_vector(f->n_rules);
    for (i=0; i<n_rules; i++) {
      f->rtp[i] = undef;
      f->imp[i] = 1.;
      f->usage[i] = 0;
    }
  }
  f->lrule = NULL;
}

void free_and_allocate_rules(fams *f, int n_rules)
{
  free_rules(f);
  allocate_rules(f, n_rules);
}

/* set_crossindex: when copying, atributes for two fams may have
   different location. tt[att_f] tells for att_f from fam f where it it's
   position for fam t */

void set_crossindex(fams *f, fams *t, char *tt)
{
  int i, j;

  for (i=0; i<f->n_in; i++) {
    for (j=0; j<t->n_in; j++)
      if (f->in[i]==t->in[j]) {
	tt[i]=j;
	break;
      }
      else tt[i]=-1;
  }
}

/* rule scanning. This is used to traverse a list of
   rules. reset_next_rule(fam) initiates the proces, while
   get_next_rule(class, rtp) returns the next rule. */

int nr_index;
rule_list *nr_list;
char *nr_att = NULL;

int anr_index;
rule_list *anr_list;
char *anr_att = NULL;

void reset_next_rule(fams *f)
{
  if (f->er == er_table) {
    FREE(nr_att);
    for (nr_index=0; nr_index<f->n_rules; nr_index++)
      if (f->rtp[nr_index] != undef) break;
    nr_att = c_vector(f->n_in);
  }
  else if (f->er == er_list)
    nr_list = f->lrule;
}

char get_next_rule(fams *f, char *class, rule_type *rtp)
{
  if (f->er == er_table) {
    if (nr_index >= f->n_rules) return FALSE;
    *class = f->rule[nr_index]; *rtp = f->rtp[nr_index];

    for (nr_index++; nr_index < f->n_rules; nr_index++)
      if (f->rtp[nr_index] != undef) break;
  }
  else if (f->er == er_list) {
    if (nr_list == NULL) return FALSE;
    *class = nr_list->class; *rtp = nr_list->rtp;
    nr_list = nr_list->next;
  }
  return TRUE;
}

char get_next_rulea(fams *f, char *class, rule_type *rtp, char **att)
{
  if (f->er == er_table) {
    if (nr_index >= f->n_rules) return FALSE;
    *class = f->rule[nr_index]; *rtp = f->rtp[nr_index];
    indx2att(f, nr_index, nr_att); *att = nr_att;
    for (nr_index++; nr_index<f->n_rules; nr_index++)
      if (f->rtp[nr_index] != undef) break;
  }
  else if (f->er == er_list) {
    if (nr_list == NULL) return FALSE;
    *class = nr_list->class; *rtp = nr_list->rtp; *att = nr_list->att;
    nr_list = nr_list->next;
  }
  return TRUE;
}

void areset_next_rule(fams *f)
{
  if (f->er == er_table) {
    FREE(anr_att);
    for (anr_index=0; anr_index<f->n_rules; anr_index++)
      if (f->rtp[anr_index] != undef) break;
    anr_att = c_vector(f->n_in);
  }
  else if (f->er == er_list)
    anr_list = f->lrule;
}

char aget_next_rule(fams *f, char *class, rule_type *rtp)
{
  if (f->er == er_table) {
    if (anr_index >= f->n_rules) return FALSE;
    *class = f->rule[anr_index]; *rtp = f->rtp[anr_index];

    for (anr_index++; anr_index < f->n_rules; anr_index++)
      if (f->rtp[anr_index] != undef) break;
  }
  else if (f->er == er_list) {
    if (anr_list == NULL) return FALSE;
    *class = anr_list->class; *rtp = anr_list->rtp;
    anr_list = anr_list->next;
  }
  return TRUE;
}

char aget_next_rulea(fams *f, char *class, rule_type *rtp, char **att)
{
  if (f->er == er_table) {
    if (anr_index >= f->n_rules) return FALSE;
    *class = f->rule[anr_index]; *rtp = f->rtp[anr_index];
    indx2att(f, anr_index, anr_att); *att = anr_att;
    for (anr_index++; anr_index<f->n_rules; anr_index++)
      if (f->rtp[anr_index] != undef) break;
  }
  else if (f->er == er_list) {
    if (anr_list == NULL) return FALSE;
    *class = anr_list->class; *rtp = anr_list->rtp; *att = anr_list->att;
    anr_list = anr_list->next;
  }
  return TRUE;
}

/* count_rules: gives the number of non-undefined rules */

int count_rules(fams *f)
{
  int i, k=0;
  rule_list *rl;

  if (f->er == er_table) {
    for (i=0; i<f->n_rules; i++) if (f->rtp[i]!=undef) k++;
  }
  else if (f->er == er_list) {
    for (rl=f->lrule; rl!=NULL; rl=rl->next, k++);
  }
  else { printf("something's wrong\n"); exit(0);}
  return k;
}

double count_rule_items(fams *f)
{
  int i, j;
  rule_list *rl;
  double d=0;

  for (rl=f->lrule; rl!=NULL; rl=rl->next)
    for (j=0; j<f->out->ndesc; j++)
      d += rl->dist[j];
  return d;
}

/* set_apriory: derives apriory information from the list of rules */

void set_apriory(fams *f)
{
  char class;
  rule_type rtp;
  int i;
  rule_list *rl;

  f->out->c_apriory = d_vector_ini(f->out->ndesc, 0.);
  if (use_distribution) {
    for (rl=f->lrule; rl!=NULL; rl=rl->next)
      for (i=0; i<f->out->ndesc; i++) 
	f->out->c_apriory[i] += rl->dist[i];
  }
  else {
    reset_next_rule(f);
    while (get_next_rule(f, &class, &rtp))
      f->out->c_apriory[class]++;
  }
  d_vector_normalize(f->out->c_apriory, f->out->ndesc);
  d_vector_max(f->out->c_apriory, f->out->ndesc, &i);
  f->out->c_default = i;
}

void set_apriory_id(char *id)
{
  fams *f;

  if ((f = find_fam_vartable(id)) == NULL) return;
  set_apriory(f);
}

/****************************************************************************
****************************************************************************/

void add_fam(Str255 fname, Str255 vname, list_of_str *inname)
{
  var_type *vout;
  fams *fam;
  list_of_fams *lf;
  list_of_str *s;
  list_of_vars *v;
  int n, nin, i;  

  vout = find_var(variables, vname);
  if (vout==NULL) vout = add_var(&variables,vname,FALSE,0)->var;
  
  for (fam=vout->famsout; fam!=NULL; fam=fam->next)
    if (!strcmp(fam->name, fname)) {
      printf("error: fam %s exists. to change, delete it first.\n", fname);
      return;
    }

  for (s=inname; s!=NULL; s=s->next) {
    if (find_var(variables, s->str)==NULL) {
      add_var(&variables,s->str,FALSE,0);
    }
    if (!strcmp(s->str, vname)) {
      printf("error: can't do this, immediate recursion found\n");
      return;			/* BBB check for cycles in general */
    }
  }

  fam = (fams *) malloc(sizeof(*fam));
  fam->hash_table = NULL;
  fam->out = vout;
  strcpy(fam->name, fname);

  for (nin=0, s=inname; s!=NULL; nin++, s=s->next);
  fam->n_in = nin;
  fam->in = (var_type **) malloc(sizeof(*fam->in)*nin);
  for (i=0, n=1, s=inname; s!=NULL; i++, s=s->next) {
    fam->in[i] = find_var(variables, s->str);
    n *= fam->in[i]->ndesc;
    lf = (list_of_fams *) malloc(sizeof(*lf));
    lf->fam = fam;
    lf->next = fam->in[i]->famsin;
    fam->in[i]->famsin = lf;
  }

  allocate_rules(fam, n);
  fam->degrees = (double *) malloc(sizeof(*fam->degrees) * vout->ndesc);

  fam->next = vout->famsout;
  vout->famsout = fam;
}

void list_fams()
{
  list_of_vars *lv;
  fams *fam;
  list_of_fams *lf;
  var_type *v;

  for (lv=variables; lv!=NULL; lv=lv->next) {
    v = lv->var;
    printf("variable %s:\n", v->name);
    printf("  fams out: ");
    for (fam=v->famsout; fam!=NULL; fam=fam->next)
      printf("%s ", fam->name);
    printf("\n  fams in: ");
    for (lf=v->famsin; lf!=NULL; lf=lf->next)
      printf("%s ", lf->fam->name);
    printf("\n");
  }
}

fams *find_fam(Str255 vname, Str255 fname) 
{
  fams *f;
  var_type *v;

  v = find_var(variables, vname);
  if (v==NULL) {
    printf("error: variable %s not found\n", vname);
    return NULL;
  }
  
  for (f=v->famsout; f!=NULL; f=f->next)
    if (!strcmp(f->name, fname))
      return f;
  printf("error: fam %s for variable %s not found\n", fname, vname);
  return NULL;
}

void list_rule(fams *f, char class, char *att)
{
  int j;
  
  printf("    ");
  for (j=0; j < f->n_in; j++)
    printf("%s ", att[j]==CUNDEF?"?":
	   strn(f->in[j]->desc[att[j]].name, print_short));
  printf("  %s\n", class==CUNDEF? "?" : 
	 strn(f->out->desc[class].name, print_short));
}

void list_rule_rl(fams *f, rule_list *rl)
{
  int i, j;

  printf("    ");
  for (j=0; j < f->n_in; j++)
    if (print_short)
      printf("%s ", rl->att[j]==CUNDEF?"?":
	     strn(f->in[j]->desc[rl->att[j]].name, print_short));
    else printf("%s\t", f->in[j]->desc[rl->att[j]].name);
    
  if (print_short) {
    printf("  %s", (rl->rtp==undef || rl->class==CUNDEF)? "?" : 
	   strn(f->out->desc[rl->class].name, print_short));
  }
  else {
    printf("\t%s\t", rl->rtp==undef || rl->class==CUNDEF ? "?" :
	   f->out->desc[rl->class].name);
    printf("%s", rl->rtp==undef?"            ":rl->rtp==man?"man         ":
	   rl->rtp==autom?"auto        ":rl->rtp==manfix?"man fixed   ":
	   "auto fixed  ");
  }

  if (use_distribution) {
    printf("  (");
    for (i=0; i<f->out->ndesc; i++) printf("%4.2lf ", rl->dist[i]);
    printf(")");
  }
  printf("\n");
}

void list_rule_table(fams *f, char print_short)
{
  int i, j, n;
  char *att, class;
  rule_type rtp;
  rule_list *r;

  if (!print_short) printf("%s has %d in variables\n", f->name, f->n_in);

  for (j=1, i=0; i<f->n_in; i++)
    j *= f->in[i]->ndesc;
  f->n_rules = j;

  printf("nrules %d\n", count_rules(f));

  printf("(%s, def %s, apri ", f->er==er_table?"table":"list",
	 f->out->desc[f->out->c_default].name);
  d_vector_print(f->out->c_apriory, f->out->ndesc);
  printf(")\n    ");
  for (i=0; i < f->n_in; i++)
    if (print_short) printf("%s ", strn(f->in[i]->name, print_short));
    else printf("%s\t", f->in[i]->name);
  if (print_short) printf("= %s\n", strn(f->out->name, print_short));
  else printf("=\t%s\n", f->out->name);

  if (f->er == er_table) {
    reset_next_rule(f);
    while (get_next_rulea(f, &class, &rtp, &att)) {
      printf("    ");
      list_rule(f, class, att);
    }
  }
  else {
    rule_list *rl;
    for (rl=f->lrule; rl!=NULL; rl=rl->next) list_rule_rl(f, rl);
  }

/*  att = c_vector(f->n_in);
  for (i=0; i < f->n_rules; i++) {
    indx2att(f, i, att);
    class = get_nrule_e(f, i, &rtp);
    printf("    ");
    list_rule(f, class, att);
  }
  FREE(att); */

}

void list_rules(Str255 fname, Str255 vname)
{
  fams *f;

  f = find_fam(vname, fname);
  if (f==NULL) return;
  list_rule_table(f, print_short);
}

void list_rules_vect(fams *f)
{
  char *att;
  int i, j, n;

  printf("%s = ", f->out->name);
  for (i=0; i<f->n_in; i++) printf("%s ", f->in[i]->name);
  printf("\n");

  for (n=1, i=0; i<f->n_in; i++) n += f->in[i]->ndesc;
  att = c_vector(f->n_in);
  for (i=0; i<f->n_in; i++) {
    for (j=0; j<n; j++) {
      indx2att(f, j, att);
      printf("%c", f->in[i]->desc[att[i]].name[0]);
    }
    printf("\n");
  }

  printf("\n");
  for (i=0; i<f->n_rules; i++) {
    j = get_nrule(f, i);
    if (j==CUNDEF) printf("-");
    else 
      printf("%c", f->out->desc[j].name[0]);
  }
  printf("\n");
  FREE(att);
}

void list_apriory(fams *f)
{
  int i, j;

  printf("(%s, def %s, apri ", f->er==er_table?"table":"list",
	 f->out->desc[f->out->c_default].name);
  for (i=0; i<f->out->ndesc; i++) {
    printf("%s/%6.3lf, ", f->out->desc[i].name,
	   f->out->c_apriory[i]);
  }
  printf(")\n    ");
}

void plot_rules(Str255 fname, Str255 vname)
{
  fams *f;
  var_type *v;
  int i, j, n;
  char fn[] = "000.rules";
  FILE *out;

  f = find_fam(vname, fname);
  if (f==NULL) return;
  v = find_var(variables, vname);

  out = fopen(fn, "w");
  if (f==NULL) {
    printf("error: can't open temporary %s file\n", fn);
    return;
  }

  fprintf(out, "%s\n%d\n", v->name, v->ndesc);
  for (i=0; i<v->ndesc; i++)
    fprintf(out, "%s\n", v->desc[i].name);

  fprintf(out, "%d\n", f->n_in);
  for (i=0; i<f->n_in; i++) {
    v = f->in[i];
    fprintf(out, "%s\n%d\n", v->name, v->ndesc);
    for (j=0; j<v->ndesc; j++)
      fprintf(out, "%s\n", v->desc[j].name);
  }

  fprintf(out, "%d\n", f->n_rules);
  for (i=0; i<f->n_rules; i++)
    fprintf(out, "%d\n", f->rtp[i]!=undef ? f->rule[i] : -1);
  fclose(out);
}

void add_rule_num(Str255 fname, Str255 vname, int n, Str255 qval)
{
  fams *f;
  int i;

  f = find_fam(vname, fname);
  if (f==NULL) return;

  if (n<0 || n>=f->n_rules) {
    printf("error: illegal rule number %d\n", n);
    return;
  }

  for (i=0; i<f->out->ndesc && strcmp(f->out->desc[i].name, qval); i++);
  if (i<f->out->ndesc)
    set_nrule(f, n, i, man);
  else
    printf("error: %d not a legal value of %s\n", qval, vname);
}

/* find_rule: returns the rule number according to the FAM and out
   variable name and list of descriptors */

int find_rule(Str255 fname, Str255 vname, list_of_str **s)
{
  fams *f;
  int i, j, n, nn=0;
  char *att;

  f = find_fam(vname, fname);
  if (f==NULL) {
    printf("error: fam %s not found\n", fname);
    return 0;
  }
  att = c_vector(f->n_in);
  for (i=0; *s!=NULL && i<f->n_in; *s=(*s)->next, i++) {
    for (j=0; j<f->in[i]->ndesc && strcmp(f->in[i]->desc[j].name,(*s)->str);
	 j++);
    if (j<f->in[i]->ndesc) att[i] = j;
    else {
      printf("error: %s is not legal qval for %s\n", (*s)->str,f->in[i]->name);
      return(-1);
    }
  }
  nn = att2indx(f, att);
  free(att);
  return nn;
}

/* find_rule_att_str: set the att for the string of descriptors of rule */

char find_rule_att_str(fams *f, char **att, list_of_str **s)
{
  int i, j, n, nn=0;

  *att = c_vector(f->n_in);
  for (i=0; *s!=NULL && i<f->n_in; *s=(*s)->next, i++) {
    for (j=0; j<f->in[i]->ndesc && strcmp(f->in[i]->desc[j].name,(*s)->str);
	 j++);
    if (j<f->in[i]->ndesc) *(*att+i) = j;
    else {
      printf("error: %s is not legal qval for %s\n", (*s)->str,f->in[i]->name);
      return FALSE;
    }
  }
  return TRUE;
}

void add_rule_str(Str255 fname, Str255 vname, list_of_str *s)
{
  fams *f;
  int i, j, n;
  char *att;

  f = find_fam(vname, fname);
  if (f==NULL) return;
  find_rule_att_str(f, &att, &s);

  if (s!=NULL && s->next==NULL) {
    for (j=0; j<f->out->ndesc && strcmp(f->out->desc[j].name, s->str); j++);
    if (j<f->out->ndesc) set_rule(f, att, j, man);
    else {
      printf("error: %s not legal qval for %s\n", s->str, f->out->name);
      return;
    }    
  }
  else
    printf("error: illegal rule\n");
}

void del_rules(fams *f)
{
  int i;
  if (f->er = er_table)
    for (i=0; i<f->n_rules; i++) f->rtp[i] = undef;
  else if (f->er = er_list) {
    free_lrule(f->lrule);
    f->lrule = NULL;
  }
}

void del_rule_str(Str255 fname, Str255 vname, list_of_str *s)
{
  fams *f;
  int i, j, n, nn=0;

  f = find_fam(vname, fname);
  if (f==NULL) return;
  nn = find_rule(fname, vname, &s);

  if (s==NULL && nn >= 0)
    f->rtp[nn] = undef;
  else
    printf("error: illegal rule\n");
}

void fix_rule(fams *f, int n, char fix)
{
  if (fix) {
    if (f->rtp[n] == man) f->rtp[n] = manfix;
    else if (f->rtp[n] == autom) f->rtp[n] = autofix;
  }
  else {
    if (f->rtp[n] == manfix) f->rtp[n] = man;
    else if (f->rtp[n] == autofix) f->rtp[n] = autom;
  }
}

void fix_rule_str(Str255 fname, Str255 vname, list_of_str *s, char fix)
{
  fams *f;
  int n;

  f = find_fam(vname, fname);
  if (f==NULL) return;
  n = find_rule(fname, vname, &s);

  if (s==NULL && n >= 0) 
    fix_rule(f, n, fix);
  else
    printf("error: illegal rule\n");
}

void fix_rule_num(Str255 fname, Str255 vname, int n, char fix)
{
  fams *f;

  f = find_fam(vname, fname);
  if (f==NULL) return;

  if (n<0 || n>=cfam->n_rules)
    printf("error: illegal rule number %d\n", n);
  else 
    fix_rule(f, n, fix);
}


/* list_informativity: computes the informativity for the node nname
   and for its children, and prints the information out */

void list_informativity(Str255 idname)
{
  fams *f;
  char deb = deb_dec;

  deb_dec = 3;
  if ((f = find_fam_vartable(idname)) == NULL) return;
  informativity_one(f);
/*  informativity_two(f); */
  deb_dec = deb;
}

/* reset_rule_usage: reset the usage field for rules of fams in the
   tree of variable vname */

void reset_usage(var_type *v) 
{
  int i;
  fams *f;

  if (v->mark) return;
  v->mark = TRUE;
  for (f=v->famsout; f!=NULL; f=f->next) {
    for (i=0; i<f->n_rules; i++) {
      f->usage[i]=0.0;
    }
 
    for (i=0; i<f->n_in; i++)
      reset_usage(f->in[i]);
  }
}

void reset_rule_usage(Str255 vname)
{
  var_type *c;			/* class */

  c = find_var(variables, vname);
  if (c==NULL) {
    printf("error: variable %s not found\n", vname);
    return;
  }

  if (c->famsout == NULL) {
    printf("error: %s is a leaf\n", c->name);
    return;
  }

  mark_var(c, FALSE);
  reset_usage(c);
}

/****************************************************************************
TABLES

Tables are used to store specific tables of rules. The structure of
the table is that of fam, but eventually only fields in, rtp, rule are
used. Tables are useful when we want to experiment with rule/structure
formation, and we want to compore results to original table.

Following rutines provide means of add, listing and deleting of tables.
****************************************************************************/

list_of_fams *find_tablel(Str255 name)
{
  list_of_fams *f;

  for (f=tables; f!=NULL; f=f->next)
    if (!strcmp(f->fam->name, name))
      return f;
  return NULL;
}

fams *find_table(Str255 name)
{
  list_of_fams *f;

  for (f=tables; f!=NULL; f=f->next)
    if (!strcmp(f->fam->name, name))
      return f->fam;
  return NULL;
}

void list_tables()
{
  list_of_fams *lf;
  int i;

  if (tables==NULL) {printf("No tables defined.\n"); return;}
  printf("Tables:\n");
  for (lf=tables; lf!=NULL; lf=lf->next) {
    printf("%s: ", lf->fam->name);
    for (i=0; i<lf->fam->n_in; i++)
      printf("%s ", lf->fam->in[i]->name);
    printf("\n");
  }
}

/****************************************************************************
COPY AND SPLIT OF RULE TABLES
****************************************************************************/

fams *create_table(fams *f, char *tname)
{
  fams *t;
  list_of_fams *tab;
  int i;

  tab = (list_of_fams *) malloc(sizeof(*tab));
  t = (fams *) malloc(sizeof(*t));
  t->hash_table = NULL;
  strcpy(t->name, tname);
  tab->fam = t;
  t->out = f->out;
  t->n_in = f->n_in;
  t->in = (var_type **) malloc(sizeof(*t->in)*t->n_in);
  for (i=0; i<f->n_in; i++) t->in[i] = f->in[i];
  
  allocate_rules(t, 0);
  
  tab->next = tables;
  if (tables!=NULL) tables->prev = tab;
  tab->prev = NULL;
  tables = tab;
  return tables->fam;
}

void copy_rules(fams *f, fams *t)
{
  int i, j;
  rule_list *rl, *rl1;
  rule_type rtp;
  char *tt, class;
  char *att, *att1;

  tt = c_vector(f->n_in);
  set_crossindex(f, t, tt);
/*printf("zzz "); for (i=0;i<f->n_in; i++)printf("%d ", tt[i]); printf("\n");*/

  i = t->n_rules;
  del_rules(t);
  allocate_rules(t, i);
  att1 = c_vector(t->n_in);
  reset_next_rule(f);
  while (get_next_rulea(f, &class, &rtp, &att)) {
    for (i=0; i<f->n_in; i++) if (tt[i]>=0) att1[tt[i]] = att[i];
    set_rule(t, att1, class, rtp);
  }
  FREE(tt); FREE(att1);
}

void copy_to_table(Str255 vname, Str255 tname)
{
  var_type *v;
  fams *t;

  v = find_var(variables, vname);
  if (v==NULL) {
    printf("error: variable %s not found\n", vname);
    return;
  }
  if (v->famsout == NULL) {
    printf("error: %s has no rules defined\n", vname);
    return;
  }

  if (find_tablel(tname) == NULL) {
    t = create_table(v->famsout, tname);
    copy_rules(v->famsout, t);
  }
  else
    printf("error: table %s is already defined\n", tname);
}

void copy_sel_rules(fams *f, fams *t, char *sel)
{
  int i, j;
  rule_list *rl, *rl1;
  rule_type rtp;
  char *tt, class;
  char *att, *att1;

  tt = c_vector(t->n_in);
  set_crossindex(f, t, tt);

  i = t->n_rules;
  del_rules(t);
  allocate_rules(t, i);
  att1 = c_vector(t->n_in);
  reset_next_rule(f); j=0;
  while (get_next_rulea(f, &class, &rtp, &att))
    if (sel[j++]) {
      for (i=0; i<f->n_in; i++) if (tt[i]>=0) att1[tt[i]] = att[i];
      set_rule(t, att1, class, rtp);
    }
  FREE(tt); FREE(att1);
}

char *sel = NULL;

/* mark_sel(): preparation for unstratified split */

void mark_sel(int count, double perc)
{
  int i, j, n;

  sel = c_vector_ini(count, FALSE);
  n = count * perc / 100.;
  for (i=0; i<n;) {
    j = count * rnd1e();
    if (!sel[j]) {
      sel[j] = TRUE;
      i++;
    }
  }
/* for (i=0; i<count;i++) printf("%d", sel[i]==TRUE?1:0); printf("\n"); */
  printf("%d %d\n", n, count);
}

/* mark_strat_sel(): preparation for stratified split */

void mark_strat_sel(fams *f, int count, double perc)
{
  int *a;			/* number of instances in each class */
  int *n;			/* how many from each class */
  int nc;			/* number of classes */
  int i, j;
  char *class, c;
  rule_type rtp;
  
  nc = f->out->ndesc;
  a = i_vector(nc);
  n = i_vector(nc);

  reset_next_rule(f);
  class = c_vector(count);
  i=0;
  while (get_next_rule(f, &c, &rtp)) {
    a[c]++;
    class[i++]=c;
  }

  for (i=0; i<nc; i++) n[i]=a[i]*perc/100.;
  i_vector_print(n, nc); printf("\ndone\n");

  sel = c_vector_ini(count, FALSE);
  for (c=0; c<nc; c++) {
    for (i=0; i<n[c];) {
      j = (int) count * rnd1e();
      if (!sel[j] && class[j]==c) {
	sel[j] = TRUE;
	i++;
      }
    }
  }

}

char copy_table(Str255 tname, Str255 vname, double perc)
{
  var_type *v;
  fams *ff, *ft;
  int i, j, count;
  double d;

  if ((v = find_var(variables, vname))==NULL) {
    printf("error: variable %s not found\n", vname);
    return FALSE;
  }
  ft=v->famsout;

  if ((ff=find_table(tname)) == NULL) { 
    printf("error: table %s node defined\n", tname);
    return FALSE;
  }

  if (!(ff->out == v)) {
    printf("error: table %s and variable %s not compatible\n", tname, vname);
    if (ff->out != v) printf("error: %s <> %s\n", ff->out->name, v->name);
    if (ff->n_rules!=ft->n_rules) printf("error: %d <> %d\n", 
					 ff->n_rules, ft->n_rules);
    fprintf(lfile, "ERRO table %s and var %s not compatible to copy\n",
	    tname, vname);
    return FALSE;
  }

  if (perc >= 100) {
    copy_rules(ff, ft);
    count =  count_rules(ff);
  }
  else {
    mark_sel(count = count_rules(ff), perc);
    copy_sel_rules(ff, ft, sel);
    FREE(sel);
  }

  printf("Copy from table %s to %s (%lf)\n", tname, vname, perc);
  fprintf(lfile,"RUCP %lf%% of table %s to %s\n", perc, tname, vname);
  return TRUE;
}

void set_class_distr(var_type *v);
char split_rules(char *id1, char *id2, char *id3, double perc)
{
  fams *f1, *f2, *f3;
  int i, count;
  var_type *v;

  if ((f1 = find_fam_vartable(id1)) == NULL) {
    printf("error: %s not found\n", id1);
    exit(0);
  }
  if ((f2 = find_fam_vartable(id2)) == NULL)
    f2 = create_table(f1, id2);
  if ((f3 = find_fam_vartable(id3)) == NULL)
    f3 = create_table(f1, id3);

/*  mark_sel(count = count_rules(f1), perc); */
  mark_strat_sel(f1, count = count_rules(f1), perc);
  copy_sel_rules(f1, f2, sel);
  for (i=0; i<count; i++) sel[i] = !sel[i];
  copy_sel_rules(f1, f3, sel);
  FREE(sel);

  if (v = find_var(variables,id2)) set_class_distr(v);
  if (v = find_var(variables,id3)) set_class_distr(v);
}

/****************************************************************************
Copy and split of rule tables
using cross-validation
****************************************************************************/

int cv_n_folds;
char *cv_fold_mark;		/* stores index of fold for each example */

void prepare_cv(char *id, int nf)
{
  fams *f;
  int n;
  int *fold_elements;		/* number of examples in each fold */
  char *fold_indx;
  int i, k, l;
  int j;			/* counts the number of folds free */

  cv_n_folds = nf;

  if ((f = find_fam_vartable(id)) == NULL) {
    printf("error: %s not found\n", id);
    exit(0);
  }

  n = count_rules(f);
  
  FREE(cv_fold_mark); cv_fold_mark = c_vector_ini(n, 0);
  fold_elements = i_vector_ini(cv_n_folds, n / cv_n_folds);
  for (i=n-((int)n/cv_n_folds)*cv_n_folds; i>0; i--)
    ++fold_elements[(int) (rnd1e() * cv_n_folds)];
  
  fold_indx = c_vector(cv_n_folds);
  for (i=0; i<cv_n_folds; i++) fold_indx[i] = i;

  j = cv_n_folds;
  for (i=0; i<n; i++) {
    cv_fold_mark[i] = fold_indx[k = (int)(rnd1e() * j)];

    --fold_elements[k];
    if (fold_elements[k] == 0) {
      for (l=k; l<cv_n_folds-1; l++) {
	fold_elements[l] = fold_elements[l+1];
	fold_indx[l] = fold_indx[l+1];
      }
      --j;
    }
  }

/*  i_vector_set(fold_elements, cv_n_folds, 0);
  for (i=0; i<n; i++) fold_elements[cv_fold_mark[i]]++;
  i_vector_print(fold_elements, cv_n_folds);
  printf("\n"); */

  FREE(fold_elements);
}

char split_using_cv(char *id1, char *id2, char *id3, int fold)
{
  fams *f1, *f2, *f3;
  int i, n;
  var_type *v;
  char *s;

  fold--;			/* indices by user are from 1 to 10 */

  if ((f1 = find_fam_vartable(id1)) == NULL) {
    printf("error: %s not found\n", id1);
    exit(0);
  }
  if ((f2 = find_fam_vartable(id2)) == NULL)
    f2 = create_table(f1, id2);
  if ((f3 = find_fam_vartable(id3)) == NULL)
    f3 = create_table(f1, id3);

  n = count_rules(f1);
  s = c_vector(n);

  for (i=0; i<n; i++) s[i] = cv_fold_mark[i] != fold;
  copy_sel_rules(f1, f2, s);

  for (i=0; i<n; i++) s[i] = !s[i];
  copy_sel_rules(f1, f3, s);
  FREE(s);

  if (v = find_var(variables,id3)) set_class_distr(v);
  if (v = find_var(variables,id2)) set_class_distr(v);
}


/****************************************************************************
Rutines for assertion, modification, and deletion of FAMs. 
****************************************************************************/

void del_table(Str255 name)
{
  list_of_fams *lt;

  if ((lt=find_tablel(name))!=NULL) {
    if (lt==tables) {
      if (lt->next != NULL) tables->next->prev = NULL;
      tables = tables->next;
      free(lt->fam);
    }
    else {
      lt->prev->next = lt->next;
      if (lt->next != NULL) lt->next->prev = lt->prev;
    }
    free(lt);
  }
  else
    printf("error: table %s not known\n", name);
}

void list_table(Str255 name)
{
  fams *f;
  if ((f=find_table(name)) != NULL) { 
    printf("listing table %s\n", f->name);
    list_rule_table(f, print_short);
  }
  else
    printf("error: table %s node defined\n", name);
}

/* copy rules from table to fam */

/* add noise to rules */

void add_noise(Str255 name, double perc)
{
  fams *f;
  char *sel;
  int i, j, n, count, nvals;
  rule_list *rl;

  if ((f = find_fam_vartable(name)) == NULL) return;

  /* make a selection array for the rules to change  */
  sel = c_vector_ini(count = count_rules(f), FALSE);
  n = count * perc / 100.;
  for (i=0; i<n;) {
    j = count * rnd1e();
    if (!sel[j]) {
      sel[j] = TRUE;
      i++;
    }
  }  

  nvals = f->out->ndesc;
  for (i=0, rl=f->lrule; rl != NULL; i++, rl=rl->next) 
    if (sel[i]) {
/*      while ((j = nvals * rnd1e()) == rl->class); */
      j = nvals * rnd1e();
      rl->class = j;
      d_vector_set(rl->dist, nvals, 0);
      rl->dist[j] = 1.;
    }
  FREE(sel);
}

/****************************************************************************
Compare table to what is derived from rule tables
****************************************************************************/

char evaltype = 1;

void compare_table_var(Str255 tname, Str255 vname)
{
  var_type *v;
  fams *f;
  int j, k, l;
  char *att, class;
  rule_type rtp;
  int n=0, ndiff=0, ndefaults=0, ndefaults_right=0;
  extern char default_used;

  FIND_VAR(v, vname);
  FIND_TABLE(f, tname);

  if (deb_dec>4) {
    for (j=0; j<f->n_in; j++) printf("%s ", f->in[j]->name);
    printf("-> %s\n", f->out->name);
  }

  build_hash_table(v->famsout);

  reset_next_rule(f);
  while (get_next_rulea(f, &class, &rtp, &att)) {
    n++;
    if (n % 100 == 0) {printf("."); fflush(stdout);}
    if (n % 5000 == 0) {printf("\n");}

    for (j=0; j<f->n_in; j++) {
      for (k=0; k<f->in[j]->ndesc; k++) f->in[j]->opt[k] = 0.;
      if (att[j]!=CUNDEF) 
	f->in[j]->opt[att[j]] = 1.;
      f->in[j]->class = att[j];
    }
    if (evaltype==0) {
      mark_var(v, FALSE);
      default_used = FALSE;
      derive_var(v);
    }
    else {
      default_used = FALSE;
      nderive_var(v);
    }
    if (default_used) {
      ndefaults++;
      if (v->class == class) ndefaults_right++;
    }
/*    if (v->class == CUNDEF) printf("zz undefined\n"); BBBB VERY DANGEROUS */
    if (v->class != class) {
      ndiff++;
/*      printf("%lf, %lf %lf\n", v->opt[class], v->opt[0], v->opt[1]); */

      if (deb_dec>4) {
	for (j=0; j<f->n_in; j++) 
	  printf("%s ", f->in[j]->desc[att[j]].name);
	if (v->class != CUNDEF)
	  printf("should be %s, but is %s\n", 
		 f->out->desc[class].name, f->out->desc[v->class].name);
	else
	  printf("should be %s, but is unknown\n", 
		 f->out->desc[class].name);
      }
    }
  }

  printf("\n");
  fprintf(lfile,"RUCM %s vs table %s, %d different from %d (%6.2lf%%)\n",
	  vname, tname, ndiff, n, (1.-ndiff/(double)n)*100.);
  printf("of %d different %d (%6.2lf%%) def %d(right %d)\n", n, ndiff, 
	 100.*ndiff/(double)n, ndefaults, ndefaults_right);
}

/****************************************************************************
Test of decomposition technique via comparison of derived rules and
table of rules stored previously. Test executes the following steps:

for i=lo to hi step s do
  for j=1 to repeat_tests do
    compose vname
    del rule
    copy i % of rules from tname to vname
    decompose vname
    compare vname to table tname

****************************************************************************/

void test_decomposition(Str255 vname, Str255 tname, int from, int to, int step)
{
  var_type *v;
  fams *f;
  int i, j;

  FIND_VAR(v, vname);
  FIND_TABLE(f, tname);

  fprintf(lfile, "COMM test decomposition %s, table %s, (%d,%d,%d), rep %d\n",
	  vname, tname, from, to, step, repeat_tests);
  if (save_tests) {
    log_tests = fopen("log_c45", "w"); 
    do_save_tests = TRUE;
    in_order = f->in;
  }

  for (i=from; i<=to; i+=step) 
    for (j=0; j<repeat_tests; j++) {
      compose(vname);
      v = find_var(variables, vname);
      del_rules(v->famsout);
      if (!copy_table(tname, vname, (double)i)) return;
      bottom_up_decompose(vname, TRUE, FALSE); 
      compare_table_var(tname, vname);
    }
  if (save_tests) {fclose(log_tests); do_save_tests = FALSE; }
}

/****************************************************************************
Count rules
****************************************************************************/

void count_id_rules(char *idname)
{
  fams *f;
  fams *find_fam_vartable(char *);

  if ((f = find_fam_vartable(idname)) == NULL) return;
  printf("Rules: %d\n", count_rules(f));
  if (use_distribution) printf("Items: %5.1lf\n", count_rule_items(f));
}

/****************************************************************************
Sort rules
****************************************************************************/

extern rule_list **rules_l;	/* sorted rules to learn and to test */
int nn;

int rule_cmp_att_s(rule_list **a, rule_list **b)
{
  int i, j;
  rule_list *r1 = *a, *r2 = *b;

  for (i=0; i<nn; i++) 
    if (r1->att[i] < r2->att[i])
      return -1;
    else if (r1->att[i] > r2->att[i])
      return 1;
  return 0;
}

void sort_rules(fams *f)
{
  rule_list *rl, *rl1, *prl;
  int i;
  int n_rules = count_rules(f);
  
  rules_l = (rule_list **) malloc(n_rules * sizeof(*rules_l));
  for (i=0, rl=f->lrule; rl!=NULL; i++, rl=rl->next)
    rules_l[i] = rl;

  nn = f->n_in;

/*  printf("looking for\n");
  for (i=0; i<n_rules; i++)
    list_rule(f, rules_l[i]->class, rules_l[i]->att); */

  qsort(rules_l, n_rules, sizeof(rules_l), rule_cmp_att_s);

/*  printf("found for\n");
  for (i=0; i<n_rules; i++)
    list_rule(f, rules_l[i]->class, rules_l[i]->att);  */

  rl = NULL;
  for (i=0; i<n_rules; i++) {
    rl1 = (rule_list *) malloc(sizeof(*rl1));
    rl1->att = rules_l[i]->att;
    rl1->class = rules_l[i]->class;
    rl1->dist = rules_l[i]->dist;

    rl1->rtp = rules_l[i]->rtp;
    rl1->imp = rules_l[i]->imp;
    rl1->n = rules_l[i]->n;
    rl1->dist = rules_l[i]->dist;

    rl1->prev = rl;
    rl = rl1;
  }


  prl = NULL;
  for (rl1=rl; rl1!=NULL; rl1=rl1->prev) {
    rl1->next = prl;
    prl = rl1;
  }
  
  f->lrule = prl;
  FREE(rules_l);
}

/****************************************************************************
Testing for duplicates
****************************************************************************/

void test_duplicates(char *id)
{
  fams *f;
  rule_list *rl;
  int n = 0, clash = 0;

  if ((f = find_fam_vartable(id)) == NULL) return;
  sort_rules(f);
  rl=f->lrule;
  if (rl!=NULL && rl->next!=NULL) 
    for (; rl->next!=NULL; rl=rl->next) {
      if (!rule_cmp_att_s(&rl, &(rl->next))) {
	list_rule(f, rl->class, rl->att);
	list_rule(f, rl->next->class, rl->next->att);
	printf("---\n");
	n++;
	if (rl->class != rl->next->class) clash++;
      }
    }
  printf("At least %d duplicate pairs, %d clashes\n", n, clash);
}

void merge_duplicates(char *id)
{
  fams *f;
  rule_list *rl, *prl=NULL;
  int i, n = 0, clash = 0;

  if ((f = find_fam_vartable(id)) == NULL) return;
  sort_rules(f);
  for (rl=f->lrule; rl->next!=NULL; rl=rl->next) {
    rl->prev = prl;
    prl = rl;
  }

  rl=f->lrule;
  if (rl!=NULL && rl->next!=NULL) 
    for (; !((rl == NULL) || (rl->next==NULL)); ) {
      if (!rule_cmp_att_s(&rl, &(rl->next))) {
	for (i=0; i<f->out->ndesc; i++)
	  rl->dist[i] += rl->next->dist[i];
	rl->next = rl->next->next; 
/*	list_rule(f, rl->class, rl->att); */
	n++;
      }
      else rl=rl->next; 
    }
  printf("%d rules merged\n", n);
  for (rl=f->lrule; rl->next!=NULL; rl=rl->next) {
    d_vector_max_max(rl->dist, f->out->c_apriory, f->out->ndesc, &i);
    rl->class = i;
  }
}

/****************************************************************************
Set the class distribution. Based on ->dist[]
****************************************************************************/

void set_class_distr(var_type *v)
{
  rule_list *rl;
  int j;
  fams *f = v->famsout;
  double *dist = d_vector_ini(v->ndesc,0);

  if (f==NULL) return;
  for (rl=f->lrule; rl!=NULL; rl=rl->next)
    for (j=0; j<v->ndesc; j++) dist[j] += rl->dist[j];
  d_vector_normalize(dist, v->ndesc);
  d_vector_max(dist, v->ndesc, &j);
  v->c_default = j;
  for (j=0; j<v->ndesc; j++) v->c_apriory[j] = dist[j];

  FREE(dist);
}

/****************************************************************************
Expand rule (with unknown values)
****************************************************************************/

void set_apriory_desc(fams *f)
{
  rule_list *rl;
  var_type *v;
  int i, j, k;
  double w, ww=0;
  int nclass = f->out->ndesc;

				/* initialize */
  for (i=0; i<f->n_in; i++) {
    v = f->in[i];
    for (j=0; j<v->ndesc; j++) {
      v->desc[j].apriory = 0;
      v->desc[j].cap = d_vector_ini(nclass, 0);
    }
  }

				/* set */
  for (rl=f->lrule; rl!=NULL; rl=rl->next) {
    for (w=0., i=0; i<f->out->ndesc; i++)
      w += rl->dist[i];
    ww += w;
    
    for (i=0; i<f->n_in; i++) 
      if (rl->att[i] != CUNDEF) {
	f->in[i]->desc[rl->att[i]].apriory += w;
	f->in[i]->desc[rl->att[i]].cap[rl->class] += w;
      }
      else {
	for (j=0; j<f->in[i]->ndesc; j++) {
	  f->in[i]->desc[j].apriory += w/(double)f->in[i]->ndesc;
	  f->in[i]->desc[j].cap[rl->class] += w / (double)f->in[i]->ndesc;
	}
      }
  }

  for (i=0; i<f->n_in; i++) {
    v = f->in[i];
    printf("%s: ", v->name);
    for (j=0; j<v->ndesc; j++) {
      v->desc[j].apriory = v->desc[j].apriory/ww;
      for (k=0; k<nclass; k++)
	v->desc[j].cap[k] = v->desc[j].cap[k]/ww;
      printf("%s/%lf ", v->desc[j].name, v->desc[j].apriory);
    }
    printf("\n");
  }
}

double min_exp_w=0.0;		/* minimum weight of the rule to be
				   used for expansion */

void expand_rules(fams *f)
{
  rule_list *rl, *prl, *nrl;
  char *c = c_vector(f->n_in);
  char *att = c_vector(f->n_in);
  char b;
  int i, p, nud, nn=1, j=0;
  double dis;
  int old_class;

  if (f==NULL) return;

  set_apriory_desc(f);
  prl = NULL;
  for (rl=f->lrule; rl!=NULL; rl=rl->next) {
    rl->prev = prl;
    prl = rl;
  }

  for (rl=f->lrule; rl!=NULL; rl=rl->next) {j++;

    for (nud=0, i=0; i<f->n_in; i++) {
      c[i] = rl->att[i]==CUNDEF; 
      if (c[i]) {
	nud++;
	nn *= f->in[i]->ndesc;
      }
    }

   if (nud>0) {		/* delete this rule */
     old_class = rl->class;
      if (rl->prev == NULL) {
	f->lrule = rl->next;
	if (rl->next != NULL)
	  rl->next->prev = NULL;
       }
      else {
	(rl->prev)->next = rl->next;
	if (rl->next != NULL) (rl->next)->prev = rl->prev;
      }

				/* expand the rule */
      for (i=0; i<f->n_in; i++) 
	if (rl->att[i]!=CUNDEF) att[i]=rl->att[i]; else att[i] = 0;
      for (i=f->n_in-1; i>=0; i--)
	if (rl->att[i]==CUNDEF) 
	  {p=i; break;}
    
      while (p>=0) {
/*	for (i=0; i<f->n_in; i++) printf("%d ", att[i]); printf("\n"); */

	/*dis = rl->dist[rl->class]/(double)nn;*/ /* first variant */
	dis = 1; 
	if (dont_care == dc_dont_know)
	  for (i=0; i<f->n_in; i++) {
	    if (rl->att[i]==CUNDEF)
				/* before re Kononko's remark */
/*	      dis *= f->in[i]->desc[att[i]].apriory;  */
	      dis *= f->in[i]->desc[att[i]].cap[old_class]; 
	  }

/*	if (dis>0) set_rule_ndis(f, att, rl->class, rl->rtp, dis); */
	if (dis>min_exp_w) set_rule_ndis(f, att, rl->class, rl->rtp, dis);
	while (p>=0) {
	  if (att[p] < f->in[p]->ndesc-1) {
	    att[p]++;
	    for (i=p+1; i<f->n_in; i++)
	      if (rl->att[i]==CUNDEF) { p = i; att[i]=0; }
	    break;
	  }
	  for (p--; p>=0; p--)
	    if (rl->att[p]==CUNDEF)
	      break;
	}
      }

    }
  }
  free(c); free(att);
}
