/****************************************************************************
io.c

Rutines that are used for saving the data (structure, fams, options,
descriptions). The data is stored in a script format file, so that
loading is done with use of script interpreter.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define GL extern
#include "sds.h"

#define APPEND_FILE out=fopen(fname,"a"); \
  if (out==NULL) {printf("error: could not open %s\n", fname); \
		    return;}

/* save_struct: saves the variables and the structure in terms of FAM
   definitions */

void save_struct(char *fname)
{
  list_of_vars *lv;
  FILE *out;
  
  printf("callled\n");
  APPEND_FILE;
  fprintf(out,"\n# VARIABLES\n");
  for (lv=variables; lv!=NULL; lv=lv->next)
    fprintf(out,"add var %s\n", lv->var->name);
}
  
void save_fam_info(char *fname)
{
  list_of_vars *lv;
  var_type *v;
  fams *f;
  int i;
  FILE *out;
  
  APPEND_FILE;
  fprintf(out,"\n# DEPENDENCIES\n");
  for (lv=variables; lv!=NULL; lv=lv->next) {
    v = lv->var;
    for (f=v->famsout; f!=NULL; f=f->next) {
/*      fprintf(out,"fam %s for %s is {", f->name, v->name); */
      fprintf(out,"%s depends on {", f->name, v->name);
      for (i=0; i<f->n_in; i++)
	fprintf(out,"%s%s", f->in[i]->name, i!=f->n_in-1 ? ", " : "}\n");
    }
  }
  fclose(out);
}

/* save_des: stores the list of qualitative descriptors for each of
   the variables, together with fuzzy membership functions - if they
   exist */

void save_des(char *fname)
{
  list_of_vars *lv;
  var_type *v;
  desc_t *d;
  int i;
  FILE *out;

  APPEND_FILE;
  fprintf(out,"\n# DESCRIPTORS\n");
  for (lv=variables; lv!=NULL; lv=lv->next) {
    v = lv->var;
    if (v->ndesc <= 0)
      fprintf(out,"# %s none\n", v->name);
    else {
      fprintf(out,"%s in {", v->name);
      for (i=0; i<v->ndesc; i++)
	fprintf(out,"%s%s", v->desc[i].name, i!=v->ndesc-1 ? ", " : "}\n");
    }
  }

/*
  fprintf(out,"\n# MEMBERSHIP FUNCTIONS\n");
  for (lv=variables; lv!=NULL; lv=lv->next) {
    v = lv->var;
    for (i=0; i<v->ndesc; i++) {
      d = &(v->desc[i]);
      if (d->tp!=none) fprintf(out,"%s of %s is ", d->name, v->name);
      switch (d->tp) {
      case left:
	fprintf(out,"[%11.7e,%11.7e,%11.7e)\n", d->a, d->b,d->c);
	break;
      case right:
	fprintf(out,"(%11.7e,%11.7e,%11.7e]\n", d->a, d->b,d->c);
	break;
      case regular:
	fprintf(out,"(%11.7e,%11.7e,%11.7e,%11.7e)\n", d->a, d->b, d->c, d->d);
	break;
      case none:
	fprintf(out,"# %s of %s is none\n", d->name, v->name);
	break;
      }
    }
    fprintf(out,"\n");
  }
*/
  fclose(out);
}

/* save_opt: saves the options. If qualitative descriptor is give,
   this is prefered. */

void save_opt(char *fname)
{
  list_of_vars *lv;
  var_type *v;
  list_of_opt *opt;
  int i, j, n;
  FILE *out;

  APPEND_FILE;
  set_opt("tmp opt");		/* save current option */
  fprintf(out,"\n# OPTIONS\n");

  for (opt=options->next; opt!=NULL; opt=opt->next) {
    
    fprintf(out,"\n# option %s\n", opt->name);
    fprintf(out,"undef opt\n");

    select_opt(opt->name);
    for (lv=variables; lv!=NULL; lv=lv->next) {
      v=lv->var;
      for (i=0, n=0; i<v->ndesc; i++)
	if (v->opt[i] > 0) n++;
      if (n>0) {
	fprintf(out,"%s = {", v->name);
	for (i=0, j=0; i<v->ndesc; i++)
	  if (v->opt[i] > 0)
	    fprintf(out,"%s/%11.7e%s",v->desc[i].name,v->opt[i],
		    j++!=n-1?", ":"}\n");

      }
      if (v->valdef)
	fprintf(out,"%s = %11.7e\n", v->name, v->val);
      if (v->expect_def)
	fprintf(out,"expect %s = %11.7e\n", v->name, v->expect);
    }

    fprintf(out,"set opt %s\n", opt->name);
  }
  select_opt("tmp opt");	/* restore current option */
  del_opt("tmp opt");		/* delete it from the list */
  fclose(out);
}

/* save_rules: save all rules for all fams within a structure */

void save_rules(char *fname)
{
  list_of_vars *lv;
  var_type *v;
  fams *f;
  int i, j, n;
  FILE *out;
  char *att, class;
  rule_type rtp;
  
  APPEND_FILE;
  fprintf(out,"\n# RULES\n");
  for (lv=variables; lv!=NULL; lv=lv->next) {
    v = lv->var;
    for (f=v->famsout; f!=NULL; f=f->next) {
/*      fprintf(out,"\n# rules for fam %s: {", f->name); */
      fprintf(out,"\n# rules for {");
      for (i=0; i<f->n_in; i++)
	fprintf(out,"%s%s", f->in[i]->name, i!=f->n_in-1 ? ", " : "} ");
/*    fprintf(out,"-> %s\nsel fam %s for %s\n", v->name, f->name, v->name); */
      fprintf(out,"-> %s\nsel %s\n", v->name, v->name);

      fprintf(out,"rule table {\n");
      reset_next_rule(f);
      while (get_next_rulea(f, &class, &rtp, &att)) {
	for (j=0; j<f->n_in; j++)
	  if (att[j]==CUNDEF) fprintf(out, "? ");
	  else fprintf(out,"%s ", f->in[j]->desc[att[j]].name);
	fprintf(out, " %s\n", f->out->desc[class].name);
      }
      fprintf(out,"}\n");      
    }
  }
  fclose(out);
}

/* save_all: */

void save_all(char *fname)
{
  Str255 s;
  sprintf(s, "rm %s", fname);
  system(s);
/*  save_struct(fname); */
  save_des(fname);
  save_fam_info(fname);
  save_rules(fname);
  save_opt(fname);
}


/****************************************************************************
Saving to C4.5 format.
To be used also with MLC++ utilities.
****************************************************************************/

fams *find_fam_vartable(char *idname)
{
  fams * f;
  var_type *v;

  v = find_var(variables, idname);
  if ((v = find_var(variables, idname)) != NULL) f=v->famsout;
  else
    if ((f = find_table(idname)) == NULL) {
/*      printf("error: variable or table %s not found\n"); */
      return NULL;
    }
  return f;
}

void do_save_c45_rules(FILE *fout, fams *f, int *x)
{
  char class, *att;
  rule_type rtp;
  int i;

/*  for (i=0; i<f->n_in; i++) printf("%d %d %s\n", i, x[i], f->in[i]->name); */
  att = c_vector(f->n_in);
  reset_next_rule(f);
  while (get_next_rulea(f, &class, &rtp, &att)) {
    for (i=0; i<f->n_in; i++)
      if (att[x[i]]!=CUNDEF)
	fprintf(fout, "%s,", f->in[x[i]]->desc[att[x[i]]].name);
      else 
	fprintf(fout, "?,");
    fprintf(fout, "%s.\n", f->out->desc[class].name);
  }
  FREE(att);
}

void save_c45_rules(char *idname, char *fname)
{
  fams *f;
  FILE *fout;
  Str255 s;
  int *x;
  int i, j;

  printf("Save rules for %s into %s\n", idname, fname);

  if ((f = find_fam_vartable(idname)) == NULL) return;
  if (strstr(fname, ".test")==NULL && strstr(fname, ".data")==NULL)
    sprintf(s, "%s.data", fname);
  else 
    sprintf(s, "%s", fname);

  x = i_vector(f->n_in);
  for (i=0; i<f->n_in; i++) x[i] = i;
  
  fout = fopen(s, "w");
  if (fout==NULL) {printf("error: can't open %s\n", s); return;}
  fprintf(fout, "| Data file for %s\n\n", f->out->name);
  do_save_c45_rules(fout, f, x);
  fclose(fout);
  FREE(x);
}

void save_c45_names(char *idname, char *fname)
{
  Str255 s;
  fams *f;
  FILE *fout;
  int i,j;

  if ((f = find_fam_vartable(idname)) == NULL) return;
  if (strstr(fname, ".names")==NULL) sprintf(s, "%s.names", fname);
  else sprintf(s, "%s", fname);

  fout = fopen(s, "w");
  if (fout==NULL) {printf("error: can't open %s\n", s); return;}

  fprintf(fout, "| Names file for %s\n\n", f->out->name);
  for (i=0; i<f->out->ndesc; i++)
    fprintf(fout,"%s%c ", f->out->desc[i].name, (i+1)==f->out->ndesc?' ':',');
  fprintf(fout, "\n\n");

  for (j=0; j<f->n_in; j++) {
    fprintf(fout, "%s: ", f->in[j]->name);
    for (i=0; i<f->in[j]->ndesc; i++)
      fprintf(fout, "%s%c ", f->in[j]->desc[i].name, 
	      (i+1)==f->in[j]->ndesc?'.':',');
    fprintf(fout, "\n");
  }
  fclose(fout);
}

/****************************************************************************
Saving to dot format.
To be used with dot graphics utility program from AT&T.
****************************************************************************/

void struct_to_dot(FILE *fout, var_type *v)
{
  fams *f;
  int i;

  if (v->mark || v->famsout==NULL || v->famsout->n_in==0) return;
  f = v->famsout;
  if (FALSE) {
    fprintf(fout, "  %c%s/%d%c -> {", '"', v->name, v->ndesc, '"');
    for (i=0; i<f->n_in; i++)
      fprintf(fout, "%c%s/%d%c%s ", '"', f->in[i]->name, f->in[i]->ndesc, '"',
	      (i+1)==f->n_in?"}\n":";"); 
  }
  else {
    fprintf(fout, "  %c%s%c -> {", '"', v->name, '"');
    for (i=0; i<f->n_in; i++)
      fprintf(fout, "%c%s%c%s ", '"', f->in[i]->name, '"',
	      (i+1)==f->n_in?"}\n":";");
  }
  v->mark = TRUE;
  for (i=0; i<f->n_in; i++) {
    struct_to_dot(fout, f->in[i]);
/*    if (f->in[i]->famsout == NULL) {
      if (TRUE) fprintf(fout, "  %s [shape=box];\n", f->in[i]->name);
      else fprintf(fout, "  %c%s/%d%c [shape=box];\n", '"', f->in[i]->name, 
		   f->in[i]->ndesc, '"');
    } */
  } 
}

void save_dot_struct(char *idname, char *fname)
{
  Str255 s;
  var_type *v;
  FILE *fout;

  v = find_var(variables, idname);
  if (v==NULL) {printf("error: variable %s not found\n", idname); return;}
  if (strstr(fname, ".")==NULL) sprintf(s,"%s.dot", fname);
  fout = fopen(s, "w");
  if (fout==NULL) {printf("error: can't open %s\n", s); return;}

  fprintf(fout, "digraph G {\n");
  fprintf(fout, "  edge [dir=back];\n");
  mark_var(v, FALSE);
  struct_to_dot(fout, v);
  fprintf(fout, "}\n");
  fclose(fout);
}

/****************************************************************************
Loading indexed rules directly for the file.
****************************************************************************/

void load_irules(char *fname, fams *f, int offset)
{
  FILE *fin;
  int att_num = 0;
  char *att = c_vector(f->n_in);
  int val, n=0;

  if (f==NULL) { printf("error: no FAM current\n"); return; }

  fin = fopen(fname, "r");
  if (fin==NULL) {printf("error: can't open %s\n", fname); return;}
  
  while (fscanf(fin, "%d", &val)!=EOF) {
    val -= offset;
    if (att_num == f->n_in) {
      if (val<0 || val>=f->out->ndesc) {
	printf("error_ %d: index %d not legal for %s\n", 
	       yylineno, val, f->out->name);
	exit(1);
      }
      n++;
      if (n % 100 == 0) {printf("."); fflush(stdout);}
      if (n % 5000 == 0) {printf(" %7d\n", n);}
      set_rule(f, att, val, man);
      att_num = 0;
    }
    else {
      if (val<0 || val>=f->in[att_num]->ndesc) {
	printf("error_ %d: index %d not legal for %s\n", 
	       yylineno, val, f->in[att_num]->name);
	exit(0);
      }
      else
	att[att_num] = val;
      att_num++;
    }
  }
  printf("\n");
  fclose(fin);
  if (att_num != 0) {
    printf("error: premature end of %s\n");
  }
}

