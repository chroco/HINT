/****************************************************************************
desc.c

This file holds the rutines that create and modify the describtion of
a variables. I.e., changing the qualitative description, and creating
and changing membership functions of these qualitative terms.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#define GL extern
#include "sds.h"
#include <string.h>


/* add_qdesc_var: assignes a qualitative description (list of strings
   s) to a variable named name */

void add_qdesc_var(Str255 name, list_of_str *sroot)
{
  var_type *v, *tv;
  int old_ndesc, n, i;
  list_of_str *s;
  list_of_fams *lf;
  fams *f;

  v = find_var(variables, name);

  if (v==NULL)
    v = add_var(&variables, name, FALSE,0)->var;
  else {
    FREE(v->desc);
    FREE(v->opt);
    old_ndesc = v->ndesc;
  }

  /* BBB free sroot */

  for (n=0, s=sroot; s!= NULL; s=s->next, n++);
  add_var_ndesc(v, n);

  for (i=0, s=sroot; s!= NULL; s=s->next, i++)
    strcpy(v->desc[i].name, s->str);

  /* have to reset the FAMs this variable was in, if the number of
     qualitative descriptors changed */

				/* BBBBB this needs a complete redifinition */
  if (old_ndesc!=v->ndesc) {
    if (v->famsout!=NULL) {
      f=v->famsout;
      del_rules(f);
      FREE(f->degrees);
      f->degrees = (double *) malloc(sizeof(*f->degrees) * v->ndesc);
    }
    for (lf=v->famsin; lf!=NULL; lf=lf->next) {
      f = lf->fam;
      for (i=0, n=1; i<f->n_in; i++) n *= f->in[i]->ndesc;
      f->n_rules = n;
      free_rules(f);
      allocate_rules(f, n);
/*      printf("%d\n", n); print this shit BBBBB */
    }
  }
}

void add_qdesc_varlist(list_of_str *s1, list_of_str *sroot)
{
  list_of_str *s=s1;
  for (; s!=NULL; s=s->next) 
    add_qdesc_var(s->str, sroot);
  free_str_list(s1);
}

/* list_var_desc: lists all the variables and their descriptions. */

void list_var_desc(list_of_vars *root)
{
  list_of_vars *v;
  desc_t *d;
  int i;

  for (v=root; v!=NULL; v=v->next) {
    printf("   %s: ", strn(v->var->name, print_short));
    if (v->var->ctype == ct_nominal) {
      for (i=0; i<v->var->ndesc; i++) {
	d = &(v->var->desc[i]);
	printf("%s ", strn(d->name, print_short));
      }
      printf("\n");
    }
    else {
      printf("\n      ");
      for (i=0; i<v->var->ndesc; i++) {
	d = &(v->var->desc[i]);
	printf("%s ", strn(d->name, print_short));
	
	if (eval_method == e_fuzzy) {
	  switch (d->tp) {
	  case left:
	    printf("[%5.2lf,%5.2lf,%5.2lf)\n", d->a, d->b,d->c);
	    break;
	  case right:
	    printf("(%5.2lf,%5.2lf,%5.2lf]\n", d->a, d->b,d->c);
	    break;
	  case regular:
	    printf("(%5.2lf,%5.2lf,%5.2lf,%5.2lf)\n", d->a, d->b, d->c, d->d);
	    break;
	  case none:
	    printf("none\n");
	    break;
	  }
	}
	else if (eval_method == e_interval) {
	  printf("%5.2lf .. %5.2lf\n      ", d->start, d->start+d->delta);
	}
      }
    }
    printf("\n");
  }
}

/* free_list_of_vars: frees mem space occupied by a list of
   variables. If dest, then frees variables pointed to as well, else
   it frees just a list (usual case). */

void free_list_of_vars(list_of_vars *lv, char dest)
{
  list_of_vars *tmp;
  
  for (; lv != NULL; lv = tmp) {
    tmp = lv->next;
    /* BBB if (dest) ... free variable desc */
    free(lv);
  }
}

/* rename_desc: rename descriptor a of vname to b */

void rename_desc(Str255 vname, Str255 a, Str255 b)
{
  var_type *v;
  int i,j,k;
  char bb=FALSE;

  printf("rename %s of %s to %s\n", a, vname, b);

  FIND_VAR(v, vname);
  for (i=0; i<v->ndesc; i++) 
    if (!strcmp(b, v->desc[i].name)) {
      printf("error: descriptor %s of %s is already defined\n", b, v->name);
      return;
    } 

  for (i=0; i<v->ndesc && !bb; i++) bb = !strcmp(a, v->desc[i].name);
  if (bb) {
    strcpy(v->desc[--i].name, b);
  }
  else
    printf("error: no descriptor %s of %s \n", a, v->name);
}

/****************************************************************************
Interval logic stuff
****************************************************************************/

/* add_idesc_var: adds intervals info to a variable. If descriptors
   are not named yet, it names them, else, it just attaches the
   interval info. */

void add_idesc_var(Str255 vname, list_of_num *ln)
{
  int i, n;
  list_of_num *l;
  var_type *v;

  v = find_var(variables, vname);
  if (v==NULL) {
    printf("error: variable %s not defined\n", vname);
    return;
  }

  for (n=0, l=ln; l!=NULL; n++, l=l->next);

  if (v->ndesc==0) {
    v->desc = (desc_t *) malloc(sizeof(*v->desc)*(n-1));
    v->ndesc = n-1;
    for (i=0; i<n-1; i++) sprintf(v->desc[i].name, "i%d", i);
  }
  if (v->ndesc==n-1) {
    for (l=ln, i=0; i<v->ndesc; i++, l=l->next) {
      v->desc[i].start = l->num;
      v->desc[i].delta = l->next->num - l->num;
    }
  }
  else
    printf("error: number of descriptors and intervals must match\n");
}

/****************************************************************************
Fuzzy logic stuff
****************************************************************************/

/* add_fdesc_var: adds a fuzzy description to a qualitative descriptor
   dname of variable vname. Rectangular descriptions are assumed. a,
   b, c, d state the four basic points of rectangular. For left or
   right one, just a and b are used. */

void add_fdesc(desc_t *des, trapz_type tp,
	       double a, double b, double c, double d)
{
  int i;
      
  des->a = a; des->b = b; des->c = c; des->d = d; 
  /* computation of trapezoid's slope */
  des->tp = tp;      
  switch (des->tp) {
  case left:
    des->r_slope = 1.0/(b - c);
    des->l_slope = 0.0;
    des->m = (b+c)/2.-a;
    des->cent = (b*b + b*c - 3.*a*a + c*c) / (6. * des->m);
    break;
  case regular:
    des->l_slope = 1.0/(b - a);
    des->r_slope = 1.0/(c - d);
    des->m = (c+d-b-a)/2.;
    des->cent = (-(a*a + a*b + b*b) + c*c + c*d + d*d) / (6. * des->m);
    break;
  case right:
    des->l_slope = 1.0/(b - a);
    des->r_slope = 0.0;
    des->m = c-(a+b)/2.;
    des->cent = (-a*a - a*b - b*b + 3.*c*c) / (6. * des->m);
    break;
  }
}

void add_fdesc_var(Str255 dname, Str255 vname, trapz_type tp,
		   double a, double b, double c, double d)
{
  var_type *v;
  desc_t *des;
  int i;
  
  v = find_var(variables, vname);
  if (v!=NULL) {
    
    for (i=0; i<v->ndesc && strcmp(v->desc[i].name,dname); i++);
    if (i<v->ndesc) {
      /*      printf("found %s for %s\n", dname, vname); */
      add_fdesc(&(v->desc[i]), tp, a, b, c, d);
    }
    else
      printf("error: unknown descriptor %s of variable %s\n", dname, vname);
  }
  else
    printf("error: variable %s not defined\n", vname);
}

/* plot_var_desc: plots out the membership function for the variables
   in given list.  for now, uses latex, which is odd but is still
   works. should shift to xplot, or something standard. */

/* #define N(x) (((x)-min)/(max-min)*100) */
#define N(x) (x)

void plot_var_desc_latex(list_of_vars *root)
{
  double min, max;
  list_of_vars *lv;
  var_type *v;
  desc_t *d;
  int i;
  FILE *f;
  char fname[] = "tmp000";
  char s[100];

  sprintf(s,"%s.tex", fname);
  f = fopen(s, "w");
  if (f==NULL) {
    printf("error: can't open tmp file\n");
    return;
  }

  fprintf(f, "\\documentstyle[piclatex]{article}\n\\pagestyle{empty}\n");
  fprintf(f, "\\begin{document}\n\n");
  for (lv=root; lv!=NULL; lv=lv->next) {
    v = lv->var;
    fprintf(f, "\\begin{figure}\n$$\n\\beginpicture\n");

    min = v->desc[0].a;
    max = v->desc[v->ndesc-1].c;

    fprintf(f, "\\setcoordinatesystem units <%lfmm,4cm>\n", 100*1/max-min);
    fprintf(f, "\\setplotarea x from %lf to %lf, y from 0 to 1\n", min, max);
    fprintf(f, "\\axis bottom ticks numbered out from %lf to %lf by %lf /\n",
	    min, max, (max-min)/4.);
    fprintf(f, "\\axis left ticks out from 0 to 1 by 0.5 /\n");

    for (i=0; i<v->ndesc; i++) {
      d = &(v->desc[i]);
      switch (d->tp) {
      case left:
	fprintf(f, "\\plot %lf 1 %lf 1 %lf 0 /\n", N(d->a), N(d->b), N(d->c));
	break;
      case right:
	fprintf(f, "\\plot %lf 0 %lf 1 %lf 1 /\n", N(d->a), N(d->b), N(d->c));
	break;
      case regular:
	fprintf(f, "\\plot %lf 0 %lf 1 %lf 1 %lf 0 /\n",
		N(d->a), N(d->b), N(d->c), N(d->d));
	break;
      case none:
	break;
      }
    }
    fprintf(f, "\\endpicture\n$$\n");
    fprintf(f, "\\caption{Membership function for %s}\n",v->name);
    fprintf(f, "\\end{figure}\n\n");
  }
  fprintf(f, "\\end{document}\n");
  fclose(f);
  sprintf(s, "latex %s > tmp000.log", fname);
/*   printf("%s\n", s); */
  printf("latex ..."); fflush(stdout);
  system(s);
  sprintf(s, "xdvi %s &", fname);
/*   printf("%s\n", s); */
  printf("\nxdvi ...");
  system(s);
  printf("\n");
}

/* desc_defined: returns TRUE if all descriptors for v have their
   membership functions defined */

char desc_defined(var_type *v) {
  int i;
  for (i=0; i<v->ndesc; i++)
    if (v->desc[i].tp == none)
      return FALSE;
  return TRUE;
}

void plot_var_desc_graph(list_of_vars *root)
{
  list_of_vars *lv;
  var_type *v;
  desc_t *d;
  char defined;			/* is description defined */
  int i, ndes;
  FILE *f;
  char fdname[] = "tmp_plot.data";
  char fpname[] = "tmp_plot.plot";
  Str255 s;
  double p_up_off=0.1, p_down_off=0.1, p_m_off=0.1;
  double p_off, p_width, p_offm;

				/* count the number of descriptions */

  for (ndes=0, lv=root; lv!=NULL; lv=lv->next)
    if (desc_defined(lv->var))
      ndes++;
  
  p_width = ((1.0-p_up_off-p_down_off)-(ndes-1.0)*p_m_off) / ndes;
  p_off = p_up_off;

  sprintf(s, "rm -f %s", fpname);
  system(s);

  for (lv=root; lv!=NULL; lv=lv->next) 
    if (desc_defined(lv->var)) {
      p_off += p_width;
      p_offm = 1 - p_off;

      f = fopen(fdname, "w");
      if (f==NULL) {
	printf("error: can't open temporary %s file\n", fdname);
	return;
      }
      
      v = lv->var;
      
      for (i=0; i<v->ndesc; i++) {
	d = &(v->desc[i]);
	switch (d->tp) {
	case left:
	  fprintf(f, "%lf 0.5 \"%s\"\n", N(d->a)+(N(d->b)-N(d->a))/2.,
		  v->desc[i].name);
	  fprintf(f, "%lf 1\n%lf 1\n%lf 0\n\"\"\n", N(d->a), N(d->b), N(d->c));
	  break;
	case right:
	  fprintf(f, "%lf 0.5 \"%s\"\n", N(d->b)+(N(d->c)-N(d->b))/2., 
		  v->desc[i].name);
	  fprintf(f, "%lf 0\n%lf 1\n%lf 1\n", N(d->a), N(d->b), N(d->c));
	  break;
	case regular:
	  fprintf(f, "%lf 0.5 \"%s\"\n", N(d->b)+(N(d->c)-N(d->b))/2.,
		  v->desc[i].name);
	  fprintf(f, "%lf 0\n%lf 1\n%lf 1\n%lf 0\n\"\"\n",
		  N(d->a), N(d->b), N(d->c), N(d->d));
	  break;
	case none:
	  break;
	}
      }
      fclose(f);
      sprintf(s, "gnugraph -s -u %lf -h %lf -y 0 1.2 -b -L \"Membership functions for %s\" < %s >> %s",
	      p_offm, p_width, v->name, fdname, fpname);
/*      printf("xx %lf %lf\n", p_offm, p_width);
      printf("%s\n", s); */
      p_off += p_m_off;
      system(s);
    }
  sprintf(s, "xplot -f 7x13 -geometry 500x%d < %s", 150*ndes, fpname);
  system(s);
}

void plot_var_desc(list_of_vars *root)
{
  if (plot_type == latex)
    plot_var_desc_latex(root);
  else if (plot_type == graph)
    plot_var_desc_graph(root);
}

