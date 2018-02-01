/****************************************************************************
struct.c

This file holds the rutines that create and modify the describtion of
a decision structure. I.e., adding and deleting the variables, and
chaning the inter-variables dependecies.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define GL extern
#include "sds.h"

void list_vars(list_of_vars *root);

/* find_var: returns a pointer to a variable named name. This is
   actually a pointer to the entry in the list. To get a pointer to a
   variable record itself, one has to refer to var entry. */

list_of_vars *find_varl(list_of_vars *root, Str255 name)
{
  list_of_vars *v;

  for (v=root; v!=NULL; v=v->next)
    if (!strcmp(v->var->name, name))
      return v;
  return NULL;
}


var_type *find_var(list_of_vars *root, Str255 name)
{
  list_of_vars *v;

  for (v=root; v!=NULL; v=v->next)
    if (!strcmp(v->var->name, name))
      return v->var;
  return NULL;
}

var_type *find_var_num(list_of_vars *root, Str255 name, int *i)
{
  list_of_vars *v;

  for (*i=0, v=root; v!=NULL; (*i)++, v=v->next)
    if (!strcmp(v->var->name, name))
      return v->var;
  return NULL;
}

/* add_var: checks if variable is already in the list of variables. If
   not, adds variable to a list and initialises some other fields. */

void add_var_ndesc(var_type *nv, int ndesc)
{
  int i;
  Str255 s;

  nv->ndesc = ndesc;
  nv->desc = (desc_t *) malloc(sizeof(*nv->desc) * ndesc);
  for (i=0; i<ndesc; i++) {
    sprintf(s, "%d", i);
    strcpy(nv->desc[i].name, s);
    nv->desc[i].tp = none;
  }
  nv->opt = d_vector(ndesc);
  nv->c_apriory = d_vector_ini(ndesc, 0.);
}

list_of_vars *add_var(list_of_vars **root, Str255 name, char b, int ndesc)
{
  list_of_vars *v;

  if ((v=find_varl(*root,name))==NULL) {
    v = (list_of_vars *) malloc(sizeof(*v));
    v->var = (var_type *) malloc(sizeof(*(v->var)));
    strcpy(v->var->name, name);
    v->var->famsin = NULL;
    v->var->famsout = NULL;
    v->var->desc = NULL;
    v->var->valdef = FALSE;
    v->var->expect_def = FALSE;
    if (ndesc>0) add_var_ndesc(v->var, ndesc);

    v->next = (*root);

    if ((*root)!=NULL) (*root)->prev = v;
    v->prev = NULL;
    (*root) = v;
  }
  else if (b)
    printf("error: variable %s already known\n", name);

  return v;
}


/* del_var: deletes a variable from a list, and handles the rules, and
   those depended on and to */

list_of_fams *del_fam_from_list(fams *f, list_of_fams *lf)
{
  list_of_fams *l;

  if (lf->fam == f) {
    if (lf->next != NULL) lf->next->prev = NULL;
    l = lf->next;
    free(lf);
    return l;
  }
  for (l=lf->next; l!=NULL; l=l->next) 
    if (l->fam == f) {
      if (l->next != NULL) l->next->prev = l->prev;
      l->prev->next = l->next;
      free(l);
    }
  return lf;
}

void del_fam(fams *f)
{
  myfree(f->in); myfree(f->rule); myfree(f->rtp);
  myfree(f->usage); myfree(f->degrees); myfree(f->imp);
}

void del_var(list_of_vars **root, Str255 name)
{
  list_of_vars *lv;
  var_type *v;
  fams *f, *fo, *fi, *fr, *ff;	/* v=fo(A,fi(B)) */
  list_of_fams *lf, *mlf, *lfprev;
  int i, j, k, l, n, nn;
  char *att, *att1, *att2;
  char class_i, class_o, b;
  rule_type rtp_i, rtp_o;

  int nis_both;
  char *is_both;		/* which of A is also in B*/
  char *pos_both;		/* where of B is i-th of A, is_both[i] */
  

  if ((lv=find_varl(*root,name))==NULL) {
    printf("error: variable %s not known\n", name);
    return;
  }

  v =lv->var;

				/* CASE 1: variable is at the root */
  if (v->famsin == NULL && v->famsout != NULL) {
    f = v->famsout;
    for (i=0; i<f->n_in; i++)
      f->in[i]->famsin = del_fam_from_list(f, f->in[i]->famsin);
    del_fam(f);
    myfree(v->desc); myfree(v->opt);
  }
				/* CASE 2: variable is a leaf */
  else if (v->famsin != NULL && v->famsout == NULL) {
    for (lf=v->famsin; lf!=NULL; lf=lf->next) {
      f = lf->fam->out->famsout;
      for (i=0; f->in[i] != v; i++); /*printf("xxx %d\n", i); */
      for (; i<f->n_in-1; i++) f->in[i] = f->in[i+1];
      f->n_in--;
      if (f->n_in == 0) {
	del_fam(lf->fam->out->famsout);
	lf->fam->out->famsout = NULL;
      }
      else {
	f->n_rules = f->n_rules / v->ndesc;
	for (i=0; i<f->n_rules; i++) f->rtp[i] = undef;
      }
    }
    myfree(v->desc); myfree(v->opt);
  }
				/* CASE 3: var is an internal node */
				/* this fails for y=f(a,g(b,a)) */
  else if (v->famsin != NULL && v->famsout != NULL) {
          
    for (mlf=v->famsin; mlf!=NULL; mlf=mlf->next) {
      fo = mlf->fam;
      fi = v->famsout;

      is_both = c_vector_ini(fi->n_in, FALSE);
      pos_both = c_vector(fi->n_in);
      nis_both = 0;
      for (i=0; i<fo->n_in; i++)
	for (j=0; j<fi->n_in; j++)
	  if (fo->in[i] == fi->in[j]) {
	    is_both[j] = TRUE;
	    pos_both[j] = i;
	    nis_both++;
	  }

      f = (fams *) malloc(sizeof(*f));
      f->hash_table = NULL;
      f->next = NULL;
      f->n_in = fo->n_in + fi->n_in - 1 - nis_both;
      
      f->out = fo->out;
      strcpy(f->name, fo->out->name);
      f->in = (var_type **) malloc(sizeof(*f->in) * f->n_in);

/* printf("xx %s %s\n", fo->name, v->name);
for (i=0; i<fo->n_in; i++) printf("%s ", fo->in[i]->name); printf("\n"); */
      for (i=0; fo->in[i]!=v; i++) f->in[i] = fo->in[i];
      for (; i<fo->n_in-1; i++) f->in[i] = fo->in[i+1];
      for (j=0, k=0; j<fi->n_in; j++)
	if (!is_both[j]) f->in[i+(k++)] = fi->in[j];
      
/*      printf("IN: %d %s\n", f->n_in, f->in[0]->name);
      for (i=0; i<f->n_in; i++) printf(" %s\n", f->in[i]->name); */
      
      for (n=1, i=0; i<f->n_in; i++) n *= f->in[i]->ndesc;
      allocate_rules(f, n);
      f->degrees = d_vector(f->out->ndesc);
      
      for (n=0; fo->in[n] != v; n++);

      att2 = c_vector(f->n_in+1);
      reset_next_rule(fo);
      while (get_next_rulea(fo, &class_o, &rtp_o, &att)) {
	nn = att[n];
      	areset_next_rule(fi);
	while (aget_next_rulea(fi, &class_i, &rtp_i, &att1)) {
	  for (b=TRUE, i=0; i<fi->n_in && b; i++)
	    if (is_both[i] && att1[i] != att[pos_both[i]]) b = FALSE;

/*	  printf("zz %s %s %s\n", b?"yes":"no ", 
		 (class_i==nn)?"yes":"no ",f->out->name); */

	  if (class_i == nn && b) {
	    for (i=0; i<fo->n_in; i++) att2[i] = att[i];
	    for (k=n; k<fo->n_in; k++) att2[k] = att2[k+1];

	    for (j=0, i=0; i<fi->n_in; i++)
	      if (!is_both[i]) att2[(j++)+fo->n_in-1] = att1[i];
	    set_rule(f, att2, class_o, rtp_o);
	  }
	}
      }
      fr = mlf->fam;
				/* update var on the same level */
      for (i=0; i<fo->n_in; i++) {
	for (lf=fo->in[i]->famsin; lf->fam!=fr; lf=lf->next);
	lf->fam = f;
      } 

				/* update vars beneath the deleted */
      for (i=0; i<fi->n_in; i++) {
/*	printf("%s %s (", fi->in[i]->name, fi->name);
	for (lf=fi->in[i]->famsin; lf!=NULL; lf=lf->next)
	  printf("%s ", lf->fam->name);
	printf(")\n"); */

	for (lfprev=NULL, lf=fi->in[i]->famsin; lf->fam!=fi; lf=lf->next)
	  lfprev = lf;
	if (is_both[i]) {
	  if (lfprev != NULL) lfprev->next = lf->next;
	}
	else
	  lf->fam = f;
      }

      fo->out->famsout = f; /* BBBBB problematicen */
/*      mlf->fam->out->famsout = f; */ /* BBBBB problematicen */
/*      del_fam(fr); */
      FREE(is_both); FREE(pos_both);
    }
  }

				/* del variable from list of variables */
  if (lv==(*root)) {
    if (lv->next != NULL) (*root)->next->prev = NULL;
    (*root) = (*root)->next;
    free(lv->var);
  }
  else {
    lv->prev->next = lv->next;
    if (lv->next != NULL) lv->next->prev = lv->prev;
  }
  free(lv);
}


/* del_list_vars: deletes a list of variables, but leaves the entries
   of variables untouched. */

void del_list_vars(list_of_vars **root)
{
  list_of_vars *v;

  for (; *root != NULL; *root = (*root)->next) {
    v=(*root)->next;
    free(*root);
  }
}

/****************************************************************************
  INPUT / OUTPUT RUTINES
  screen listing
****************************************************************************/

/* list_vars: lists all the variables enetered so far. */

void list_vars(list_of_vars *root)
{
  list_of_vars *v,*vd;

  for (v=root; v!=NULL; v=v->next)
    printf("%s ", v->var->name);
  printf("\n");
}

/* list_struct: takes a list of variables, figures out for which
   variable thre is no up dependecies, and for these variables lists
   out a tree of variables that influence them. */

void list_struct_tree(var_type *v, int n, ltype lt)
{
  int i, nspace;
  list_of_vars *lv;
  char b = TRUE;
  fams *f;

/*  if (v->mark) return; */
  v->mark = TRUE;

  for (i=0; i<n; i++) printf("  ");
/*  printf("%s (%s)/%d", v->name, v->ctype==ct_nominal?"n":"q", v->ndesc); */
  printf("%s/%d", v->name, v->ndesc);

  nspace = DESC_COL - n*2 - strlen(v->name);
  for (i=0; i<nspace; i++) printf(" ");

  switch (lt) {
  case desc:
    for (i=0; i<v->ndesc; i++)
      printf("%s\t", v->desc[i].name);
    break;
  case curropt:
    if (eval_method==e_crisp) {
      for (b=TRUE, i=0; i<v->ndesc && b; i++) {
	b = !(v->opt[i] > 0);
	if (!b) printf("%s", v->desc[i].name);
      }
      if (b) printf("?");
    }
    else {
      if (v->valdef) printf("%6.3lf\t", v->val); else printf("      \t");
      if (v->expect_def) printf("%6.3lf\t", v->expect);else printf("      \t");
      for (i=0; i<v->ndesc; i++)
	if (v->opt[i] > 0) {
	  printf("%s/%4.2lf\t", v->desc[i].name, v->opt[i]);
	  b = FALSE;
	}
      if (b) printf("?");
    }
    break;
  case infor:
    printf("%lf (%lf)", v->infor, v->gain);
    break;
  }

  printf("\n"); 

  for (f=v->famsout; f!=NULL; f=f->next)
    for (i=0; i<f->n_in; i++)
      list_struct_tree(f->in[i], n+1, lt);
}

void list_struct(list_of_vars *root, ltype lt)
{
  list_of_vars *v;

  for (v=root; v!=NULL; v=v->next)
    if (v->var->famsin == NULL) {
      mark_var(v->var, FALSE);
      list_struct_tree(v->var, 0, lt);
    }
}

/****************************************************************************
deriving the difference between two structures
****************************************************************************/

/*
get the leaf variables of A
get the leaf variables of B
subsection
get distance matrix for A
get distance matrix for B
compare
*/

int node_distance(var_type *v1, var_type *v2)
{
  int i, n;
  fams *f;
  list_of_fams *lf;

  v1->mark = FALSE;

  if (v1 == v2) return 0;

  f = v1->famsout;		/* search down */
  if (f!=NULL) 
    for (i=0; i<f->n_in; i++)
      if (f->in[i]->mark) {
	n = node_distance(f->in[i], v2);
	if (n>=0) {
	  if (f->n_in>1) return n+1;
	  else return n;	/* in case of redundancy link */
	}
      }

  for (lf=v1->famsin; lf!=NULL; lf=lf->next) /* search up */
    if (lf->fam->out->mark) {
      n = node_distance(lf->fam->out, v2);
      if (n>=0) return n+1;
    }

  return -1;			/* not found */
}

var_type **leaves = NULL, **leaves1 = NULL;
int **distance = NULL;
int nleaves, nleaves1, tmp_nleaves;

int get_n_leaves(var_type *v, char mark)
{
  int i, n=0;
  fams *f;

  if (mark==-1 || v->mark == mark) {
    if (v->famsout == NULL) return 1;
    else {
      f = v->famsout;
      for (i=0; i<f->n_in; i++)
	n += get_n_leaves(f->in[i], mark);
    }
  }
  return n;
}

int ndes_inodes;
int get_n_inodes(var_type *v, char mark)
{
  int i, n=0;
  fams *f;

  if (mark==-1 || v->mark == mark) {
    if (v->famsout == NULL) return 0;
    else {
      f = v->famsout;
      for (i=0; i<f->n_in; i++)
	n += get_n_inodes(f->in[i], mark);
      printf("yy %s %d in %d\n", v->name, v->ndesc, f->n_in);
      if (f->n_in>1) { ndes_inodes += v->ndesc; return n+1;}
      else return n;
    }
  }
  return 0;
}

void build_leaves_array(var_type *v, var_type **leaves)
{
  int i;

  if (!v->mark) return;
  if (v->mark && v->famsout == NULL) leaves[tmp_nleaves++] = v;
  else
    for (i=0; i<v->famsout->n_in; i++)
      build_leaves_array(v->famsout->in[i], leaves);
}

/* derive_distance_matrix: derives the distance matrix for all vars
   rooted in v that are marked TRUE */

void derive_distance_matrix(var_type *v)
{
  int i, j;
  list_of_vars *lv;

/*printf("FOR %s\n", v->name);*/
  distance = i_matrix(nleaves, nleaves);

  for (lv=variables; lv!=NULL; lv=lv->next) mark_var(lv->var, FALSE);
  mark_var(v, TRUE);

/*  for (i=0; i<nleaves; i++) printf("%s ", leaves[i]->name); printf("\n");*/

  for (i=0; i<nleaves-1; i++) 
    for (j=i+1; j<nleaves; j++) {
      mark_var(v, TRUE);
      distance[i][j] = node_distance(leaves[i], leaves[j]) - 1;
    }

/*  for (i=0; i<nleaves-1; i++) {
    for (j=i+1; j<nleaves; j++)
      printf("%2d ", distance[i][j]);
    printf("\n");	
  } */
}

/* it is expected that id2 is that with a subset of leaves of id1 */

int set_leaves_array(var_type *v, var_type ***leaves)
{
  list_of_vars *lv;

  for (lv=variables; lv!=NULL; lv=lv->next) mark_var(lv->var, FALSE);
  mark_var(v, TRUE);
  nleaves = get_n_leaves(v, TRUE);
  FREE(*leaves); 
  *leaves = (var_type **) malloc(sizeof(**leaves) * nleaves);
  tmp_nleaves = 0; 
  build_leaves_array(v, *leaves);
  return tmp_nleaves;
}

double compare_struct_dist(char *id1, char *id2)
{
  var_type *v1, *v2;
  list_of_vars *lv;
  int n, **distance1, i, j, k;
  double dist;

  FIND_VAR(v1, id1);
  FIND_VAR(v2, id2);
  
  nleaves1 = set_leaves_array(v1, &leaves1);
  nleaves = set_leaves_array(v2, &leaves);
/*  for (i=0; i<nleaves1; i++) printf("%s ", leaves1[i]->name); printf("\n");
  for (i=0; i<nleaves; i++) printf("%s ", leaves[i]->name); printf("\n"); */
  
  derive_distance_matrix(v1);
  distance1 = distance;
  derive_distance_matrix(v2);

  dist = 0.;
  for (i=0; i<nleaves-1; i++) 
    for (j=i+1; j<nleaves; j++)
      dist += ABS(distance1[i][j] - distance[i][j]);
  dist = dist * 2 / (double) (nleaves * (nleaves - 1));
/*  printf("%lf\n", dist); */
  return dist;
}

int count_id_leaves(char *id)
{
  var_type *v;
  list_of_vars *lv;

  FIND_VAR(v, id);
  for (lv=variables; lv!=NULL; lv=lv->next) mark_var(lv->var, FALSE);
  mark_var(v, TRUE);
  return get_n_leaves(v, TRUE);
}

int count_id_inodes(char *id, double *ndes)
{
  var_type *v;
  list_of_vars *lv;
  int n;

  FIND_VAR(v, id);
  for (lv=variables; lv!=NULL; lv=lv->next) mark_var(lv->var, FALSE);
  mark_var(v, TRUE);
  ndes_inodes = 0;
  n = get_n_inodes(v, TRUE);
  printf("xxx %d\n", ndes_inodes);
  *ndes = ndes_inodes / (double) n;
  return n;
}
