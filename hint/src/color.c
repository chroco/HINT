/****************************************************************************
color.c

Graph coloring rutines
****************************************************************************/

#include <stdio.h>
#define GL extern
#include "sds.h"

int num_nodes;			/* number of nodes to color */
int maxc;			/* max number of colors to use */
char *vcr=NULL, *colors=NULL;	/* tmp and final colors of nodes */
extern int nvcol;
extern char nondis_dec;
double **comp_evidence;		/* evidence for compatibility */

/****************************************************************************
Backtracking coloring algorithm for a fixed number of colors
see JA McHugh: Algorithmic Graph Theory, p 231, Prentice-Hall, 1990
****************************************************************************/

void next(char **im, int k)
{
  char b;
  int i;

  vcr[k]++;
  do {
    for (b=TRUE, i=0; i<num_nodes && b; i++)
      if (im[k][i])
	b = vcr[i]!=vcr[k];
    if (!b) vcr[k]++;
/*    printf("xx %d %d (%d %d)\n", vcr[k], b, vcr[im[k][i]], vcr[k]); */
  } while (!b);
}

int color_backtrack(char **im, int M)
{
  int i, j, k = 0;

  if (deb_dec>2) printf("Coloring for %d ->", M);
  for (i=0; i<num_nodes; i++) vcr[i] = -1;
  do {
    next(im, k);
    if (vcr[k]<M && k<num_nodes-1) 
      k++;			/* advance */
    else if (vcr[k]<M && k==num_nodes-1) { /* found */
      for (j=0, i=0; i<num_nodes; i++) if (j<vcr[i]) j=vcr[i];
      if (deb_dec>2) printf(" yes\n");
      return j+1;
    }
    else if (vcr[k]>=M && k>0)	/* backup */
      vcr[k--]=0;
    else if (vcr[k]>=M && k==0) {
      if (deb_dec>2) printf(" no\n");
      return FALSE;
    }

/*    for (i=0; i<num_nodes; i++) printf("%d ", vcr[i]);
    printf("   %d\n", k); */
  } while (TRUE);
}

int color_backtrack_set(char **im, int M)
{
  int i, j, k = 0;
  int n=0;

  if (deb_dec>2) printf("Coloring for %d\n", M);
  for (i=0; i<num_nodes; i++) vcr[i] = -1;
  do {
    next(im, k);
    if (vcr[k]<M && k<num_nodes-1) 
      k++;			/* advance */
    else if (vcr[k]<M && k==num_nodes-1) { /* found */
      for (j=0, i=0; i<num_nodes; i++) if (j<vcr[i]) j=vcr[i];

      printf("%d: ", n++);
      for (i=0; i<=M; i++) printf("%d ", vcr[i]); printf("\n");

/*      return j+1; */
    }
    else if (vcr[k]>=M && k>0)	/* backup */
      vcr[k--]=0;
    else if (vcr[k]>=M && k==0)
      return FALSE;

/*    for (i=0; i<num_nodes; i++) printf("%d ", vcr[i]);
    printf("   %d\n", k); */
  } while (TRUE);
}

/****************************************************************************
Bisection used to find a minimum number of colors to color a graph.
Start with (1,maxc).
****************************************************************************/

#define COPY_VCR  {int ii; for (ii=0; ii<num_nodes; ii++) colors[ii]=vcr[ii];}

int find_minc_bisect(char **im, int hi)
{
  int i, lo=1, med;

  if (color_backtrack(im, lo)) {
    COPY_VCR;
    return lo;
  }
  if (!(hi = color_backtrack(im, hi))) {
    printf("warning: could not color");
    return FALSE;
  }
  COPY_VCR;
  
  while (lo + 1 < hi) {
    med = (lo+hi)/2;
    if (i = color_backtrack(im, med)) {
      hi = i;
      COPY_VCR;
    }
    else 
      lo = med;
  }
  if (deb_dec>2) {
    printf("Coloring returns #colors  %d\n", hi+1);
    pspace(nvcol); for (i=0; i<num_nodes; i++) printf("%d", colors[i]);
    printf("\n");
  }
  return hi;
}

/****************************************************************************
Perkowski's coloring
****************************************************************************/


#define SWITCH(i,j,k) {k=i; i=j; j=k;}

int *nedges;			/* number of edges for each node */
int *nsort;			/* indices of nodes sorted by num of edges */
char **freec;			/* is color j free for node i */

int perkowski(char **im)
{
  int i, j, k, indx, best;
  char b;
  double best_ev, ev;
  extern char use_evidence;
  void new_color_analysis(char **im, char *colors);

  for (i=0; i<num_nodes; i++) colors[i]=-1; /* BZ */
  nedges = i_vector_ini(num_nodes, 0);
  nsort = i_vector(num_nodes);
  freec = c_matrix_ini(num_nodes, maxc, TRUE);

  for (i=0; i<num_nodes; i++) nsort[i] = i;
  for (i=0; i<num_nodes; i++)
    for (j=0; j<num_nodes; j++) {
      if (im[i][j]) nedges[i]++;
      im[i][j] = im[i][j];
    }

  for (j=0; j<num_nodes-1; j++)
    for (i=0; i<num_nodes-1; i++)
      if (nedges[i] < nedges[i+1]) {
	SWITCH(nedges[i], nedges[i+1], k);
	SWITCH(nsort[i], nsort[i+1], k);
      }

  for (i=0; i<num_nodes; i++) {
    indx = nsort[i];


    if (use_evidence) {
      best_ev=-1;
      for (j=0; j<maxc; j++)
	if (freec[indx][j]) {
	  ev = 0;
	  for (b=FALSE, k=0; k<num_nodes; k++)
	    if (colors[k]==j) {b=TRUE; ev += comp_evidence[indx][k];}
	  if (b && ev>best_ev) { best_ev = ev; best = j;}
	}

      if (best_ev>0) {
	colors[indx] = best;
	j = best;
      }
      else {
	for (j=0; !freec[indx][j]; j++);
	colors[indx] = j;
      }


    }
    else {
      for (j=0; !freec[indx][j]; j++);
      colors[indx] = j;
    }

    for (k=0; k<num_nodes; k++)
      if (im[k][indx])
	freec[k][j] = FALSE;
  }

  for (j=0, i=0; i<num_nodes; i++) if (colors[i]>j) j=colors[i];

  if (deb_dec>2) {
    printf("Coloring returns #colors %d (%d)\n", j+1, nvcol);
    pspace(nvcol); for (i=0; i<num_nodes; i++) printf("%d", colors[i]);
    printf("\n");
  }

/* special for AIM J */
/*  printf("\nColor analysis:\n"); */
/*  new_color_analysis(im, colors); */
/*   printf("End Analysis\n\n"); */

  free(nedges); free(nsort);
  return j+1;
}

/****************************************************************************
Set-up and main rutine
****************************************************************************/

void print_im(char *s, char **im, int num_nodes)
{
  int i, j;

  printf("%s:\n", s);
  for (i=0; i<num_nodes; i++) {
    pspace(nvcol);
    for (j=0; j<num_nodes; j++)
      printf("%c", im[i][j] ? '1' : '0');
    printf("\n");
  }
}

void ini_test_color(char **im)
{
  FILE *f;
  int i, j, ndata;

  f = fopen("coltmp", "r");
  if (f==NULL) {printf("Can't open coltmp\n"); exit(0);}
  
  fscanf(f, "%d", &num_nodes);
  im = c_matrix_ini(num_nodes, num_nodes, FALSE);
  while (fscanf(f, "%d %d", &i, &j) != EOF) {
    im[i][j] = im[j][i] = TRUE;
  }
  if (deb_dec>2) print_im("Incompat matrix", im, num_nodes);
}

int count_nnodes(char **im)
{
  int i, j, k, max=0;

  for (i=0; i<num_nodes; i++) {
    k=0;
    for (j=0; j<num_nodes; j++) if (im[i][j]) k++;
    if (k>max) max=k;
  }
  return max;
}

#define TESTING FALSE

int save_no = 0;

void save_im_to_file(char **im, int ncol)
{
  Str255 s;
  int i, j;

  FILE *f;
  sprintf(s, "%s%d", im_fname, save_no);
  f = fopen(s, "w");
  if (f==NULL) {printf("error: can't open im file %s\n", s); return;}
  fprintf(f, "%d\n", num_nodes);
  for (i=0; i<num_nodes; i++)
    for (j=0; j<num_nodes; j++)
      if (im[i][j]) fprintf(f, "%d %d\n", i, j);
  fclose(f);
  f = fopen("log_im", "a");
  fprintf(f, "%s -> %d (%s)\n", s, ncol,
	  coloring==col_optimal?"opt":"sort");
  fclose(f);
  save_no++;
}

/****************************************************************************
MAIN RUTINE TO BE CALLED FOR COLORING
****************************************************************************/


int color_graph(char **im, int nn)
{
  int i;

  num_nodes = nn;
  if (TESTING) {
    ini_test_color(im); 
    deb_dec = 4;
  }

  if (deb_dec>2) print_im("Incompat matrix:",im, nn);
  maxc = count_nnodes(im) + 1;
  FREE(vcr); FREE(colors);
  vcr = c_vector(num_nodes);
  colors = c_vector(num_nodes);
  if (deb_dec>2) printf("Start with %d colors\n", maxc);

  if (coloring == col_optimal)
    i = find_minc_bisect(im, perkowski(im));
/*    i = find_minc_bisect(im, maxc); */
  else if (coloring == col_heuristic || coloring == col_soft)
    i = perkowski(im);
  else if (coloring == col_set) {
    i = find_minc_bisect(im, perkowski(im));
    printf("Start with %d colors\n", i);
    color_backtrack_set(im, i);    
  }

  if (save_im) save_im_to_file(im, i);
  free(vcr);
  return i;
}


/****************************************************************************
ANALYSIS OF COLORING SOLUTION AND ITERACTIVE COLORING
****************************************************************************/

extern Str255 joined_var_name;

char nth_digit(int k, int n)
{
  int i;
  for (i=0; i<n; i++) k = k/10;
  return k % 10;
}

extern char **dm, **im, **ndm, **nim, *dcol, *tdcol, *ndcol;
extern int nrow, ncol, ncn, nboth;
extern char nc, nnc;
extern char *colors;
extern fams *sf;
extern char *is_col, *is_row, *is_both;
extern Str255 joined_var_name;		/* global storage for new node */

char *sure_color = NULL;	/* the colors that are trivial to define */
char *needed_color = NULL;	/* the colors that would be needed */
int *nex = NULL;		/* number of examples for each colon */
char *ccol;			/* current color */

void list_color(char *s, char *color)
{
  int i, max=0;

  for (i=0; i<ncol; i++) if (color[i]!=CUNDEF) max = MMAX(max,color[i]);

  printf("%s (%d colors)\n", s, max+1);
  pspace(nvcol);
  for (i=0; i<ncol; i++) 
    if (color[i]!=CUNDEF) printf("%d", color[i]); else printf("-");
  printf("\n");
}

void list_ndm()
{
  int i, j;
  printf("ndm:\n");
  for (j=0; j<nboth; j++) {
    for (i=0; i<nnc; i++)
      if (ndm[j][i]==CUNDEF) printf("- "); else printf("%d ", ndm[j][i]);
    printf("\n");
  }
}

void list_current_color()
{
  list_color("Current color", ccol);
  list_color("Example map color", dcol);
}

void color_analyze_start()
{
  c_vector_copy(tdcol, &ccol, ncol);
}

/* derive sure colors */

char *nsure_color;

void new_derive_sure(char ncolors, char ncol, char **im, char *col)
{
  int i, j, k, indx, best;
  char b;
  double best_ev, ev;
  extern char use_evidence;
  char *poss;			/* possible color to assign */
  char *my_col;

  nsure_color = c_vector_ini(num_nodes, CUNDEF);

  nedges = i_vector_ini(num_nodes, 0);
  nsort = i_vector(num_nodes);
  freec = c_matrix_ini(num_nodes, maxc, TRUE);
  my_col = c_vector_ini(num_nodes, CUNDEF);

				/* count connectivity */
  for (i=0; i<num_nodes; i++) nsort[i] = i;
  for (i=0; i<num_nodes; i++)
    for (j=0; j<num_nodes; j++) {
      if (im[i][j]) nedges[i]++;
      im[i][j] = im[i][j];
    }

				/* sorting */
  for (j=0; j<num_nodes-1; j++)
    for (i=0; i<num_nodes-1; i++)
      if (nedges[i] < nedges[i+1]) {
	SWITCH(nedges[i], nedges[i+1], k);
	SWITCH(nsort[i], nsort[i+1], k);
      }

  poss = c_vector(ncolors);

  for (i=0; i<num_nodes; i++) {
    indx = nsort[i];

    for (j=0; !freec[indx][j]; j++);
    my_col[indx] = j;

    /* are there nodes that are not adjacent but of different color? */
    c_vector_set(poss, ncolors, FALSE);
    for (k=0; k<num_nodes; k++)
      if ((indx!=k) && (!im[k][indx]) && (my_col[k] != CUNDEF)) {
	poss[my_col[k]] = TRUE;
      }

    nsure_color[indx] = 0;
    for (k=0; k<ncolors; k++) if (poss[k]) nsure_color[indx]++;
    
    for (k=0; k<num_nodes; k++)
      if (im[k][indx])
	freec[k][j] = FALSE;
  }

  list_color("sure sure", nsure_color);

}

void derive_sure(char ncolors, char ncol, char **my_im, char *my_col,
		 char *sure_color)
{
  char *sure = c_vector(ncolors);
  int i, j;
  char b, stop;
				/* first step, to see if some has less */
				/* colors in my_im than needed */
  for (i=0; i<ncol; i++) {
    for (j=0; j<ncolors; j++) sure[j]=FALSE;
    for (j=0; j<ncol; j++) if (my_im[i][j])
      if (my_col[j]!=CUNDEF)
	sure[my_col[j]] = TRUE;
    for (b=TRUE, j=0; j<ncolors && b; j++) if (my_col[i]!=j) b=sure[j];
    sure_color[i] = b;
  }
  
  /* iterative steps, colors needed should be those that have sure colors */

  /* for each column, it should not be compatible to those columns that are
     sure and are colored with different color */

  stop = FALSE;
  while (!stop) {
    stop = TRUE;
    for (i=0; i<ncol; i++) 
      if (sure_color[i]) {
/*	for (j=0; j<ncolors; j++) sure[j]=FALSE;
	for (j=0; j<ncol; j++) 
	  if (my_im[i][j] && sure_color[j]) sure[my_col[j]] = TRUE; */

	for (j=0; j<ncolors; j++) sure[j]=TRUE;
	for (j=0; j<ncol; j++) 
	  if (!my_im[i][j] && sure_color[j]) sure[my_col[j]] = FALSE;

	for (b=TRUE, j=0; j<ncolors && b; j++) 
	  if (my_col[i]!=j) b=sure[j];

	if (!b) {
	  sure_color[i] = FALSE;
	  stop = FALSE;
	}
      }
  }
  FREE(sure);
}

void analyze_nd_count_entries()
{
  char *old_dcol = dcol;
  FREE(dcol); c_vector_copy(ccol, &dcol, ncol);
  nnc = nd_count_entries(nc);
  map_nd_colors(dcol, ndcol);
}

/* analyze_color: for disjoint or non-disjoint case computes the redundant and
   sure colors. */

void analyze_color(char *my_color)
{
  int i, j;
  char class, *att, b, *tmp;
  rule_type rtp;

  analyze_nd_count_entries();
  print_im("R Incompat matrix:", im, ncol);
  free(needed_color); free(nex);
  needed_color = c_vector_ini(ncol, FALSE);
  tmp = c_vector(ncol);
  nex = i_vector_ini(ncol, 0);

  list_color("Start color", ccol);
				/* count #examples in each colon */
  nex = i_vector_ini(ncol, 0);
  reset_next_rule(sf);
  while (get_next_rulea(sf, &class, &rtp, &att)) {
    i = satt2indx(sf, att, is_col);
    nex[i]++;
  }

				/* report the table attributes */
  printf("of %s/%d dec %s = ", sf->out->name, sf->n_in, joined_var_name);
  for (i=0; i<sf->n_in; i++) if (is_col[i]) printf("%s ", sf->in[i]->name);
  printf("\n");

  att = c_vector(sf->n_in);
  for (i=0; i<sf->n_in; i++) if (is_col[i]) {
    pspace(nvcol);
    for (j=0; j<ncol; j++) {
      indx2satt(sf, j, att, is_col);
      printf("%c", sf->in[i]->desc[att[i]].name[0]);
    }
    printf("\n");
  }
  FREE(att);

  reset_next_rule(sf); i=0;
  while (get_next_rulea(sf, &class, &rtp, &att)) i++;
				/* report the number of col elements */
  printf("Num examples (%d total)\n", i);
  for (j=5; j>=0; j--) {
    for (b=FALSE, i=0; i<ncol && !b; i++) b = nth_digit(nex[i],j) != 0;
    if (b) {
      pspace(nvcol);
      for (i=0; i<ncol; i++) printf("%d", nth_digit(nex[i],j)); printf("\n");
    }
  }
				/* derive needed colors */
  needed_color = c_vector_ini(ncol, FALSE);
  reset_next_rule(sf);
  while (get_next_rulea(sf, &class, &rtp, &att)) {
    i = satt2indx(sf, att, is_col);
    needed_color[i] = TRUE;
  }
/*  for (i=0; i<ncol; i++) if (!needed_color[i]) ccol[i] = CUNDEF; */

  for (i=0; i<ncol; i++) tmp[i] = needed_color[i] ? ccol[i] : CUNDEF;
  list_color("Needed color", tmp);

  FREE(sure_color);
  sure_color = c_vector_ini(ncol, FALSE);
  derive_sure(nc, ncol, im, ccol, sure_color);
  
  for (i=0; i<ncol; i++) tmp[i] = sure_color[i] ? ccol[i] : CUNDEF;
  list_color("Sure color", tmp);
}

void analyze_current_color()
{
  analyze_color(ccol);
}

void set_color(Str255 s)
{
  int i, max;
  for (i=0; i<ncol; i++) 
    if (s[i] != '-') ccol[i] = s[i]-'0'; else ccol[i]=CUNDEF;
  list_color("Color set to", ccol);
  for (i=0; i<ncol; i++) dcol[i] = ccol[i];
  for (i=0, max=0; i<ncol; i++) max = MMAX(max, dcol[i]);
  nc = nnc = max + 1;
}

void set_sure_color()
{
  int i;
  for (i=0; i<ncol; i++) if (!sure_color[i]) ccol[i] = CUNDEF;
/*  list_color("Color set to", ccol); */
}

void set_need_color()
{
  int i;
  void set_nim();

  for (i=0; i<ncol; i++) if (!needed_color[i]) {
    ccol[i] = CUNDEF; dcol[i] = CUNDEF;
  }
  list_color("Color set to", ccol);
  analyze_nd_count_entries();
/*  set_nim(); */
}

/* list_color_groups: for the color[i] lists, for every color, those
   substructure rules that belong to that group. Also, it lists a
   number of examples that belong to that rule. */

char *sure_groups=NULL;

void list_color_groups()
{
  int i, j, k, max=0;
  char *att;

  for (i=0; i<ncol; i++) if (ccol[i]!=CUNDEF) max = MMAX(max,ccol[i]);

  att = c_vector(sf->n_in);
  for (i=0; i<=max; i++) {
    printf("Group %d:\n", i);
    for (j=0; j<ncol; j++)
      if (ccol[j]==i && sure_color[j]) {
	indx2satt(sf, j, att, is_col);
	printf("  %3d:  ", j);
	for (k=0; k<sf->n_in; k++)
	  if (is_col[k]) printf("%s ", sf->in[k]->desc[att[k]].name);
	printf("\n");
      }
  }
  printf("\n");
  FREE(att);

  if (!nondis_dec) return;

  printf("Map to %d sets, i.e. ", nnc);
  for (i=0; i<nnc; i++) {
    printf("%d (", i);
    for (j=0; j<nc; j++) if (ndcol[j]==i) printf("%d ", j);
    printf(") ");
  }
  printf("\n");

  FREE(sure_groups);
  sure_groups = c_vector(nc);
  derive_sure(nnc, nc, nim, ndcol, sure_groups);
  printf("Sure groups: ");
  for (i=0; i<nc; i++) 
    if (!sure_groups[i] || ndcol[i]==CUNDEF) printf("-");
    else printf("%d", i);
  printf("\n");
}

void list_unsure_instances()
{
  int i, j, k;
  char *att, *compat;

  att = c_vector(sf->n_in);
  compat = c_vector(nc);

  for (i=0; i<ncol; i++)
    if ((ccol[i]==CUNDEF || !sure_color[i]) && needed_color[i]) {
      indx2satt(sf, i, att, is_col);
      printf("  %3d:  ", i);
      for (k=0; k<sf->n_in; k++)
	if (is_col[k]) printf("%s ", sf->in[k]->desc[att[k]].name);
      printf(" --> ");
      
      /* find all colors that this instance can got to */
      for (j=0; j<nc; j++) compat[j]=TRUE;
      for (j=0; j<ncol; j++) 
	if (j!=i && sure_color[j] && im[j][i]) compat[ccol[j]]=FALSE;
      for (j=0; j<nc; j++) if (compat[j]) printf("%d ", j);
      printf("\n");
    }
  FREE(att);
}

void set_diff_color(int a, int b)
{
}




/************************************************************************
new color analysis as from July 1997 (for AIM J)

*************************************************************************/

/* new_color_analysis: for AIM J. */

void new_color_analysis(char **im, char *colors)
{
  void print_im(char *s, char **im, int num_nodes);
  int i, j;
  char class, *att, b, *tmp;
  rule_type rtp;
  char trace = TRUE;

  for (nc=0, i=0; i<ncol; i++) if (colors[i]>nc) nc=colors[i];
  nc++;

  if (trace) print_im("IM", im, ncol);

  free(needed_color); free(nex);
  needed_color = c_vector_ini(ncol, FALSE);
  tmp = c_vector(ncol);
  ccol = c_vector(ncol);
  nex = i_vector_ini(ncol, 0);

  for (i=0; i<ncol; i++) ccol[i] = colors[i];

  if (trace) list_color("Start color", ccol);

				/* count #examples in each colon */
  nex = i_vector_ini(ncol, 0);
  reset_next_rule(sf);
  while (get_next_rulea(sf, &class, &rtp, &att)) {
    i = satt2indx(sf, att, is_col);
    nex[i]++;
  }

				/* report the table attributes */
  if (trace) {
    printf("of %s/%d dec %s = ", sf->out->name, sf->n_in, joined_var_name);
    for (i=0; i<sf->n_in; i++) if (is_col[i]) printf("%s ", sf->in[i]->name);
    printf("\n");
  

    att = c_vector(sf->n_in);
    for (i=0; i<sf->n_in; i++) if (is_col[i]) {
      pspace(nvcol);
      for (j=0; j<ncol; j++) {
	indx2satt(sf, j, att, is_col);
	printf("%c", sf->in[i]->desc[att[i]].name[0]);
      }
      printf("\n");
    }
    FREE(att);
  }

  reset_next_rule(sf); i=0;
  while (get_next_rulea(sf, &class, &rtp, &att)) i++;
				/* report the number of col elements */
  if (trace) {
  printf("Num examples (%d total)\n", i);
  for (j=5; j>=0; j--) {
    for (b=FALSE, i=0; i<ncol && !b; i++) b = nth_digit(nex[i],j) != 0;
    if (b) {
      pspace(nvcol);
      for (i=0; i<ncol; i++) printf("%d", nth_digit(nex[i],j)); printf("\n");
    }
  }
}
				/* derive needed colors */
  needed_color = c_vector_ini(ncol, FALSE);
  reset_next_rule(sf);
  while (get_next_rulea(sf, &class, &rtp, &att)) {
    i = satt2indx(sf, att, is_col);
    needed_color[i] = TRUE;
  }
/*  for (i=0; i<ncol; i++) if (!needed_color[i]) ccol[i] = CUNDEF; */

  for (i=0; i<ncol; i++) tmp[i] = needed_color[i] ? ccol[i] : CUNDEF;
  list_color("Needed color", tmp);

  FREE(sure_color);
  sure_color = c_vector_ini(ncol, FALSE);

  new_derive_sure(nc, ncol, im, ccol);
  for (i=0; i<ncol; i++) tmp[i] = nsure_color[i]>1 ? CUNDEF : ccol[i];
  if (trace) list_color("Sure color", tmp);

				/* report on sure */
  printf("TYPICAL FOR ");
  for (i=0; i<sf->n_in; i++) if (is_col[i]) printf("%s ", sf->in[i]->name);
  printf("\n");

  att = c_vector(sf->n_in);
  for (j=0; j<ncol; j++)
    if (nsure_color[j]<=1) {
    for (i=0; i<sf->n_in; i++) if (is_col[i]) {
      indx2satt(sf, j, att, is_col);
      printf("%s ", sf->in[i]->desc[att[i]].name);
    }
    printf("\n");
  }
  FREE(att);
  printf("END TYPICAL\n");

}

/****************************************************************************/

/* set_instance_to_color: sets the instance to specific color. If
   sure_color[ins], this can only be done for non-disjoint
   decomposition if the old and the new colors are not in conflict in
   nim. For disjoint decomposition, this should not be in conflict
   with any of sure colors. */

char set_im_instance_to_color(int ins, int new_col)
{
  int j, old_col = ccol[ins];

  if (old_col == CUNDEF) {
    for (j=0; j<ncol; j++)
      if (j!=ins && sure_color[j]) {
	if (ccol[j]==old_col) im[ins][j] = im[j][ins] = FALSE;
	else im[ins][j] = im[j][ins] = TRUE;
      }
  }
  else if (nondis_dec) {
    for (j=0; j<ncol; j++) {
      if (j!=ins) {
	if (sure_color[j]) {
	  if (ccol[j]==new_col) { im[ins][j] = im[j][ins] = FALSE; }
	  else { im[ins][j] = im[j][ins] = TRUE; }
	  }
	else {
/*BUG*/	
	  if (ccol[j]!=CUNDEF && nim[new_col][ccol[j]])
	    im[ins][j]=im[j][ins]=TRUE;
	  else im[ins][j]=im[j][ins]=FALSE; 
	}
      }
    }
  }
  ccol[ins] = new_col;
  return TRUE;
}

void set_nim()
{
  int i;

/*  print_im("RN Incompat matrix:", im, ncol); */
  nc = color_graph(im, ncol);
  for (i=0; i<ncol; i++) if (ccol[i]!=CUNDEF) ccol[i] = colors[i]; 
  analyze_color(ccol);

  if (nondis_dec) {
    analyze_nd_count_entries();
    list_color_groups();
  }
}

void set_instance_to_color(int ins, int new_col)
{
  printf("set instance %d %d\n", ins, new_col);

  if (set_im_instance_to_color(ins, new_col)) {
    if (nondis_dec) set_nim();
    else analyze_color(ccol);
/*    set_sure_color(); */
  }
}

/* set_same_color: for all in im, all that are compatible with a,
   should be compatible with b, and vice versa. Basically removes
   color b. */

void set_same_color(int a, int b)
{
  int i;

  if (!nim[a][b] || TRUE) {
    for (i=0; i<ncol; i++)
      if (ccol[i]==b) set_im_instance_to_color(i, a);
    set_nim();
  }
  else
    printf("Groups %d and %d are not compatible.\n", a, b);
}

/****************************************************************************
TRY TO USE BETTER COLORING this measures the "evidence" of the
incompatibility. For the coloring, if there are several candidates for
colors, it uses the one that is most evident.
****************************************************************************/

extern int nn_rules;
extern rule_list **rules;

derive_inc_evidence()
{
  int s1, e1, s2, e2;		/* limits of two subgroups within the
				   same group */
  int i, i1, i2, j, k, l;
  char c1, c2;

  extern double **d_matrix_ini(int nr, int nc, double ini);
  comp_evidence = d_matrix_ini(ncol, ncol, 0.);
/*  ppr(); */


  e1 = 0;
  for (i=e1; i<nn_rules; i++) {
    s1 = i;
    for (; rules[i]->unused; i++);
    e1 = i;
    printf("%d %d\n", s1, e1);
    
    if (rules[i]->mark) {
      j=e1+1;
      for (;; j++) {
	s2 = j;
	for (; rules[j]->unused; j++);
	e2 = j;

	printf("%d %d -- %d %d, cr=%d %d\n", s1, e1, s2, e2, c1, c2);


	if (c1>=ncol || c2>=ncol || c1<0 || c2<0) {
	  printf("ssx %d (%d-%d) (%d %d)\n", nn_rules, j, k, c1, c2); 
	  exit(0);
	}

	for (k=s1; k<=e1; k++)
	  for (l=s2; l<=e2; l++) {
	    c1 = rules[k]->cr;
	    c2 = rules[l]->cr;

	    comp_evidence[c1][c2] += 
	      rules[k]->dist[rules[k]->class] * rules[l]->dist[rules[l]->class];
	    comp_evidence[c2][c1] = comp_evidence[c1][c2];
	  }

	if (!rules[j]->mark) break;
      }
      i = e2;
    }
  }

  for (c1=0; c1<ncol; c1++) {
    for (c2=0; c2<ncol; c2++)
      printf("%5.2lf ", comp_evidence[c1][c2]);
    printf("\n");
  }
}

derive_evidence()
{
  int s1, e1, s2, e2;		/* limits of two subgroups within the
				   same group */
  int i, i1, i2, j, k, l;
  int c1, c2;

  extern double **d_matrix_ini(int nr, int nc, double ini);
  comp_evidence = d_matrix_ini(ncol, ncol, 0.);
/*  ppr(); */


  e1 = 0;
  for (i=e1; i<nn_rules; i++) {
    s1 = i;
    for (; rules[i]->unused; i++);
    e1 = i;
/*    printf("%d %d\n", s1, e1); */

    for (k=s1; k<=e1-1; k++)
      for (l=k+1; l<=e1; l++) {
	c1 = rules[k]->cr;
	c2 = rules[l]->cr;

	comp_evidence[c1][c2] += 
	  rules[k]->dist[rules[k]->class] * rules[l]->dist[rules[l]->class];
	comp_evidence[c2][c1] = comp_evidence[c1][c2];
      }
  }    

/*  for (c1=0; c1<ncol; c1++) {
    for (c2=0; c2<ncol; c2++)
      printf("%5.2lf ", comp_evidence[c1][c2]);
    printf("\n");
  } */
}
