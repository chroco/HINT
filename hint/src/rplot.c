#include <stdio.h>
#include <stdlib.h>
#include <forms.h>
#include "rplotf.h"
#define GL
#include "sds.h"

#define MAX_LN 10		/* max number of child properties */
#define POINT_SIZE 8

typedef struct {
  FL_FORM *mainform;
  FL_OBJECT *doshowrules, *doexit, *slide1, *pos1, *slide2, *pos2;
  FL_OBJECT *ygroup, *y1, *y2, *xgroup, *x1, *x2, *sel_property;
  void *vdata;
  long ldata;
} FD_mainform;

FD_mainform *fd_main;		/* main form, selects the property and
				   child properties values */
FD_brform *fd_br;

var_type *sel_var;		/* variable that is to observe */
fams *fam;			/* selected fam */
var_type **main_v;		/* list of output and interm properties */
int n_main_v;
int n_leaves = 0;		/* number of property childreen */
int v_sel = 0;

int xaxis = 0;			/* x axis */
int yaxis = 1;			/* y axis */
int selected=0;			/* selected point from browser */

FL_OBJECT *x[MAX_LN], *y[MAX_LN];
FL_OBJECT *slide[MAX_LN];
FL_OBJECT *bleft[MAX_LN], *value[MAX_LN], *bright[MAX_LN];
FL_OBJECT *txt_x, *txt_y, *txt_a, *txt_b;
int sel[MAX_LN];		/* selected qval */

char showhide = FALSE;		/* FALSE = show, TRUE = hide */

/****************************************************************************
Rule browser
****************************************************************************/

void show_rules()
{
  int i, j, n;
  Str255 s, s1, s2;
  char b;

  fl_freeze_form(fd_br->brform);
  fl_clear_browser(fd_br->br);
  for (i=0; i<n_leaves; i++)
    if ((i!=xaxis) && (i!=yaxis)) {
      sprintf(s, "%s = %s", fam->in[i]->name, fam->in[i]->desc[sel[i]].name);
      fl_addto_browser(fd_br->br, s);
    }
  fl_addto_browser(fd_br->br, "");
  sprintf(s, "%-10s %-10s   %-10s", fam->in[xaxis]->name,
	  fam->in[yaxis]->name, sel_var->name);
  fl_addto_browser(fd_br->br, s);
  fl_addto_browser(fd_br->br, " ");

  for (i=0; i < fam->n_rules; i++) {
    b = TRUE;
    for (n=1, j=0; j < n_leaves; n*=fam->in[j]->ndesc, j++)
      if (j!=xaxis && j!=yaxis)
	b = b && sel[j] == ((i / n) % fam->in[j]->ndesc);
    if (b) {
      sprintf(s, "");
      for (n=1, j=0; j < n_leaves; n*=fam->in[j]->ndesc, j++)
	if (j==xaxis)
	  sprintf(s, "%s%-10s ", s,
		  fam->in[j]->desc[(i / n) % fam->in[j]->ndesc].name);
      for (n=1, j=0; j < n_leaves; n*=fam->in[j]->ndesc, j++)
	if (j==yaxis)
	  sprintf(s, "%s%-10s ", s,
		  fam->in[j]->desc[(i / n) % fam->in[j]->ndesc].name);
      sprintf(s, "%s  %-10s", s, sel_var->desc[fam->rule[i]].name);
      fl_addto_browser(fd_br->br, s);
    }
  }
  fl_select_browser_line(fd_br->br, selected+(n_leaves-2)+3);
  fl_unfreeze_form(fd_br->brform);
}

/****************************************************************************
Plot rules
****************************************************************************/

void geomview_rules()
{
  int i, j, n;
  char b;

  printf("(read geometry { define myrules \n");
  printf("LIST {\n");

  printf("{appearance {+face +edge linewidth 3}\n");
  printf("{\n");
  printf("{\ndefine mymesh\nZMESH\n%d %d\n", 
	 fam->in[xaxis]->ndesc, fam->in[yaxis]->ndesc);
  for (i=0; i < fam->n_rules; i++) {
    b = TRUE;
    for (n=1, j=0; j < n_leaves; n*=fam->in[j]->ndesc, j++)
      if (j!=xaxis && j!=yaxis)
	b = b && sel[j] == ((i / n) % fam->in[j]->ndesc);
    if (b)
      printf("%d\n", fam->rule[i]);
  }
  printf("}\n");		/* ZMESH ends here */
  printf("}}\n");		/* its appearance definitions ends here */

  printf("{appearance\n");
  printf("{ +edge linewidth %d}\n", POINT_SIZE);
  printf("{\n");
  printf("define coord\n");
  printf("VECT\n3 3 3\n1 1 1\n1 1 1\n");
  printf("%d 0 0\n0 %d 0\n0 0 %d\n", fam->in[xaxis]->ndesc-1,
	 fam->in[yaxis]->ndesc-1, sel_var->ndesc-1);
  printf("0 1 0 1\n1 0 0 1\n0 1 1 1\n");
  printf("}}\n");		/* axis definition end here */

  printf("define reference\n");
  printf("{appearance\n");
  printf("{ +edge linewidth %d}\n", POINT_SIZE);
  printf("{\n");
  printf("VECT\n1 1 1\n1\n1\n");
  printf("0 0 0\n");
  printf("1 1 1 1\n");
  printf("}}\n");

  printf("}\n");		/* LIST ends here */
  printf("})\n");		/* (read geometry ends here */

				/* defining a background */
  printf("(read geometry { define myback \n");
  printf("{appearance {-face +edge linewidth 1}\n");
  printf("{\n");
  printf("LIST {\n");
  printf("{\nZMESH\n%d %d\n", fam->in[xaxis]->ndesc, fam->in[yaxis]->ndesc);
  for (i=0; i < fam->in[xaxis]->ndesc; i++) {
    for (j=0; j < fam->in[yaxis]->ndesc; j++)
      printf("%d ", 0);
    printf("\n", 0);
  }
  printf("}\n");		/* ZMESH ends here */
  printf("{\nuMESH\n%d %d\n", fam->in[xaxis]->ndesc, sel_var->ndesc);
  for (j=0; j < sel_var->ndesc; j++) {
  for (i=0; i < fam->in[xaxis]->ndesc; i++)
      printf("%d %d %d ", i, 0, j);
    printf("\n", 0);
  }
  printf("}\n");		/* ZMESH ends here */
  printf("{\nvMESH\n%d %d\n", fam->in[yaxis]->ndesc, sel_var->ndesc);
  for (j=0; j < sel_var->ndesc; j++) {
  for (i=0; i < fam->in[yaxis]->ndesc; i++)
      printf("%d %d %d ", 0, i, j);
    printf("\n", 0);
  }
  printf("}\n");		/* ZMESH ends here */
  printf("}}\n");		/* its appearance definitions ends here */
  printf("}\n");		/* LIST ends here */
  printf("})\n");		/* (read geometry ends here */

  fflush(stdout);
}


/****************************************************************************
Call-back functions
****************************************************************************/

void adjust_br_val(int x, int y);
void set_mainform();

void do_y(FL_OBJECT *obj, long data)
{
  if (data == xaxis) {
    fl_set_button(obj, 0);
    fl_set_button(y[yaxis], 1);
  }
  else {
    yaxis = data;
    selected = sel[xaxis] + fam->in[xaxis]->ndesc * sel[yaxis];
    show_rules(); 
    geomview_rules();
  }
}

void do_x(FL_OBJECT *obj, long data)
{
  if (data == yaxis) {
    fl_set_button(obj, 0);
    fl_set_button(x[xaxis], 1);
  }
  else {
    xaxis = data;
    selected = sel[xaxis] + fam->in[xaxis]->ndesc * sel[yaxis];
    show_rules(); 
    geomview_rules();
  }
}

void do_slide(FL_OBJECT *obj, long i)
{
  int j;

  j = (int) fl_get_slider_value(obj);
  sel[i] = j;
  fl_set_object_label(value[i], fam->in[i]->desc[sel[i]].name);
  if ((i!=xaxis) && (i!=yaxis)) {
    show_rules(); 
    geomview_rules();
  }
  else {
    adjust_br_val(sel[xaxis], sel[yaxis]);
    selected = sel[xaxis] + fam->in[xaxis]->ndesc * sel[yaxis];
    fl_select_browser_line(fd_br->br, selected+(n_leaves-2)+3);
  }
}

void do_left(FL_OBJECT *obj, long i)
{
  if (sel[i] != 0) {
    sel[i]--;
    fl_set_object_label(value[i], fam->in[i]->desc[sel[i]].name);
    fl_set_slider_value(slide[i], sel[i]);
    if ((i!=xaxis) && (i!=yaxis)) {
      show_rules(); 
      geomview_rules();
    }
    else {
      adjust_br_val(sel[xaxis], sel[yaxis]);
      selected = sel[xaxis] + fam->in[xaxis]->ndesc * sel[yaxis];
      fl_select_browser_line(fd_br->br, selected+(n_leaves-2)+3);
    }
  }
}

void do_right(FL_OBJECT *obj, long i)
{
  if (sel[i] < fam->in[i]->ndesc-1) {
    sel[i]++;
    fl_set_object_label(value[i], fam->in[i]->desc[sel[i]].name);
    fl_set_slider_value(slide[i], sel[i]);
    if ((i!=xaxis) && (i!=yaxis)) {
      show_rules(); 
      geomview_rules();
    }
    else {
      adjust_br_val(sel[xaxis], sel[yaxis]);
      selected = sel[xaxis] + fam->in[xaxis]->ndesc * sel[yaxis];
      fl_select_browser_line(fd_br->br, selected+(n_leaves-2)+3);
    }
  }
}

void do_showhide(FL_OBJECT *obj, long i)
{
  if (!showhide) {
    fl_set_object_label(obj, "Hide");
    show_rules(); 
    fl_show_form(fd_br->brform,FL_PLACE_CENTER,FL_FULLBORDER,"Rule Browser");
  }
  else {
    fl_set_object_label(obj, "Show");
    fl_hide_form(fd_br->brform);
  }
  showhide = !showhide;
}

FD_mainform *create_form_mainform(void);

void do_sel_property()
{
  int i;
  i = fl_get_choice(fd_main->sel_property) - 1;
  if (main_v[i] != sel_var) {
    sel_var = main_v[i];
    fam = sel_var->famsout;
    n_leaves = fam->n_in;

    fl_hide_form(fd_main->mainform);
    fl_free_form(fd_main->mainform);
    
    fd_main = create_form_mainform();
    set_mainform();
    fl_show_form(fd_main->mainform,FL_PLACE_CENTER,FL_FULLBORDER,"Plot Rules");
    v_sel = i;
    fl_set_choice(fd_main->sel_property, v_sel+1);
    show_rules();
    geomview_rules();
  }
}

void adjust_br_val(int x, int y) 
{
  int n, nn, i;
  sel[xaxis] = x;
  sel[yaxis] = y;
  fl_set_slider_value(slide[xaxis], x);
  fl_set_object_label(value[xaxis], fam->in[xaxis]->desc[sel[xaxis]].name);
  fl_set_slider_value(slide[yaxis], y);
  fl_set_object_label(value[yaxis], fam->in[yaxis]->desc[sel[yaxis]].name);

/*      fprintf(stderr, "selected %d %d\n", x, y); */

  for (nn=0, i=0, n=1; i<n_leaves; n*=fam->in[i]->ndesc, i++)
    nn += n*sel[i];

  printf("(read geometry { \n");
  printf("define reference\n");
  printf("{appearance\n");
  printf("{ +edge linewidth %d}\n", POINT_SIZE);
  printf("{\n");
  printf("VECT\n1 1 1\n1\n1\n");
  printf("%d %d %d\n", x, y, fam->rule[nn]);
  printf("1 1 1 1\n");
  printf("}}\n");
  printf("})\n");
  fflush(stdout);
}

void do_br(FL_OBJECT *obj, long i)
{
  int x, y;
  selected =  fl_get_browser(obj) - (n_leaves-2) - 3;
/*    fprintf(stderr, "selected %d %d\n", fl_get_browser(obj), selected); */
  x = selected % fam->in[xaxis]->ndesc;
  y = selected / fam->in[xaxis]->ndesc;
  adjust_br_val(x, y);
}

/****************************************************************************
Setting-up the interface
****************************************************************************/

FD_mainform *create_form_mainform(void)
{
  FL_OBJECT *obj;
  FD_mainform *fdui = (FD_mainform *) fl_calloc(1, sizeof(FD_mainform));
  int off;

  off = n_leaves * 40;

  fdui->mainform = fl_bgn_form(FL_NO_BOX, 440, 120+off);
  obj = fl_add_box(FL_UP_BOX,0,0,440,120+off,"");
    fl_set_object_lsize(obj,11);
  fdui->doshowrules = obj = fl_add_button(FL_NORMAL_BUTTON,235,20,90,30,
					  "Show Rules");
  fl_set_object_lsize(obj,FL_NORMAL_SIZE);
  fl_set_object_lstyle(obj,FL_BOLD_STYLE);
  fl_set_object_callback(obj, (void *) do_showhide, NULL);  

  fdui->doexit = obj = fl_add_button(FL_NORMAL_BUTTON,330,20,90,30,"Exit");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_lstyle(obj,FL_BOLD_STYLE);
  fdui->sel_property = obj = fl_add_choice(FL_DROPLIST_CHOICE,75,20,150,30,
					   "Property");
    fl_set_object_boxtype(obj,FL_DOWN_BOX);
  fl_end_form();

  return fdui;
}

void create_the_forms()
{
  fd_main = create_form_mainform();
  fd_br = create_form_brform();
}

/* process_sds: sds related part of prules. reads the input command
   file and builds the property structure */

void process_sds(int argc, char *argv[], char *fname)
{
/*  if (argc<=1) {
    fprintf(stderr, "usage: %s fname\n", argv[0]);
    exit(0);
  }
*/
  init_init(argv, argc);
  yyin = fopen(fname, "r");
  if (yyin == NULL) {
    fprintf(stderr, "%s: could not open cmd file %s\n", argv[0], argv[1]);
    exit(0);
  }
/*  load_cmds(FALSE); */

  interact_mode = FALSE;
  yyparse();
  fclose(yyin);
}

/* init_prules: finds the output and interm properties and sets the
   default property to observe */

void init_prules()
{
  int i;
  list_of_vars *l;
  for(n_main_v=0, l = variables; l!=NULL; l=l->next)
    if (l->var->famsout != NULL) n_main_v++;
  main_v = (var_type **) calloc(n_main_v, sizeof(*main_v));
  for(i=0, l = variables; l!=NULL; l=l->next)
    if (l->var->famsout != NULL)
      main_v[i++] = l->var;
  sel_var = main_v[0];
  fam = sel_var->famsout;
  n_leaves = fam->n_in;
}

void set_mainform()
{
  int off;
  int top;
  FL_OBJECT *obj;
  int i;

  for (i=0; i<n_main_v; i++) {
    fl_addto_choice(fd_main->sel_property, main_v[i]->name);
  }

  fl_addto_form(fd_main->mainform);
  top = 65+(n_leaves-1)*40;
  off = n_leaves * 40;
  fl_set_form_size(fd_main->mainform, 440, 120+off);

  obj = fl_add_text(FL_NORMAL_TEXT,340,85+off,80,20,"prules V1.0");
  fl_set_object_boxtype(obj,FL_NO_BOX);
  fl_set_object_color(obj,FL_COL1,FL_COL1);
  fl_set_object_lsize(obj,FL_TINY_SIZE);
  fl_set_object_lalign(obj,FL_ALIGN_RIGHT);
  obj = fl_add_text(FL_NORMAL_TEXT,310,70+off,110,20,"(C) 1995 Blaz Zupan");
  fl_set_object_boxtype(obj,FL_NO_BOX);
  fl_set_object_color(obj,FL_COL1,FL_COL1);
  fl_set_object_lsize(obj,FL_TINY_SIZE);
  fl_set_object_lalign(obj,FL_ALIGN_RIGHT);
  obj = fl_add_text(FL_NORMAL_TEXT,10,50+off,30,30,"X");
  fl_set_object_boxtype(obj,FL_NO_BOX);
  fl_set_object_color(obj,FL_COL1,FL_COL1);
  obj = fl_add_text(FL_NORMAL_TEXT,35,50+off,30,30,"Y");
  fl_set_object_boxtype(obj,FL_NO_BOX);
  fl_set_object_color(obj,FL_COL1,FL_COL1);

  fl_bgn_group();
  for (i=0; i<n_leaves; i++) {
    y[i] = obj = fl_add_lightbutton(FL_RADIO_BUTTON,40,top-i*40,
				    35,25, fam->in[i]->name);
    fl_set_object_boxtype(obj,FL_NO_BOX);
    fl_set_object_boxtype(obj,FL_NO_BOX);
    fl_set_object_color(obj,FL_LEFT_BCOL,FL_RED);
    fl_set_object_callback(obj, (void *) do_y, (long) i);
  }
  fl_set_button(y[yaxis], TRUE);
  fl_end_group();

  fl_bgn_group();
  for (i=0; i<n_leaves; i++) {
    x[i] = obj = fl_add_lightbutton(FL_RADIO_BUTTON,10.0,top-i*40,
				    35,25,"");
    fl_set_object_boxtype(obj,FL_NO_BOX);
    fl_set_object_color(obj,FL_LEFT_BCOL,FL_GREEN);
    fl_set_object_lsize(obj,11);
    fl_set_object_callback(obj, (void *) do_x, (long) i);
  }
  fl_set_button(x[xaxis], TRUE);
  fl_end_group();

  for (i=0; i<n_leaves; i++) {
    slide[i] = obj = fl_add_slider(FL_HOR_SLIDER,325,top-i*40,95,25,"");
    fl_set_object_callback(obj, (void *) do_slide, (long) i);
    fl_set_object_lsize(obj,11);
 
    fl_set_slider_bounds(slide[i], 0, fam->in[i]->ndesc-1);
    fl_set_slider_value(slide[i], 0);
    fl_set_slider_step(slide[i], 1.0);

    bleft[i] = obj = fl_add_button(FL_TOUCH_BUTTON,125.0,top-i*40,25,25,"@<");
    fl_set_object_callback(obj, (void *) do_left, (long) i);
    sel[i] = 0;
    value[i] = obj = fl_add_box(FL_BORDER_BOX,160,top-i*40,120,25,"");
      bright[i] = obj = fl_add_button(FL_TOUCH_BUTTON,290,top-i*40,25,25,"@>");
    fl_set_object_callback(obj, (void *) do_right, (long) i);
    fl_set_object_label(value[i], fam->in[i]->desc[sel[i]].name);
  }
  fl_end_form();
}


/****************************************************************************
Main rutine
****************************************************************************/

main(int argc, char *argv[])
{
  FL_OBJECT *obj;
  char *fname;

  fl_initialize(argv[0], "Plot Rules", 0, 0 ,&argc, argv);

  fname = (char *) 
    fl_show_fselector("Select property strucuture file", 0, "*.s*", 0);
  if (fname == NULL) exit(0);
  process_sds(argc, argv, fname);
  init_prules();
  create_the_forms();
  set_mainform();

  printf("(geometry rules { : myrules})\n"); /* why is this here? */
  printf("(geometry background { : myback})\n"); /* why is this here? */

  geomview_rules();
  printf("(bbox-draw World no)\n");
  printf("(camera-reset Camera)");
  printf("(transform rules rules focus rotate -1.5708 0 0)\n");
  printf("(transform rules rules focus rotate 0 -1.5708 0)\n");
  printf("(transform background background focus rotate -1.5708 0 0)\n");
  printf("(transform background background focus rotate 0 -1.5708 0)\n");

  fl_show_form(fd_main->mainform,FL_PLACE_CENTER,FL_FULLBORDER,"Plot Rules");
/*  fl_show_form(fd_br->brform,FL_PLACE_CENTER,FL_FULLBORDER,"Rule Browser"); */
  fl_set_browser_fontstyle(fd_br->br, FL_FIXED_STYLE);
  fl_set_browser_fontsize(fd_br->br, FL_SMALL_SIZE);
  fl_set_object_callback(fd_br->br, (void *) do_br, NULL);
  show_rules();
  do {
    obj = fl_do_forms();
    if (obj == fd_main->sel_property)
      do_sel_property();
  }
  while (obj != fd_main->doexit);
  fl_hide_form(fd_main->mainform);
}
