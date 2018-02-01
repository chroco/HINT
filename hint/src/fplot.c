/*
c89 -c -O -D_HPUX_SOURCE +z -I../FORMS -I/usr/include/X11R5 demo15.c
c89 -O -s  demo15.o  -o demo15 -L../FORMS -lforms -lXpm -L/usr/lib/X11R5 -lX11 -lm 
*/

#include <stdio.h>
#include <stdlib.h>
#include "forms.h"
#define GL
#include "sds.h"
#include "colormap.h"

#define MAX_OUT 50
#define MAX_IN 50
#define MIN_MESH 2
#define MAX_MESH 50

#define POINT_SIZE 8

typedef struct {
  FL_FORM *fform;

  FL_OBJECT *exit, *evaluate, *redraw, *read_data;

  FL_OBJECT *grp_update, *upd_man, *upd_auto;
  FL_OBJECT *grp_plot, *plot_man, *plot_auto;
  FL_OBJECT *grp_lplot, *lplot_man, *lplot_auto;

  FL_OBJECT *xgrid, *ygrid;
  FL_OBJECT *colormap;
  FL_OBJECT *plot;

  FL_OBJECT *xmin, *xmax, *ymin, *ymax, *xslide, *yslide;
  FL_OBJECT *grp_x, *grp_y;
  FL_OBJECT *x[MAX_IN], *y[MAX_IN], *val[MAX_IN],
  *slide[MAX_IN], *qval[MAX_IN];
  FL_OBJECT *grp_out, *out[MAX_OUT], *oval[MAX_OUT], *oqval[MAX_OUT];

  void *vdata;
  long ldata;
} FD_fform;

extern FD_fform * create_form_fform(void);
typedef struct {
	FL_FORM *plot;
	FL_OBJECT *lplot;
	void *vdata;
	long ldata;
} FD_pform;

FD_fform *fform;
FD_pform *pform;

var_type *out_var[MAX_OUT];	/* variables with fuzzy mem functions */
list_of_vars *lout;		/* same, but as list of variables */
var_type *in_var[MAX_IN];	/* variables with fuzzy mem functions */
int nout, nin;			/* number of them */
int zaxis, xaxis, yaxis;	/* indx's of vars to plot */
char auto_upd;			/* update of out values automatically? */
char auto_plot;			/* update visualisation automatically? */
char auto_lplot;		/* update of line plot automatically? */
double xmin, xmax, ymin, ymax;	/* plot limits */
int xgrid = 5, ygrid = 5;	/* grid size */
char line_plot = FALSE;		/* line plot is shown? */

void update_out_vars();
double *mycolor;

/****************************************************************************
Geomview related stuff - visualisation
****************************************************************************/

struct {
  double x, y, z;
} xyz[100][100];

				/* scaling */
#define SX(x) (((x)-xmin)/(xmax-xmin))
#define SY(y) (((y)-ymin)/(ymax-ymin))
#define SZ(z) (((z)-zmin)/(zmax-zmin))

void geomview_func()
{
  int i, j, n;
  char b;
  list_of_vars plotl;
  double zmin, zmax;
  double xval, yval, zval;
  double *d;

  plotl.next = NULL;
  plotl.var = out_var[zaxis];

  xval = in_var[xaxis]->val; yval = in_var[yaxis]->val;
  derive_lvar(&plotl);
  zval = out_var[zaxis]->val;

  zmax = -1e100; zmin = 1e100;
  for (j=0; j<ygrid; j++) {
    in_var[yaxis]->val = ymin + (ymax-ymin)*j/(ygrid-1.);
    real_to_fuzzy(in_var[yaxis]);
    
    for (i=0; i<xgrid; i++) {
      xyz[i][j].x = in_var[xaxis]->val = xmin + (xmax-xmin)*i/(xgrid-1.);
      xyz[i][j].y = in_var[yaxis]->val;
      real_to_fuzzy(in_var[xaxis]);
      derive_lvar(&plotl);
      if (zmax<plotl.var->val) zmax=plotl.var->val;
      if (zmin>plotl.var->val) zmin=plotl.var->val;
      xyz[i][j].z = plotl.var->val;
    }
  }

  for (i=0; i<xgrid; i++) 
    for (j=0; j<ygrid; j++) {
      xyz[i][j].x = SX(xyz[i][j].x);
      xyz[i][j].y = SY(xyz[i][j].y);
      xyz[i][j].z = SZ(xyz[i][j].z);
    }

  printf("(read geometry { define myfunct \n");
  printf("LIST {\n");

  printf("{appearance {+face +edge linewidth 1}\n");
  printf("{\n");
  printf("{\ndefine mymesh\nCMESH\n%d %d\n", xgrid, ygrid);
  for (j=0; j<ygrid; j++) {
    for (i=0; i<xgrid; i++) {
/*      printf("%10.5e %10.5e %10.5e\n", xyz[i][j].x, xyz[i][j].y, xyz[i][j].z); */
      printf("%5.3lf %5.3lf %5.3lf\n", xyz[i][j].x, xyz[i][j].y, xyz[i][j].z);
/*      n = (xyz[i][j].z-zmin)/(zmax-zmin) * N_COPPER; */
      n = xyz[i][j].z * N_COPPER;
      d = mycolor+n*3;
/*      fprintf(stderr, "xxx %d %d %d", i, j, n);
      fprintf(stderr, "(%5.3lf %5.3lf %5.3lf)\n", *(d), *(d+1), *(d+2));  */
      printf("%5.3lf %5.3lf %5.3lf 1\n",
	     mycolor[n*3], mycolor[n*3+1], mycolor[n*3+2]);
/*      printf("%5.3lf %5.3lf %5.3lf 1\n", *(d), *(d+1), *(d+2));   */
/*   printf("%d %d %d 1\n", c_copper[n][0], c_copper[n][1], c_copper[n][2]); */
    }
  }
  printf("}\n");		/* ZMESH ends here */
  printf("}}\n");		/* its appearance definitions ends here */

  printf("define myballs\n{appearance\n");
  printf("{ +edge linewidth %d}\n", POINT_SIZE);
  printf("{\n");
  printf("define coord\n");
  printf("VECT\n4 4 4\n1 1 1 1\n1 1 1 1\n");
  printf("%10.5e %10.5e %10.5e \n%10.5e %10.5e %10.5e\n%10.5e %10.5e %10.5e\n%10.5e %10.5e %10.5e\n",
	 SX(xmax), SY(ymin), SZ(zmin),
	 SX(xmin), SY(ymax), SZ(zmin),
	 SX(xmin), SY(ymin), SZ(zmax),
	 SX(xmin), SY(ymin), SZ(zmin)
	 );
  printf("0 1 0 1\n1 0 0 1\n0 1 1 1\n0 0 0 1");
  printf("}}\n");		/* axis definition end here */

  printf("define myball\n");
  printf("{appearance\n");
  printf("{ +edge linewidth %d}\n", POINT_SIZE);
  printf("{\n");
  printf("VECT\n1 1 1\n1\n1\n");

  
  printf("%10.5e %10.5e %10.5e\n", SX(xval), SY(yval), SZ(zval));
  printf("1 1 1 1\n");
  printf("}}\n");

  printf("}\n");		/* LIST ends here */
  printf("})\n");		/* (read geometry ends here */

				/* defining a background */
  printf("(read geometry { define myback \n");
  printf("{appearance {-face +edge linewidth 1}\n");
  printf("{\n");
  printf("LIST {\n");
  printf("{\nMESH\n%d %d\n", xgrid, ygrid);
  for (j=0; j < ygrid; j++)
    for (i=0; i < xgrid; i++) {
      printf("%8.3e %8.3e %8.3e\n", 
	     SX(xmin + (xmax-xmin)*i/(xgrid-1.)),
	     SY(ymin + (ymax-ymin)*j/(ygrid-1.)),
	     SZ(zmin));

/*      printf("%8.3e %8.3e %8.3e\n", xyz[i][j].x, xyz[i][j].y, SZ(zmin)); */
  }
  printf("}\n");		/* ZMESH ends here */

  printf("{\nuMESH\n%d %d\n", xgrid, out_var[zaxis]->ndesc);
  for (j=0; j < out_var[zaxis]->ndesc; j++) {
  for (i=0; i < xgrid; i++)
    printf("%8.3e %8.3e %8.3e\n", SX(xmin + (xmax-xmin)*i/(xgrid-1.)), 
	   SY(ymin), 
	   SZ(zmin + (zmax-zmin)*j/(out_var[zaxis]->ndesc-1.)));

    printf("\n");
  }
  printf("}\n");

  printf("{\nvMESH\n%d %d\n", ygrid, out_var[zaxis]->ndesc);
  for (j=0; j < out_var[zaxis]->ndesc; j++) {
    for (i=0; i < ygrid; i++) 
      printf("%3.8e %8.3e %8.3e\n", SX(xmin), 
	     SY(ymin + (ymax-ymin)*i/(ygrid-1.)),
	     SZ(zmin + (zmax-zmin)*j/(out_var[zaxis]->ndesc-1.))); 
    printf("\n");
  }
  printf("}\n");

  printf("}}\n");		/* its appearance definitions ends here */
  printf("}\n");		/* LIST ends here */
  printf("})\n");		/* (read geometry ends here */

  fflush(stdout);
  in_var[xaxis]->val = xval; in_var[yaxis]->val = yval;
  real_to_fuzzy(in_var[xaxis]); real_to_fuzzy(in_var[yaxis]); 
  out_var[zaxis]->val = zval;
}

/****************************************************************************
Line plot related stuff - plotting
****************************************************************************/

float xdata[MAX_MESH], ydata[MAX_MESH];
float px[2], py[2];
char did_lplot=FALSE;

void draw_line_plot()
{
  int i;
  list_of_vars plotl;
  Str255 str;
  double xval, zval;

  plotl.next = NULL;
  plotl.var = out_var[zaxis];

  xval = in_var[xaxis]->val;
  derive_lvar(&plotl);
  zval = out_var[zaxis]->val;

  for (i=0; i<xgrid; i++) {
    xdata[i] = in_var[xaxis]->val = xmin + (xmax-xmin)*i/(xgrid-1.);
    real_to_fuzzy(in_var[xaxis]);
    derive_lvar(&plotl);
    ydata[i] = plotl.var->val;
  }

  fl_freeze_form(pform->plot);
  fl_set_xyplot_xbounds(pform->lplot, xdata[0], xdata[xgrid-1]);
  fl_set_xyplot_ybounds(pform->lplot,
			out_var[zaxis]->vmin, out_var[zaxis]->vmax);
  sprintf(str, "%s(%s)", out_var[zaxis]->name, in_var[xaxis]->name);
  fl_set_xyplot_data(pform->lplot, xdata, ydata, xgrid, str,
		     in_var[xaxis]->name, out_var[zaxis]->name);

  in_var[xaxis]->val = xval;
  out_var[zaxis]->val = zval;

/*  px[0] = in_var[xaxis]->val;
  px[1] = px[0]+0.01;
  py[0] = out_var[zaxis]->vmin;
  py[1] = out_var[zaxis]->vmax;
  if (did_lplot) fl_delete_xyplot_overlay(pform->lplot, 1);
  if (!did_lplot) fl_add_xyplot_overlay(pform->lplot, 1, px, py, 2, FL_RED);
  did_lplot = TRUE; */
  fl_unfreeze_form(pform->plot);
}

/****************************************************************************
Utility functions, called from callback functions
****************************************************************************/

void update_out_vars() 
{
  static Str255 s;
  int i;

  derive_lvar(lout);
  for (i=0; i<nout; i++) {
    sprintf(s, "%10.5e", out_var[i]->val);
/*    printf("%s\n", s); */
    fl_set_object_label(fform->oval[i], s);
    fl_set_object_label(fform->oqval[i], get_q_val(out_var[i]));
  }
}

void set_plot_limits(FL_OBJECT *cmin, FL_OBJECT *cmax, 
		       FL_OBJECT *s, int indx, double min, double max)
{
  double v, w, x;
  fl_set_counter_value(cmax, max);
  fl_set_counter_value(cmin, min);
  w = (max-min)/(in_var[indx]->vmax -in_var[indx]->vmin);
  fl_set_slider_size(s,w);
/*  v = (min+max)/2; */

  x = (min+max)/2;
  w = max - min;
  if ((in_var[indx]->vmax -in_var[indx]->vmin - w) == 0)
    v = (in_var[indx]->vmax + in_var[indx]->vmin)/2.0;
  else
    v = (x-in_var[indx]->vmin-w/2.) /
      (in_var[indx]->vmax -in_var[indx]->vmin - w) *
	(in_var[indx]->vmax -in_var[indx]->vmin) + in_var[indx]->vmin;

/*  printf("%lf %lf - %lf %lf - %lf\n", in_var[indx]->vmin, in_var[indx]->vmax, min, max, v); */
  fl_set_slider_value(s, v);
}

set_ini_plot_limits(FL_OBJECT *cmin, FL_OBJECT *cmax, FL_OBJECT *s, int indx)
{
  fl_set_counter_bounds(cmin, in_var[indx]->vmin, in_var[indx]->vmax);
  fl_set_counter_bounds(cmax, in_var[indx]->vmin, in_var[indx]->vmax);
  fl_set_counter_step(cmin,
		      (in_var[indx]->vmax -in_var[indx]->vmin)/500.0,
		      (in_var[indx]->vmax -in_var[indx]->vmin)/50.0);
  fl_set_counter_step(cmax,
		      (in_var[indx]->vmax -in_var[indx]->vmin)/500.0,
		      (in_var[indx]->vmax -in_var[indx]->vmin)/50.0);

  fl_set_slider_bounds(s, in_var[indx]->vmin, in_var[indx]->vmax);
  fl_set_slider_step(s, (in_var[indx]->vmax -in_var[indx]->vmin)/100.0);
}

void set_ini_x_plot_limits()
{
  set_ini_plot_limits(fform->xmin, fform->xmax, fform->xslide, xaxis);
  set_plot_limits(fform->xmin, fform->xmax, fform->xslide, xaxis, xmin, xmax);
}

void set_ini_y_plot_limits()
{
  set_ini_plot_limits(fform->ymin, fform->ymax, fform->yslide, yaxis);
  set_plot_limits(fform->ymin, fform->ymax, fform->yslide, yaxis, ymin, ymax);
}

/****************************************************************************
Callback functions
****************************************************************************/

void do_plot(FL_OBJECT *obj, long data)
{
  if (line_plot) {
    fl_set_object_label(obj, "Show Plot");
    fl_hide_form(pform->plot);
  }
  else {
    fl_set_object_label(obj, "Hide Plot");
    draw_line_plot();
    fl_show_form(pform->plot,FL_PLACE_CENTER,FL_FULLBORDER,"Line plot");
  }
  line_plot = !line_plot;
}

void do_colormap(FL_OBJECT *obj, long data)
{
  int i;
  i = fl_get_choice(obj);
/*  fprintf(stderr, "yyy %d\n", i); */
  switch (i) {
  case 1:
    mycolor = (double *) c_hot;
    break;
  case 2:
    mycolor = (double *) c_cool;
    break;
  case 3:
    mycolor = (double *) c_copper;
    break;
  case 4:
    mycolor = (double *) c_hsv;
    break;
  case 5:
    mycolor = (double *) c_gray;
    break;
  }
  if (auto_plot) geomview_func();
}

void do_grid_size(FL_OBJECT *obj, long i)
{
  if (i==0)
    xgrid = fl_get_counter_value(obj);
  else
    ygrid = fl_get_counter_value(obj);
  if (auto_plot) geomview_func();
  if (auto_lplot) draw_line_plot();
}
  

void do_redraw(FL_OBJECT *obj, long i)
{
  geomview_func();
}

void do_evaluate(FL_OBJECT *obj, long i)
{
  update_out_vars();
}

void do_i_slider(FL_OBJECT *obj, long i)
{
  double zmax, zmin; int ii,j;
  in_var[i]->val = fl_get_slider_value(obj);
  real_to_fuzzy(in_var[i]);

  fl_set_counter_value(fform->val[i], in_var[i]->val);
  fl_set_object_label(fform->qval[i], get_q_val(in_var[i]));
  if (auto_upd) update_out_vars();
  if (auto_plot) geomview_func();
  if (auto_lplot) draw_line_plot();
}

void do_i_counter(FL_OBJECT *obj, long i)
{
  in_var[i]->val = fl_get_counter_value(obj);
  real_to_fuzzy(in_var[i]);
  fl_set_slider_value(fform->slide[i], in_var[i]->val);
  fl_set_object_label(fform->qval[i], get_q_val(in_var[i]));
  if (auto_upd) update_out_vars();
  if (auto_plot) geomview_func();
  if (auto_lplot) draw_line_plot();
}

void do_x(FL_OBJECT *obj, long data)
{
  if (data == yaxis) {
    fl_set_button(obj, 0);
    fl_set_button(fform->x[xaxis], 1);
  }
  else {
    xaxis = data;
    xmin = in_var[xaxis]->vmin;
    xmax = in_var[xaxis]->vmax;
    set_ini_x_plot_limits();
    if (auto_plot) geomview_func();
    if (auto_lplot) draw_line_plot();
  }
}

void do_out(FL_OBJECT *obj, long data)
{
  zaxis = data;
  if (auto_plot) geomview_func();
  if (auto_lplot) draw_line_plot();
}

void do_y(FL_OBJECT *obj, long data)
{
  if (data == xaxis) {
    fl_set_button(obj, 0);
    fl_set_button(fform->y[yaxis], 1);
  }
  else {
    yaxis = data;
    ymin = in_var[yaxis]->vmin;
    ymax = in_var[yaxis]->vmax;
    set_ini_y_plot_limits();
    if (auto_plot) geomview_func();
    if (auto_lplot) draw_line_plot();
  }
}

void do_set_update(FL_OBJECT *obj, long data)
{
  auto_upd = data;
}

void do_set_plot(FL_OBJECT *obj, long data)
{
  auto_plot = data;
}

void do_set_lplot(FL_OBJECT *obj, long data)
{
  auto_lplot = data;
}

void do_xmin(FL_OBJECT *obj, long data)
{
  xmin = fl_get_counter_value(obj);
  if (xmin < xmax)
    set_plot_limits(fform->xmin,fform->xmax,fform->xslide, xaxis, xmin, xmax);
  else {
    xmin -= (in_var[xaxis]->vmax -in_var[xaxis]->vmin)/100.0;
    fl_set_counter_value(fform->xmin, xmin);
  }
  if (auto_plot) geomview_func();
  if (auto_lplot) draw_line_plot();
}

void do_xmax(FL_OBJECT *obj, long data)
{
  xmax = fl_get_counter_value(obj);
  if (xmin < xmax)
    set_plot_limits(fform->xmin,fform->xmax,fform->xslide, xaxis, xmin, xmax);
  else {
    xmax += (in_var[xaxis]->vmax -in_var[xaxis]->vmin)/100.0;
    fl_set_counter_value(fform->xmax, xmax);
  }
  if (auto_plot) geomview_func();
  if (auto_lplot) draw_line_plot();
}

void do_ymin(FL_OBJECT *obj, long data)
{
  ymin = fl_get_counter_value(obj);
  if (ymin < ymax)
    set_plot_limits(fform->ymin,fform->ymax,fform->yslide, yaxis, ymin, ymax);
  else {
    ymin -= (in_var[yaxis]->vmax -in_var[yaxis]->vmin)/100.0;
    fl_set_counter_value(fform->ymin, ymin);
  }
  if (auto_plot) geomview_func();
}

void do_ymax(FL_OBJECT *obj, long data)
{
  ymax = fl_get_counter_value(obj);
  if (ymin < ymax)
    set_plot_limits(fform->ymin,fform->ymax,fform->yslide, yaxis, ymin, ymax);
  else {
    ymax += (in_var[yaxis]->vmax -in_var[yaxis]->vmin)/100.0;
    fl_set_counter_value(fform->ymax, ymax);
  }
  if (auto_plot) geomview_func();
}

/****************************************************************************
Interface Setup
****************************************************************************/

FD_fform *create_form_fform(void)
{
  int off;
  int top, otop;
  FL_OBJECT *obj;
  int i;
  FD_fform *fdui = (FD_fform *) fl_calloc(1, sizeof(FD_fform));

  top = 155 + nin * 30 + nout * 35 + 30;

  fdui->fform = fl_bgn_form(FL_NO_BOX, 690, top+5);
  obj = fl_add_box(FL_UP_BOX,0,0,690,top+5,"");

				/* main and simple buttons */
  fdui->redraw = obj = fl_add_button(FL_NORMAL_BUTTON,595,45,75,25,"Redraw");
  fl_set_object_callback(obj, do_redraw, NULL);
  fdui->exit = obj = fl_add_button(FL_NORMAL_BUTTON,595,15,75,25,"Exit");
  fdui->evaluate=obj = fl_add_button(FL_NORMAL_BUTTON,420,45,75,25,"Evaluate");
  fl_set_object_callback(obj, do_evaluate, NULL);
  fdui->read_data=obj=fl_add_button(FL_NORMAL_BUTTON,510,45,75,25,"Read Data");
  fdui->plot = obj = fl_add_button(FL_NORMAL_BUTTON,510,15,75,25,"Show Plot");
  fl_set_object_callback(obj,do_plot,0);

				/* colormap */
  fdui->colormap = obj = fl_add_choice(FL_DROPLIST_CHOICE,420,15,75,25,"");
  fl_set_object_boxtype(obj,FL_DOWN_BOX);
  fl_addto_choice(obj, "Hot");
  fl_addto_choice(obj, "Cool");
  fl_addto_choice(obj, "Copper");
  fl_addto_choice(obj, "Hsv");
  fl_addto_choice(obj, "Gray");
  fl_set_object_callback(obj,do_colormap,0);

				/* mesh size */
  fdui->xgrid = obj = fl_add_counter(FL_SIMPLE_COUNTER,85,50,100,25,"");
  fl_set_counter_bounds(obj,MIN_MESH,MAX_MESH);
  fl_set_counter_step(obj, 1, 1);
  fl_set_counter_value(obj, xgrid);
  fl_set_counter_return(obj, TRUE);
  fl_set_counter_precision(obj, 0);
  fl_set_object_callback(obj, do_grid_size, 0);
  fdui->ygrid = obj = fl_add_counter(FL_SIMPLE_COUNTER,85,20,100,25,"");
  fl_set_counter_bounds(obj,MIN_MESH,MAX_MESH);
  fl_set_counter_step(obj, 1, 1);
  fl_set_counter_return(obj, TRUE);
  fl_set_counter_precision(obj, 0);
  fl_set_counter_value(obj, ygrid);
  fl_set_object_callback(obj, do_grid_size, 1);
  obj = fl_add_text(FL_NORMAL_TEXT,10,50,70,25,"Mesh Size");

				/* update type */
    fl_set_object_callback(obj,do_set_update,0);

  obj = fl_add_text(FL_NORMAL_TEXT,210,50,45,25,"Update");
  fdui->grp_update = fl_bgn_group();
  fdui->upd_man=obj= fl_add_checkbutton(FL_RADIO_BUTTON,260,50,75,25,"Manual");
  fl_set_object_callback(obj, do_set_update, 0);
  fdui->upd_auto=obj= fl_add_checkbutton(FL_RADIO_BUTTON,330,50,70,25,"Auto");
  fl_set_object_callback(obj, do_set_update, 1);
  fl_set_button(obj, 1);
  fl_end_group();

				/* plot type */
  obj = fl_add_text(FL_NORMAL_TEXT,210,30,50,25,"Visualize");
  fdui->grp_plot= fl_bgn_group();
  fdui->plot_man=obj=fl_add_checkbutton(FL_RADIO_BUTTON,260,30,75,25,"Manual");
  fl_set_object_callback(obj, do_set_plot, 0);
  fdui->plot_auto=obj=fl_add_checkbutton(FL_RADIO_BUTTON,330,30,70,25,"Auto");
  fl_set_object_callback(obj, do_set_plot, 1);
  fl_set_button(obj, 1);
  fl_end_group();

				/* lplot type */
  obj = fl_add_text(FL_NORMAL_TEXT,210,10,45,25,"Plot");
  fdui->grp_lplot = fl_bgn_group();
  fdui->lplot_man=obj=fl_add_checkbutton(FL_RADIO_BUTTON,260,10,75,25,"Manual");
  fl_set_object_callback(obj,do_set_plot,0);
  fdui->lplot_auto=obj=fl_add_checkbutton(FL_RADIO_BUTTON,330,10,70,25,"Auto");
  fl_set_object_callback(obj,do_set_plot,1);
  fl_set_button(obj, 1);
  fl_end_group();


  obj = fl_add_text(FL_NORMAL_TEXT,5,75,680,15,"@DnLine");

				/* plot limits */
  obj = fl_add_text(FL_NORMAL_TEXT,10,125,90,25,"X Plot Limits");
  fdui->xmin = obj = fl_add_counter(FL_NORMAL_COUNTER,100,125,175,25,"");
  fl_set_counter_precision(obj, 3);
  fl_set_object_callback(obj, do_xmin, 0);
  fl_set_counter_return(obj, TRUE);
  fdui->xmax = obj = fl_add_counter(FL_NORMAL_COUNTER,285,125,195,25,"");
  fl_set_counter_precision(obj, 3);
  fl_set_object_callback(obj, do_xmax, 0);
  fl_set_counter_return(obj, TRUE);
  fdui->xslide = obj = fl_add_slider(FL_HOR_SLIDER,490,125,180,25,"");
  fl_deactivate_object(obj);

  obj = fl_add_text(FL_NORMAL_TEXT,10,95,90,25,"Y Plot Limits");
  fdui->ymin = obj = fl_add_counter(FL_NORMAL_COUNTER,100,95,175,25,"");
  fl_set_counter_precision(obj, 3);
  fl_set_object_callback(obj, do_ymin, 0);
  fl_set_counter_return(obj, TRUE);
  fdui->ymax = obj = fl_add_counter(FL_NORMAL_COUNTER,285,95,195,25,"");
  fl_set_counter_precision(obj, 3);
  fl_set_object_callback(obj, do_ymax, 0);
  fl_set_counter_return(obj, TRUE);
  fdui->yslide = obj = fl_add_slider(FL_HOR_SLIDER,490,95,180,25,"");
  fl_deactivate_object(obj);

  obj = fl_add_text(FL_NORMAL_TEXT,5,155,675,15,"@DnLine");

				/* output variables */
  otop = top;
  fdui->grp_out = fl_bgn_group();
  for (i=0; i<nout; i++) {
    top -= 35;
    fdui->out[i] =obj = fl_add_checkbutton(FL_RADIO_BUTTON,15,top,125,30,
				     out_var[i]->name);
    fl_set_object_callback(obj, do_out, (long) i);
  }
  fl_end_group();
  fl_set_button(fdui->out[zaxis], 1);

  for (top=otop, i=0; i<nout; i++) {
    top -= 35;
    fdui->oval[i] = obj = fl_add_text(FL_NORMAL_TEXT,140,top,145,30,"oval");
    fl_set_object_boxtype(obj,FL_DOWN_BOX);
    fdui->oqval[i] = obj = fl_add_text(FL_NORMAL_TEXT,290,top,380,30,
					"oqval0");
    fl_set_object_boxtype(obj,FL_DOWN_BOX);
  }

  top -= 15;
  obj = fl_add_text(FL_NORMAL_TEXT,5,top,675,15,"@DnLine");
  
				/* input variables */
  otop = top;
  fdui->grp_x = fl_bgn_group();
  for (i=0; i<nin; i++) {
    top -= 30;
    fdui->x[i] = obj = fl_add_checkbutton(FL_RADIO_BUTTON,15,top,25,25,"");
    fl_set_object_color(obj,FL_LEFT_BCOL,FL_GREEN);
    fl_set_object_callback(obj, do_x, (long) i);
  }
  fl_end_group();
  fl_set_button(fdui->x[xaxis], 1);

  fdui->grp_y = fl_bgn_group();
  for (top=otop, i=0; i<nin; i++) {
    top -= 30;
    fdui->y[i] = obj = fl_add_checkbutton(FL_RADIO_BUTTON,45,top,80,25,
				    in_var[i]->name);
    fl_set_object_color(obj,FL_LEFT_BCOL,FL_RED);
    fl_set_object_callback(obj, do_y, (long) i);
  }
  fl_end_group();
  fl_set_button(fdui->y[yaxis], 1);

  for (top=otop, i=0; i<nin; i++) {
    top -= 30;
    fdui->val[i] = obj = fl_add_counter(FL_NORMAL_COUNTER,135,top,195,25,"");
    fl_set_counter_return(obj, TRUE);
    fl_set_counter_bounds(obj, in_var[i]->vmin, in_var[i]->vmax);
    fl_set_counter_value(obj, in_var[i]->val);
    fl_set_counter_step(obj,
			(in_var[i]->vmax -in_var[i]->vmin)/500.0,
			(in_var[i]->vmax -in_var[i]->vmin)/50.0);
    fl_set_counter_precision(obj, 3);
    fl_set_object_callback(obj, do_i_counter, i);
    fdui->slide[i] = obj = fl_add_slider(FL_HOR_SLIDER,340,top,90,25,"");
    fl_set_slider_bounds(obj, in_var[i]->vmin, in_var[i]->vmax);
    fl_set_slider_value(obj, in_var[i]->val);
    fl_set_slider_step(obj,(in_var[i]->vmax -in_var[i]->vmin)/100.0);
    fl_set_object_callback(obj, do_i_slider, i);

    fdui->qval[i] = obj = fl_add_text(FL_NORMAL_TEXT,445,top,225,25,"qval0");
    fl_set_object_boxtype(obj,FL_DOWN_BOX);
  }

  fl_end_form();

  return fdui;
}

FD_pform *create_form_pform(void)
{
  FL_OBJECT *obj;
  FD_pform *fdui = (FD_pform *) fl_calloc(1, sizeof(FD_pform));

  fdui->plot = fl_bgn_form(FL_NO_BOX, 380, 210);
  obj = fl_add_box(FL_UP_BOX,0,0,380,210,"");
  fdui->lplot = obj = fl_add_xyplot(FL_NORMAL_XYPLOT,10,10,360,190,"");
  fl_end_form();

  return fdui;
}
void create_the_forms()
{
  fform = create_form_fform();
  pform = create_form_pform();
  set_ini_x_plot_limits();
  set_ini_y_plot_limits();
}

/****************************************************************************
Initialization & Main rutine
****************************************************************************/

/* process_sds: sds related part of prules. reads the input command
   file and builds the property structure */

void process_sds(int argc, char *argv[], char *fname)
{
  init_init(argv, argc);
  yyin = fopen(fname, "r");
  if (yyin == NULL) {
    fprintf(stderr, "%s: could not open cmd file %s\n", argv[0], argv[1]);
    exit(0);
  }

  interact_mode = FALSE;
  yyparse();
  fclose(yyin);
}

/* sets the vmin and vmax for the property according to given
   descriptions, assumes at least two descriptions of type left and
   right. BBB some more caution should be paid here */

void set_vmin_max(var_type *v)
{
  v->vmin = v->desc[0].a;
  v->vmax = v->desc[v->ndesc-1].c;
}

/* init_fplot: finds the output and interm properties and sets the
   default property to observe */

void init_fplot()
{
  int i;
  list_of_vars *l, *tlv;

  for(lout=NULL, nout=0, l = variables; l!=NULL; l=l->next)
    if ((l->var->ndesc > 0 && l->var->desc[0].tp!=none)
	&& (l->var->famsout != NULL)) {
      out_var[nout++] = l->var;
      tlv = (list_of_vars*)malloc(sizeof(*tlv));
      tlv->var = l->var;
      tlv->next = lout;
      lout = tlv;
      set_vmin_max(l->var);
    }

  for(nin=0, l = variables; l!=NULL; l=l->next)
    if ((l->var->ndesc > 0 && l->var->desc[0].tp!=none)
	&& (l->var->famsout == NULL)) {
      in_var[nin] = l->var;
      set_vmin_max(l->var);
      in_var[nin]->val = (in_var[nin]->vmin + in_var[nin]->vmax)/2;
      real_to_fuzzy(in_var[nin]);
      l->var->valdef = TRUE;
      nin++;
    }

/*  printf("in:\n");
  for (i=0; i<nin; i++) printf("%s\n", in_var[i]->name);
  printf("out:\n");
  for (i=0; i<nout; i++) printf("%s\n", out_var[i]->name); */
/*  for (tlv=out; tlv!=NULL; tlv=tlv->next) printf("%s\n", tlv->var->name); */

  xaxis = 0; yaxis = 1; zaxis = 0;
  xmin = in_var[xaxis]->vmin;
  xmax = in_var[xaxis]->vmax;
  ymin = in_var[yaxis]->vmin;
  ymax = in_var[yaxis]->vmax;
  auto_upd = TRUE;
  auto_plot = TRUE;
  auto_lplot = TRUE;
  mycolor = (double *) c_cool;
}

/****************************************************************************
Main rutine
****************************************************************************/

main(int argc, char *argv[])
{
  char *fname;
  FL_OBJECT *obj;
  int i;

  fl_initialize(argv[0], "FormDemo", 0, 0 ,&argc, argv);

/*  strcpy(fname, "bbb.s"); */
  fname = (char *) 
    fl_show_fselector("Select property strucuture file", 0, "*.s*", 0);
  if (fname == NULL) exit(0);
  process_sds(argc, argv, fname);
  init_fplot();
  create_the_forms();
  for (i=0; i<nin; i++)
    fl_set_object_label(fform->qval[i], get_q_val(in_var[i]));
  update_out_vars();

  printf("(geometry funct { : myfunct})\n"); 
  printf("(geometry background { : myback})\n");
/*  printf("(geometry balls { : myballs})\n");
  printf("(geometry ball { : myball})\n");*/

  geomview_func();
  printf("(bbox-draw World no)\n");
  printf("(camera-reset Camera)");
/*  printf("(transform funct funct focus rotate -1.5708 0 0)\n");
  printf("(transform funct funct focus rotate 0 -1.5708 0)\n");
  printf("(transform background background focus rotate -1.5708 0 0)\n");
  printf("(transform background background focus rotate 0 -1.5708 0)\n"); */

  fl_show_form(fform->fform,FL_PLACE_CENTER,FL_FULLBORDER,"Plot Functions");
  do obj = fl_do_forms();  while (obj != fform->exit);
  fl_hide_form(fform->fform);
}
