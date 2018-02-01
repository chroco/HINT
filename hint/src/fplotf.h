#ifndef FD_fform_h_
#define FD_fform_h_
/* Header file generated with fdesign. */

/**** Callback routines ****/

extern void do_colormap(FL_OBJECT *, long);
extern void do_plot(FL_OBJECT *, long);
extern void do_set_lplot(FL_OBJECT *, long);
extern void do_set_plot(FL_OBJECT *, long);
extern void do_set_update(FL_OBJECT *, long);



/**** Forms and Objects ****/

typedef struct {
	FL_FORM *fform;
	FL_OBJECT *exit;
	FL_OBJECT *y1;
	FL_OBJECT *x1;
	FL_OBJECT *slide1;
	FL_OBJECT *xmin;
	FL_OBJECT *xmax;
	FL_OBJECT *ymin;
	FL_OBJECT *ymax;
	FL_OBJECT *redraw;
	FL_OBJECT *y0;
	FL_OBJECT *x0;
	FL_OBJECT *slide0;
	FL_OBJECT *slide0;
	FL_OBJECT *xslide;
	FL_OBJECT *qval0;
	FL_OBJECT *qval1;
	FL_OBJECT *yslide;
	FL_OBJECT *evaluate;
	FL_OBJECT *out0;
	FL_OBJECT *oqval0;
	FL_OBJECT *out1;
	FL_OBJECT *oqval1;
	FL_OBJECT *oval0;
	FL_OBJECT *oval0;
	FL_OBJECT *read_data;
	FL_OBJECT *xgrid;
	FL_OBJECT *xgrid;
	FL_OBJECT *colormap;
	FL_OBJECT *plot;
	FL_OBJECT *grp_lplot;
	FL_OBJECT *plot_man;
	FL_OBJECT *plot_auto;
	FL_OBJECT *grp_plot;
	FL_OBJECT *plot_man;
	FL_OBJECT *plot_auto;
	FL_OBJECT *grp_update;
	FL_OBJECT *upd_man;
	FL_OBJECT *upt_auto;
	void *vdata;
	long ldata;
} FD_fform;

extern FD_fform * create_form_fform(void);
typedef struct {
	FL_FORM *pform;
	FL_OBJECT *lplot;
	void *vdata;
	long ldata;
} FD_pform;

extern FD_pform * create_form_pform(void);

#endif /* FD_fform_h_ */
