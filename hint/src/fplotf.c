/* Form definition file generated with fdesign. */

#include "/usr/local/include/forms.h"
#include "fplotf.h"

FD_fform *create_form_fform(void)
{
  FL_OBJECT *obj;
  FD_fform *fdui = (FD_fform *) fl_calloc(1, sizeof(FD_fform));

  fdui->fform = fl_bgn_form(FL_NO_BOX, 690, 445);
  obj = fl_add_box(FL_UP_BOX,0,0,690,445,"");
  fdui->exit = obj = fl_add_button(FL_NORMAL_BUTTON,595,15,75,25,"Exit");
  fdui->y1 = obj = fl_add_checkbutton(FL_RADIO_BUTTON,45,170,80,25,"var1");
  fdui->x1 = obj = fl_add_checkbutton(FL_PUSH_BUTTON,15,170,25,25,"");
  obj = fl_add_counter(FL_NORMAL_COUNTER,135,170,195,25,"val1");
  fdui->slide1 = obj = fl_add_slider(FL_HOR_SLIDER,340,170,90,25,"");
  fdui->xmin = obj = fl_add_counter(FL_NORMAL_COUNTER,100,125,175,25,"");
  fdui->xmax = obj = fl_add_counter(FL_NORMAL_COUNTER,285,125,195,25,"");
  fdui->ymin = obj = fl_add_counter(FL_NORMAL_COUNTER,100,95,175,25,"");
  fdui->ymax = obj = fl_add_counter(FL_NORMAL_COUNTER,285,95,195,25,"");
  obj = fl_add_text(FL_NORMAL_TEXT,10,125,90,25,"X Plot Limits");
  obj = fl_add_text(FL_NORMAL_TEXT,10,95,90,25,"Y Plot Limits");
  obj = fl_add_text(FL_NORMAL_TEXT,5,75,680,15,"@DnLine");
  fdui->redraw = obj = fl_add_button(FL_NORMAL_BUTTON,595,45,75,25,"Redraw");
  fdui->y0 = obj = fl_add_checkbutton(FL_RADIO_BUTTON,45,200,80,25,"var0");
  fdui->x0 = obj = fl_add_checkbutton(FL_PUSH_BUTTON,15,200,25,25,"");
  fdui->slide0 = obj = fl_add_counter(FL_NORMAL_COUNTER,135,200,195,25,"");
  fdui->slide0 = obj = fl_add_slider(FL_HOR_SLIDER,340,200,90,25,"");
  obj = fl_add_text(FL_NORMAL_TEXT,5,155,675,15,"@DnLine");
  obj = fl_add_text(FL_NORMAL_TEXT,5,230,675,15,"@DnLine");
  fdui->xslide = obj = fl_add_slider(FL_HOR_SLIDER,490,125,180,25,"");
  fdui->qval0 = obj = fl_add_text(FL_NORMAL_TEXT,445,200,225,25,"qval0");
    fl_set_object_boxtype(obj,FL_DOWN_BOX);
  fdui->qval1 = obj = fl_add_text(FL_NORMAL_TEXT,445,170,225,25,"qval0");
    fl_set_object_boxtype(obj,FL_DOWN_BOX);
  fdui->yslide = obj = fl_add_slider(FL_HOR_SLIDER,490,95,180,25,"");
  fdui->evaluate = obj = fl_add_button(FL_NORMAL_BUTTON,420,45,75,25,"Evaluate");
  fdui->out0 = obj = fl_add_checkbutton(FL_PUSH_BUTTON,15,280,125,30,"Output0");
  fdui->oqval0 = obj = fl_add_text(FL_NORMAL_TEXT,290,280,380,30,"oqval0");
    fl_set_object_boxtype(obj,FL_DOWN_BOX);
  fdui->out1 = obj = fl_add_checkbutton(FL_PUSH_BUTTON,15,245,125,30,"Output1");
  fdui->oqval1 = obj = fl_add_text(FL_NORMAL_TEXT,290,245,380,30,"oqval1");
    fl_set_object_boxtype(obj,FL_DOWN_BOX);
  fdui->oval0 = obj = fl_add_text(FL_NORMAL_TEXT,140,280,145,30,"oval0");
    fl_set_object_boxtype(obj,FL_DOWN_BOX);
  fdui->oval0 = obj = fl_add_text(FL_NORMAL_TEXT,140,245,145,30,"oval0");
    fl_set_object_boxtype(obj,FL_DOWN_BOX);
  fdui->read_data = obj = fl_add_button(FL_NORMAL_BUTTON,510,45,75,25,"Read Data");
  fdui->xgrid = obj = fl_add_counter(FL_SIMPLE_COUNTER,85,50,100,25,"");
  fdui->xgrid = obj = fl_add_counter(FL_SIMPLE_COUNTER,85,20,100,25,"");
  obj = fl_add_text(FL_NORMAL_TEXT,10,50,70,25,"Mesh Size");
  fdui->colormap = obj = fl_add_choice(FL_DROPLIST_CHOICE,420,15,75,25,"");
    fl_set_object_boxtype(obj,FL_DOWN_BOX);
    fl_set_object_callback(obj,do_colormap,0);
  fdui->plot = obj = fl_add_button(FL_NORMAL_BUTTON,510,15,75,25,"Show Plot");
    fl_set_object_callback(obj,do_plot,0);

  fdui->grp_lplot = fl_bgn_group();
  fdui->plot_man = obj = fl_add_checkbutton(FL_RADIO_BUTTON,260,10,75,25,"Manual");
    fl_set_object_callback(obj,do_set_lplot,0);
  fdui->plot_auto = obj = fl_add_checkbutton(FL_RADIO_BUTTON,330,10,70,25,"Auto");
    fl_set_object_callback(obj,do_set_lplot,1);
  obj = fl_add_text(FL_NORMAL_TEXT,210,10,45,25,"Plot");
  fl_end_group();


  fdui->grp_plot = fl_bgn_group();
  fdui->plot_man = obj = fl_add_checkbutton(FL_RADIO_BUTTON,260,30,75,25,"Manual");
    fl_set_object_callback(obj,do_set_plot,0);
  fdui->plot_auto = obj = fl_add_checkbutton(FL_RADIO_BUTTON,330,30,70,25,"Auto");
    fl_set_object_callback(obj,do_set_plot,1);
  obj = fl_add_text(FL_NORMAL_TEXT,210,30,50,25,"Visualize");
  fl_end_group();


  fdui->grp_update = fl_bgn_group();
  fdui->upd_man = obj = fl_add_checkbutton(FL_RADIO_BUTTON,260,50,75,25,"Manual");
    fl_set_object_callback(obj,do_set_update,0);
  fdui->upt_auto = obj = fl_add_checkbutton(FL_RADIO_BUTTON,330,50,70,25,"Auto");
    fl_set_object_callback(obj,do_set_update,1);
  obj = fl_add_text(FL_NORMAL_TEXT,210,50,45,25,"Update");
  fl_end_group();

  fl_end_form();

  return fdui;
}
/*---------------------------------------*/

FD_pform *create_form_pform(void)
{
  FL_OBJECT *obj;
  FD_pform *fdui = (FD_pform *) fl_calloc(1, sizeof(FD_pform));

  fdui->pform = fl_bgn_form(FL_NO_BOX, 380, 210);
  obj = fl_add_box(FL_UP_BOX,0,0,380,210,"");
  fdui->lplot = obj = fl_add_xyplot(FL_NORMAL_XYPLOT,10,10,360,190,"");
  fl_end_form();

  return fdui;
}
/*---------------------------------------*/

