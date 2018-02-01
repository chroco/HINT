/* Form definition file generated with fdesign. */

#include "forms.h"
#include "rplotf.h"

FD_brform *create_form_brform(void)
{
  FL_OBJECT *obj;
  FD_brform *fdui = (FD_brform *) fl_calloc(1, sizeof(FD_brform));

  fdui->brform = fl_bgn_form(FL_NO_BOX, 430, 320);
  obj = fl_add_box(FL_UP_BOX,0,0,430,320,"");
  fdui->br = obj = fl_add_browser(FL_HOLD_BROWSER,20,20,390,280,"");
    fl_set_object_lstyle(obj,FL_FIXED_STYLE);
  fl_end_form();

  return fdui;
}
/*---------------------------------------*/

