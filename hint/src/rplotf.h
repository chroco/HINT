#ifndef FD_brform_h_
#define FD_brform_h_
/* Header file generated with fdesign. */

/**** Callback routines ****/



/**** Forms and Objects ****/

typedef struct {
	FL_FORM *brform;
	FL_OBJECT *br;
	void *vdata;
	long ldata;
} FD_brform;

extern FD_brform * create_form_brform(void);

#endif /* FD_brform_h_ */
