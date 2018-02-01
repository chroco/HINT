#include <stdio.h>
#include <stdlib.h>
#define GL
#include "sds.h"

void main(int argc, char *argv[])
{
  init_init(argc, argv);
  
  if (argc < 2) {
    printf(">> ");
    load_cmds(TRUE);		/* commands from keyboard */
  }

  else {
    yyin = NULL;
    read_options(argc, argv);
    if (yyin != NULL)		/* commands from the file */
      load_cmds(FALSE);
    if (!n_errors) {		/* commands from keyboard */
      yyin = (FILE *) fdopen(0,"rt");
      ungetc('\n',yyin);
      load_cmds(TRUE);
    }
  }
  fclose(lfile);
}
