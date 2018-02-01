/*****************************************************************************/
/* DECLARATIONS                                                              */
/*****************************************************************************/

%{
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#define GL extern
#include "sds.h"

extern yylineno;		/* line number in screen def. source */
char *att = NULL;		/* attrib of property, used to enter rules */
int att_num = 0;

list_of_vars *tmp_vars=NULL;

%}

%union {
  int intV;
  double doubV;
  char chV;
  Str255 strV;
  char *chpV;
  list_of_vars *varL;
  list_of_str *strL;
  list_of_num *numL;
  list_of_q *strQ;
  fa_type faV;
  ga_error_type gaeV;
  ga_policy_type gapV;
  eval_method_type evalV;
  coloring_type colV;
  er_type erV;
}

%token <intV>	ADD	AND	VAR	QUIT	DEL	LIST
%token <intV>	VARS	DEPENDS	ON	IN	IS	STRUCTURE
%token <intV>	OF	FAM	FOR	RULE	ALL	DESCRIPTIONS
%token <intV>	SET	CURRENT	SHOW	OPT	NODE
%token <intV>	DERIVE	SEL	LOAD	UNDEF	TEST	EVALUATE
%token <intV>	EXPECT	COMMENT	LEARN	STAT	GA	POPULATION
%token <intV>	POINTS	PLOT	MAN	WAIT	DEBUG	ITERATIONS
%token <intV>	LEX	PRINT	FREQ	ERROR	MAX	YECHO
%token <intV>	SAVE	IMP	WEIGHT	MIN	MULT	INFORMATIVITY
%token <intV>	GAABS	YES	NO	PREFER MUTATION
%token <intV>	POLICY	CUSTOM	OR	LSQR	RESET	CROSSOVER
%token <intV>	USAGE	NORM	WIDTH	TABLE	SPLIT	DECOMPOSE
%token <intV>	REDO	GAMMA	TS	USING	CL	SIGM
%token <intV>   KNN	COLOR	COMPUTE	ONE	COMPOSE	HEURISTIC
%token <intV>   COPY	TO	QUALITATIVE	FROM	QUANTITATIVE
%token <intV>   COMPARE STEP	REPEAT	YLOG	GR	GINI
%token <intV>   RELIEFF FUZZY	CRISP	PREFIX	JOIN	INTERVAL INSTANCE
%token <intV>	TREE	DM	PJOIN	FJOIN	SURE	NEED	REDUNDANT
%token <intV>	VECT	SAME	GROUP	ANALYZE	NAMES	WRITE	DIFFERENT
%token <intV>	SOFT	MNUMBER	DOT	COUNT	LOCAL	GLOBAL
%token <intV>	SDTIC	DFC    	DTIC	CM
%token <intV>	MERGE	PAIRS	LAPLACE	CLUSTER	NOISE	SEED
%token <intV>	DONTCARE	DONTKNOW	APRIORY	DISTRIBUTION
%token <intV>	CV	SAMPLE	SORT	MDL	LEAF	NOD
%token <intV>	DUPLICATE	EXPAND

%token <doubV>	NUM
%token <intV>	INUM
%token <strV>	STRING	FNAME	ID

%type  <strL>	str_list	str_list_rev
%type  <numL>	num_list	num_list_rev
%type  <varL>   var_list	var_list_rev
%type  <strQ>	q_list		q_list_el
%type  <faV>	fand
%type  <gaeV>   ga_error_entry;
%type  <chV>	yesno
%type  <gapV>   ga_policy_entry
%type  <evalV>  eval_type;
%type  <colV>	coloring_entry;
%type  <intV>	att_entry, att_id;
%type  <erV>	encode_rules_type;
%type  <strV>	id_num
%type  <doubV>  num

%start	commands

/*****************************************************************************/
/* GRAMMAR RULES AND ACTIONS                                                 */
/*****************************************************************************/

%%

commands:	commands command '\n'
{ if (interact_mode) printf(">> "); fflush(stdout); }
	|	command '\n'	
{ if (interact_mode) printf(">> "); fflush(stdout); }
	|	commands COMMENT	{}
	|	COMMENT			{}
	|	error '\n'
{ 
  if (interact_mode) {
    printf("syntax error\n>> ");
    fflush(stdout); /* MyRefresh(); */
  }
  else {
    yyerror();
    exit(1);
  }
}
;

/****************************************************************************
  COMMAND LINE HANDLING
****************************************************************************/

command: {}

/****************************************************************************
  manipulation with variable lists and dependecies
****************************************************************************/

	|	ADD VAR ID		{ add_var(&variables,$3,TRUE,0); }
	|	DEL VARS '{' new_vars '}'	{ del_list_vars(&tmp_vars); }
	|	DEL '{' new_vars '}'	{ del_list_vars(&tmp_vars); }
	|	DEL VAR ID		{ del_var(&variables,$3); }
	|	DEL ID			{ del_var(&variables,$2); }
	|	LIST VARS		{ list_vars(variables); }
	|	ID DEPENDS ON '{' str_list '}' { add_fam($1, $1, $5); }
	|	LIST STRUCTURE		{ list_struct(variables, empty); }
	|	INFORMATIVITY ID	{ list_informativity($2); }
	|	SET PRINT INUM		{ print_short = $3; }
	|	QUALITATIVE ID		
{ var_type *v; FIND_VAR(v,$2); if (v!=NULL) v->ctype=ct_nominal; }
	|	QUANTITATIVE ID		
{ var_type *v; FIND_VAR(v,$2); if (v!=NULL) v->ctype=ct_contin; }

	|	COMPARE STRUCTURE ID ID
{ printf("Struct %s %s diff %6.3lf\n", $3, $4, compare_struct_dist($3, $4)); }
	|	COUNT LEAF ID		
{ printf("Leaves %s %d\n", $3, count_id_leaves($3)); }
	|	COUNT NODE ID		
{ double d;
  int n = count_id_inodes($3,&d);
  printf("INodes %s %d(%5.2lf)\n", $3, n, d); }

/****************************************************************************
  changing the description of variables
****************************************************************************/

	|	str_list IN '{' str_list '}'	{ add_qdesc_varlist($1, $4); }
	|	ID '/' id_num TO id_num		{ rename_desc($1, $3, $5); }
	|	id_num TO id_num
{ rename_desc(cfam->out->name,$1,$3); }
	|	ID OF ID IS '[' NUM ',' NUM ',' NUM ')'
{ add_fdesc_var($1,$3,left,$6,$8,$10,0); }
	|	ID OF ID IS '(' NUM ',' NUM ',' NUM']'
{ add_fdesc_var($1,$3,right,$6,$8,$10,0); }
	|	ID OF ID IS '(' NUM ',' NUM ',' NUM ',' NUM ')'
{ add_fdesc_var($1,$3,regular,$6,$8,$10,$12); }
	|	INTERVAL OF ID IS '{' num_list '}'   { add_idesc_var($3, $6); }
	|	DERIVE INTERVAL FOR var_list USING ID
{ derive_int_using_itable($6, $4); }
	|	LEARN INTERVAL FOR ID USING ID	
{
#ifdef CBIG
  learn_intervals($4, $6); 
#endif
}
	|	LIST DESCRIPTIONS	{ list_var_desc(variables); }
	|	PLOT DESCRIPTIONS	{ plot_var_desc(variables); }

/****************************************************************************
  defining and changing RULES
****************************************************************************/

	|	SET RULE encode_rules_type	{ encode_rules = $3; }
	|	FAM ID FOR ID IS '{' str_list '}' { add_fam($2, $4, $7); }
	|	LIST ALL FAM		{ list_fams(); }
	|	LIST FAM ID FOR ID	{ list_rules($3, $5); }
	|	LIST FAM
{
  if (cfam==NULL)
    printf("error: no FAM current\n");
  else
    list_rules(cfam->name, cfam->out->name);
}
	|	LIST
{
  if (cfam==NULL)
    printf("error: no FAM current\n");
  else
    list_rules(cfam->name, cfam->out->name);
}
	|	VECT
{
  if (cfam==NULL)
    printf("error: no FAM current\n");
  else
    list_rules_vect(cfam);
}
	|	LIST RULE
{
  if (cfam==NULL)
    printf("error: no FAM current\n");
  else
    list_rules(cfam->name, cfam->out->name);
}
	|	SORT RULE
{
  if (cfam==NULL)
    printf("error: no FAM current\n");
  else
    sort_rules(cfam);
}
	|	PLOT FAM
{
  if (cfam==NULL)
    printf("error: no FAM current\n");
  else
    plot_rules(cfam->name, cfam->out->name);
}
	|	SEL FAM ID FOR ID
{
  cfam = find_fam($5, $3);
}
	|	SHOW CURRENT FAM
{
  if (cfam==NULL)
    printf("no FAM current\n");
  else
    printf("FAM %s for %s\n", cfam->name, cfam->out->name);
}
	|	FAM ID FOR ID RULE NUM '=' ID
{
  add_rule_num($2, $4, $6, $8);
}
	|	RULE NUM '=' ID
{
  add_rule_num(cfam->name, cfam->out->name, (int) $2, $4);
}
	|	RULE NUM '=' NUM
{
  if (cfam==NULL)
    printf("error: no FAM current\n");
  else {
    Str255 s;
    sprintf(s, "%d", (int)$4);
    add_rule_num(cfam->name, cfam->out->name, (int) $2, s);
  }
}
	|	RULE '(' str_list ')'
{
  if (cfam==NULL)
    printf("error: no FAM current\n");
  else
    add_rule_str(cfam->name, cfam->out->name, $3);
}
	|	RULE TABLE prepare_rt '{' poss_nl rule_table poss_nl '}'
{
  set_apriory(cfam);
}
	|	SET APRIORY ID 
{
  set_apriory_id($3);
}
	|	RULE RULE TABLE prepare_rt '{' poss_nl i_rule_table poss_nl '}'
{
  set_apriory(cfam);
}
	|	DEL RULE INUM
{
  if (cfam==NULL)
    printf("error: no FAM current\n");
  else {
    int n = (int) $3;
    if (cfam->er == er_table) {
      if (n<0 || n>=cfam->n_rules)
	printf("error: illegal rule number %d\n", n);
      else 
	cfam->rtp[n] = undef;
    }
    else if (cfam->er == er_list) {
      rule_list *rl, *lrl;
      int i;
      if (n==0) cfam->lrule = cfam->lrule->next;
      else {
	rl = lrl = cfam->lrule;
	for (i=0; i<n; i++, rl=rl->next)
	  lrl =rl;
	lrl->next = rl->next;
      }
    }
  }
}
	|	DEL RULE '(' str_list ')'
{
  if (cfam==NULL)
    printf("error: no FAM current\n");
  else
    del_rule_str(cfam->name, cfam->out->name, $4);
}
	|	DEL RULE
{
  if (cfam==NULL)
    printf("error: no FAM current\n");
  else
    del_rules(cfam);
}
	|	IMP NUM '=' NUM
{
  if (cfam==NULL)
    printf("error: no FAM current\n");
  else
    if ($2<cfam->n_rules && $2>=0 && $4>=0 && $4<=1)
      cfam->imp[(int)$2] = $4;
    else
      printf("error: parameters out of range\n");
}
	|	RESET RULE USAGE FOR ID { reset_rule_usage($5);}
	|	RESET NUM 		{ printf("Number %lf\n", $2); }
	|	RESET ID		{ printf("Id %s\n", $2); }
	|	TEST DUPLICATE ID	{ test_duplicates($3); }
	|	MERGE DUPLICATE ID	{ merge_duplicates($3); }
	|	EXPAND RULE 		{ expand_rules(cfam); }
	|	EXPAND RULE num		{ extern double min_exp_w;
					  min_exp_w=$3;
					  expand_rules(cfam); }

/****************************************************************************
  RULE tables
****************************************************************************/

	|	LIST TABLE		{ list_tables(); }

	|	COPY ID TO TABLE ID	{ copy_to_table($2, $5); }
	|	SAVE RULE IN TABLE ID	
{
  if (cfam==NULL)
    printf("error: no FAM current\n");
  else
    copy_to_table(cfam->out->name, $5);
}
	|	DEL TABLE ID		{ del_table($3); }
	|	LIST TABLE ID		{ list_table($3); }
	|	COPY num '%' OF TABLE ID TO ID
{ copy_table($6, $8, $2); }
	|	SPLIT num '%' OF ID TO ID AND ID
{ split_rules($5, $7, $9, $2); }
	|	COPY RULE FROM ID TO ID	{ copy_table($4, $6, 100); }
	|	COMPARE ID TO TABLE ID	{ compare_table_var($5, $2); }
	|	COUNT RULE ID		{ count_id_rules($3); }
	|	LIST APRIORY		{ list_apriory(cfam); }
/****************************************************************************
  decision trees for RULEs
****************************************************************************/

	|	TREE ID		 	{ derive_dec_tree($2,0); }
	|	TREE ID INFORMATIVITY 	{ derive_dec_tree($2,0); }
	|	TREE ID GINI 		{ derive_dec_tree($2,1); }
	|	TREE ID GR 		{ derive_dec_tree($2,2); }
	|	TREE ID RELIEFF 	{ derive_dec_tree($2,3); }

	|	LIST TREE ID		{ list_dec_tree($3); }

/****************************************************************************
  RULE decomposition
****************************************************************************/

	|	SET REPEAT TEST NUM	{ repeat_tests = (int) $4; }
	|	SET SAVE TEST yesno     { save_tests = $4; }
	|	TEST DECOMPOSE ID FOR NUM TO NUM STEP NUM RULE OF ID
{
  test_decomposition($3, $12, (int)$5, (int)$7, (int)$9);
}
	|	SET DECOMPOSE LOCAL	{ dec_global = FALSE; }
	|	SET DECOMPOSE GLOBAL	{ dec_global = TRUE; }

	|	SET DECOMPOSE CM	{ dec_crit = c_cm; }
	|	SET DECOMPOSE DFC	{ dec_crit = c_dfc; }
	|	SET DECOMPOSE ERROR	{ dec_crit = c_error; }
	|	SET DECOMPOSE MNUMBER	{ dec_crit = c_m; }
	|	SET DECOMPOSE SDTIC	{ dec_crit = c_r; }
	|	SET DECOMPOSE DTIC	{ dec_crit = c_cr; }


	|	SET DECOMPOSE INUM	{ n_dis_lo = n_dis_hi = $3; }
	|	SET DECOMPOSE INUM TO INUM  { n_dis_lo = $3; n_dis_hi = $5; }
	|	SET DECOMPOSE INUM AND INUM
{ n_dis_lo = n_dis_hi = $3; n_ndis_lo = n_ndis_hi = $5; }
	|	SET DECOMPOSE INUM TO INUM AND INUM TO INUM
{ n_dis_lo = $3; n_dis_hi = $5; n_ndis_lo = $7; n_ndis_hi = $9; }
	|	TEST REDUNDANT ID	{ check_redundant($3, FALSE, NULL); }
	|	DEL REDUNDANT ID	{ check_redundant($3, TRUE, NULL); }
	|	DEL REDUNDANT INFORMATIVITY ID
{ del_redundant_sorted($4, M_INFORM, NULL); }
	|	DEL REDUNDANT GINI ID
{ del_redundant_sorted($4, M_GINI, NULL); }
	|	DEL REDUNDANT GR ID
{ del_redundant_sorted($4, M_GR, NULL); }
	|	DEL REDUNDANT RELIEFF ID
{ del_redundant_sorted($4, M_RELIEFF, NULL); }
	|	TEST REDUNDANT ID FOR var_list
					{ check_redundant($3, FALSE, $5); }
	|	DEL REDUNDANT ID FOR var_list
					{ check_redundant($3, TRUE, $5); }
	|	YLOG INFORMATIVITY 	{ log_inform = !log_inform; }
	|	YLOG RELIEFF 		{ log_relieff = !log_relieff; }
	|	YLOG GINI 		{ log_gini = !log_gini; }
	|	YLOG GR 		{ log_gr = !log_gr; }
	|	YLOG MDL 		{ log_mdl = !log_mdl; }
	|	SET DEBUG INFORMATIVITY yesno	{ deb_inform = $4; }
	|	SET DEBUG RELIEFF yesno	{ deb_relieff = $4; }
	|	SET DEBUG GINI yesno	{ deb_gini = $4; }
	|	SET DEBUG GR yesno	{ deb_gr = $4; }
	|	SET DM yesno		{ use_dm = $3; }

	|	DFC ID { printf("dfc %s is %d\n", $2, get_dfc_measure($2)); }
	|	SDTIC ID 
{ printf("sdtic %s is %5.3lf\n", $2, get_code_measure($2)); }
	|	DTIC ID 
{
  extern double get_dtic_measure(Str255 vname);
  printf("dtic %s is %5.3lf\n", $2, get_dtic_measure($2));
}

	|	DECOMPOSE INFORMATIVITY ID INUM { fast_decompose($3, M_INFORM, (int)$4); }
	|	DECOMPOSE GINI ID INUM { fast_decompose($3, M_GINI, $4); }
	|	DECOMPOSE GR ID INUM { fast_decompose($3, M_GR, $4); }
	|	DECOMPOSE RELIEFF ID INUM { fast_decompose($3, M_RELIEFF, $4); }
	|	DECOMPOSE MDL ID INUM { fast_decompose($3, M_MDL, $4); }


	|	SET DECOMPOSE DECOMPOSE INUM 
{
  /* this is to use new IM derivation procedure */
  extern char dec_speed;
  dec_speed = $4;
}

/****************************************************************************
  changing values of variables, option manipulation
****************************************************************************/

	|	LIST OPT		{ list_struct(variables, curropt); }
	|	LIST OPT ALL		{ list_options(); }
	|	DEL OPT ID		{ del_opt($3); }
	|	SET OPT ID		{ set_opt($3); }
	|	SEL OPT ID		{ select_opt($3); }

	|	ID '=' NUM
{
  Str255 s;
  if (pref_qual) {
    sprintf(s, "%d", (int)$3);
    set_var_val_q($1,s);
  }
  else set_var_val_num($1, $3);
}
	|	ID '=' '{' q_list '}'	{ set_var_val_qlist($1, $4); }
	|	ID '=' ID		{ set_var_val_q($1,$3); }
	|	EXPECT ID '=' NUM	{ set_var_expect($2, $4); }
	|	EVALUATE ID		{ evaluate($2); }
	|	EVALUATE ID FOR ID
{
  if (find_opt($4, TRUE) != NULL) {
    select_opt($4);
    evaluate($2);
    set_opt($4);
  }
}
	|	DERIVE var_list		{ derive_lvar($2); }
	|	DERIVE var_list FOR ID
{
  if (find_opt($4, TRUE) != NULL) {
    select_opt($4);
    derive_lvar($2);
    set_opt($4);
  }
}
  	|	DERIVE var_list FOR ALL
{
  list_of_opt *opt;
  for (opt=options; opt!=NULL; opt=opt->next) {
    select_opt(opt->name);
    derive_lvar($2);
    set_opt(opt->name);
  }
}
  	|	DERIVE var_list FOR '{' str_list '}'
{
  list_of_str *s;
  for (s=$5; s!=NULL; s=s->next)
    if (find_opt(s->str, TRUE) != NULL) {
      select_opt(s->str);
      derive_lvar($2);
      set_opt(s->str);
    }
}
	|	LIST ID FOR ALL		{ list_var_opt_all($2); }
	|	PLOT ID			{ plot_var_opt_all($2); }
	|	LIST ID FOR '{' str_list '}'	{ list_var_opt_str($2, $5); }
	|	UNDEF OPT
{
				/* BBB try not only to undefine the
                                   whole option, but also just
                                   internal nodes */
  undefine_opt();
}
	|	EVALUATE INUM {extern char evaltype; evaltype=$2;}
	|	SET EVALUATE eval_type	{ eval_method = $3; }
	|	INSTANCE TABLE ID FOR ID '{' poss_nl str_list '}'
					{ add_itable($3, $5, $8); }
	|	LIST INSTANCE TABLE	{ list_instance_tables(); }
	|	LIST INSTANCE TABLE ID	{ list_instance_table($4,FALSE); }
	|	DERIVE INSTANCE TABLE ID { list_instance_table($4,TRUE); }
	|	SEL INSTANCE NUM FROM ID
				      { load_instance_from_table((int)$3,$5);}
	|	LOAD INSTANCE FROM ID	{ load_instances_from_table($4); }

	|	DERIVE RULE TABLE ID FROM INSTANCE TABLE ID
				{ derive_rtable_from_itable($4, $8, TRUE); }
	|	DERIVE RULE FOR ID FROM INSTANCE TABLE ID
				{ derive_rtable_from_itable($4, $8, FALSE); }

	|	ID FROM INSTANCE TABLE ID
				{ discretize($1, $5); }

	|	SPLIT num '%' INSTANCE TABLE ID TO ID AND ID
{ split_instance_table($6, $8, $10, $2); }

/****************************************************************************
  learning fuzzy rules and descriptions
****************************************************************************/

	|	SET GA POPULATION INUM	{ ga_population_size = $4; }
	|	SET GA ITERATIONS INUM	{ ga_maxiter = $4; }
	|	SET GA DESCRIPTIONS POINTS '=' NUM	{ ga_desc_pts = $6; }
	|	SET GA WEIGHT POINTS '=' NUM	{ ga_weight_pts = $6; }
	|	SET GA PRINT FREQ '=' NUM	{ ga_print_freq = $6; }
	|	SET GA WIDTH  '=' NUM		{ ga_width = $5; }
	|	SET MAX PLOT ERROR '=' NUM	{ ga_max_error = $6; }
	|	SET ERROR ga_error_entry	{ ga_error_method = $3; }
	|	SET AND '=' fand		{ fa_method = $4; }
	|	SET GA MUTATION FAM '=' NUM	{ ga_mut_fams = $6; }
	|	SET GA MUTATION DESCRIPTIONS '=' NUM	{ ga_mut_desc = $6; }
	|	SET GA MUTATION WEIGHT '=' NUM	{ ga_mut_w = $6; }
	|	SET GA CROSSOVER FAM '=' NUM	{ ga_cross_fams = $6; }
	|	SET GA CROSSOVER DESCRIPTIONS '=' NUM	{ ga_cross_desc = $6; }
	|	SET GA CROSSOVER WEIGHT '=' NUM	{ ga_cross_w = $6; }
	|	SET GA LEARN FAM '=' yesno	{ ga_learn_fams = $6; }
	|	SET GA LEARN DESCRIPTIONS '=' yesno	{ ga_learn_desc = $6; }
	|	SET GA LEARN WEIGHT '=' yesno	{ ga_learn_w = $6; }
	|	SET GA POLICY '=' ga_policy_entry	{ ga_policy = $5; }
	|	LEARN FAM FOR ID		{}
	|	LEARN WEIGHT FOR ID		{}
	|	LEARN DESCRIPTIONS FOR ID	{}
	|	LEARN FOR ID			{}
	|	LEARN var_list			
{
#ifdef CBIG
  ga_learn($2);
#endif
}
	|	LEARN FAM ID FOR ID     	{}
	|	STAT FOR ID			{ stat_for_var($3); }
	|	STAT FOR INSTANCE TABLE ID	
{
  printf("Error: %10.5lf\n", stat_for_itable($5));
}
	|	PLOT GA ERROR			
{
#ifdef CBIG
  plot_ga_error();
#endif
}

/****************************************************************************
  node decomposition - constructive induction
****************************************************************************/

	|	TEST DECOMPOSE ID	{ bottom_up_decompose($3,TRUE,TRUE);}
	|	DECOMPOSE ID		
{ 
  if (dec_global) global_decompose($2,TRUE,FALSE);
  else bottom_up_decompose($2,TRUE,FALSE);
}
/*	|	TEST ID		{ global_decompose($2,TRUE,FALSE);} */
	|	COMPOSE ID		{ compose($2); }
	|	ONE DECOMPOSE ID	{ bottom_up_decompose($3,FALSE,FALSE);}
	|	DECOMPOSE TABLE ID ON var_list AND var_list
{ decomposition_table($3, $5, $7); }

	|	JOIN var_list OF ID TO ID
{ join_nodes($2,NULL,$4,$6,0); }
	|	JOIN var_list AND var_list OF ID TO ID
						{ join_nodes($2,$4,$6,$8,0); }
	|	PJOIN var_list OF ID TO ID
{ join_nodes($2,NULL,$4,$6,1); }
	|	PJOIN var_list AND var_list OF ID TO ID
						{ join_nodes($2,$4,$6,$8,1); }
	|	FJOIN
{ join_nodes(NULL,NULL,NULL,NULL,2); }

	|	SET HEURISTIC SPLIT yesno	{ heur_split = $4; }
	|	SET HEURISTIC COMPUTE yesno	{ heur_comp = $4; }
	|	SET DECOMPOSE DEBUG INUM	{ deb_dec = $4; }
	|	SET DECOMPOSE KNN NUM		{ krelieff = $4; }

	|	SET COLOR coloring_entry	{ coloring = $3; }
/*	|	SET MNUMBER NUM		{ mcriteria = $3; } */
	|	SET SAVE COLOR yesno		{ save_im = $4; }
	|	SET SAVE COLOR ID		{ strcpy(im_fname, $4); }
	|	ANALYZE COLOR			{ analyze_current_color(); }
	|	LIST COLOR			{ list_current_color(); }
	|	LIST GROUP		{ list_color_groups(); }
	|	LIST INSTANCE		{ list_unsure_instances(); }
	|	LIST DM			{ list_ndm(); }
	|	INSTANCE NUM TO NUM	{ set_instance_to_color((int)$2, (int)$4); }

	|	SET COLOR STRING	{ set_color($3); }
	|	SURE COLOR		{ set_sure_color(); }
	|	NEED COLOR		{ set_need_color(); }
	|	SAME NUM NUM		{ set_same_color((int)$2,(int)$3); }
	|	DIFFERENT NUM NUM	{ set_diff_color((int)$2,(int)$3); }

	|	SET COLOR ON		{ extern char use_evidence;
					  use_evidence = TRUE;
					}

/****************************************************************************
  noise handling experiments
****************************************************************************/

	|	num '%' NOISE TO ID     { add_noise($5, $1); }
	|	SET MNUMBER num		{ m_param = $3; }

	|	TEST ERROR	{ get_class_error();}
	|	TEST num	{ test($2,0);}
	|	TEST num num	{ test($2,$3);}
/*	|	TEST EVALUATE   { test_classify(); } */
	|	TEST EVALUATE   { double e, m; 
				  int c=get_partition_error(&e,&m);
				printf("FIN c=%d e=%6.2lf %5.2lf\n", c, e,m);}
	|	MERGE INUM INUM	{ merge_colors($2, $3, -1.); }
	|	MERGE		{ repeat_merge(-1); }
	|	MERGE INUM	{ repeat_merge($2); }
	|	TEST MERGE num num STEP num	{ test_merge($3, $4, $6); }
	|	MERGE EVALUATE  { printf("Error: %7.3lf\n", 
					 get_class_error()); }
/*	|	CLUSTER NUM	{ cluster_fixed_colors((int)$2); } */
	|	PAIRS		{ int i, j; double d;
				  select_pair(&i, &j, &d);
				}
	|	SET NOISE MNUMBER       { noise_handling = n_mprob; }
	|	SET NOISE MNUMBER CM    { noise_handling = n_mprob;
					  m_dec = m_cm; }
	|	SET NOISE MNUMBER ERROR { noise_handling = n_mprob;
					  m_dec = m_error; }
	|	SET NOISE MNUMBER CM ERROR { noise_handling = n_mprob;
					  m_dec = m_cm_error; }
	|	SET NOISE INFORMATIVITY	{ noise_handling = n_entropy; }
	|	SET NOISE LAPLACE	{ noise_handling = n_laplace; }
	| 	SET NOISE CLUSTER	{ noise_handling = n_cluster; }

	|	SET DONTCARE	{ dont_care = dc_dont_care; }
	|	SET DONTKNOW	{ dont_care = dc_dont_know; }
	|	SET APRIORY	{ dont_care = dc_apriory; }

	|	SEED num	{ srand((int)$2); }

	|	SET DISTRIBUTION	{ use_distribution=!use_distribution; }
	|	SET CV 			{ test_type = t_cv; }
	|	SET CV INUM		{ int extern n_folds;
					  test_type = t_cv;
					  n_folds = (int)$3;
					}
	|	SET SAMPLE		{ test_type = t_sample; }
	|	SET SAMPLE num num INUM	{ extern int n_folds;
					  extern double p_sample_l, p_sample_mut;
					  test_type = t_sample;
					  p_sample_l = $3; p_sample_mut = $4;
					  n_folds = (int)$5;
					}

	|	SET CV INUM FOR ID	{ prepare_cv($5, $3); }
	|	SPLIT ID TO ID AND ID USING INUM
				{ split_using_cv($2, $4, $6, $8); }
	|	SET NOISE LOCAL { char extern sel_mfreq; sel_mfreq=TRUE; }
	|	SET NOISE GLOBAL { char extern sel_mfreq; sel_mfreq=FALSE; }


/****************************************************************************
  input / output
****************************************************************************/

	|	LOAD ID
{
  FILE *f;

  if (myno_files == MY_MAX_FILES) {
    printf("error: no of files to be loaded exceeded the limit\n");
    exit(0);
  }
  printf("load %s\n", $2);
  f = fopen($2, "r");
  if (f == NULL) {
    printf("error: could not open %s as cmd file\n", $2);
  }
  else {
    myno_files++;
    myfiles[myno_files].f = yyin;
    myfiles[myno_files].lineno = yylineno;
    yyin = f;
    yylineno = 1;
    interact_mode = FALSE;
  }
}

/* fast loading of rules which opens the file ID and loads the rules
   that are in the index form. */

	|	LOAD RULE ID		
{ load_irules($3, cfam, 0); set_apriory(cfam); }
	|	LOAD RULE ID USING INUM	
{ load_irules($3, cfam, $5); set_apriory(cfam); }

/****************************************************************************
  store commands
****************************************************************************/
	|	SAVE STRUCTURE ID	{ save_struct($3); }
	|	SAVE DESCRIPTIONS ID	{ save_des($3); }
	|	SAVE RULE ID		{ save_rules($3); }
	|	SAVE OPT ID		{ save_opt($3); }
	|	SAVE ID			{ save_all($2); }

	|	WRITE ID TO ID		{ save_c45_rules($2,$4);
					  save_c45_names($2,$4);
					}
	|	WRITE ID ID TO ID	
{
  Str255 s1, s2;
  save_c45_names($2,$5);
  strcpy(s1,$5);
  save_c45_rules($2,strcat($5,".data"));
  save_c45_rules($3,strcat(s1,".test"));
}
	|	WRITE RULE ID TO ID	{ save_c45_rules($3,$5); }
	|	WRITE NAMES ID TO ID	{ save_c45_names($3,$5); }
	|	WRITE DOT ID TO ID	{ save_dot_struct($3,$5); }

/****************************************************************************
  other commands
****************************************************************************/
	|	SEL ID
{				/* select fam or option -> FIFO  */
  var_type *v;
  list_of_opt *opt;

  if (((v = find_var(variables,$2)) != NULL) &&
      find_fam($2, v->name) != NULL)
    cfam = find_fam($2, v->name);
  else if (find_opt($2,TRUE)) {
    opt = find_opt($2,TRUE);
    if (opt != NULL)
      select_opt($2);
    else
      printf("error: not found, neither a FAM nor an option\n");
  }
}
	|	SHOW
{
  printf("GENETIC ALGORITHM:\n");
  printf("population size    %d\n", ga_population_size);
  printf("max iterations     %d\n", ga_maxiter);
  printf("description points %d\n", ga_desc_pts);
  printf("weight points      %d\n", ga_weight_pts);
  printf("print freq         %d\n", ga_print_freq);
  printf("width (desc enlar) %lf\n", ga_width);
  printf("policy (cross,mut) %s\n",
	 ga_policy==gap_or?"or":(ga_policy==gap_and?"and":"custom"));
  printf("max plot error     %lf\n", ga_max_error);
  printf("error measure      %s\n",
	 (ga_error_method==gae_abs)?"abs":
	 (ga_error_method==gae_norm)?"norm":
	 (ga_error_method==gae_sqr)?"sqr":
	 (ga_error_method==gae_perc)?"perc":
	 (ga_error_method==gae_sqrperc)?"sqr perc":"max perc");
  printf("mutation probab    %5.2lf/L (fams) %5.2lf/L (desc) %5.2lf/L (w)\n",
	 ga_mut_fams, ga_mut_desc, ga_mut_w);
  printf("crossover probab   %5.2lf   (fams) %5.2lf   (desc) %5.2lf   (w)\n",
	 ga_cross_fams, ga_cross_desc, ga_cross_w);
  printf("learn              %s   (fams) %s   (desc) %s   (w)\n",
	 ga_learn_fams?" yes ":" no  ",
	 ga_learn_desc?" yes ":" no  ",
	 ga_learn_w?" yes ":" no  ");

  printf("\nDEBUGGING:\n");
  printf("evaluate           %s\n", debug_e ? "yes" : "no");
  printf("lexical analysis   %s\n", debug_l ? "yes" : "no");
  printf("ga learning        %s\n", debug_g ? "yes" : "no");
  printf("descriptor constr  %s\n", debug_g ? "yes" : "no");
  printf("print length       %d\n", print_short);

  printf("\nEVALUATION:\n");
  printf("fuzzy and          %s\n", (fa_method==famin)?"min" :
	 ((fa_method==famax)?"max":"mult"));

/*  printf("\nFUZZY IDENTIFICATION:\n");
  printf("init num clusters  %d\n", cl_K);
  printf("fuzziness          %5.3lf\n", cl_m);
  printf("merging treshold   %5.3lf\n", cl_gamma); 
  printf("clustering e       %7.5lf\n", cl_e);
  printf("redo the loop      %d\n", cl_redo);
  printf("sigm memberships   %s\n", cl_sigm ? "yes" : "no");
  printf("deb clust proto    %s\n", deb_cl_cp ? "yes" : "no");
  printf("data and partit    %s\n", deb_cl_pm ? "yes" : "no");
  printf("identif matrices   %s\n", deb_cl_f ? "yes" : "no"); */

  printf("\nDECOMPOSITION:\n");
  printf("coloring           %s\n", coloring==col_optimal ? "opt" :
	 coloring==col_heuristic ? "heuristic" : "ga");
  printf("m-criteria         %lf\n", mcriteria);
  printf("heuristic compute  %s\n", heur_comp ? "yes" : "no");
  printf("heuristic split    %s\n", heur_split ? "yes" : "no");
  printf("dec debug          %d\n", deb_dec);
  printf("k for NN relieff   %d\n", krelieff);
  printf("log informativity  %s\n", log_inform ? "yes" : "no");
  printf("log relieff        %s\n", log_relieff ? "yes" : "no");
  printf("log gini           %s\n", log_gini ? "yes" : "no");
  printf("log gain-ratio     %s\n", log_gr ? "yes" : "no");
  printf("debug inform       %s\n", deb_inform ? "yes" : "no");
  printf("debug relieff      %s\n", deb_relieff ? "yes" : "no");
  printf("debug gini         %s\n", deb_gini ? "yes" : "no");
  printf("debug gain-ratio   %s\n", deb_gr ? "yes" : "no");
  printf("save im            %s\n", save_im ? "yes" : "no");
  printf("save im file name  %s\n", im_fname);
  printf("bound              %d to %d\n", n_dis_lo, n_dis_hi);

  printf("\nEVALUATION / DATA REPRESENTATION:\n");
  printf("encoding rules     %s\n", encode_rules==er_table? "as table" :
	 encode_rules==er_list? "as list" : "as tree");
  printf("use dm             %s\n", use_dm?"yes":"no");
  printf("eval method        %s\n", eval_method==e_crisp? "crisp" :
	 eval_method==e_fuzzy? "fuzzy" :
	 eval_method==e_numerical? "numerical" : "interval");
  printf("entry preference   %s\n", pref_qual?"qualitative":"quantitative");
  printf("test repetition    %d\n", repeat_tests);
  printf("var prefix         %s\n", prefix_var);
  printf("opt prefix         %s\n", prefix_opt);
}
	|	DEBUG EVALUATE 	{ debug_e = TRUE;}
/*	|	DEBUG LEX	{ debug_l = TRUE;} */
	|	DEBUG GA	{ debug_g = TRUE;}
	|	DEBUG DESCRIPTIONS	{ debug_d = TRUE;}
	|	PREFER QUALITATIVE { pref_qual = TRUE; }
	|	PREFER QUANTITATIVE { pref_qual = FALSE; }
	|	VAR PREFIX ID	{ strcpy(prefix_var, $3); iprefix_var=1; }
	|	OPT PREFIX ID	{ strcpy(prefix_opt, $3); iprefix_opt=1; }
	|	WAIT NUM
{
/*  struct tms bbuf, ebuf;
  times(&bbuf);
  times(&ebuf);
  while (((float)(ebuf.tms_utime-bbuf.tms_utime)/60) < $2)
    times(&ebuf); */
}
	|	YECHO STRING		{ printf("%s\n", $2); }
	|	YLOG STRING		{ fprintf(lfile, "%s\n", $2); }
	|	COLOR			{ color_graph(); }
	|	QUIT			{ exit(1); }
;

yesno:		YES	{ $$=TRUE; }
	|	NO	{ $$=FALSE; }
;

fand:		MIN	{ $$=famin; }
	|	MAX	{ $$=famax; }
	|	MULT	{ $$=famult; }
;

eval_type:	FUZZY	{ $$=e_fuzzy; }
	|	CRISP	{ $$=e_crisp; }
	|	INTERVAL { $$=e_interval; }
;

coloring_entry:	OPT	 	{ $$=col_optimal; }
	|	HEURISTIC	{ $$=col_heuristic; }
	|	GA		{ $$=col_ga; }
	|	SET		{ $$=col_set; }
/*	|	SOFT		{ $$=col_soft; } */

ga_error_entry:	GAABS	{ $$=gae_abs; }
	|	LSQR	{ $$=gae_sqr; }
	|	'%'	{ $$=gae_perc; }
	|	LSQR '%'	{ $$=gae_sqrperc; }
	|	MAX '%'	{ $$=gae_maxperc; }
	|	NORM	{ $$=gae_norm; }
;

ga_policy_entry:	AND	{ $$=gap_and; }
	|		OR	{ $$=gap_or; }
	|		CUSTOM	{ $$=gap_custom; }
;

new_vars:	new_vars ',' new_var	{}
	|	new_var			{}
;

new_var:	ID
{
  add_var(&tmp_vars, $1, FALSE, 0);
  add_var(&variables, $1, FALSE, 0);
}
;


str_list:	str_list_rev
{
  list_of_str *sl, *slt;

  slt = NULL;
  for (sl=$1; sl!=NULL; sl=sl->prev) {
    sl->next = slt;
    slt = sl;
  }
  $$ = slt;
}
;

str_list_rev:	str_list_rev separator_nl ID
{
  $$ = (list_of_str *) malloc(sizeof(*$$));
  if ($$==NULL) printf("out of mem\n");
  strcpy($$->str, $3);
  $$->prev = $1;
}
	|	str_list_rev separator_nl NUM
{
  Str255 s;
  $$ = (list_of_str *) malloc(sizeof(*$$));
  if ($$==NULL) printf("out of mem\n");
  sprintf(s, "%lf", $3);
  strcpy($$->str, s);
  $$->prev = $1;
}
	|	str_list_rev separator_nl INUM
{
  Str255 s;
  $$ = (list_of_str *) malloc(sizeof(*$$));
  if ($$==NULL) printf("out of mem\n");
  sprintf(s, "%d", (int)$3);
  strcpy($$->str, s);
  $$->prev = $1;
}
	|	str_list_rev separator_nl '?'
{
  $$ = (list_of_str *) malloc(sizeof(*$$));
  if ($$==NULL) printf("out of mem\n");
  strcpy($$->str, "?");
  $$->prev = $1;
}
	|	ID
{
  $$ = (list_of_str *) malloc(sizeof(*$$));
  if ($$==NULL) printf("out of mem\n");
  strcpy($$->str, $1);
  $$->prev = NULL;
}
	|	'?'
{
  $$ = (list_of_str *) malloc(sizeof(*$$));
  if ($$==NULL) printf("out of mem\n");
  strcpy($$->str, "?");
  $$->prev = NULL;
}
	|	NUM
{
  Str255 s;
  $$ = (list_of_str *) malloc(sizeof(*$$));
  if ($$==NULL) printf("out of mem\n");
  sprintf(s, "%d", (int)$1);
  strcpy($$->str, s);
  $$->prev = NULL;
}
	|	INUM
{
  Str255 s;
  $$ = (list_of_str *) malloc(sizeof(*$$));
  if ($$==NULL) printf("out of mem\n");
  sprintf(s, "%d", (int)$1);
  strcpy($$->str, s);
  $$->prev = NULL;
}
;

num_list:	num_list_rev
{
  list_of_num *sl, *slt;

  slt = NULL;
  for (sl=$1; sl!=NULL; sl=sl->prev) {
    sl->next = slt;
    slt = sl;
  }
  $$ = slt;
}
;

num_list_rev:	num_list_rev separator_nl num
{
  $$ = (list_of_num *) malloc(sizeof(*$$));
  $$->num = $3;
  $$->prev = $1;
}
	|	num
{
  $$ = (list_of_num *) malloc(sizeof(*$$));
  $$->num = $1;
  $$->prev = NULL;
}
;

separator_nl:	',' | ';' | '\n' | ;
separator:	',' | ';' | ;
poss_nl:	'\n' | ;

var_list:	var_list_rev
{
  list_of_vars *sl, *slt;

  slt = NULL;
  for (sl=$1; sl!=NULL; sl=sl->prev) {
    sl->next = slt;
    slt = sl;
  }
  $$ = slt;
}
;

var_list_rev:	var_list_rev ',' ID
{
  var_type *v;
  v = find_var(variables,$3);
  if (v==NULL) {
    printf("error: %s not found\n", $3);
    $$ = $1;
  }
  else {
    $$ = (list_of_vars *) malloc(sizeof(*$$));
    $$->var = v;
    $$->prev = $1;
  }
}
	|	ID
{
  var_type *v;
  v = find_var(variables,$1);
  if (v==NULL) {
    printf("error: %s not found\n", $1);
    $$ = NULL;
  }
  else {
    $$ = (list_of_vars *) malloc(sizeof(*$$));
    $$->var = v;
    $$->prev = NULL;;
  }
}

q_list:		q_list ',' q_list_el
{
  $3->next = $1;
  $$ = $3;
}
	|	q_list_el
{
  $$ = $1;
}
;

q_list_el:	ID '/' NUM
{
  $$ = (list_of_q *) malloc(sizeof(*$$));
  strcpy($$->desc, $1);
  $$->degree = $3;
  $$->next = NULL;
}
	|	NUM '/' NUM
{
  Str255 s;
  $$ = (list_of_q *) malloc(sizeof(*$$));
  sprintf(s, "%d", (int)$1);
  strcpy($$->desc, s);
  $$->degree = $3;
  $$->next = NULL;
}
;

/*****************************************************************************
entering rules as a table
*****************************************************************************/

                                     /****************************************
                                     using the names
                                     ****************************************/


rule_table:	rule_table '\n' att_entry
	|	att_entry
;

att_entry:	att_entry separator att_id
	|	att_id
;

att_id:		id_num
{
  int j;

  if (att_num == cfam->n_in) {
    for (j=0; j<cfam->out->ndesc 
	 && strcmp(cfam->out->desc[j].name, $1);
	 j++);
    if (j<cfam->out->ndesc) set_rule(cfam, att, j, man);
    else {
      printf("error %d: %s not legal qval for %s\n", 
	     yylineno, $1, cfam->out->name);
      exit(1);
    }

    att_num = 0;
  }
  else {
    if (!strcmp($1,"?")) {
      att[att_num] = CUNDEF;
    }
    else {
      for (j=0; j<cfam->in[att_num]->ndesc 
	   && strcmp(cfam->in[att_num]->desc[j].name, $1);
	   j++);
      if (j<cfam->in[att_num]->ndesc) att[att_num] = j;
      else {
	printf("error_ %d: %s not legal qval for %s\n", 
	       yylineno, $1, cfam->in[att_num]->name);
      }
    }
    att_num++;
  }
}
;

                                     /****************************************
                                     using the indices
                                     ****************************************/

i_rule_table:	i_rule_table '\n' i_att_entry
	|	i_att_entry
;

i_att_entry:	i_att_entry separator i_att_id
	|	i_att_id
;

i_att_id:	INUM
{
  int j;

  if (att_num == cfam->n_in) {
    if ($1<0 || $1>=cfam->out->ndesc) {
      printf("error i %d: index %d not legal for %s\n", 
	     yylineno, $1, cfam->out->name);
      exit(1);
    }
    set_rule(cfam, att, $1, man);
    att_num = 0;
  }
  else {
/*    if ($1==CUNDEF)
      att[att_num] = $1; */
    if ($1<0 || $1>=cfam->in[att_num]->ndesc) {
      printf("error j %d: index %d not legal for %s\n", 
	     yylineno, $1, cfam->in[att_num]->name);
      exit(0);
    }
    else
      att[att_num] = $1;
    att_num++;
  }
}
;

                                     /****************************************
                                     utilities
                                     ****************************************/

id_num:	ID		{ strcpy($$,$1); }
	|	NUM	{ Str255 s; sprintf(s, "%lf", (int)$1); strcpy($$,s);}
	|	INUM	{ Str255 s; sprintf(s, "%d", (int)$1); strcpy($$,s);}
	|	'?'	{ strcpy($$,"?"); }
;

num:	NUM		{ $$ = $1; }
	|	INUM	{ $$ = $1; }

prepare_rt:
{
  if (cfam==NULL)
    printf("error: no FAM current\n");
  else {
    if (att!=NULL) free(att);
    att = c_vector(cfam->n_in);
  }
}
;

encode_rules_type:	TABLE	{ $$ = er_table; }
	|		LIST	{ $$ = er_list; }
	|		TREE	{ $$ = er_tree; }

/*****************************************************************************
PROGRAMS SECTION
*****************************************************************************/

%%

yyerror()
{
  if (!interact_mode) {
    printf("Syntax error detected in line %d\n",yylineno);
    ++n_errors;
  }
}
