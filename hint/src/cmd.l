/****************************************************************************
cmd.l

Lexical analyzer
****************************************************************************/


%p 6000
%a 8000
%o 15000

%{
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define GL extern
#include "sds.h"
#include "y.tab.h"

#define ECH if (debug_l) ECHO
#define RET(x) {if (debug_l) ECHO; return(x);}

%}

delim    [ \t]
ws       {delim}+
letter   [A-Za-z]
digit    [0-9]
id       ({letter}|(\/))({letter}|{digit}|(\.|\+|\-)|'|(_)|(\/))*
string   \"(.)*\"
comment	 #((#)|{id}|(.))*\n
echo	 echo{delim}+((#)|{id}|(.))*\n

inumber  (\-)?{digit}+
number   (\-)?{digit}+(\.{digit}+)?([eE][+\-]?{digit}+)?

%%

{ws}		{ ECH; }
\n		RET('\n');

";"		RET(';');
","		RET(',');
":"		RET(':');
"("		RET('(');
")"		RET(')');
"["		RET('[');
"]"		RET(']');
"{"		RET('{');
"}"		RET('}');
"\""            RET('"');
"="		RET('=');
"/"		RET('/');
"%"		RET('%');
"?"		RET('?');


add		RET(ADD);
all		RET(ALL);
and		RET(AND);
apriory		RET(APRIORY);
cm		RET(CM);
cross		RET(CROSSOVER);
comp(are)?	RET(COMPARE);
copy		RET(COPY);
count		RET(COUNT);
color           RET(COLOR);
cv		RET(CV);
deb(ug)?	RET(DEBUG);
dec(ompose)?	RET(DECOMPOSE);
(rm|del(ete)?)	RET(DEL);
dist(ribution)?	RET(DISTRIBUTION);
duplicate(s)?	RET(DUPLICATE);
dfc		RET(DFC);
dot		RET(DOT);
dontcare	RET(DONTCARE);
dontknow	RET(DONTKNOW);
dot		RET(DOT);
depends		RET(DEPENDS);
diff(erent)	RET(DIFFERENT);
dm		RET(DM);
dtic		RET(DTIC);
sdtic		RET(SDTIC);
echo		RET(YECHO);
error		RET(ERROR);
fam(s)?		RET(FAM);
for		RET(FOR);
from		RET(FROM);
gini		RET(GINI);
global		RET(GLOBAL);
gain		RET(GR);
in		RET(IN);
inf(or(m(ativity)?)?)?		RET(INFORMATIVITY);
is		RET(IS);
int(erval)?	RET(INTERVAL);
ins(tance)?	RET(INSTANCE);
join		RET(JOIN);
pjoin		RET(PJOIN);
fjoin		RET(FJOIN);
(ls|list)	RET(LIST);
learn		RET(LEARN);
load		RET(LOAD);
log		RET(YLOG);
m		RET(MNUMBER);
max		RET(MAX);
mdl		RET(MDL);
merge		RET(MERGE);
min		RET(MIN);
mult		RET(MULT);
noise		RET(NOISE);
name(s)?	RET(NAMES);
of		RET(OF);
on		RET(ON);
print		RET(PRINT);
quit		RET(QUIT);
nominal		RET(QUALITATIVE);
quant(itative)?	RET(QUANTITATIVE);
red(undant)?	RET(REDUNDANT);
relieff		RET(RELIEFF);
reset		RET(RESET);
rule(s)?	RET(RULE);
save		RET(SAVE);
sel(ect)?	RET(SEL);
seed		RET(SEED);
set		RET(SET);
split		RET(SPLIT);
stat(istics)?	RET(STAT);
show		RET(SHOW);
sure		RET(SURE);
struct(ure)?	RET(STRUCTURE);
to		RET(TO);
table		RET(TABLE);
test		RET(TEST);
tree		RET(TREE);
undef		RET(UNDEF);
using		RET(USING);
var(iable)?	RET(VAR);
var(s|iables)	RET(VARS);
wait		RET(WAIT);
weight(s)?	RET(WEIGHT);
write?		RET(WRITE);
yes		RET(YES);

{id}            {
  		  ECH;
                  strcpy(yylval.strV,yytext);
                  return(ID);
                }

{inumber}        {
                  ECH;
                  yylval.intV = atoi((const char *) yytext);
                  return(INUM);
                }

{number}        {
                  ECH;
                  yylval.doubV = atof((const char *) yytext);
                  return(NUM);
                }

{string}        {
  		  ECH;
                  strcpy(yylval.strV,yytext+1);
		  yylval.strV[strlen(yylval.strV)-1] = '\0';
                  return(STRING);
                }

{echo}        {
  		  ECH;
                  printf("%s", &yytext[5]);
                  return(COMMENT);
                }

{comment}	{
		  ECH;
/*		  printf("comment <%s>\n", yytext); */
/*		  return(COMMENT); */
		}

%%
