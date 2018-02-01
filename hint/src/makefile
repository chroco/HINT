#########################################################################
# Makefile for hint, (C) 1995, 96 by Blaz Zupan

#########################################################################
# HINTVER defines whether to include all the code (CBIG) or just a part
# that is directly related to decomposition (CSMALL)

HINTSIZE=small
#HINTSIZE=big

#########################################################################
# Compile on hp or sun machines, usign bound checking or not

#MACHINE=hp
MACHINE=sun
#MACHINE=sgi

BOUNDS_CHECK=no
#BOUNDS_CHECK=yes

#########################################################################

ifeq ($(HINTSIZE),small)
HINTVER=-DCSMALL
else
HINTVER=-DCBIG
endif

ifeq ($(MACHINE),hp)
  NR = /users/blaz/lib/nrc/double/
  ifeq ($(BOUNDS_CHECK),no)
#    CC     = cc -O +DS770
#    CC     = cc -g
#    CC     = cc +O2 +DS770
#    CC     = cc 
    CC = cc -O
    OPTS = -Aa -Dhp -DOPTIMIZE $(HINTVER)
#    OPTS = -Dhp -DOPTIMIZE $(HINTVER)
  else
    CC = gbcc 
    OPTS = -Dhp -DOPTIMIZE $(HINTVER)
  endif

else 
ifeq ($(MACHINE),sun)
  OPTS = -Dsun -DOPTIMIZE $(HINTVER)
  ifeq ($(BOUNDS_CHECK),no)
    CC = gcc -O
  else
    CC = /usr/local/bin/bgcc/xgcc -B/usr/local/bin/bgcc/ \
         -fbounds-checking -fno-builtin
  endif

else
ifeq ($(MACHINE),sgi)
  NR = /user/blaz/lib/nrc/double/
  OPTS = -Dsgi -DOPTIMIZE $(HINTVER)
  ifeq ($(BOUNDS_CHECK),no)
#    CC = gcc -O3
    CC = cc -O2
#    CC = cc
  else
    CC = gcc -fbounds-checking -fno-builtin
  endif
endif
endif
endif

#########################################################################
# libraries, objects

ifeq ($(HINTSIZE),small)
  LIBS = -ll -ly -lm
else
  LIBS = -ll -ly libpga.a -lm
endif

FDLIB = -L/usr/local/lib -lforms -lXpm -L/usr/lib/X11R5 -lX11 -lm

ifeq ($(HINTSIZE),small)
  OBJ_HINT = struct.o desc.o rules.o opt.o eval.o \
	misc.o io.o decomp.o color.o heuristics.o noise.o id3.o \
	my_mnbrak.o y.tab.o lex.yy.o
else
  OBJ_HINT = struct.o desc.o rules.o opt.o eval.o \
	misc.o learn.o io.o fuzzy.o decomp.o color.o heuristics.o noise.o \
	my_mnbrak.o \
	$(NR)golden.o \
	id3.o y.tab.o lex.yy.o
endif

OBJS = hint.o $(OBJ_HINT)
OBJS_PR = rplotf.o rplot.o $(OBJ_HINT)
OBJS_PF = fplot.o $(OBJ_HINT)

ifeq ($(HINTSIZE),small)
  NROBJ =
else
  NR  = ../lib/nrc/double/
  NROBJ = $(NR)nrutil.o $(NR)ludcmp.o $(NR)lubksb.o $(NR)jacobi.o
endif

#########################################################################
# dependencies

hint: $(OBJS)
	$(CC) $(OPTS)  -o hint $(OBJS) $(NROBJ) $(LIBS)

y: cmd.y
	yacc -v -d cmd.y
	$(CC) -c $(OPTS) y.tab.c
	$(CC) $(OPTS)  -o hint $(OBJS) $(NROBJ) $(LIBS)

rplot:	$(OBJS_PR)
	$(CC) $(OPTS) -s -o rplot $(OBJS_PR) $(LIBS) $(FDLIB)

fplot:	$(OBJS_PF)
	$(CC) $(OPTS) -o fplot $(OBJS_PF) $(NROBJ) $(LIBS) $(FDLIB)

fplot.o: fplot.c
	$(CC) -I/usr/local/include -c $(OPTS) $<

.c.o:
	$(CC) -c $(OPTS) $<

y.tab.o: cmd.y
	yacc -d cmd.y
	$(CC) -c $(OPTS) y.tab.c

# lex.yy.o: sds.h cmd.l y.tab.h
lex.yy.o: sds.h cmd.l
	lex cmd.l
	$(CC) -c $(OPTS) lex.yy.c

y.tab.h: cmd.y sds.h
	yacc -d cmd.y
	$(CC) -c $(OPTS) y.tab.c

hint.o:		sds.h
struct.o:	sds.h
desc.o:		sds.h
rules.o:	sds.h
opt.o:		sds.h
eval.o:		sds.h
misc.o:		sds.h
learn.o:	sds.h
io.o:		sds.h
fuzzy.o:	sds.h	
fplot.o:	sds.h colormap.h
decomp.o:	sds.h
color.o:	sds.h
heuristics.o:	sds.h
id3.o:		sds.h
noise:		sds.h

cp:
	cp *.c *.l *.y *.h *.s *.pl *.fd tmp* todo gen extract makefile bac

hp:
	cp *.c *.l *.y *.h 00 hpstuff
	tar cvf hpstuff.tar hpstuff
	gzip hpstuff.tar

wc:
	wc -l hint.c struct.c desc.c rules.c opt.c eval.c learn.c \
              misc.c io.c fuzzy.c decomp.c color.c heuristics.c noise.c

clean:
	rm *.o hint
