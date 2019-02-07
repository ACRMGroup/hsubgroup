#COPT	= -ansi -Wall -pedantic -O3
COPT	= -ansi -Wall -pedantic -g
LINK1	= 
LINK2	=
CC	= cc -I$(HOME)/include -L$(HOME)/lib

EXE	= hsubgroup
OFILES	= hsubgroup.o sophie.o fullmatrix.o

$(EXE) : $(OFILES) $(LFILES)
	$(CC) $(COPT) -o $(EXE) $(OFILES) $(LFILES) -lbiop -lgen -lm -lxml2

.c.o :
	$(CC) $(COPT) -o $@ -c $<

clean :
	/bin/rm -f $(EXE) $(OFILES) $(LFILES)

test : $(EXE)
	(cd t; ./test.sh)

