COPT  = -ansi -Wall -pedantic -L$(HOME)/lib -I$(HOME)/include
LIBS  = -lbiop -lgen -lm -lxml2
CC    = gcc

OFILES = hsubgroup.o sophie.o

hsubgroup : $(OFILES)
	$(CC) $(COPT) -o hsubgroup $(OFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -o $@ -c $<

clean :
	rm -f $(OFILES)


