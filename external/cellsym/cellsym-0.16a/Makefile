# Adjust the next two lines so that libspg can be found
#LDFLAGS=-L/usr/local/shared/spglib/1.9.4/lib
CPPFLAGS=-I/usr/local/spglib/1.9.4/include

# Warning: dependences on .h files are not given. If you change
# these, type "make clean; make"

OBJS=cellsym.o cell_read.o periodic_table.o pdb_write.o \
     cell_write.o pdb_read.o basis.o ident_sym.o shelx_read.o \
     shelx_write.o


# Linux / IA32 / gcc
# Do not define SPGLIBVER unless your version of spglib is >=1.8.3
CFLAGS=-Wall -D_FILE_OFFSET_BITS=64 -DSPGLIBVER

CC=cc

cellsym: cellsym.h $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o cellsym $(OBJS) -lsymspg -lm

cellsym.o: cellsym.h


clean:
	rm cellsym $(OBJS)

