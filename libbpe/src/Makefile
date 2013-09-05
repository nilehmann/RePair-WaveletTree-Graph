CC=g++
LIBCDS=../../libcds
LIB=$(LIBCDS)/lib/libcds.a

OBJ=obj

DIRS+=$(OBJ)

TARGETS+=$(OBJ)/bpe.o
TARGETS+=$(OBJ)/mylist.o
TARGETS+=$(OBJ)/heap.o
TARGETS+=$(OBJ)/records.o
TARGETS+=$(OBJ)/hash.o
TARGETS+=$(OBJ)/array.o
TARGETS+=$(OBJ)/maloc.o
TARGETS+=$(OBJ)/factorization.o
TARGETS+=$(OBJ)/factorization_var.o

CFLAGS= -ggdb -m32 -fpermissive  -I$(LIBCDS)/includes/
#CFLAGS= -O9 -DNDEBUG -m32 -fpermissive  -I$(LIBCDS)/includes/

BPE=lib_bpe.a

all: lib

# pattern rule for all objects files
$(OBJ)/%.o: %.c *.h
	@echo " [C++] Compiling $<"
	@$(CC) -c $(CFLAGS) $< -o $@

$(OBJ)/%.o: %.cpp *.h
	@echo " [C++] Compiling $<"
	@$(CC) -c $(CFLAGS) $< -o $@


lib: ${DIRS} ${TARGETS}
	@echo " [LIB] Packing the object files"
	@ar -rvcs ../lib/lib_bpe.a ${TARGETS}
	@echo " [INCL] Moving includes files"
	@cp -f *.h ../includes/

${DIRS}:
	@echo " [MSG] Creating directories"
	@mkdir -p ${DIRS}

clean:
	@echo " [CLN] Removing object files and binaries"
	@rm -f *~ $(OBJ)/*.o ../bin/test;

