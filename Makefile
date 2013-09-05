CC=g++
OBJ=obj
DIRS+=$(OBJ)

LIBCDS=libcds
LIBBPE=libbpe
LIB= $(LIBBPE)/lib/lib_bpe.a $(LIBCDS)/lib/libcds.a 

CFLAGS= -ggdb -m32 -fpermissive  -I$(LIBCDS)/includes/ -I$(LIBBPE)/includes/
#CFLAGS= -O9 -DNDEBUG -m32 -fpermissive  -I$(LIBCDS)/includes/

TARGETS = $(OBJ)/RePairWaveletTree.o $(OBJ)/RePairWaveletTreeBuilder.o $(OBJ)/GraphSequence.o $(OBJ)/GraphSequenceTest.o

.PHONY: build check

all: ${DIRS} pre build check test


test: pre $(OBJ)/test.o ${TARGETS}
	$(CC) -o test  $(OBJ)/test.o ${TARGETS} ${LIB} $(CFLAGS)

check: pre $(OBJ)/check.o ${TARGETS}
	$(CC) -o check  $(OBJ)/check.o ${TARGETS} ${LIB} $(CFLAGS)

build: pre $(OBJ)/build.o ${TARGETS}
	$(CC) -o build  $(OBJ)/build.o ${TARGETS} ${LIB} $(CFLAGS)

pre:
	@echo " [BLD] Building BPE"
	@make --no-print-directory -C $(LIBBPE)


$(OBJ)/%.o: %.cpp *.h
	$(CC) -c $(CFLAGS) $< -o $@


${DIRS}:
	@echo " [MSG] Creating directories"
	@mkdir -p ${DIRS}


clean:
	@echo " [CLN] Cleaning objects and binaries"
	@rm -f *~ $(OBJ)/*.o  build check;

cleanall: clean
	@echo " [CLN] Cleaning BPE"
	@make --no-print-directory -C $(LIBBPE) cleanall
	@echo " [CLN] Out of cleaning BPE"
	
