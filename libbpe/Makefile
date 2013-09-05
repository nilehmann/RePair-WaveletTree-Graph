DIRS=includes bin lib
LIBCDS=../libcds


all: ${DIRS} pre libcompact 

${DIRS}:
	@echo " [MSG] Creating directories"
	@mkdir -p ${DIRS}

libcompact:
	@echo " [MSG] Entering directory src"
	@make --no-print-directory -C src 

pre:
	@echo " [BLD] Building libcds"
	@make --no-print-directory -C $(LIBCDS) libcompact

clean:
	@echo " [MSG] Entering directory src"
	@make --no-print-directory -C src clean
	@echo " [MSG] Out directory src"
	@echo " [CLN] Cleaning lib folder"
	@rm -f lib/*
	@echo " [CLN] Cleaning includes folder"
	@rm -f includes/*
	@rm -f *~


cleanall: clean
	@echo " [CLN] Cleaning libcds"
	@make --no-print-directory -C $(LIBCDS) clean
