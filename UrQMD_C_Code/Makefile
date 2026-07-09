#A Makefile using with root system and oracle grid
#'make help' for detail information

SHELL := /bin/bash

ROOTP=YES
NUM=20
CLONE_LIST=filelist
SUBMIT=./submit

LIST_PREV_NUM = 30

PROJ_INCS =

C_SRC=$(wildcard *.cc)
SHELL_SRC=$(wildcard *.shell)
SRC=$(C_SRC) $(SHELL_SRC)
HEAD=$(wildcard *.h)
HEAD+=$(wildcard inc/*.h)
PROG_OBJ=$(patsubst %.cc, _build/%.o, $(C_SRC))
HC_SRC=$(patsubst %.h, %.cpp, $(HEAD))
HC_SRC+=$(wildcard inc/*.cpp)
HC_OBJ=$(patsubst %.h, _build/*.o, $(HEAD))
BC_SRC=$(wildcard *.cpp)
LIBS_SRC=$(foreach cfile, $(HC_SRC), $(findstring $(cfile), $(BC_SRC)))
DLIBS_OBJ=
ifneq (, $(LIBS_SRC))
LIBS_OBJ=$(DLIBS_OBJ) $(patsubst %.cpp, _build/%.o, $(LIBS_SRC))
LIBS_OBJ+=$(patsubst inc/%.cpp, _build/%.o, $(wildcard inc/*.cpp))
endif
C_PROGS=$(patsubst %.cc, _build/%, $(C_SRC))
SHELL_PROGS=$(patsubst %.shell, _build/%, $(SHELL_SRC))
PROGS=$(C_PROGS) $(SHELL_PROGS)
PROGS_SH=$(patsubst %, _build/%.sh, $(PROGS))
PROGS_N=$(patsubst _build/%, %, $(PROGS))
PROGS_RUN=$(patsubst _build/%, r_%, $(PROGS))

RUN=$(word 1, $(PROGS_N))
RUNNAME=
MERGEROOT=
TRY=0

CXX ?= g++

MORECFLAGS=
MORELDFLAGS=
MORELIBS=

.PHONY: all clean run local-run local-clone local-parallel rebuild-local dist-merge merge clone grid old-clone declone cpdata exdata update update_back dist-ext ext runname name submit_script
.DEFAULT_GOAL := all

RUN_DEP =
COMP_DEP =
CLONE_DEP =
EXTCLEAN =
OTHERCLEAN =

# define by conf.mk
ifneq ($(wildcard conf.mk), )
include conf.mk
endif

ifdef ROOTP
EXTCFLAGS=$(MORECFLAGS) $(shell root-config --cflags)
EXTLDFLAGS=$(MORELDFLAGS) $(shell root-config --ldflags)
LIBS=$(MORELIBS) $(shell root-config --libs)
GLIBS=$(shell root-config --glibs)
else
EXTCFLAGS=$(MORECFLAGS)
EXTLDFLAGS=$(MORELDFLAGS)
LIBS=$(MORELIBS)
endif

ifdef DB
CFLAGS  = -g -Wall $(EXTCFLAGS)
LDFLAGS = -g -Wall $(EXTLDFLAGS)
else
CFLAGS  = -O2 $(EXTCFLAGS)
LDFLAGS = -O2 $(EXTLDFLAGS)
endif

all: $(C_PROGS) $(RUN_DEP)
	@echo === Programs placed at _build/ ===

#$(C_PROGS): $(PROG_OBJ) $(LIBS_OBJ)
$(C_PROGS): % : %.o $(LIBS_OBJ)

#$(PROG_OBJ): $(C_SRC) $(HEAD)

ifdef FINE_RULE
$(LIBS_OBJ): $(LIBS_SRC) $(HEAD)
endif

run: $(PROGS_RUN)
	echo "herehrerherhehr"
	@if [ $(TRY) -ne 0 ]; then \
		head -n $(TRY) $(CLONE_LIST) > r_$(RUN)/$(CLONE_LIST); \
	fi
	cd r_$(RUN) && ./$(RUN)

local-run: _build/$(RUN) $(CLONE_DEP)
	@cdir=`pwd`; \
	rdir="r_$(RUN)"; \
	mkdir -p "$$rdir"; \
	ln -sf "$$cdir/_build/$(RUN)" "$$rdir/$(RUN)"; \
	if [[ -e "$(CLONE_LIST)" ]]; then \
		head -n 1 "$(CLONE_LIST)" > "$$rdir/$(CLONE_LIST)"; \
	else \
		echo 12345 > "$$rdir/$(CLONE_LIST)"; \
	fi; \
	if [[ -n "$(PROJ_INCS)" ]]; then \
		for f in $(PROJ_INCS); do \
			ln -sf "$$cdir/$$f" "$$rdir/"; \
		done; \
	fi; \
	echo "Running one local job in $$rdir"; \
	cd "$$rdir" && ./$(RUN)

local-clone: _build/$(RUN) Data $(CLONE_DEP)
	@cdir=`pwd`; \
	rm -rf Data/*; \
	for i in `seq -w $(NUM)`; do \
		rdir="Data/$(RUN)_$${i}_run"; \
		mkdir -p "$$rdir"; \
		ln -sf "$$cdir/_build/$(RUN)" "$$rdir/$(RUN)"; \
		seed_line=$$((10#$$i)); \
		if [[ -e "$(CLONE_LIST)" ]]; then \
			sed -n "$${seed_line}p" "$(CLONE_LIST)" > "$$rdir/$(CLONE_LIST)"; \
		else \
			echo "$${seed_line}" > "$$rdir/$(CLONE_LIST)"; \
		fi; \
		if [[ -n "$(PROJ_INCS)" ]]; then \
			for f in $(PROJ_INCS); do \
				ln -sf "$$cdir/$$f" "$$rdir/"; \
			done; \
		fi; \
	done; \
	echo $(NUM) > Data/total; \
	echo "Created $(NUM) local run directories under Data/"

local-parallel: local-clone
	@echo "Running up to $(JOBS) local UrQMD jobs in parallel"; \
	find Data -maxdepth 1 -type d -name '$(RUN)_*_run' -print0 | \
		xargs -0 -P $(JOBS) -I {} bash -c 'cd "$$1" && ./$(RUN) > run.log 2>&1' _ {}; \
	echo "Finished local jobs. Use: bash pull_data.sh"

rebuild-local:
	$(MAKE) -C urqmd-4.0 distclean
	$(MAKE) -C src/mysrc clean
	$(MAKE)

r_%: _build/% $(COMP_DEP) submit
	@for item in $(PROGS_N); do \
		pwd=`pwd`; \
		if ! [ -d "r_$${item}" ]; then \
			mkdir r_$${item}; \
		fi; \
		ln -sf $${pwd}/_build/$$item $${pwd}/r_$${item}/$${item}; \
		./submit tmpl $${item} "./$${item}" "$${item}.sh"; \
		mv $${item}.sh r_$${item}/; \
		mv subq.con r_$${item}/; \
		if [[ -e $(CLONE_LIST) ]]; then \
			head -n $(LIST_PREV_NUM) $(CLONE_LIST) > r_$${item}/$(CLONE_LIST); \
		fi; \
		if [[ -n "$(PROJ_INCS)" ]]; then \
			for f in $(PROJ_INCS); do \
				ln -sf $${pwd}/$$f $${pwd}/r_$${item}/; \
			done; \
		fi; \
	done; \

_build/%: _build/%.o
	@if ! [ -e _build ]; then \
		mkdir _build; \
	fi
	$(CXX) -o $@ $< $(LIBS_OBJ) $(LDFLAGS) $(LIBS)

_build/%.o: %.cc $(HEAD)
	@if ! [ -e _build ]; then \
		mkdir _build; \
	fi
	$(CXX) -c $(CFLAGS) $< -o $@

_build/%.o: %.cpp
	@if ! [ -e _build ]; then \
		mkdir _build; \
	fi
	$(CXX) -c $(CFLAGS) $< -o $@

_build/%.o: inc/%.cpp
	@if ! [ -e _build ]; then \
		mkdir _build; \
	fi
	$(CXX) -c $(CFLAGS) $< -o $@

_build/%: %.shell
	@if ! [ -e _build ]; then \
		mkdir _build; \
	fi
	cp $< $@; chmod +x $@

update: update_back

update_back:
	rm -rf _build && make

clean: declone $(OTHERCLEAN)
	rm -f *.o
	rm -rf r_*
	rm -rf _build
	rm -rf _merge
	@if [[ -x ./submit ]]; then ./submit clean; fi
	rm -rf $(EXTCLEAN)
	@if [ -e proj.conf ]; then \
		for f in `cat proj.conf`; do \
			rm -rf $$f; \
		done; \
	fi

ifeq (, $(RUNNAME))
submit_script: _build/$(RUN).sh
else
submit_script: _build/$(RUNNAME).sh
endif

_build/$(RUN).sh: submit
	@./submit tmpl $(RUN) './$(RUN)' '$(RUN).sh'; \
	mv $(RUN).sh _build; \
	mv subq.con _build

_build/$(RUNNAME).sh: submit
	@./submit tmpl $(RUNNAME) './$(RUN)' '$(RUNNAME).sh'; \
	mv $(RUNNAME).sh _build; \
	mv subq.con _build

#$(CLONE_LIST):
	#/bin/ls $(DIST)/*.root > $(CLONE_LIST)

Data:
	@mkdir Data

grid: _build/$(RUN) Data  _build/$(RUN).sh $(CLONE_DEP)
	@cdir=`pwd`;\
	echo wait a while; \
	exesh=$(RUN).sh; \
	if [[ -n "$(RUNNAME)" ]]; then \
		./submit tmpl $(RUNNAME) './$(RUN)' '$(RUNNAME).sh'; \
		rm -f _build/$(RUNNAME).sh && mv $(RUNNAME).sh _build; \
		rm -f _build/subq.con && mv subq.con _build; \
		exesh=$(RUNNAME).sh; \
	fi; \
	rdir=''; \
	if [ -e Data/total ]; then \
		echo ' >>> Data/total exists'; \
		printf 'Continue to delete Data first? (*y*/n/ctrl-c): '; \
		read chos; \
		if [[ $$chos == "n" ]]; then \
			echo 'do noting...'; \
			exit 0; \
		fi; \
		echo 'removing and clone again...'; \
	fi; \
	rm -rf Data/*; \
	if ! [ -e $(CLONE_LIST) ]; then \
		do_list=0; \
		finn=$(NUM); \
	else \
		do_list=1; \
		linen=$$(wc -l $(CLONE_LIST)|cut -d' ' -f1); \
		eachf=`python -c 'print(int(round('"$${linen}"'.0/$(NUM))))' `; \
		tdir=__clone_tmp.dir; \
		rm -rf $$tdir; \
		mkdir $$tdir; \
		ln -s $$PWD/$(CLONE_LIST) $$tdir/cf; \
		cd $$tdir; \
		split -a 5 -l $$eachf cf; \
		rm -f cf; \
		cd ..; \
		cfs=(`/bin/ls $$tdir/*`); \
		acn=`/bin/ls $$tdir/|wc -l`; \
		if [ $$acn -lt 2 ]; then \
			echo Too less number to run; \
			exit 0; \
		fi; \
		if [ $$acn -gt $(NUM) ]; then \
			cat $${cfs[((acn-1))]} >> $${cfs[((acn-2))]} ;\
			cfs=($${cfs[@]:0:((acn-1))}) ;\
		fi; \
		finn=$(NUM); \
		if [ $$acn -lt $(NUM) ]; then \
			finn=$$acn; \
		fi; \
	fi; \
	((c=0));\
	for i in `seq -w $$finn`; do \
		rdir=Data/$(RUN)_$${i}_run; \
		if ! [ -d "$${rdir}" ]; then \
			mkdir $${rdir}; \
		fi; \
		ln -fs $${cdir}/_build/$$exesh $${cdir}/$${rdir}/$$exesh; \
		ln -fs $${cdir}/_build/$(RUN) $${cdir}/$${rdir}/$(RUN); \
		ln -fs $${cdir}/_build/subq.con $${cdir}/$${rdir}/subq.con; \
		if [ $$do_list -eq 1 ]; then \
			mv $${cfs[$${c}]} $${rdir}/$(CLONE_LIST); \
		fi; \
		if [[ -n "$(PROJ_INCS)" ]]; then \
			for f in $(PROJ_INCS); do \
                fullname=`dirname "$$f"`; \
				if [ $${fullname} == "." ]; then \
					fullname=""; \
				fi; \
                nodirname=$${f/$$fullname/}; \
                ln -sf $${cdir}/$$f $${rdir}/$$nodirname; \
			done; \
		fi; \
		((c=c+1)); \
	done; \
	if [ $$do_list -eq 1 ]; then \
		rm -rf $$tdir; \
	fi; \
	echo $$finn > Data/total; \

clone: grid

# submit:
#	@ln -s $(SUBMIT) ./submit

## _merge:
## 	@mkdir _merge
## 
## _extract:
## 	@mkdir _extract
## 
## merge.sh: submit
## 	@if [ -z '$(RUNNAME)' ]; then \
## 		./submit tmpl merge.sh "make dist-merge" merge.sh; \
## 	else \
## 		./submit tmpl '$(RUNNAME)' 'make dist-merge' merge.sh; \
## 	fi
## 
## merge: dist-merge
## 
## dist-ext: _extract
## 	find Data -name '*.root' -exec mv {} _extract/ \;
## 
## ext: dist-ext
## 	
## dist-merge: _merge
## 	@cdir=`pwd`; \
## 	roots=(); \
## 	if [ -n "$(MERGEROOT)" ]; then \
## 		echo Merge at $(MERGEROOT); \
## 		merge_root="$(MERGEROOT)"; \
## 	else \
## 		merge_root=$$cdir; \
## 	fi; \
## 	if [ -e $$merge_root/merge.config ]; then \
## 		rootsn=(`cat $$merge_root/merge.config`); \
## 	else \
## 		for i in $$merge_root/Data/*_run; do \
## 			cd $$i; \
## 			rootsn=(`find * -name '*.root'`); \
## 			cd $$cdir; \
## 			rsync -av --include='*/' --exclude='*' $$i/ _merge/; \
## 			break; \
## 		done; \
## 	fi; \
## 	if [ -n "$(LIST)" ]; then \
## 		rm -rf to.merge.list; \
## 	fi; \
## 	for rf in $${rootsn[@]}; do \
## 		roots=(); \
## 		for i in $$merge_root/Data/*_run; do \
## 			if [ -e $$merge_root/error ]; then \
## 				grep -q $$i $$merge_root/error 2>/dev/null; \
## 				if [ $$? -eq 0 ]; then \
## 					continue; \
## 				fi; \
## 			fi; \
## 			roots=($${roots[@]} $${i}/$${rf}); \
## 			continue; \
## 			rts=`/bin/ls $${i}/$${rf}`; \
## 			if [ $$? -ne 0 ]; then \
## 				exit 1; \
## 			fi; \
## 			roots=($${roots[@]} $${rts}); \
## 		done; \
## 		if [ -n "$(LIST)" ]; then \
## 			echo hadd -f _merge/$${rf} $${roots[@]} >> to.merge.list; \
## 			echo hadd -f _merge/$${rf}; \
## 			continue; \
## 		fi; \
## 		hadd -f _merge/$${rf} $${roots[@]}; \
## 	done; \


	
