SHELL := /bin/bash # Use bash syntax
FFLAGS	= -O3 -mcmodel=medium
LDFLAGS = -O3 -mcmodel=medium

.DEFAULT_GOAL := default

include ../../../conf.mk

ifeq ($(TOOLCHAIN), intel)
F_CPL=source /softwares/ifort/intel/composer_xe_2013.5.192/bin/compilervars.sh intel64 && ifort
URQMD_OBJDIR=$(URQMD_SRC)/obj_x86_64.Intel
else
F_CPL=$(FC)
URQMD_OBJDIR=$(URQMD_SRC)/obj_x86_64
endif

URQMD_SRC=../..

default: libcurqmd.a


libcurqmd.a: curqmd_evt.f cmain.f curqmd_out.f cinit.f curqmd_input.f $(URQMD_OBJDIR)/dir
	$(F_CPL) $(FFLAGS) -c cmain.f -o cmain.o -I$(URQMD_SRC)
	$(F_CPL) $(FFLAGS) -c cinit.f -o cinit.o -I$(URQMD_SRC)
	$(F_CPL) $(FFLAGS) -c curqmd_input.f -o curqmd_input.o -I$(URQMD_SRC)
	$(F_CPL) $(FFLAGS) -c curqmd_evt.f -o curqmd_evt.o -I$(URQMD_SRC)
	$(F_CPL) $(FFLAGS) -c curqmd_out.f -o curqmd_out.o -I$(URQMD_SRC)
	objs=(`/bin/ls $(URQMD_OBJDIR)/*.o | grep -v '/urqmd.o'`); \
	ar r libcurqmd.a $${objs[@]} cmain.o cinit.o curqmd_input.o curqmd_evt.o curqmd_out.o

$(URQMD_OBJDIR)/dir:
	 -d $(URQMD_OBJDIR) || { echo "Missing $(URQMD_OBJDIR); build vanilla UrQMD first"; exit 1; }
	 $@

clean:
	rm -f *.a *.o
