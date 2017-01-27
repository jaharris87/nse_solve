## Build search path
VPATH = .:./src

## Inlcude options for build configurations defined in Makefile.opt
include ./Makefile.opt

## Configure Libraries and files to link appropriate EOS
ifeq ($(EOS),BAHCALL)
  EOS_OBJ = eos_bahcall.o
else ifeq ($(EOS),HELMHOLTZ)
  EOS_OBJ = eos_helm.o helmholtz.o
  VPATH += $(HELMHOLTZ_PATH)
endif

.DEFAULT_GOAL := $(EXE)

## Inlcude system configurations and architecture/installation specific MACROS
include ./Makefile.internal
-include ./Makefile.dev

$(EXE): control.o data.o ffn.o net_preprocess.o nse_module.o $(EOS_OBJ) nse_slice.o
	$(LDR) $(LDRFLAGS) -o xnse \
	    control.o data.o ffn.o net_preprocess.o nse_module.o $(EOS_OBJ) nse_slice.o \
	    $(notdir $(EXTRA_OBJ)) $(EXTRA_LINK) $(LIBS)

#
# Rules for compiling individual files.
#
eos_helm.o: eos_helm.f90
	$(FC) $(FFLAGS) -I$(HELMHOLTZ_PATH) -c $< -o $@
nse_module.o: nse_module.f90
	$(FC) $(FFLAGS) $(LAPACK_INC) -c $< -o $@
common.o: common.f90
	$(FC) $(FFLAGS) $(LAPACK_INC) -c $< -o $@

$(LAPACK_OBJ_F90): %.o: %.f90
	$(FC) $(FFLAGS) $(LAPACK_INC) -c $< -o $(notdir $@)
$(LAPACK_OBJ_F): %.o: %.f
	$(FC) $(FFLAGS) $(LAPACK_INC) -c $< -o $(notdir $@)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f core *.o *.oo *.mod *.lst *.cub *.ptx *.i *.T *.diag xref.db
	rm -rf $(INLINE_DB)
