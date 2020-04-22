SHELL := /bin/bash
# possible values for f90comp are 'ifort' or 'gfortran'
f90comp = gfortran

# only change to yes when developing
debug = no

# different switches for different compilers
fastintel = -mkl -ipo -O3 -static -fp-model fast=2 -xAVX -axCORE-AVX2
debugintel = -mkl -static -C -debug all -check all -g
fastgnu = -O2 -lblas -llapack #-static-libgfortran
debuggnu =  -Wall -Wrealloc-lhs -fcheck=all -lblas -llapack 

ifeq ($(f90comp),ifort)
ifeq ($(debug),yes)
switch=$(debugintel)
else
switch=$(fastintel)
endif
else
ifeq ($(debug),yes)
switch=$(debuggnu)
else
switch=$(fastgnu)
endif
endif

target := rrreml.exe

objects = getZGZMats.o dsptrf_Ldet.o countNumberLines.o askFilename.o getEffects.o reml_module.o RRReml.o

$(target): $(objects)
	$(f90comp) $(switch) -o $(target) $(objects)
	@mkdir -p obj
	@mv *.o obj/
	@mv *.mod obj/

reml_module.mod: reml_module.o reml_module.f90
	$(f90comp) $(switch) -c reml_module.f90
reml_module.o: reml_module.f90
	$(f90comp) $(switch) -c reml_module.f90
RRReml.o: dsptrf_Ldet.o countNumberLines.o askFilename.o reml_module.mod RRReml.f90
	$(f90comp) $(switch) -c RRReml.f90
%.o: %.f90
	$(f90comp) -c $(switch) $<
clean:
	@command rm -vf obj/*.o obj/*.mod $(target) testdata/*.*.*

test:
	@cd testdata; command rm -f {ranEff,var,fixEff}.*.$(f90comp) ; echo -n "running the case with correlation..."; ../$(target) < input1 > /dev/null  && echo "test correlation successful" ; mv fixedEffects fixEff.cor.$(f90comp); mv randomEffects ranEff.cor.$(f90comp); mv variances var.cor.$(f90comp); echo -n "running the case without correlation..." ; ../$(target) < input2 > /dev/null && echo "test zero correlation successful"; mv randomEffects ranEff.uncor.$(f90comp) ; mv fixedEffects fixEff.uncor.$(f90comp); mv variances var.uncor.$(f90comp); cd ..
