# Only works on 64-bit.
# Edit these variables.
MEXSUFFIX   =  mex
# wangyu mexa64
MEX         = mex
CXX         = g++
PARDISOHOME = C:\WorkSpace\MATLAB\bbw_demo\my_toolbox\external\pardiso\pardiso-matlab\libpardiso500-WIN-X86-64
# wangyu $(HOME)/pardiso/lib

# Do not edit anything past here (unless you know what you're doing).
MEXFLAGS   = -cxx CXX=$(CXX) CC=$(CXX) LD=$(CXX) -largeArrayDims \
             -L$(PARDISOHOME) -lpardiso -lmwlapack -lmwblas \
              /libpardiso500-WIN-X86-64/libpardiso500-WIN-X86-64.lib  -lgfortran -lpthread -lm -Wl,--rpath,$(PARDISOHOME)
# /usr/lib64/libgomp.so.1
# \WorkSpace\MATLAB\bbw_demo\my_toolbox\external\pardiso\pardiso-matlab\libpardiso500-WIN-X86-64\libpardiso500-WIN-X86-64.lib
TARGETS_WO_MEXSUFFIX = pardisoinit pardisoreorder \
          pardisofactor pardisosolve \
          pardisofree

TARGETS = pardisoinit.$(MEXSUFFIX) pardisoreorder.$(MEXSUFFIX) \
          pardisofactor.$(MEXSUFFIX) pardisosolve.$(MEXSUFFIX) \
          pardisofree.$(MEXSUFFIX)

all: $(TARGETS_WO_MEXSUFFIX)

%: %.cpp common.cpp matlabmatrix.cpp sparsematrix.cpp \
                pardisoinfo.cpp
	$(MEX) $(MEXFLAGS) -output $@ $^

clean:
	rm -f $(TARGETS) 
