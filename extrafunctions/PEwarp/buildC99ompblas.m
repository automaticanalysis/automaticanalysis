function buildC99ompblas(mexfile)

eval(['mex -v -largeArrayDims CFLAGS="\$CFLAGS -std=c99 -O3 -fopenmp" LDFLAGS="\$LDFLAGS -lgomp -lmwblas" ' mexfile]);

