

gcc pazm.c  psrFits.c psrPalett.c  -o pazm -O3 -fopenmp -DOPENMP -lcfitsio -lcpgplot -lpgplot -lX11 -lgfortran -lm -lgsl -lgslcblas -L/home/zhoudejiang/soft/libpng/lib -lpng -lz -Wall
