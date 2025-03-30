# PAZM: A Lawn Mowing Tool for .ar Files

## Description
PAZM is a specialized software tool designed to "mow the lawn" for `.ar` files. These single pulse files are FITS files folded by `dspsr`. The tool helps clean and process these folded pulsar data files.

## Dependencies
PAZM requires the following libraries to be installed:

- `cfitsio` - For FITS file handling
- `gsl` (version 2.7 or higher) - GNU Scientific Library
- `pgplot` - PGPLOT graphics library

## Installation
To compile PAZM, you have two options:

1. **Simple method**: Run the provided shell script:
   ```bash
   ./sh.sh
   ```
2. **Manual compilation**: Use the following gcc command:
   ```bash
   gcc pazm.c psrFits.c psrPalett.c -o pazm -O3 -fopenmp -DOPENMP \
   -lcfitsio -lcpgplot -lpgplot -lX11 -lgfortran -lm -lgsl -lgslcblas \
   -L/home/zhoudejiang/soft/libpng/lib -lpng -lz -Wall
   ```

## Usage
After successful compilation, you can run the program with:

   ```bash
   ./pazm [options] [input_file]
   ```

## Author

Dejiang Zhou @ NAOC