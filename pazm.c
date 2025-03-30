#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include "psrFits.h"

// Use OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

void help();
int main(int argc,char **argv)
{
    int i;
    OutputPs=0;
    fs = 1;
    bs = 1;
    begin_phase = 0.0f;
    end_phase = 1.0f;
    colorType = 1;

    subintnum=1;
    on_mark=0;
    off_mark=0;
    no_plt=0;
    plotI=0;
    sub_BIN_N=8;

    impulse_filter_plot=31;

    sigmat=4.0;
    plotIQUV=0;

    zapbin_mark=0;
    alls=0;

    if(argc < 2) help();
    for(i=0;i<argc;i++)
    {
        if (strcmp(argv[i],"-h")==0) help();
        else if(strcmp(argv[i],"-t")==0)
        {
            sscanf(argv[++i],"%d",&ncpus);
        }
        else if (strcmp(argv[i],"-g")==0) OutputPs=1;
        // else if (strcmp(argv[i],"-d")==0) 
        // {
        //     sscanf(argv[++i],"%f,%f",&begin_phase,&end_phase);
        // }
        // else if (strcmp(argv[i],"-f")==0) 
        // {
        //     sscanf(argv[++i],"%d",&fs);
        // }
        // else if (strcmp(argv[i],"-b")==0) 
        // {
        //     sscanf(argv[++i],"%d",&bs);
        // }
        else if(strcmp(argv[i],"-type")==0)
            sscanf(argv[++i],"%d",&colorType);
        else if (strcmp(argv[i],"-s")==0) 
        {
            sscanf(argv[++i],"%d",&subintnum);
        }
        // else if (strcmp(argv[i],"-on")==0) 
        // {
        //     sscanf(argv[++i],"%f,%f",&on_x_1[0], &on_x_2[0]);
        //     on_x_1[1] = on_x_1[0];
        //     on_x_2[1] = on_x_2[0];
        //     on_mark=1;
        // }
        // else if (strcmp(argv[i],"-off")==0) 
        // {
        //     sscanf(argv[++i],"%f,%f",&off_x_1[0], &off_x_2[0]);
        //     off_x_1[1] = off_x_1[0];
        //     off_x_2[1] = off_x_2[0];
        //     off_mark=1;
        // }
        // else if (strcmp(argv[i],"-np")==0) 
        //     no_plt=1;
        else if (strcmp(argv[i],"-alls")==0) 
            alls=1;
        
        else if(strcmp(argv[i],"-ipp")==0)
        {
            sscanf(argv[++i],"%d",&impulse_filter_plot);
        }
        else if(strcmp(argv[i],"-sig")==0)
        {
            sscanf(argv[++i],"%lf",&sigmat);
        }
        else if(strcmp(argv[i],"-zapb")==0)
        {
            sscanf(argv[++i],"%f,%f",&zapbin_s[zapbin_mark], &zapbin_e[zapbin_mark]);
            zapbin_mark++;
        }
        
    }
    strcpy(fitsname, argv[argc-1]);
    printf("Read file: %s\n",fitsname);
    // fitsname=argv[argc-1]; 

    if(zapbin_mark==0 && impulse_filter_plot!=0)
    {
        zapbin_mark=1;
        zapbin_s[0] = 0.0;
        zapbin_e[0] = 1.0;
    }

    if (ncpus > 1) {
#ifdef _OPENMP
        int maxcpus = omp_get_num_procs();
        int openmp_numthreads = (ncpus <= maxcpus) ? ncpus : maxcpus;
        ncpus = openmp_numthreads;
        // Make sure we are not dynamically setting the number of threads
        omp_set_dynamic(0);
        omp_set_num_threads(openmp_numthreads);
        printf("Using %d threads with OpenMP\n\n", ncpus);
#endif
    } else {
#ifdef _OPENMP
        omp_set_num_threads(1); // Explicitly turn off OpenMP
#endif
    }


    ReadFits();


    return 0;
}

void help()
{
    printf("\n");
    printf("A software for Mow the lawn for a .ar file.\n");
    printf("Usage: \n");
    printf("       pazm [...] file.ar\n");
    printf("\n");

    printf("Selection options:\n");
#ifdef _OPENMP
    printf("  -t n           umber of processors to use with OpenMP (default 1)\n");
#endif
    // printf("  -d a,b              Set plot range (default 0,0 -- 1.0)\n");
    // printf("  -f fs               F-scrunch factor just for plot\n");
    // printf("  -b bs               B-scrunch factor just for plot\n");
    // printf("  -on  a,b            Mark on pulse range\n");
    printf("  -s s                Set view subint num (default 1)\n");
    // printf("  -np                 Do not plot\n");
    printf("  -g                  output as a ps file\n");
    printf("\nImpulse filter: \n");
    printf("  -ipp  k        Add plot line of profile after Impulse filter and set the symmetric centered moving window of size 'k' (0 is not do this, surggest 31 for FAST 49.152us samp)\n");
    printf("  -zapb s,e      Set sep bin phase (default 0.0 to 1.0)\n");
    printf("  -sig  s        Set the sigma for Impulse filter (do with -ipp, and default 4.0)\n");
    printf("  -alls          Do Impulse filt for all subint ( must do with -ipp, can with -zapb and -sig)\n");

    printf("\n");
    printf("examples:\n");
    printf("pazm -ipp 31 -zapb 0.0,0.5 -zapb 0.6,1.0 -sig 4.0 -alls  ***.ar     /     pazm -s 2 \n");
    printf("\n");
    printf("  -h                   plot this help\n");

    exit(0);
}