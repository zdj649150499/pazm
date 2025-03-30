#ifndef PSRFITS_H_
#define PSRFITS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"
#include "psrPalett.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_filter.h"

// #include "fftw3.h"


#ifdef __cplusplus
extern "C" {
#endif
    
    int ncpus;
    int OutputPs;
    int fs;
    int bs;
    float begin_phase;
    float end_phase;
    char fitsname[256];
    char Ffitsname[256];



    int status;
    int nulval;
    int anynul;

    int NCHAN;
    int NBIN;
    int NSUBINT;
    int NCHAN_f;
    int NBIN_b;
    double         TBIN;                    //获得采样时间

    char POL_TYPE[16];
    int NPOL;

    double *FREQArrayback;
    float *Freq_ori_Array;
    float *Data_ori_Array_I;
    float *Data_ori_Array_I_bk;
    float *Data_ori_profile_I;
    float *Phase_ori_Array;
    
    float *WTSArray;
    float *OFFSArray;
    float *SCLArray;
    short *DATAArrayback;

    int subintnum;
    int on_mark, off_mark;
    float on_x_1[2];
    float on_y_1[2];
    float off_x_1[2];
    float off_y_1[2];

    float on_x_2[2];
    float on_y_2[2];
    float off_x_2[2];
    float off_y_2[2];

    float no_plt;

    int plotIQUV;

    int plotI;

    float SIGMA;

    int sub_BIN_N;
    int sub_BIN_off_mark;


    int impulse_filter_plot;
    double sigmat;

    float zapbin_s[128], zapbin_e[128];
    int zapbin_mark;


    int alls;  /* Do impulse detection*/

    void iniFitsParameter();
    void ReadFits();
    void debasechannels_0(float *array,int x, int y);
    float median(float arr[], int n);
    float maxInFloatArray(float *array, long start, long end);
    void GetSNRfromSum(float *array, int N, float *SNR);
    void subBaseOfImageArray(float *a,int nx,int ny);
    void subBaseOfImageArrayIQUV(float *I, float *Q, float *U, float *V, int nx,int ny);
    void Getsigmaofprofile_1(float *array, int n, float *baselineprofile, float *sigmaprofile);
    void Getsigmaofprofile(float *array, int n, float *baselineprofile, float *sigmaprofile);

    void GetProfile(float *inarray, float *profile, float  *profile_x, int x, int y);
    void Ploat_PF_I(float *array, float *profile, float  *profile_x, int x, int y, int plt_xl, int plt_xr);


    // fftwf_plan plan_transpose(int rows, int cols, float *in, float *out);
    void Impulse_Detection_X(float *array, float *wt, int sN, int eN, int X, int Y);
    void Impulse_Detection_Y(float *array, int N, int M);


    void CompressCounts( float *data, short *target_data, int size,float *scale, float *offset);
    void minmax(float *data, int size, float *min, float *max);
    
#ifdef __cplusplus
}
#endif



#endif
