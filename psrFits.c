#define _FILE_OFFSET_BITS 64
#include "psrFits.h"


void iniFitsParameter()
{
  status=0;
  nulval=-1;
  anynul=0;
}

void ReadFits()
{
    int i,j, k;
    
    int col;
    int NAXIS2history;

    fitsfile *fp;

    iniFitsParameter();

    if(alls==0)
        fits_open_file(&fp,fitsname,READONLY,&status);
    else
    {
        char shell[128];
        sprintf(shell, "chmod +xw %s", fitsname);
        system(shell);
        fits_open_file(&fp,fitsname,READWRITE,&status);
    }
        

    fits_movabs_hdu(fp, 1, NULL, &status);
    fits_movnam_hdu(fp, BINARY_TBL, "HISTORY", 0,         &status);
    fits_read_key  (fp, TINT,       "NAXIS2",   &NAXIS2history, NULL, &status);
    fits_get_colnum(fp, CASESEN,    "TBIN",     &col, &status);
    fits_read_col(fp, TDOUBLE, col,  NAXIS2history, 1, 1, &nulval,&TBIN, &anynul, &status);
    fits_movnam_hdu(fp, BINARY_TBL, "SUBINT",   0,              &status);
    fits_read_key  (fp, TINT,       "NCHAN",    &NCHAN,   NULL, &status);   // 获得频率通道数
    fits_read_key  (fp, TINT,       "NAXIS2",   &NSUBINT, NULL, &status);   // 获得共有多少个观测SUB
    // fits_read_key  (fp, TINT,       "NAXIS1",   &NAXIS1,  NULL, &status);
    // fits_read_key  (fp, TINT,       "NSBLK",    &NSBLK,   NULL, &status);   // 获得每行数据的采样点数
    fits_read_key  (fp, TINT,       "NPOL",     &NPOL,    NULL, &status);   // 获得极化数
    fits_read_key  (fp, TINT,       "NBIN",     &NBIN,    NULL, &status);   //一个脉冲周期内的点数
    fits_read_key  (fp, TSTRING,    "POL_TYPE", POL_TYPE, NULL, &status);


    FREQArrayback = malloc(sizeof(double)*NCHAN);
    Freq_ori_Array = malloc(sizeof(float)*NCHAN);
    WTSArray = malloc(sizeof(float)*NCHAN);
    OFFSArray = malloc(sizeof(float)*NCHAN*NPOL);
    SCLArray = malloc(sizeof(float)*NCHAN*NPOL);
    DATAArrayback = malloc(sizeof(short)*NBIN*NCHAN*NPOL);

    if(alls==0)
        Data_ori_Array_I = malloc(sizeof(float)*NBIN*NCHAN);
    else
        Data_ori_Array_I = malloc(sizeof(float)*NBIN*NCHAN*NPOL);
    // Data_ori_Array_I_bk = malloc(sizeof(float)*NBIN*NCHAN);
    Data_ori_profile_I = malloc(sizeof(float)*NBIN);
    Phase_ori_Array = malloc(sizeof(float)*NBIN);

    if(strcmp(POL_TYPE, "IQUV") == 0 && NPOL==4)
        plotIQUV=1;

    fits_get_colnum(fp, CASESEN,  "DAT_FREQ",  &col,  &status);
    fits_read_col(fp, TDOUBLE, col, subintnum, 1, NCHAN, &nulval,FREQArrayback, &anynul, &status);


    
    // fftwf_plan tplan1, tplan2;

    // if(impulse_filter_plot)
    // {
    //     tplan1 = plan_transpose(NBIN, NCHAN, Data_ori_Array_I_bk, Data_ori_Array_I_bk);  // for fits to ar
    //     tplan2 = plan_transpose(NCHAN, NBIN, Data_ori_Array_I_bk, Data_ori_Array_I_bk);  // for ar to fits
    // }
    

    for(k=1;k<=NSUBINT; k++)
    {
        if(alls==1)
            subintnum=k;

        fits_get_colnum(fp, CASESEN,  "DAT_WTS",  &col,  &status);
        fits_read_col(fp, TFLOAT, col,  subintnum, 1, NCHAN, &nulval,WTSArray, &anynul, &status);
        fits_get_colnum(fp, CASESEN,  "DAT_OFFS", &col, &status);
        fits_read_col(fp, TFLOAT, col,  subintnum, 1, NCHAN*NPOL, &nulval,OFFSArray, &anynul, &status);
        fits_get_colnum(fp, CASESEN,  "DAT_SCL",  &col,  &status);
        fits_read_col(fp, TFLOAT, col,  subintnum, 1, NCHAN*NPOL, &nulval,SCLArray, &anynul, &status);
        fits_get_colnum(fp, CASESEN,  "DATA",     &col,          &status);
        fits_read_col(fp, TSHORT, col,  subintnum, 1, NBIN*NCHAN*NPOL, &nulval,DATAArrayback, &anynul, &status);

        
        if(alls==0)
        {
            if(strcmp(POL_TYPE, "IQUV") == 0 && NPOL==4)
                for(i=0;i<NCHAN;i++)
                {
                    Freq_ori_Array[i]=(float)FREQArrayback[i];
                    for(j=0;j<NBIN;j++)
                    {
                        Data_ori_Array_I[i*NBIN+j]=WTSArray[i]*(DATAArrayback[i*NBIN+j]*SCLArray[i]+OFFSArray[i]);
                    }
                }
            else if(NPOL==2)
                for(i=0;i<NCHAN;i++)
                {
                    Freq_ori_Array[i]=(float)FREQArrayback[i];
                    for(j=0;j<NBIN;j++)
                    {
                        Data_ori_Array_I[i*NBIN+j]=WTSArray[i]*((DATAArrayback[i*NBIN+j]*SCLArray[i]+OFFSArray[i]) + (DATAArrayback[NBIN*NCHAN+i*NBIN+j]*SCLArray[i+NCHAN]+OFFSArray[i+NCHAN]));
                    }
                }
            else
                for(i=0;i<NCHAN;i++)
                {
                    Freq_ori_Array[i]=(float)FREQArrayback[i];
                    for(j=0;j<NBIN;j++)
                    {
                        Data_ori_Array_I[i*NBIN+j]=WTSArray[i]*(DATAArrayback[i*NBIN+j]*SCLArray[i]+OFFSArray[i]);
                    }
                }
            subBaseOfImageArray(Data_ori_Array_I, NBIN, NCHAN);

            // if(impulse_filter_plot)
            // {
            //     for(i=0; i<zapbin_mark; i++)
            //     {
            //         printf("Detect imppulse   %f -- %f\n", zapbin_s[i], zapbin_e[i]);
            //         Impulse_Detection_X(Data_ori_Array_I, WTSArray, NBIN*zapbin_s[i], NBIN*zapbin_e[i], NBIN, NCHAN);
            //     }
                
            // }
            GetProfile(Data_ori_Array_I, Data_ori_profile_I,Phase_ori_Array, NBIN, NCHAN);
            Ploat_PF_I(Data_ori_Array_I, Data_ori_profile_I, Phase_ori_Array, NBIN, NCHAN, begin_phase*NBIN, NBIN*end_phase);
        }
        else
        {
            int a, b;
            float scl, offs;

            for(i=0;i<NPOL; i++)
            {
#ifdef _OPENMP
#pragma omp parallel for num_threads(ncpus)  schedule(dynamic) private(a, b, scl, offs) shared(i, Data_ori_Array_I, DATAArrayback, SCLArray, OFFSArray)
#endif
                for(a=0; a<NCHAN; a++)
                {
                    scl=SCLArray[a+i*NCHAN];
                    offs=OFFSArray[a+i*NCHAN];
                    for(b=0;b<NBIN; b++)
                    {
                        Data_ori_Array_I[b+a*NBIN+NBIN*NCHAN*i] = (DATAArrayback[b+a*NBIN+NBIN*NCHAN*i]*scl+offs);
                    }
                }
            }

            if(ncpus==1)
                for(j=0; j<zapbin_mark; j++)
                {
                    int size = impulse_filter_plot;
                    double t = sigmat;
                    size_t noutlier=0.0;
                    int sb=NBIN*zapbin_s[j];
                    int eb=NBIN*zapbin_e[j];
                    int N = eb-sb;
                    gsl_vector *xx=gsl_vector_alloc(N);
                    gsl_vector *yy=gsl_vector_alloc(N);
                    gsl_vector *xmedian = gsl_vector_alloc(N);
                    gsl_vector *xsigma = gsl_vector_alloc(N);
                    gsl_vector_int *ioutlier = gsl_vector_int_alloc(N);
                    gsl_filter_impulse_workspace *impulse_p = gsl_filter_impulse_alloc(size);
                    
                    for(i=0;i<NPOL; i++)
                    {
                        for(a=0; a<NCHAN; a++)
                        {
                            if(WTSArray[a] == 0.0f)
                                continue;
                            for(b=0;b<N;b++)
                            {
                                gsl_vector_set(xx,b,Data_ori_Array_I[sb+b+a*NBIN+i*NCHAN*NBIN]);
                            }
                            gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_MAD, t, xx, yy, xmedian, xsigma, &noutlier, ioutlier, impulse_p);

                            for(b=0;b<N;b++)
                            {
                                Data_ori_Array_I[sb+b+a*NBIN+NBIN*NCHAN*i] = gsl_vector_get(yy,b);
                            }
                        }
                    }
                    gsl_vector_free(xx);
                    gsl_vector_free(yy);
                    gsl_vector_free(xmedian);
                    gsl_vector_free(xsigma);
                    gsl_vector_int_free(ioutlier);
                    gsl_filter_impulse_free(impulse_p);
                }
            else
            {
#ifdef _OPENMP
#pragma omp parallel for num_threads(ncpus)  schedule(dynamic) private(a, b, i, j) shared(NCHAN, NBIN, NPOL, impulse_filter_plot,sigmat, WTSArray, zapbin_s, zapbin_e, Data_ori_Array_I)
#endif
                for(a=0; a<NCHAN; a++)
                {
                    if(WTSArray[a] == 0.0f)
                        continue;
                    for(j=0; j<zapbin_mark; j++)
                    {
                        int size = impulse_filter_plot;
                        double t = sigmat;
                        size_t noutlier=0.0;
                        int sb=NBIN*zapbin_s[j];
                        int eb=NBIN*zapbin_e[j];
                        int N = eb-sb;
                        gsl_vector *xx=gsl_vector_alloc(N);
                        gsl_vector *yy=gsl_vector_alloc(N);
                        gsl_vector *xmedian = gsl_vector_alloc(N);
                        gsl_vector *xsigma = gsl_vector_alloc(N);
                        gsl_vector_int *ioutlier = gsl_vector_int_alloc(N);
                        gsl_filter_impulse_workspace *impulse_p = gsl_filter_impulse_alloc(size);
                        
                        for(i=0;i<NPOL; i++)
                        {
                                for(b=0;b<N;b++)
                                {
                                    gsl_vector_set(xx,b,Data_ori_Array_I[sb+b+a*NBIN+i*NCHAN*NBIN]);
                                }
                                gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_MAD, t, xx, yy, xmedian, xsigma, &noutlier, ioutlier, impulse_p);

                                for(b=0;b<N;b++)
                                {
                                    Data_ori_Array_I[sb+b+a*NBIN+NBIN*NCHAN*i] = gsl_vector_get(yy,b);
                                }
                        }
                        gsl_vector_free(xx);
                        gsl_vector_free(yy);
                        gsl_vector_free(xmedian);
                        gsl_vector_free(xsigma);
                        gsl_vector_int_free(ioutlier);
                        gsl_filter_impulse_free(impulse_p);
                    }
                }
            }


            for(i=0; i<NPOL; i++)
            {
#ifdef _OPENMP
#pragma omp parallel for num_threads(ncpus)  schedule(dynamic) private(j) shared(i, NCHAN, NBIN, WTSArray, Data_ori_Array_I, DATAArrayback, SCLArray)
#endif
                for(j=0;j<NCHAN;j++)
                {
                    if(WTSArray[j] == 0.0f)
                        continue;
                    CompressCounts(Data_ori_Array_I+j*NBIN+i*NBIN*NCHAN, DATAArrayback+j*NBIN+i*NBIN*NCHAN, NBIN, &SCLArray[j+i*NCHAN], &OFFSArray[j+i*NCHAN]);
                }
            }
            fits_get_colnum(fp, CASESEN,  "DAT_OFFS",    &col,      &status); 
            fits_write_col (fp, TFLOAT,col, subintnum, 1, NCHAN*NPOL,  OFFSArray, &status);
            fits_get_colnum(fp, CASESEN,  "DAT_SCL",    &col,      &status); 
            fits_write_col (fp, TFLOAT,col, subintnum, 1, NCHAN*NPOL,  SCLArray, &status);
            fits_get_colnum(fp, CASESEN,  "DATA",    &col,      &status); 
            fits_write_col (fp, TSHORT,col, subintnum, 1, NBIN*NCHAN*NPOL,  DATAArrayback, &status);
        }
        if(alls==0)
            break;
    }
    fits_close_file(fp,&status);

    free(FREQArrayback);
    free(Freq_ori_Array);
    free(WTSArray);
    free(OFFSArray);
    free(SCLArray);
    free(DATAArrayback);
    free(Data_ori_Array_I);
    free(Data_ori_profile_I);

}


void GetProfile(float *inarray, float *profile, float  *profile_x, int x, int y)
{
    int i, j;

    for(i=0; i<x; i++)
    {
        profile[i] = 0.0f;
    }

    for(j=0; j<y; j++)
    {
        for(i=0; i<x; i++)
        {
            profile[i] += inarray[i+j*x];
        }
    }

    for(i=0; i<x; i++)
    {
        profile[i] /= y;
        profile_x[i] = 1.0f*i/x;
    }
}


void Ploat_PF_I(float *array, float *profile, float  *profile_x, int x, int y, int plt_xl, int plt_xr)
{
    int i,j;
    float min, max,med=0.0;
    float min_profile, max_profile;
    float data_bk;
    float sig=0.0;
    float *array_bk;
    float  bright=0.5;
    float  contra=1.0;
    float  Ti[6];


    float *profile_bk;
    
    profile_bk = profile;
    // profile_bk=array+6*x;

    Ti[0]=0.0;
    Ti[1]=1.0;
    Ti[2]=0.0;
    Ti[3]=0.0;
    Ti[4]=0.0;
    Ti[5]=1.0;

    array_bk = malloc(sizeof(float)*(plt_xr-plt_xl)*y);

    min = max = array[plt_xl];
    for(j=0; j<y; j++)
        for(i=plt_xl; i<plt_xr; i++)
        {
            data_bk = array[i+ j*x];
            array_bk[i-plt_xl + j*(plt_xr-plt_xl)] = data_bk;
            min = (min > data_bk) ? data_bk : min;
            max = (max < data_bk) ? data_bk : max;
        }
    
    min_profile = max_profile = profile_bk[plt_xl];
    for(i=plt_xl; i<plt_xr; i++)
    {
        min_profile = (min_profile > profile_bk[i]) ? profile_bk[i] : min_profile;
        max_profile = (max_profile < profile_bk[i]) ? profile_bk[i] : max_profile;
    }


    if(colorType==1) 
    {
        sig=sqrt(gsl_stats_float_variance(array_bk,1,(plt_xr-plt_xl)*y));
        med = median(array_bk, (plt_xr-plt_xl)*y);
    }

    free(array_bk);


    if(OutputPs==0)
        cpgbeg(0,"/xs",1,1);
    else
    {
        char outps[1000];
        sprintf(outps,"%s_%05d",fitsname,subintnum);
        strcat(outps,".png");
        printf("Output png file: %s\n",outps);
        strcat(outps,"/png");
        cpgbeg(0,outps,1,1);
        cpgscr(0,1,1,1);
        cpgslw(2);
    }
    cpgscf(2);
    cpgpap(12, 0.8);
    cpgsvp(0.08,0.98,0.25,0.9);
    palett(colorType,contra,bright);
    cpgswin(plt_xl-0.5,plt_xr+0.5,0.5,y-0.5);


    if(colorType!=1)
        cpgimag(array,x,y,1,x,1,y,min,max,Ti);
    else
        cpgimag(array,x,y,1,x,1,y,med-sig*0.5,med+sig*2.0,Ti);

    cpgswin(profile_x[plt_xl], profile_x[plt_xr-1],Freq_ori_Array[0],Freq_ori_Array[NCHAN-1]);
    if(OutputPs==1)
    {
        cpgsci(14);
    }
    cpgbox("bcist",0,0,"bcnist",0,0);
    cpglab("","Frequency (MHz)","");

    cpgsvp(0.08,0.98,0.1,0.25);
    cpgswin(profile_x[plt_xl], profile_x[plt_xr-1], min_profile-0.05*(max_profile-min_profile), max_profile+0.05*(max_profile-min_profile));
    cpgline(x, profile_x,profile_bk);
    // cpgline(x, profile_x,array+5*x);
    cpgbox("bcnst",0,0,"bct",0,0);
    cpglab("Phase","","");

    if(impulse_filter_plot)
    {

        printf("Plot impulse filted line: \n");
        int k;
        int size = impulse_filter_plot;
        int N=x;
        double t = sigmat;
        size_t noutlier=0.0;


        float *baseline = malloc(sizeof(float)*N);
        gsl_vector *xx=gsl_vector_alloc(N);
        gsl_vector *yy=gsl_vector_alloc(N);
        gsl_vector *xmedian = gsl_vector_alloc(N);
        gsl_vector *xsigma = gsl_vector_alloc(N);
        gsl_vector_int *ioutlier = gsl_vector_int_alloc(N);

        gsl_filter_impulse_workspace *impulse_p = gsl_filter_impulse_alloc(size);
        for(k=0;k<N;k++)
        {
            gsl_vector_set(xx,k,profile_bk[k]);
        }
        // gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_QN, t, xx, yy, xmedian, xsigma, &noutlier, ioutlier, impulse_p);
        // gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_MAD, t, xx, yy, xmedian, xsigma, &noutlier, ioutlier, impulse_p);
        gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_IQR, t, xx, yy, xmedian, xsigma, &noutlier, ioutlier, impulse_p);
        // gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_SN, t, xx, yy, xmedian, xsigma, &noutlier, ioutlier, impulse_p);
        for(k=0;k<N;k++)
        {
            baseline[k] = gsl_vector_get(yy,k);
            cpgsci(2);
            if(baseline[k] != profile_bk[k])
                cpgpt(1, profile_x+k, profile_bk+k, 20);
            // baseline[k] = gsl_vector_get(xmedian, k);

            
        }

        gsl_vector_free(xx);
        gsl_vector_free(yy);
        gsl_vector_free(xmedian);
        gsl_vector_free(xsigma);
        gsl_vector_int_free(ioutlier);
        gsl_filter_impulse_free(impulse_p);
        cpgsci(3);
        cpgline(N,profile_x,baseline);
        cpgsci(1);

        free(baseline);
    }
    
    
    cpgsvp(0.08,0.98,0.9,0.999);
    cpgswin(0,1,0,1);

    char txt[1024];
    float SNR;
    GetSNRfromSum(profile_bk, x, &SNR);
    cpgsci(11);
    sprintf(txt,"%s",fitsname);
    cpgptxt(0,0.8,0,0,txt);
    sprintf(txt,"ori-NBIN: %-9d ori-NCHAN: %-8d ori-NSUB: %-8d",NBIN, NCHAN, NSUBINT);
    cpgptxt(0,0.52,0,0,txt);
    sprintf(txt,"plot-NBIN: %-8d plot-NCHAN: %-7d plot-SUB: %-7d S/N: %.1f",x, y,subintnum, SNR);
    cpgptxt(0,0.25,0,0,txt);
    

    cpgend();
}




void Impulse_Detection_X(float *array, float *wt, int sN, int eN, int X, int Y)
{
    int i,k;
    int size = impulse_filter_plot;
    double t = sigmat;
    size_t noutlier=0.0;
    int N = eN-sN+1;



    gsl_vector *xx=gsl_vector_alloc(N);
    gsl_vector *yy=gsl_vector_alloc(N);
    gsl_vector *xmedian = gsl_vector_alloc(N);
    gsl_vector *xsigma = gsl_vector_alloc(N);
    gsl_vector_int *ioutlier = gsl_vector_int_alloc(N);
    gsl_filter_impulse_workspace *impulse_p = gsl_filter_impulse_alloc(size);

    for(i=0; i<Y; i++)
    {
        if(wt[i] == 0.0f)
            continue;
        for(k=0;k<N;k++)
        {
            gsl_vector_set(xx,k,array[sN+k+i*X]);
        }
        // gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_QN, t, xx, yy, xmedian, xsigma, &noutlier, ioutlier, impulse_p);
        // gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_MAD, t, xx, yy, xmedian, xsigma, &noutlier, ioutlier, impulse_p);
        gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_IQR, t, xx, yy, xmedian, xsigma, &noutlier, ioutlier, impulse_p);
        // gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_SN, t, xx, yy, xmedian, xsigma, &noutlier, ioutlier, impulse_p);
        for(k=0;k<N;k++)
        {
            array[sN+k+i*X] = gsl_vector_get(yy,k);
        }
    }

    gsl_vector_free(xx);
    gsl_vector_free(yy);
    gsl_vector_free(xmedian);
    gsl_vector_free(xsigma);
    gsl_vector_int_free(ioutlier);
    gsl_filter_impulse_free(impulse_p);
}

void Impulse_Detection_Y(float *array, int N, int M)
{
    int i,k;
    int size = impulse_filter_plot;
    double t = sigmat;
    size_t noutlier=0.0;


    gsl_vector *xx=gsl_vector_alloc(N);
    gsl_vector *yy=gsl_vector_alloc(N);
    gsl_vector *xmedian = gsl_vector_alloc(N);
    gsl_vector *xsigma = gsl_vector_alloc(N);
    gsl_vector_int *ioutlier = gsl_vector_int_alloc(N);
    gsl_filter_impulse_workspace *impulse_p = gsl_filter_impulse_alloc(size);

    for(i=0; i<M; i++)
    {
        for(k=0;k<N;k++)
        {
            gsl_vector_set(xx,k,array[k+i*N]);
        }
        // gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_QN, t, xx, yy, xmedian, xsigma, &noutlier, ioutlier, impulse_p);
        // gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_MAD, t, xx, yy, xmedian, xsigma, &noutlier, ioutlier, impulse_p);
        gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_IQR, t, xx, yy, xmedian, xsigma, &noutlier, ioutlier, impulse_p);
        // gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_SN, t, xx, yy, xmedian, xsigma, &noutlier, ioutlier, impulse_p);
        for(k=0;k<N;k++)
        {
            array[k+i*N] = gsl_vector_get(yy,k);
        }
    }

    gsl_vector_free(xx);
    gsl_vector_free(yy);
    gsl_vector_free(xmedian);
    gsl_vector_free(xsigma);
    gsl_vector_int_free(ioutlier);
    gsl_filter_impulse_free(impulse_p);
}



float maxInFloatArray(float *array, long start, long end)
{
  int i;
  float m;
  m=array[start];
  for(i=start+1;i<=end;i++)
    if(m<array[i])
      m=array[i];
  return m;
}


void GetSNRfromSum(float *array, int N, float *SNR)
{
    float meanprofile,sigmaprofile;
    int i,num=0;
    float SUM=0;

    //Debaseline(array,N);
    // Getsigmaofprofile(array,N,&meanprofile,&sigmaprofile);
    Getsigmaofprofile_1(array,N,&meanprofile,&sigmaprofile);
    for(i=0;i<N;i++)
    {
        if(array[i] - meanprofile >3*sigmaprofile)
        {
            num++;
            SUM+=(array[i]-meanprofile);
        }
    }
    
    //*SNR=SUM*sqrt(num)/sigmaprofile;
    *SNR=SUM/(sigmaprofile*sqrt(num));

    if(num==0)
    {
        float max;
        max=maxInFloatArray(array,0,N-1);
        *SNR=(max-meanprofile)/sigmaprofile;
    }
    SIGMA = sigmaprofile;
}

void debasechannels_0(float *array,int x, int y)
{
  int i,j;
  float med;//rms;
  float *backuparray;

  backuparray=malloc(sizeof(float)*x);
  for(i=0;i<y;i++)
  {
    for(j=0;j<x;j++)
    {
      backuparray[j]=array[j+i*x];
    }
    med=median(backuparray,x);
    for(j=0;j<x;j++)
    {
      array[j+i*x]-=med;
    }
  }
  free(backuparray);
}

#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

float median(float arr[], int n)
{
    int low, high;
    int median;
    int middle, ll, hh;

    low = 0;
    high = n - 1;
    median = (low + high) / 2;
    for (;;) {
        if (high <= low)        /* One element only */
            return arr[median];

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]);
            return arr[median];
        }

        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if (arr[middle] > arr[high])
            ELEM_SWAP(arr[middle], arr[high]);
        if (arr[low] > arr[high])
            ELEM_SWAP(arr[low], arr[high]);
        if (arr[middle] > arr[low])
            ELEM_SWAP(arr[middle], arr[low]);

        /* Swap low item (now in position middle) into position (low+1) */
        ELEM_SWAP(arr[middle], arr[low + 1]);

        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        for (;;) {
            do
                ll++;
            while (arr[low] > arr[ll]);
            do
                hh--;
            while (arr[hh] > arr[low]);

            if (hh < ll)
                break;

            ELEM_SWAP(arr[ll], arr[hh]);
        }

        /* Swap middle item (in position low) back into correct position */
        ELEM_SWAP(arr[low], arr[hh]);

        /* Re-set active partition */
        if (hh <= median)
            low = ll;
        if (hh >= median)
            high = hh - 1;
    }
}

#undef ELEM_SWAP


void subBaseOfImageArray(float *a,int nx,int ny)
{
  int i,j;
  float axave,axmax;

//   ax=creatOneDimFloatArray(0,nx-1);

  for(i=0;i<ny;i++)
  {
    // for(j=0;j<nx;j++)
    //   ax[j]=a[i*nx+j];
    // Getsigmaofprofile(a+i*nx,nx,&axave,&axmax);
    Getsigmaofprofile_1(a+i*nx,nx,&axave,&axmax);
    for(j=0;j<nx;j++)
      a[i*nx+j]-=axave;
  }
//   releaseOneDimFloatArray(ax,0);
}

void subBaseOfImageArrayIQUV(float *I, float *Q, float *U, float *V, int nx,int ny)
{
  int i,j;
  float axave[4],axmax[4];

  for(i=0;i<ny;i++)
  {
    // Getsigmaofprofile(I+i*nx,nx,&axave[0],&axmax[0]);
    // Getsigmaofprofile(Q+i*nx,nx,&axave[1],&axmax[1]);
    // Getsigmaofprofile(U+i*nx,nx,&axave[2],&axmax[2]);
    // Getsigmaofprofile(V+i*nx,nx,&axave[3],&axmax[3]);

    Getsigmaofprofile_1(I+i*nx,nx,&axave[0],&axmax[0]);
    Getsigmaofprofile_1(Q+i*nx,nx,&axave[1],&axmax[1]);
    Getsigmaofprofile_1(U+i*nx,nx,&axave[2],&axmax[2]);
    Getsigmaofprofile_1(V+i*nx,nx,&axave[3],&axmax[3]);
    for(j=0;j<nx;j++)
    {
        I[i*nx+j]-=axave[0];
        Q[i*nx+j]-=axave[1];
        U[i*nx+j]-=axave[2];
        V[i*nx+j]-=axave[3];
    }
  }
}

void Getsigmaofprofile_1(float *array, int n, float *baselineprofile, float *sigmaprofile)
{
	int i,j;
	int N=sub_BIN_N;
	float mean[N];
	float sigma[N];

	for(j=0;j<N;j++)
	{
		mean[j]=0;
		for(i=j*n/N;i<(j+1)*n/N;i++)
		{
			mean[j]+=array[i];
		}
		mean[j]/=(n/N);
		sigma[j]=0;
		for(i=j*n/N;i<(j+1)*n/N;i++)
		{
			sigma[j]+=pow(array[i]-mean[j],2);
		}
		sigma[j]=sqrt(sigma[j]/(1.0*n/N));
	}
	float min_mean,min_sigma;
	min_mean=mean[0];
	min_sigma=sigma[0];
    sub_BIN_off_mark=0;
	for(i=1;i<N;i++)
	{
		if(min_mean>mean[i])
		{
			min_mean=mean[i];
			min_sigma=sigma[i];
            sub_BIN_off_mark=1;
		}
	}

	*baselineprofile=min_mean;
	*sigmaprofile=min_sigma;
}


void Getsigmaofprofile(float *array, int n, float *baselineprofile, float *sigmaprofile)
{
  int i;
  float mean=0;
  float sigma=0;
  float sigma1=0;
  float mean1=mean;
  int N;

  for(i=0;i<n;i++)
  {
    mean+=array[i];
  }
  mean/=n;

  for(i=0;i<n;i++)
  {
    sigma+=pow(array[i]-mean,2);
  }
  sigma=sqrt(sigma/n);

  while (fabs(sigma1-sigma)>(0.05*sigma))
  {
    sigma1=sigma;
    mean1=mean;
    N=0;
    mean=0;
    sigma=0;

    for(i=0;i<n;i++)
    {
      if(array[i]-mean1<(3*sigma1))
      {
        N++;
        mean+=array[i];
      }
    }
    mean/=N;

    for(i=0;i<n;i++)
    {
      if(array[i]-mean1<(3*sigma1))
      {
        sigma+=pow(array[i]-mean,2);
      }
    }
    sigma=sqrt(sigma/N);

    if(sigma==0)
    {
      mean=mean1;
      sigma=sigma1;
      break;
    }
  }

  *baselineprofile=mean;
  *sigmaprofile=sigma;
}




// fftwf_plan plan_transpose(int rows, int cols, float *in, float *out)
// {
// // FFTW can be tricked into doing *very* fast transposes
// // (initial testing showed 6-7x faster than TOMs!)
// // http://agentzlerich.blogspot.com/2010/01/using-fftw-for-in-place-matrix.html
// // http://www.fftw.org/faq/section3.html#transpose
//     const unsigned flags = FFTW_MEASURE;        /* other flags are possible */
//     fftwf_iodim howmany_dims[2];
//     howmany_dims[0].n = rows;
//     howmany_dims[0].is = cols;
//     howmany_dims[0].os = 1;
//     howmany_dims[1].n = cols;
//     howmany_dims[1].is = 1;
//     howmany_dims[1].os = rows;
//     return fftwf_plan_guru_r2r( /*rank= */ 0, /*dims= */ NULL,
//                                /*howmany_rank= */ 2, howmany_dims,
//                                in, out, /*kind= */ NULL, flags);
// }


void CompressCounts( float *data, short *target_data, int size,float *scale, float *offset)
{
    // get the range of vlues
    float  min, max;
    float scale_bk, offset_bk;

    min = max = data[0];

    minmax( data, size, &min, &max);

    // determine the scale and offset
    float range = max - min;
    if (max==0 && min==0)
    {
        scale_bk = 1.0;
        offset_bk = 0.0;
    }
    else {
        if (range==0)
        scale_bk = max / 65534.0f;
        else
        scale_bk = range / 65534.0f;
        offset_bk = min - ( -32768.0f * scale_bk );
    }
    
    int d;
    for( d = 0; d < size; d++ )
    {
        target_data[d] = (int) floor( ( data[d] - offset_bk ) / scale_bk + 0.5f);
    }

    *scale = scale_bk;
    *offset = offset_bk;
}

void minmax(float *data, int size, float *min, float *max)
{
    int i;
    float min_bk, max_bk;
    min_bk = max_bk = data[0];
    for(i=1; i<size; i++)
    {
        min_bk = (min_bk>data[i]) ? data[i]:min_bk;
        max_bk = (max_bk<data[i]) ? data[i]:max_bk;
    }
    *min = min_bk;
    *max = max_bk;
}


