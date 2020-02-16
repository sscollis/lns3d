#include <stdio.h>

#ifdef CRAY
  #define _cptofcd(c, l)        ((char *)((((long)(c))&0xfc000000ffffffff) | \
                                         ((((long)(l))&0x7fffff)<<35)))
  #define _fcdtocp(f)           ((char *)(((long)(f))&0xfc000000ffffffff))
  #define _fcdlen(f)            ((unsigned)((((long)(f))>>35)&0x7fffff))
#else
  #define vfork                 (fork)
  #define _cptofcd(c, l)        ((char *)c)
  #define _fcdtocp(f)           ((char *)f)
  #define _fcdlen(f)            (strlen(f))
#endif

/* This subroutine writes a binary PLOT3D data file.  It works both
   on the iris and cray, but it is most useful on the irises */

#ifdef CRAY
void WDATA( double *q1, double *q2, double *q3, double *q4, double *q5, 
            int *nx, int *ny, int *nz, int *lnx, int *lny, char* name)
#else
void wdata_( double *q1, double *q2, double *q3, double *q4, double *q5, 
             int *nx, int *ny, int *nz, int *lnx, int *lny, char* name)
#endif
{
  int i, j, k, l;
  float fsmach=0.0, alpha=0.0, re=0.0, time=0.0;
  float tmp;
  FILE *fp;
    
  fp = fopen(_fcdtocp(name),"w+b");
  
  fwrite( nx,      sizeof(int), 1, fp);
  fwrite( ny,      sizeof(int), 1, fp);
  fwrite( nz,      sizeof(int), 1, fp);
  fwrite( &fsmach, sizeof(float), 1, fp);
  fwrite( &alpha,  sizeof(float), 1, fp);
  fwrite( &re,     sizeof(float), 1, fp);
  fwrite( &time,   sizeof(float), 1, fp);

  for (k=0; k<*nz; k++) {
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
        l = i + *lnx * j + *lnx * *lny * k;
        tmp = (float)(*(q1+l));
        fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
  }
  for (k=0; k<*nz; k++) {
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
        l = i + *lnx * j + *lnx * *lny * k;
        tmp = (float)(*(q2+l));
        fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
  }
  for (k=0; k<*nz; k++) {
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
        l = i + *lnx * j + *lnx * *lny * k;
        tmp = (float)(*(q3+l));
        fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
  }
  for (k=0; k<*nz; k++) {
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
        l = i + *lnx * j + *lnx * *lny * k;
        tmp = (float)(*(q4+l));
        fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
  }
  for (k=0; k<*nz; k++) {
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
        l = i + *lnx * j + *lnx * *lny * k;
        tmp = (float)(*(q5+l));
        fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
  }
  fclose(fp);
}
