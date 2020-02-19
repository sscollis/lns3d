#include <stdio.h>

#ifdef CRAY
  #define _cptofcd(c, l)  	((char *)((((long)(c))&0xfc000000ffffffff) | \
				  	 ((((long)(l))&0x7fffff)<<35)))
  #define _fcdtocp(f)     	((char *)(((long)(f))&0xfc000000ffffffff))
  #define _fcdlen(f)      	((unsigned)((((long)(f))>>35)&0x7fffff))
#else
  #define vfork			(fork)
  #define _cptofcd(c, l)  	((char *)c)
  #define _fcdtocp(f)     	((char *)f)
  #define _fcdlen(f)      	(strlen(f))
#endif

/* This subroutine writes a binary PLOT3D data file.  It works both
   on the iris and cray, but it is most useful on the irises */

#ifdef CRAY
void WDATA( double *fsmach, double *alpha, double *re, double *time,
            double *q1, double *q2, double *q3, double *q4, double *q5, 
            int *nx, int *ny, int *nz, int *lnx, int *lny, char* name)
#else
void wdata_( double *fsmach, double *alpha, double *re, double *time,
             double *q1, double *q2, double *q3, double *q4, double *q5, 
             int *nx, int *ny, int *nz, int *lnx, int *lny, char* name)
#endif
{
  int i, j, k, l;
  float tmp;
  FILE *fp;
    
  fp = fopen(_fcdtocp(name),"w+b");
  
  fwrite( nx, sizeof(int), 1, fp);
  fwrite( ny, sizeof(int), 1, fp);
  fwrite( nz, sizeof(int), 1, fp);
  tmp = (float)(*fsmach);
  fwrite( &tmp, sizeof(float), 1, fp);
  tmp = (float)(*alpha);
  fwrite( &tmp, sizeof(float), 1, fp);
  tmp = (float)(*re);
  fwrite( &tmp, sizeof(float), 1, fp);
  tmp = (float)(*time);
  fwrite( &tmp, sizeof(float), 1, fp);

  for (k=0; k<*nz; k++) {
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
	l = i + *lnx * j + *lnx * *lny * k;
	tmp = (float)(*(q1+l));
/* 	printf( "%d %d %d  %e\n", i, j, k, tmp); */
	fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
  }
  for (k=0; k<*nz; k++) {
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
	l = i + *lnx * j + *lnx * *lny * k;
	tmp = (float)(*(q2+l));
/* 	printf( "%d %d %d  %e\n", i, j, k, tmp); */
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

#ifdef CRAY
void WDATA2D( double *fsmach, double *alpha, double *re, double *time,
              double *q1, double *q2, double *q3, double *q4, 
              int *nx, int *ny, int *lnx, char* name)
#else
void wdata2d_( double *fsmach, double *alpha, double *re, double *time,
               double *q1, double *q2, double *q3, double *q4, 
               int *nx, int *ny, int *lnx, char* name)
#endif
{
  int i, j, l;
  float tmp;
  FILE *fp;
    
  fp = fopen(_fcdtocp(name),"w+b");
  
  fwrite( nx, sizeof(int), 1, fp);
  fwrite( ny, sizeof(int), 1, fp);
  tmp = (float)(*fsmach);
  fwrite( &tmp, sizeof(float), 1, fp);
  tmp = (float)(*alpha);
  fwrite( &tmp, sizeof(float), 1, fp);
  tmp = (float)(*re);
  fwrite( &tmp, sizeof(float), 1, fp);
  tmp = (float)(*time);
  fwrite( &tmp, sizeof(float), 1, fp);

    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
	l = i + *lnx * j ;
	tmp = (float)(*(q1+l));
	fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
	l = i + *lnx * j ;
	tmp = (float)(*(q2+l));
	fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
	l = i + *lnx * j ;
	tmp = (float)(*(q3+l));
	fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
	l = i + *lnx * j ;
	tmp = (float)(*(q4+l));
	fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
  fclose(fp);
}

