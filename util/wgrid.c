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

/* This subroutine writes a binary PLOT3D grid file.  It works on both 
   on the iris and cray, but it is most useful on the irises */

#ifdef CRAY
void WGRID( double *x, double *y, double *z, int *nx, int *ny, int *nz,
            int *lnx, int *lny, char *name)
#else
void wgrid_( double *x, double *y, double *z, int *nx, int *ny, int *nz,
             int *lnx, int *lny, char *name)
#endif
{
  int i, j, k, l;
  float tmp;
  FILE *fp;

/*  All doubles are cast to floats before writing out  */
      
  fp = fopen(_fcdtocp(name),"w+b");
  
  fwrite( nx, sizeof(int), 1, fp);
  fwrite( ny, sizeof(int), 1, fp);
  fwrite( nz, sizeof(int), 1, fp);
  for (k=0; k<*nz; k++) {
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
	l = i + *lnx * j + *lnx * *lny * k;
	tmp = (float)(*(x+l));
	fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
  }
  for (k=0; k<*nz; k++) {
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
	l = i + *lnx * j + *lnx * *lny * k;
	tmp = (float)(*(y+l));
	fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
  }
  for (k=0; k<*nz; k++) {
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
	l = i + *lnx * j + *lnx * *lny * k;
	tmp = (float)(*(z+l));
	fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
  }
  fclose(fp);
}


#ifdef CRAY
void WGRID2D( double *x, double *y, int *nx, int *ny, int *lnx, char *name)
#else
void wgrid2d_( double *x, double *y, int *nx, int *ny, int *lnx, char *name)
#endif
{
  int i, j, l;
  float tmp;
  FILE *fp;

/*  All doubles are cast to floats before writing out  */
      
  fp = fopen(_fcdtocp(name),"w+b");
  
  fwrite( nx, sizeof(int), 1, fp);
  fwrite( ny, sizeof(int), 1, fp);
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
	l = i + *lnx * j;
	tmp = (float)(*(x+l));
	fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
	l = i + *lnx * j;
	tmp = (float)(*(y+l));
	fwrite( &tmp, sizeof(float), 1, fp );
      }
    }
  fclose(fp);
}

#ifdef CRAY
void WGRIDf( float *x, float *y, int *nx, int *ny, int *nz,
             int *lnx, int *lny, char *name)
#else
void wgridf_( float *x, float *y, int *nx, int *ny, int *nz,
              int *lnx, int *lny, char *name)
#endif
{
  int i, j, k, l;
  float z=0.0;
  FILE *fp;
  
  fp = fopen(_fcdtocp(name),"w+b");
  
  fwrite( nx, sizeof(int), 1, fp);
  fwrite( ny, sizeof(int), 1, fp);
  fwrite( nz, sizeof(int), 1, fp);
  for (k=0; k<*nz; k++) {
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
	l = i + *lnx * j + *lnx * *lny * k;
	fwrite( (x+l), sizeof(float), 1, fp );
      }
    }
  }
  for (k=0; k<*nz; k++) {
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
	l = i + *lnx * j + *lnx * *lny * k;
	fwrite( (y+l), sizeof(float), 1, fp );
      }
    }
  }
  for (k=0; k<*nz; k++) {
    for (j=0; j<*ny; j++) {
      for (i=0; i<*nx; i++) {
	l = i + *lnx * j + *lnx * *lny * k;
	fwrite( &z, sizeof(float), 1, fp );
      }
    }
  }
  fclose(fp);
}
