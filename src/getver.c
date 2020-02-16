#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <sys/wait.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

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
  
int archpid = 0;

#ifdef CRAY
void ADDVER( char *base, char *suf, int *iver, char *filename, int *length )
#else
void addver_( char *base, char *suf, int *iver, char *filename, int *length )
#endif
{
  if ( sprintf(filename,"%s.%s.%d",_fcdtocp(base),_fcdtocp(suf),*iver) > 
       *length-1 ) {
    printf("Error in addver...\n");
    printf("  Needed string length greater that %d.\n",*length);
    printf("  Increase string length in calling routine\n\n");
    exit(1);
  }
  filename = _cptofcd(filename,*length);
}

#ifdef CRAY
int IGETVER( char *base, char *suf )
#else
int igetver_( char *base, char *suf )
#endif
{
  FILE *istr;
  char command[256], line[256];
  int iver=0;
  
  sprintf(command,"verlist %s.%s 0 0\n", _fcdtocp(base), _fcdtocp(suf) );
  istr = popen(command,"r");
  fgets(line,256,istr);
  pclose(istr);
  sscanf(line,"%d",&iver);
  if ( strlen(line) < 1 ) iver = 0;
  return iver;
}
