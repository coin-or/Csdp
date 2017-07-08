/*
 * User exit routine for psd.  This version of the routine simply returns
 * 0, so that psd will not stop.
 */

#include "declarations.h"


#ifdef USESIGTERM

#include <unistd.h>
#include <sys/types.h>
#include <signal.h>

int sigterm_signaled=0;

/*
 * This sets sigterm_signaled to 1.  The next time
 * user_exit runs, it will see this and return 1, so CSDP will stop.
 */

void catch_sigterm(signal)
     int signal;
{
  sigterm_signaled=1;
}

int user_exit(n,k,C,a,dobj,pobj,constant_offset,constraints,X,y,Z,params)
     int n;
     int k;
     struct blockmatrix C;
     double *a;
     double dobj;
     double pobj;
     double constant_offset;
     struct constraintmatrix *constraints;
     struct blockmatrix X;
     double *y;
     struct blockmatrix Z;     
     struct paramstruc params;
{
  /*
   * Stop on any of the following signals.  SIGTERM, SIGXCPU, SIGINT, SIGQUIT
   */
  
  signal(SIGTERM,catch_sigterm);
  signal(SIGXCPU,catch_sigterm);
  signal(SIGINT,catch_sigterm);
  signal(SIGQUIT,catch_sigterm);
  
  /*
   * If a signal to terminate has been raised, then quit.
   */
  
  if (sigterm_signaled==1)
    {
      /*
       * This will stop CSDP with a return code of 10.
       */
      return(1);
    }
  else
    {
      /*
       * This will tell CSDP to continue.
       */
      return(0);
    };
  
  /*
   * Any other positive value >= 2 returned will stop CSDP with a success.
   */
}

#else

int user_exit(n,k,C,a,dobj,pobj,constant_offset,constraints,X,y,Z,params)
     int n;
     int k;
     struct blockmatrix C;
     double *a;
     double dobj;
     double pobj;
     double constant_offset;
     struct constraintmatrix *constraints;
     struct blockmatrix X;
     double *y;
     struct blockmatrix Z;     
     struct paramstruc params;
{
  /*
   * This will tell CSDP to continue.
   */
  return(0);

  /*
   * A returned value of 1 stops CSDP with failure (10)
   * Any other positive value >= 2 returned will stop CSDP with a success.
   */

}

#endif
