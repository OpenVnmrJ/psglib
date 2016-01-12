/* cpmgt2 -
  carr-purcell meiboom-gill t2 sequence
  
  contact- G. Gray (palo alto)   revision -      */

#include <standard.h>
pulsesequence()
{   
    double n,bt,r;

/* calculate 'big tau' values */
    bt = getval("bt");
    n =  bt/(2.0*d2);
    n = (int)((n/2.0)+0.5)*2.0;
    initval(n,v3);

/* equilibration period */
    status(A);
    hsdelay(d1);
/* calculate exact delay and phases */
    r = d2-p1/2.0-rof2;   /* correct delay for pulse width */
    mod2(oph,v2);   /* 0,1,0,1 */
    incr(v2);   /* 1,2,1,2 = y,y,-y,-y */

/* spin-echo loop */
    status(B);
    rgpulse(pw,oph,rof1,rof2);
    starthardloop(v3);
    	delay(r);
    	rgpulse(p1,v2,rof2,rof2); 
    	delay(r);
    endhardloop();
/* observation period */
    status(C);
} 

