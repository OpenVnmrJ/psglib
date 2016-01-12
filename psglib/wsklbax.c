/*  wsklbax-
      water suppression using soft-hard-lock pulse sequence
      uses waveform generator for soft pulse
      all parameters and usage the same as sklbax.c except
      that "shape" defines the shapelib entry for the soft pulse.
      see Sklenar and Bax, JMR 75,378(1987)
*/

#include <standard.h>

pulsesequence()
{
   double satdly,satpwr,trim,phaseinc;
    char shape[MAXSTR];
    getstr("shape",shape);
    if (shape[0] == '\0') 
    {
      text_error("no shape entry? ABORT");
      abort(1);
    }
   stepsize(0.5,TODEV);
   phaseinc = getval("phaseinc");
   phaseinc =2.0*phaseinc;
    if (phaseinc < 0.0) phaseinc=720+phaseinc;
   initval(phaseinc,v1);
     trim   = getval("trim");
     satdly = getval("satdly");    
     satpwr = getval("satpwr");
     tpwr = getval("tpwr");

     status(A);
     rlpower(satpwr,TODEV);           /* reduced xmtr power for soft pulse */
     xmtrphase(v1);         /* correction for phase error from power switch */
     hsdelay(d1);
     rcvroff();
     shaped_pulse(shape,satdly*1e-3,t1,4.0e-6,0.0);    /*satdly entered in ms*/
     rlpower(tpwr,TODEV);
     xmtrphase(zero);                      /* return to normal phase control */
     rgpulse(pw,t2,5e-6,0.0);
     rgpulse(trim,t3,0.0,rof2);
     setreceiver(t4);
     status(B);
}

/* sklbax 
t1 = 3 3 3 3 1 1 1 1 
t2 = 0 2 0 2 0 2 0 2
t3 = 1 1 3 3 1 1 3 3
t4 = 0 1 0 1 0 1 0 1
*/
