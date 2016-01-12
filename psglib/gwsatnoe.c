/*  gwsatnoe-
               d1...p1x4....pw....at......

    noe difference: dof is cycled alternately between dof and control
     receiver is alternated; dm should be yyn
         G.Gray Palo Alto, 5 apr 94    */

#include <standard.h>
pulsesequence()
{
   double control,cycles,gzlvl0,gt0,p1pwr,pwpwr;
   char p1shape[MAXSTR], pwshape[MAXSTR];
     cycles = getval("cycles");
    control = getval("control");
      p1pwr = getval("p1pwr"); 
      pwpwr = getval("pwpwr");
     gzlvl0 = getval("gzlvl0");
        gt0 = getval("gt0");
     getstr("p1shape",p1shape);
     getstr("pwshape",pwshape);
     initval(cycles,v1);
     mod2(ct,v2);                   /* 01010101.....*/ 
     loadtable("gwsatnoe");
    /* equilibrium period */
    status(A);
       ifzero(v2);
        offset(dof,DODEV);
       elsenz(v2);
        offset(control,DODEV);
       endif(v2);
       hsdelay(d1);
       lk_hold();
    status(B);
       rlpower(p1pwr,TODEV); 
       if (cycles>0.0)
       {
       starthardloop(v1);
       shaped_pulse(p1shape,p1,zero,20.0e-6,rof2);
         zgradpulse(0.6*gzlvl0,gt0);
       shaped_pulse(p1shape,p1,one,20.0e-6,rof2);
         zgradpulse(0.7*gzlvl0,gt0);
       shaped_pulse(p1shape,p1,two,20.0e-6,rof2);
         zgradpulse(0.8*gzlvl0,gt0);
       shaped_pulse(p1shape,p1,three,20.0e-6,rof2);
         zgradpulse(gzlvl0,gt0);
       endhardloop();
       }
        hsdelay(d2);
        rlpower(pwpwr,TODEV);
     status(C);
       shaped_pulse(pwshape,pw,t1,20.0e-6,rof2);
     setreceiver(t2);
     delay(alfa+1/fb);
     acquire(np,1/sw);
     lk_sample();
}
