
/* hmqcnoesyhmqc3d3rf - X-filtered 2D cross relaxation experiment 
               USES THIRD CHANNEL FOR X
               USES SS "read" pulse in final hmqc step

   Based on Frienkel, et.al, JMR, 90,420(1990)


   coded by J. Lee (Harvard Med. Sch.) ; revised G. Gray (varian)
 */

#include <standard.h>

static int ph1[1] = {0},
           ph2[2] = {0,2},
           ph3[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
           ph4[8] = {0,0,0,0,2,2,2,2},
           ph5[4] = {0,0,2,2},
           ph6[16] = {0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2},
           ph7[8] = {0,2,2,0,2,0,0,2},
           ph9[4] = {1,3,3,1};

pulsesequence()
{
   double phase,phase2,gzlvl0,gt0,gzlvl1,gt1,j,pwx2,pwx2lvl,
          SSpwr,SSpw,satfrq,satdly,satpwr,mix,tau; 
   int iphase,iphase2;
   char SSshape[MAXSTR],sspul[MAXSTR],satmode[MAXSTR];

/* LOAD VARIABLES */
   satpwr=getval("satpwr"); satfrq=getval("satfrq");
   satdly=getval("satdly"); tau=getval("tau");
   gzlvl0=getval("gzlvl0"); gt0=getval("gt0");
   j=getval("j"); pwx2=getval("pwx2"); pwx2lvl=getval("pwx2lvl");
   SSpw=getval("SSpw"); SSpwr=getval("SSpwr");
   phase=getval("phase"); phase2=getval("phase2");
   gzlvl1=getval("gzlvl1"); gt1=getval("gt1");
   iphase = (int) (phase + 0.5); mix=getval("mix");
   iphase2 = (int) (phase2 + 0.5); 
   if (j>0.0) tau=1.0/(2.0*j); 
   getstr("sspul",sspul); getstr("satmode",satmode);
   getstr("SSshape",SSshape);

   settable(t1,1,ph1);
   settable(t2,2,ph2);
   settable(t3,32,ph3);
   settable(t4,8,ph4);
   settable(t5,4,ph5);
   settable(t6,16,ph6);
   settable(t7,8,ph7);
   settable(t9,4,ph9);

/* states-tppi */
   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v14);
   initval(2.0*(double)(((int)(d3*getval("sw2")+0.5)%2)),v13);

   getelem(t2,ct,v2); getelem(t5,ct,v5); getelem(t7,ct,oph);
   if (iphase == 2) incr(v5); if (iphase2 == 2) incr(v2);

   add(v14,v5,v5); add(v14,oph,oph); add(v13,v2,v2); add(v13,oph,oph);

 status(A);
      rlpower(pwx2lvl,DO2DEV);
     if (sspul[A] == 'y')
     {
      rlpower(tpwr,TODEV); 
      rgpulse(15.0*pw,zero,rof1,0.0);
      rgpulse(15.3*pw,one,0.0,rof2);
     }
      zgradpulse(gzlvl0,gt0);             /* homospoil  */
      delay(d1);
     if (satmode[A] == 'y')
     {
      if (fabs(satfrq-tof)>0.0) offset(satfrq,TODEV);
      rlpower(satpwr,TODEV);
      rgpulse(satdly,t9,rof1,rof1);
      rlpower(tpwr,TODEV); 
      if (fabs(satfrq-tof)>0.0) offset(tof,TODEV);
     }
      rcvroff();
 status(B);
      rgpulse(pw, t1, rof1, 0.0);
      delay(tau);
      dec2rgpulse(pwx2,v2, 0.0,0.0);
      if ((d3/2.0)>0)
       {
         delay(d3/2 -  pw -2.0*pw/3.1416);
         rgpulse(2*pw, t3,0.0,0.0);
         delay(d3/2 -  pw -2.0*pw/3.1416);
       }
       else rgpulse(2*pw,t3,0.0,0.0);
      dec2rgpulse(pwx2, t1,0.0,0.0);
      delay(tau);
 status(C);
      rgpulse(pw,t1,0.0,0.0);
      if (satmode[C] == 'y')
       {
        rlpower(satpwr,TODEV);
        delay(hst);
        zgradpulse(gzlvl0,gt0);
        rgpulse(mix-hst-gt0,t9,0.0,0.0);
        rlpower(SSpwr,TODEV);
        delay(rof1);
       }
      else
       {
        delay(hst);
        zgradpulse(gzlvl0,gt0);
        rlpower(SSpwr,TODEV);
        delay(mix-hst-gt0);
       }
      shaped_pulse(SSshape,SSpw,t4,0.0,0.0);
      rlpower(SSpwr+6.0,TODEV);            /* 6db for SS as 180 */
      zgradpulse(gzlvl1,gt1);
      delay(tau-POWER_DELAY-gt1);
      dec2rgpulse(pwx2,v5,0.0,0.0);
    if ((d2/2.0) > (SSpw/2.0 + pwx2))
     { 
      delay((d2/2.0) -SSpw/2.0 -(2.0*pwx2/3.14159));
      shaped_pulse(SSshape,SSpw,t6,0.0,0.0);
      delay((d2/2.0) -SSpw/2.0 -(2.0*pwx2/3.14159));
     }
      else shaped_pulse(SSshape,SSpw,t6,0.0,0.0);
      dec2rgpulse(pwx2,t1,0.0,0.0);
      zgradpulse(gzlvl1,gt1);
      delay(tau-POWER_DELAY-gt1);
      rlpower(dpwr2,DO2DEV);
 status(D);
}

