/*   Heteronuclear Polarization Transfer by Isotropic Mixing

    ref.: M. Ernst, C. Griesinger, R. R. Ernst, and W. Bermel, 
          Mol. Phys. 74, 212-252 (1991). (sequence of Fig. 1c.)

   Parameters:

      satmode  = flag for optional solvent presaturation
      satdly   = presaturation delay
      satfrq   = presaturation frequency
      satpwr   = power level for presaturation
      sspul    = optional GRAD-90(proton)-GRAD sequence at the beginning of d1
      gt0      = gradient duration for sspul
      gzlvl0   = gradient power for sspul
      pw       = 90 deg. proton pulse at tpwr
      pwx      = 90 deg. x pulse on channel 2 (at pwxlvl)
      pwxlvl   = power lovel for pwx on channel 2
      pwx2     = 90. deg. x pulse on channel 3 (at pwx2lvl)
      pwx2lvl  = power level for pwx2 om channel 3
      mix      = first istropic mixing period
      mix2     = second isotropic mixing period
      slpw     = 90 deg. proton pulse for mixing (at slpwr)
      slpwr    = power level for slpw on channel 1
      slpwx    = 90 deg. x pulse for mixing on channel 2 (at slpwrx)
      slpwrx   = power level for mixing on channel 2 
      slpwx2   = 90 deg. x pulse for mixing on channel 3 (at slpwrx2)
      slpwrx2  = power level for mixing on channel 3
      Htrim    = proton trim pulse at tpwr-6 at end of 2nd mix
      Xtrim    = X trim pulse at pwxlvl-6 at start of 2nd mix
      X2trim   = X2 trim pulse at pwx2lvl-6 at start of 2nd mix

      E. Hoffmann and P. Sandor (Varian Darmstadt)
*/

#include <standard.h>

dipsy2hx(ph1,ph2,ph3,ph4)
codeint ph1,ph2,ph3,ph4;
{
   double slpw = getval("slpw"),
          slpwx = getval("slpwx"),
          slpwx2 = getval("slpwx2");
   sim3pulse(3.556*slpw,3.556*slpwx,3.556*slpwx2,ph1,ph3,ph3,rof1,rof2);
   sim3pulse(4.556*slpw,4.556*slpwx,4.556*slpwx2,ph2,ph4,ph4,rof1,rof2);
   sim3pulse(3.222*slpw,3.222*slpwx,3.222*slpwx2,ph1,ph3,ph3,rof1,rof2);
   sim3pulse(3.167*slpw,3.167*slpwx,3.167*slpwx2,ph2,ph4,ph4,rof1,rof2);
   sim3pulse(0.333*slpw,0.333*slpwx,0.333*slpwx2,ph1,ph3,ph3,rof1,rof2);
   sim3pulse(2.722*slpw,2.722*slpwx,2.722*slpwx2,ph2,ph4,ph4,rof1,rof2);
   sim3pulse(4.167*slpw,4.167*slpwx,4.167*slpwx2,ph1,ph3,ph3,rof1,rof2);
   sim3pulse(2.944*slpw,2.944*slpwx,2.944*slpwx2,ph2,ph4,ph4,rof1,rof2);
   sim3pulse(4.111*slpw,4.111*slpwx,4.111*slpwx2,ph1,ph3,ph3,rof1,rof2); 
}

DIPSY2hx(x,mx,x1,mx1)
codeint x,mx,x1,mx1;
{
   dipsy2hx(x,mx,x1,mx1); dipsy2hx(mx,x,mx1,x1); dipsy2hx(mx,x,mx1,x1); dipsy2hx(x,mx,x1,mx1);
}

pulsesequence()
{
/* VARIABLE DECLARATION */
   double 
          pwx2 = getval("pwx2"),
          pwx2lvl = getval("pwx2lvl"),
          mix = getval("mix"),
          mix2 = getval("mix2"), 
          slpw = getval("slpw"),
          slpwr = getval("slpwr"),
          slpwx = getval("slpwx"),
          slpwrx = getval("slpwrx"),
          slpwx2 = getval("slpwx2"),
          slpwrx2 = getval("slpwrx2"), 
          Htrim = getval("Htrim"),
          Xtrim = getval("Xtrim"),
          X2trim = getval("X2trim"),  
          gt0 = getval("gt0"),
          gzlvl0 = getval("gzlvl0"),
          cycles, cycles2; 
   char   sspul[MAXSTR]; 
   int    iphase = (int) (getval("phase") + 0.5);

   getstr("sspul",sspul);

   loadtable("heteroHH");
   sub(ct,ssctr,v12);
   getelem(t1,v12,v1);
   getelem(t3,v12,v4);
   getelem(t4,v12,v6);
   getelem(t5,v12,oph);
   assign(zero,v14);
   if (iphase == 2)
      {incr(v4); incr(v14);}
 
/* HYPERCOMPLEX MODE USES REDFIELD TRICK TO MOVE AXIAL PEAKS TO EDGE */
   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v13);
   if ((iphase==1)||(iphase==2))
      {add(v4,v13,v4); add(v14,v13,v14); add(oph,v13,oph);}

/* CHECK CONDITIONS */

   if ((slpw > slpwx) && (slpw > slpwx2))
       cycles = mix/(115.112*slpw);
   else
      {
      if (slpwx > slpwx2)
         cycles = mix/(115.112*slpwx);
      else cycles = mix/(115.112*slpwx2);
      } 

   if ((slpw > slpwx) && (slpw > slpwx2))
       cycles2 = mix2/(115.112*slpw);
   else
      {
      if (slpwx > slpwx2)
         cycles2 = mix2/(115.112*slpwx);
      else cycles2 = mix2/(115.112*slpwx2);
      } 

/* Check for correct DM settings */
   if ((dm[A] == 'y') || (dm[B] == 'y') || (dm2[A] == 'y') || (dm2[B] == 'y'))
    {
      text_error("DM and DM2 must be set to either 'nny' or 'nnn'.");
      abort(1);
    }
 
/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      initval(cycles+0.5,v10);
      initval(cycles2+0.5,v11);
      rlpower(tpwr,TODEV);
      rlpower(pwxlvl,DODEV); rlpower(pwx2lvl,DO2DEV);
      if (sspul[A] == 'y')
        {
          zgradpulse(gzlvl0,gt0);
          delay(5.0e-5);
          rgpulse(pw,zero,rof1,rof1);
          zgradpulse(gzlvl0,gt0);
        }
      if (d1>satdly+hst) hsdelay(d1-satdly);
      if (satmode[A] == 'y')
        {
          rlpower(satpwr,TODEV);
          if (satfrq != tof) offset(satfrq,TODEV);
          rgpulse(satdly,v14,rof1,rof1);
          if (satfrq != tof) offset(tof,TODEV);
          rlpower(tpwr,TODEV);
        }
          
   status(B);
      rgpulse(pw,v1,rof1,rof1);
      if (cycles>0.6)
      {
         rlpower(slpwr,TODEV); 
         rlpower(slpwrx,DODEV); rlpower(slpwrx2,DO2DEV);
         getelem(t2,v12,v2);
         add(v2,two,v3);
         add(v4,two,v5);
/*         loop(v10,v7);
            DIPSY2hx(v2,v3,v4,v5);
         endloop(v7); */
         starthardloop(v10);
            DIPSY2hx(v2,v3,v4,v5); 
         endhardloop(); 
         delay(2.0e-5);
         rlpower(tpwr,TODEV);
         rlpower(pwxlvl,DODEV); rlpower(pwx2lvl,DO2DEV);
      }
      delay(d2/2.0);
      sim3pulse(2.0*pw,2.0*pwx,2.0*pwx2,v6,v6,v6,rof1,rof1);
      delay(d2/2.0);                                             
      rlpower(pwxlvl-6.0,DODEV); rlpower(pwx2lvl-6.0,DO2DEV);
      sim3pulse((double)0.0,Xtrim,X2trim,zero,v6,v6,rof1,rof2); 
      if (cycles2>0.6)
      {
         rlpower(slpwr,TODEV); 
         rlpower(slpwrx,DODEV); rlpower(slpwrx2,DO2DEV);
   /*   loop(v10,v8); 
            DIPSY2hx(v2,v3,v4,v5); 
         endloop(v8);  */
         starthardloop(v11);
            DIPSY2hx(v2,v3,v2,v3); 
         endhardloop();  
         delay(2.0e-5);
      } 
      rlpower(tpwr-6.0,TODEV);
      rgpulse(Htrim,v6,rof1,rof2);
      rlpower(dpwr,DODEV); rlpower(dpwr2,DO2DEV);
      decphase(zero); dec2phase(zero);
   status(C);
}

/* phase table
t1 = 1 1 3 3 2 2 0 0 3 3 1 1 0 0 2 2
t2 = 0 0 2 2 1 1 3 3 2 2 0 0 3 3 1 1
t3 = 0 2 2 0 1 3 3 1 2 0 0 2 3 1 1 3
t4 = 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
t5 = 0 2 2 0 1 3 3 1 2 0 0 2 3 1 1 3
*/           
