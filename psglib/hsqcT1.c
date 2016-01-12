/* hscT1 -          heteronuclear Overbodenhausen experiment using
                    refocused INEPT and broadband decoupling of H1
                    during t1;  uses 3rd channel for X (2nd decoupler)
                    X magnetization put on Z axis for relaxation period
                    then put back in XY plane before reverse inept
                 
                    can use proton scrambling pulses at end of mix period
                    and trim pulse for water dephasing during first inept.
   Parameters:

     satmode = 'y':  presaturation during relaxation period (satdly) with xmtr
              'n':  no presaturation during relaxation period (satdly)
     satfrq = presaturation frequency
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
       tpwr = power level for H1 transmitter pulses
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
    hdshape = H1 programmable decoupling shape
      hdpwr = power level for H1 programmable decoupling
       hd90 = 90 degree pulse length for H1 programmable decoupling shape
              at power level 'hdpwr'
      hdres = tip-angle resolution for the H1 programmable decoupling
              shape
    scramble= time for (trim(x)-trim(y)-delay) loop
       trim = water spoiling trim pulse (typ. 1ms)
      mix = relaxation delay following inept transfer 
     pwx2lvl = power level for X decoupler pulses
        pwx2 = 90 degree decoupler pulse length for X
        jxh = one-bond heteronuclear coupling constant to X (in Hz)
        dm  = 'nnnn':  no broadband decoupling of X during acquisition
              'nnny':  broadband heteronuclear decoupling of X during
		        acquisition


      phase = 1,2:  hypercomplex experiment with F1 quadrature (complex F1-FT)
*/


#include <standard.h>

static int	phs1[8]  = {0,0,0,0,2,2,2,2},
		phs2[32] = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,
                            2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3},
		phs3[2]  = {1,3},
                phs4[4]  = {0,0,2,2},
		phs5[32] = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,
                            2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3},
                phs6[8]  = {0,2,2,0,2,0,0,2},
		phs7[8]  = {0,0,0,0,1,1,1,1};
                phs8[8]  = {1,1,1,1,0,0,0,0};

pulsesequence()
{
/* VARIABLE DECLARATION */
   char         satmode[MAXSTR],
		hdshape[MAXSTR];
   int          phase;
   double  
                scramble,
                pwx2lvl,
                cycles,
                pwx2,
                trim,
                mix,
                jxh,
                ni,
                sw1,
		bird,
                satfrq,
                satdly,
                satpwr,
		hdpwr,
		hd90,
		hdres;

/* Load variables */
   scramble = getval("scramble");
     trim = getval("trim");
    mix = getval("mix");
   satfrq = getval("satfrq");
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   pwx2lvl = getval("pwx2lvl");
   pwx2 = getval("pwx2");
    ni = getval("ni");
   sw1 = getval("sw1");
   jxh = getval("jxh");
   hd90 = getval("hd90");
   hdres = getval("hdres");
   hdpwr = getval("hdpwr");
   phase = (int) (getval("phase") + 0.5);
   getstr("satmode", satmode);
   getstr("hdshape", hdshape);	/* H1 t1 decoupling pattern	*/


/* Load phase tables */
   settable(t1, 8, phs1);
   settable(t2, 32, phs2);
   settable(t3, 2, phs3);
   settable(t4, 4, phs4);
   settable(t5, 32, phs5);
   settable(t6, 8, phs6);
   settable(t7, 8, phs7);
   settable(t8, 8, phs8);

      bird = 1/(2*jxh);
     cycles= scramble/(32*pw);
     cycles= 2.0*(double)(int)(cycles/2.0);
     scramble=cycles*(32*pw);
     initval(cycles,v1);

/* Phase incrementation for hypercomplex 2D data */
   if (phase == 2)
      tssub(t1, 1, 4);

/* check for too much proton power   */

  if (trim>0.03)
   {
      printf("The proton trim pulse is too long. Use < 0.003.\n");
      abort(1);
   }
  if (scramble>0.02)
   {
      printf("The proton scrambing sequence is too long. Use < 0.02.\n");
      abort(1);
   }
   if (hdpwr>40)
   {
      printf("The proton decoupling power is too high. Use < 40.\n");
      abort(1);
   }

/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      rlpower(pwx2lvl, DO2DEV);
      rlpower(tpwr, TODEV);
      if (ni>1)
      {
       rlpower(hdpwr, TODEV);
       obsprgon(hdshape, hd90, hdres);
       xmtron();
       delay((ni/sw1)-d2);  /* total 1H decoupling time independent of d2 */
       xmtroff();
       obsprgoff();
       rlpower(tpwr, TODEV);
      }
      hsdelay(d1);

/* selective saturation period */
      if (satmode[A] == 'y')
      {
         if (fabs(tof-satfrq)>0.1)   offset(satfrq, TODEV);
         rlpower(satpwr, TODEV);
         rgpulse(satdly, zero, 4.0e-5, 0.2e-6);
         if (fabs(tof-satfrq)>0.1)   offset(tof, TODEV);
         rlpower(tpwr,TODEV);
      }

   status(B);
      rcvroff();
      rgpulse(pw, zero, rof1, 0.0);
      delay(0.4*bird);
      sim3pulse(2*pw, 0.0,2*pwx2,zero, zero, zero, 0.0, 0.0);
      delay(0.4*bird);
      rgpulse(trim,zero,0.0,0.0);
      sim3pulse(pw, 0.0,pwx2,one ,zero, t1, 0.0, 0.0);
      delay( bird/2);
      sim3pulse(2*pw, 0.0,2*pwx2,zero, zero,t2, 0.0, 0.0);
      rlpower(hdpwr, TODEV);
      delay(bird/2);
      obsprgon(hdshape, hd90, hdres);
      xmtron();
      delay(d2);
      dec2rgpulse(pwx2,t3,0.0,0.0);
      delay(mix);
      xmtroff();
      obsprgoff();
      rlpower(tpwr,TODEV);
     if (scramble>0.0)
       {
       initval(50.0,v1);
       starthardloop(v1);
        rgpulse(8*pw, t7, 0.0, 0.0);
        rgpulse(8*pw, t8, 0.0, 0.0);
        delay(16*pw);
       endhardloop();
       }
      dec2rgpulse(pwx2,t4,0.0,0.0);

   status(C);
      delay(bird/2);
      sim3pulse(2*pw, 0.0,2*pwx2, zero,zero, zero, 0.0, 0.0);
      delay( bird/2 );
      sim3pulse(pw,0.0, pwx2, zero,zero, one, 0.0, 0.0);
      delay(0.4*bird);
      sim3pulse(2*pw, 0.0, 2*pwx2, t5,zero, zero, 0.0, 0.0);
      rlpower(dpwr2,DO2DEV);
      rcvron();
      delay(0.4*bird);

   status(D);
      setreceiver(t6);
}
