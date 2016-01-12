/* slphncoca.c      -  HN(CO)CA 3D experiment by Kay et al.; uses a
	   	    frequency-shifted pulse on the Ca spins;
                  Krish - Nov. 27 1991
                          NOTE : decoupler = N15
                                 decoupler2= C13
                  uses only FAD in F1 and F2
           (SLP presaturation included )
   Parameters:

      sspul = 'y':  selects for Trim(x)-Trim(y) sequence at the start
		    of the pulse sequence
     satflg = 'y':  H1 presaturation during relaxtion delay
     satfrq = frequency of 1H presaturation for all periods
     satdly = saturation time during the relaxation period
     satpwr = saturation power for all periods of presaturation with xmtr
     slpflg = 'y' : off-resonance presaturation with slp
     satshape = slpsatd (tof is downfield to water)
                slpsatu (tof is upfield to water)
     h2off = water off-resonance frequency (tof - satfrq)
     mult = multiplication factor shape_pw calculation (should be 10 for the two
             shapes noted above)
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
       tpwr = power level for 1H transmitter pulses
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
       dof2 = should be set to C0 frequency
    pwcalvl = power level for Ca decoupler pulses
       pwca = 90 degree decoupler pulse length for Ca at `pwcalvl`
    cashape = frequency-shifted pulse for Ca excitation
    pwcolvl = power level for C0 decoupler pulses
       pwco = 90 degree decoupler pulse length for C0 at `pwcolvl`
   pwco2lvl = power level for shorter CO pulse
      pwco2 = 90 degree decoupler pulse length for C0 at `pwco2lvl`
   pwn15lvl = power level for N15 decoupler pulses
      pwn15 = 90 degree decoupler pulse length for N15
       dpwr = power level for N15 broadband decoupling
       tau1 = first delay in sequence (NH scalar coupling)
       tau2 = second delay in sequence
       tau3 = third delay in sequence
       tau4 = fourth delay in sequence
         dm = 'nnnnn':  no broadband decoupling of N15 during acquisition
              'nnnny':  broadband heteronuclear decoupling of N15
		        during acquisition
        dm2 = 'n':  no broadband decoupling of C13 during acquisition

      phase = 1,2:  hypercomplex experiment with F1 quadrature (complex F1-FT)
     phase2 = 1,2:  hypercomplex experiment with F2 quadrature (complex F2-FT)
*/




#include <standard.h>

pulsesequence()
{
/* VARIABLE DECLARATION */
   char         satflg[MAXSTR],
		sspul[MAXSTR],
                slpflg[MAXSTR],
                satshape[MAXSTR],
		cashape[MAXSTR];
   int          phase,
		phase2,
		satmove;
   double       pwcalvl,
                pwca,
                pwcolvl,
                pwco,
                pwco2lvl,
                pwco2,
                pwn15lvl,
                pwn15,
                satfrq,
                satdly,
                satpwr,
                cycle,
                shape_pw,
                mult,
                h2off,
		tau1,
		tau2,
		tau3,
		tau4;


/* Load variables */
   satfrq = getval("satfrq");
   satdly = getval("satdly");
   satpwr = getval("satpwr");
   mult = getval("mult");
   h2off = getval("h2off");

   pwcalvl = getval("pwcalvl");
   pwca = getval("pwca");
   pwcolvl = getval("pwcolvl");
   pwco = getval("pwco");
   pwco2lvl = getval("pwco2lvl");
   pwco2 = getval("pwco2");
   pwn15lvl = getval("pwn15lvl");
   pwn15 = getval("pwn15");
   tau1 = getval("tau1");
   tau2 = getval("tau2");
   tau3 = getval("tau3");
   tau4 = getval("tau4");
   phase = (int) (getval("phase") + 0.5);
   phase2 = (int) (getval("phase2") + 0.5);

   getstr("cashape", cashape);
   getstr("sspul", sspul);
   getstr("satflg", satflg);
   getstr("slpflg",slpflg);
   getstr("satshape",satshape);

   shape_pw = mult/h2off;
   if (shape_pw < 0.0)
     shape_pw = -shape_pw;
   cycle = (satdly/(shape_pw + 15.4e-6));
   cycle = 2.0*(double)(int)(cycle/2.0);
   initval(cycle,v14);

/* Load phase cycles */

   loadtable("slphncoca");                /* t1 = 1 1 3 3
                                          t2 = 0 2
                                          t3 = 0 0 0 0 2 2 2 2
                                          t4 = 0 2 2 0 2 0 0 2   */



/* Check conditions */
   satmove = ( fabs(tof - satfrq) >= 0.1 );


/* Add in States-Haberkorn element */
   getelem(t2,ct,v2);
   getelem(t3,ct,v3);
   getelem(t4,ct,oph);

   if (phase2 == 2)		/* N15 t1 element */
      incr(v2);

   if (phase == 2)		/* C13 t2 element */
      incr(v3);


/* Add in FAD */

   initval(2.0*(double)(((int)(d3*getval("sw2")+0.5)%2)),v12);
   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v13);

     add(v2,v12,v2);
     add(oph,v12,oph);
     add(v3,v13,v3);
     add(oph,v13,oph);



/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      rlpower(tpwr, TODEV);		/*  H1 hard-pulse power level	*/
      rlpower(pwn15lvl, DODEV);  	/* N15 hard-pulse power level	*/
      rlpower(pwcolvl, DO2DEV); 	/* CO hard-pulse power		*/

      if (sspul[A] == 'y')
      {
         hsdelay(hst);
         rgpulse(pw, zero, rof1, rof2);
         hsdelay(hst + d1);
      }
      else
      {
         hsdelay(d1);
      }

/* selective saturation period */
      if (satflg[A] == 'y')
      {
        rlpower(satpwr, TODEV);
        if (slpflg[0] == 'n')
        {
         if (satmove)
            offset(satfrq, TODEV);
         rgpulse(satdly, zero, 4.0e-5, 0.2e-6);
         if (satmove)
            offset(tof, TODEV);
        }
        if (slpflg[0] == 'y')
        {
         rcvroff();
         delay(rof1);
         starthardloop(v14);
          shaped_pulse(satshape,shape_pw,zero,0.0,1e-6);
         endhardloop();
         delay(rof2);
         rcvron();
        }
         rlpower(tpwr, TODEV);
         delay(1.0e-5);
      }


/* Pulse train */
   status(B);
      rcvroff();

      rgpulse(pw, zero, rof1, 0.0);
      decphase(t1);	
      dec2phase(zero);
      delay(tau1);

      simpulse(2*pw, 2*pwn15, zero, t1, 0.0, 0.0);
      txphase(t1);
      decphase(v2);
      delay(tau1);

      simpulse(pw, pwn15, t1, v2, 0.0, 0.0);
      txphase(zero);
      decphase(zero);
      delay(tau2);

      delay(d3/2);
      rgpulse(2*pw, zero, 0.0, 0.0);
      delay(d3/2);


   status(C);
      delay(tau3);
      gensim2pulse(2*pwn15, 2*pwco, zero, zero, 0.0, 0.0, DODEV, DO2DEV);
      delay(tau2 + tau3);
      gensim2pulse(pwn15, pwco, zero, zero, 0.0, 0.0, DODEV, DO2DEV);

      rlpower(pwcalvl, DO2DEV);
      dec2phase(v3);
      delay(tau4 - 14.2e-6);
      dec2shaped_pulse(cashape,pwca,v3,0.0,0.0);


      dec2phase(zero);
      rlpower(pwco2lvl, DO2DEV);
      if (d2 > 0.0)
      delay(d2/2 - 2*pwca/3.1416 - pwco2 - 18.6e-6);
      else
      delay(d2/2);
      gensim2pulse(2*pw, 2*pwco2, zero, zero, 0.0, 0.0, TODEV, DO2DEV);

      rlpower(pwcalvl, DO2DEV);
      if (d2 > 0.0)
      delay(d2/2 - 2*pwca/3.1416 - pwco -  14.2e-6);
      else
      delay(d2/2);

      dec2shaped_pulse(cashape, pwca, zero, 0.0, 0.0);
      rlpower(pwcolvl, DO2DEV);
      delay(tau4 - 8.4e-6);


   status(D);
      gensim2pulse(pwn15, pwco, zero, zero, 0.0, 0.0, DODEV, DO2DEV);
      delay(tau3 + tau2);
      gensim2pulse(2*pwn15, 2*pwco, zero, zero, 0.0, 0.0, DODEV, DO2DEV);
      delay(tau3);

      rgpulse(2*pw, zero, 0.0, 0.0);
      delay(tau2);
      simpulse(pw, pwn15, zero, zero, 0.0, 0.0);
      delay(tau1);
      simpulse(2*pw, 2*pwn15, zero, zero, 0.0, 0.0);

      rlpower(dpwr, DODEV);
      rlpower(dpwr2, DO2DEV);
      delay(rof2);
      rcvron();
      delay(tau1 - rof2 - 8.4e-6);


   status(E);
}
