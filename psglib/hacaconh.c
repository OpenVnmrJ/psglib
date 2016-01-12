/* ha(caco)nh.c  Grzesiek and Bax (J. Biol Molecular NMR 185,204)
                      (modified to give ha correlations only)

                 Water suppression with spinlock purge pulses and gradients only

                 uses shifted laminar pulses to do the co hard pulses
                  
                 Carbon   - DODEV  (1st decoupler)
                 Nitrogen - DO2DEV (2nd decoupler)

   Parameters:

      sspul = 'y':  selects for Trim(x)-Trim(y) sequence at the start
		    of the pulse sequence
         hs = 'yn':  homospoil pulse (hst) during the d1 relaxation delay
      f2180 = 'n':  standard t2 timing
              'y':  modified t2 timing for t2(1) = half the dwell time
       fad1 = 'y':  TPPI axial-peak displacement along t1
              'n':  standard phasecycle
       fad2 = 'y':  TPPI axial-peak displacement along t2 (3D experiment)
              'n':  standard phasecycle
        nni = set equal to final ni
      hdlvl = power level for H1 broadband decoupling during t1, and t2
    hdshape = H1 broadband decoupling sequence at TOF
      hdres = tip-angle resolution for H1 broadband decoupling sequence
       hd90 = 90 degree pulse length for H1 broadband decoupling
       tpwr = power level for H1 transmitter pulses
         pw = 90 degree xmtr pulse length for protons (the observed nucleus)
    pwcalvl = power level for Ca hard pulses
       pwca = 90 degree decoupler pulse length for aliphatic C at `pwcalvl`
      cafrq = frequency for Carbon alpha spins
    pwcolvl = power level for shaped C0 pulses (applied dof=cabfrq)
       pwco = 180 degree pulse length for shaped C0 pulse at `pwcolvl`
    coshape = decoupling pattern for selective CO decoupling
   pwcolvlh = power level for hard C0 pulses (dof=cofrq)
      pwcoh = 90 degree pulse length for hard C0 pulse at `pwcolvlh`
      cofrq = freq of co for hard pulses
       padj = phase 7 on co pulse
   pwn15lvl = power level for N15 hard pulses
      pwn15 = 90 degree pulse length for N15
        jnh = one-bond heteronuclear coupling constant to NH (in Hz) (set to 92)
        jch = one-bond heteronuclear coupling constant to CH (in Hz) (set to 160)
      delta = 1/(2*JCH)
    epsilon = refocus time for HC, set to 2.1 msec 
      kappa = 1/(2*JNH) (refocus of antiphase nitrogen proton magnetization)
     lambda = 1/(4*JNH)
       zeta = polarization transfer cb to ca, ca to co (3.7 msec)
      theta = co to n polarization transfer (11.4 msec)
    tndelay = nitrogen constant time period (11.2 msec)
       tau2 = t2 evolution time
      phase = 1,2:  hypercomplex experiment with F1 quadrature (complex F1-FT)
     phase2 = 1,2:  hypercomplex experiment with F2 quadrature (complex F2-FT)
        zap = length of spin lock purge pulse
   amidefrq = amide frequency
     gzlvl1 = size of gradient zz filter after 1st inept
     gzlvl2 = size of two gradients to clean up pi pulses
     gzlvl3 = zz filter
     gzlvl4 = zz filter
     gzlvl5 = purge gradient before detection
        grt = gradient recovery time (500us)

     Gordon Rule, U. Virginia
*/

#include <standard.h>
#include <math.h>
static int      phs1[1] = {1},
		phs2[2] = {0,2},
		phs3[1] = {1},
                phs5[8] = {0,0,0,0,2,2,2,2},
                phs6[4] = {0,0,2,2},
                phs8[8] = {0,0,0,0,2,2,2,2},
                phs9[16]= {0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2},
                phs10[8]= {0,2,2,0,2,0,0,2};

static double	sim3_delay = 15e-6, pi = 3.1415926;

/*-------------------------------
|				|
|       pulsesequence()/0	|
|				|
+------------------------------*/
pulsesequence()
{
/* VARIABLE DECLARATION */
   char 	sspul[MAXSTR],
		fad1[MAXSTR],
		fad2[MAXSTR],
		f1180[MAXSTR],
		f2180[MAXSTR],
		coshape[MAXSTR],
		hdshape[MAXSTR];
   int          iphase,iphase2,t1index;
   double       ss,epsilon,kappa,lambda,eta,tndelay,delta,
		sstrim,theta,zeta, padj,acqt,t1a,t1b,t1c,ni,nni,
		sw1,sw2,tau2,
		t2dly1,
                pwcalvl,pwca,
		cafrq,cofrq,
                pwcolvl,pwcolvlh,
                pwco,pwcoh,
                pwn15lvl,
                pwn15,
                jnh,jch,
		hd90,
		hdlvl,
		hdres,
                amidefrq,zap,gzlvl1,gzlvl2,gzlvl3,gzlvl4,gzlvl5,
                gt1,gt2,gt3,gt4,gt5,grt;

/* LOAD VARIABLES */
   hdlvl = getval("hdlvl");
   hd90 = getval("hd90");
   hdres = getval("hdres");
   pwcalvl = getval("pwcalvl");
   pwca = getval("pwca");
   pwcolvl = getval("pwcolvl");
   pwco = getval("pwco");
   pwcolvlh = getval("pwcolvlh");
   pwcoh = getval("pwcoh");
   cafrq = getval("cafrq");
   cofrq = getval("cofrq");
   padj = getval("padj");
   amidefrq = getval("amidefrq");
   pwn15lvl = getval("pwn15lvl");
   pwn15 = getval("pwn15");
   tndelay = getval("tndelay");
   jnh = getval("jnh");
   jch = getval("jch");
   eta = getval("eta");
   epsilon = getval("epsilon");
   theta = getval("theta");
   zeta = getval("zeta");
   ss = getval("ss");
   sstrim = getval("sstrim");
   sw1 = getval("sw1");
   sw2 = getval("sw2");
   ni = getval("ni");
   nni = getval("nni");
   iphase = (int) (getval("phase") + 0.5);
   iphase2 = (int) (getval("phase2") + 0.5);
   zap = getval("zap");
   gzlvl1 = getval("gzlvl1");
   gzlvl2 = getval("gzlvl2");
   gzlvl3 = getval("gzlvl3");
   gzlvl4 = getval("gzlvl4");
   gzlvl5 = getval("gzlvl5");
      gt1 = getval("gt1");
      gt2 = getval("gt2");
      gt3 = getval("gt3");
      gt4 = getval("gt4");
      gt5 = getval("gt5");
      grt = getval("grt");
   getstr("coshape", coshape);
   getstr("hdshape", hdshape);
   getstr("sspul", sspul);
   getstr("fad1", fad1);
   getstr("fad2", fad2);
   getstr("f1180", f1180);
   getstr("f2180", f2180);

/* INITIALIZE VARIABLES */
   delta = 1.0/(2.0*jch);
   kappa = 1.0/(2.0*jnh);
   lambda = 1.0/(4.0*jnh);
   tau2=d3;
   if(f2180[0] == 'y')tau2=(d3+(1/(2*sw2)));
   initval(1.0,v7);     /*small angle phase shift for co pulse*/
   stepsize(padj,DODEV);
/* LOAD PHASE TABLE */
   settable(t1, 1, phs1);
   settable(t2, 2, phs2);
   settable(t3, 1, phs3);
   settable(t5, 8, phs5);
   settable(t6, 4, phs6);
   settable(t8, 8, phs8);
   settable(t9,16, phs9);
   settable(t10,8, phs10);
   getelem(t1,ct,v1);
   getelem(t8,ct,v8);
   getelem(t10,ct,oph);

/* CHECK CONDITIONS */
   if(gzlvl1 <0)gzlvl1=gzlvl1*-1.0;
   if(gzlvl2 >0)gzlvl2=gzlvl2*-1.0;
   if(gzlvl3 <0)gzlvl3=gzlvl3*-1.0;
   if(gzlvl4 <0)gzlvl4=gzlvl4*-1.0;
   if(gzlvl5 <0)gzlvl5=gzlvl5*-1.0;
   if((gzlvl1>20000)||(gzlvl2>15000)||(gzlvl3>23000)||(gzlvl4>10000)||(gzlvl5>10000))
   {
      printf("Gradients too large.\n");
      abort(1);
   }
   if((gt1>.002)||(gt2>.002)||(gt3>.003)||(gt4>.002)||(gt5>.002))
   {
     printf("Gradients too long.\n");
     abort(1);
   }
   if ( (strcmp(dm2, "n") != 0) && (strcmp(dm2, "nn") != 0) &&
	(strcmp(dm2, "ny") != 0) )
   {
      text_error("`dm(N15)` must be to 'n', 'nn', or 'ny'\n");
      abort(1);
   }

   if ( (strcmp(dm, "n") != 0) && (strcmp(dm, "nn") != 0) )
   {
      text_error("`dm(C13)` must be to 'n' or 'nn'\n");
      abort(1);
   }
   if (dpwr2 >48)
   {
	printf("decoupler power too large\n");
        abort(1);
   }
   if(pwcolvlh > 54)
   {
       printf("carbon power too high for pi pulses at co\n");
       abort(1);
   }
   if(pwcalvl > 54)
   {
       printf("carbon power too high for pi pulses at ca\n");
       abort(1);
   }
   if( getval("ni2")*(1/sw2)-pwco<0)
   {
       printf("nitrogen sweep width is too wide for pwco\n");
       abort(1);
   }
   if( zap > .0015)
   {
       printf("water purge pulses too long, reduce zap to 1.5msec\n");
       abort(1);
   }
   if( nni != ni)printf("WARNING - nni not equal to ni\n");
/* DETERMINE STEADY-STATE MODE */
   if (ss < 0)
   {
      ss *= (-1);
      initval(ss, ssval);
      initval(ss, ssctr);
   }

/* ADD IN STATES-HABERKORN ELEMENT & FAD */
   if (iphase == 2)incr(v1);
   if (iphase2 == 2)incr(v8);
   if (fad1[A] == 'y')
   {
	initval(2.0*(double)((int)(d2*sw1+0.5)%2),v12);
        add(v1,v12,v1); add(oph,v12,oph);
   }
   if (fad2[A] == 'y')
   {
	initval(2.0*(double)((int)(d3*sw2+0.5)%2),v14);
        add(v8,v14,v8); add(oph,v14,oph);
   }
/* BEGIN ACTUAL PULSE SEQUENCE CODE */
   status(A);
      rlpower(tpwr, TODEV);		/*  H1 hard-pulse power level	*/
      txphase(v1);
      rlpower(pwn15lvl, DO2DEV);   	/* N15 pulse power level	*/
      rlpower(pwcalvl+6, DODEV);  	/* Cab pulse power+6 for pi	*/
      offset(cafrq, DODEV); 		/* Carbon frequency		*/
      offset(tof, TODEV);		/* Set proton on water		*/

/*Calculate delays for t1 evolution */
      acqt=(nni-1)/sw1;
      t1index=(int)(d2*sw1+0.5);
      t1a=(delta/2.0) + ( (acqt/2.0)*t1index/(nni-1) );
      t1b=((acqt-delta)/2.0)*t1index/(nni-1);
      t1c=(delta/2.0) - ( (delta/2.0)*t1index/(nni-1) );

      rcvroff();
      if (sspul[A] == 'y')
      {
         rgpulse(sstrim, zero, rof1, 1.0e-6);
         rgpulse(sstrim, one, rof1, rof2);
         hsdelay(d1);
      }
      else
      {
         hsdelay(d1);
      }
/* INEPT transfer from H1 to C13 */
      rgpulse(pw, v1, rof1, rof1);
      delay(t1a - (pwca/2.0) - rof1 - (2*pw/pi));
      decrgpulse(pwca, t2, 0.0, rof1);
      rlpower(pwcalvl,DODEV);
      delay(t1b - (pwca/2.0) - 2*rof1 - pw);
      rgpulse(2*pw, zero, rof1, rof1);
      delay(t1c - pw - 2*rof1);
      rgpulse(pw, t2, rof1, rof1);
      rlpower(hdlvl, TODEV);
      txphase(zero);
      zgradpulse(gzlvl1,gt1);
      delay(grt);

/*Transfer from ca to co */
      decrgpulse(pwca, zero, rof1, rof1);
      rlpower(pwcolvl,DODEV);	
      decshaped_pulse(coshape,pwco,zero,rof1,0.0);
      delay(epsilon);
      obsprgon(hdshape, hd90, hdres);
      xmtron();
      rlpower(pwcalvl+6,DODEV);
      delay( zeta - epsilon - (pwca/2.0) );
      decrgpulse(pwca,t5,0.0,rof1);
      rlpower(pwcolvl,DODEV);
      decshaped_pulse(coshape,pwco,zero,rof1,0.0);
      rlpower(pwcalvl,DODEV);
      delay(zeta);
      decrgpulse(pwca,zero,0.0,rof1);

/*Transfer from co to n */

        rlpower(pwcolvlh,DODEV);
        offset(cofrq,DODEV);
        xmtroff();
        obsprgoff();
        zgradpulse(gzlvl2,gt2);
        delay(grt);
        obsprgon(hdshape, hd90, hdres);
        xmtron();
        decrgpulse(pwcoh,t6,rof1,0.0);
        
        rlpower(pwcalvl+6,DODEV);	/*refocus co-ca coupling*/
        offset(cafrq,DODEV);
        delay(eta-(pwca/2.0));
        decrgpulse(pwca,zero,0.0,0.0);
        rlpower(pwcolvlh+6,DODEV);
        offset(cofrq,DODEV);
        delay(theta-eta-(pwca/2.0)-(pwcoh/2.0));
        xmtroff();
        obsprgoff();
        sim3pulse(0.0,pwcoh,2*pwn15,zero,zero,zero,0.0,0.0);
        obsprgon(hdshape, hd90, hdres);
        xmtron();
        rlpower(pwcolvlh,DODEV);
        
        delay(theta-(pwcoh/2.0));
        
        dcplrphase(v7);
        decrgpulse(pwcoh,zero,0.0,rof1);
        dcplrphase(zero);
        rlpower(pwcolvlh+6,DODEV);
        xmtroff(); 
        obsprgoff();
        rlpower(tpwr-8,TODEV);
        rgpulse(3*zap,zero,rof1,rof1);
        rgpulse(2*zap,one,rof1,rof1);
        zgradpulse(gzlvl3,gt3*1.8);
        delay(grt);
        rlpower(tpwr,TODEV);
        rgpulse(pw,zero,rof1,rof1);
        rlpower(hdlvl,TODEV);
        if(amidefrq != tof)offset(amidefrq,TODEV);
        zgradpulse(gzlvl3,gt3);
        delay(grt);
        obsprgon(hdshape, hd90, hdres);
        xmtron();
        dec2rgpulse(pwn15,v8,rof1,rof1);
        t2dly1=tndelay-(tau2/2.0)-2*rof1-sim3_delay;
        if(t2dly1 < 0.0)t2dly1=0.0;
        delay(t2dly1);
        xmtroff();
        obsprgoff();
        sim3pulse(0.0,pwcoh,2*pwn15,zero,zero,t9,rof1,rof1);
        obsprgon(hdshape, hd90, hdres);
        xmtron();
        rlpower(pwcalvl+6,DODEV);
        offset(cafrq,DODEV);
        if((tau2/2.0)+(pwca/2.0)<kappa)
        {
        
        delay(tndelay +(tau2/2.0)-kappa);
        xmtroff();
        obsprgoff();
        rlpower(tpwr, TODEV);
        if(kappa-(tau2/2.0)-(pwca/2.0)>0)delay(kappa-(tau2/2.0)-(pwca/2.0));
        decrgpulse(pwca,zero,rof1,rof1);
        delay((tau2/2.0)-(pwca/2.0));
	}
        else
        {
        
        if(tndelay-(pwca/2.0)>0)delay(tndelay-(pwca/2.0));
        decrgpulse(pwca,zero,rof1,rof1);
        if((tau2/2.0)-(pwca/2.0)-kappa>0)delay((tau2/2.0)-(pwca/2.0)-kappa);
        xmtroff();
        obsprgoff();
        rlpower(tpwr, TODEV);
        delay(kappa);
	}

        dec2rgpulse(pwn15, zero, rof1, rof1);

        zgradpulse(gzlvl4,gt4);
        delay(grt);

       rgpulse(pw,zero, rof1, rof1);
       delay(lambda-2*rof1-(2*pw/pi)-pwn15);
       sim3pulse(2*pw, 0.0, 2*pwn15, zero, zero, zero, rof1, rof1);
       rlpower(dpwr2, DO2DEV);
       dec2phase(zero);
       txphase(one);
       delay(lambda-2*rof1-(2*pw/pi)-pwn15);

       rgpulse(pw, one, rof1, rof1);
       zgradpulse(gzlvl5,gt5);
       delay(grt);
       txphase(zero);
       rgpulse(pw, zero, rof1, rof1);	
       rcvron();
   status(B);
}
