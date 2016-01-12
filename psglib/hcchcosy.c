
/* hcch_cosy.c


   See  Bax et al., J.Mag.Reson.87, 620-627 (1990).
        Bax et al., J.Mag.Reson.88, 425-431 (1990).
        Clore et al., Biochemistry 29, 8172-8184 (1990).

    Various flags need to be set to indicate which experiment is
    to be run:

    Use programmed decoupling of carbonyls during t2.  Use
    herm180 for co refocusing

    Set dm = 'nny', dmm = 'ccp' [13C decoupling during acquisition]

    Must set phase = 1,2 and phase2 = 1,2 for hypercomplex in t1 and t2.

    Typical acquisition times are 28 ms [t1], 10 ms [t2], and 47 ms [t3]
       with 128 complex [t1], 32 complex [t2], and 512 real [t3]. 

    Written by Mike Wittekind 1/7/91 [hcch_16.c]
    Rwritten to constant-time experiment on April 20, 1992 by L.M..
    converted into SLP format June 15, 1992 GG(varian)
*/

#include <standard.h>
#include <math.h>

/* Define phase-table */
   static int  	psi1[1] = {0},
		phi1[2] = {0,1},
		phi2[1] = {0},
		phi3[1] = {0},
		phi4[1] = {1},
		phi5[8] = {0,0,1,1,2,2,3,3},
		phi6[8] = {0,0,0,0,2,2,2,2},
		phi7[1] = {0},

		rec[8]  = {0,2,2,0,0,2,2,0};

   static double d2_init = 0.0, d3_init = 0.0;


pulsesequence()
{
/* DECLARE VARIABLES */

 char        satmode[MAXSTR],
             ca90sh[MAXSTR],   /* SLP Ca pulse (dof at CO freq.) */
             ca180sh[MAXSTR],   /* SLP Ca pulse (dof at CO freq.) */
             ca180sh2[MAXSTR],   /* SLP Ca softer pulse (dof at CO freq.) */
             coshape[MAXSTR],   /* dec pattern for CO decoupling  */
             fhfdwt1[MAXSTR],   /* Flag to indicate half dwell start in t1 */
             fhfdwt2[MAXSTR];   /* Flag to indicate half dwell start in t2 */

 int       
             iphase,
	     satmove,
             iphase2,
	     t1_counter,
	     t2_counter;

 double      tauhc,            /* 1 / 4*J[13C-H]     [~1.5ms]          */
	     tau1,		/* t1-evolution variable */
	     tau2,		/* t2-evolutionvariable */
	     satfrq,		/* presaturation frequency */
             pwca,             /* PW90 for 13C nucleus                 */
             pwca180,          /* PW180 for 13C nucleus w. null at CO*/
             pwmlev,           /* PW90 for 1H during mlev16             */
             pwco,             /* PW90 for 13C carbonyl decoupling      */
             satpwr,          /* low level 1H trans.power for presat    */
             pwcalvl,           /* power level for 13C pulses on dec1   */
             mlevpwr,         /*  power level for 1H mlev decoupling     */
             cores,          /* pattern resol. for C=O decoupling  */
	     h1resol,		/* mlev decoupler pattern resol. = 90o */
             pwcolvl,          /* power level for C=O pulses on dec1   */
	     s180lvl,		/* power level for 180 with null at CO */
             jch,              /* coupling for C-C  (set to 40 Hz   )  */
             jcc,              /* coupling for C-C  (set to 40 Hz   )  */
             delta1,delta2,sw1,sw2;

/* LOAD VARIABLES */
  sw1 = getval("sw1");
  sw2 = getval("sw2");
  pwca = getval("pwca");
  pwca180 = getval("pwca180");
  pwmlev = getval("pwmlev");  
  pwco = getval("pwco");
  cores = getval("cores");
  satfrq = getval("satfrq");
  satpwr = getval("satpwr");
  pwcalvl = getval("pwcalvl");
  s180lvl = getval("s180lvl");
  mlevpwr = getval("mlevpwr"); 
  pwcolvl = getval("pwcolvl");
  jch = getval("jch");       /*   Use 152   [overestimate; true J ~140 ] */
  jcc = getval("jcc");       /*   Use 152   [overestimate; true J ~140 ] */
  getstr("ca180sh",ca180sh);
  getstr("ca180sh2",ca180sh2);
  getstr("ca90sh",ca90sh);
  getstr("coshape",coshape);
  getstr("satmode",satmode);
  getstr("fhfdwt1",fhfdwt1);
  getstr("fhfdwt2",fhfdwt2);
  iphase = (int)(getval("phase") + 0.5);
  iphase2 = (int)(getval("phase2") + 0.5);


/* calculate delays */
  h1resol = (90.0);
  tauhc = (1.0/(4*jch)); 
  delta1 = (1.0/(3.0*jch));
  delta2 = (1.0/(3.3*jcc));
  satmove = (fabs(tof - satfrq) >= 0.2);
/* CHECK VALIDITY OF PARAMETER RANGE */

    if((dm[A] != 'n') || (dm[B] != 'n'))
    {
        printf("dm should be 'nnn' or 'nny' ");
        abort(1);
    }

    if((dm2[A] != 'n') || (dm2[B] != 'n') || (dm2[C] != 'n'))
    {
        printf("dm2 should be 'nnn'  ");
        abort(1);
    }

    if( satpwr > 20 )
    {
        printf("SATPWR too large !!!  ");
        abort(1);
    }


    if( pwcolvl > 42 )
    {
        printf("DCOPWR too large!");
        abort(1);
    }

    if( dpwr > 48 )
    {
        printf("don't fry the probe, DCPWR too large!  ");
        abort(1);
    }

    if(( pwca > 2.5e-5 || pw > 2.5e-5 ))
    {
        printf("Pulsewidths [pw, pwca] must be shorter than 25 us. Abort."); 
        abort(1);
    }


    if ( pwco > 2.0e-3 )
    {
        printf("PWCO must be less than 2 ms. Abort.");
        abort(1);
    }


/* INITIALIZE VARIABLES */

  if(fhfdwt1[0] == 'y')
    tau1 = (d2 + 1/(2*sw1) );
  else
      tau1 = d2;

  if(fhfdwt2[0] == 'y')
    tau2 = (d3 + 1/(2*sw2));
  else
       tau2 = d3;


	settable(t21,1,psi1);
	settable(t11,2,phi1);
	settable(t12,1,phi2);
	settable(t13,1,phi3);
	settable(t14,1,phi4);
	settable(t15,8,phi5);
	settable(t16,8,phi6);
	settable(t17,1,phi7);
	settable(t10,8,rec);
                                  /* Phase table:

                                     phi1 = t11 = 0 1


                                     phi2 = t12 = 0

                                     phi3 = t13 = 0 0 0 0 0 0 0 0

                                     phi3'= t14 = 1 1 1 1 1 1 1 1 


                                     phi5 = t15 = 0 0 1 1 2 2 3 3

                                     phi6 = t16 = 0 0 0 0 2 2 2 2

   THESE TO BE SOFTWARE-MODIFIED BASED ON t1, t2 VALUES:
                                     rec  = t10 = 0 2 2 0 0 2 2 0
				     psi1 = t21 = 0

                                     psi2 = t22 = 0
*/

 if(iphase==2) tsadd(t21,1,4);
 if(iphase2==2) tsadd(t12,1,4);
 
 if(ix==1)
 {
	d2_init = 0.0;
	d3_init = 0.0;
 }

 t1_counter = (int) ((d2-d2_init)*sw1 + 0.5);
 t2_counter = (int) ((d3-d3_init)*sw2 + 0.5);

 if(t1_counter%2)
 {
	tsadd(t21,2,4);
	tsadd(t10,2,4);
 }
 if(t2_counter%2)
 {
	tsadd(t12,2,4);
	tsadd(t10,2,4);
 }


/* BEGIN ACTUAL PULSE SEQUENCE */


status(A);
   rcvroff();
   if (satmode[A] == 'y')
   {
	if(satmove)
		offset(satfrq,TODEV);
	rlpower(satpwr,TODEV);
    	rgpulse(d1,zero,rof1,rof1);
	if(satmove)
		offset(tof,TODEV);
   }
   else
   {
    	delay(d1);
   }
   decphase(t13);
   rlpower(tpwr, TODEV);   /* Set transmitter power for hard 1H pulses */
   rlpower(pwcalvl, DODEV);
status(B);
   rgpulse(pw,t21,1.0e-5,0.0);                /* First 1H 90 degree pulse     */
   txphase(t11);
   if(tau1 < (2*pwca + pw))	/* first t1-value */
   {
	delay(tauhc - 2*pwca - 1.0e-6);  /* assuming pwca > pw */
     	decshaped_pulse(ca90sh,pwca,t13,0.0,0.0); /* 13C composite 180 deg pulse  */
	sim3shaped_pulse("hard",ca180sh,"hard",2*pw, 2*pwca,0.0,t11, t14,zero, 1.0e-6, 0.0);
     	decshaped_pulse(ca90sh,pwca,t13,0.0,0.0); /* 13C composite 180 deg pulse  */
	txphase(one);
	decphase(t12);
	delay(tauhc - 2*pwca - 1.0e-6);  /* assuming pwca > pw */
   }
   else
   {
	delay(tauhc + tau1/2 - 2*pwca -1.0e-6); /* t1-evol. plus pol. trans. */
     	decshaped_pulse(ca90sh,pwca,t13,0.0,0.0); /* 13C composite 180 deg pulse  */
     	decshaped_pulse(ca180sh,2*pwca,t14,1.0e-6,0.0); /* 13C composite 180 deg pulse  */
     	decshaped_pulse(ca90sh,pwca,t13,0.0,0.0); /* 13C composite 180 deg pulse  */
	delay(tau1/2 - 2*pwca - pw);       /* continued t1-evol. */
        rgpulse(2*pw,t11,0.0,0.0);      /* proton echo pulse */
	txphase(one);
	decphase(t12);
	delay(tauhc - pw);
   }

   sim3shaped_pulse("hard",ca90sh,"hard",pw,pwca,0.0,one,t12,zero,1.0e-6,0.2e-6);

	decphase(t15);
      	rlpower(pwcolvl, DODEV);
      	rlpower(mlevpwr, TODEV);
      	decprgon(coshape, pwco, cores);
	 if (pwcolvl>0) decon();
      	txphase(zero);
      	delay(delta1 - 3*POWER_DELAY - PRG_START_DELAY);
	obsprgon("mlev16", pwmlev, h1resol);
	xmtron();
	delay(delta2/2 - delta1 + tau2/2 - PRG_STOP_DELAY - pwca);
	if (pwcolvl>0) decoff();
	decprgoff();
      	rlpower(pwcalvl, DODEV);
      	decshaped_pulse(ca180sh, 2*pwca,t15,0.0,0.0);
	decphase(t17);
        txphase(zero);
      	rlpower(pwcolvl, DODEV);
        decprgon(coshape, pwco, cores);
	 if (pwcolvl>0) decon();
	delay(delta2/2 - tau2/2 - 2*POWER_DELAY - PRG_START_DELAY - PRG_STOP_DELAY);
	if (pwcolvl>0) decoff();
	decprgoff();
      	rlpower(pwcalvl, DODEV);
        decshaped_pulse(ca90sh,pwca,t17,0.0,0.0); /*90 deg 13C pulse*/
	rlpower(s180lvl, DODEV);    /* 180 with null at CO  */
	decphase(zero);
	delay(delta2/2 - POWER_DELAY - pwca180);

        decshaped_pulse(ca180sh2,pwca180,zero,0.0,0.0); /*sel. 180 deg 13C pulse*/
	rlpower(pwcalvl, DODEV);
        delay(delta2/2 - delta1 - 2*POWER_DELAY - pwca180 - PRG_STOP_DELAY);
	xmtroff();
        obsprgoff();
	rlpower(tpwr, TODEV);
	delay(delta1);
     /* start reversed INEPT */
        sim3shaped_pulse("hard",ca90sh,"hard",pw,pwca,0.0,t13,zero,zero,0.0,0.0);
   decphase(zero);
   txphase(zero);
   delay(tauhc - 2*pwca - 1.0e-6);    /* delay = 1/4J     */
   decshaped_pulse(ca90sh,pwca,zero,0.0,0.0);
   sim3shaped_pulse("hard",ca180sh,"hard",2*pw,2*pwca,0.0,zero,one,zero,1.0e-6,0.0);
                                                  /* 1H + 13C 180 deg pulses      */
   decshaped_pulse(ca90sh,pwca,zero,1.0e-6,rof2);
   rcvron();
   delay(tauhc - 2*pwca - 1.0e-6 - rof2);     /* delay = 1/4J       */
   decshaped_pulse(ca90sh,pwca,zero,0.0,0.0);  /* Filter out IySz terms        */
   decshaped_pulse(ca90sh,pwca,t16,1.0e-6,0.0);
   rlpower(dpwr, DODEV);           /* Set power for decoupling     */

/* BEGIN ACQUISITION */
   setreceiver(t10);
status(C);
}


