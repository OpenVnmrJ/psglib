/* ghsqcse
	heteronuclear single quantum with gradient selection.
        Kay,Keifer and Saarinen,JACS,114,10663(1992)
        includes water flipback selective pulse
        options for presat and scuba
        valid for vnmr5.1 or later
        modified phase cycling GG 12jan95 palo alto
*/

#include <standard.h>

static double d2_init = 0.0;

static int phi1[1] = {0},
           phi2[1] = {1},
           phi3[8] = {0,0,1,1,2,2,3,3}, 
           phi4[1] = {0},
           phi5[1] = {0},
           phi7[2] = {0,2},
           rec[4] = {0,2,2,0};

pulsesequence()
{
/* DECLARE VARIABLES */

 char        fscuba[MAXSTR],f1180[MAXSTR],flipshap[MAXSTR],
             flipback[MAXSTR],sspul[MAXSTR],c180_flg[MAXSTR];

 int	     icosel,t1_counter;

 double      hscuba,                /* length of 1/2 scuba delay */
             tauxh,                 /* 1 / 4J(XH)                */
             pwxN15,                /* PW90 for N-nuc            */
             pwxC13,                /* PW90 for C-nuc            */
             pwxNlvl,               /* power level for N hard pulses */
             pwxClvl,               /* power level for C hard pulses */
             jxh,                   /* coupling for XH           */
	     tau1,	      	    /* t1/2 */
	     sw1,
             flippw,                 /* pw for selective pulse at flippwr  */
             flippwr,

             gzlvl0,               /* level of grad. */
             gt0,                  /* grad time */
             gzlvl1,               /* level of grad. */
             gt1,                  /* grad time */
             gzlvl2,
             gt2,
             gzlvl3,
             gt3,
             gzlvl4,
             gt4,
             gzlvl5,
             gt5,
             BigT,                  /* constant time for N evolution */
             BigT1;


/* LOAD VARIABLES */

  jxh = getval("jxh");
  pwxN15 = getval("pwxN15");
  pwxC13 = getval("pwxC13");
  pwxNlvl = getval("pwxNlvl"); 
  pwxClvl = getval("pwxClvl"); 
  hscuba = getval("hscuba");
  sw1 = getval("sw1");
  BigT = getval("BigT");
  BigT1 = getval("BigT1");
  flippw = getval("flippw");
  flippwr = getval("flippwr");

  gt0 = getval("gt0");
  gt1 = getval("gt1");
  gt2 = getval("gt2");
  gt3 = getval("gt3");
  gt4 = getval("gt4");
  gt5 = getval("gt5");
  gzlvl0 = getval("gzlvl0");
  gzlvl1 = getval("gzlvl1");
  gzlvl2 = getval("gzlvl2");
  gzlvl3 = getval("gzlvl3");
  gzlvl4 = getval("gzlvl4");
  gzlvl5 = getval("gzlvl5");

  getstr("fscuba",fscuba); 
  getstr("f1180",f1180); 
  getstr("sspul",sspul);
  getstr("c180_flg",c180_flg);
  getstr("flipshap",flipshap);
  getstr("flipback",flipback);

/* check validity of parameter range */

    if((dm[A] == 'y' || dm[B] == 'y' || dm[C] == 'y'))
	{
	printf("incorrect Dec1 decoupler flags!  ");
	abort(1);
    } 

    if((dmm[A] == 'y' || dmm[B] == 'y' || dmm[C] == 'y'))
	{
	printf("incorrect Dec2 decoupler flags!  ");
	abort(1);
    } 

    if( satpwr > 8 )
    {
	printf("satpwr too large !!!  ");
	abort(1);
    }

    if( dpwr > 50 )
    {
	printf("don't fry the probe, dpwr too large!  ");
	abort(1);
    }

    if( dpwr2 > 50 )
    {
	printf("don't fry the probe, dpwr2 too large!  ");
	abort(1);
    }

    if(gt0 > 15.0e-3 || gt1 > 15.0e-3 || gt2 > 15.0e-3 || gt3 > 15.0e-3 || gt4 > 15.0e-3 || gt5 > 15.0e-3)
    {
        printf("gti must be less than 15 ms \n");
        abort(1);
    }


/* LOAD VARIABLES */

  settable(t1, 1, phi1);
  settable(t2, 1, phi2);
  settable(t3, 8, phi3);
  settable(t4, 1, phi4);
  settable(t5, 1, phi5);
  settable(t7, 2, phi7);
  settable(t6, 4, rec);

/* INITIALIZE VARIABLES */

  tauxh = ((jxh != 0.0) ? 1/(4*(jxh)) : 2.25e-3);

/* Phase incrementation for hypercomplex data */

   if ( phase1 == 1 )     /* Hypercomplex in t1 */
   {
        tsadd(t5, 2, 4); 
        icosel = 1; /*  and reverse the sign of the gradient  */
   } 
   else icosel = -1;


/* calculate modification to phases based on current t1 values
   to achieve States-TPPI acquisition */
 
 
   if(ix == 1)
      d2_init = d2;
      t1_counter = (int) ( (d2-d2_init)*sw1 + 0.5);

      if(t1_counter %2) {
        tsadd(t7,2,4);
        tsadd(t6,2,4);
      }


/* set up so that get (-90,180) phase corrects in F1 if f1180 flag is y */

   tau1 = d2;
   if(f1180[A] == 'y')  tau1 += ( 1.0/(2.0*sw1) );
   tau1 = tau1/2.0;
   

/* BEGIN ACTUAL PULSE SEQUENCE */

status(A);
   obspower(tpwr);               /* Set power for pulses  */
   decpower(pwxClvl);              /* Set decoupler1 power to pwxClvl */
   dec2power(pwxNlvl);            /* Set decoupler2 power to pwxNlvl */
status(B);
  if(sspul[A] == 'y')
  {
    rgpulse(200*pw,zero,rof1,2.0e-6);
    rgpulse(200*pw/1.62,one,2.0e-6,rof1);
  }
 if(satmode[A] == 'y')
{
  /* Presaturation Period */
  obspower(satpwr);
  txphase(zero);
  rgpulse(satdly,zero,2.0e-6,2.0e-6);
  obspower(tpwr);                /* Set power for hard pulses  */
    if (fscuba[0] == 'y')            /* Scuba pulse sequence */
    {  
      hsdelay(hscuba);
      rgpulse(pw,zero,1.0e-6,0.0);	/* 90x180y90x */
      rgpulse(2.0*pw,one,1.0e-6,0.0);
      rgpulse(pw,zero,1.0e-6,0.0);
      txphase(zero);
      delay(hscuba);        
    }
}
 delay(d1);
status(C);
  rcvroff();
  delay(20.0e-6);

  /* eliminate all magnetization originating on 15N */

  dec2rgpulse(pwxN15,zero,0.0,0.0);
  delay(2.0e-6);
  zgradpulse(gzlvl0,gt0);
  delay(150.0e-6);

  rgpulse(pw,zero,1.0e-6,0.0);    

  txphase(zero);
  dec2phase(zero);

  delay(2.0e-6);
  zgradpulse(gzlvl1,gt1);
  delay(2.0e-6);

  delay(tauxh - gt1 - 4.0e-6);               /* delay=1/4J(XH)   */

  sim3pulse(2.0*pw,0.0e-6,2*pwxN15,zero,zero,zero,0.0,0.0);

  txphase(t2);
  dec2phase(t7);

  delay(tauxh - gt1 - 202.0e-6);               /* delay=1/4J(XH)   */

  delay(2.0e-6);
  zgradpulse(gzlvl1,gt1);
  delay(200.0e-6);

  rgpulse(pw,t2,0.0,0.0);

  /* shaped pulse */
  if (flipback[A]=='y')
  {
    obspower(flippwr);
    shaped_pulse(flipshap,flippw,two,2.0e-6,0.0);
    delay(2.0e-6);
    obspower(tpwr);
  }

  delay(2.0e-6);
  zgradpulse(gzlvl2,gt2);
  delay(150.0e-6);

  dec2rgpulse(pwxN15,t7,0.0,0.0);

  txphase(t4); dec2phase(t3);
/* assumption is that pwxC13>pw */ 
  if(c180_flg[A] == 'y')
    {
     if(d2>0.0) delay(tau1-pwxC13);
     simpulse(2*pw,2*pwxC13,t4,zero,0.0,0.0);
     if(d2>0.0) delay(tau1-pwxC13);
    }
  else
    {
     if(d2>0.0) delay(tau1-pw/2.0);
     rgpulse(2.0*pw,t4,0.0,0.0);
     if(d2>0.0) delay(tau1-pw/2.0);
    }
  txphase(zero);

  delay(BigT);

  dec2rgpulse(2*pwxN15,t3,0.0,0.0);

  dec2phase(t5);

  delay(2.0e-6);
  zgradpulse(gzlvl3,gt3);
  delay(2.0e-6);
 
  if( (c180_flg[A] == 'y') && (pwxC13 > pw) )
     delay(BigT + 2*pwxC13 - gt3 - 2*GRADIENT_DELAY - 4.0e-6);  
  else 
     delay(BigT + 2.0*pw - gt3 - 2*GRADIENT_DELAY - 4.0e-6);  

  sim3pulse(pw,0.0e-6,pwxN15,zero,zero,t5,0.0,0.0);   /* X read pulse */

  txphase(zero);
  dec2phase(zero);

  delay(2.0e-6);
  zgradpulse(gzlvl5,gt5);
  delay(200.0e-6);

  delay(tauxh - gt5 - 202.0e-6);           /* delay=1/4J (XH)  */

  sim3pulse(2.0*pw,0.0e-6,2*pwxN15,zero,zero,zero,0.0,0.0);  

  txphase(one); dec2phase(one);

  delay(tauxh - gt5 - 202.0e-6);            /* delay=1/4J (XH)  */

  delay(2.0e-6);
  zgradpulse(gzlvl5,gt5);
  delay(200.0e-6);

  sim3pulse(pw,0.0e-6,pwxN15,one,zero,one,0.0,0.0);

  dec2phase(zero); txphase(zero);

  delay(2.0e-6);
  zgradpulse(gzlvl5,gt5);
  delay(200.0e-6);

  delay(tauxh - gt5 - 202.0e-6);            /* delay=1/4J (XH)  */

  sim3pulse(2.0*pw,0.0e-6,2*pwxN15,zero,zero,zero,0.0,0.0);

  dec2power(dpwr2);  /* decoupler power for decoupling on decouper channel 2 */
  decpower(dpwr);

  delay(tauxh - 2.0*POWER_DELAY - gt5 - 202.0e-6);  /* delay=1/4J(XH)   */

  delay(2.0e-6);
  zgradpulse(gzlvl5,gt5);
  delay(200.0e-6);

  rgpulse(pw,zero,0.0,0.0);

  delay(BigT1);

  rgpulse(2.0*pw,zero,0.0,0.0);

  delay(2.0e-6);
  zgradpulse(icosel*gzlvl4,gt4);
  delay(2.0e-6);

  delay(BigT1 - gt4 - 4.0e-6 - 2.0*GRADIENT_DELAY);

/* acquire data */

status(D);
     setreceiver(t6);
}
