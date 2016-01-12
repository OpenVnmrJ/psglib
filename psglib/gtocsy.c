/* gtocsy     
	2D tocsy 
       Uses pfg to help suppress water
       Use about a 30 us H 90 for the dipsi portion of the sequence
       Use about a 10us H 90 for the hard proton pulses
           includes the super scuba sequence of Brown, Weber and Mueller.  Set
           it equal to 'n' to do standard presaturation.  Note: presaturation
           is done using the transmitter with the power level set by 'satpwr'

       Written by L.E. Kay 10-15-92
       Modified by G.Gray 5-13-93

*/

#include <standard.h>
#define PI 3.1416

static double d2_init = 0.0;

static int ph1[8] = {2,2,2,2,0,0,0,0},
           ph2[16] = {0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2},
           ph3[4] = {0,0,2,2}, 
           ph4[16]  = {1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3},
           ph5[4] = {0,1,2,3},
           ph6[1] = {0},
           rec[16] = {1,0,3,2,3,2,1,0,3,2,1,0,1,0,3,2};

pulsesequence()
{
/* DECLARE VARIABLES */

 char        satmode[MAXSTR], sspul[MAXSTR],
             Sshape[MAXSTR],SSshape[MAXSTR]; 

 int	     phase,t1_counter;

 double      scuba,                /* length of 1/2 scuba delay */
             satpwr,               /* low power level for presat*/
             satfrq,               /* saturation frequency  */
             satdly,               /* pre saturation time */
             p1lvl,                /* power level for hard pulses */
             trim,                  /* delay for trims on either side of dipsi
                                      use trim = 1.0e-3  */
	     sw1,                  /* spectral width in 1H dimension */
             dipsitime,             /* total length of dipsi sequence */
             mix,                  /* dipsi mix time: pw90=pw at tpwr */
             SSpw,             /* PW90 for ss shape pulse  */
             Spw,              /* PW180 for s shape pulse  */
             SSpwr,             /* Power level for ss pulse */
             Spwr,              /* Power level for s pulse */
             p10,                   /* 320 proton pulse width at tpwr */
             p11,                   /* 410 proton pulse width at tpwr */
             p12,                   /* 290 proton pulse width at tpwr */
             p13,                   /* 285 proton pulse width at tpwr */
             p14,                   /* 30 "" */
             p15,                   /* 245 proton pulse width at tpwr */
             p16,                   /* 375 */ 
             p17,                   /* 265 */
             p18,                   /* 370 */
             p19,                   /* 60 */ 
             delay_clean,          /* delay = 0.5*dipsi mix to get "clean"
                                      spectra            */
             ncyc,                 /* number of cycles of dipsi sequence
                                      a typical number is 10 */
             gt1,
             gzlvl1,
             gt2,
             gzlvl2,
             BigT;


/* LOAD VARIABLES */

  SSpw = getval("SSpw");
  Spw = getval("Spw");
  SSpwr = getval("SSpwr");
  Spwr = getval("Spwr");
  gt1 = getval("gt1");
  gzlvl1 = getval("gzlvl1");
  gt2 = getval("gt2");
  gzlvl2 = getval("gzlvl2");
  BigT = getval("BigT");
   mix = getval("mix");
  satdly = getval("satdly");
  pw = getval("pw");
  p1lvl = getval("p1lvl");
  satpwr = getval("satpwr");
  satfrq = getval("satfrq");
  scuba = getval("scuba");
  trim = getval("trim");
  phase = (int) (getval("phase") + 0.5);
  sw1 = getval("sw1");

  getstr("sspul",sspul);
  getstr("satmode",satmode); 
  getstr("Sshape",Sshape);
  getstr("SSshape",SSshape);

/* check validity of parameter range */

    if((dm[A] == 'y' )||( dm[C] == 'y'))
	{
	printf("incorrect Dec1 decoupler flags!  ");
	abort(1);
    } 

    if((dm2[A] == 'y')||(dm2[C] == 'y'))
	{
	printf("incorrect Dec2 decoupler flags!  ");
	abort(1);
    } 

    if( satpwr > 8 )
    {
	printf("satpwr too large !!!  ");
	abort(1);
    }

    if( tpwr > 57 )
    {
	printf("tpwr too large !!!  ");
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

    if( mix > 0.1 )
    {
	printf("don't fry the probe, mix too large !  ");
	abort(1);
    }

    if ( trim > 0.004 )
    {
        printf("dont fry the probe, trim is too long ! " );
        abort(1);
     }


/* LOAD VARIABLES */

  settable(t1, 8,  ph1);
  settable(t2, 16,  ph2);
  settable(t3, 4, ph3);
  settable(t4, 16,  ph4);
  settable(t5, 4,  ph5);
  settable(t6, 1, ph6);
  settable(t7,16,rec);

/* INITIALIZE VARIABLES */
  p10 = (320.0/90.0)*pw;
  p11 = (410.0/90.0)*pw;
  p12 = (290.0/90.0)*pw;
  p13 = (285.0/90.0)*pw;
  p14 = (30.0/90.0)*pw;
  p15 = (245.0/90.0)*pw;
  p16 = (375.0/90.0)*pw;
  p17 = (265.0/90.0)*pw;
  p18 = (370.0/90.0)*pw;
  p19 = (60.0/90.0)*pw;

  dipsitime = (4.0*(p10+p11+p12+p13+p14+p15+p16+p17+p18) +p19);
  ncyc=mix/dipsitime;
  ncyc=2.0*(double)(int)(ncyc/2.0);
    if( ncyc > 20 )
    {
	printf("don't fry the probe, mix too large !  ");
	abort(1);
    }
  initval(ncyc,v8);

  delay_clean = ncyc*pw*115.778;
  delay_clean = delay_clean/2.0;  /* divide by two since roe is -2X the rate 
                                     of noe  */
  

/* Phase incrementation for hypercomplex data */

   if ( phase == 2 )     /* Hypercomplex in t1 */
        tsadd(t1, 1, 4);

/* calculate modifications to phases bsed on current t1 vaules to 
   achieve States-TPPI acquisition */

   if( ix == 1)
      d2_init = d2;
   t1_counter = (int) ( (d2-d2_init)*sw1 + 0.5);
   if(t1_counter % 2) {
      tsadd(t1,2,4);
      tsadd(t7,2,4);
   } 


/* BEGIN ACTUAL PULSE SEQUENCE */

status(A);
   rlpower(p1lvl,TODEV);

/* kill any remaining steady state magnetization */
 if (sspul[A]=='y')
  {
   rgpulse(11.2*p1,zero,0.0,0.0);
   rgpulse(16.4*p1,one,2.0e-6,0.0);
  }
/* Presaturation Period */

 if(satmode[A] == 'y')
{
  if (fabs((satfrq-tof)>0.0)) offset(satfrq,TODEV);
  rlpower(satpwr,TODEV);
  rgpulse(satdly,zero,rof1,0.0);
  if (fabs((satfrq-tof)>0.0)) offset(tof,TODEV);
  rlpower(p1lvl,TODEV);                /* Set power for hard pulses  */

    if (scuba>0)                       /* Scuba pulse sequence */
    {  
      hsdelay(scuba);

      rgpulse(p1,zero,1.0e-6,0.0);	/* 90x180y90x */
      rgpulse(2*p1,one,1.0e-6,0.0);
      rgpulse(p1,zero,1.0e-6,0.0);
 
      txphase(zero);
      delay(scuba);        
    }
}
 delay(d1);
 rcvroff();
status(B);

  initval(1.0,v2);
  stepsize(45.0,TODEV);
  xmtrphase(v2);   

  rgpulse(p1,t1,2.0e-6,0.0);        /* Proton excitation pulse */
  xmtrphase(zero);  

  txphase(t3);
  dec2phase(t2);
  if (d2>0)
    delay(d2 - 4.0*p1/PI);
  rgpulse(p1,t3,2.0e-6,0.0); 
 status(C);
  rlpower(tpwr,TODEV);              /* dipsi power level */
  txphase(zero);
  rgpulse(trim,zero,0.0,0.0);
  
  if(mix > 0.0) {

     starthardloop(v8);

      rgpulse(p10,one,0.0,0.0);      
      rgpulse(p11,three,0.0,0.0);      
      rgpulse(p12,one,0.0,0.0);      
      rgpulse(p13,three,0.0,0.0);      
      rgpulse(p14,one,0.0,0.0);      
      rgpulse(p15,three,0.0,0.0);      
      rgpulse(p16,one,0.0,0.0);      
      rgpulse(p17,three,0.0,0.0);      
      rgpulse(p18,one,0.0,0.0);      

      rgpulse(p10,three,0.0,0.0);      
      rgpulse(p11,one,0.0,0.0);      
      rgpulse(p12,three,0.0,0.0);      
      rgpulse(p13,one,0.0,0.0);      
      rgpulse(p14,three,0.0,0.0);      
      rgpulse(p15,one,0.0,0.0);      
      rgpulse(p16,three,0.0,0.0);      
      rgpulse(p17,one,0.0,0.0);      
      rgpulse(p18,three,0.0,0.0);      

      rgpulse(p10,three,0.0,0.0);      
      rgpulse(p11,one,0.0,0.0);      
      rgpulse(p12,three,0.0,0.0);      
      rgpulse(p13,one,0.0,0.0);      
      rgpulse(p14,three,0.0,0.0);      
      rgpulse(p15,one,0.0,0.0);      
      rgpulse(p16,three,0.0,0.0);      
      rgpulse(p17,one,0.0,0.0);      
      rgpulse(p18,three,0.0,0.0);      

      rgpulse(p10,one,0.0,0.0);      
      rgpulse(p11,three,0.0,0.0);      
      rgpulse(p12,one,0.0,0.0);      
      rgpulse(p13,three,0.0,0.0);      
      rgpulse(p14,one,0.0,0.0);      
      rgpulse(p15,three,0.0,0.0);      
      rgpulse(p16,one,0.0,0.0);      
      rgpulse(p17,three,0.0,0.0);      
      rgpulse(p18,one,0.0,0.0);      

      rgpulse(p19,zero,0.0,0.0);

     endhardloop();
  } 
  txphase(zero);
  rgpulse(trim,zero,0.0,0.0);

  rlpower(p1lvl,TODEV);                /*  power for  high power pulses */

  /* flip back pulse followed by period to get clean magnetization */

  rgpulse(p1,t4,2.0e-6,2.0e-6);

  txphase(t5);
  rgradient('z',gzlvl1);
  delay(gt1);
  rgradient('z',0.0);
  delay(delay_clean - gt1);
  rlpower(SSpwr,TODEV);
  shaped_pulse(SSshape,SSpw,t5,2.0e-6,2.0e-6);
  
  rgradient('z',gzlvl2);
  delay(gt2);
  rgradient('z',0.0);

  delay(BigT - gt2 - POWER_DELAY);

  rlpower(Spwr,TODEV);
  shaped_pulse(Sshape,Spw,zero,2.0e-6,2.0e-6);

  rgradient('z',gzlvl2);
  delay(gt2);
  rgradient('z',0.0);

  delay(BigT - gt2 + SSpw/2.0);

status(D);
  setreceiver(t7);
  rcvron();
}
