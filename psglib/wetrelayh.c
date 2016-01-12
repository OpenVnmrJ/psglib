/*  wetrelayh.c - made from relayh.c (9.2 9/29/94) - 
relayed cosy, including regular cosy and delayed cosy

   Parameters:

         tau = delay to allow the relay of magnetization via
	       scalar coupling
       relay = the number of relays to be employed; if relay = 0,
	       a delayed COSY is performed for tau > 0
          nt = 2**(relay+2), e.g., nt = 8 for relay = 1

  f. sauriol july 1986; s. patt june/july 1987
  P.A.Keifer 950406 - made wetrelayh.c
  P.A.Keifer 950920 - updated wet
  P.A.Keifer 960116 - added tpwrf control to wet4


  */


#include <standard.h>

pulsesequence()
{
   int          i,
		relay;

/* GET NEW PARAMETERS */
   relay = (int) (getval("relay") + 0.5);

/* CHECK CONDITIONS */
   if (p1 == 0.0)
      p1 = pw;


/* CALCULATE PHASES */
   sub(ct, ssctr, v11);		     /* v11 = ct-ss */
   initval(256.0, v12);
   add(v11, v12, v11);		     /* v11 = ct-ss+256 */
   hlv(v11, v1);		     /* v1 = cyclops = 00112233 */
   for (i = 0; i < relay + 1; i++)
      hlv(v1, v1);		     /* cyclops = 2**(relay+2) 0's, 1's, 2's,
					3's */

   mod4(v11, v2);		     /* v2 = 0123 0123... */
   dbl(v2, oph);		     /* oph = 0202 0202 */
   add(v2, v1, v2);		     /* v2 = 0123 0123 + cyclops */
   add(oph, v1, oph);		     /* oph = 0202 0202 + cyclops */
   hlv(v11, v3);		     /* v3 = 0011 2233 4455 ... */


/* PULSE SEQUENCE */
   status(A);
      if (rof1 == 0.0)
      {
         rof1 = 1.0e-6;		     /* phase switching time */
         rcvroff();
      }
      hsdelay(d1);		     /* preparation period */
     if (getflag("wet")) wet4(zero,one);

   status(B);
      rgpulse(pw, v1, rof1, rof1);
      delay(d2);		     /* evolution period */
   status(A);
      if (relay == 0)
         delay(tau);		     /* for delayed cosy */
      rgpulse(p1, v2, rof1, rof1);   /* start of mixing period */
      if (relay == 0)
         delay(tau);		     /* for delayed cosy */

      for (i = 0; i < relay; i++)    /* relay coherence */
      {
         hlv(v3, v3);		     /* v3=2**(relay+1) 0's, 1's, 2's, 3's */
         dbl(v3, v4);		     /* v4=2**(relay+1) 0's, 2's, 0's, 2's */
         add(v2, v4, v5);	     /* v5=v4+v2 (including cyclops) */
         delay(tau/2);
         rgpulse(2.0*pw, v2, rof1, rof1);
         delay(tau/2);
         rgpulse(pw, v5, rof1, rof1);
      }

   status(C);
      delay(rof2);
      rcvron();
}

/* wet4 - Water Elimination */
wet4(phaseA,phaseB)
  codeint phaseA,phaseB;
 
{
  double finepwr,gzlvlw,gtw,gswet,dmfwet,dpwrwet,dofwet,wetpwr,pwwet,dz;
  int c13wet;
  char   wetshape[MAXSTR];
  c13wet=getflag("c13wet");             /* Water suppression flag        */  
  getstr("wetshape",wetshape);    /* Selective pulse shape (base)  */
  wetpwr=getval("wetpwr");        /* User enters power for 90 deg. */
  pwwet=getval("pwwet");        /* User enters power for 90 deg. */
  dmfwet=getval("dmfwet");
  dpwrwet=getval("dpwrwet");
  dofwet=getval("dofwet");
  dz=getval("dz");
  finepwr=wetpwr-(int)wetpwr;     /* Adjust power to 152 deg. pulse*/
  wetpwr=(double)((int)wetpwr);
  if (finepwr==0.0) {wetpwr=wetpwr+5; finepwr=4095.0; }
  else {wetpwr=wetpwr+6; finepwr=4095.0*(1-((1.0-finepwr)*0.12)); }
  rcvroff();
  if (c13wet)
    {
    setstatus(DECch,FALSE,'w',FALSE,dmfwet);
    decoffset(dofwet);
    decpower(dpwrwet);
    }
  obspower(wetpwr);         /* Set to low power level        */
  gzlvlw=getval("gzlvlw");      /* Z-Gradient level              */
  gtw=getval("gtw");            /* Z-Gradient duration           */
  gswet=getval("gswet");        /* Post-gradient stability delay */
  chess(finepwr*0.5059,wetshape,pwwet,phaseA,20.0e-6,rof2,gzlvlw,gtw,gswet,c13wet);
  chess(finepwr*0.6298,wetshape,pwwet,phaseB,20.0e-6,rof2,gzlvlw/2.0,gtw,gswet,c13wet);
  chess(finepwr*0.4304,wetshape,pwwet,phaseB,20.0e-6,rof2,gzlvlw/4.0,gtw,gswet,c13wet);
  chess(finepwr*1.00,wetshape,pwwet,phaseB,20.0e-6,rof2,gzlvlw/8.0,gtw,gswet,c13wet);
  if (c13wet)
    {
    setstatus(DECch,FALSE,'c',FALSE,dmf);
    decoffset(dof);
    decpower(dpwr);
    }
  obspower(tpwr);               /* Reset to normal power level   */
  obspwrf(tpwrf);
  rcvron();
  delay(dz);
}
/* chess - CHEmical Shift Selective Suppression */
chess(pulsepower,pulseshape,duration,phase,rx1,rx2,gzlvlw,gtw,gswet,c13wet)  double pulsepower,duration,rx1,rx2,gzlvlw,gtw,gswet;
  int c13wet;
  codeint phase;
  char* pulseshape;
{
  obspwrf(pulsepower);
  if (c13wet) decon();
  shaped_pulse(pulseshape,duration,phase,rx1,rx2);
  if (c13wet) decoff();
  zgradpulse(gzlvlw,gtw);

  delay(gswet);
}
 
int getflag(str)
char str[MAXSTR];
{
   char strval[MAXSTR];
 
   getstr(str,strval);
   if ((strval[0]=='y') || (strval[0]=='Y')) return(TRUE);
     else                                    return(FALSE);
}


