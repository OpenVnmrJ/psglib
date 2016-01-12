/* wetgcosy.c - made from gcosy (9.1 9/28/93)
   gcosy.c - pulsed field gradient COSY (magnitude mode)

	gzlvl1 = gradient amplitude (-32768 to +32768)
	gt1 = gradient duration in seconds (0.001)
	grise = gradient rise and fall time (in seconds; 0.00001)
	qlvl = 1 (quantum selection level)
	gstab = optional delay for stability
	phase = 1 (selects echo N-type coherence selection; default)
              = 2 (selects antiecho P-type coherence selection) 

        process N-type data with wft2d(1,0,0,1)
	process P-type data with wft2d(1,0,0,-1)
                   the ('t2dc') arguement to wft2d may be useful 

	PAK 920420 - 
	PAK 920429 - gcosy.c from pgcosy.c
	PAK 920514 - corrected and verified P vs. N-type aspects
	P.A.Keifer 950406 - made wetgcosy.c
	P.A.Keifer 950920 - updated wet
        P.A.Keifer 960116 - added tpwrf control to wet4
 */


#include <standard.h>

pulsesequence()
{
  double gzlvl1,gt1,grise,qlvl,gstab,phase;
  int icosel;

  gzlvl1 = getval("gzlvl1");
  gt1 = getval("gt1");
  grise = getval("grise");
  qlvl = getval("qlvl");
  gstab = getval("gstab");
  phase = getval("phase");


  if (phase == 2.0) 
    { 
     icosel=-1;
     if (ix==1) 
     fprintf(stdout,"P-type COSY\n");
    }
  else
    { 
     icosel=1;    /* Default to N-type experiment */
     if (ix==1) 
     fprintf(stdout,"N-type COSY\n");
    }

/* BEGIN ACTUAL PULSE SEQUENCE */
  status(A);
     delay(d1);
   if (getflag("wet")) wet4(zero,one);

  status(B);
     rgpulse(pw,oph,rof1,rof2);
     delay(d2);
     rgradient('z',gzlvl1);
     delay(gt1+grise);
     rgradient('z',0.0);
     delay(grise);

     rgpulse(pw,oph,rof1,rof2);
     rgradient('z',gzlvl1*icosel);
     delay(qlvl*(gt1+grise));
     rgradient('z',0.0);
     delay(grise);
     delay(gstab);

  status(C);
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


