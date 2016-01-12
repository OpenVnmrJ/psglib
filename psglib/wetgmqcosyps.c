/* wetgmqcosyps - made from gmqcosyps (from gmqcosy V1.1 10/28/92)
	phase-sensitive pulsed gradient enhanced mq cosy
	(currently (930326) this is set up with the gmqcosy macro)
   
   Parameters:
	glvl - set to 2 for a DQFCOSY
	gzlvl1 - try 15000 (must be < 32767/qlvl)
	gt1 - try 0.0025 seconds for organics;
		at least 0.008-0.012 sec for proteins
	phase - 1,2 (run and process as a standard hypercomplex dataset)
	pwwet
	wetpwr
	wetshape

   Processing:
	if phase=1,2: use wft2da
		(after running the gmqcosy macro to setup, set rp1=90
		and rp=rp+90 as compared to the phase 1D spectrum
		(the gmqcosy macro may do part of this))
	 
	pak 920506 -
	eh 9301 - tested and working (was ehgmqcosy)
	pak 930326 - added more documentation
	P.A.Keifer 930806 - added rcvr control
	SHS & PAK 9502 - made wetgmqcosyps.c
	P.A.Keifer 950407 - standardized
	P.A.Keifer 950920 - updated wet
        P.A.Keifer 960116 - added tpwrf control to wet4
*/


#include <standard.h>

pulsesequence()
{
        double gzlvl1,qlvl,grise,gstab,gt1,ss,phase;
	int icosel,iphase;
  	ss = getval("ss");
        grise=getval("grise");
        gstab=getval("gstab"); 
        gt1=getval("gt1");
        gzlvl1=getval("gzlvl1");
        qlvl=getval("qlvl");
        phase=getval("phase");
        iphase = (int) (phase + 0.5);



/* DETERMINE STEADY-STATE MODE */
   if (ss < 0) ss = (-1) * ss;
   else
      if ((ss > 0) && (ix == 1)) ss = ss;
      else
         ss = 0;
   initval(ss, ssctr);
   initval(ss, ssval);

   assign(oph,v2);
   assign(oph,v1);
   if (iphase == 2)
      incr(v1);

/* HYPERCOMPLEX MODE */
   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v6);
   if ((iphase==1) || (iphase==2))
      {add(v1,v6,v1); add(oph,v6,oph);} 
	

     status(A);
        rcvroff();
	delay(d1);
     if (getflag("wet")) wet4(zero,one);
	rlpower(tpwr,TODEV);

	rgpulse(pw,v1,rof1,rof2);
     status(C);
        if (d2 > rof1 + 4.0*pw/3.1416)
           delay(d2 - rof1 - 4.0*pw/3.1416);
     status(C);
        rgpulse(pw,v2,rof1,rof2);
      
        delay(gt1+2.0*grise+24.4e-6);
        rgpulse(2.0*pw,v2,rof1,rof2);
        icosel=-1;
        rgradient('z',gzlvl1*(double)icosel);
	delay(gt1+grise);
	rgradient('z',0.0);
        txphase(oph);
	delay(grise);

	rgpulse(pw,v2,rof1,rof2);

        rgradient('z',gzlvl1*qlvl);

	delay(gt1+grise);
	rgradient('z',0.0);
	delay(grise);
       rcvron();
	delay(gstab);  
     status(D);
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

 

