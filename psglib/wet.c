/* wet */

#include <standard.h>
extern double getvalnwarn();

pulsesequence()
{
  /* DECLARE & READ IN NEW PARAMETERS */

  char p1shape[MAXSTR];
  double p1pwr,dz;
  int composite, wet;

  loadtable("wet");              /* Phase table                   */
  composite=getflag("composit");  /* Composite observe pulse flag  */
  wet=getflag("wet");             /* Water suppression flag        */

  if (wet) {                      /* Get WET parameters            */
    getstr("p1shape",p1shape);    /* Selective pulse shape (base)  */
    p1pwr=getval("p1pwr");        /* User enters power for 90 deg. */
    dz=getval("dz");              /* delay to allow return to zero */
   }


  /* PULSE SEQUENCE */

  status(A);

    hsdelay(d1);

  status(B);

    if (wet) { wet4(p1pwr,p1shape,p1,t1,t2); delay(dz); }

  status(C); 

    if (composite) composite_pulse(pw,t3,rof1,rof2,v1);
      else rgpulse(pw,t3,rof1,rof2);
    setreceiver(t4);
}





/* wet4 - Water Elimination */
wet4(pulsepower,shape,duration,phaseA,phaseB)
  double pulsepower,duration;
  codeint phaseA,phaseB;
  char* shape;
{
  double finepwr,gzlvlw,gtw,gstab;
  finepwr=pulsepower-(int)pulsepower;     /* Adjust power to 152 deg. pulse*/
  pulsepower=(double)((int)pulsepower);             
  if (finepwr==0.0) {pulsepower=pulsepower+5; finepwr=4095.0; }
  else {pulsepower=pulsepower+6; finepwr=4095.0*(1-((1.0-finepwr)*0.12)); }
  rcvroff();
  obspower(pulsepower);         /* Set to low power level        */
  gzlvlw=getval("gzlvlw");      /* Z-Gradient level              */
  gtw=getval("gtw");            /* Z-Gradient duration           */
  gstab=getval("gstab");        /* Post-gradient stability delay */
  chess(finepwr*0.6452,shape,p1,phaseA,20.0e-6,rof2,gzlvlw,gtw,gstab);
  chess(finepwr*0.5256,shape,p1,phaseB,20.0e-6,rof2,gzlvlw/2.0,gtw,gstab);
  chess(finepwr*0.4928,shape,p1,phaseB,20.0e-6,rof2,gzlvlw/4.0,gtw,gstab);
  chess(finepwr*1.00,shape,p1,phaseB,20.0e-6,rof2,gzlvlw/8.0,gtw,gstab);
  obspower(tpwr); obspwrf(tpwrf);  /* Reset to normal power level   */
  rcvron();
}

/* chess - CHEmical Shift Selective Suppression */
chess(pulsepower,pulseshape,duration,phase,rx1,rx2,gzlvlw,gtw,gstab)
  double pulsepower,duration,rx1,rx2,gzlvlw,gtw,gstab;
  codeint phase;
  char* pulseshape;
{
  status(B);
  obspwrf(pulsepower);
  shaped_pulse(pulseshape,duration,phase,rx1,rx2);
  status(C);
  zgradpulse(gzlvlw,gtw);                                                   
  delay(gstab);
}

composite_pulse(width,phasetable,rx1,rx2,phase)
  double width,rx1,rx2;
  codeint phasetable,phase;
{
  getelem(phasetable,ct,phase); /* Extract observe phase from table */
  incr(phase); rgpulse(width,phase,rx1,rx1);  /*  Y  */
  incr(phase); rgpulse(width,phase,rx1,rx1);  /* -X  */
  incr(phase); rgpulse(width,phase,rx1,rx1);  /* -Y  */
  incr(phase); rgpulse(width,phase,rx1,rx2);  /*  X  */
}

int getflag(str)
char str[MAXSTR];
{
   char strval[MAXSTR];

   getstr(str,strval);
   if ((strval[0]=='y') || (strval[0]=='Y')) return(TRUE);
     else                                    return(FALSE);
}

