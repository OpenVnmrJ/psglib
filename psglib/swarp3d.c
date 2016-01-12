/* 
 *	A multi-slice3D spin warp experiment with optional inversion recovery
 */
#include <standard.h>
#include <stdio.h>

static   double gro,grof,gpe1,gpe2,xinc,yinc,te,trise,ti,tr,dtr;
static   double tpe,ni,at,dx,gpeoffset,gpex1,gpex2,dy,dz;
static   double tpwr1,tpwri,pi,sw1,sw2,ni2,d3;
static   char gread,gphase,gslice,ir[MAXSTR];
static   char p1pat[MAXSTR],pwpat[MAXSTR],pipat[MAXSTR];

pulsesequence()
{
/*  


Authors:
 Evan Williams
*/


   get_parameters();
   status(A);
/* 
   program loop 
*/
  
    hsdelay(dtr);
    status(B);
    if (ir[0] == 'y')
     inversion_pulse();
    selective_90();
    phase_n_compensate();
    delay(dx);
    selective_180();
    delay(dy);
    read_out();
}

get_parameters()
{ 
   if (getorientation(&gread,&gphase,&gslice,"orient") < 0) 
     exit(1);
   gro = getval("gro");
   grof = getval("grof");
   gpe1 = getval("gpe1");
   gpe2 = getval("gpe2");	/*note - what was slice select is now 2nd pe*/
   at   = getval("at");
   d3   = getval("d3");
   sw1   = getval("sw1");
   sw2   = getval("sw2");
   ni =  (int) getval("ni");  
   ni2 =  (int) getval("ni2");  
   trise = getval("trise");
   te = getval("te");
   ti = getval("ti");
   tr = getval("tr");
   tpe = getval("tpe");
   tpwr1 = getval("tpwr1");
   tpwri = getval("tpwri");
   pi = getval("pi");
   getstr("p1pat",p1pat);
   getstr("pwpat",pwpat);
   getstr("pipat",pipat);
   getstr("ir",ir);

/* now we need to find out what increment we are at*/

    xinc=(d2*sw1)+1.0;
    yinc=(d3*sw2)+1.0;

/*
   do some error checking on string variables 
*/
   if (p1pat[0] == '\0') 
   { 
     text_error("p1pat file spec error\n"); 
     exit(2); 
   }
   if (pwpat[0] == '\0') 
   { 
     text_error("pwpat file spec error\n"); 
     exit(2); 
   }
   if ((ir[0] == 'y') && (pipat[0] == '\0')) 
   { 
     text_error("pipat file spec error\n"); 
     exit(2); 
   }
/*  */
/* calculate dx: te = p1/2 + tpe + dx + pw + dy + aq/2 + 3*trise */
/*  */

   if (tpe >= p1)
     dx = (te/2) - tpe - trise - 0.5*(p1+pw);
     else
     dx = (te/2) - p1 - trise - 0.5*(p1+pw);
   if (dx < 1.0e-6) 
   { 
     text_error("cannot fit tpe trise p1 pw into te... ABORT");
     exit(2);
   }
   dy = (te/2) - at/2 - trise - pw/2;
   if (dy < 1.0e-6) 
     dy = 1.0e-6;
/* 
 *   calculate dz: ti = 0.5*(p1+pi) + 2*trise + dz
 */
   dz = ti - 0.5*(p1+pi) - trise;
   if (dz < 1.0e-6) 
     dz=1.0e-6;
/* 
   calculate dtr: tr = dtr + nfi*(ti + te + at/2 + trise + pi/2)
*/
   if (ir[0] =='y')
     dtr = tr - (ti+te+at/2+trise+pi/2);
   else
     dtr = tr - (te+at/2+trise+p1/2);

   if (dtr < 1.0e-6)
     dtr=1.0e-6;
  
   if (ni > 1.0 ) 
   { 				/* use gpe1 as step size with ni */
     gpeoffset = ni * gpe1/ 2.0; /* commented out so to collect 1/2 echo */
     gpex1 = xinc*gpe1 - gpeoffset;
   }
   else 
     gpex1 = gpe1;

   if (ni2 > 1.0 ) 
   { 				/* use gpe2 as step size with ni2 */
/*     gpeoffset = ni2 * gpe2/ 2.0; *//* 1/2 echo again */
     gpex2 = yinc*gpe2/* - gpeoffset*/;
   }
   else 
      gpex2 = gpe2;
/* set the powers */
   initval(tpwr,v1);
   initval(tpwr1,v2);
   initval(tpwri,v4);
   assign(oph,v3);
/* set the phase cycle */
   if (((int) xinc % 2) == 0) 
      add(one,v3,v3); 
   else 
      add(three,v3,v3);
}



inversion_pulse()
{
  power(v4,TODEV);
  shaped_pulse(pipat,pi,oph,rof1,rof2);
  delay(dz);                          /* inversion period ends */
}

selective_90()
{
   power(v2,TODEV);
   shaped_pulse(p1pat,p1,oph,rof1,rof2);
   rgradient(gread,grof*gro);
   rgradient(gphase,gpex1);
   rgradient(gslice,gpex2);
}

phase_n_compensate()
{
  if (tpe >= p1)
  { 				/* calculate delays for all gradients  */
    delay(p1);
    if ((tpe-p1) == 0.0)
    {
      rgradient(gread,0.0); 
      rgradient(gphase,0.0);
      rgradient(gslice,0.0);
   }
    else
    {
      delay(tpe-p1);
      rgradient(gread,0.0);
      rgradient(gphase,0.0);
      rgradient(gslice,0.0);
    }
  }
  else
  {
    delay(tpe);
    rgradient(gread,0.0);
    rgradient(gphase,0.0);
    rgradient(gslice,0.0);
    delay(p1-tpe);
  }
}

selective_180()
{
   power(v1,TODEV);
   status(C);
   shaped_pulse(pwpat,pw,v3,rof1,rof2);
   offset(tof,TODEV);		/*  go back to observe position */
}

read_out()
{
  rgradient(gread,gro);
  delay(trise);
  acquire(np,1.0/sw);		/* Acquire  */
  rgradient(gread,0.0);
}
