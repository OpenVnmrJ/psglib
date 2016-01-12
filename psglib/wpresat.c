 
/*  wpresat-
               d1...satdly...p1...d2(sat. opt.)...pw....at......

	two-pulse sequence with WFG control of satdly(s) and pulses 

   satfrq = frequency for saturation
    p1frq = frequency for p1 pulse
      tof = frequency for pw pulse
    p1pwr = power level for p1 pulse
    pwpwr = power level for pw pulse
   satpwr = power level for saturation pulse(s)
  p1shape = 'name'   (file found in shapelib)
  pwshape = 'name'   (file found in shapelib)
 satshape = 'name'   (file found in shapelib) 
 satmode  = determines when saturation occurs. Use as "dm" in s2pul

         G.Gray Palo Alto, Sept 1991     */

#include <standard.h>
pulsesequence()
{
   double satdly,satpwr,p1frq,satfrq,p1pwr,pwpwr;
   char satmode[MAXSTR], satshape[MAXSTR], p1shape[MAXSTR], pwshape[MAXSTR];
     satdly = getval("satdly");
     satfrq = getval("satfrq");
     satpwr = getval("satpwr"); initval(satpwr,v7);
      p1frq = getval("p1frq");
      p1pwr = getval("p1pwr");  initval(p1pwr,v8);
      pwpwr = getval("pwpwr");  initval(pwpwr,v9);
     getstr("p1shape",p1shape);
     getstr("pwshape",pwshape);
     getstr("satshape",satshape);
     getstr("satmode",satmode);

    /* equilibrium period */
    status(A);
     if (satmode[A]=='y') 
     {          
       if (fabs(satfrq - tof)>0.1)  offset(satfrq,TODEV);
       hsdelay(d1);
       power(v7,TODEV); 
                             dbl(ct,v1);      /*v1=02020202	*/
                             hlv(ct,v2);      /* v2=00112233	*/
                             mod2(v2,v3);     /* v3=00110011	*/
                             add(v1,v3,v4);   /* v4=02130213	*/
                             add(v2,one,v5);  /* v5=11223300	*/
                             mod4(v4,oph);    /* oph=02130213	*/
       shaped_pulse(satshape,satdly,v5,rof1,rof2); 
     }
     else
       hsdelay(d1);
    status(B);
     if (p1>0.0)
      {
       if (fabs(satfrq - p1frq)>0.1)  offset(p1frq,TODEV);
       power(v8,TODEV); 
       shaped_pulse(p1shape,p1,zero,20.0e-6,rof2);
       if (satmode[B]=='y')
        {
         if (fabs(satfrq-p1frq)>0.0) offset(satfrq,TODEV); 
         power(v7,TODEV); shaped_pulse(satshape,d2,v5,rof1,rof2);
         if (fabs(satfrq-tof)>0.0) offset(tof,TODEV);
         power(v9,TODEV);
        }
       else
        {
         if (fabs(p1frq-tof)>0.0) offset(tof,TODEV); 
         power(v9,TODEV); hsdelay(d2);
        }
       }
      else
       {
        if (fabs(satfrq-tof)>0.1) offset(tof,TODEV); power(v9,TODEV);
       }
     status(C);
       shaped_pulse(pwshape,pw,oph,20.0e-6,rof2);
}
