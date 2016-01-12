/* heteronuclear spin echo  non-c13 proton suppression */
/*    requires linear amplifiers                       */
/*   contact- G. Gray (palo alto)  revision-           */

#include <standard.h>

pulsesequence()
{
   double ss,j,pwx,pwxlvl,satpwr,satfrq,satdly; tpwr=getval("tpwr");
   satpwr=getval("satpwr"); satfrq=getval("satfrq"); satdly=getval("satdly");
   pwx=getval("pwx"); pwxlvl=getval("pwxlvl"); initval(pwxlvl,v10); initval(dpwr,v11);
   initval(satpwr,v8); initval(tpwr,v9); j=getval("j");
   
   ss=getval("ss");
   if (ss<0)
    {
       ss = (-1)*ss;
       initval(ss,ssval);
       initval(ss,ssctr);
     }

   mod2(ct,v2); dbl(v2,v2);               /* v2=0202020202020202 */
   assign(zero,v1);                       /* v1=0000000000000000 */
   hlv(ct,v3);                            /* v3=0011223344556677 */
   add(v3,v2,v2); add(v3,v1,v1);            /* add cyclops  */
   assign(v2,oph);  add(v1,one,v4);
   status(A);
     power(v8,TODEV);
      power(v10,DODEV);
      offset(satfrq,TODEV);
      hsdelay(d1);  
      rgpulse(satdly,zero,rof1,0.0);
      power(v9,TODEV); 
      offset(tof,TODEV);
   status(B);
      rgpulse(pw,v1,rof1+18.0e-6,10.0e-6); 
      delay(1.0/(2.0*j));  
      decpulse(pwx,v2);
      rgpulse(2.0*pw,v4,0.0,0.0);
      delay(10.0e-6);
      decpulse(pwx,v1);
      power(v11,DODEV);
      delay(1.0/(2.0*j));
   status(C);
      txphase(zero); decphase(zero);
}
