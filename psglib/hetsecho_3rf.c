/* heteronuclear spin echo using second decoupler   */

#include <standard.h>

pulsesequence()
{
   double ss,j,pwx2,pwx2lvl,satpwr,satfrq,satdly; 
   satpwr=getval("satpwr"); satfrq=getval("satfrq"); satdly=getval("satdly");
   pwx2=getval("pwx2"); pwx2lvl=getval("pwx2lvl"); 
   j=getval("j");
   
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
      rlpower(satpwr,TODEV);
      rlpower(pwx2lvl,DO2DEV);
      if ((satfrq-tof)!=0.0) offset(satfrq,TODEV);
      hsdelay(d1);  
      rcvroff();
      pulse(satdly,zero);
      rlpower(tpwr,TODEV); 
      if ((satfrq-tof)!=0.0) offset(tof,TODEV);
   status(B);
      rgpulse(pw,v1,5.0e-6,0.0); 
      delay(1.0/(2.0*j));  
      dec2rgpulse(pwx2,v2,0.0,0.0);
      delay(d2/2);
      rgpulse(2.0*pw,v4,0.0,0.0);
      delay(d2/2);
      dec2rgpulse(pwx2,v1,0.0,0.0);
      rlpower(dpwr2,DO2DEV);
      delay(1.0/(2.0*j));
   status(C);
      rcvron();
}
