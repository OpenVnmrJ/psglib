/* dante -   selective excitation via pulse train  */
/*       contact- G. Gray  (palo alto)  revision-  */
#include <standard.h>

pulsesequence()
{
    double count,offset;
    offset = getval("offset");
    count = getval("count");
     d2=1/offset;
    initval(count,v1);
    /* equilibrium period */
     status(A);
     hsdelay(d1);
     if (count == 1.0) pulse(pw,oph);
     else
     {
     starthardloop(v1);
      rgpulse(pw,oph,rof1,rof2);
      delay(d2-rof2-rof1);
     endhardloop();
     }
}
