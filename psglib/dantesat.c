/*  dantesat-
              {trim(x)trim(y)} ......d1...satdly..pw....at......
        sspul = 'y' does trim(x)trim(y) to destroy all magnetization
               G.Gray, Palo Alto  Sept. 1991   version 4.1 vnmr
*/
#include <standard.h>
pulsesequence()
{
  double tau,offset,satdly,satpwr;
  char sspul[MAXSTR];
  satdly = getval("satdly");  satpwr = getval("satpwr");
  getstr("sspul",sspul); 
  offset = getval("offset");
       dbl(ct,v1);      /*v1=02020202	*/
       hlv(ct,v2);      /* v2=00112233	*/
       mod2(v2,v3);     /* v3=00110011	*/
       add(v1,v3,v4);   /* v4=02130213	*/
       add(v2,one,v5);  /* v5=11223300	*/
       mod4(v4,oph);    /* oph=02130213	*/

  initval(satdly/(1/offset),v13);
  tau=(1/offset)-p1-2*rof1;
  if (satpwr > 45)
   { printf("satpwr too large - acquisition aborted./n"); abort(1); }
  status(A);
      if (sspul[A]=='y')
       { rgpulse(1000*1e-6,zero,rof1,0.0); rgpulse(1000*1e-6,one,0.0,rof1);}
        hsdelay(d1);
        rcvroff();
        rlpower(satpwr,TODEV); txphase(v5);
        starthardloop(v13);
         rgpulse(p1,v5,rof1,rof1);
         delay(tau);
        endhardloop();
        rcvron();
        rlpower(tpwr,TODEV);
  status(B);
       rgpulse(pw,oph,rof1,rof2);
}
