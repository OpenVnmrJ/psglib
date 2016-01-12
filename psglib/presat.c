/*  presat-
              {trim(x)trim(y)} ......d1...satdly...p1...d2...pw....at......
	standard two-pulse sequence with optional composite observe pulse 
        satmod = 'y' : obs xmtr saturation at satfrq with power satpwr
                     (use like dm, i.e. satmode='yyn' or 'ynn')
        sspul = 'y' does trim(x)trim(y) to destroy all magnetization
        composit='y': uses composite 90 for pw (discriminates relative to B1)
               G.Gray, Palo Alto  Sept. 1991   version 4.1 vnmr
*/
#include <standard.h>
pulsesequence()
{
  double satdly,satpwr,satfrq;
  char composit[MAXSTR],satmode[MAXSTR],sspul[MAXSTR];
  satdly = getval("satdly"); satfrq = getval("satfrq"); satpwr = getval("satpwr");
  getstr("composit",composit); getstr("satmode",satmode); getstr("sspul",sspul); 
  initval(satpwr,v9); initval(tpwr,v10);
       dbl(ct,v1);      /*v1=02020202	*/
       hlv(ct,v2);      /* v2=00112233	*/
       mod2(v2,v3);     /* v3=00110011	*/
       add(v1,v3,v4);   /* v4=02130213	*/
       add(v2,one,v5);  /* v5=11223300	*/
       mod4(v4,oph);    /* oph=02130213	*/
  if (satpwr > 45)
   { printf("satpwr too large - acquisition aborted./n"); abort(1); }
  status(A);
      if (sspul[A]=='y')
       { rgpulse(1000*1e-6,zero,rof1,0.0); rgpulse(1000*1e-6,one,0.0,rof1);}
      hsdelay(d1);
      if (satmode[A] == 'y') 
       {
        if (fabs(tof-satfrq)>0.0) offset(satfrq, TODEV);
        power(v9,TODEV); txphase(v5);
        rgpulse(satdly,v5,rof1,rof1); 
        power(v10,TODEV);
        if (fabs(tof-satfrq)>0.0)
         {  offset(tof,TODEV); delay(40.0e-6); }
       }
  status(B);
      pulse(p1,zero);
      if (satmode[B] == 'y') 
       {
        if (fabs(tof-satfrq)>0.0) offset(satfrq,TODEV);
        power(v9,TODEV); txphase(v5); 
        rgpulse(d2,zero,rof1,rof1); 
        power(v10,TODEV);
        if (fabs(tof-satfrq)>0.0)
         {  offset(tof,TODEV); delay(40.0e-6); }
       }
      else 
        hsdelay(d2);
     status(C);
      if (composit[A] == 'y')
       {
       add(oph,one,v2); add(oph,two,v3); add(oph,three,v4);
       rgpulse(pw,v2,rof1,1.0e-6);   /* 90(+y)90(-x)90(-y)90(x) */   
       rgpulse(pw,v3,0.0,1.0e-6);
       rgpulse(pw,v4,0.0,1.0e-6);
       rgpulse(pw,oph,0.0,rof2);
       }
     else
       rgpulse(pw,oph,rof1,rof2);
}
