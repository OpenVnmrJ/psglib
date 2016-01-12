/* t1cp - T1 measurement by cp for solids
     D.A. Torchia, JMR ,v30, 613, (1978)

normal xpolar parameters are used, with the following additions

 cp90 - carbon 90 deg pulse
 d2   - delay for T1 measurement

E. Williams (palo alto)  ,Oct 1987   */

#include <standard.h>

pulsesequence()

{
double level1,level2,p2,cp90,level1f,level2f;
extern double getval();

level1=getval("level1");
level2=getval("level2");
level1f=getval("level1f");
level2f=getval("level2f");
cp90=getval("cp90");
p2=getval("p2");

initval(level2,v14);
initval(level1,v13);
initval(level2f,v12);
initval(level1f,v11);
hlv(ct,v1);         /* 0,0,1,1,2,2,3,3...*/
assign(v1,v2);
incr(v2);           /* 1,1,2,2,3,3,0,0...*/
mod2(ct,v3);        /* 0,1,0,1...........*/
dbl(v3,v3);         /* 0,2,0,2,0,2,......*/
add(v1,v3,oph);     /* 0,2,1,3,2,0,3,1...*/
mod4(v1,v1);
mod4(v2,v2);
mod4(oph,oph);

decphase(zero);
txphase(zero);      /* to remind us*/

status(A);
power(v13,DODEV);
pwrf(v11,DODEV);

delay(d1);

decpulse(pw,oph);   /* 1H polarization*/
decphase(v2);
status(C);          /* dec spin lock*/
rcvroff();
rgpulse(p2,v1,rof1,0.0); /* contact*/
status(A);          /* off with 1H*/
rgpulse(cp90,v2,1.0e-6,0.0);
rcvron();

delay(d2-rof1);          /*relaxation delay*/

rcvroff();
rgpulse(cp90,v1,rof1,0.0);  /* return 13C to xy */
power(v14,DODEV);     /* DD decouple */
pwrf(v12,DODEV);
status(C);
delay(rof2);
rcvron();
}

