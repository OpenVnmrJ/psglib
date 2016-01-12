/*  1-1 echo:   Sklenar and Bax, J. Magn. Reson., 74, 469(1987)  
 try 55us for tau ;  set hst to .002 hs='yyn'
 p1,p2 can be used to make small correction for radiation-damping recovery  */
#include <standard.h>

pulsesequence()
{
 double tau,p2;
     p2=getval("p2");
     tau=getval("tau");

    /* phase cycle calculation */
       mod2(ct,v1);               /* 0101010101010101....*/
       dbl(v1,v2);                /* 02020202020202020202....*/
       hlv(ct,v3);                /* 0011223344556677....*/
       hlv(v3,v4);                /* 00001111222233334444....*/
       mod4(v4,v5);               /* [0][1][2][3][0][1][2][3]..  []=repeat 4 times */
       add(two,v4,v6);            /* [2][3][0][1][2][3][0][1]..*/
       hlv(v4,v7);                /* [0][0][1][1][2][2][3][3][0][0]..*/
       hlv(v7,v7);                /* 16-0,16-1,16-2,16-3,16-0........*/
       hlv(v7,v8);                /* 32-0,32-1,32-2,32-3,32-0........*/
       mod2(v7,v9);               /* 16-0,16-1,16-0,16-1.............*/
       dbl(v9,v9);                /* 16-0,16-2,16-0,16-2.............*/
       mod2(v5,v10);              /* [0][1][0][1][0][1]..............*/
       dbl(v10,v10);              /* [0][2][0][2][0][2]..............*/
       add(v10,v9,v10);           /* [0][2][0][2][2][0][2][0]........*/
       add(v10,v8,v10);           /* increment v10 every 32 scans....*/
       add(oph,v5,v1);            /* 0123 1230 2301 3012 ............*/
       add(v1,two,v3);            /* inverse of v1*/
       add(v2,v5,v4);             /* 0202 1313 2020 3131*/     
       assign(v4,oph);
 
    /* equilibrium period */
    status(A);
     hsdelay(hst);                /* hst should be 2ms for later use */
     hsdelay(hst);                /* 4 homospoils used prior to d1   */
     hsdelay(hst);
     hsdelay(hst+d1-4*hst);

    status(B);
     rcvroff();
     rgpulse(pw,v5,rof1,rof1);
     delay(tau-2.0*rof1);
     rgpulse(pw-p1,v6,rof1,rof2);
     if (hs[B]=='y'){ hsdelay(hst); delay(hst);}
     else delay(2.0*hst);                
     rgpulse(pw,v1,rof1,rof1);
     delay(2.0*tau-2.0*rof1);
     rgpulse(pw-p2,v3,rof1,rof2);
     if (hs[B]=='y'){ hsdelay(hst); delay(hst);}
     else delay(2.0*hst);             
     rcvron();

    status(C);
}
