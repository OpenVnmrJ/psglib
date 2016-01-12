/* FLOCK  long range 13C/1H correlation  13C-detected 
	from Reynolds, et.al. Mag. Res. Chem. 27, p162-169 (1989).
	Version 1.0 by R.C. Crouch Burroughs Wellcome Co.
			RTP, NC, 27709 USA (919) 248 - 3840
			email  crouch@bwco.com
	Adapted from D. Lankin's (Searle) pascal version   */

#include <standard.h>

pulsesequence()
{
	double j1xh,jnxh,pp,pplvl,d7,d5,d6;
	char	presat[MAXSTR], flokoff[MAXSTR];

pp = getval("pp");
pplvl = getval("pplvl");
j1xh = getval("j1xh");
jnxh = getval("jnxh");
getstr("flokoff",flokoff);
getstr("presat",presat);

if (jnxh > 0.0)
	{
		d7 = 1.0/(2*jnxh);
		d5 = 1.0/(4*jnxh);
	}
else
	{
		d7 = getval("d7");
		d5 = getval("d5");
	}
d6 = 1.0/(2*j1xh);

/* PHASE CYCLING */
hlv(ct,v2); hlv(v2,v2);   /* v2 = 0000111122223333 */
add(two,v2,v3);
mod2(v2,v1); dbl(v1,v1);  /* v1 = 0000222200002222 */
sub(three,oph,v4);
add(one,v4,v4);
add(v1,v4,v4);		  /* v4 = 0321210303212103 */
add(one,v1,v6);
hlv(ct,v5);  hlv(v5,v5);
hlv(v5,v5); hlv(v5,v5);
add(v5,v1,v1);
add(v5,v2,v2);
add(v5,v6,v6);
add(v5,v3,v3);
add(v5,v4,v4);
mod4(ct,oph);
add(v5,oph,oph);

status(A);
      if (dm[0] == 'y')
      {
         fprintf(stdout, "decoupler must be set as dm=nny\n");
         abort(1);
      }
 
      if (newdecamp)
      {
	rlpower(pplvl, DODEV);
	rlpower(tpwr, TODEV);
      }
      else
      {
         declvlon();		/* sets maximum dhp value */
      }
hsdelay(d1);
if (presat[0] == 'y')
	{
		rgpulse(pw,zero,rof1,rof1);
		hsdelay(0.05);
	}
status(B);
rcvroff();
decpulse(pp,v1);
delay(d2/2.0);
if (flokoff[0] == 'y')			/* simple composite 180 */
	{
		rgpulse(pw,v1,rof1,0.0);
		rgpulse(2*pw,v6,0.0,0.0);
		rgpulse(pw,v1,0.0,rof1);
	}
else					/* bird composite 180 */
	{
		decpulse(pp,zero);
		delay(d6-rof1-1.5*pp);
		rgpulse(pw,v1,rof1,0.0);
		simpulse(2*pw,2*pp,v6,one,0.0,0.0);
		rgpulse(pw,v1,0.0,rof1);
		delay(d6-rof1-1.5*pp);
		decpulse(pp,two);
	}
delay(d2/2.0);
delay(d7/2.0);
decpulse(pp,v2);
delay(d6-rof1-1.5*pp);
rgpulse(pw,v1,rof1,0.0);
simpulse(2*pw,2*pp,v6,v2,0.0,0.0);
rgpulse(pw,v1,0.0,rof1);
delay(d6-rof1-1.5*pp);
decpulse(pp,v3);
delay(d7/2.0 - rof1);
simpulse(pw,pp,v1,v4,rof1,rof1);
delay(d5/2.0 - rof1);
decpulse(pp,v1);
delay(d6-rof1-pw-2.5*pp);
rgpulse(pw,v1,rof1,0.0);
simpulse(2*pw,2*pp,v6,v6,0.0,0.0);
rgpulse(pw,v1,0.0,rof1);
delay(d6-rof1-pw-1.5*pp);
decpulse(pp,v1);
if (newdecamp)
	{
	 rlpower(dpwr, DODEV);
	}
else
	{
	 declvloff();
	}
delay(d5/2.0);
rcvron(); decphase(zero);
status(C);
}

