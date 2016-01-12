/* makeSS - C code to make S or SS shaped pulses.
            These are amplitude-only modulated shaped pulses
             with either a full (SS) or half (S) cycle of cosine
             amplitude.  
     Ref:  Stephen H. Smallcombe, JACS 115, 4776-4785 (1993)

            Valid options are:
              -g   add 2 lines before and after shape that are
                   gated off to permit WFG amplitude to settle.
              -s   Make S pulse instead of default SS pulse.
              -f   Specify an excitation frequency instead of pulse width.
              -h   Display help screen
              -r   Specify resolution

            Required arguement is length of pulse in microseconds, or
               frequency of excitation maximum with -f option.
*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>

void displayhelp()
{
printf("makeSS - Create S and SS cosine modulated shaped pulse file\n");
printf("Usage: makeSS [-sgrfhd] [resolution] pulse_width|frequency\n\n");
printf("Options:  (May be upper or lower case, contrary to normal unix.)\n");
printf("      -s  Create S rather than SS pulse shape (half cosine).\n");
printf("      -g  Added gating statements to top and bottom of shape.\n");
printf("      -r  Specify resolution in microseconds.  Default=0.2\n");
printf("      -f  Interpret argument as frequency max instead of pulse width.\n");
printf("      -h  Display this help screen.\n");
printf("      -d  Display extra debugging output.\n\n");
printf("Examples:\n");
printf("   makeSS -gf 1500  : Make SS pulse with excitation maxima at +/-1500 Hz\n");
printf("                      from carrier, add gating elements, use 200 nsec res.\n");
printf("   makeSS 600       : Make SS pulse with 600 microsecond width, no gating,\n");
printf("                      default 200 nsec resolution.\n");
printf("   makeSS -sr 2 500 : Make S pulse of 500 microsecond width, no gating,\n");
printf("                      2 microsecond resolution.\n");
printf(" \nCAUTION: No check is made to verify a unique file name.  If a file with\n");
printf("  the same name exists, it will be overwritten.  There is no provision here\n");
printf("  for writing the output to anything other than the current directory.\n"); 
}

int main(argc, argv)
int argc;
char *argv[];
{
int i,j,n,step;
int gateflag,freqflag,typeflag,helpflag,resflag;
int debug;
FILE *output;
char fname[24],flag;
float length,fmax,cycle,pwfactor,amplitude,resolution;
struct timeval tp;
struct timezone tzp;
struct tm *time;
char *username,*timenow;

if (gettimeofday(&tp, &tzp)!=0) { 
    printf("Error getting current system time\n");
    return(-1); }
time=localtime(&tp.tv_sec);
timenow=asctime(time);
username=getenv("LOGNAME");
printf("MakeSS - Creating a cosine-modulated shaped pulse file\n");
printf("   for %s on %s",username,timenow);
if (argc<2) {printf("You must supply at least the desired pulse width!\n");
             return(1); }
length=0.0;
resolution=0.0;
freqflag=0; typeflag=0; debug=0; gateflag=0; resflag=0; helpflag=0;
for (i=1;i<argc;i++) {
       if (*argv[i]=='-') 
          for (j=1;j<strlen(argv[i]);j++){
          flag=toupper(*(argv[i]+j));
          if (flag=='G') gateflag=1;
            else if (flag=='F') freqflag=1;
            else if (flag=='S') typeflag=1;
            else if (flag=='D') debug=1;
            else if (flag=='R') resflag=1;
            else if (flag=='H') helpflag=1;
            else printf("Ignoring unknown option %c\n",*(argv[i]+j)); }
        else 
          if (length!=0.0) {printf("Too many arguments supplied!\n");
                            return(2); }
          else if ((resflag==1)&&(resolution==0.0))
                         { if (sscanf(argv[i],"%e",&resolution)==0){
                         printf("Error with %s as numeric argument\n",argv[i]);
                         return(3); } }
                else if (sscanf(argv[i],"%e",&length)==0) {
                        printf("Error with %s as numeric argument\n",argv[i]);
                        return(3); }
     }
if (length==0.0) { printf("You must supply the desired pulse width!\n\n");
                   displayhelp();
                   return(4); }
if (helpflag==1) displayhelp();
if (typeflag==1) {pwfactor=666667.0;  cycle=0.5;} 
            else {pwfactor=1100000.0;  cycle=1.0;}
if (debug) printf("Requested pulse length=%f\n",length);
if (freqflag==1) {fmax=length; length=pwfactor/fmax;}
if (debug) printf("Pulse = %f, fmax= %f\n",length,fmax);
if (resolution==0.0) {
   resolution=0.2;
   step=1; }
 else {
   step=(int)(resolution/0.2+0.5);
   resolution=0.2*(float)step; }
n=(int)(length/resolution+.5);
length=(float)n*resolution;
fmax=pwfactor/length;
if (gateflag==1) length+=0.85;
i=1;
fname[0]='S';
if (typeflag!=1) fname[i++]='S';
fname[i]='\0';
if (debug) printf("Partial name= %s\n",fname);
sprintf(fname+i,"%d",(int)(length*10));
for (;fname[i]!='\0';i++);
if (debug) printf("Partial name= %s, index= %i\n",fname,i);
if (gateflag==1)   sprintf(fname+i,"g.RF\0"); 
              else sprintf(fname+i,".RF\0");
if (debug) printf("Output file name will be %s\n",fname);
if ((output=fopen(fname,"w"))==NULL)
  { printf("Unable to open %s as output file\n",fname);
    return(5); }
length=(float)((int)(length*10))/10;
if (typeflag==1)
fprintf(output,"# S-pulse, shape analog of 1-1 binomial pulse.\n");
else
fprintf(output,"# SS pulse, shaped pulse analog of 1-2-1 binomial pulse.\n");
fprintf(output,"# Reference: Stephen H. Smallcombe, JACS 115, 4776-4785 (1993).\n");
fprintf(output,"# Shaped pulse file generated with makeSS program.\n");
fprintf(output,"# Created by user %s\n",username);
fprintf(output,"# Creation date: %s",timenow);
fprintf(output,"# To use properly, Pulse width should be %5.1f microseconds.\n",length);
fprintf(output,"# Excitation maxima occur +/- %5.0f Hertz from carrier.\n",fmax);
if (gateflag==1) {
    fprintf(output,"  0.0   1023.0   1.0   0\n");
    fprintf(output,"  0.0   1023.0   1.0   0\n"); }
for (i=0;i<n;i++) {
    amplitude=1023*cos(2*3.14159*cycle*i/n);
    if (amplitude<0.0) fprintf(output,"180.0   %6.1f   %d \n",0.0-amplitude,step);
                  else fprintf(output,"  0.0   %6.1f   %d \n",amplitude,step); }
if (gateflag==1) if (typeflag==1) {
    fprintf(output,"180.0   1023.0   1.0   0\n");
    fprintf(output,"180.0   1023.0   1.0   0\n"); }
  else {
    fprintf(output,"  0.0   1023.0   1.0   0\n");
    fprintf(output,"  0.0   1023.0   1.0   0\n"); }
fclose(output);
printf("Wrote file %s.\n",fname);
printf("Exitation maxima occur +/- %5.0f Hertz from carrier\n",fmax);
printf("Proper pulse width is %5.1f microseconds\n",length);
return(0);
}
