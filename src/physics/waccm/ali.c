#include <stdio.h>
#include <stdlib.h>

int ali_(double *H,
         double *P,
         double *T,
         double *CO2,
         double *O,
         double *N2,
         double *O2,
         double *CH_SUM,
         int    *N)

{
    int ii, ND;
    
    ND=N[0];
    
    for(ii=0;ii<ND-1;ii++){
        CH_SUM[ii] = -0.0005;
    }
    
    ii = 20;
    CH_SUM[ii] = H[ii];
    ii = 10;
    CH_SUM[ii] = CO2[ii];
    
    return 0;
}
