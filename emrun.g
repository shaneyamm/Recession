                   output file=icont1.out reset;
format /m1 /ro 14,6;


/* This program evaluates a specified start-up value with different
   means but the same variance by iterating on normal equations in a
    vector system

    Variable meanings:
            n = number of observations
            kc = 1 suppresses output print
                  kc = 2 prints probabilities and smoothed probabilities
            kb =1  just evaluates likelihood function
                  kb = 3 calculates gradients
                  kb = 4 iterates on normal equations */

   nv = 1;                     @ nv is the number of variables analyzed @
   nvmu2 = nv + 1;             @ nvmu2 is the index at which mu2 paramters
                                 begin @
   nvp = (2*nv) + 1;           @ nvp is the index for p @
   nvq  = (2*nv)+2;            @ nvq is the index for q @
   nvsig1 = (2*nv) + 3;        @ nvsig1 is the index at which sig1 parameters
                                 begin @
   nvsig2 = (2*nv) + 3 + (nv^2);  @ nvsig2 minus 1 is the last element @

icontrol = 1;                 @ icontrol = 0 for full-sample analysis @
                                     @ icontrol = 1 to implement loop using real-time data @

if icontrol == 0;

	nx = 233;
	ny = 159;
        load yall[nx,ny] = real_gdp_real_time.csv;
            "Data analyzed:  ";;
                     "US GDP growth";
	yq = yall[2:nx,ny];
	n2 = rows(yq);
	y = 400*ln(yq[2:n2,1] ./ yq[1:n2-1,1]);     @ growth rates @
	n = rows(y);
	iend = 2004.75;
  else;
	"Real time data used through date";;icontrol;
	nx = 233;
	ny = 159;
        load yall[nx,ny] = real_gdp_real_time.csv;
        iend = 1968.0;
endif;
ixx = zeros(nx-2,3);                             @ ixx will have the real-time smoothed probs of specified lag @
do until iend > 2004.75;
	if icontrol > 0;
	        i = 4*(iend - 1947) + 2;
		j = 4*(iend - 1965.5) + 2;
        	"----------------";"";"check: data lists this data set as known as of date";;yall[1,j];
                "and running through date";;yall[i,1];
        	yq = yall[2:i,j];
	        n2 = rows(yq);
        	y = 400*ln(yq[2:n2,1] ./ yq[1:n2-1,1]);     @ growth rates @
	        n = rows(y);
	endif;
                
    let th[5,1] =
     4.73 -1.19 0.94 0.80 12.04;
      alpha = 0; beta = .5; nu = 0;

   "";"";":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::";"";
   "Parameters for Bayesian prior (alpha, beta, nu)";
    alpha;;beta;;nu;"";

kb = 4;
kc = 1;
it = 1;
ierr = 0;             @ a value of 1 indicates an 80287 exception encountered @

   dh = 1.e-8;         @ definition of small step for gradient and hessian @
   deltalik = 10;      @ change in likelihood function @
   deltach = zeros(6,1); @ change in parameter vector over previous iteration @
   deltath  = 10;      @ biggest element of deltach @
   thtol  = 1.e-5;     @ when change in parameters is less than this,
                        convergence is declared to be obtained @
   gradz = 10+zeros(6,1);  @ gradient vector @


  vof = 0;
  kc = 1; kb = 4;

#include smooth.g;

          "starting parameter values";
          "   Means for state 1";
              th[1:(nvmu2-1),1]';
           "   Means for state 2";
               th[nvmu2:(nvp-1),1]';
           "   Prob from state 1 to state 1 (p)";
               th[nvp,1];
            "   Prob from state 2 to state 2 (q)";
                th[nvq,1];
            "   Variance-covariance matrix ";;
                reshape(th[nvsig1:(nvsig2-1),1],nv,nv);

     kb =1; {vof,qax,ind} = ofn(th); kb = 4;
     deltath = 10;  ierr = 0;
     mxit =200;    @number of iterations on normal equations @

  @ iterate on normal equations @

     it = 1;

     do until (it > mxit) or (deltath < thtol) or (ierr == 1);
                if kc == 2;
                 "iteration";;it;
                 "   Means for state 1";
                     th[1:(nvmu2-1),1]';
                 "   Means for state 2";
                     th[nvmu2:(nvp-1),1]';
                 "   Prob from state 1 to state 1 (p)";
                     th[nvp,1];
                  "   Prob from state 2 to state 2 (q)";
                      th[nvq,1];
                  "   Variance-covariance matrix ";;
                      reshape(th[nvsig1:(nvsig2-1),1],nv,nv);
                 endif;

           {th,qax,ind} = ofn(th);
     it = it+1;
     endo;

     "";"      Number of iterations required:";;it;
        "      Change in theta vector";;deltath;


     @ print out final results @

     kb=1;"";
     "";"             FINAL RESULTS";"";

"Estimated theta vector:";th';"";

"Estimated parameters:";
   "   Means for state 1";
        th[1:(nvmu2-1),1]';
   "   Means for state 2";
        th[nvmu2:(nvp-1),1]';
   "   Prob from state 1 to state 1 (p)";
        th[nvp,1];
   "   Prob from state 2 to state 2 (q)";
        th[nvq,1];
   "   Variance-covariance matrix";;
        reshape(th[nvsig1:(nvsig2-1),1],nv,nv);
     "Value of log likelihood function:";;{vof,qax,ind}=ofn(th);vof;


kc=2;kb=1;{vof,qax,ind} = ofn(th);

@"post-call check:";
ind~qax;@

if icontrol >0;
    ixx[i-2,1] = qax[n,6];  @ this stores the 1-qtr lag smoothed prob @
    ixx[i-2,2] = yall[1,j];    @ this is the quarter at the start of which data were observed @
    ixx[i-2,3] = ind[n-1,1];   @ this is the quarter about which an inference is drawn @
endif;

iend = iend + 0.25;
endo;

if icontrol == 1;
	ixx;
endif;

