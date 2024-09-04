/* SMOOTH4

   This proc evaluates the smoothed probabilities necessary to
   iterate on normal equations when variances are the same for
   the two states.  It returns:

       value of likelihood when kb = 1 or 2
       new parameter vector when kb = 3      */

proc(3)=ofn(th);
   local mu1,mu2,p,q,sig1,sig2,rho,pa,pax,qax,pfx,its,
         p1,f,f0,fy,ind,thn,yxx,fk,f00,
         t,itt,xsig1,xsig2,detsig1,detsig2;


@ First read in initial parameter values @

   ind = seqa(1947.25,0.25,n);  
   qax = 0;                                          @ prevents tripping error if not calculating smoothed probs @    
   mu1 = th[1:(nvmu2-1),1];               @ mean for state 1 @
   mu2 = th[nvmu2:(nvp-1),1];             @ mean for state 2 @
   p = th[nvp,1];         @ prob from 1 to 1 @
   q = th[nvq,1];         @ prob from 2 to 2 @
   sig1 = reshape(th[nvsig1:(nvsig2-1),1],nv,nv);
                                          @ variance for state 1 @
   sig2 = sig1;
                                          @ variance for state 2 @

   rho = (1-q)/(2-p-q);                   @ prob. that s0 = 1 @

   pa = zeros(2,1);
      pa[1,1] = rho;
      pa[2,1] = 1 - rho;
   p1 = zeros(4,1);

   pax = zeros(n,4);  /* pax denotes probabilities from filter
                          col. 1: st = 1, st-1 = 1
                          col. 2: st = 2, st-1 = 1
                          col. 3: st = 1, st-1 = 2
                          col. 4: st = 2, st-1 = 2
                          col. 5: st = 1  [added later]
                          col. 6: st-1 = 1  */
      /*   qax denotes probabilities from smoother, same
                          structure as pax */
   pfx = zeros(n,1);    @ likelihoods from filter @


@ Calculate probability weighted likelihoods for each obs @

   xsig1 = invpd(sig1);
       detsig1 = detl;
   xsig2 = invpd(sig2);
       detsig2 = detl;

   yxx = zeros(n,2);
   itt = 1;

   do until itt > n;
   yxx[itt,.] =
            exp(-(1/2)*(((y[itt,.]-mu1')*xsig1)*(y[itt,.]-mu1')')) ~
            exp(-(1/2)*(((y[itt,.]-mu2')*xsig2)*(y[itt,.]-mu2')'));
   itt = itt+1;
   endo;

   yxx = ((p/(detsig1^.5))*yxx[.,1])~(((1-p)/(detsig2^.5))*yxx[.,2])~
         (((1-q)/(detsig1^.5))*yxx[.,1])~((q/(detsig2^.5))*yxx[.,2]);

   if ndpchk(24);
      ierr = 1;
      "Underflow or overflow encountered; estimation terminated";
      "theta vector";th'; thn = th;
      ndpclex;
      goto stopnow;
   endif;


@ Next call basic filter, store results in pax, pfx  @

   its = 1;
   f = 0;
   do until its > n;
      p1[1,1] = pa[1,1]*yxx[its,1];
      p1[2,1] = pa[1,1]*yxx[its,2];
      p1[3,1] = pa[2,1]*yxx[its,3];
      p1[4,1] = pa[2,1]*yxx[its,4];
      pfx[its,1] = sumc(p1);
      f = f + ln(pfx[its,1]);

      p1 = p1/pfx[its,1];
      pax[its,.] = p1';
      pa[1,1] = p1[1,1] + p1[3,1];
      pa[2,1] = p1[2,1] + p1[4,1];

   its = its+1;
   endo;

   f00 = f;
   f0 = f - (nu/2)*(mu1'*xsig1*mu1)  - (nu/2)*(mu2'*xsig2*mu2)
                -  (alpha*ln(detsig1))
          - (alpha*ln(detsig2)) - sumc(beta*diag(xsig1))
            - sumc(beta*diag(xsig2));

@ Now calculate smoother, if desired  @

   If (kb==1) ;
       vof = f0;
       thn = vof;
   endif;
   if ((kb>1) or (kc==2)) ;

      qax = pax[1,.]~pax[1,.];
         qax[1,2]=0; qax[1,4]=0; qax[1,5]=0; qax[1,7]=0;
      t = 2;
      do until t > n;
          qax = (((yxx[t,1]*qax[.,1:4])+(yxx[t,3]*qax[.,5:8]))~
                ((yxx[t,2]*qax[.,1:4])+(yxx[t,4]*qax[.,5:8])))/
                pfx[t,1];
          qax = qax | (pax[t,.] ~ pax[t,.]);
             qax[t,2]=0; qax[t,4]=0; qax[t,5]=0; qax[t,7]=0;
      t=t+1;
      endo;

      qax = qax[.,1:4] + qax[.,5:8];


         @ Calculate filter and smoother probs that st=1
            (col. 5) and st-1=1 (col. 6) @
              pax = pax ~ (pax[.,1] + pax[.,3])
                        ~ (pax[.,1] + pax[.,2]);
         qax = qax ~ (qax[.,1] + qax[.,3])
                        ~ (qax[.,1] + qax[.,2]);


   endif;


@ Print output, if desired @

   if kc == 2;
      format /rd /m1 6,4;
      "";"filter probabilities";"";
      "Obs   st=1     st=2    st=1     st=2     st=1    st-1=1";
      "      st-1=1   st-1=1  st-1=2   st-1=2  ";
  
      ind~pax;

      "";"full sample smoother probabilities";"";
      "Obs   st=1     st=2     st=1     st=2     st=1     st-1=1";
      "      st-1=1   st-1=1   st-1=2   st-1=2  ";
      ind~qax;

      format /ro 14,6;

   endif;


@ Calculate gradient (the global vector gradz), if desired @

   if (kb == 3) ;
      vof = f0;
      gradz[1,1] = xsig1*(-nu*mu1 + sumc((y-mu1).*(qax[.,5])));
      gradz[2,1] = xsig2*(-nu*mu2 + sumc((y-mu2).*(1-qax[.,5])));
      gradz[3,1] = (1/p)*sumc(qax[2:n,1]) - (1/(1-p))*sumc(qax[2:n,2])
                  +  (1/(1-p))*(qax[1,5] - rho);
      gradz[4,1] = (1/q)*sumc(qax[2:n,4]) - (1/(1-q))*sumc(qax[2:n,3])
                  +  (1/(1-q))*(rho - qax[1,5]);
      gradz[5,1] = vec((1/2)*(sig1*sumc(qax[.,5])
                  - ((y-mu1)'*(qax[.,5].*(y-mu1)))))
                 + alpha*sig1 - beta*eye(nv) - (1/2)*nu*mu1*mu1';
        if nv == 1; gradz[5,1] = (-1/(sig1^2))*gradz[5,1]; endif;
      gradz[6,1] = vec((1/2)*(sig2*sumc(1-qax[.,5])
                    -((y-mu2)'*((1-qax[.,5]).*(y-mu2)))))
                    +  alpha*sig2 - beta*eye(nv) - (1/2)*nu*mu2*mu2';
        if nv == 1; gradz[6,1] = (-1/(sig2^2))*gradz[6,1]; endif;

      thn = vof;
    endif;


@ Produce a new iteration on normal equations, if desired  @

   if kb == 4;
      thn = zeros(nvsig2-1,1);
      p = sumc(qax[2:n,1])/
          (sumc(qax[2:n,6]) + rho - qax[1,5]);
      q = sumc(qax[2:n,4])/
          (sumc(1-qax[2:n,6]) + qax[1,5] - rho);
      mu1 = sumc(y.*qax[.,5])/
          (sumc(qax[.,5]) + nu);
      mu2 = sumc(y.*(1-qax[.,5]))/
          (sumc(1-qax[.,5])+ nu);
      sig1 = ((y'-mu1)*(qax[.,5].*(y-mu1')) +
                (y'-mu2)*((1-qax[.,5]).*(y-mu2')))/n;
      thn[1:(nvmu2-1),1] = mu1;
      thn[nvmu2:(nvp-1),1] = mu2;
      thn[nvp,1] = p;
      thn[nvq,1] = q;
      thn[nvsig1:(nvsig2-1),1] = reshape(sig1,nv*nv,1);

      deltalik = f0 - vof;
      deltach = thn-th;
      deltath = maxc(abs(deltach));

/*    "value of likelihood function:";;f0;
          "    change over previous value:";;deltalik;

          "";"iteration";;it;
       @   "new theta vector"; thn';  @
           "biggest change in param. vector";;deltath; */


      vof = f0;

   endif;

goto stopnow;


return;

stopnow:
   retp(thn,qax,ind);
   endp;



