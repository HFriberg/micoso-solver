function stepdir(self, prob, centercorrector)
  if nargin ~= 3
    centercorrector = true;
  end
  
  pdnorm = min(norm(self.xx, 2), norm(self.ss, 2));
  %self.xx = self.xx / pdnorm;
  self.ss = self.ss / pdnorm;
  self.yy = self.yy / pdnorm;
  
  %
  % Compute NT scaling
  %
  [OWm1,OWp1] = self.stepdir_ntscaling(self.xx, self.ss);
  TinvarwX = self.T * self.invarw(self.T*OWp1*self.xx);
  %Test.Nodesolve.verifyntscaling(self, self.xx, self.ss, OWm1, OWp1);
  
  %
  % Scale and factorize system
  %
  PTyy = self.sysP'*self.yy;
  scalPTA = self.PTA*OWm1;
  sysLT = self.cholinf( scalPTA*scalPTA' );
  sysL = sysLT';
  %Test.Nodesolve.verifyfactorization(sysLT, sysL, scalPTA*scalPTA');
  
  %
  % Compute subdirection (gx and gy)
  %
  scalc = OWm1*prob.c;
  PTgy = sysLT \ (sysL \ (self.PTb + scalPTA*scalc));
  OWp1gx = scalPTA'*PTgy - scalc;  
  %Test.Nodesolve.verifysubdir(scalPTA, OWp1gx, PTgy, self.PTb, scalc);

  %
  % Mehrotra's predictor direction
  %
  PTr1 = self.tau*self.PTb - self.PTA*self.xx;
  r2 = self.tau*prob.c - self.PTA'*PTyy - self.ss;
  r3 = self.kap + prob.c'*self.xx - self.PTb'*PTyy;
  
  if ~self.wpdir
    
    TinvarwXr4 = -OWm1*self.ss;
    r5 = -self.kap*self.tau;

    [dxx1, dtau1, PTdyy1, dss1, dkap1] = ...
    self.stepdir_compute(OWp1,OWm1,scalPTA,sysL,sysLT,scalc,...
          OWp1gx,PTgy,...
          PTr1,OWm1*r2,r3,TinvarwXr4,r5);
    %Test.Nodesolve.verifydir(self, prob, OWp1, OWm1, dxx1, dtau1, PTdyy1, dss1, dkap1, PTr1, r2, r3, TinvarwXr4, r5, TinvarwX);
    
  else
    
    dxx1   = self.dxx;
    dtau1  = self.dtau;
    PTdyy1 = self.sysP'*self.dyy;
    dss1   = self.dss;
    dkap1  = self.dkap;
    
  end
  
  %
  % Analyze progress and choose aggresiveness
  %
  self.prstat = dtau1/self.tau - dkap1/self.kap;
  self.mu = (self.xx'*self.ss + self.tau*self.kap)/(self.nk+1);
  self.gam = 0;   % also used in steplength computation
  if centercorrector
    prediction = 1 - self.steplength(dxx1, dtau1, dss1, dkap1);
    self.gam = prediction*min(0.3, prediction^2);
  end
  
  % Experimental
  self.gam = 0.5;
  
  %
  % Mehrotra's center-corrector direction
  %
  PTr1 = -self.gam*PTr1;
  r2 = -self.gam*r2;
  r3 = -self.gam*r3;
  r4 = self.gam*self.mu*self.e1;
  r5 = self.gam*self.mu;
  
  if centercorrector && ~self.wpdir
    % Add corrector terms
    r4 = r4 - self.arw(self.T*OWp1*dxx1) * (self.T*OWm1*dss1);
    r5 = r5 - dkap1*dtau1;
  end
  
  [dxx2, dtau2, PTdyy2, dss2, dkap2] = ...
  self.stepdir_compute(OWp1,OWm1,scalPTA,sysL,sysLT,scalc,...
        OWp1gx,PTgy,...
        PTr1,OWm1*r2,r3,TinvarwX*r4,r5);
  %Test.Nodesolve.verifydir(self, prob, OWp1, OWm1, dxx2, dtau2, PTdyy2, dss2, dkap2, PTr1, r2, r3, r4, r5);
  
  %
  % Final step direction
  %
  self.dxx  =  dxx1 +  dxx2;
  self.dtau = dtau1 + dtau2;
  self.dyy  = self.sysP*(PTdyy1 + PTdyy2);
  self.dss  =  dss1 +  dss2;
  self.dkap = dkap1 + dkap2;
  
end
