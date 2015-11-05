function [xx,tau,yy,ss,kap,pfeas,dfeas,gfeas] = step(self, prob, dxx, dtau, dyy, dss, dkap)
  if nargin() == 2
    dxx  = self.dxx;
    dtau = self.dtau;
    dyy  = self.dyy;
    dss  = self.dss;
    dkap = self.dkap;
  end

  [ssp,ssd] = self.steplength(dxx, dtau, dss, dkap);
  
  xx  = self.step_compute(self.xx , dxx , ssp);
  if any(xx(1:prob.K.l) <= 0)
    disp('ohh no!');
  end
  tau = self.step_compute(self.tau, dtau, ssp);
  yy  = self.step_compute(self.yy , dyy , ssd);
  ss  = self.step_compute(self.ss , dss , ssd);
  kap = self.step_compute(self.kap, dkap, ssd);
  
  % Measures of progress
  pfeas = norm(prob.A *xx - prob.b*tau, Inf);
  dfeas = norm(prob.A'*yy + ss - prob.c*tau, Inf);
  gfeas = norm(-prob.c'*xx + prob.b'*yy - kap, Inf);
  
  if nargout == 0
    if max([pfeas, dfeas, gfeas]) <= max([self.pfeas, self.dfeas, self.gfeas])
      self.xx  = xx;
      self.tau = tau;
      self.yy  = yy;
      self.ss  = ss;
      self.kap = kap;
      self.pfeas = pfeas;
      self.dfeas = dfeas;
      self.gfeas = gfeas;
    end
  end
end
