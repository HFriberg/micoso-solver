function [SSprim, SSdual] = steplength(self, dxx, dtau, dss, dkap)
  
  SSprim = self.steplength_compute(self.xx, dxx, self.tau, dtau);
  SSprim = self.steplength_safeguard(self.xx, dxx, SSprim);

  SSdual = self.steplength_compute(self.ss, dss, self.kap, dkap);
  SSdual = self.steplength_safeguard(self.ss, dss, SSdual);

  if self.param_stepcommonprimalduallength
    SSprim = min(SSprim, SSdual);
    SSdual = SSprim;
  end
  
%   % Scaled safeguard
%   n = 0;
%   unsafe = true;
%   
%   while unsafe
%     n = n + 1;
%     SSprim = self.param_stepreduction(n) * SSprim;
%     SSdual = self.param_stepreduction(n) * SSdual;
%     
%     xxstep = self.step_compute(self.xx, dxx, SSprim);
%     ssstep = self.step_compute(self.ss, dss, SSdual);
%     [OWm1,OWp1] = self.stepdir_ntscaling(xxstep, ssstep);
% 
%     xxmat  = sparse(1:self.nx, self.kidx, OWp1*xxstep);
%     ssmat  = sparse(1:self.nx, self.kidx, OWm1*ssstep);
%     unsafe = any(xxstep(self.posvaridx) <= 0) || any(diag(xxmat'*self.Q*xxmat) <= 0);
%     unsafe = unsafe || any(ssstep(self.posvaridx) <= 0) || any(diag(ssmat'*self.Q*ssmat) <= 0);
%     
%     xxmat  = sparse(1:self.nx, self.kidx, xxstep);
%     ssmat  = sparse(1:self.nx, self.kidx, ssstep);
%     unsafe = unsafe || any(xxstep(self.posvaridx) <= 0) || any(diag(xxmat'*self.Q*xxmat) <= 0);
%     unsafe = unsafe || any(ssstep(self.posvaridx) <= 0) || any(diag(ssmat'*self.Q*ssmat) <= 0);
%   end
  
%     
%     plotx = [linspace(0, 1, 200), SSprim];
%     ploty = zeros(size(plotx));
%     for i = 1:length(plotx);
%       xxmat = sparse(1:self.nx, self.kidx, self.xx + 0.99*plotx(i)*dx);
%       ploty(i) = min(diag(xxmat'*self.Q*xxmat));
%     end
%     plot(plotx(1:end-1),ploty(1:end-1),'-',plotx(end),ploty(end),'o')

%   end
  %disp({SSprim,SSdual})
end
