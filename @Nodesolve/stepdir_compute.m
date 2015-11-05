function [dx,dt,PTdy,ds,dk] = stepdir_compute(self,OWp1,OWm1,scalPTA,sysL,sysLT,scalc,...
  OWp1gx,PTgy,...
  PTr1,scalr2,r3,TinvarwXr4,r5)
%
% Definitions:
% scalPTA  = sysPT*A*NTm1
% scalc  = NTm1*c
% scalr2 = NTm1*r2
%

  PThy = sysLT \ (sysL \ (PTr1 + scalPTA*(scalr2 - TinvarwXr4)));
  OWp1hx = scalPTA'*PThy - (scalr2 - TinvarwXr4);
  %verifysubdir(scalA, scalAT, OWp1hx, hy, r1, scalr2 - TinvarwXr4);
  
  dt = ((r5 + self.tau*r3) + self.tau*(scalc'*OWp1hx - self.PTb'*PThy)) ...
     / (      self.kap     - self.tau*(scalc'*OWp1gx - self.PTb'*PTgy));
  
  PTdy   =     PThy +     PTgy*dt;
  OWp1dx = OWp1hx + OWp1gx*dt;
  dx     = OWm1*OWp1dx;
  
  dk = (r5 - self.kap*dt) / self.tau;
  ds = OWp1*(TinvarwXr4 - OWp1dx);
end
