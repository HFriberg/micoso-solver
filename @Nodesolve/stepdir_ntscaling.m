function [OWm1,OWp1] = stepdir_ntscaling(self, xx, ss)
  smat = sparse(1:self.nx, self.kidx, ss);  
  xmat = sparse(1:self.nx, self.kidx, xx);
  diagsQs = sqrt(full(diag(smat'*self.Q*smat)));
  diagxQx = sqrt(full(diag(xmat'*self.Q*xmat)));
  
  theta = sqrt(diagsQs ./ diagxQx);
  sc = sqrt(full(diag(xmat'*smat)) + diagsQs.*diagxQx);
  th1 = (theta .* sc);
  th2 = (theta ./ sc);
  omega = sparse(1:self.nx,self.kidx,ss./th1(self.kidx)) + sparse(1:self.nx,1:self.nx,th2(self.kidx)) * (self.Q * xmat);
  
%   WW = omega*omega' - self.Q;
%   OWp2 = sparse(1:self.nx,1:self.nx,theta(self.kidx).^2) * WW;
%   OWm2 = sparse(1:self.nx,1:self.nx,theta(self.kidx).^(-2)) * self.Q*WW*self.Q;
  
  omega = omega / sqrt(2);
  omega = (self.Te1mat+omega) * spfun(@(x) 1./sqrt(1+x), self.Te1mat'*omega);
  
  WW = omega*omega' - self.Q;
  OWp1 = sparse(1:self.nx,1:self.nx,theta(self.kidx)) * WW;
  OWm1 = sparse(1:self.nx,1:self.nx,theta(self.kidx).^(-1)) * self.Q*WW*self.Q;
end