function SS = steplength_compute(self, xx, dxx, tau, dtau)
  
  SStau = -tau ./ dtau; 
  SStau = SStau(SStau >= 0);

%   lst = 1:nx;
%   c = full(diag(kmat1 * sparse(lst, kidx, xx.*xx)));
%   b = full(diag(kmat1 * sparse(lst, kidx, xx.*dx)));
%   a = full(diag(kmat1 * sparse(lst, kidx, dx.*dx)));
  
  xxmat = sparse(1:self.nx, self.kidx, self.T*xx);
  dxmat = sparse(1:self.nx, self.kidx, self.T*dxx);
  c = full(diag(self.kmat1 * (xxmat.*xxmat)));
  b = full(diag(self.kmat1 * (dxmat.*xxmat)));
  a = full(diag(self.kmat1 * (dxmat.*dxmat)));
  
  % Minimum positive root of ax^2 + 2bx + c = 0
  d = real(sqrt(b.^2 - a.*c));
  SSxx = c ./ (-b + d);
  SSxx = SSxx(SSxx >= 0);
  
  SS = min([1; SStau; SSxx]);
  
end

