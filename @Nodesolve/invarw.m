function X = invarw(self,x)
  xmat = sparse(1:self.nx, self.kidx, x);
  
  sc1 = full(diag(self.kmat1 * (xmat.*xmat)));
  sc0 = full(sqrt(diag(self.kmat0 * xmat)).^(-1));
  sc0 = sc0(self.kidx);
  sc0(self.kbegs) = -sc0(self.kbegs);
  xsc = x .* sc0 ./ sqrt(abs(sc1(self.kidx)));
  
  sc0(self.kbegs) = 0;
  X = sparse(1:self.nx, self.kidx, sign(sc1(self.kidx)).*xsc) * ...
      sparse(self.kidx, 1:self.nx, xsc) + sparse(1:self.nx, 1:self.nx, sc0.^2);
end