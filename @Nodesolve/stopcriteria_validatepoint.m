function [code,msg] = stopcriteria_validatepoint(self, prob)
  
  %
  % Validate internal state
  %
  if ~all([isfinite(self.tau); isfinite(self.kap); isfinite(self.xx); isfinite(self.ss); isfinite(self.yy)]);
    code = self.stopcodes.ERR;
    msg = 'NaNs or Infs produced in solution vectors';
    return;
  end
  if (self.tau <= 0)
    code = self.stopcodes.ERR;
    msg = 'tau is not positve';
    return;
  end
  if (self.kap <= 0)
    code = self.stopcodes.ERR;
    msg = 'kappa is not positve';
    return;
  end
  if any(self.xx(1:prob.K.l) <= 0)
    code = self.stopcodes.ERR;
    msg = 'x in K.l is not positve';
    return;
  end
  if any(self.ss(1:prob.K.l) <= 0)
    code = self.stopcodes.ERR;
    msg = 's in K.l is not positve';
    return;
  end
  xmat = sparse(1:self.nx, self.kidx, self.xx);
  if any(diag(xmat'*self.Q*xmat) <= 0)
    code = self.stopcodes.ERR;
    msg = 'x is not in cone interior';
    return;
  end
  smat = sparse(1:self.nx, self.kidx, self.ss);
  if any(diag(smat'*self.Q*smat) <= 0)
    code = self.stopcodes.ERR;
    msg = 's is not in cone interior';
    return;
  end

  code = self.stopcodes.ITER;
  msg = '';
  
end