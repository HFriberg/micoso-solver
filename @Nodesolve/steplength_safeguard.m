function SS = steplength_safeguard(self, xx, dxx, SSunsafe)
  n = 0;
  unsafe = true;  
  
  while unsafe
    n = n + 1;
    SS = self.param_stepreduction(n)*SSunsafe;
    xxstep = self.step_compute(xx, dxx, SS);
    xxmat  = sparse(1:self.nx, self.kidx, xxstep);
    unsafe = any(xxstep(self.posvaridx) <= 0) || any(diag(xxmat'*self.Q*xxmat) <= 0);
  end
end

% function SS = steplength_safeguard2(self, nx, kidx, Q, xx, dxx, stepdamping, stepreduction, SSunsafe)
%   n = length(stepreduction);
%   safe = true;
%   
%   while safe && (n >= 2)
%     n = n - 1;
%     SS = stepreduction(n)*SSunsafe;
%     xxmat = sparse(1:nx, kidx, xx + stepdamping*SS*dxx);
%     safe = all(diag(xxmat'*Q*xxmat) > 0);
%   end
%   
%   if ~safe
%     SS = stepreduction(n+1)*SSunsafe;
%   end
% end