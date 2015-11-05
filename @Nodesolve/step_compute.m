function xx = step_compute(self, xx, dxx, SS)

% Best accuracy near optimum (always keeps distance to boundary)
%   xx = xx + self.param_stepdamping * SS * dxx;
% 
% Speedy ending near optimum (jumps to boundary when centering is near-zero)
  xx = xx + max(self.param_stepdamping, 1-self.gam>0.99) * SS * dxx;
end
