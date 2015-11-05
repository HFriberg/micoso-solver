function progressreport(self, prob)  
%PROGRESSREPORT Summary of this function goes here
%   Detailed explanation goes here

  pobj = full(prob.c' * (self.xx / self.tau));
  dobj = full(prob.b' * (self.yy / self.tau));
  
  fprintf('%d\t%.1e\t%.1e\t%.1e\t%.2e\t%e\t%e\t%.1e\n', self.iter, self.pfeas, self.dfeas, self.gfeas, self.prstat, pobj, dobj, self.mu);
end

