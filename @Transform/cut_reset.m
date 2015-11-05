function cut_reset(self)

  offset = size(self.prob.A,1) - self.ncutmax;
  
  self.prob.A(offset+1:end,:) = 0;
  self.prob.b(offset+1:end) = 0;
  self.prob.transform.ncut = 0;
  
 