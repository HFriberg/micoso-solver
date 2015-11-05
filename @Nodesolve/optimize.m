function [stopcode,res] = optimize(self, prob, warmpoint)
  if nargin == 2
    self.init(prob);
  else
    self.init(prob, warmpoint);
  end
  
  % Progress report header
  fprintf('iter\tPFEAS\tDFEAS\tGFEAS\tPRSTAT\t\tPOBJ\t\tDOBJ\t\tMU\n');
  
  self.iter = 0;
  while self.iter < 100 && self.stopcriteria(prob) == self.stopcodes.ITER
    self.stepdir(prob);
    self.progressreport(prob);
    
    self.step(prob);
    self.wpdir = false;
    self.iter  = self.iter + 1;
  end
  [stopcode,stopmsg,res] = self.stopcriteria(prob);
  self.progressreport(prob);
  disp(stopmsg);
end
