function [code,msg,res] = stopcriteria(self, prob)
%STOPCRITERIA Summary of this function goes here
%   Detailed explanation goes here
  
  res = NaN;

  %
  % Validate point
  %
  [code,msg] = self.stopcriteria_validatepoint(prob);
  if code ~= self.stopcodes.ITER
    return
  end

  %
  % Validate progress
  %  
  if max([self.pfeas, self.dfeas, self.gfeas]) + eps < self.termstat
    self.termstat = max([self.pfeas, self.dfeas, self.gfeas]);
    code = self.stopcodes.ITER;
    msg = '';
    
  else
    code = self.stopcodes.STALL;
    msg = 'stalled';
  end
  
  if nargout == 3
    [code,msg,res] = self.stopcriteria_classifypoint(prob);
  end

end