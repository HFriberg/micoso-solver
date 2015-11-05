function wp = getwarmpoint(self)
%GETWARMPOINT Summary of this function goes here
%   Detailed explanation goes here

  wp = struct('tau', self.tau, 'kap', self.kap, 'xx', self.xx, 'ss', self.ss, 'yy', self.yy);
  
end

