function [code,msg,res] = stopcriteria_classifypoint(self, prob)
  
  % Check for primal-dual feasible solution
  xscal = self.xx / self.tau;
  yscal = self.yy / self.tau;
  sscal = self.ss / self.tau;
  
  err = max([...
    norm(prob.A *xscal - prob.b, Inf) / (1 + norm(prob.b,Inf)), ...
    norm(prob.A'*yscal + sscal - prob.c, Inf) / (1 + norm(prob.c,Inf)), ...
    min(xscal'*sscal, abs(prob.b'*yscal - prob.c'*xscal)) / max(1, min(abs([prob.c'*xscal, prob.b'*yscal]))) ...
  ]);

  code = self.stopcodes.OPT;
  msg = 'primal-dual optimal solution';
  res = struct('x',xscal,'y',yscal,'s',sscal,'err',err);

  
%   % Check for primal infeasibility
%   if prob.b'*self.yy > 0
%     
%     err = (norm(prob.A'*self.yy + self.ss, Inf) * norm(prob.b,Inf)) / ...
%           ((prob.b'*self.yy) * max(1,norm(prob.c,Inf)));
% 
%     if (err < res.err)
%       code = self.stopcodes.PINF;
%       msg = 'primal infeasibility certificate';
%       
%       if err == 0
%         scal = 1;
%       else
%         scal = err * max(1,norm(prob.c,Inf)) / (norm(prob.A'*self.yy + self.ss, Inf) * norm(prob.b,Inf));
%       end
%       res = struct('x',NaN,'y',scal*self.yy,'s',scal*self.ss,'err',err);
%     end
%   end
%   
%   
%   % Check for dual infeasibility
%   if -prob.c'*self.xx > 0
%     
%     err = (norm(prob.A*self.xx, Inf) * norm(prob.c,Inf)) / ...
%           ((-prob.c'*self.xx) * max(1,norm(prob.b,Inf)));
% 
%     if (err < res.err)
%       code = self.stopcodes.DINF;
%       msg = 'dual infeasibility certificate';
%       
%       if err == 0
%         scal = 1;
%       else
%         scal = err * max(1,norm(prob.b,Inf)) / (norm(prob.A*self.xx, Inf) * norm(prob.c,Inf));
%       end
%       res = struct('x',scal*self.xx,'y',NaN,'s',NaN,'err',err);
%     end
%   end

  
  % Check for ill-posed problem
  err = max([...
    norm(prob.A'*self.yy + self.ss, Inf) * norm(prob.b,Inf) / max(1,norm(prob.c,Inf)), ...
    norm(prob.A*self.xx, Inf) * norm(prob.c,Inf) / max(1,norm(prob.b,Inf)) ...
  ]);
        
  if (err < res.err)
    code = self.stopcodes.ILL;
    msg = 'ill-posed problem certificate';
    res = struct('x',self.xx,'y',self.yy,'s',self.ss,'tau',self.tau,'kap',self.kap,'err',err);
  end
  
end