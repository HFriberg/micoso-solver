classdef Transform < matlab.mixin.Copyable
  %CBFDATA Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = public)
    A
    b
    ncutmax
    ncut
  end
  
  methods (Access = public)
    presolve(self, eliminate_dummy_variables, eliminate_dummy_variables_2, eliminate_free_variables)
    sol = postsolve(self, sol)
    
    cut_init(self, ncutmax, ncutnz)
    cut_reset(self)
    cut_bounds(self, blj, blv, buj, buv)
    cut_support(self, xval, scale_up)
  end
  
  %
  % Construction
  %
  methods (Access = public)
    function self = Transform(prob, varnum)
      self.A = speye(varnum);
      self.b = zeros(varnum, 1);
      self.ncutmax = 0;
      self.ncut = 0;
      self.anchor(prob);
    end
    
    function anchor(self, prob)
      self.prob = prob;
    end
  end
  
  properties (Access = private)
    prob
  end
  
end

