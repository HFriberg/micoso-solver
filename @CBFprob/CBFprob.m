classdef CBFprob < matlab.mixin.Copyable
  %CBFDATA Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
      K
      J
      c
      c0
      A
      b
      transform
  end
  
  methods (Access = public)
    self = parse(self, filepath)
    write(self, filepath)
  end
  
  methods(Access = protected)
    % Overwrites the behavior of copy() to not just copy
    % the 'transform' class object by reference
    function out = copyElement(in)
       out = copyElement@matlab.mixin.Copyable(in);
       out.transform = copy(in.transform);
       out.transform.anchor(out);
    end
  end
  
end

