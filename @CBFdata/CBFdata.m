classdef CBFdata < matlab.mixin.Copyable
  %CBFDATA Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    ver
    objsense
    varnum
    varstackdomain
    varstackdim
    intvar
    mapnum
    mapstackdomain
    mapstackdim
    c
    c0
    A
    b
  end
  
  methods (Access = public)
    self = parse(self, filepath)
    write(self, filepath)
  end
  
end

