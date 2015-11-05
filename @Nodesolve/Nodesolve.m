classdef Nodesolve < matlab.mixin.Copyable
  %Nodesolve Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = public)
    % Current point in homogeneous model
    xx
    tau
    yy
    ss
    kap
    
    % Solver parameters
    param_safemode
    param_stepdamping
    param_stepreduction
    param_stepcommonprimalduallength
    
    % Stop codes
    stopcodes
  end
  
  methods (Access = public)
    function self = Nodesolve()
      % Default parameters
      self.param_safemode = false;
      self.param_stepdamping = 0.99;
      self.param_stepreduction = [1, 1-eps, exp(2.^((1:15)/3)-1).^-(1), 0];
      self.param_stepcommonprimalduallength = true;
      self.stopcodes = struct('ERR',-1, 'ITER',0, 'STALL',1, 'OPT',2, 'PINF',3, 'DINF',4, 'ILL',5);
    end
    
    [stopcode,res] = optimize(self, prob, warmpoint)
    wp = getwarmpoint(self)
  end
  
  %
  % Internal methods and datastructures
  %
  properties (SetAccess = private, GetAccess = ?Test.Nodesolve)
    iter
    mu
    gam
    pfeas
    dfeas
    gfeas
    prstat
    termstat
    
    wpdir
    dxx
    dtau
    dyy
    dss
    dkap
    
    ksiz
    ktyp
    
    nx
    nc
    nk
    posvaridx
    kbegs
    kidx
    e1
    kmat0
    kmat1
    
    sysP
    PTA
    PTb
    sys_cholerr
    sys_ldlerr
    
    Q
    T
    Te1vec
    Te1mat
    arwi
    arwj
    arwv
    
    % For the aggregated vector [x; s]
    COMMONnx
    COMMONkidx
    COMMONkmat1
    COMMONQ
    COMMONT
  end
  
end