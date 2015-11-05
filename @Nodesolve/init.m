function init(self, prob, wp)

  % Solver datastructures
  if prob.K.f >= 1
    error(['Free variables have to be presolved away. The "transform" member of the ', class(prob), ' class instance can do this.'])
  end
  self.ksiz = [ones(1,prob.K.l), prob.K.q, prob.K.r];  
  self.ktyp = [repmat('l',1,prob.K.l), repmat('q',1,length(prob.K.q)), repmat('r',1,length(prob.K.r))];
  
  self.nx = length(prob.c);
  self.nc = length(prob.b);
  self.nk = length(self.ksiz);
  self.posvaridx = getposvaridx(prob);
  self.kbegs = find(getkbegvec(self.ksiz, self.nx, 1, 0));
  self.kidx = getkidx(self.kbegs, self.nx);
  self.e1 = getkbegvec(self.ksiz, self.nx, 1, 0);
  self.kmat0 = sparse(self.kidx, 1:self.nx, self.e1);
  self.kmat1 = sparse(self.kidx, 1:self.nx, getkbegvec(self.ksiz, self.nx, 1,-1));
  
  self.Q = getQmat(self.ksiz, self.ktyp, self.kbegs, self.nx);
  self.T = getTmat(self.ksiz, self.ktyp, self.kbegs, self.nx);
  self.Te1vec = self.T * self.e1;
  self.Te1mat = self.T * self.kmat0';
  [self.arwi,self.arwj,self.arwv] = getarwpattern(self.kidx, self.kbegs, self.nx);
  
  if self.param_stepcommonprimalduallength
    self.COMMONnx    = 2*self.nx;
    self.COMMONkidx  = [self.kidx, self.nk + self.kidx];
    self.COMMONkmat1 = blkdiag(self.kmat1, self.kmat1);
    self.COMMONQ     = blkdiag(self.Q, self.Q);
    self.COMMONT     = blkdiag(self.T, self.T);
  end

  % System ordering
  scalpatA  = prob.A*self.kmat1';
  sysperm   = symamd(tril(scalpatA*scalpatA'));
  self.sysP = sparse(sysperm, 1:self.nc, 1);
  self.PTA  = self.sysP' * prob.A;
  self.PTb  = self.sysP' * prob.b;
  
  % Initialize execution flags
  self.sys_cholerr = false;
  self.sys_ldlerr = false;
  
  % Compute coldpoint
  cp = self.getcoldpoint();
  self.wpdir = false;
  
  if nargin == 3

    disp('Computing search direction from warmpoint')
    self.tau = wp.tau;
    self.kap = wp.kap;
    self.xx  = wp.xx;
    self.ss  = wp.ss;
    self.yy  = wp.yy;
    
    if self.stopcriteria_verifypoint(prob) ~= 0
      error('Warmpoint needs to be feasible')
    end
    
    self.heuristic_repair(prob);
    self.stepdir(prob, false);
    
    self.tau = self.tau + self.dtau;
    self.kap = self.kap + self.dkap;
    self.xx  = self.xx + self.dxx;
    self.ss  = self.ss + self.dss;
    self.yy  = self.yy + self.dyy;
    
    self.pfeas = norm(prob.A *self.xx - prob.b*self.tau, Inf);
    self.dfeas = norm(prob.A'*self.yy + self.ss - prob.c*self.tau, Inf);
    self.gfeas = norm(-prob.c'*self.xx + prob.b'*self.yy - self.kap, Inf);

    self.stopcriteria_verifypoint(prob)
    
    self.dtau = (wp.tau + self.dtau) - cp.tau;
    self.dkap = (wp.kap + self.dkap) - cp.kap;
    self.dxx  = (wp.xx + self.dxx) - cp.xx;
    self.dss  = (wp.ss + self.dss) - cp.ss;
    self.dyy  = (wp.yy + self.dyy) - cp.yy;
    
    self.wpdir = true;
    
    % Reset execution flags
    self.sys_cholerr = false;
    
  end
  
  % Set initial point
  self.tau = cp.tau;
  self.kap = cp.kap;
  self.xx  = cp.xx;
  self.ss  = cp.ss;
  self.yy  = cp.yy;
  
  self.pfeas = norm(prob.A *self.xx - prob.b*self.tau, Inf);
  self.dfeas = norm(prob.A'*self.yy + self.ss - prob.c*self.tau, Inf);
  self.gfeas = norm(-prob.c'*self.xx + prob.b'*self.yy - self.kap, Inf);
  self.termstat = Inf;
  
  warning('off','MATLAB:nearlySingularMatrix')
  %warning('off','MATLAB:singularMatrix')
end

function posvaridx = getposvaridx(prob)
  qbegs = cumsum([1,prob.K.q]);
  rbegs = cumsum([1,prob.K.r]);
  posvaridx = [1:prob.K.l, prob.K.l + ...
              [qbegs(1:end-1), (qbegs(end)-1) + ...
              [rbegs(1:end-1), 1+rbegs(1:end-1)]]];
end

function kbegvec = getkbegvec(ksiz, nx, hitval, missval)
  kbegs = [1, 1+cumsum(ksiz(1:end-1))];
  kbegvec = ones(nx,1)*missval;
  kbegvec(kbegs) = hitval;
end

function kidx = getkidx(kbegs, nx)
  kidx = zeros(nx,1);
  kidx(kbegs) = 1;
  kidx = cumsum(kidx);
end

function Q = getQmat(ksiz, ktyp, kbegs, nx)
  j = 1:nx;
  v = getkbegvec(ksiz, nx, 1,-1);
  rotcones = kbegs(ktyp == 'r');
  
  j(rotcones)   = j(rotcones) + 1;
  j(rotcones+1) = j(rotcones+1) - 1;
  v(rotcones+1) = 1;
  
  Q = sparse(1:nx,j,v);
end

function T = getTmat(ksiz, ktyp, kbegs, nx)
  rotcones = kbegs(ktyp == 'r')'; 
  
  v = getkbegvec(ksiz, nx, 1,-1);
  v(rotcones) = v(rotcones)/sqrt(2);
  v(rotcones+1) = v(rotcones+1)/sqrt(2);
  v = [v; repmat(1/sqrt(2), 2*length(rotcones), 1)];
  
  T = sparse([1:nx, rotcones  , rotcones+1],...
             [1:nx, rotcones+1, rotcones  ],...
             v);
end

function [i,j,v] = getarwpattern(kidx, kbegs, nx)

  limkidx = kidx;
  limkidx(kbegs) = [];
  
  limrange = (1:nx)';
  limrange(kbegs) = [];
  
  %    vert            horz            diag
  i = [kbegs(limkidx); limrange;       (1:nx)'    ];
  j = [limrange;       kbegs(limkidx); (1:nx)'    ];
  v = [limrange;       limrange;       kbegs(kidx)];
  
  % This seems like the fastest ordering
  [v,p] = sort(v);
  i = i(p);
  j = j(p);
end

