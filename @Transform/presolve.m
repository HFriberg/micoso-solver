function presolve(self, eliminate_dummy_variables, ...
                        eliminate_dummy_variables_2, ...
                        eliminate_free_variables)
  prob = self.prob;
  Jvar = sparse(prob.J,ones(size(prob.J)),1,size(prob.A,2),1);
                                   
  %
  % PRESOLVE STEP 1
  % Search for constraints:  x = y,
  % where x is free, and substitute out the variable x.
  % (x is integer => y is integer)
  %
  if (eliminate_dummy_variables)
    tA = prob.A';

    % Find x=y constraints
    r = find(sum(logical(tA))==2 & sum(tA)==0.0 & prob.b'==0.0);
    [c,~] = find(tA(:,r));
    cr = c(2:2:end)';
    cl = c(1:2:end)';

    % Limit to those where cl is free
    r  = r(cl <= prob.K.f);
    cr = cr(cl <= prob.K.f);
    cl = cl(cl <= prob.K.f);

    % Replace cl by cr
    prob.A(:,cr) = prob.A(:,cr) + prob.A(:,cl);
    prob.c(cr)   = prob.c(cr)   + prob.c(cl);
    self.A(:,cr) = self.A(:,cr) + self.A(:,cl);
    Jvar(cr)     = Jvar(cr)     | Jvar(cl);

%     % Move integer requirements
%     [hitcl,hitidx] = ismember(cl, prob.J);
%     hitcr = ismember(cr, prob.J);
%     prob.J(hitidx(hitcl & ~hitcr)) = cr(hitcl);
%     prob.J(hitidx(hitcl &  hitcr)) = [];
    
    % Eliminate cl
    prob.A(:,cl) = [];
    prob.c(cl) = [];
    self.A(:,cl) = [];
    Jvar(cl) = [];
    prob.K.f = prob.K.f - length(cl);

    % Eliminate r
    prob.A(r,:) = [];
    prob.b(r) = [];
  end

  %
  % PRESOLVE STEP 2
  % Search for constraints:  y = x,
  % where x >= 0 and y is on the '>=' side of a second-order cone, 
  % and substitute out the variable x.
  % (x is integer => y is integer)
  %
  if (eliminate_dummy_variables_2)
    tA = prob.A';

    % Find x=y constraints
    r = find(sum(logical(tA))==2 & sum(tA)==0.0 & prob.b'==0.0);
    [c,~] = find(tA(:,r));
    cr = c(2:2:end)';
    cl = c(1:2:end)';

    % Limit to those where cl is linear
    r  = r(prob.K.f < cl&cl <= prob.K.f+prob.K.l);
    cr = cr(prob.K.f < cl&cl <= prob.K.f+prob.K.l);
    cl = cl(prob.K.f < cl&cl <= prob.K.f+prob.K.l);
    
    % Limit to those where cr is on the '>=' side of a second-order cone
    qbegs = cumsum([1,prob.K.q]);
    rbegs = cumsum([1,prob.K.r]);
    grcvar = prob.K.f + prob.K.l + [qbegs(1:end-1), ...
                  (qbegs(end)-1) + [rbegs(1:end-1), 1+rbegs(1:end-1)]];
    r  = r(ismember(cr, grcvar));
    cr = cr(ismember(cr, grcvar));
    cl = cl(ismember(cr, grcvar));

    % Replace cl by cr
    prob.A(:,cr) = prob.A(:,cr) + prob.A(:,cl);
    prob.c(cr)   = prob.c(cr)   + prob.c(cl);
    self.A(:,cr) = self.A(:,cr) + self.A(:,cl);
    Jvar(cr)     = Jvar(cr)     | Jvar(cl);
    

    % Move integer requirements
%     [hitcl,hitidx] = ismember(cl, prob.J);
%     hitcr = ismember(cr, prob.J);
%     prob.J(hitidx(hitcl & ~hitcr)) = cr(hitcl);
%     prob.J(hitidx(hitcl &  hitcr)) = [];
    
    % Eliminate cl
    prob.A(:,cl) = [];
    prob.c(cl) = [];
    self.A(:,cl) = [];
    Jvar(cl) = [];
    prob.K.l = prob.K.l - length(cl);

    % Eliminate r
    prob.A(r,:) = [];
    prob.b(r) = [];
  end
  
  %
  % PRESOLVE STEP 3
  % Search for all variables:  x free,
  % and replace by (x',x) in conic quadratic domain.
  %
  if (eliminate_free_variables && prob.K.f >= 1)
    
    
    prob.A = [prob.A(:,prob.K.f + (1:prob.K.l)),  ...
              spalloc(size(prob.A,1),1,0),        ...
              prob.A(:,1:prob.K.f),               ...
              prob.A(:,(prob.K.f + prob.K.l + 1):end)];
            
    prob.c = [prob.c(prob.K.f + (1:prob.K.l));  ...
              0;                                ...
              prob.c(1:prob.K.f);               ...
              prob.c((prob.K.f + prob.K.l + 1):end)];
            
    self.A = [...
      self.A(:,prob.K.f + (1:prob.K.l)),  ...
      spalloc(size(self.A,1),1,0),        ...
      self.A(:,1:prob.K.f),               ...
      self.A(:,(prob.K.f + prob.K.l + 1):end)];
    
    prob.K.q = [1+prob.K.f, prob.K.q];
    prob.K.f = 0;
    
  end
  
  prob.J = find(Jvar);
end
