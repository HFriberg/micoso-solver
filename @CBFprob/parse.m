function self = parse(self, filepath)
  p = CBFdata().parse(filepath);
  
  % Record transformations
  self.transform = Transform(self, p.varnum);

  % Sedumi style (discards objective constant and changes sign for maximization problems)
  self.A = p.A;
  self.b = -p.b;
  self.c0 = p.c0;
  if strcmp(p.objsense, 'MIN')
    self.c = p.c;
  else
    self.c = -p.c;
  end
  
  self.K = struct('f',0,'l',0,'q',[],'r',[],'s',[]);
  
  % Remove unconstrained maps
  kbegs = find(getkbegvec(p.mapstackdim', p.mapnum, 1, 0));
  kidx = getkidx(kbegs, p.mapnum);
  rowdomain = p.mapstackdomain(kidx);
  
  freedoms = strcmp(rowdomain,'F');
  self.A(freedoms,:) = [];
  self.b(freedoms) = [];
  rowdomain(freedoms) = [];
  p.mapstackdomain(freedoms) = [];
  p.mapstackdim(freedoms) = [];
  
  % Convert all domains to 'L=' by introducing slack
  slackdoms = find(~strcmp(rowdomain,'L='));
  self.c((end+1):(end+length(slackdoms))) = 0;
  self.A(slackdoms, (end+1):(end+length(slackdoms))) = -speye(length(slackdoms));
  self.transform.A = [self.transform.A, spalloc(size(self.transform.A,1),length(slackdoms),0)];
  p.varnum = p.varnum + length(slackdoms);
  
  slackdoms = find(~strcmp(p.mapstackdomain,'L='));
  p.varstackdomain = [p.varstackdomain; p.mapstackdomain(slackdoms)]; %{end+1:end+length(slackdoms)} = p.mapstackdomain{slackdoms};
  p.varstackdim = [p.varstackdim; p.mapstackdim(slackdoms)]; %  (end+1:end+length(slackdoms)) = p.mapstackdim(slackdoms);
  
  % Setup variable domains
  vartype = zeros(p.varnum,1);
  curvar = 0;
  for i = 1:length(p.varstackdomain)
    switch(p.varstackdomain{i})
      case 'F'
        self.K.f = self.K.f + p.varstackdim(i);
        vartype(curvar + (1:p.varstackdim(i))) = 1;
      case 'L+'
        self.K.l = self.K.l + p.varstackdim(i);
        vartype(curvar + (1:p.varstackdim(i))) = 2;
      case 'Q'
        self.K.q(end+1) = p.varstackdim(i);
        vartype(curvar + (1:p.varstackdim(i))) = 3;
      case 'QR'
        self.K.r(end+1) = p.varstackdim(i);
        vartype(curvar + (1:p.varstackdim(i))) = 4;
      case 'L-'
        % Change sign of coefficients
        self.A(:, curvar + (1:p.varstackdim(i))) = -self.A(:, curvar + (1:p.varstackdim(i)));
        self.c(curvar + (1:p.varstackdim(i))) = -self.c(curvar + (1:p.varstackdim(i)));
        self.transform.A(:, curvar + (1:p.varstackdim(i))) = -self.transform.A(:, curvar + (1:p.varstackdim(i)));
        
        self.K.l = self.K.l + p.varstackdim(i);
        vartype(curvar + (1:p.varstackdim(i))) = 2;
      case 'L='
        % Remove coefficients
        self.A(:, curvar + (1:p.varstackdim(i))) = [];
        self.c(curvar + (1:p.varstackdim(i))) = [];
        self.transform.A(:,curvar + (1:p.varstackdim(i))) = [];
        
        vartype(curvar + (1:p.varstackdim(i))) = [];
        curvar = curvar - p.varstackdim(i);
    end
    
    curvar = curvar + p.varstackdim(i);
  end
  
  % Stable sort
  [~,varperm] = sort(vartype);
  self.A = self.A(:,varperm);
  self.c = self.c(varperm);
  self.transform.A = self.transform.A(:,varperm);
  
  % Integer variables
  if isprop(p,'intvar')
    Jvar = sparse(p.intvar,ones(size(p.intvar)),1,size(self.A,2),1);
    Jvar = Jvar(varperm);
    self.J = find(Jvar);
  else
    self.J = [];
  end
  
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
