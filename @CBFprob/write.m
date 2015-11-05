function write(self, filepath)

  p = CBFdata();

  p.ver = 1;
  p.objsense = 'MIN';  
  p.varnum = length(self.c);
  p.varstackdomain = [repmat({'F'},sign(self.K.f),1); repmat({'L+'},sign(self.K.l),1); repmat({'Q'},length(self.K.q),1); repmat({'QR'},length(self.K.r),1)];
  p.varstackdim = [repmat(self.K.f,sign(self.K.f),1); repmat(self.K.l,sign(self.K.l),1); self.K.q'; self.K.r'];
  p.intvar = self.J;
  p.mapnum = length(self.b);
  p.mapstackdomain = {'L='};
  p.mapstackdim = p.mapnum;
  p.c = self.c;
  p.c0 = self.c0;
  p.A = self.A;
  p.b = -self.b;
  
  p.write(filepath);
  