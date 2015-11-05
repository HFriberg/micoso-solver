function cut_bounds(self, blj, blv, buj, buv)

  prob = self.prob;
  blcuts = length(blj);
  bucuts = length(buj);

  if self.ncut + blcuts + bucuts > self.ncutmax
    error('Not enough room for more cuts')
  end

  offset = size(prob.A,1) - self.ncutmax + self.ncut;
  
  prob.A(offset+(1:blcuts), self.ncut+(1:blcuts)) = -speye(blcuts);
  prob.A(offset+(1:blcuts), prob.J(blj))          =  speye(blcuts);
  prob.A(offset+blcuts+(1:bucuts), self.ncut+blcuts+(1:bucuts)) = speye(bucuts);
  prob.A(offset+blcuts+(1:bucuts), prob.J(buj))                 = speye(bucuts);
  prob.b(offset+(1:blcuts+bucuts)) = [blv; buv];
  
  self.ncut = self.ncut + blcuts + bucuts;
  