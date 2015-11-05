function cut_init(self, ncutmax, ncutnz)

prob = self.prob;

if prob.K.f ~= 0
  error('cuts not supported given free variable, try presolving first')
end

prob.K.l = ncutmax + prob.K.l;
prob.J   = ncutmax + prob.J;
prob.c   = [spalloc(ncutmax,1,0); prob.c];
prob.A   = [spalloc(size(prob.A,1), ncutmax, 0), prob.A;...
            spalloc(ncutmax, ncutmax+size(prob.A,2), ncutnz)];
prob.b   = [prob.b; spalloc(ncutmax, 1, ncutmax)];

self.ncutmax = ncutmax;

% TODO: update self.A and self.b
