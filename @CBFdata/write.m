function write(self, filepath)

  f = fopen(filepath, 'wt');
  if (f == -1)
    error(['fopen failed: ', filepath])
  end
  
  fprintf(f, 'VER\n%i\n\n', self.ver);
  
  fprintf(f, 'OBJSENSE\n%s\n\n', self.objsense);
  
  if self.varnum >= 1
    fprintf(f, 'VAR\n%u %u\n', self.varnum, length(self.varstackdim));
    data = [self.varstackdomain'; num2cell(self.varstackdim')];
    fprintf(f, '%s %u\n', data{:});
    fprintf(f, '\n');
  end
  
  if length(self.intvar) >= 1
    fprintf(f, 'INT\n%u\n', length(self.intvar));
    fprintf(f, '%u\n', self.intvar-1);
    fprintf(f, '\n');
  end
  
  if self.mapnum >= 1
    fprintf(f, 'CON\n%u %u\n', self.mapnum, length(self.mapstackdim));
    data = [self.mapstackdomain'; num2cell(self.mapstackdim')];
    fprintf(f, '%s %u\n', data{:});
    fprintf(f, '\n');
  end
  
  if nnz(self.c) >= 1
    fprintf(f, 'OBJACOORD\n%u\n', nnz(self.c));
    [i,~,v] = find(self.c);
    fprintf(f, '%u %.16g\n', [i-1,v]');
    fprintf(f, '\n');
  end
  
  if self.c0 ~= 0
    fprintf(f, 'OBJBCOORD\n%.16g\n\n', self.c0);
  end

  if nnz(self.A) >= 1
    fprintf(f, 'ACOORD\n%u\n', nnz(self.A));
    [i,j,v] = find(self.A);
    fprintf(f, '%u %u %.16g\n', [i-1,j-1,v]');
    fprintf(f, '\n');
  end
  
  if nnz(self.b) >= 1
    fprintf(f, 'BCOORD\n%u\n', nnz(self.b));
    [i,~,v] = find(self.b);
    fprintf(f, '%u %.16g\n', [i-1,v]');
    fprintf(f, '\n');
  end
  
  fclose(f);