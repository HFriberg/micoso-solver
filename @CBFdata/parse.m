function self = parse(self, filepath)
  self.ver = 0;

  f = fopen(filepath, 'rt');
  if (f == -1)
    error(['fopen failed: ', filepath])
  end
  
  TXT = fscanf(f, '%c');
  NEWLINES = find(TXT == 10);
  ST = [1, NEWLINES(1:end-1)+1];
  EN = NEWLINES(1:end)-1;
  
  DATALINES = (ST<=EN & TXT(ST)~='#');
  ST = ST(DATALINES);
  EN = EN(DATALINES);
  
  fclose(f);

  li = 0;
  while li < length(ST)
    [li,line] = next(TXT,ST,EN,li,1);
    
    if (self.ver == 0)  
      if strcmp(line, 'VER')
        [li,line] = next(TXT,ST,EN,li,1);
        self.ver = sscanf(line, '%i');
        if (self.ver ~= 1)
          error('The version of the file format is not supported.')
        end
      else
        error('First keyword should be VER.')
      end
      
    else
      
      switch(line)
        case 'OBJSENSE'
          [li,line] = next(TXT,ST,EN,li,1);
          self.objsense = sscanf(line, '%s');
          
        case 'VAR'
          [li,line] = next(TXT,ST,EN,li,1);
          buf = sscanf(line, '%u %u');
          self.varnum = buf(1); 
          varstacknum = buf(2);
          
          [li,line] = next(TXT,ST,EN,li,varstacknum);
          buf = textscan(line, '%s %u');
          self.varstackdomain = buf{1};
          self.varstackdim = double(buf{2});
          
        case 'INT'
          [li,line] = next(TXT,ST,EN,li,1);
          intvarnum = sscanf(line, '%u');
          self.intvar = zeros(intvarnum, 1);
          
          [li,line] = next(TXT,ST,EN,li,intvarnum);
          buf = textscan(line, '%u');
          self.intvar = double(buf{1}) + 1;
          
        case 'CON'
          [li,line] = next(TXT,ST,EN,li,1);
          buf = sscanf(line, '%u %u');
          self.mapnum = buf(1); 
          mapstacknum = buf(2);
          
          [li,line] = next(TXT,ST,EN,li,mapstacknum);
          buf = textscan(line, '%s %u');
          self.mapstackdomain = buf{1};
          self.mapstackdim = double(buf{2});
          
        case 'PSDVAR'
          error('NO SUPPORT')
          
        case 'PSDCON'
          error('NO SUPPORT')

        case 'OBJFCOORD'
          error('NO SUPPORT')
          
        case 'OBJACOORD'
          [li,line] = next(TXT,ST,EN,li,1);
          objannz = sscanf(line, '%u');
          
          [li,line] = next(TXT,ST,EN,li,objannz);
          buf = textscan(line, '%u %f');
          
          self.c = sparse(double(buf{1})+1, 1, buf{2}, ...
                          self.varnum, 1, objannz);
          
        case 'OBJBCOORD'
          [li,line] = next(TXT,ST,EN,li,1);
          self.c0 = sscanf(line, '%f');
          
        case 'FCOORD'
          error('NO SUPPORT')
          
        case 'ACOORD'
          [li,line] = next(TXT,ST,EN,li,1);
          annz = sscanf(line, '%u');
          
          [li,line] = next(TXT,ST,EN,li,annz);
          buf = textscan(line, '%u %u %f');
          
          self.A = sparse(double(buf{1})+1, double(buf{2})+1, buf{3}, ... 
                          self.mapnum, self.varnum, annz);
          
        case 'BCOORD'
          [li,line] = next(TXT,ST,EN,li,1);
          bnnz = sscanf(line, '%u');
          
          [li,line] = next(TXT,ST,EN,li,bnnz);
          buf = textscan(line, '%u %f');
          
          self.b = sparse(double(buf{1})+1, 1, buf{2}, ...
                          self.mapnum, 1, bnnz);
          
        case 'HCOORD'
          error('NO SUPPORT')
          
        case 'DCOORD'
          error('NO SUPPORT')
          
        otherwise
          error(['Keyword "', line, '" not recognized!'])
      end
    end
  end
  
  if (isempty(self.c) && self.varnum >= 1)
    self.c = sparse(self.varnum, 1);
  end
  
  if (isempty(self.b) && self.mapnum >= 1)
    self.b = sparse(self.mapnum, 1);
  end
end

function [li,line] = next(TXT,ST,EN,li,num)
  line = TXT(ST(li+1):EN(li+num));
  li = li + num;
end

