function X = arw(self,x)
  X = sparse(self.arwi, self.arwj, x(self.arwv));
end