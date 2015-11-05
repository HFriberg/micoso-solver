function sol = postsolve(self, sol)
  sol = self.A * sol + self.b;
end
