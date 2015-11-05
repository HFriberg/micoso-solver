function R = cholinf(self, X)
%cholinf A replacement written by Michael Grant for the CVX project (copied 
% from functions/@cvx/quad_form.m and specialized to this project), for 
% the deprecated built-in MATLAB function:
%   cholinc(X, 'inf').
%
  tolLDL = 2*eps;
  
  % Any symmetric matrix can be L*D*L' factorized into a
  % lower-triangular L and block-diagonal D with block-sizes 
  % of 1x1 and 2x2 only. Fixing all 2x2 blocks, and negative 
  % 1x1 blocks, of D to zero, we obtain the approximated
  % Cholesky factorization: (L*sqrt(D))*(L*sqrt(D))'.
  %
  if ~self.sys_cholerr
    [R,self.sys_cholerr] = chol(X + eps*speye(size(X)));
    if self.sys_cholerr
      disp('Changing to LDL factorization')
    end
  end
  
  if self.sys_cholerr
    if ~self.sys_ldlerr
      try
        [ R, DD, prm ] = ldl(X, 'upper', 'vector');
      catch
        self.sys_ldlerr = true;
        disp('Changing to unstable LDL factorization');
      end
    end
    
    if self.sys_ldlerr
      [ R, DD, prm ] = ldl(X, eps, 'upper', 'vector');
    end
    
    tt = diag(DD,1) == 0;
    DD = diag(DD);
    tt = [ tt ; true ] & [ true ; tt ] & ( DD > tolLDL*sum(diag(X)) );
    R  = bsxfun( @times, sqrt( DD(tt) ), R(tt,:) );
    if any( diff(prm) ~= 1 )
      R(:, prm) = R;
    end
  end
end
