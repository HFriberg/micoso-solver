function cut_support(self, xval, scale_up)
% xval:     The supported points of a cut are combined rowwise. These
%           matrices are combined blockdiagonally in xval.
% scale_up: The support cut should allow xval to scale up (as oppose to 
%           scale down). Defined per cut.

nx = size(self.prob.A,2);
nc = size(xval,2) / nx;
if abs(nc-round(nc)) > eps || nc ~= length(scale_up)
  error('inconsistent dimensions')
end

if self.ncut + nc > self.ncutmax
  error('Not enough room for more cuts')
end

% compute and set coefficients
rows = size(self.prob.A,1) - self.ncutmax + self.ncut + (1:nc);
self.prob.A(rows, 1:end) = reshape(xval \ ones(size(xval,1),1), [nx, nc])';
self.prob.b(rows) = ones(nc,1);

% set slack variables
slackcols = self.ncut + (1:nc);
self.prob.A(rows(scale_up),  slackcols(scale_up))  = -speye(sum(scale_up));
self.prob.A(rows(~scale_up), slackcols(~scale_up)) =  speye(sum(~scale_up));

self.ncut = self.ncut + nc;
