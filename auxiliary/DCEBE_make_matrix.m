function [X, nablaK] = DCEBE_make_matrix(N, t, k, deriv_pattern)
%DCEBE_make_matrix Make spline matrix

% number of elements in derivative pattern
l = k+1;

% timepoint must be within valid bounds
if (t < 1) || (t > N-l)
    error('t out of bounds')
end

% decompose time into integer and fractional part
t_int = floor(t);
t_frac = t - t_int;
h = 1 - t_frac;

% make spline matrices
s = N - t_int + 1; % number of elements in the spline
% data fidelity matrix
XU = sparse(t_int-1, s);
XU(:,1) = 1;
XD = speye(s);
X = [XU; XD];
% derivative matrix
colsNablaK = ones(s-l+1, 1) * deriv_pattern;
nablaK = spdiags( colsNablaK, 0:k, s-l+1, s );

% replace line in nabla matrix to account for fractional change point
nablaK(1,1:l) = DCEBE_make_deriv_pattern(k, l, h) * sqrt(h);


end

