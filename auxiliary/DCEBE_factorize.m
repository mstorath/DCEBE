function [U, d, p, q, r] = DCEBE_factorize(X, C, k)
%DCEBE_factorize Performs a matrix factorization to speed up computation of
%the GCV score when multiple model parameters are evaluated

% init size
p = size(C, 1);
N = size(X, 1);

%%% Step 1: Decompose C
[Gamma, Delta] = eig(C);

% revert order
Gamma = Gamma(:,end:-1:1);
Delta = Delta(end:-1:1,end:-1:1);
delta = diag(Delta);
% extraxt positive eigenvalues (r is the rank)
r = min(sum(delta > sqrt(eps)), p-k);
%r = p-k;
delta_1 = delta(1:r);
%delta_1(delta_1<0) = 0;

delta_1 = max(0, delta_1); % saves negative values due to numerical errors
P_1 = Gamma * diag([delta_1.^(-1/2); ones(p -r ,1)]);  % P_1^T C P_1 = [I 0; 0 0]
X_2 = X * P_1;
X_2_1 = X_2(:,1:r);
X_2_2 = X_2(:,r+1:end);

%%% Step 2
B = (X_2_2' * X_2_2)^(-1);
ReBSqrt = real(B^(1/2)); % can be complex due to machine arithmetic
P_2 = [eye(r), zeros(r, p -r); -B * X_2_2' * X_2_1, ReBSqrt];
X_3 = X_2 * P_2;

%%% Step 3: SVD of X_3
[U,S,~] = svd(X_3);

% compute the indices of the p-q 1-eigenvalues
%q = r;  % (for smoothing splines, see paper by Kent)
q = max(1, min(sum(diag(S) > eps), r));

[~,min_perm] = sort(abs(diag(S)-1));
idx_ones = min_perm(1:p-q);

% resort SVD result such that 1 eigenvalues are in last position
permutation = [1:min(idx_ones)-1, max(idx_ones)+1:p, idx_ones(:)'];
permutation2 = [permutation, numel(permutation)+1:N];
U = U(:,permutation2);
S = S(permutation2,permutation);
S_diag = diag(S);

% copmute d: S_diag = sqrt(d), and take first q elements
d = S_diag(1:q).^2;

end

