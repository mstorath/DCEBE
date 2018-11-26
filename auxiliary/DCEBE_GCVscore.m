function score = DCEBE_GCVscore(y, t_arr, beta_arr, k, factorize, minBeta )
%DCEBE_GCVscore Computes the GCV scoore for all values in t_arr and beta_arr

% init
tol = 1e-6;
alpha = beta_arr.^(2*k);
I = numel(t_arr);
J = numel(alpha);
[N, M] = size(y);

score = zeros(M, I, J);
rss = zeros(M, I, J);
df =zeros(M, I, J);
deriv_pattern = DCEBE_make_deriv_pattern(k);

% loop over all time points
for i = 1:I
    
    % changepoint index candidate
    t = t_arr(i);
    
    % prevent h from getting too small
    t_frac = t - floor(t);
    h = 1 - t_frac;
    if abs(h) < tol
        t = ceil(t) - tol;
    end
    
    % only t between 1 and N-k-1 is valid
    t_inbound = max(min(t, N-k-1), 1);
    
    [X, nablaK] = DCEBE_make_matrix(N, t_inbound, k, deriv_pattern);
    
    % two variants for computing the GCV score are available
    if factorize
        % compute GCV score by matrix factorization
        % see paper by J. T. Kent, M. Mohammadzadeh
        % "Global optimization of the generalized cross-validation
        % criterion", 2000
        
        NTN = nablaK' * nablaK;
        [U, d, p, q, r] = DCEBE_factorize(X, full(NTN), k);
        
        %%%  apply U' to y
        yStar = U' * y;
        yStar1 = yStar(1:q,:);
        %yStar2 = yStar((q+1):q+r-q,:); % is empty
        yStar2 = [];
        yStar4 = yStar((p+1):end,:);
        z = yStar1;
        w = [yStar2; yStar4];
        w_sum = sum(w.^2, 1); % row vector (size m)
        
        % compute a constant
        %npqr = N - p; % actually %npqr = N - p - q + r; but  q = r
        npqr = N - p - q + r;
        
        % check for imag
        if any(imag(z) ~= 0)
            warning('Imaginary parts occurred. May be due to round-off.')
        end
        
        %%% calculate scores for all alpha
        v = sum( alpha./ (d + alpha), 1) + npqr;
        for m = 1:M
            u = sum( z(:,m).^2 .* alpha.^2 ./ (d + alpha).^2) + w_sum(m);
            score(m,i,:) = N * u ./ v.^2;
        end
        
    else
        % direct variant using QR decomposition, 
        % better conditioned than factorization variant
        for j=1:J
            Z = [speye(N); sparse(size(nablaK, 1), N)];
            hat_matrix = X * ([ X; sqrt(alpha(j)) * nablaK ]  \ Z);
            y_hat = hat_matrix * y;
            tr = trace(hat_matrix);
            df(:,i,j) = tr;
            for m = 1:M
                rss(m,i,j) = sum((y_hat(:,m) - y(:,m)).^2);
                score(m,i,j) = (1/N) * rss(m,i,j)./ (1 - df(m,i,j)/N)^2;
            end
        end
        
    end
end

% add barrier functions to bound the search space to feasible values (alpha > alphaMin, t \in (1, N-k))
score = score + barrier(shiftdim(beta_arr(:), -2), minBeta-1) + barrier(shiftdim(t_arr(:), -1), 1) + barrier(-shiftdim(t_arr(:), -1), -(N-k));

end

% auxiliary barrier function to constrain the search space to a feasible interval  
function y = barrier(x, b)
y = -log( sin((x-b)*pi/2) );
y(x < b) = Inf;
y(x > (b + 1)) = 0;
end

