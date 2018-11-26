function y_hat = DCEBE_hat_fun( alpha, y, X, nabla )
% returns the model estimate y_hat of y

% Naive approach can get ill conditioned for large signals: y_hat = X * ((XTX + alpha * NTN) \ (X' * y));
% so we solve problem in a stable way using least squares solver based on
% QR decomposition:
y_hat = X * ([X; sqrt(alpha) * nabla] \ [y; zeros(size(nabla, 1), 1)]);

end

