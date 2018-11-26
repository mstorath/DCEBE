function deriv_pattern = DCEBE_make_deriv_pattern(k, l, h)
%DCEBE_make_deriv_pattern Computes the finite difference weights for the
%nodes [0, h, h+1, ...h+l-2]

% l is number of elements in derivtive pattern
if ~exist('l', 'var')
    l = k+1;
end

% create derivative pattern for order k
if ~exist('h', 'var')
    h = 1;
end

x = [0, h + (0:l-2)];
deriv_pattern = fdweights(0,x,k);

end

