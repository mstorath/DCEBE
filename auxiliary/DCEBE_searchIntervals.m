function interval = DCEBE_searchIntervals(y, method)
% determines a reasonable search interval for a signal y

[N, M] = size(y);
interval = zeros(2, M);

minb = 2;
maxb = N/2;

switch method
    % search interval between min and max of the signal
    case 'min-max'
        for m=1:M
            [~,lb] = min(y(:,m)); % left bound
            [~,rb] = max(y(:,m)); % right bound
            interval(:, m) = [max(minb,lb), min(maxb,rb)];
        end
    
    % search interval between min and mean between min and max of the signal, 
    % and left hand side extended by a few indices
    case 'min-max-refined'
        for m=1:M
            [~,lb] = min(y(:,m)); % left bound
            [~,rb] = max(y(:,m)); % max right bound
            tr = mean([y(lb,m),y(rb,m)]); 
            [~,rb_aux] = min(abs(y(lb:rb,m) - tr));
            rb = rb_aux + lb -1;
            lb = lb - 4; % give left slightly more space
            interval(:, m) = [max(minb,lb), min(maxb,rb)];
        end
    
    otherwise
        error('This method does not exist.')
end

