function [BAT, output] = DCEBE_estimateBAT( y, varargin )
% The method estimates the bolus arrival time (BAT) of DCE-MRI signals.
% It is particularly intended to work with signals that do not have a fast upslope
% as it is often the case for data of small animals.
% The proposed method employs a spline-based approximation model.
% Parameter estimation is done by generalized cross validation.
% 
% More details are described in the paper
% A. Bendinger, C. Debus, C. Glowa, C. Karger, J. Peter, M. Storath. 
% Bolus arrival time estimation in dynamic contrast-enhanced MRI of small animals based on spline models. 2018
%
% Input
% y: NxM matrix, M DCE signals of length N
%
% Optional input
% 'search_interval': search interval for coarse search of BAT
% 'coarse_res': time resolution for coarse search (0.25 corresponds to 4 samples between two samples)
% 'beta_coarse_search': coarse search in stiffness parameter alpha = beta^(2*k)
% 'orders': tested orders of the spline
% 'verbosity': quantity of information shown during computation (0 no output, 1 a bit, 2 more)
% 'common_BAT': set true if a common BAT parameter is assumed for all signals
% 'solver_type': 'fminunc' for Newton method, 'fminsearch' for derivative free method
%
% Output
% BAT: estimated BAT (for each signal)
%
% Optional output
% output: struct of additional outputs
% 
% Copyright 2018 Alina Bendinger, Martin Storath

% init
[N, M] = size(y); % (M DCE signals each of length N)

% parse optional parameters
ip = inputParser;
ip.addParameter('search_interval', []);     % search interval for coarse search of BAT
ip.addParameter('coarse_res', 0.25);        % time resolution for coarse search (0.25 corresponds to 4 samples between two samples)
ip.addParameter('beta_coarse_search', (linspace(1, 25, 50)));   % coarse search in stiffness parameter alpha = beta^(2*k)
ip.addParameter('orders', 3:6);             % tested orders of the spline
ip.addParameter('verbosity', 1);            % show some intermediate output?
ip.addParameter('common_BAT', false);       % estimate a common BAT for all M signals?
ip.addParameter('solver_type', 'fminunc');  % specify solver type ('fminunc' or 'fminsearch')
ip.parse(varargin{:});
input = ip.Results;

% determine search intervals
if isempty(input.search_interval)
    interval_cand = DCEBE_searchIntervals(y, 'min-max-refined');
    input.search_interval(1) = min(interval_cand(1));
    input.search_interval(2) = max(interval_cand(2));
end

% options for optimizer
switch input.solver_type
    case 'fminsearch'
        options_fminsearch =optimset('Display', 'notify', 'MaxFunEvals', 100000);
        solver = @(f, param) fminsearch(f, param, options_fminsearch);
        
    case 'fminunc'
        options_fminunc = optimoptions(@fminunc,...
            'Display','notify',...
            'Algorithm','quasi-newton',...
            'OptimalityTolerance', 1e-6, ...
            'FiniteDifferenceType', 'central', ...
            'MaxIterations', 1000,...
            'MaxFunctionEvaluations', 10000 ...
            );
        solver = @(f, param) fminunc(f, param, options_fminunc);
        
    otherwise
        error('Solver mist be fminsearch or fminunc.')
end

% init
coarse_t = input.search_interval(1) : input.coarse_res : input.search_interval(2);
T = numel(coarse_t);
L = numel(input.orders);

% containers for the optimal estimates
score_opt = Inf(M,1);
BAT = NaN(M,1);
beta_opt = NaN(M,1);
k_opt = NaN(M,1);

% main loop
for l=1:L % loop over orders
    k = input.orders(l);
    
    if input.verbosity > 0
        fprintf('Order %i (%i of %i):\n  Coarse search... \n', k, l, L)
    end
    
    % coarse search
    A = numel(input.beta_coarse_search);
    minBeta = input.beta_coarse_search(1);
    
    % coarse search scores (dim M x T x A)
    score_coarse =  log(DCEBE_GCVscore(y, coarse_t, input.beta_coarse_search, k, false, minBeta));
    
    % begin fine search
    if input.verbosity > 0
        fprintf('  Fine search... \n')
    end
    
    if input.common_BAT % case of single BAT for all signals
        score_coarse_sum = squeeze(mean(score_coarse, 1)); % sum up scores over all signals
        [min_score_coarse, min_idx] = min(score_coarse_sum(:));
        [idx_BAT, idx_beta] = ind2sub([T, A], min_idx);
        param = [coarse_t(idx_BAT), input.beta_coarse_search(idx_beta)];
        f = @(par) squeeze(mean(log(DCEBE_GCVscore(y, par(1), par(2), k, false, minBeta)), 1));
        [param_fine, score_fine] = solver(f, param);
        if score_fine > f(param)
                warning('Fine search not successfull.')
        end
        param = param_fine;
        score = score_fine;
        if score < score_opt(1)
            score_opt(:) = score;
            BAT(:) = param(1);
            beta_opt(:) = param(2);
            k_opt(:) = k;
        end
        
    else % case of multiple BATs for all signals
        for m=1:M
            % determine best coarse search parameters
            [~,min_idx] = min(score_coarse(m, :));
            [idx_BAT, idx_beta] = ind2sub([T, A], min_idx);
            m_param = [coarse_t(idx_BAT), input.beta_coarse_search(idx_beta)];
            f = @(par) log(DCEBE_GCVscore(y(:,m), par(1), par(2), k, false, minBeta));
            
            [m_param_fine, m_score_fine] = solver(f, m_param);
            if m_score_fine > f(m_param)
                warning('Fine search not successfull.')
            end
            m_param = m_param_fine;
            m_score = m_score_fine;
            
            if m_score < score_opt(m)
                score_opt(m) = m_score;
                BAT(m) = m_param(1);
                beta_opt(m) = m_param(2);
                k_opt(m) = k;
            end
            
            if (input.verbosity > 1) || ((input.verbosity > 0) && (mod(m, 50) == 0))
                fprintf('  %i of %i, LogGCV: %f, beta: %f, k_opt %i, CP: %.3f \n', m, M, score_opt(m), beta_opt(m), k_opt(m), BAT(m))
            end
        end
    end
    % end fine search
end

% get some additional output
output.k_opt = k_opt;
output.beta_opt = beta_opt;
output.score_opt = score_opt;
output.cp_int = floor(BAT);
for m=1:M
    beta = beta_opt(m);
    k = k_opt(m);
    b = BAT(m);
    deriv_pattern = DCEBE_make_deriv_pattern(k);
    [X, nablaK] = DCEBE_make_matrix(N, b, k, deriv_pattern);
    output.y_hat(:,m) = DCEBE_hat_fun( beta.^(2*k), y(:,m), X,  nablaK );
    % add the x-axis location of the BAT (replacing its lefthand neighbor)
    output.y_hat_s_x(:,m) = (1:N);
    output.y_hat_s_x(output.cp_int(m), m) = b;
end

end % end of function
