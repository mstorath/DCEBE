% the script estimates the common bolus arrival time of a series of DCE-MRI
% (estimating a single BAT from multiple signal increases stability)

% load simulated Rat data (high temporal resolution delta_t_hr = 0.25 sec)
load('DCEBE_Rat_ETM.mat')
%load('DCEBE_Rat_2CXM.mat')

curve_type = 2 ; %(1, 2 or 3)
concentration_hr = input_data(:,curve_type);
delta_t_hr = 0.25;
[N,~] = size(concentration_hr);
t_hr = (0:N-1) * delta_t_hr;
BAT_true = (find(abs(diff(concentration_hr)), 1, 'first') - 1) * delta_t_hr; % true BAT is defined as last point with baseline concentration of high resolution signal

% simulate lower temporal resolution of delta_t_lr sec by downsampling
delta_t_lr = 2;
ratio = delta_t_lr/delta_t_hr;
concentration_lr = concentration_hr(1:ratio:end,:);
t_lr = t_hr(1:ratio:end);

% add noise to the signals (M realizations)
rng(123) % random seed for reproducibility
SNR = 10; % signal to noise ratio
M = 20; % number of signals
y = DCEBE_add_gaussian_noise_SNR(repmat(concentration_lr, [1,M]), SNR);

% estimate a common BAT index for all signals simultaneously (columns of y)
BAT_est_idx = DCEBE_estimateBAT(y, 'common_BAT', true);

% convert from index to time
BAT_est = (BAT_est_idx - 1) * delta_t_lr;
   
%% show output
plotopts = {'MarkerSize', 10,'Linewidth', 2};
plot(t_lr, y, '.-')
hold on
plot(BAT_est, zeros(size(BAT_est)), 'xg', plotopts{:}, 'DisplayName', 'Estimated BAT' )
plot(BAT_true, 0, 'ok', plotopts{:}, 'DisplayName', 'True BAT')
xlabel('Time [s]')
ylabel('Concentration')
hold off
legend('show')