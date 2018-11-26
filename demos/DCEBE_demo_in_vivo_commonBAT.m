% load in vivo rat data (94 signals of 94 voxels withing one tumor)
load('DCEBE_Rat_in_vivo.mat')
y = input_data;

% estimation time with the default search interval can take up to 30 min
si = [];
% to lower estimation time, restrict search space to a reasonable interval, 
% e.g. to si = [40, 60] by uncommenting the following line
%si = [40, 60];

% estimate a common BAT index for all signals (columns of y)
tic
BAT_est_idx = DCEBE_estimateBAT(y, 'search_interval', si, 'common_BAT', true);
toc

% convert from index to time
delta_t = 0.75;
BAT_est = (BAT_est_idx - 1) * delta_t;
t = (0:size(y, 1)-1) * delta_t;   

%% show output
plotopts = {'MarkerSize', 10,'Linewidth', 2};
plot(t, y, '.-')
hold on
plot(BAT_est, zeros(size(BAT_est)), 'xg', plotopts{:}, 'DisplayName', 'Estimated BATs' )
xlabel('Time [s]')
ylabel('Concentration')
hold off
legend('show')