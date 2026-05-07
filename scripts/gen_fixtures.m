% Generate MATLAB-side reference fixtures for the Python parity tests.
%
% Run from the repo root: `matlab -nodisplay -batch "run('scripts/gen_fixtures.m')"`
%
% Each fixture stores:
%   y                — N x M input signals
%   bat_matlab       — M x 1 BAT estimates (1-based, sample units)
%   k_opt            — M x 1 selected spline orders
%   beta_opt         — M x 1 stiffness parameters
%   score_opt        — M x 1 final log-GCV scores
%   y_hat            — N x M reconstructed signals
%   meta             — string describing args used (re-applied in Python)
%
% RNG seed pinned to 123 to match the existing demos. The fixtures are
% small enough (M = 2..5, N = 30..50) for the parity tests to run in
% under a minute total.

repo_root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(repo_root);
addpath(fullfile(repo_root, 'auxiliary'));
addpath(fullfile(repo_root, 'external', 'fdweights'));
addpath(fullfile(repo_root, 'demos'));

fixture_dir = fullfile(repo_root, 'matlab_fixtures');
if ~exist(fixture_dir, 'dir')
    mkdir(fixture_dir);
end

%% Fixture 1: small synthetic signals, default search interval
fprintf('--- fixture_small ---\n');
rng(123);
N = 30;
M = 2;
true_bat_idx = 12;
y_clean = zeros(N, 1);
y_clean(true_bat_idx:end) = linspace(0, 1, N - true_bat_idx + 1);
y = repmat(y_clean, 1, M);
SNR = 30;
y = y + (1 / SNR) * randn(size(y));

[bat_matlab, output] = DCEBE_estimateBAT(y, 'verbosity', 0);
k_opt = output.k_opt;
beta_opt = output.beta_opt;
score_opt = output.score_opt;
y_hat = output.y_hat;
meta = 'default-search-interval (envelope across signals after bug fix); per-signal mode';
save(fullfile(fixture_dir, 'fixture_small.mat'), ...
    'y', 'bat_matlab', 'k_opt', 'beta_opt', 'score_opt', 'y_hat', 'meta');
fprintf('  saved %s\n', fullfile(fixture_dir, 'fixture_small.mat'));

%% Fixture 2: explicit search interval, per-signal mode
fprintf('--- fixture_explicit_si ---\n');
rng(123);
N = 40;
M = 3;
true_bats = [15, 20, 18];
y = zeros(N, M);
for m = 1:M
    bi = true_bats(m);
    sig = zeros(N, 1);
    sig(bi:end) = linspace(0, 1, N - bi + 1);
    y(:, m) = sig;
end
SNR = 25;
y = y + (1 / SNR) * randn(size(y));

search_interval = [10, 30];
[bat_matlab, output] = DCEBE_estimateBAT(y, ...
    'search_interval', search_interval, 'verbosity', 0);
k_opt = output.k_opt;
beta_opt = output.beta_opt;
score_opt = output.score_opt;
y_hat = output.y_hat;
meta = 'search_interval=[10, 30]; per-signal mode';
save(fullfile(fixture_dir, 'fixture_explicit_si.mat'), ...
    'y', 'bat_matlab', 'k_opt', 'beta_opt', 'score_opt', 'y_hat', ...
    'search_interval', 'meta');
fprintf('  saved %s\n', fullfile(fixture_dir, 'fixture_explicit_si.mat'));

%% Fixture 3: common BAT mode
fprintf('--- fixture_common ---\n');
rng(123);
N = 40;
M = 3;
true_bat = 18;
y_clean = zeros(N, 1);
y_clean(true_bat:end) = linspace(0, 1, N - true_bat + 1);
y = repmat(y_clean, 1, M);
SNR = 25;
y = y + (1 / SNR) * randn(size(y));

search_interval = [10, 30];
[bat_matlab, output] = DCEBE_estimateBAT(y, ...
    'search_interval', search_interval, 'common_BAT', true, 'verbosity', 0);
k_opt = output.k_opt;
beta_opt = output.beta_opt;
score_opt = output.score_opt;
y_hat = output.y_hat;
meta = 'search_interval=[10, 30]; common_bat=true';
save(fullfile(fixture_dir, 'fixture_common.mat'), ...
    'y', 'bat_matlab', 'k_opt', 'beta_opt', 'score_opt', 'y_hat', ...
    'search_interval', 'meta');
fprintf('  saved %s\n', fullfile(fixture_dir, 'fixture_common.mat'));

fprintf('\nAll fixtures generated in %s\n', fixture_dir);
