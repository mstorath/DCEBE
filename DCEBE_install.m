% show license
disp('----------')
lic = fileread('LICENSE');
disp(lic);
disp('----------')
clear('lic')

disp('Setting paths for DCEBE...')
cd(fileparts(which(mfilename)))
addpath('.')
addpath('auxiliary')
addpath(fullfile('external', 'fdweights'));
savepath;
disp('Done.')
