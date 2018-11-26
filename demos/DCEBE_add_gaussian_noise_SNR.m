function [y_noisy, sigma] = add_gaussian_noise_SNR(y, SNR)
%ADD_GAUSSIAN_NOISE_SNR Adds Gaussian noise such that the output signal has
%a signal to noise ratio of SNR (here SNR = max(signal) / std(noise))

if SNR == Inf
    y_noisy = y;
    sigma = 0;
else
    [~,M] = size(y);
    y_noisy = zeros(size(y));
    sigma = zeros(M,1);
    
    for m=1:M
        signal = y(:, m);
        signal_max = max(signal);
        noise = randn(size(signal));
        sigma(m) = signal_max / (SNR * std(noise) ) ;
        y_noisy(:, m) = signal + sigma(m) * noise;
    end
    
end

