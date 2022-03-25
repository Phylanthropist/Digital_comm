% Note: n is the ray at a particular instance is moves from n to N
function [omega_mTau,T_k] = jakes_model(N_equal_strenght_Ray, carrier_freq, velocity_of_mobile, speed_of_light, symbol_rate, Length_generated_signal, NumberWave)
   No = N_equal_strenght_Ray/4;  % where No is the number of oscillators
   K = sqrt(2/No); % Normalization factor
   n = (1:No);
   arrival_angle = (2*pi*(n-0.5))/N_equal_strenght_Ray;
   Beta = (pi*n)/No;
   Maximum_doppler_shift = (2*pi*carrier_freq*velocity_of_mobile)/ speed_of_light;
   A_j = hadamard(No); % The hadamard matrix use for the auto and cross correlation
   % The angle theta has to be randomize in order to produce various waveform.
   rng('default')
   theta_n = rand(1, length(n))*2*pi; % Oscillator phase
   t=(1/(symbol_rate*1000):1/(symbol_rate*1000):1/(symbol_rate*1000) * Length_generated_signal);
   % USING THE VECTORIZED ALGEBRA FORMAT
   T_k = K .* (A_j(1:NumberWave, n).* (exp(-1i*Beta)))* cos(((Maximum_doppler_shift.* t).' * cos(arrival_angle)) + theta_n).';% we use the number of  wave when we want to limit the number of wave output.
   omega_mTau = (1/(symbol_rate*1000)) * (Maximum_doppler_shift/2*pi);
   [C, Lags] = autocorr(T_k(1,:), length(t)-1);
   length(t);
   tau = Lags * omega_mTau;
   plot(tau, C, 'ob');
   xlabel('Normalized Time Delay'); ylabel('Normalized Autocorrelation'); 
   title('Autocorrelation of the first waveform k=1');
   grid on
end


