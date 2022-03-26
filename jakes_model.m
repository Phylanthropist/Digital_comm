            % __________________________________________
            % MODIFICATION OF JAKES MODEL (DENTS MODEL)
            % __________________________________________

% N_equal_strenght_Ray =  Number of equal strength signal arriving at the
% receiver, using a multiple of two will simulate fine

% carrier_freq =  Carrier frequency of the signal from the transmitter. 
% Normally 800 or 1900 MHz for mobile comms

% velocity_of_mobile = velocity of the mobile that is to receive the signal
% expressed in meters per second

% speed_of_light =  3*10^8 m/s
% SymbolRate - scalar power of 2 and is in kilo-symbols-per-sec is used to
% specify what should be the transmission data rate. Slower rates will
% provide slowly fading channels. Normal voice and soem data rates are
% 64-256 ksps

% Length of the generated signal seqence 
% NumberWave = Number of wave form to be selected in order not to use the
% entire wave generated (these waves are uncorrelated since they are formed
% using the Hadamard matrix


function [omega_mTau,T_k] = jakes_model(N_equal_strenght_Ray, carrier_freq, velocity_of_mobile, speed_of_light, symbol_rate, Length_generated_signal, NumberWave)

   No = N_equal_strenght_Ray/4;  % No =  is the number of oscillators
   K = sqrt(2/No); % Normalization factor
   n = 1:No;
   arrival_angle = (2*pi*(n-0.5))/N_equal_strenght_Ray;
   Beta = (pi*n)/No;
   Maximum_doppler_shift = (2*pi*carrier_freq*velocity_of_mobile)/ speed_of_light;
   A_j = hadamard(No);
   rng('default')
   theta_n = rand(1, length(n))*2*pi; % Random oscillator phase 
   t=(1/(symbol_rate*1000):1/(symbol_rate*1000):1/(symbol_rate*1000) * Length_generated_signal);

   % APPLYING VECTORIZED ALGEBRA 
   T_k = K .* (A_j(1:NumberWave, n).* (exp(-1i*Beta)))* cos(((Maximum_doppler_shift.* t).' * cos(arrival_angle)) + theta_n).';
   omega_mTau = (1/(symbol_rate*1000)) * (Maximum_doppler_shift/(2*pi));
   [C, Lags] = autocorr(T_k(1,:), length(t)-1);
   tau = Lags * omega_mTau;
   plot(tau, C, '-b');
   xlabel('Normalized times'); ylabel('Normalized Autocorrelation'); 
   title('Autocorrelation of the first waveform k=1');
   grid on

   % _____________________________________________________________________
   % THEORETICAL COMPUTATION OF AUTOCORRELATION USING THE BESSEL FUNCTION
   % _____________________________________________________________________
  
   hold on
   z = Maximum_doppler_shift*tau;
   r = besselj(0,z);
   plot(tau, r, '+r');
end


