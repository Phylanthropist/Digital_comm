% Note: n is the ray at a particular instance is moves from n to N
function T_k = jakes_model(N_equal_strenght_Ray, carrier_freq, velocity_of_mobile, speed_of_light, Number_of_wave )
   No = N_equal_strenght_Ray/4;  % where No is the number of oscillators
   K = sqrt(2/No); % Normalization factor
   n = (1:No);
   arrival_angle = (2*pi*(n-0.5))/N_equal_strenght_Ray;
   Beta = pi*n/No;
   Maximum_doppler_shift = (2*pi*carrier_freq*velocity_of_mobile)/ speed_of_light;
   doppler_shift = Maximum_doppler_shift * cos(arrival_angle);
   A_j = hadamard(No); % The hadamard matrix use for the auto and cross correlation
   % The angle theta has to be randomize in order to produce various waveform.
   rng('shuffle')
   theta_n = rand(1, length(n))*2*pi; % Oscillator phase
   % USING THE VECTORIZED ALGEBRA FORMAT
   T_k = K* (A_j(Number_of_wave, 1:No) .* (exp(-1i*Beta).*cos(doppler_shift + theta_n)));
end


