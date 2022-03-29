%%%%%%%%%%%%%%%%%%%%%%% JAKES MODEL REVISITED (DENT MODEL) %%%%%%%%%%%%%%%%
%                  이름: ONYEKWELU MICHAEL CHISOM ( 마이클 치솜) 
%                  학번: 2022151222
%                  전공: 융합전자공학과
%                  과정 제목: 무선 통신
%                  주제: FADING CHANNEL (MULTIPLE UNCORRELATED SIGNAL)
%                  교수: 문의찬 교수님
%                  날짜: 2022년 03월 30일
%                  과정 번호: ENE 1015

%%%%%%%%%%%%%%%%%%%%%% DECLARATION OF VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%
N_equal_strength_Ray = 1024;
Number_of_wave  = 2;
sampling_interval = 383.5e-6;
Number_of_samples = 1e6;
Number_of_oscillators = N_equal_strength_Ray/4;  
Normalization_factor = sqrt(2/Number_of_oscillators); 
n = 1:Number_of_oscillators;
arrival_angle = (2*pi*(n-0.5))/N_equal_strength_Ray;
Beta = (pi*n)/Number_of_oscillators;
Maximum_doppler_shift = 2*pi*83;
Hadamard_matrices = hadamard(Number_of_oscillators);
theta_n = rand(1, length(n))*2*pi; % Random oscillator phase 
t=0:sampling_interval:(Number_of_samples-1)*sampling_interval;
% APPLYING VECTORIZED ALGEBRA 
T_k = Normalization_factor .* Hadamard_matrices(Number_of_wave,:).* (exp(-1i*Beta))* transpose(cos((transpose(Maximum_doppler_shift.* t) * cos(arrival_angle)) + theta_n));

%%%%%%%%%%%%%%%%%%%% PLOTS AND DISPLAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
[Auto_Corr_Value, tau] = xcorr(T_k, 'normalized');
Auto_Corr_Value = real(Auto_Corr_Value);
idx0 = find(tau==0);
tau = tau * sampling_interval; 
plot(tau, Auto_Corr_Value, 'ob');
xlim([tau(idx0) tau(idx0+50)]);
xlabel('\omega_m\tau'); 
ylabel('Re[Autocorrelation]'); 
title('Autocorrelation graph for a fading channel');
grid on
hold on
J0 = besselj(0, Maximum_doppler_shift*tau);
plot(tau, J0, 'r')
hold off
legend('Simulated', 'theoretical')

figure(2)
Abs_T = abs(T_k);
Histogram = histogram(Abs_T, 'Normalization','pdf', 'DisplayStyle','bar', 'EdgeColor','auto');
hold on;
ax = gca;
x = linspace(ax.XLim(1), ax.XLim(2), 1000);
plot(x, Histogram, 'LineWidth', 2)
xlabel('x'); 
ylabel('probability density'); 
title('Probability density function (PDF)');
axis tight
grid on




