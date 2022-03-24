tic
close all; clear var; clc;

% INPUTS AND CONVERSION
Number_of_input_bit = 10^6;
EbN0_dB = -8:2:10;
% ____________________
% BPSK's BER ANALYSIS
% ____________________


Simulated_BER_BPSK_AWGN = zeros(1, length(EbN0_dB));
Simulated_BER_BPSK_RAY = zeros(1, length(EbN0_dB));

for k = 1:length(EbN0_dB)
    data = rand(1, Number_of_input_bit)>=0.5; % A logical operation. 
    transmitted_BPSK = 2 * data -1; % BPSK mapping (antipodal points)
    % Convert the noise in dB to linear scale for every iteration
    h = 1/sqrt(2)* (randn(1, Number_of_input_bit) + 1i* randn(1, Number_of_input_bit));
    noise = 1/sqrt(2)* (randn(1, Number_of_input_bit) + 1i* randn(1, Number_of_input_bit));
    Generated_noise = noise * 10^(-EbN0_dB(k)/20); % for every iteration the noise changes. This value need not be static
    Recieved_BPSK_AWGN = Generated_noise + transmitted_BPSK;
    Recieved_BPSK_RAY = (h.* transmitted_BPSK) + Generated_noise;
    Recieved_BPSK_RAY_1 = Recieved_BPSK_RAY./h;

    % ____________________________________________________________________
    % THRESHOLD DETECTOR OF RECEIVED BIT FOR ERROR/ BIT ERROR CALCULATION
    % ____________________________________________________________________
    estimated_bit_AWGN = real(Recieved_BPSK_AWGN) >=0;
    estimated_bit_RAY = real(Recieved_BPSK_RAY_1) >0;
    
    Simulated_BER_BPSK_AWGN(k) = sum(xor(data, estimated_bit_AWGN))/length(data);
    Simulated_BER_BPSK_RAY(k) = sum(xor(data, estimated_bit_RAY))/length(data);
end

% _____________________
% QPSK's BER ANALYSIS
% _____________________
QPSK_constant = 2; % 1 for BPSK, 2 for QPSK
QPSK_M = 2^QPSK_constant;
Rm_QPSK = log2(QPSK_M); % bit per symbol, here it is = 2
Rc = 1;

BER_QPSK_AWGN = zeros(1, length(EbN0_dB));
BER_QPSK_RAYLEIGH = zeros(1, length(EbN0_dB));

for k = 1:length(EbN0_dB)
    data = randn(1, Number_of_input_bit)>=0; % A logical operation. 
    inphase_data = data(1:2:end);
    quadrature_data = data(2:2:end);
    transmitted_QPSK = sqrt(1/2)*(1i*(2*quadrature_data-1)+(2*inphase_data-1)); % QPSK mapping, generates a pair of real and imaginary components
    % Convert the noise in dB to linear scale for every iteration
    Noise_variance = sqrt(1./((2 * Rm_QPSK * Rc) * 10.^(EbN0_dB(k)/10)));
    Generated_noise = Noise_variance * (randn(1, length(transmitted_QPSK)) + 1i*randn(1, length(transmitted_QPSK)));
    h_RAY = 1/sqrt(2) *  (randn(1, length(transmitted_QPSK)) + 1i*randn(1, length(transmitted_QPSK)));
    % Just like the transmitted data, the noise should also have the inphase and quadrature points
    Recieved_QPSK_AWGN = transmitted_QPSK + Generated_noise;
    Recieved_QPSK_RAYLEIGH = h_RAY.* transmitted_QPSK + Generated_noise;
    Recieved_QPSK_RAYLEIGH_1 = Recieved_QPSK_RAYLEIGH./h_RAY;
    % ____________________________________________________________________
    % THRESHOLD DETECTOR OF RECEIVED BIT FOR ERROR/ BIT ERROR CALCULATION
    % ____________________________________________________________________
    detected_inphase_AWGN = real(Recieved_QPSK_AWGN)>=0;
    detected_quadrature_AWGN = imag(Recieved_QPSK_AWGN)>=0;
       
    detected_inphase_RAYLEIGH = real(Recieved_QPSK_RAYLEIGH_1)>=0;
    detected_quadrature_RAYLEIGH = imag(Recieved_QPSK_RAYLEIGH_1)>=0;
   
 
    Received_symbol1_AWGN = detected_inphase_AWGN + 1i * detected_quadrature_AWGN;
    estimated_bits_AWGN = reshape([detected_inphase_AWGN; detected_quadrature_AWGN],1,[]);
    
    Received_symbol1_RAYLEIGH = detected_inphase_RAYLEIGH + 1i * detected_quadrature_RAYLEIGH;
    estimated_bits_RAYLEIGH = reshape([detected_inphase_RAYLEIGH; detected_quadrature_RAYLEIGH],1,[]);

    BER_QPSK_AWGN(k) = sum(xor(data, estimated_bits_AWGN))/length(data);
    BER_QPSK_RAYLEIGH(k) = sum(xor(data, estimated_bits_RAYLEIGH))/length(data);
end


% _____________________
% 16QAM's BER ANALYSIS\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% _____________________
QAM_M = 16; % for 16-QAM modulation (16 constellation)

%      A           B     3 +      E            F
%                          +
%                          +
%                          +
%                        1 +
%      C           D       +      G            H
%                          +
%  ++++++++++++++++++++++++++++++++++++++++++++++++
%     -3           -1      +      1            3
%                          +
%      I           J    -1 +      M            N
%                          +
%                          +
%                          +
%      K           L    -3 +      O            P
% MQAM = +-(Square_root(M)), therefore array from constellation
Norm_factor = 1/sqrt(10);
% creating the real and imaginary points of the constellation
Real = [(-1*sqrt(QAM_M)/2-1):2:-1 1:2:sqrt(QAM_M)-1]; % [-3,-1,1,3] constellation points
Imaginary = [(-1*sqrt(QAM_M)/2-1):2:-1 1:2:sqrt(QAM_M)-1];
KK = log2(QAM_M);
% A symbol consist of 4 bits.
EsN0_dB = EbN0_dB + 10*log10(KK);

% Mapping to gray code
ref = (0:k-1);
map = bitxor(ref, floor(ref/2));
[value, index] = sort(map);
bitERR_AWGN = zeros(1, length(EbN0_dB));
bitERR_RAYLEIGH = zeros(1, length(EbN0_dB));

for i = 1:length(EbN0_dB)
    % ________________
    % MODULATION
    % ________________

    % Generate the data and group as a symbol of four
    % divide the generated bits into two, and apply gray coding
    Raw_data = rand(1, Number_of_input_bit*KK, 1)>0.5;
    Raw_data = reshape(Raw_data,KK, Number_of_input_bit).';
    bin2Dec = ones(Number_of_input_bit,1)*(2.^((KK/2-1):-1:0));
    % Real part computation
    Raw_data_real = Raw_data(:,(1:KK/2));
    bin2Dec_real = sum(Raw_data_real.*bin2Dec,2);
    Dec2Gray_real = bitxor(bin2Dec_real, floor(bin2Dec_real/2));
    % Imaginary part computation
    Raw_data_Imag = Raw_data(:,(KK/2+1:KK));
    bin2Dec_Imag = sum(Raw_data_Imag.*bin2Dec,2);
    Dec2Gray_Imag = bitxor(bin2Dec_Imag, floor(bin2Dec_Imag/2));
    % Map the Real and Imaginary Gray code into the constellation
    Map_Re = Real(Dec2Gray_real+1);
    Map_Imag = Imaginary(Dec2Gray_Imag+1);
    map_data = Map_Re + 1i*Map_Imag;
    Norm_map_data = Norm_factor * map_data;

    Noise = 1/sqrt(2)*(randn(1,Number_of_input_bit) + 1i*(randn(1, Number_of_input_bit)));
    Transmitted_data_AWGN = Norm_map_data + 10^(-EsN0_dB(i)/20)* Noise;

    h_RAY_1 = 1/sqrt(2) * (randn(1, length(Raw_data)) + 1i*randn(1, length(Raw_data)));
    Transmitted_data_RAYLEIGH = h_RAY_1.* Norm_map_data + 10^(-EsN0_dB(i)/20)* Noise;
    Transmitted_data_RAYLEIGH_1 = Transmitted_data_RAYLEIGH./h_RAY_1;

    % ________________
    % DEMODULATION
    % ________________

    Transmitted_data_Real_AWGN = real(Transmitted_data_AWGN)/Norm_factor;
    Transmitted_data_Imag_AWGN = imag(Transmitted_data_AWGN)/Norm_factor;

    Transmitted_data_Real_RAYLEIGH = real(Transmitted_data_RAYLEIGH_1)/Norm_factor;
    Transmitted_data_Imag_RAYLEIGH = imag(Transmitted_data_RAYLEIGH_1)/Norm_factor;


    % Real part demodulation computation for AWGN
    Tx_Real_AWGN = 2* floor(Transmitted_data_Real_AWGN/2)+1;
    Tx_Real_AWGN(Tx_Real_AWGN>max(Real)) = max(Real);
    Tx_Real_AWGN(Tx_Real_AWGN<min(Real)) = min(Real);
    % Real part demodulation computation for AWGN
    Tx_imag_AWGN = 2*floor(Transmitted_data_Imag_AWGN/2) +1;
    Tx_imag_AWGN(Tx_imag_AWGN>max(Imaginary)) = max(Imaginary);
    Tx_imag_AWGN(Tx_imag_AWGN<min(Imaginary)) = min(Imaginary);

    % Real part demodulation computation for RAYLEIGH
    Tx_Real_RAYLEIGH = 2* floor(Transmitted_data_Real_RAYLEIGH/2)+1;
    Tx_Real_RAYLEIGH(Tx_Real_RAYLEIGH>max(Real)) = max(Real);
    Tx_Real_RAYLEIGH(Tx_Real_RAYLEIGH<min(Real)) = min(Real);
    % Real part demodulation computation for RAYLEIGH
    Tx_imag_RAYLEIGH = 2*floor(Transmitted_data_Imag_RAYLEIGH/2) +1;
    Tx_imag_RAYLEIGH(Tx_imag_RAYLEIGH>max(Imaginary)) = max(Imaginary);
    Tx_imag_RAYLEIGH(Tx_imag_RAYLEIGH<min(Imaginary)) = min(Imaginary);


    % converting constellation to decimal for AWGN
    Tx_Real_Dec_AWGN = index(floor((Tx_Real_AWGN+4)/2+1))-1; % 4 is added because the lowest number is -3, 
    Tx_Imag_Dec_AWGN = index(floor((Tx_imag_AWGN+4)/2+1))-1;
    % converting decimal constellation to binary for AWGN
    Tx_Real_bin_AWGN = dec2bin(Tx_Real_Dec_AWGN, KK/2) - '0';
    Tx_imag_bin_AWGN = dec2bin(Tx_Imag_Dec_AWGN, KK/2) - '0';

    % converting constellation to decimal for RAYLEIGH
    Tx_Real_Dec_RAYLEIGH = index(floor((Tx_Real_RAYLEIGH+4)/2+1))-1; % 4 is added because the lowest number is -3, 
    Tx_Imag_Dec_RAYLEIGH = index(floor((Tx_imag_RAYLEIGH+4)/2+1))-1;
    % converting decimal constellation to binary for RAYLEIGH
    Tx_Real_bin_RAYLEIGH = dec2bin(Tx_Real_Dec_RAYLEIGH, KK/2) - '0';
    Tx_imag_bin_RAYLEIGH = dec2bin(Tx_Imag_Dec_RAYLEIGH, KK/2) - '0';

    % finally demodualated bit, inorder to calculate the bit error
    bitERR_AWGN(i) = size(find(Raw_data_real - Tx_Real_bin_AWGN ),1) + size(find(Raw_data_Imag - Tx_imag_bin_AWGN ),1);
    bitERR_RAYLEIGH(i) = size(find(Raw_data_real - Tx_Real_bin_RAYLEIGH ),1) + size(find(Raw_data_Imag - Tx_imag_bin_RAYLEIGH ),1);
end

BER_16QAM_AWGN = bitERR_AWGN/(Number_of_input_bit * KK);
BER_16QAM_RAYLEIGH = bitERR_RAYLEIGH/(Number_of_input_bit * KK);

% ____________________
%         GRAPH PLOT
% ____________________
subplot(2,1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATED BPSK OVER AWGN  %%%%%%%%%%%%%%%%%%
semilogy(EbN0_dB, Simulated_BER_BPSK_AWGN, 'b--', 'LineWidth',1);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATED QPSK OVER AWGN %%%%%%%%%%%%%%%%%%
semilogy(EbN0_dB, BER_QPSK_AWGN, 'K--', 'LineWidth',1);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATED 16QAM OVER AWGN %%%%%%%%%%%%%%%%%%
semilogy(EbN0_dB, BER_16QAM_AWGN, 'r-', 'LineWidth',1);
hold on;
% ___________________
% Theoretical BERs AWGN
% ___________________
%%%%%%%%%%%%%%%%%%%%%%%% BPSK OVER AWGN %%%%%%%%%%%%%%%%%%%%%%%%%%%
semilogy(EbN0_dB, 0.5*(erfc(sqrt(10.^(EbN0_dB/10)))), 'y*', 'LineWidth',0.9);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%% QPSK OVER AWGN %%%%%%%%%%%%%%%%%%%%%%%%%%%
semilogy(EbN0_dB, 0.5*(erfc(sqrt(10.^(EbN0_dB/10)))), 'g*', 'LineWidth',0.9);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%% 16QAM OVER AWGN %%%%%%%%%%%%%%%%%%%%%%%%%%%
semilogy(EbN0_dB, (1/KK)*3/2*erfc(sqrt(KK*0.1*(10.^(EbN0_dB/10)))), 'bs', 'MarkerFacecolor', 'b', 'LineWidth',0.9);

title("BPSK, QPSK, 16QAM's BER over AWGN channel");
xlabel('Normalized Signal to Noise ratio (Eb/N0)');
ylabel('Bit Error rate (BER) in dB');

legend('Simluated BPSK ','Simluated QPSK', 'Simluated 16QAM','Theoretical BPSK', ...
    'Theoretical QPSK', 'Theoretical 16QAM', Location='southwest');


subplot(2,1,2)
%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATED BPSK OVER RAYLEIGH  %%%%%%%%%%%%%%%%%%
semilogy(EbN0_dB, Simulated_BER_BPSK_RAY, 'm--', 'LineWidth',1);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATED QPSK OVER RAYLEIGH  %%%%%%%%%%%%%%%%%%
semilogy(EbN0_dB, BER_QPSK_RAYLEIGH, 'yo', 'LineWidth',2);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATED 16QAM OVER RAYLEIGH  %%%%%%%%%%%%%%%%%%
semilogy(EbN0_dB, BER_16QAM_RAYLEIGH, 'c*-', 'LineWidth',1);
hold on;

title("BPSK, QPSK, 16QAM's BER over RAYLEIGH channel");
xlabel('Normalized Signal to Noise ratio (Eb/N0)');
ylabel('Bit Error rate (BER) in dB');
% _________________________
% Theoretical BERs RAYLEIGH
% _________________________

%%%%%%%%%%%%%%%%%%%%%%%% BPSK OVER RAYLEIGH %%%%%%%%%%%%%%%%%%%%%%%%
semilogy(EbN0_dB, 0.5*(1-sqrt(10.^(EbN0_dB/10)./(1+(10.^(EbN0_dB/10))))), 'k*', 'LineWidth',0.9);
hold on; 
%%%%%%%%%%%%%%%%%%%%%%%% QPSK OVER RAYLEIGH %%%%%%%%%%%%%%%%%%%%%%%%
semilogy(EbN0_dB, 0.5.*(1-sqrt(10.^(EbN0_dB/10)./(10.^(EbN0_dB/10)+1))), 'gd', 'LineWidth',1.5);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%% 16QAM OVER RAYLEIGH %%%%%%%%%%%%%%%%%%%%%%%
theoretical_QAM_RAY = 3/8 * ( 1 - sqrt(2/5*10.^(EbN0_dB/10)*log2(QAM_M)/log2(QAM_M)./(1+2/5*10.^(EbN0_dB/10)*log2(QAM_M)/log2(QAM_M))) );
semilogy(EbN0_dB, theoretical_QAM_RAY, 'o', 'MarkerFacecolor', 'r', 'LineWidth',0.9);

title("BPSK, QPSK, 16QAM's BER over RAYLEIGH channel");
xlabel('Normalized Signal to Noise ratio (Eb/N0)');
ylabel('Bit Error rate (BER) in dB');

legend('Simluated BPSK ','Simluated QPSK', 'Simluated 16QAM','Theoretical BPSK', ...
    'Theoretical QPSK', 'Theoretical 16QAM', Location='southwest');
hold on; 
axis tight
grid on
toc
