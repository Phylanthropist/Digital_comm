clc, clear var, close all;
% ____________________________________________________________________
% ENTER THE NUMBER OF BITS, REQUIRED PHASE OFFSET AND FREQUENCY OFFSET
% ____________________________________________________________________
N = input("Enter the number of  input bit needed to be generated = ");
disp('Binary Input Information at Transmitter: ');
disp(N);
% Phase offset input
Phase_off = input("Enter the phase offset = ");
disp('Phase offset: ');
disp(Phase_off);
% frequency offset input
frequency_off = input("Enter the frequency offset = ");
disp('frequency offset: ' );
disp(frequency_off);


BPSK_constellation(N, Phase_off, frequency_off);
QPSK_constellation(N, Phase_off, frequency_off);
QAM_constellation(N, Phase_off, frequency_off)

% generate the data and create the BPSK modulator
function BPSK_constellation(Number_of_input, phaseOffset, frequencyOffset)
    % data must be a column vector
    data = transpose(randn(1,Number_of_input)>0);
    bpskModulator = comm.BPSKModulator;
    BPSK_data = bpskModulator(data);
    awgn_BPSK_data = awgn(BPSK_data, 10, 'measured'); % adding white guassian noise to the signal.
    phase_Offset = comm.PhaseFrequencyOffset('PhaseOffset', phaseOffset);
    frequency_Offset = comm.PhaseNoise("FrequencyOffset",frequencyOffset, 'Level', -50);
    BPSK_data_po = phase_Offset(awgn_BPSK_data);
    BPSK_data_fo = frequency_Offset(BPSK_data_po);
    scatterplot(BPSK_data_fo)
    title('BPSK constellation');
    hold on;
end

function QPSK_constellation(Number_of_input, phaseOffset, frequencyOffset)
    % data must be a column vector
    data = (randn(1,Number_of_input)>0).';
    qpskModulator = comm.QPSKModulator('BitInput',true);
    QPSK_data = qpskModulator(data);
    awgn_QPSK_data = awgn(QPSK_data, 10, 'measured');
    phase_freq_Offset_1 = comm.PhaseFrequencyOffset('PhaseOffset', phaseOffset);
    frequency_Offset = comm.PhaseNoise("FrequencyOffset",frequencyOffset, 'Level', -50);
    QPSK_data_po = phase_freq_Offset_1(awgn_QPSK_data);
    QPSK_data_fo= frequency_Offset(QPSK_data_po);
    scatterplot(QPSK_data_fo)
    title('QPSK constellation');
    hold on;
end

function QAM_constellation(Number_of_input, phaseOffset, frequencyOffset)
    M = 16;
    data = randi([0 1], Number_of_input*log2(M),1);
    Qam_mod = qammod(data, M, 'Inputtype', 'bit', 'UnitAveragePower', true );
    phaseoff_QAM = comm.PhaseFrequencyOffset('PhaseOffset', phaseOffset);
    frequency_Offset = comm.PhaseNoise("FrequencyOffset",frequencyOffset, 'Level', -50);
    awgn_Qam_data = awgn(Qam_mod, 10, 'measured');
    QAM_mod_po = phaseoff_QAM(awgn_Qam_data);
    QAM_mod_fo = frequency_Offset(QAM_mod_po);
    scatterplot(QAM_mod_fo)
    title('QAM constellation');
    hold on;
end