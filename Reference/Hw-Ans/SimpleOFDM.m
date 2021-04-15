%% ************** Preparation part ********************
clear all; clc;
% system parameters
ml = 2;                      % Modulation level: 2--4QAM; 4--16QAM; 6--64QAM
NormFactor = sqrt(2/3*(ml.^2-1));
fftlen = 64;
bits_per_sym = fftlen*ml;
NumSyms = 50;
TotalNumBits = bits_per_sym*NumSyms;
%% ************** Loop start***************************
snr = [0:2:10];
ber = zeros(1,length(snr));
for snr_index = 1:length(snr)
    
%% *********************** Transmitter ******************************
    % Generate the information bits
    inf_bits = randn(1,TotalNumBits)>0; % standard normal distribution (u = 0 and sigma = 1)

    %Modulate
    paradata = reshape(inf_bits,length(inf_bits)/ml,ml);
    mod_ofdm_syms = qammod(bi2de(paradata),2^ml)/NormFactor;
    mod_ofdm_syms = reshape(mod_ofdm_syms,fftlen,NumSyms);

    % IFFT
    tx_blks = sqrt(fftlen)*ifft(mod_ofdm_syms); % help(fft) for scaler

    % P/S
    tx = reshape(tx_blks,NumSyms*fftlen,1);
    % SNR = 10lg(S/N) --> N = 1/10^(SNR/10) --> /2 for I/Q
    noise_var = 1/(10^(snr(snr_index)/10))/2;
    len = length(tx);
    noise = sqrt(noise_var) * (randn(len,1) + j*randn(len,1));
    % add noise
    rx_signal = tx + noise;

    rx_signal = reshape(rx_signal, fftlen, NumSyms);
    freq_data = fft(rx_signal)/sqrt(fftlen);     
    Data_seq = reshape(freq_data,fftlen*NumSyms,1)*NormFactor;

    % demodulate
    DemodSeq = de2bi(qamdemod(Data_seq,2^ml),ml);
    SerialBits = reshape(DemodSeq,size(DemodSeq,1)*ml,1).';
    err = sum(abs(SerialBits-inf_bits));
    ber(snr_index) = err/TotalNumBits;
end
semilogy(snr,ber,'-b.');hold on;