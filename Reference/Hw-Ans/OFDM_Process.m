clear all; clc;

fftlen = 16;
NumSyms = 8;
total_len = fftlen*NumSyms;
sym_seq = [1:total_len];
mod_ofdm_syms = reshape(sym_seq,fftlen,NumSyms);

% IFFT
tx_blks = sqrt(fftlen)*ifft(mod_ofdm_syms);

% P/S
tx = reshape(tx_blks,NumSyms*fftlen,1);

rx_signal = tx;
rx_signal = reshape(rx_signal, fftlen, NumSyms);
freq_data = fft(rx_signal)/sqrt(fftlen);     
Data_seq = reshape(freq_data,fftlen*NumSyms,1);