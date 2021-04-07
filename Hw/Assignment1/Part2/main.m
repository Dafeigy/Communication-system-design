clc
clear all;
close all;
Nt_carr=256;      %子载波数=FFT点数---256
Np_carr=Nt_carr/2-1; %实际子载波数---127
Sig_per_carr=100; %每子载波含符号数---100
bits_per_symbol=4;      %每符号含比特数,16QAM调制为log_2^16=4
SNR=15; %信噪比,经实验验证信噪比越大的情况下单个子波误码越少，20的时候已经没有误码
%==============================================================
%======================信号产生===================================

baseband_out_length=Np_carr * Sig_per_carr * bits_per_symbol;  %所输入的比特数目  127*500*4
baseband_out=round(rand(1,baseband_out_length));%输出待调制的二进制比特流
%===================16QAM调制====================================

complex_carrier_matrix=qammod(baseband_out',16,'InputType','bit');%列向量
figure;
plot(complex_carrier_matrix,'*r');%16QAM调制后星座图
title('16QAM调制星座图');
axis([-5,5,-5,5]);
grid on   %显示网格线
axis square
%saveas(gcf,'16QAM-100.jpg')
%=============================串并转换================================

complex_carrier_matrix1=reshape(complex_carrier_matrix',Np_carr,Sig_per_carr)';%串并变换Sig_per_carr*Nt_carr 矩阵
carrier=[1:127];%选定载波
%==============================================================
%========================IFFT===================================

time_wave_matrix=ifft(complex_carrier_matrix1,Nt_carr,2);%OFDM调制 即IFFT行变换
%时域波形矩阵，行为每载波所含符号数，列IFFT点数，N个子载波映射在其内，每一行即为一个OFDM符号
%====================经过高斯加性白噪声信道===========================
received_time_wave_sequence=awgn(time_wave_matrix,SNR,'measured');
%=======================FFT====================================
receive_sequence=fft(received_time_wave_sequence,Nt_carr,2);
receive_sequence=receive_sequence(:,carrier);
%=====================并串转换=============================

received_complex_carrier_matrix1=reshape(receive_sequence',Np_carr*Sig_per_carr,1)';

%=================画图===========================
figure;
plot(received_complex_carrier_matrix1,'*r');%接收端星座图
title('接收端16QAM【解调】前星座图');
axis([-5,5,-5,5]);
grid on   %显示网格线
axis square
%saveas(gcf,'16QAMr-100.jpg')

%======================16QAM解调================
demodu_baseband_out=qamdemod(received_complex_carrier_matrix1,16,'OutputType','bit');
demodu_baseband_out=reshape(demodu_baseband_out,Np_carr * Sig_per_carr * bits_per_symbol,1);
[~,ber]=symerr(demodu_baseband_out,baseband_out');
ber_carriers=zeros(1,Np_carr);
for j=1:Np_carr
    for i=1:Sig_per_carr
        for k=1:bits_per_symbol
            if demodu_baseband_out((i-1)*Np_carr*bits_per_symbol+(j-1)*bits_per_symbol+k)~=baseband_out((i-1)*Np_carr*bits_per_symbol+(j-1)*bits_per_symbol+k)
                ber_carriers(j)=ber_carriers(j)+1;
            end
        end
    end
end
figure;
plot(1:Np_carr,ber_carriers,'--r*');
title('各子载波误码个数')
ylabel('误码个数');
xlabel('子载波编号');
xlim([1 Np_carr]);
grid on;
%saveas(gcf,'BER-100.jpg')


