clc
clear all;
close all;
Nt_carr=256;      %���ز���=FFT����---256
Np_carr=Nt_carr/2-1; %ʵ�����ز���---127
Sig_per_carr=100; %ÿ���ز���������---100
bits_per_symbol=4;      %ÿ���ź�������,16QAM����Ϊlog_2^16=4
SNR=15; %�����,��ʵ����֤�����Խ�������µ����Ӳ�����Խ�٣�20��ʱ���Ѿ�û������
%==============================================================
%======================�źŲ���===================================

baseband_out_length=Np_carr * Sig_per_carr * bits_per_symbol;  %������ı�����Ŀ  127*500*4
baseband_out=round(rand(1,baseband_out_length));%��������ƵĶ����Ʊ�����
%===================16QAM����====================================

complex_carrier_matrix=qammod(baseband_out',16,'InputType','bit');%������
figure;
plot(complex_carrier_matrix,'*r');%16QAM���ƺ�����ͼ
title('16QAM��������ͼ');
axis([-5,5,-5,5]);
grid on   %��ʾ������
axis square
%saveas(gcf,'16QAM-100.jpg')
%=============================����ת��================================

complex_carrier_matrix1=reshape(complex_carrier_matrix',Np_carr,Sig_per_carr)';%�����任Sig_per_carr*Nt_carr ����
carrier=[1:127];%ѡ���ز�
%==============================================================
%========================IFFT===================================

time_wave_matrix=ifft(complex_carrier_matrix1,Nt_carr,2);%OFDM���� ��IFFT�б任
%ʱ���ξ�����Ϊÿ�ز���������������IFFT������N�����ز�ӳ�������ڣ�ÿһ�м�Ϊһ��OFDM����
%====================������˹���԰������ŵ�===========================
received_time_wave_sequence=awgn(time_wave_matrix,SNR,'measured');
%=======================FFT====================================
receive_sequence=fft(received_time_wave_sequence,Nt_carr,2);
receive_sequence=receive_sequence(:,carrier);
%=====================����ת��=============================

received_complex_carrier_matrix1=reshape(receive_sequence',Np_carr*Sig_per_carr,1)';

%=================��ͼ===========================
figure;
plot(received_complex_carrier_matrix1,'*r');%���ն�����ͼ
title('���ն�16QAM�������ǰ����ͼ');
axis([-5,5,-5,5]);
grid on   %��ʾ������
axis square
%saveas(gcf,'16QAMr-100.jpg')

%======================16QAM���================
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
title('�����ز��������')
ylabel('�������');
xlabel('���ز����');
xlim([1 Np_carr]);
grid on;
%saveas(gcf,'BER-100.jpg')


