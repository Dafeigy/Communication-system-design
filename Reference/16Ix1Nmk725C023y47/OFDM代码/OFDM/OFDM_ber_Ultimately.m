clc
clear all;
close all;
Nt_carr=256;      %���ز���=FFT����---256
Np_carr=Nt_carr/2-1; %ʵ�����ز���---127
Sig_per_carr=500; %ÿ���ز���������---500
bits_per_symbol=4;      %ÿ���ź�������,16QAM����---4
CP_Ratio=1/8;
CP_length=CP_Ratio*Nt_carr;%ѭ��ǰ׺����
TS_length=Nt_carr+CP_length;%ѵ�����г���
SNR_max=25;

%==================================================
%================�źŲ���===================================
baseband_out_length=Np_carr * Sig_per_carr * bits_per_symbol;  %������ı�����Ŀ  500*256*4
rand( 'twister',0);
baseband_out=round(rand(1,baseband_out_length));%��������ƵĶ����Ʊ�����

%==============16QAM����====================================
complex_carrier_matrix=qam16(baseband_out);%������

%==============�����任====================================
%����ת��ʱ������ȡ����
complex_carrier_matrix1=reshape(complex_carrier_matrix',Np_carr,Sig_per_carr)';%�����任symbols_per_carrier*carrier_count ����
%complex_carrier_matrix1=conj(reshape(complex_carrier_matrix',carrier_count,symbols_per_carrier)')%����ת��

%==============��������ӳ��====================================
carriers=(1:Np_carr)+1;%����Գ����ز�ӳ��  �������ݶ�Ӧ��IFFT������
conjugate_carriers=Nt_carr-carriers+2;%����Գ����ز�ӳ��  �������Ӧ��IFFT������
IFFT_modulation=zeros(Sig_per_carr,Nt_carr);%��0���IFFT_bin_length IFFT ����
IFFT_modulation(:,carriers )=complex_carrier_matrix1 ;%δ��ӵ�Ƶ�ź� �����ز�ӳ���ڴ˴�
IFFT_modulation(:,conjugate_carriers )=conj(complex_carrier_matrix1);%�����ӳ��

%=================IFFT===========================
time_wave_matrix=ifft(IFFT_modulation,Nt_carr,2);%OFDM���� ��IFFT�б任

%=====================���ѭ��ǰ׺CP====================================
CP=time_wave_matrix(:,Nt_carr-CP_length+(1:CP_length));
time_wave_matrix_add_CP=[CP,time_wave_matrix];

%=========================����ת��======================================
time_wave_sequence=reshape(time_wave_matrix_add_CP',(Nt_carr+CP_length)*Sig_per_carr,1)';

%=====================���ѵ������TS====================================
% BPSKTable=[-1-8i,1+8i];
% TS_0=zeros(1,TS_length);
% TS_0=2*sqrt(65)*randi([0,1],1,TS_length);
BPSKTable=[-1,1];
TS_0=sqrt(2)*randi([0,1],1,TS_length);
% TS_0=2*randi([0,1],1,TS_length)-1;
TS_Schmidl_1=BPSKTable(randi([0,1],1,Nt_carr/2)+1);
TS_Schmidl=[TS_Schmidl_1(Nt_carr/2-CP_length+1:Nt_carr/2) TS_Schmidl_1 TS_Schmidl_1];
TS_Minn_1=BPSKTable(randi([0,1],1,Nt_carr/4)+1);
TS_Minn=[-TS_Minn_1(Nt_carr/4-CP_length+1:Nt_carr/4) TS_Minn_1 TS_Minn_1 -TS_Minn_1 -TS_Minn_1];
TS_Park_1=BPSKTable(randi([0,1],1,Nt_carr/4)+1);
TS_Park_2=zeros(1,Nt_carr/4);
for i=1:Nt_carr/4
    TS_Park_2(i)=conj(TS_Park_1(Nt_carr/4-i+1));
end
TS_Park=[conj(TS_Park_2(Nt_carr/4-CP_length+1:Nt_carr/4)) TS_Park_1 TS_Park_2 conj(TS_Park_1) conj(TS_Park_2)];
TS=[TS_0 TS_Schmidl TS_0 TS_Minn TS_0 TS_Park TS_0];
time_wave_sequence_add_TS=[TS time_wave_sequence];

ber_equilibrium=zeros(1,SNR_max+1);
ber=zeros(1,SNR_max+1);
for SNR=0:SNR_max
    %=====================������˹�����ŵ�====================================
    received_time_wave_sequence_add_TS=awgn(time_wave_sequence_add_TS,SNR,'measured');
    received_TS=received_time_wave_sequence_add_TS(1:(7*TS_length));
    received_TS_Schmidl=received_TS((TS_length+1):(2*TS_length));
    received_TS_Minn=received_TS((3*TS_length+1):(4*TS_length));
    received_TS_Park=received_TS((5*TS_length+1):(6*TS_length));

    %===========================����ת��========================================
    received_time_wave_matrix_add_CP=reshape(received_time_wave_sequence_add_TS((7*TS_length+1):end)',Nt_carr+CP_length,Sig_per_carr)';%ȥ��ѵ������

    %=====================ȥ��ѭ��ǰ׺CP====================================
    received_time_wave_matrix=received_time_wave_matrix_add_CP(:,CP_length+1:CP_length+Nt_carr);%ȥCP

    %=====================�ŵ�������ʱ�����====================================
    % TS_FFT=fft(TS_Schmidl((CP_length+1):TS_length));
    % received_TS_FFT=fft(received_TS_Schmidl((CP_length+1):TS_length));
    % H=received_TS_FFT./TS_FFT;

    H=received_TS_Schmidl((CP_length+1):TS_length)./TS_Schmidl((CP_length+1):TS_length);

    % H=fft(received_time_wave_matrix(1,:))./fft(time_wave_matrix(1,:));

    received_time_wave_matrix_FFT=fft(received_time_wave_matrix,Nt_carr,2);
    received_time_wave_matrix_equilibrium_FFT=received_time_wave_matrix_FFT./H;

    received_time_wave_matrix_equilibrium=ifft(received_time_wave_matrix_equilibrium_FFT,Nt_carr,2);

    %=================================����======================================
    %=================================FFT======================================
    received_equilibrium_Hermite=fft(received_time_wave_matrix_equilibrium,Nt_carr,2);

    %===========================�Ⱓ������ӳ��==================================
    received_equilibrium=received_equilibrium_Hermite(:,carriers);

    %=============================����ת��=====================================
    received_complex_carrier_matrix1=reshape(received_equilibrium',Np_carr*Sig_per_carr,1)';

    %===========================16QAM���==================================
    demodu_baseband_out=deqam16(received_complex_carrier_matrix1);
    ber_equilibrium(SNR+1)=symerr(demodu_baseband_out,baseband_out)/baseband_out_length;
    
	%===============================δ����======================================
    %=================================FFT======================================
    received_Hermite=fft(received_time_wave_matrix,Nt_carr,2);

    %===========================�Ⱓ������ӳ��==================================
    received=received_Hermite(:,carriers);

    %=============================����ת��=====================================
    received_complex_carrier_matrix1_without_equilibrium=reshape(received',Np_carr*Sig_per_carr,1)';

    %===========================16QAM���==================================
    demodu_baseband_out_without_equilibrium=deqam16(received_complex_carrier_matrix1_without_equilibrium);
    ber(SNR+1)=symerr(demodu_baseband_out_without_equilibrium,baseband_out)/baseband_out_length;
end
semilogy(0:SNR_max,ber_equilibrium,'--r*');
hold on;
semilogy(0:SNR_max,ber,'--g*');
legend('�������','δ����');
title('�źž�����˹�����ŵ������������ߣ���˥�䣩')
xlabel('SNR/dB')
ylabel('Pe')
grid on;