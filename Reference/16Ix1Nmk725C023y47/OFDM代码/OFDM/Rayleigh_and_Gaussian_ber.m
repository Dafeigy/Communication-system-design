clc
clear all;
close all;
Nt_carr=256;      %子载波数=FFT点数---256
Np_carr=Nt_carr/2-1; %实际子载波数---127
Sig_per_carr=500; %每子载波含符号数---500
bits_per_symbol=4;      %每符号含比特数,16QAM调制---4
CP_Ratio=1/8;
CP_length=CP_Ratio*Nt_carr;%循环前缀长度
TS_length=Nt_carr+CP_length;%训练序列长度
SNR_max=25;

%==================================================
%================信号产生===================================
baseband_out_length=Np_carr * Sig_per_carr * bits_per_symbol;  %所输入的比特数目  500*256*4
rand( 'twister',0);
baseband_out=round(rand(1,baseband_out_length));%输出待调制的二进制比特流

%==============16QAM调制====================================
complex_carrier_matrix=qam16(baseband_out);%列向量

%==============串并变换====================================
%矩阵转置时附加有取共轭
complex_carrier_matrix1=reshape(complex_carrier_matrix',Np_carr,Sig_per_carr)';%串并变换Sig_per_carr*Nt_carr 矩阵
%complex_carrier_matrix1=conj(reshape(complex_carrier_matrix',Nt_carr,Sig_per_carr)')%两次转置

%==============埃尔米特映射====================================
carriers=(1:Np_carr)+1;%共轭对称子载波映射  复数数据对应的IFFT点坐标
conjugate_carriers=Nt_carr-carriers+2;%共轭对称子载波映射  共轭复数对应的IFFT点坐标
IFFT_modulation=zeros(Sig_per_carr,Nt_carr);%添0组成IFFT_bin_length IFFT 运算
IFFT_modulation(:,carriers )=complex_carrier_matrix1 ;%未添加导频信号 ，子载波映射在此处
IFFT_modulation(:,conjugate_carriers )=conj(complex_carrier_matrix1);%共轭复数映射

%=================IFFT===========================
time_wave_matrix=ifft(IFFT_modulation,Nt_carr,2);%OFDM调制 即IFFT行变换

%=====================添加循环前缀CP====================================
CP=time_wave_matrix(:,Nt_carr-CP_length+(1:CP_length));
time_wave_matrix_add_CP=[CP,time_wave_matrix];

%=========================并串转换======================================
time_wave_sequence=reshape(time_wave_matrix_add_CP',(Nt_carr+CP_length)*Sig_per_carr,1)';

%=====================添加训练序列TS====================================
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
    %===================经过瑞利信道+高斯噪声信道================================
    H_Rayleigh=raylrnd(0.1,1,(Sig_per_carr+7)*TS_length);
    received_time_wave_sequence_add_TS=time_wave_sequence_add_TS.*H_Rayleigh;

    received_time_wave_sequence_add_TS=awgn(received_time_wave_sequence_add_TS,SNR,'measured');
    received_TS=received_time_wave_sequence_add_TS(1:(7*TS_length));
    received_TS_Schmidl=received_TS((TS_length+1):(2*TS_length));
    received_TS_Minn=received_TS((3*TS_length+1):(4*TS_length));
    received_TS_Park=received_TS((5*TS_length+1):(6*TS_length));

    %=============================均衡==========================================
    %=====================信道估计与时域均衡====================================
    received_time_wave_sequence_add_TS_equilibrium=received_time_wave_sequence_add_TS./H_Rayleigh;

    %===========================串并转换========================================
    received_time_wave_matrix_add_CP_equilibrium=reshape(received_time_wave_sequence_add_TS_equilibrium((7*TS_length+1):end)',Nt_carr+CP_length,Sig_per_carr)';%去除训练序列
    
    %=====================去除循环前缀CP====================================
    received_time_wave_matrix_equilibrium=received_time_wave_matrix_add_CP_equilibrium(:,CP_length+1:CP_length+Nt_carr);%去CP

    %=================================FFT======================================
    received_equilibrium_Hermite=fft(received_time_wave_matrix_equilibrium,Nt_carr,2);

    %===========================解埃尔米特映射==================================
    received_equilibrium=received_equilibrium_Hermite(:,carriers);

    %=============================并串转换=====================================
    received_complex_carrier_matrix1=reshape(received_equilibrium',Np_carr*Sig_per_carr,1)';

    %===========================16QAM解调==================================
    demodu_baseband_out=deqam16(received_complex_carrier_matrix1);
    [~,ber_equilibrium(SNR+1)]=symerr(demodu_baseband_out,baseband_out);
    
	%===============================未均衡======================================
    %===========================串并转换========================================
    received_time_wave_matrix_add_CP=reshape(received_time_wave_sequence_add_TS((7*TS_length+1):end)',Nt_carr+CP_length,Sig_per_carr)';%去除训练序列
    
    %=====================去除循环前缀CP====================================
    received_time_wave_matrix=received_time_wave_matrix_add_CP(:,CP_length+1:CP_length+Nt_carr);%去CP
    
    %=================================FFT======================================
    received_Hermite=fft(received_time_wave_matrix,Nt_carr,2);

    %===========================解埃尔米特映射==================================
    received=received_Hermite(:,carriers);

    %=============================并串转换=====================================
    received_complex_carrier_matrix1_without_equilibrium=reshape(received',Np_carr*Sig_per_carr,1)';

    %===========================16QAM解调==================================
    demodu_baseband_out_without_equilibrium=deqam16(received_complex_carrier_matrix1_without_equilibrium);
    [~,ber(SNR+1)]=symerr(demodu_baseband_out_without_equilibrium,baseband_out);
end
semilogy(0:SNR_max,ber_equilibrium,'--r*');
hold on;
semilogy(0:SNR_max,ber,'--g*');
legend('迫零均衡','未均衡');
title('信号经过高斯噪声信道+瑞利衰落信道的误码率曲线（有衰落）')
xlabel('SNR/dB')
ylabel('Pe')
grid on;