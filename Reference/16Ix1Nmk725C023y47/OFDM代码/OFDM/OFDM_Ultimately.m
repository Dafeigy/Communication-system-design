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
SNR=20;

%==================================================
%================信号产生===================================
baseband_out_length=Np_carr * Sig_per_carr * bits_per_symbol;  %所输入的比特数目  500*256*4
% rand( 'twister',0);
baseband_out=round(rand(1,baseband_out_length));%输出待调制的二进制比特流

%==============16QAM调制====================================
complex_carrier_matrix=qam16(baseband_out);%列向量

figure;
plot(complex_carrier_matrix,'*r');%16QAM调制后星座图
title('16QAM调制星座图');
axis([-5,5,-5,5]);
grid on   %显示网格线
axis square

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
%时域波形矩阵，行为每载波所含符号数，列IFFT点数，N个子载波映射在其内，每一行即为一个OFDM符号

%=====================计算峰均比PARR====================================
PAPR=10*log10(Sig_per_carr*Nt_carr*(max(max(time_wave_matrix)))^2/sum(sum(time_wave_matrix.^2)));

%==============降PAPR预操作====================================
b=[exp(pi/4*1i) exp(5*pi/4*1i)];
for i=3:Np_carr
    complex_carrier_matrix1(:,i)=b(1)*complex_carrier_matrix1(:,i);
end
PAPR_temp=zeros(1,4);
for i=1:2
    for j=1:2
        complex_carrier_matrix1_temp=[b(i)*complex_carrier_matrix1(:,1) b(j)*complex_carrier_matrix1(:,2) complex_carrier_matrix1(:,3:Np_carr)];
        carriers=(1:Np_carr)+1;
        conjugate_carriers=Nt_carr-carriers+2;
        IFFT_modulation=zeros(Sig_per_carr,Nt_carr);
        IFFT_modulation(:,carriers )=complex_carrier_matrix1_temp ;
        IFFT_modulation(:,conjugate_carriers )=conj(complex_carrier_matrix1_temp);
        time_wave_matrix=ifft(IFFT_modulation,Nt_carr,2);
        PAPR_temp((i-1)*2+j)=10*log10(Sig_per_carr*Nt_carr*(max(max(time_wave_matrix)))^2/sum(sum(time_wave_matrix.^2)));
    end
end
[~,index]=min(PAPR_temp);
complex_carrier_matrix1_temp=[b(ceil(index/2))*complex_carrier_matrix1(:,1) b(index-2*(ceil(index/2)-1))*complex_carrier_matrix1(:,2) complex_carrier_matrix1(:,3:Np_carr)];
carriers=(1:Np_carr)+1;
conjugate_carriers=Nt_carr-carriers+2;
IFFT_modulation=zeros(Sig_per_carr,Nt_carr);
IFFT_modulation(:,carriers )=complex_carrier_matrix1_temp ;
IFFT_modulation(:,conjugate_carriers )=conj(complex_carrier_matrix1_temp);
time_wave_matrix=ifft(IFFT_modulation,Nt_carr,2);
        
figure;
stem(0:Nt_carr-1, abs(IFFT_modulation(2,1:Nt_carr)),'b*-')%第一个OFDM符号的频谱
grid on
axis ([0 Nt_carr -0.5 4.5]);
ylabel('Magnitude');
xlabel('IFFT Bin');
title('第一个OFDM符号各频率对应的幅度（幅频响应）');

figure;
stem(0:Nt_carr-1, (180/pi)*angle(IFFT_modulation(2,1:Nt_carr)), 'b*-')
hold on
plot(0:Nt_carr-1, (180/pi)*angle(IFFT_modulation(2,1:Nt_carr)), 'go')
axis ([0 Nt_carr -200 +200])
grid on
ylabel('Phase (degrees)')
xlabel('IFFT Bin')
title('第一个OFDM符号各频率对应的相位（相频响应）')

figure;
subplot(211)
plot(0:Nt_carr-1,time_wave_matrix(1,:));%第一个OFDM符号的时域波形
axis([0,Nt_carr-1+CP_length, -0.6, 0.6]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('发送端第一个OFDM符号的时域波形');

%=====================添加循环前缀CP====================================
CP=time_wave_matrix(:,Nt_carr-CP_length+(1:CP_length));
time_wave_matrix_add_CP=[CP,time_wave_matrix];
subplot(212)
plot(0:Nt_carr-1+CP_length,time_wave_matrix_add_CP(1,:));%第一个OFDM符号的时域波形
axis([0,Nt_carr-1+CP_length, -0.6, 0.6]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('发送端第一个OFDM符号的时域波形（加CP）');

%=========================并串转换======================================
time_wave_sequence=reshape(time_wave_matrix_add_CP',(Nt_carr+CP_length)*Sig_per_carr,1)';

%=====================添加训练序列TS====================================
% TS_0=zeros(1,TS_length);
% BPSKTable=[-1-8i,1+8i];
% TS_0=2*sqrt(65)*randi([0,1],1,TS_length);
% TS_0=5*sqrt(2)*randi([0,1],1,TS_length);%干扰序列平均功率为TS序列的五倍
% TS_0=2*randi([0,1],1,TS_length)-1;%干扰序列平均功率与TS序列相当
TS_0=sqrt(2)*randi([0,1],1,TS_length);
BPSKTable=[-1,1];
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

figure;
subplot(311)
plot(0:Nt_carr-1+CP_length,TS_Schmidl);%训练序列波形
axis([0,Nt_carr-1+CP_length, -1.5, 1.5]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('Schmidl训练序列波形');

subplot(312)
plot(0:Nt_carr-1+CP_length,TS_Minn);%训练序列波形
axis([0,Nt_carr-1+CP_length, -1.5, 1.5]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('Minn训练序列波形');

subplot(313)
plot(0:Nt_carr-1+CP_length,TS_Park);%训练序列波形
axis([0,Nt_carr-1+CP_length, -1.5, 1.5]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('Park训练序列波形');

%=====================经过高斯噪声信道====================================
figure;
received_time_wave_sequence_add_TS=awgn(time_wave_sequence_add_TS,SNR,'measured');
received_TS=received_time_wave_sequence_add_TS(1:(7*TS_length));
received_TS_Schmidl=received_TS((TS_length+1):(2*TS_length));
received_TS_Minn=received_TS((3*TS_length+1):(4*TS_length));
received_TS_Park=received_TS((5*TS_length+1):(6*TS_length));
subplot(311)
plot(0:Nt_carr-1+CP_length,received_TS_Schmidl);%接收端Schmidl训练序列波形
axis([0,Nt_carr-1+CP_length, -1.5, 1.5]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('接收端Schmidl训练序列波形');

subplot(312)
plot(0:Nt_carr-1+CP_length,received_TS_Minn);%接收端Minn训练序列波形
axis([0,Nt_carr-1+CP_length, -1.5, 1.5]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('接收端Minn训练序列波形');

subplot(313)
plot(0:Nt_carr-1+CP_length,received_TS_Park);%接收端Park训练序列波形
axis([0,Nt_carr-1+CP_length, -1.5, 1.5]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('接收端Park训练序列波形');

%=====================Schmidl方法进行符号定时====================================
Syn_length=3*TS_length/2;
p_Schmidl=zeros(1,Syn_length);
R_Schmidl=zeros(1,Syn_length);
M_Schmidl=zeros(1,Syn_length);
Syn_Schmidl=received_TS(1:(3*TS_length));
for m=1+CP_length:Syn_length+CP_length
    for k=0:Nt_carr/2-1
        p_Schmidl(m-CP_length)=p_Schmidl(m-CP_length)+conj(Syn_Schmidl(m+k))*Syn_Schmidl(m+k+Nt_carr/2);
        R_Schmidl(m-CP_length)=R_Schmidl(m-CP_length)+abs(Syn_Schmidl(m+k+Nt_carr/2))^2;
    end
    M_Schmidl(m-CP_length)=abs(p_Schmidl(m-CP_length))^2/R_Schmidl(m-CP_length)^2;
end
[maximum_Schmidl,max_index_Schmidl]=max(M_Schmidl);

%=====================Minn方法进行符号定时====================================
p_Minn=zeros(1,Syn_length);
R_Minn=zeros(1,Syn_length);
M_Minn=zeros(1,Syn_length);
Syn_Minn=received_TS((2*TS_length+1):(5*TS_length));
for m=CP_length+1:Syn_length+CP_length
    for l=0:1
        for k=0:Nt_carr/4-1
            p_Minn(m-CP_length)=p_Minn(m-CP_length)+conj(Syn_Minn(m+k+l*Nt_carr/2))*Syn_Minn(m+k+l*Nt_carr/2+Nt_carr/4);
            R_Minn(m-CP_length)=R_Minn(m-CP_length)+abs(Syn_Minn(m+k+l*Nt_carr/2+Nt_carr/4))^2;
        end
    end
    M_Minn(m-CP_length)=abs(p_Minn(m-CP_length))^2/R_Minn(m-CP_length)^2;
end
[maximum_Minn, max_index_Minn]=max(M_Minn);

%=====================Park方法进行符号定时====================================
p_Park=zeros(1,Syn_length);
R_Park=zeros(1,Syn_length);
M_Park=zeros(1,Syn_length);
Syn_Park=received_TS((4*TS_length+1):(7*TS_length));
for m=CP_length+1:3*TS_length/2+CP_length
    for k=0:Nt_carr/2-1
        p_Park(m-CP_length)=p_Park(m-CP_length)+Syn_Park(m+k)*Syn_Park(m-k+Nt_carr-1);
        R_Park(m-CP_length)=R_Park(m-CP_length)+abs(Syn_Park(m+k))^2;
    end
    M_Park(m-CP_length)=abs(p_Park(m-CP_length))^2/R_Park(m-CP_length)^2;
end
[maximum_Park, max_index_Park]=max(M_Park);

figure;
plot(1:Syn_length,M_Schmidl,'r')
hold on
plot(1:Syn_length,M_Minn,'g')
hold on
plot(1:Syn_length,M_Park,'b')
grid on; 
% axis([0,400,0,1.1]); 
title('三种方法得到的定时量度曲线'); 
legend('Schmidl算法','Minn算法','Park算法')
xlabel('Time(sample)'); 
ylabel('Timing Metric');

%===========================串并转换========================================
received_time_wave_matrix_add_CP=reshape(received_time_wave_sequence_add_TS((7*TS_length+1):end)',Nt_carr+CP_length,Sig_per_carr)';%去除训练序列
figure;
subplot(311)
plot(0:Nt_carr-1+CP_length,received_time_wave_matrix_add_CP(1,:));%接收端第一个OFDM符号的时域波形（含CP）
axis([0,Nt_carr-1+CP_length, -0.6, 0.6]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('接收端第一个OFDM符号的时域波形（含CP）');

%=====================去除循环前缀CP====================================
received_time_wave_matrix=received_time_wave_matrix_add_CP(:,CP_length+1:CP_length+Nt_carr);%去CP
subplot(312)
plot(0:Nt_carr-1,real(received_time_wave_matrix(1,:)));%接收端第一个OFDM符号的时域波形（含CP）
axis([0,Nt_carr-1+CP_length, -0.6, 0.6]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('接收端去除CP后的第一个OFDM符号的时域波形');

%=====================信道估计与均衡====================================
% TS_FFT=fft(TS_Schmidl((CP_length+1):TS_length));
% received_TS_FFT=fft(received_TS_Schmidl((CP_length+1):TS_length));
% H=received_TS_FFT./TS_FFT;

H=received_TS_Schmidl((CP_length+1):TS_length)./TS_Schmidl((CP_length+1):TS_length);

% H=fft(received_time_wave_matrix(1,:))./fft(time_wave_matrix(1,:));

received_time_wave_matrix_FFT=fft(received_time_wave_matrix,Nt_carr,2);
received_time_wave_matrix_equilibrium_FFT=received_time_wave_matrix_FFT./H;

received_time_wave_matrix_equilibrium=ifft(received_time_wave_matrix_equilibrium_FFT,Nt_carr,2);
subplot(313)
plot(0:Nt_carr-1,real(received_time_wave_matrix_equilibrium(1,:)));%接收端第一个OFDM符号的时域波形（含CP）
axis([0,Nt_carr-1+CP_length, -0.6, 0.6]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('接收端均衡后第一个OFDM符号的时域波形');

%=================================FFT======================================
received_equilibrium_Hermite=fft(received_time_wave_matrix_equilibrium,Nt_carr,2);

%===========================解埃尔米特映射==================================
received_equilibrium=received_equilibrium_Hermite(:,carriers);

%===========================降低PAPR后处理==================================
received_equilibrium(:,3:Np_carr)=received_equilibrium(:,3:Np_carr)/b(1);
received_equilibrium(:,1)=received_equilibrium(:,1)/b(ceil(index/2));
received_equilibrium(:,2)=received_equilibrium(:,2)/b(index-2*(ceil(index/2)-1));

%=============================并串转换=====================================
received_complex_carrier_matrix1=reshape(received_equilibrium',Np_carr*Sig_per_carr,1)';

figure;
plot(received_complex_carrier_matrix1,'*r');%接收端星座图
title('接收端(16QAM解调制前)星座图');
axis([-5,5,-5,5]);
grid on   %显示网格线
axis square

%===========================16QAM解调==================================
demodu_baseband_out=deqam16(received_complex_carrier_matrix1);
[~,ber]=symerr(demodu_baseband_out,baseband_out);

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