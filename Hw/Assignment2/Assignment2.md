# 16QAM调制下的OFDM的循环前缀拓展

### 任务内容：

实现如图框架，并于无循环码状态对比：

![image-20210410210206059](https://i.loli.net/2021/04/10/CcwGo3xK4dNurzL.png)



## 实现方式

对于循环前缀，可以对IFFT后的矩阵进行分块选择，选择矩阵的前L列作为循环前缀，并通过矩阵合并的方式将其放置于信码之前；去除循环前缀可采用循环卷积的方法，循环卷积的实现在`circonv.m`中有定义。

## 运行结果：

![SNR](https://i.loli.net/2021/04/15/5hIBHVjaQzCfEkU.jpg)

## 代码

主函数：

```matlab
clc
clear all;
close all;
Nt_carr=256;  %子载波数=FFT点数---256
Np_carr=127;
Sig_per_carr=100; %每子载波含符号数---100
ml=4;
bits_per_symbol=4*Nt_carr;      
NormFactor = sqrt(2/3*(ml.^2-1));
y_data=[];
CP=[];%初始化前缀CP
CP_length=10;%CP
snr=[5:20];
carrier=[1:Nt_carr/2-1];
raw_input_length=Sig_per_carr * bits_per_symbol;  %所输入的比特数目
for SNR=1:length(snr)
%SNR=15;
%信噪比,经实验验证信噪比越大的情况下单个子波误码越少，20的时候已经没有误码
%==============================================================
%======================信号产生===================================


raw_input=round(rand(1,raw_input_length));%输出待调制的二进制比特流
%===================16QAM调制====================================
paradata = reshape(raw_input,length(raw_input)/ml,ml);
complex_carrier_matrix=qammod(bi2de(paradata),2^ml)/NormFactor;%列向量
% figure;
% plot(complex_carrier_matrix,'*r');%16QAM调制后星座图
% title('16QAM调制星座图');
% axis([-5,5,-5,5]);
% grid on   %显示网格线
% axis square
% %saveas(gcf,'16QAM-100.jpg')
%=============================串并转换================================

complex_carrier_matrix1=reshape(complex_carrier_matrix,Nt_carr,Sig_per_carr);%串并变换Sig_per_carr*Nt_carr 矩阵
%==============================================================
%========================IFFT===================================

time_wave_matrix=sqrt(Nt_carr)*ifft(complex_carrier_matrix1);%OFDM调制 即IFFT行变换
%时域波形矩阵，行为每载波所含符号数，列IFFT点数，N个子载波映射在其内，每一行即为一个OFDM符号

%=======================================================================
%===================================加循环前缀=============================
%=======================================================================


CP=time_wave_matrix(:,Sig_per_carr-CP_length+(1:CP_length));%取矩阵所有行中的第x列到x+CP_length的元素作为循环前缀
time_wave_matrix_with_CP=[CP time_wave_matrix];
%figure;
% plot(0:Nt_carr-1+CP_length,time_wave_matrix_add_CP(1,:));%第一个OFDM符号的时域波形
% axis([0,Nt_carr-1+CP_length, -0.6, 0.6]);
% grid on;
% ylabel('Amplitude');
% xlabel('Time');
% title('发送端第一个OFDM符号的时域波形（加CP）');
%====================经过信道加噪声===========================
h=zeros(1,CP_length)+0.1;
h(1)=1;
cp_snr=reshape(time_wave_matrix_with_CP,Nt_carr*size(time_wave_matrix_with_CP,2),1);
noise_var = 1/(10^(snr(SNR)/10))/2;
len = length(cp_snr);
noise = sqrt(noise_var) * (randn(len,1)*1 +j*randn(len,1));
rx_signal = cp_snr + noise;
rx_signal = reshape(rx_signal,Nt_carr, size(time_wave_matrix_with_CP,2));
%===================去除循环前缀==================
af_channel=zeros(Nt_carr,size(time_wave_matrix_with_CP,2));
  for n=1:Nt_carr
        af_channel(n,:)=circonv(rx_signal(n,:),h,size(time_wave_matrix_with_CP,2))';
  end
  % Remove CP
    R_cp=zeros(Nt_carr,Sig_per_carr);
    for m=1:Sig_per_carr
        R_cp(:,m)=af_channel(:,CP_length+m);
    end
    
    % FFT
    freq_data = fft(R_cp)/sqrt(Nt_carr);     
    
    %除H（k），得出X(k)
    h_k=fft(h,Nt_carr);
    H=diag(h_k);
    H_inv=inv(H);  % 求逆矩阵，左乘恢复X（k）
    X=zeros(Nt_carr,Sig_per_carr);
    X=H_inv*freq_data;
   % figure;
% plot(received_complex_carrier_matrix1,'*r');%接收端星座图
% title('接收端16QAM【解调】前星座图');
% axis([-5,5,-5,5]);
% grid on   %显示网格线
% axis square
 %====================解调================
Data_seq = reshape(X,Nt_carr*Sig_per_carr,1)*NormFactor;
DemodSeq = de2bi(qamdemod(Data_seq,2^4),4);
SerialBits = reshape(DemodSeq,size(DemodSeq,1)*4,1).';
err = sum(abs(SerialBits-raw_input));
ber(SNR) = err/raw_input_length;
% %==================================================
% %==================================================
% %====================经过高斯加性白噪声信道===========================
% received_time_wave_sequence=awgn(time_wave_matrix',snr(SNR),'measured');
% %=======================FFT====================================
% receive_sequence=fft(received_time_wave_sequence,Nt_carr,2);
% %=====================并串转换=============================
received_complex_carrier_matrix1=reshape(receive_sequence',Nt_carr*Sig_per_carr,1)';
% %======================16QAM解调================
% demodu_baseband_out=qamdemod(received_complex_carrier_matrix1,16,'OutputType','bit');
% demodu_baseband_out=reshape(demodu_baseband_out,Sig_per_carr * bits_per_symbol,1);
% [~,bre]=symerr(demodu_baseband_out,raw_input');
% ber_carriers=zeros(1,Np_carr);
% for j=1:Np_carr
%     for i=1:Sig_per_carr
%         for k=1:bits_per_symbol/Nt_carr
%             if demodu_baseband_out((i-1)*Np_carr*bits_per_symbol/Nt_carr+(j-1)*bits_per_symbol/Nt_carr+k)~=raw_input((i-1)*Np_carr*bits_per_symbol/Nt_carr+(j-1)*bits_per_symbol/Nt_carr+k)
%                 ber_carriers(j)=ber_carriers(j)+1;
%             end
%         end
%     end
% end
% y_data(end+1)=sum(ber_carriers)/raw_input_length;
end
% plot(snr,y_data,'-b.');
xlabel('信噪比SNR');
ylabel('误比特率ber');
title('各信噪比的误比特率');
hold on;
plot(snr,ber,'-r.');
saveas(gcf,'SNR.jpg')
```

2.循环卷积

```matlab
function y = circonv(x1,x2,n)
%计算循环卷积  自己通过矩阵相乘得出的结果
% input  : 需要做循环卷积的两个序列以及点数
% output : 循环卷积的结果 
length_x2 = length(x2);
length_x1 = length(x1);%求序列的长度
x2n = zeros(1,n); x1n = zeros(1,n);%定义好拓展序列
x2n_f = zeros(n,n); %定义x2拓展后循环位移后的矩阵
if(length_x2<n ||length_x1<n )
    x2n = [x2,zeros(1,n-length_x2)];
    x1n = [x1,zeros(1,n-length_x1)];
    x2n_v = [x2n(1), fliplr(x2n)];
    x2n = x2n_v(1,1:n);
%     x1n_v = [x1n(1), fliplr(x1n)];
%     x1n = x1n_v(1,1:n);
end   %将x1和x2进行n点拓展
i=1;
while (i<=n)
    if(i==1)
        x2n_f(i,:) = x2n(1,:);
    else   
        x2n_f(i,:) = [x2n(n),x2n(1,1:n-1)];
    end
    x2n = x2n_f(i,:);
    i=i+1;
end  %构建矩阵
y = x2n_f*x1n';   %相乘即可
end

```

