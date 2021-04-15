clc
clear all;
close all;
Nt_carr=256;  %���ز���=FFT����---256
Np_carr=127;
Sig_per_carr=100; %ÿ���ز���������---100
ml=4;
bits_per_symbol=4*Nt_carr;      
NormFactor = sqrt(2/3*(ml.^2-1));
y_data=[];
CP=[];%��ʼ��ǰ׺CP
CP_length=10;%CP
snr=[5:20];
carrier=[1:Nt_carr/2-1];
raw_input_length=Sig_per_carr * bits_per_symbol;  %������ı�����Ŀ
for SNR=1:length(snr)
%SNR=15;
%�����,��ʵ����֤�����Խ�������µ����Ӳ�����Խ�٣�20��ʱ���Ѿ�û������
%==============================================================
%======================�źŲ���===================================


raw_input=round(rand(1,raw_input_length));%��������ƵĶ����Ʊ�����
%===================16QAM����====================================
paradata = reshape(raw_input,length(raw_input)/ml,ml);
complex_carrier_matrix=qammod(bi2de(paradata),2^ml)/NormFactor;%������
% figure;
% plot(complex_carrier_matrix,'*r');%16QAM���ƺ�����ͼ
% title('16QAM��������ͼ');
% axis([-5,5,-5,5]);
% grid on   %��ʾ������
% axis square
% %saveas(gcf,'16QAM-100.jpg')
%=============================����ת��================================

complex_carrier_matrix1=reshape(complex_carrier_matrix,Nt_carr,Sig_per_carr);%�����任Sig_per_carr*Nt_carr ����
%==============================================================
%========================IFFT===================================

time_wave_matrix=sqrt(Nt_carr)*ifft(complex_carrier_matrix1);%OFDM���� ��IFFT�б任
%ʱ���ξ�����Ϊÿ�ز���������������IFFT������N�����ز�ӳ�������ڣ�ÿһ�м�Ϊһ��OFDM����

%=======================================================================
%===================================��ѭ��ǰ׺=============================
%=======================================================================


CP=time_wave_matrix(:,Sig_per_carr-CP_length+(1:CP_length));%ȡ�����������еĵ�x�е�x+CP_length��Ԫ����Ϊѭ��ǰ׺
time_wave_matrix_with_CP=[CP time_wave_matrix];
%figure;
% plot(0:Nt_carr-1+CP_length,time_wave_matrix_add_CP(1,:));%��һ��OFDM���ŵ�ʱ����
% axis([0,Nt_carr-1+CP_length, -0.6, 0.6]);
% grid on;
% ylabel('Amplitude');
% xlabel('Time');
% title('���Ͷ˵�һ��OFDM���ŵ�ʱ���Σ���CP��');
%====================�����ŵ�������===========================
h=zeros(1,CP_length)+0.1;
h(1)=1;
cp_snr=reshape(time_wave_matrix_with_CP,Nt_carr*size(time_wave_matrix_with_CP,2),1);
noise_var = 1/(10^(snr(SNR)/10))/2;
len = length(cp_snr);
noise = sqrt(noise_var) * (randn(len,1)*1 +j*randn(len,1));
rx_signal = cp_snr + noise;
rx_signal = reshape(rx_signal,Nt_carr, size(time_wave_matrix_with_CP,2));
%===================ȥ��ѭ��ǰ׺==================
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
    
    %��H��k�����ó�X(k)
    h_k=fft(h,Nt_carr);
    H=diag(h_k);
    H_inv=inv(H);  % ���������˻ָ�X��k��
    X=zeros(Nt_carr,Sig_per_carr);
    X=H_inv*freq_data;
   % figure;
% plot(received_complex_carrier_matrix1,'*r');%���ն�����ͼ
% title('���ն�16QAM�������ǰ����ͼ');
% axis([-5,5,-5,5]);
% grid on   %��ʾ������
% axis square
 %====================���================
Data_seq = reshape(X,Nt_carr*Sig_per_carr,1)*NormFactor;
DemodSeq = de2bi(qamdemod(Data_seq,2^4),4);
SerialBits = reshape(DemodSeq,size(DemodSeq,1)*4,1).';
err = sum(abs(SerialBits-raw_input));
ber(SNR) = err/raw_input_length;
% %==================================================
% %==================================================
% %====================������˹���԰������ŵ�===========================
% received_time_wave_sequence=awgn(time_wave_matrix',snr(SNR),'measured');
% %=======================FFT====================================
% receive_sequence=fft(received_time_wave_sequence,Nt_carr,2);
% %=====================����ת��=============================
% 
% received_complex_carrier_matrix1=reshape(receive_sequence',Nt_carr*Sig_per_carr,1)';
% %======================16QAM���================
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
xlabel('�����SNR');
ylabel('�������ber');
title('������ȵ��������');
hold on;
plot(snr,ber,'-r.');
saveas(gcf,'SNR.jpg')