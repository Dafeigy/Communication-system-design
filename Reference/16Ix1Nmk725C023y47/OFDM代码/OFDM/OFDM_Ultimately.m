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
SNR=20;

%==================================================
%================�źŲ���===================================
baseband_out_length=Np_carr * Sig_per_carr * bits_per_symbol;  %������ı�����Ŀ  500*256*4
% rand( 'twister',0);
baseband_out=round(rand(1,baseband_out_length));%��������ƵĶ����Ʊ�����

%==============16QAM����====================================
complex_carrier_matrix=qam16(baseband_out);%������

figure;
plot(complex_carrier_matrix,'*r');%16QAM���ƺ�����ͼ
title('16QAM��������ͼ');
axis([-5,5,-5,5]);
grid on   %��ʾ������
axis square

%==============�����任====================================
%����ת��ʱ������ȡ����
complex_carrier_matrix1=reshape(complex_carrier_matrix',Np_carr,Sig_per_carr)';%�����任Sig_per_carr*Nt_carr ����
%complex_carrier_matrix1=conj(reshape(complex_carrier_matrix',Nt_carr,Sig_per_carr)')%����ת��

%==============��������ӳ��====================================
carriers=(1:Np_carr)+1;%����Գ����ز�ӳ��  �������ݶ�Ӧ��IFFT������
conjugate_carriers=Nt_carr-carriers+2;%����Գ����ز�ӳ��  �������Ӧ��IFFT������
IFFT_modulation=zeros(Sig_per_carr,Nt_carr);%��0���IFFT_bin_length IFFT ����
IFFT_modulation(:,carriers )=complex_carrier_matrix1 ;%δ��ӵ�Ƶ�ź� �����ز�ӳ���ڴ˴�
IFFT_modulation(:,conjugate_carriers )=conj(complex_carrier_matrix1);%�����ӳ��

%=================IFFT===========================
time_wave_matrix=ifft(IFFT_modulation,Nt_carr,2);%OFDM���� ��IFFT�б任
%ʱ���ξ�����Ϊÿ�ز���������������IFFT������N�����ز�ӳ�������ڣ�ÿһ�м�Ϊһ��OFDM����

%=====================��������PARR====================================
PAPR=10*log10(Sig_per_carr*Nt_carr*(max(max(time_wave_matrix)))^2/sum(sum(time_wave_matrix.^2)));

%==============��PAPRԤ����====================================
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
stem(0:Nt_carr-1, abs(IFFT_modulation(2,1:Nt_carr)),'b*-')%��һ��OFDM���ŵ�Ƶ��
grid on
axis ([0 Nt_carr -0.5 4.5]);
ylabel('Magnitude');
xlabel('IFFT Bin');
title('��һ��OFDM���Ÿ�Ƶ�ʶ�Ӧ�ķ��ȣ���Ƶ��Ӧ��');

figure;
stem(0:Nt_carr-1, (180/pi)*angle(IFFT_modulation(2,1:Nt_carr)), 'b*-')
hold on
plot(0:Nt_carr-1, (180/pi)*angle(IFFT_modulation(2,1:Nt_carr)), 'go')
axis ([0 Nt_carr -200 +200])
grid on
ylabel('Phase (degrees)')
xlabel('IFFT Bin')
title('��һ��OFDM���Ÿ�Ƶ�ʶ�Ӧ����λ����Ƶ��Ӧ��')

figure;
subplot(211)
plot(0:Nt_carr-1,time_wave_matrix(1,:));%��һ��OFDM���ŵ�ʱ����
axis([0,Nt_carr-1+CP_length, -0.6, 0.6]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('���Ͷ˵�һ��OFDM���ŵ�ʱ����');

%=====================���ѭ��ǰ׺CP====================================
CP=time_wave_matrix(:,Nt_carr-CP_length+(1:CP_length));
time_wave_matrix_add_CP=[CP,time_wave_matrix];
subplot(212)
plot(0:Nt_carr-1+CP_length,time_wave_matrix_add_CP(1,:));%��һ��OFDM���ŵ�ʱ����
axis([0,Nt_carr-1+CP_length, -0.6, 0.6]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('���Ͷ˵�һ��OFDM���ŵ�ʱ���Σ���CP��');

%=========================����ת��======================================
time_wave_sequence=reshape(time_wave_matrix_add_CP',(Nt_carr+CP_length)*Sig_per_carr,1)';

%=====================���ѵ������TS====================================
% TS_0=zeros(1,TS_length);
% BPSKTable=[-1-8i,1+8i];
% TS_0=2*sqrt(65)*randi([0,1],1,TS_length);
% TS_0=5*sqrt(2)*randi([0,1],1,TS_length);%��������ƽ������ΪTS���е��屶
% TS_0=2*randi([0,1],1,TS_length)-1;%��������ƽ��������TS�����൱
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
plot(0:Nt_carr-1+CP_length,TS_Schmidl);%ѵ�����в���
axis([0,Nt_carr-1+CP_length, -1.5, 1.5]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('Schmidlѵ�����в���');

subplot(312)
plot(0:Nt_carr-1+CP_length,TS_Minn);%ѵ�����в���
axis([0,Nt_carr-1+CP_length, -1.5, 1.5]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('Minnѵ�����в���');

subplot(313)
plot(0:Nt_carr-1+CP_length,TS_Park);%ѵ�����в���
axis([0,Nt_carr-1+CP_length, -1.5, 1.5]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('Parkѵ�����в���');

%=====================������˹�����ŵ�====================================
figure;
received_time_wave_sequence_add_TS=awgn(time_wave_sequence_add_TS,SNR,'measured');
received_TS=received_time_wave_sequence_add_TS(1:(7*TS_length));
received_TS_Schmidl=received_TS((TS_length+1):(2*TS_length));
received_TS_Minn=received_TS((3*TS_length+1):(4*TS_length));
received_TS_Park=received_TS((5*TS_length+1):(6*TS_length));
subplot(311)
plot(0:Nt_carr-1+CP_length,received_TS_Schmidl);%���ն�Schmidlѵ�����в���
axis([0,Nt_carr-1+CP_length, -1.5, 1.5]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('���ն�Schmidlѵ�����в���');

subplot(312)
plot(0:Nt_carr-1+CP_length,received_TS_Minn);%���ն�Minnѵ�����в���
axis([0,Nt_carr-1+CP_length, -1.5, 1.5]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('���ն�Minnѵ�����в���');

subplot(313)
plot(0:Nt_carr-1+CP_length,received_TS_Park);%���ն�Parkѵ�����в���
axis([0,Nt_carr-1+CP_length, -1.5, 1.5]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('���ն�Parkѵ�����в���');

%=====================Schmidl�������з��Ŷ�ʱ====================================
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

%=====================Minn�������з��Ŷ�ʱ====================================
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

%=====================Park�������з��Ŷ�ʱ====================================
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
title('���ַ����õ��Ķ�ʱ��������'); 
legend('Schmidl�㷨','Minn�㷨','Park�㷨')
xlabel('Time(sample)'); 
ylabel('Timing Metric');

%===========================����ת��========================================
received_time_wave_matrix_add_CP=reshape(received_time_wave_sequence_add_TS((7*TS_length+1):end)',Nt_carr+CP_length,Sig_per_carr)';%ȥ��ѵ������
figure;
subplot(311)
plot(0:Nt_carr-1+CP_length,received_time_wave_matrix_add_CP(1,:));%���ն˵�һ��OFDM���ŵ�ʱ���Σ���CP��
axis([0,Nt_carr-1+CP_length, -0.6, 0.6]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('���ն˵�һ��OFDM���ŵ�ʱ���Σ���CP��');

%=====================ȥ��ѭ��ǰ׺CP====================================
received_time_wave_matrix=received_time_wave_matrix_add_CP(:,CP_length+1:CP_length+Nt_carr);%ȥCP
subplot(312)
plot(0:Nt_carr-1,real(received_time_wave_matrix(1,:)));%���ն˵�һ��OFDM���ŵ�ʱ���Σ���CP��
axis([0,Nt_carr-1+CP_length, -0.6, 0.6]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('���ն�ȥ��CP��ĵ�һ��OFDM���ŵ�ʱ����');

%=====================�ŵ����������====================================
% TS_FFT=fft(TS_Schmidl((CP_length+1):TS_length));
% received_TS_FFT=fft(received_TS_Schmidl((CP_length+1):TS_length));
% H=received_TS_FFT./TS_FFT;

H=received_TS_Schmidl((CP_length+1):TS_length)./TS_Schmidl((CP_length+1):TS_length);

% H=fft(received_time_wave_matrix(1,:))./fft(time_wave_matrix(1,:));

received_time_wave_matrix_FFT=fft(received_time_wave_matrix,Nt_carr,2);
received_time_wave_matrix_equilibrium_FFT=received_time_wave_matrix_FFT./H;

received_time_wave_matrix_equilibrium=ifft(received_time_wave_matrix_equilibrium_FFT,Nt_carr,2);
subplot(313)
plot(0:Nt_carr-1,real(received_time_wave_matrix_equilibrium(1,:)));%���ն˵�һ��OFDM���ŵ�ʱ���Σ���CP��
axis([0,Nt_carr-1+CP_length, -0.6, 0.6]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('���ն˾�����һ��OFDM���ŵ�ʱ����');

%=================================FFT======================================
received_equilibrium_Hermite=fft(received_time_wave_matrix_equilibrium,Nt_carr,2);

%===========================�Ⱓ������ӳ��==================================
received_equilibrium=received_equilibrium_Hermite(:,carriers);

%===========================����PAPR����==================================
received_equilibrium(:,3:Np_carr)=received_equilibrium(:,3:Np_carr)/b(1);
received_equilibrium(:,1)=received_equilibrium(:,1)/b(ceil(index/2));
received_equilibrium(:,2)=received_equilibrium(:,2)/b(index-2*(ceil(index/2)-1));

%=============================����ת��=====================================
received_complex_carrier_matrix1=reshape(received_equilibrium',Np_carr*Sig_per_carr,1)';

figure;
plot(received_complex_carrier_matrix1,'*r');%���ն�����ͼ
title('���ն�(16QAM�����ǰ)����ͼ');
axis([-5,5,-5,5]);
grid on   %��ʾ������
axis square

%===========================16QAM���==================================
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
title('�����ز��������')
ylabel('�������');
xlabel('���ز����');
xlim([1 Np_carr]);
grid on;