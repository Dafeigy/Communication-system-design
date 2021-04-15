%------------------原数列--------------------
origin=[1 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%------------------改变位置------------------
location1=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5 0];
location2=[1 0 0 0 0 0 0.5 0 0 0 0 0 0 0 0 0];
%------------------改变幅度------------------
amplitude1=[1 1.4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
amplitude2=[1 0.3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%------------------改变相位------------------
phase1=[1 0.5*exp(1i*pi/4) 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
phase2=[1 0.5*exp(1i*pi/8) 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%------------------改变个数------------------
num1=[1 0.5 0 0.8 0 0 0 0 0 0 0 0 0 0 0 0];
num2=[1 0.5 0 0.8 0 0.3 0 0 0 0 0 0 0 0 0 0];
%------------------2048点FFT-----------------
Origin=abs(fft(origin,2048));

Location1=abs(fft(location1,2048));
Location2=abs(fft(location2,2048));

Amp1=abs(fft(amplitude1,2048));
Amp2=abs(fft(amplitude2,2048));

Phase1=abs(fft(phase1,2048));
Phase2=abs(fft(phase2,2048));

Num1=abs(fft(num1,2048));
Num2=abs(fft(num2,2048));


figure(1);
plot(Origin,'r');
hold on;
plot(Location1,'g');
hold on;
plot(Location2,'b');
legend('x1n','x2n1','x2n2');
title('改变位置');
xlabel('频率');
ylabel('幅度');


figure(2);
plot(Origin,'r');
hold on;
plot(Amp1,'g');
hold on;
plot(Amp2,'b');
legend('origin','Amp1','Amp2');
title('改变幅度');
xlabel('频率');
ylabel('幅度');
saveas(gcf,'改变幅度.jpg')
figure(3);
plot(Origin,'r');
hold on;
plot(Phase1,'g');
hold on;
plot(Phase2,'b');
legend('origin','phase1','phase2');
title('改变相位');
xlabel('频率');
ylabel('幅度');
saveas(gcf,'改变相位.jpg')
figure(4);
plot(Origin,'r');
hold on;
plot(Num1,'g');
hold on;
plot(Num2,'b');
legend('origin','num1','num2');
title('改变个数');
xlabel('频率');
ylabel('幅度');
saveas(gcf,'改变个数.jpg')