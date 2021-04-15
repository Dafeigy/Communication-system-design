h = zeros(1,64);
h(1) = 1;
h(5) = 1;
H = fft(h);

h(2) = 0.5; h(3) = 0.4; h(4) = 0.2;
H = fft(h);
hold on; plot(abs(H),'-r');