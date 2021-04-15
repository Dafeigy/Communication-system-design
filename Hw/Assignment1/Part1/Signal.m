function [E1,E2,E3] = Signal(f,v,c,d,L)

t=0.1:0.001:12;

E1=cos(2*pi*f*((1-v/c).*t-d/c))./(d+v.*t); 
E2=cos(2*pi*f*((1+v/c).*t-(2*L-d)/c))./(2*L-d-v.*t);
E3=2*sin(2*pi*f*(v*t/c+(d-L)/c)).*sin(2*pi*f*(t-L/c))./(d+v*t);