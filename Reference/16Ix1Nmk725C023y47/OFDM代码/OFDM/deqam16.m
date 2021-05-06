function [demodu_bit_symble]=deqam16(Rx_serial_complex_symbols)
%将得到的串行16QAM数据解调成二进制比特流
    N=length(Rx_serial_complex_symbols);
    complex_symbols=reshape(Rx_serial_complex_symbols,length(Rx_serial_complex_symbols),1);
    d=1;

%         mapping=[-3*d  3*d;         %7
%                -d  3*d;         %6
%                 d  3*d;         %2
%               3*d  3*d;         %3
%              -3*d  d;           %5
%                -d  d;           %4
%                 d  d;           %0
%               3*d  d;           %1
%              -3*d  -d;          %13
%                -d  -d;          %12
%                 d  -d;          %8
%               3*d  -d;          %9
%              -3*d  -3*d;        %15
%                -d  -3*d;        %14
%                 d  -3*d;        %10
%               3*d  -3*d];       %11

    mapping=[   d  d;           %0
              3*d  d;           %1
                d  3*d;         %2
              3*d  3*d;         %3
               -d  d;           %4
             -3*d  d;           %5
               -d  3*d;         %6
             -3*d  3*d;         %7
                d  -d;          %8
              3*d  -d;          %9
                d  -3*d;        %10
              3*d  -3*d;        %11
               -d  -d;          %12
             -3*d  -d;          %13
               -d  -3*d;        %14
             -3*d  -3*d];       %15
         
%             0111 0110 0010 0011        7  6  2  3
%             0101 0100 0000 0001        5  4  0  1
%             1101 1100 1000 1001       13 12  8  9
%             1111 1110 1010 1011       15 14 10 11
         
%     d=2;
%     k=sqrt(2);
%     mapping=[   d  0;           %0
%             k*d/2  k*d/2;       %1
%                 0  d;           %3
%            -k*d/2  k*d/2;       %2
%                -d  0;           %6
%            -k*d/2  -k*d/2;      %7
%                 0  -d;          %5
%             k*d/2  -k*d/2;      %4
%               2*d  0;           %8
%               k*d  k*d;         %9
%                 0  2*d;         %11
%              -k*d  k*d;         %10
%              -2*d  0;           %14
%              -k*d  -k*d;        %15
%                 0  -2*d;        %13
%               k*d  -k*d];       %12

%     mapping=[   d  0;           %0
%             k*d/2  k*d/2;       %1
%            -k*d/2  k*d/2;       %2
%                 0  d;           %3
%             k*d/2  -k*d/2;      %4
%                 0  -d;          %5
%                -d  0;           %6
%            -k*d/2  -k*d/2;      %7
%               2*d  0;           %8
%               k*d  k*d;         %9
%              -k*d  k*d;         %10
%                 0  2*d;         %11
%               k*d  -k*d;        %12
%                 0  -2*d;        %13
%              -2*d  0;           %14
%              -k*d  -k*d];       %15

	complex_mapping=complex(mapping(:,1),mapping(:,2));
    metrics=zeros(1,16);
    decode_symble=zeros(1,N);
	for i=1:N
        for j=1:16
            metrics(j)=abs(complex_symbols(i,1)-complex_mapping(j,1));
        end
        [~,decode_symble(i)]=min(metrics);%将离某星座点最近的值赋给decode_symble(i)
	end
    decode_bit_symble=de2bi((decode_symble-1)','left-msb');
    demodu_bit_symble=reshape(decode_bit_symble',1,N*4);
end