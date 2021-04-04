function [complex_qam_data]=qam16(bitdata)
%modulation of 16QAM,modulate bitdata to 16QAM complex signal
    X1=reshape(bitdata,4,length(bitdata)/4)';
    d=1;%min distance of symble
    source=bi2de(X1,'left-msb')+1;
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
    for i=1:length(bitdata)/4
    	qam_data(i,:)=mapping(source(i),:);%data mapping
    end
    complex_qam_data=complex(qam_data(:,1),qam_data(:,2));
end