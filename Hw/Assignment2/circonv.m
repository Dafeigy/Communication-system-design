function y = circonv(x1,x2,n)
%����ѭ�����  �Լ�ͨ��������˵ó��Ľ��
% input  : ��Ҫ��ѭ����������������Լ�����
% output : ѭ������Ľ�� 
length_x2 = length(x2);
length_x1 = length(x1);%�����еĳ���
x2n = zeros(1,n); x1n = zeros(1,n);%�������չ����
x2n_f = zeros(n,n); %����x2��չ��ѭ��λ�ƺ�ľ���
if(length_x2<n ||length_x1<n )
    x2n = [x2,zeros(1,n-length_x2)];
    x1n = [x1,zeros(1,n-length_x1)];
    x2n_v = [x2n(1), fliplr(x2n)];
    x2n = x2n_v(1,1:n);
%     x1n_v = [x1n(1), fliplr(x1n)];
%     x1n = x1n_v(1,1:n);
end   %��x1��x2����n����չ
i=1;
while (i<=n)
    if(i==1)
        x2n_f(i,:) = x2n(1,:);
    else   
        x2n_f(i,:) = [x2n(n),x2n(1,1:n-1)];
    end
    x2n = x2n_f(i,:);
    i=i+1;
end  %��������
y = x2n_f*x1n';   %��˼���
 
end
