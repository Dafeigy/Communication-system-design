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
