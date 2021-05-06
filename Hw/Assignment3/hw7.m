clc;
clear all;
ShortTraining = sqrt(13/6)*[0,0,1+j,0,0,0,-1-j,0,0,0,1+j,0,0,0,-1-j,0,0,0,-1-j,0,0,0,1+j,0,0,0,0,0,0,-1-j,0,0,0,-1-j,0,0,0,1+j,0,0,0,1+j,0,0,0,1+j,0,0,0,1+j,0,0]';
LongTraining=[1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1]';
UsedSubcIdx = [7:32 34:59];
reorder = [33:64 1:32];  
Shortmap = zeros(64, 1);
Shortmap(UsedSubcIdx,:)=ShortTraining;
Shortmap(reorder,:) = Shortmap;
ShortTrain=sqrt(64)*ifft(sqrt(64/52)*Shortmap);
ShortTrain=ShortTrain(1:16);%16*10;%每个短序列长度16
ShortTrain_Final=repmat(ShortTrain,10,1);%一共有10个，进行叠加；
disp(ShortTrain_Final');
Longmap = zeros(64,1);
Longmap(UsedSubcIdx,:)=LongTraining;
Longmap(reorder,:) = Longmap;
LongTrain=sqrt(64)*ifft(sqrt(64/52)*Longmap);
GI=LongTrain(33:64,:);%GI2,时域长度为32，提取LongTrain中的后32长度作为循环前缀
LongTrain_Final=[GI;LongTrain;LongTrain];
disp(LongTrain_Final');
final=[ShortTrain_Final;LongTrain_Final]';%训练序列
disp(final);