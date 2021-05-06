

# 作业四 信道编码与解码

在最基本的OFDM系统中，加入信道编码。在信号传输过程中将出现差错，故对数字信号必须采用纠、检错技术，即纠、检错编码技术，以增强数据在信道中传输时抵御各种干扰的能力，提高系统的可靠性。

## 卷积码的编码

卷积码是一种差错控制编码，由P.Elias于1955年发明。因为数据与二进制多项式滑动相关故称卷积码。卷积码在通信系统中应用广泛，如IS-95，TD-SCDMA，WCDMA，IEEE 802.11及卫星等系统中均使用了卷积码。卷积码是一种有记忆的纠错码，编码规则是将k个信息比特编码形成n个比特，编码后的n个码元不但与当前输入的k个信息有关，仍与之前的L-1组的信息有关。IEEE 802.11中的卷积码编码器实现原理如下：

![image-20210505164816235](https://i.loli.net/2021/05/05/uidvpcSosIZeKgY.png)

输入的比特流经过六个移位寄存器和两个模2加法器组成。编码前每一个字寄存器清零，信息比特流顺序送入编码器，每输入一个比特就会交替从A或B端输出，将A端与B端的结果拼接在一起即为信道编码的结果。在IEEE 802.11中经由循环码编码后比特流长度会变为原来的两倍。在`Matlab`中利用`convenc`函数完成卷积码的编码操作。

## 卷积码的译码

卷积码的译码方式有两种——代数译码和概率译码。这里采用维特比译码：

* 在j=L-1个时刻前，计算每一个状态单个路径分支度量。
* 第x-1个时刻开始，对进入每一个状态的部分路径进行计算，这样的路径有$2^k$条，挑选具有最大部分路径值的部分路径为幸存路径，删去进入该状态的其他路径，然后，幸存路径向前延长一个分支。
* 重复上一步的计算、比较和判决过程。若输入接收序列长为（x+L-1）k，其中，后L-1段是人为加入的全0段，则译码进行至（x+ L-1）时刻为止。

在`Matlab`中采用`vitdec()`函数对循环码进行译码。

## 加信道编码后的系统模拟

实验成果如图：

![](https://i.loli.net/2021/05/05/zC9NFdbKlia1ORL.png)

代码如下：

```matlab
%% ************** Preparation part ********************
clear all; clc;
% system parameters
ml = 2;                      % Modulation level: 2--4QAM; 4--16QAM; 6--64QAM
NormFactor = sqrt(2/3*(ml.^2-1));
fftlen = 64;
bits_per_sym = fftlen*ml;
NumSyms = 50;
TotalNumBits = bits_per_sym*NumSyms;
%% ************** Loop start***************************
snr = [0:1:10];
ber = zeros(1,length(snr));
trellis=poly2trellis(7,[133 171]);
for snr_index = 1:length(snr)
    
%% *********************** Transmitter ******************************
    % Generate the information bits
    inf_bits = randn(1,TotalNumBits)>0; % standard normal distribution (u = 0 and sigma = 1)
    af_cbits = convenc(inf_bits,trellis);
    NumSyms=2*NumSyms;
    %Modulate
    paradata = reshape(af_cbits,length(af_cbits)/ml,ml);
    mod_ofdm_syms = qammod(bi2de(paradata),2^ml)/NormFactor;
    mod_ofdm_syms = reshape(mod_ofdm_syms,fftlen,NumSyms);

    % IFFT
    tx_blks = sqrt(fftlen)*ifft(mod_ofdm_syms); % help(fft) for scaler

    % P/S
    tx = reshape(tx_blks,NumSyms*fftlen,1);
    % SNR = 10lg(S/N) --> N = 1/10^(SNR/10) --> /2 for I/Q
    noise_var = 1/(10^(snr(snr_index)/10))/2;
    len = length(tx);
    noise = sqrt(noise_var) * (randn(len,1) + j*randn(len,1));
    % add noise
    rx_signal = tx + noise;

    rx_signal = reshape(rx_signal, fftlen, NumSyms);
    freq_data = fft(rx_signal)/sqrt(fftlen);     
    Data_seq = reshape(freq_data,fftlen*NumSyms,1)*NormFactor;

    % demodulate
    NumSyms=0.5*NumSyms;
    DemodSeq = de2bi(qamdemod(Data_seq,2^ml),ml);
    DemodSeq = reshape(DemodSeq,1,[]);
    SerialBits = vitdec(DemodSeq,trellis,NumSyms,'trunc','hard');
    
    %SerialBits = reshape(DemodSeq,size(DemodSeq,1)*ml,1).';
    err = sum(abs(SerialBits-inf_bits));
    ber(snr_index) = err/TotalNumBits;
end
semilogy(snr,ber,'-b.');hold on;
```



# 作业五 交织

交织是为了在时域或者频域，或者同时在时域、频域上分布传输的信息比特，使信道的突发错误在时间上得以扩散，从而使得译码器可以将他们当做随机错误处理。交织放在信道编码之后，接收端先解交织将信道引起的突发错误转换为随机错误，然后再信道译码纠正随机错误，因为信道编解码可纠正随机错误，难纠正连续很多的突发错误。

## 802.11中的交织的实现

* 所有编码数据比特应由块交织器交织，其块大小对应于一个OFDM符号。
* **交织器由两步置换定义**。第一步确保相邻的编码比特被映射到不相邻的子载波上。第二步确保相邻的编码比特交替映射到星座图中较低和较高的有效比特上，从而避免低可靠性比特的长时间运行。如图所示：



![image-20210505171919843](https://i.loli.net/2021/05/05/6sjurSoPQULpVGH.png)

在`Matlab`中可以采用矩阵行列的变换达到交织的实现。给定矩阵：
$$
a=\begin{pmatrix}
50  &43 &36 &29	&22 &15	&8	&1\\
51  &44	&37	&30 &23	&16	&9	&2\\
52	&45	&38	&31	&24 &17 &10	&3\\
53	&46	&39	&32	&25	&18	&11	&4\\
54	&47	&40	&33	&26	&19	&12	&5\\
55	&48	&41	&34	&27	&20	&13	&6\\
56	&49	&42	&35	&28	&21	&14	&7\\
\end{pmatrix}
$$
令其进行列的倒序后再进行行的倒序，得到我们的交织结果：
$$
a’=\begin{pmatrix}
7	&14	&21	&28	&35	&42	&49	&56\\
6	&13	&20	&27	&34	&41	&48	&55\\
5	&12	&19	&26	&33	&40	&47	&54\\
4	&11	&18	&25	&32	&39	&46	&53\\
3	&10	&17	&24	&31	&38	&45	&52\\
2	&9	&16	&23	&30	&37	&44	&51\\
1	&8	&15	&22	&29	&36	&43	&50\\
\end{pmatrix}
$$

## 802.11中的解交织的实现

解交织的方法与交织的方法相同，即对结果进行行倒序与列倒序即可：
$$
a=\begin{pmatrix}
50  &43 &36 &29	&22 &15	&8	&1\\
51  &44	&37	&30 &23	&16	&9	&2\\
52	&45	&38	&31	&24 &17 &10	&3\\
53	&46	&39	&32	&25	&18	&11	&4\\
54	&47	&40	&33	&26	&19	&12	&5\\
55	&48	&41	&34	&27	&20	&13	&6\\
56	&49	&42	&35	&28	&21	&14	&7\\
\end{pmatrix}
$$

## 代码

```matlab
a=fliplr(reshape(1:56,7,8));
af_map=mapping(a);
de_map=demapping(af_map);
function map=mapping(a)
map=flip(fliplr(a));
end
function demap=demapping(map)
demap=flip(fliplr(map));
end
```



# 作业7 前导训练序列

训练序列是在发送比特流的某些固定位置插入一些一直的符号或序列，在接收端利用这些符号或序列按某些算法来进行信道估计。接收机的同步定时、载波频偏估计以及信道估计等都由前置的两个训练符号完成。训练符号由10个周期重复的短序列，和2个周期重复的长训练序两部分组成，通过与接收信号的相关运算，来估计定时、相位与时延，具有很强抗快衰落的能力。

### 前导训练的生成过程

前导训练序列由10个短序列、保护间隔与2个长序列组成，其中保护间隔是防止叠加干扰。前导序列的组成如图所示：

![image-20210505202554783](https://i.loli.net/2021/05/05/GLnMOi6jkADHWKt.png)

### 短序列的生成

* 首先给定能量归一化的短序列调制因子：

$$
S_{-26,26}=\sqrt{\frac{13}{6}}\times{[0 ,0 ,1+j ,0 ,0 ,0 ,-1-j ,0 ,\\0 ,0 ,1+j ,0 ,0 ,0 ,-1-j,\\ 
                       0 ,0 ,0 -1-j ,0 ,0 ,0 ,1+j ,0 ,\\0 ,0 ,0 ,0 ,0 ,-1-j ,0 ,0 ,0, \\  
                       -1-j ,0 ,0 ,0 ,1+j ,0 ,0 ,0 ,1+j \\,0 ,0 ,0 ,1+j ,0 ,0 ,0 ,1+j ,0 ,0]}
$$

前面的$\sqrt{\frac{13}{5}} $是为了将52个子载波中的12个（标号为**-24，-20，-16，-12，-8，-4，4，8，12，16，20，24**）进行能量归一化，使得时域长度为160。

* 其次对提取出来的序列进行映射操作，
* 接着进行64点的逆IFFT操作。

映射操作如图所示：

<img src="https://i.loli.net/2021/05/05/RNu89H7IGQeZL6w.png" alt="image-20210505205801662" style="zoom:50%;" />

要点如下：

* 将子载波子载波标号为（1~26）作为IFFT的1~26作为输入，标号为（-26~-1）作为IFFT的38-63作为输入，其余的为0。
* 原子载波需要经过数据处理：将数据子载波依次标号为1:5 7:19 21:26 27:32 34:46 48:52，故需要在6、20、33、47处插入导频1
* 训练序列对应的子载波则不需要这个插入导频，仅需要经过相应的映射经过IFFT变换即可。



## 长序列的生成

* 首先，同理给出长序列的调制因子$S_L$:
  $$
  S_L=\{
  1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,\\-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1
  \}
  $$

* 进行64点IFFT运算

* 生成保护间隔$G_I $,保护间隔取IFFT运算的后32位

* 生成长训练序列

## 结果展示

短序列实部：

![image-20210506153047144](https://i.loli.net/2021/05/06/hoTyIXDjqYNVdAW.png)



长序列实部：

![image-20210506154031059](https://i.loli.net/2021/05/06/mhnYCUJEylZSAqt.png)

## 代码

```matlab
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
```



