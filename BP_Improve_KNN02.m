function [Err]=BP_Improve_KNN02(Knn_Id)
format long
h=3;
i=7;
j=2;
Alpha=0.9; %一种改进的算法，设置学习率动量参数校正率实现神经网络的加速收敛
Beta=0.5; % 校正率
Gamma=0.85; % 动量参数
Tor=0.0005;                 %要求达到的精度
Maxepoch=200;              %定义迭代次数
Accuracy=0;                 %定义准确率
Ntrain=90;                 %训练集个数
Ntest=10;                   %测试集个数
%随机赋值 [-1, +1]
V=2*(rand(h,i)-2.5);   %给w1赋值
W=2*(rand(i,j)-0.5);   %w2
Pi=2*(rand(1,i)-0.5);   %b1
Tau=2*(rand(1,j)-0.5);   %b2
DeltaWOld(i,j)=0;       %dw2
DeltaVOld(h,i)=0;        %dw1
DeltaPiOld(i)=0;         %db1
DeltaTauOld(j)=0;         %db2
% the learning process
Epoch=1;
Error=10; % 初始误差赋值定义变量
load data_all.txt
Odesired=data_all(:,2);
switch Knn_Id
    case  1115
        datanew=data_all(:,3:5);
    case  1317
        datanew=data_all(:,4:6);
    case  1519
        datanew=data_all(:,5:7);
        
end
maxv=max(max(datanew));         %取数据集中最大值
minv=min(min(datanew));         %取数据集中最小值
datanorm=2*((datanew-minv)/(maxv-minv)-0.1); %将数据及进行归一化
while Error>Tor
    Err(Epoch)=0;
    for k=1:Ntrain %  k是训练集数据的索引
        a=datanorm(k,:); %取第k组数据
        % set the desired output ck[j]
        if data_all(k,2)==0    %按实际结果转化成三维矩阵，2代表第二列，也就是数据中的实际结果
            ck=[1 0];
        elseif data_all(k,2)==1
            ck=[0 1];
        end
        %  %计算隐藏层激活
        for ki=1:i %这里的i是和定值，为隐藏层节点数
            b(ki)=logsig(a*V(:,ki)+Pi(ki));
            %logsig是sigmoid，进行节点激活，隐藏层三个节点因此ki=3
            %b=sigmoid(xw+b)
        end
        % 输出层及节点激活
        for kj=1:j %j=3
            c(kj)=logsig(b*W(:,kj)+Tau(kj));
        end
        %  （求输出层w的误差
        d=c.*(1-c).*(ck-c);
        %  (求输入层w的误差
        e=b.*(1-b).*(d*W');
        
        for ki=1:i  
            for kj=1:j  
                DeltaW(ki,kj)=Alpha*b(ki)*d(kj)+Gamma*DeltaWOld(ki,kj);
            end
        end
        W=W+DeltaW;
        DeltaWOld=DeltaW;
         
        for kh=1:h % 4
            for ki=1:i % 3
                DeltaV(kh,ki)=Beta*a(kh)*e(ki);
            end
        end
        V=V+DeltaV;
        DeltaVold=DeltaV;
        
        DeltaPi=Beta*e+Gamma*DeltaPiOld;              %db1更新
        Pi=Pi+DeltaPi;                   %b1
        DeltaPiold=DeltaPi;              %保存更新值
        DeltaTau=Alpha*d+Gamma*DeltaTauOld;
        Tau=Tau+DeltaTau;
        DeltaTauold=DeltaTau;
        %   误差是d1d2中的最大值
        Err(Epoch)=Err(Epoch)+0.5*(d(1)*d(1)+d(2)*d(2));
    end %for k=1:Ntrain
    Err(Epoch)=Err(Epoch)/Ntrain;
    Error=Err(Epoch);
    
    %迭代次数过多停止训练当先练次数达到了设置的训练上限停止训练
    if Epoch > Maxepoch
        break;
    end
    Epoch = Epoch +1; %  如果没超出最大的训练次数迭代次数加一
end
for k=1:Ntest 
    a=datanorm(Ntrain+k,:);
    %  计算隐藏层节点b
    for ki=1:i
        b(ki)=logsig(a*V(:,ki)+Pi(ki));
    end
    %  计算输出层节点c
    for kj=1:j
        c(kj)=logsig(b*W(:,kj)+Tau(kj));
    end
    if (c(1)> 0.9)
        Otest(k)=0;
    elseif (c(2)> 0.9)
        Otest(k)=1;
%     elseif (c(3)> 0.9)
%         Otest(k)=2;
    else
        Otest(k)=3;
    end
    % calculate the accuracy of test sets
    if Otest(k)==Odesired(Ntrain+k)
        Accuracy=Accuracy+1;
    end
end % k=1:Ntest
% disp(Err)
% plot(Err,'b-')
disp(datanorm)
end