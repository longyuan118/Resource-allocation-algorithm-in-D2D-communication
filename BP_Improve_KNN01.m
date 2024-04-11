function Err=BP_Improve_KNN01(filename)
format long
% 输入层h、隐含层i、输出层j等一些初始化
h=3;    % 一个特征变量数据，取决于k值选取个数
i=2;
j=2;    % 两种类别

Alpha=0.9; %改进的BP算法，设置学习率动量参数校正率实现神经网络的加速收敛
Beta=0.5; % 校正率
Gamma=0.85; % 动量参数
Tor=0.0001;                 %要求达到的精度
Maxepoch=20;              %定义迭代次数
Accuracy=0;                 %定义准确率
Ntrain=90;                 %训练集个数
Ntest=10;                   %测试集个数
% 随机给w、b赋值 [-1, +1]
V=2*(rand(1,i)-2.5);   %给w1赋值
W=2*(rand(i,j)-0.5);   %w2
Pi=2*(rand(1,i)-0.5);   %b1
Tau=2*(rand(1,j)-0.5);   %b2
DeltaWOld(i,j)=0;       %dw2
DeltaVOld(h,i)=0;        %dw1
DeltaPiOld(i)=0;         %db1
DeltaTauOld(j)=0;         %db2
ck=zeros(1,2);
% 学习过程
Epoch=1;
Error=10; % 初始误差赋值定义变量

%% 训练模型：加载knn分类结果的数据集,使用switch-case语句控制k值的选取结果数据集。第一列特征属性和第三列特征数据（占比）
switch filename
    case 11     % k=11
        data11 = sprintf('data%d.txt',filename);
        load(data11)
        Odesired=data11(:,1);    % 取特征属性
        datanew=data11(:,3);     % 特征数据集
        maxv=max(max(datanew));         %取数据集中最大值
        minv=min(min(datanew));         %取数据集中最小值
        datanorm=2*((datanew-minv)/(maxv-minv)+0.5); %将数据及进行[-1,1]归一化
        while Error>Tor
            Err(Epoch)=0;
            for k=1:Ntrain
                a=datanorm(k,:);
                if data11(k,1)==0
                    ck=[1 0];
                elseif data11(k,1)==1
                    ck=[0 1];
                end
                % 计算隐藏层、输出层及节点激活
                for ki=1:i  %i是和定值，为隐藏层节点数
                    b(ki)=logsig(a*V(:,ki)+Pi(ki));         % 3x2,1x2
                end
                for kj=1:j %j=2
                    c(kj)=logsig(b*W(:,kj)+Tau(kj));
                end
                % 计算输出层、隐藏层误差error、权值w
                d=c.*(1-c).*(ck-c);
                e=b.*(1-b).*(d*W');
                % 调整隐藏层和输出层间的权值DeltaW
                for ki=1:i
                    for kj=1:j
                        DeltaW(ki,kj)=Alpha*b(ki)*d(kj)+Gamma*DeltaWOld(ki,kj);     % 引入Alpha,神经网络的加速收敛
                    end
                end
                W=W+DeltaW;
                DeltaWOld=DeltaW;
                % 调整输入层和隐藏层间的权值DeltaV
                
                for kh=1:1  % 3
                    for ki=1:i % 2
                        DeltaV(kh,ki)=Beta*a(kh)*e(ki);
                    end
                end
                V=V+DeltaV;
                DeltaVold=DeltaV;
                % 调整Pi和Tau的门限
                DeltaPi=Beta*e+Gamma*DeltaPiOld;              %db1更新
                Pi=Pi+DeltaPi;                   %b1
                DeltaPiold=DeltaPi;              %保存更新值
                DeltaTau=Alpha*d+Gamma*DeltaTauOld;
                Tau=Tau+DeltaTau;
                DeltaTauold=DeltaTau;
                % 误差是d1d2d3中的最大值???
                Err(Epoch)=Err(Epoch)+0.5*(d(1)*d(1)+d(2)*d(2));
                %         Err(Epoch)=Err(Epoch)+0.5*(d(1)*d(1)+d(2)*d(2));
            end
            Err(Epoch)=Err(Epoch)/Ntrain;
            Error=Err(Epoch);
            % 迭代次数过多停止训练,当先练次数达到了设置的训练上限停止训练
            if Epoch > Maxepoch
                break;
            end
            Epoch = Epoch +1;   % 如果没超出最大的训练次数迭代次数加一
        end
    case 13     % k=13
        data13 = sprintf('data%d.txt',filename);
        load(data13)
        Odesired=data13(:,1);    % 取特征属性
        datanew=data13(:,3);     % 特征数据集
        maxv=max(max(datanew));         %取数据集中最大值
        minv=min(min(datanew));         %取数据集中最小值
        datanorm=2*((datanew-minv)/(maxv-minv)+0.5); %将数据及进行[-1,1]归一化
        while Error>Tor
            Err(Epoch)=0;
            for k=1:Ntrain
                a=datanorm(k,:);
                if data13(k,1)==0
                    ck=[1 0];
                elseif data13(k,1)==1
                    ck=[0 1];
                end
                % 计算隐藏层、输出层及节点激活
                for ki=1:i  %i是和定值，为隐藏层节点数
                    b(ki)=logsig(a*V(:,ki)+Pi(ki));         % 3x2,1x2
                end
                for kj=1:j %j=2
                    c(kj)=logsig(b*W(:,kj)+Tau(kj));
                end
                % 计算输出层、隐藏层误差error、权值w
                d=c.*(1-c).*(ck-c);
                e=b.*(1-b).*(d*W');
                % 调整隐藏层和输出层间的权值DeltaW
                for ki=1:i
                    for kj=1:j
                        DeltaW(ki,kj)=Alpha*b(ki)*d(kj)+DeltaWOld(ki,kj);     % 引入Alpha,神经网络的加速收敛
                    end
                end
                W=W+DeltaW;
                DeltaWOld=DeltaW;
                % 调整输入层和隐藏层间的权值DeltaV
                
                for kh=1:1  % 3
                    for ki=1:i % 2
                        DeltaV(kh,ki)=Beta*a(kh)*e(ki);
                    end
                end
                V=V+DeltaV;
                DeltaVold=DeltaV;
                % 调整Pi和Tau的门限
                DeltaPi=Beta*e+Gamma*DeltaPiOld;              %db1更新
                Pi=Pi+DeltaPi;                   %b1
                DeltaPiold=DeltaPi;              %保存更新值
                DeltaTau=Alpha*d+Gamma*DeltaTauOld;
                Tau=Tau+DeltaTau;
                DeltaTauold=DeltaTau;
                % 误差是d1d2中的最大值???
                Err(Epoch)=Err(Epoch)+0.5*(d(1)*d(1)+d(2)*d(2));
            end
            Err(Epoch)=Err(Epoch)/Ntrain;
            Error=Err(Epoch);
            % 迭代次数过多停止训练,当先练次数达到了设置的训练上限停止训练
            if Epoch > Maxepoch
                break;
            end
            Epoch = Epoch +1;   % 如果没超出最大的训练次数迭代次数加一
        end
    case 15     % k=15
        data15 = sprintf('data%d.txt',filename);
        load(data15)
        Odesired=data15(:,1);    % 取特征属性
        datanew=data15(:,3);     % 特征数据集
        maxv=max(max(datanew));         %取数据集中最大值
        minv=min(min(datanew));         %取数据集中最小值
        datanorm=2*((datanew-minv)/(maxv-minv)+0.5); %将数据及进行归一化
        while Error>Tor
            Err(Epoch)=0;
            for k=1:Ntrain
                a=datanorm(k,:);
                if data15(k,1)==0
                    ck=[1 0];
                elseif data15(k,1)==1
                    ck=[0 1];
                end
                % 计算隐藏层、输出层及节点激活
                for ki=1:i  %i是和定值，为隐藏层节点数
                    b(ki)=logsig(a*V(:,ki)+Pi(ki));         % 3x2,1x2
                end
                for kj=1:j %j=2
                    c(kj)=logsig(b*W(:,kj)+Tau(kj));
                end
                % 计算输出层、隐藏层误差error、权值w
                d=c.*(1-c).*(ck-c);
                e=b.*(1-b).*(d*W');
                % 调整隐藏层和输出层间的权值DeltaW
                for ki=1:i
                    for kj=1:j
                        DeltaW(ki,kj)=Alpha*b(ki)*d(kj)+DeltaWOld(ki,kj);     % 引入Alpha,神经网络的加速收敛
                    end
                end
                W=W+DeltaW;
                DeltaWOld=DeltaW;
                % 调整输入层和隐藏层间的权值DeltaV
                
                for kh=1:1  % 3
                    for ki=1:i % 2
                        DeltaV(kh,ki)=Beta*a(kh)*e(ki);
                    end
                end
                V=V+DeltaV;
                DeltaVold=DeltaV;
                % 调整Pi和Tau的门限
                DeltaPi=Beta*e+Gamma*DeltaPiOld;              %db1更新
                Pi=Pi+DeltaPi;                   %b1
                DeltaPiold=DeltaPi;              %保存更新值
                DeltaTau=Alpha*d+Gamma*DeltaTauOld;
                Tau=Tau+DeltaTau;
                DeltaTauold=DeltaTau;
                % 误差是d1d2d3中的最大值???
                Err(Epoch)=Err(Epoch)+0.5*(d(1)*d(1)+d(2)*d(2));
                %         Err(Epoch)=Err(Epoch)+0.5*(d(1)*d(1)+d(2)*d(2));
            end
            Err(Epoch)=Err(Epoch)/Ntrain;
            Error=Err(Epoch);
            % 迭代次数过多停止训练,当先练次数达到了设置的训练上限停止训练
            if Epoch > Maxepoch
                break;
            end
            Epoch = Epoch +1;   % 如果没超出最大的训练次数迭代次数加一
        end
    case 17     % k=17
        data17 = sprintf('data%d.txt',filename);
        load(data17)
        Odesired=data17(:,1);    % 取特征属性
        datanew=data17(:,3);     % 特征数据集
        maxv=max(max(datanew));         %取数据集中最大值
        minv=min(min(datanew));         %取数据集中最小值
        datanorm=2*((datanew-minv)/(maxv-minv)); %将数据及进行归一化
        while Error>Tor
            Err(Epoch)=0;
            for k=1:Ntrain
                a=datanorm(k,:);
                if data17(k,1)==0
                    ck=[1 0];
                elseif data17(k,1)==1
                    ck=[0 1];
                end
                % 计算隐藏层、输出层及节点激活
                for ki=1:i  %i是和定值，为隐藏层节点数
                    b(ki)=logsig(a*V(:,ki)+Pi(ki));         % 3x2,1x2
                end
                for kj=1:j %j=2
                    c(kj)=logsig(b*W(:,kj)+Tau(kj));
                end
                % 计算输出层、隐藏层误差error、权值w
                d=c.*(1-c).*(ck-c);
                e=b.*(1-b).*(d*W');
                % 调整隐藏层和输出层间的权值DeltaW
                for ki=1:i
                    for kj=1:j
                        DeltaW(ki,kj)=Alpha*b(ki)*d(kj)+DeltaWOld(ki,kj);     % 引入Alpha,神经网络的加速收敛
                    end
                end
                W=W+DeltaW;
                DeltaWOld=DeltaW;
                % 调整输入层和隐藏层间的权值DeltaV
                
                for kh=1:1  % 3
                    for ki=1:i % 2
                        DeltaV(kh,ki)=Beta*a(kh)*e(ki);
                    end
                end
                V=V+DeltaV;
                DeltaVold=DeltaV;
                % 调整Pi和Tau的门限
                DeltaPi=Beta*e+Gamma*DeltaPiOld;              %db1更新
                Pi=Pi+DeltaPi;                   %b1
                DeltaPiold=DeltaPi;              %保存更新值
                DeltaTau=Alpha*d+Gamma*DeltaTauOld;
                Tau=Tau+DeltaTau;
                DeltaTauold=DeltaTau;
                % 误差是d1d2d3中的最大值???
                Err(Epoch)=Err(Epoch)+0.5*(d(1)*d(1)+d(2)*d(2));
                %         Err(Epoch)=Err(Epoch)+0.5*(d(1)*d(1)+d(2)*d(2));
            end
            Err(Epoch)=Err(Epoch)/Ntrain;
            Error=Err(Epoch);
            % 迭代次数过多停止训练,当先练次数达到了设置的训练上限停止训练
            if Epoch > Maxepoch
                break;
            end
            Epoch = Epoch +1;   % 如果没超出最大的训练次数迭代次数加一
        end
    case 19     % k=19
        data19 = sprintf('data%d.txt',filename);
        load(data19)
        Odesired=data19(:,1);    % 取特征属性
        datanew=data19(:,3);     % 特征数据集
        maxv=max(max(datanew));         %取数据集中最大值
        minv=min(min(datanew));         %取数据集中最小值
        datanorm=2*((datanew-minv)/(maxv-minv)-0.5); %将数据及进行归一化
        while Error>Tor
            Err(Epoch)=0;
            for k=1:Ntrain
                a=datanorm(k,:);
                if data19(k,1)==0
                    ck=[1 0];
                elseif data19(k,1)==1
                    ck=[0 1];
                end
                % 计算隐藏层、输出层及节点激活
                for ki=1:i  %i是和定值，为隐藏层节点数
                    b(ki)=logsig(a*V(:,ki)+Pi(ki));         % 3x2,1x2
                end
                for kj=1:j %j=2
                    c(kj)=logsig(b*W(:,kj)+Tau(kj));
                end
                % 计算输出层、隐藏层误差error、权值w
                d=c.*(1-c).*(ck-c);
                e=b.*(1-b).*(d*W');
                % 调整隐藏层和输出层间的权值DeltaW
                for ki=1:i
                    for kj=1:j
                        DeltaW(ki,kj)=Alpha*b(ki)*d(kj)+Gamma*DeltaWOld(ki,kj);     % 引入Alpha,神经网络的加速收敛
                    end
                end
                W=W+DeltaW;
                DeltaWOld=DeltaW;
                % 调整输入层和隐藏层间的权值DeltaV
                
                for kh=1:1  % 3
                    for ki=1:i % 2
                        DeltaV(kh,ki)=Beta*a(kh)*e(ki);
                    end
                end
                V=V+DeltaV;
                DeltaVold=DeltaV;
                % 调整Pi和Tau的门限
                DeltaPi=Beta*e+Gamma*DeltaPiOld;              %db1更新
                Pi=Pi+DeltaPi;                   %b1
                DeltaPiold=DeltaPi;              %保存更新值
                DeltaTau=Alpha*d+Gamma*DeltaTauOld;
                Tau=Tau+DeltaTau;
                DeltaTauold=DeltaTau;
                % 误差是d1d2d3中的最大值???
                Err(Epoch)=Err(Epoch)+0.5*(d(1)*d(1)+d(2)*d(2));
                %         Err(Epoch)=Err(Epoch)+0.5*(d(1)*d(1)+d(2)*d(2));
            end
            Err(Epoch)=Err(Epoch)/Ntrain;
            Error=Err(Epoch);
            % 迭代次数过多停止训练,当先练次数达到了设置的训练上限停止训练
            if Epoch > Maxepoch
                break;
            end
            Epoch = Epoch +1;   % 如果没超出最大的训练次数迭代次数加一
        end
end

%% 测试数据
for k=1:Ntest % 测试集
    a=datanorm(Ntrain+k,:);
    % 计算隐藏层节点b
    for ki=1:i
        b(ki)=logsig(a*V(:,ki)+Pi(ki));
    end
    % 计算输出层节点c
    for kj=1:j
        c(kj)=logsig(b*W(:,kj)+Tau(kj));
    end
    % 将输出转换为字段格式
    if (c(1)> 0.9)
        Otest(k)=0;
    elseif (c(2)> 0.9)
        Otest(k)=1;
    else
        Otest(k)=2;
    end
    % 计算测试集准确度
    if Otest(k)==Odesired(Ntrain+k)
        Accuracy=Accuracy+1;
    end
end
% Err=Err*100;

% plot(Err);
% plot the NN output and desired output during test

% N=1:Ntest;
% figure; plot(N,Otest,'b-',N,Odesired(91:100),'r-');
% % display the accuracy
% Accuracy = 100*Accuracy/Ntest;
% t=['正确率: ' num2str(Accuracy) '%' ];
% disp(t);
end