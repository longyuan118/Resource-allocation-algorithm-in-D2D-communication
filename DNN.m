function [Accuracy, Err] = testKNN02(Knn_Id)
    format long;
    h = 3;
    i = 7;
    j = 2;
    Alpha = 0.01; 
    Beta = 0.9;    
    Gamma = 0.9;   
    Tor = 0.0001;  
    Maxepoch = 1000; 
    Ntrain = 90; 
    Ntest = 10;   
    V = 2 * (rand(h, i) - 2.5);
    W = 2 * (rand(i, j) - 0.5);
    Pi = 2 * (rand(1, i) - 0.5);
    Tau = 2 * (rand(1, j) - 0.5);
    DeltaWOld = zeros(i, j);
    DeltaVOld = zeros(h, i);
    DeltaPiOld = zeros(1, i);
    DeltaTauOld = zeros(1, j);
   
    data_all = load('data_all.txt');
    switch Knn_Id
        case 1115
            datanew = data_all(:, 3:5);
        case 1317
            datanew = data_all(:, 4:6);
        case 1519
            datanew = data_all(:, 5:7);
    end
    maxv = max(max(datanew));
    minv = min(min(datanew));
    datanorm = 2 * ((datanew - minv) / (maxv - minv) - 0.1);
    Odesired = data_all(:, 2);
    Epoch = 1;
    Error = 10;  
    Err = [];
    while Error > Tor
        Err(Epoch) = 0;
        for k = 1:Ntrain
            a = datanorm(k, :);
            if Odesired(k) == 0
                ck = [1 0];
            elseif Odesired(k) == 1
                ck = [0 1];
            end
            for ki = 1:i
                b(ki) = logsig(a * V(:, ki) + Pi(ki));
            end
            for kj = 1:j
                c(kj) = logsig(b * W(:, kj) + Tau(kj));
            end
            d = c .* (1 - c) .* (ck - c);
            e = b .* (1 - b) .* (d * W');
            for ki = 1:i
                for kj = 1:j
                    DeltaW(ki, kj) = Alpha * b(ki) * d(kj) + Gamma * DeltaWOld(ki, kj);
                end
            end
            W = W + DeltaW;
            DeltaWOld = DeltaW;
            for kh = 1:h
                for ki = 1:i
                    DeltaV(kh, ki) = Beta * a(kh) * e(ki) + Gamma * DeltaVOld(kh, ki);
                end
            end
            V = V + DeltaV;
            DeltaVOld = DeltaV;
            DeltaTau = Alpha * d + Gamma * DeltaTauOld;
            Tau = Tau + DeltaTau;
            DeltaTauOld = DeltaTau;
            DeltaPi = Beta * e + Gamma * DeltaPiOld;
            Pi = Pi + DeltaPi;
            DeltaPiOld = DeltaPi;
            Err(Epoch) = Err(Epoch) + 0.5 * (d(1) * d(1) + d(2) * d(2));
        end
        Err(Epoch) = Err(Epoch) / Ntrain;
        Error = Err(Epoch);
        if Epoch >= Maxepoch
            break;
        end
        Epoch = Epoch + 1;
    end
    
    % 测试神经网络
    Accuracy = 0;
    for k = 1:Ntest
        a = datanorm(Ntrain + k, :);
        for ki = 1:i
            b(ki) = logsig(a * V(:, ki) + Pi(ki));
        end
        for kj = 1:j
            c(kj) = logsig(b * W(:, kj) + Tau(kj));
        end
        if c(1) > 0.9
            Otest(k) = 0;
        elseif c(2) > 0.9
            Otest(k) = 1;
        else
            Otest(k) = 3;
        end
        if Otest(k) == Odesired(Ntrain + k)
            Accuracy = Accuracy + 1;
        end
    end
    Accuracy = Accuracy / Ntest;
end


