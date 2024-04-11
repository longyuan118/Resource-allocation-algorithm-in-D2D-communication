clear,clc
close all

%% 初始化各值
R=500;     
Ndd=50;   
Ncell=20;  
Lmin=25;  
Lmax=80;   
Pmax=24;  
P0=-78;   
PN=116;    
N_UE=9;   
N_eNB=5; 
B=180;     
a=0.8;     
nSum=500;  
D=1;     
k1_BP=13;
Dcell_new=zeros(Ncell,1);
Ddd_new=zeros(Ncell,1);
Ddd2dd_new=zeros(Ncell,1);
Ddd2dd_cell=zeros(Ncell,1);
Ddd_eNB=zeros(Ncell,1);
PLdd=zeros(Ncell,1);
PLdd_eNB=zeros(Ncell,1);
PLcell_new=zeros(Ncell,1);
PLdd_new=zeros(Ncell,1);
PLdd2dd_new=zeros(Ncell,1);
PLdd2dd_cell=zeros(Ncell,1);
G_SINRddnew=zeros(nSum,Ncell);
G_SINRnew=zeros(nSum,Ncell);
G_SINRall=zeros(nSum,2*Ncell); 
G_SSINRcell=zeros(nSum,Ncell);
G_TSINRcell=zeros(nSum,Ncell);
G_Thnew=zeros(nSum,Ncell);
G_Thddnew=zeros(nSum,Ncell);
G_Thall=zeros(nSum,2*Ncell);
% G_tnew=zeros(nSum,Ncell);
% G_tddnew=zeros(nSum,Ncell);

SDdd_new=zeros(Ncell,1);
SDcell_new=zeros(Ncell,1);
SDdd_eNB=zeros(Ncell,1);
SPLdd=zeros(Ncell,1);
SPLdd_new=zeros(Ncell,1);
SPLcell_new=zeros(Ncell,1);
SPLdd_eNB=zeros(Ncell,1);
SG_Thnew=zeros(nSum,Ncell);
TG_Thnew=zeros(nSum,Ncell);
N=zeros(Ncell,k1_BP);
HRAN=randi([0 1],Ncell+Ndd+Ndd,3);
TPcell = zeros(Ncell, 1);
global Idd
global SINRdd
for j=1:nSum
    interest=zeros(Ndd+Ndd+Ncell,1);
    for k=1:Ncell
        interest(k)=0;
    end
    for k=Ncell+1:Ncell+Ndd
        interest(k)=1;
    end
    for k=Ncell+Ndd+1:Ncell+Ndd+Ndd
        interest(k)=2;
    end
    for m=1:Ncell
        xcell=rand(Ncell,1)*R*3/2-R;
        ycell=sqrt(3)*rand(Ncell,1)*R-sqrt(3)*R/2;
        i=1;
        while i<=Ncell
            if xcell(i)<=-R/2 && ycell(i)>=sqrt(3)*(xcell(i)+R)
                xcell(i)=xcell(i)+3/2*R;
                ycell(i)=ycell(i)-sqrt(3)*R/2;
            end
            if xcell(i)<=-R/2 && ycell(i)<=-sqrt(3)*(xcell(i)+R)
                xcell(i)=xcell(i)+3/2*R;
                ycell(i)=ycell(i)+sqrt(3)*R/2;
            end
            if xcell(i)^2+ycell(i)^2<Lmin^2
                xcell(i)=rand*3/2*R-R;
                ycell(i)=sqrt(3)*rand*R-sqrt(3)*R/2;
                i=i-1;
            end
            i=i+1;
        end
        Hcell=[xcell,ycell];
        xddt=rand(Ndd,1)*3/2*R-R;
        yddt=sqrt(3)*rand(Ndd,1)*R-sqrt(3)*R/2;
        i=1;
        while i<=Ndd
            if xddt(i)<-R/2 && yddt(i)>=sqrt(3)*(xddt(i)+R)
                xddt(i)=xddt(i)+3/2*R;
                yddt(i)=yddt(i)-sqrt(3)*R/2;
            end
            if xddt(i)<-R/2 && yddt(i)<-sqrt(3)*(xddt(i)+R)
                xddt(i)=xddt(i)+3/2*R;
                yddt(i)=yddt(i)+sqrt(3)*R/2;
            end
            if xddt(i)^2+yddt(i)^2<Lmin^2
                xddt(i)=rand*3/2*R-R;
                yddt(i)=sqrt(3)*rand*R-sqrt(3)*R/2;
                i=i-1;
            end
            i=i+1;
        end
        Hddt=[xddt,yddt];       
        Ddd_dd=rand(Ndd,1)*Lmax;
        z=rand(Ndd,1)*2*pi;
        xddr=zeros(Ndd,1);
        yddr=zeros(Ndd,1);
        for i=1:Ndd
            xddr(i)=xddt(i)+Ddd_dd(i)*cos(z(i));
            yddr(i)=yddt(i)+Ddd_dd(i)*sin(z(i));
        end
        Hddr=[xddr,yddr];    
        Hall=[Hcell;Hddt;Hddr];
        

        Hall_int=[Hall,HRAN];  
        %         newpoint=[-300 -400];
        newpoint_int=[-300 -400 1 1 0];  
        Mdl1=KDTreeSearcher(Hall_int(1:Ncell+Ndd,:));
        [n]=knnsearch(Mdl1,newpoint_int,'k',k1_BP);
        N(m,:)=n;                       
        %         tabulate(interest(n));
        TABLE=tabulate(interest(n));
        maxDate=max(TABLE(:,2));
        %         TABLE(maxDate==(TABLE(:,2)))
        

        if TABLE(maxDate==(TABLE(:,2)))==0    
            for i=1:k1_BP
                if interest(N(m,i),:)==0      
                    Dcell_new(m)=sqrt((newpoint_int(1)-Hall_int(N(m,i),1))^2+(newpoint_int(2)-Hall_int(N(m,i),2))^2);
                    
                    break
                end
            end
            for i=1:k1_BP
                if interest(N(m,i),:)==1    
                    Ddd_eNB(m)=sqrt(Hall_int(N(m,i),1)^2+Hall_int(N(m,i),2)^2);
                    break
                end
            end
            for i=1:k1_BP
                if interest(N(m,i))==0
                    SDcell_new(m)=sqrt((newpoint_int(1)-Hall_int(N(m,i),1))^2+(newpoint_int(2)-Hall_int(N(m,i),2))^2);        % 计算最近的D2D发射端用户距离
                    break
                end
            end

            if Dcell_new(m)<44.5
                PLcell_new(m)=89; 
            else
                PLcell_new(m)=40.3+40*log10(Dcell_new(m));
            end
            if Ddd_dd(i)<44.2
                PLdd(i)=38.47+20*log10(Ddd_dd(i));
            else
                PLdd(i)=40.3+40*log10(Ddd_dd(i));
            end

            if Ddd_eNB(m)<35
                PLdd_eNB(m)=89.3;
            else
                PLdd_eNB(m)=35.24+35*log10(Ddd_eNB(m));
            end
            

            TPcell=min(Pmax,P0+a*PLcell_new);
            TPdd=min(Pmax,P0+a*PLdd);
            RPdd=TPdd-PLdd-normrnd(0,8,Ncell,1)-N_UE;
            RPcell=TPcell-PLcell_new-normrnd(0,8,Ncell,1)-N_eNB;
            Idd=TPdd-PLdd_eNB;      
            SINRcell=zeros(1,Ncell);
            for i=1:Ncell
                SINRcell(i)=10*log10(10^(RPcell(i)/10)/(10^(Idd(i)/10)+10^(-PN/10)));
                %                 SINRcell(i)=10*log10(10^(RPcell(i)/10)/(10^(-PN/10)));
            end
            Thcell=zeros(1,Ncell);
            for i=1:Ncell
                Thcell(i)=(B*log2(1+10^(SINRcell(i)/10)))/100+10;
            end
            

        elseif TABLE(maxDate==(TABLE(:,2)))==1
            for i=1:5
                if interest(N(m,i),:)==1 && (newpoint_int(3)==Hall_int(N(m,i),3)||newpoint_int(4)==Hall_int(N(m,i),4)||newpoint_int(5)==Hall_int(N(m,i),5) )
                    Ddd_new(m)=sqrt((newpoint_int(1)-Hall_int(N(m,i),1))^2+(newpoint_int(2)-Hall_int(N(m,i),2))^2);        
                    break
                elseif interest(N(m,i),:)==1 && Hall_int(N(m,i),3) && newpoint_int(4)~=Hall_int(N(m,i),4)&&newpoint_int(5)~=Hall_int(N(m,i),5)    
                    Icellpoint=find(interest(N(m,:))==0,1,'first');  
                    Icell_int=Hall_int(Icellpoint,:);
                    Ddd2dd_cell(m)=sqrt((Icell_int(n(i),1)-newpoint_int(1))^2+(Icell_int(n(i),2)-newpoint_int(2))^2);
                    newnewpoint=find(interest(n)==1,1,'last');
                    newnewpoint_int=Hall_int(N(m,newnewpoint),:); 
                    Mdl2=KDTreeSearcher(Hall_int(1:Ncell+Ndd,:));
                    [n1]=knnsearch(Mdl2,newnewpoint_int,'k',k1_BP);
                    TABLE=tabulate(interest(n1));
                    maxDate1=max(TABLE(:,2));
                    if TABLE(maxDate1==(TABLE(:,2)))==1 
                        for i1=1:k1_BP
                            if interest(n1(i1))==1   
                                Ddd2dd_new(m)=sqrt((newnewpoint_int(1)-Hall_int(n1(i1),1))^2+(newnewpoint_int(2)-Hall_int(n1(i1),2))^2);
                                break
                            end
                        end
                        break
                    end
                end
            end
            for i=1:k1_BP
                if interest(N(m,i),:)==1 && (newpoint_int(3)==Hall_int(N(m,i),3)||newpoint_int(4)==Hall_int(N(m,i),4)||newpoint_int(5)==Hall_int(N(m,i),5) )
                    SDdd_new(m)=sqrt((newpoint_int(1)-Hall_int(N(m,i),1))^2+(newpoint_int(2)-Hall_int(N(m,i),2))^2);        % 计算最近的D2D发射端用户距离
                    break
                end
                if interest(N(m,i),:)==1
                    SDdd_eNB(m)=sqrt(Hall_int(N(m,i),1)^2+Hall_int(N(m,i),2)^2);   
                    break
                end
            end
            if Ddd_new(m)<44.5
                PLdd_new(m)=89;  
            else
                PLdd_new(m)=40.3+40*log10(Ddd_new(m));
            end
            if Ddd2dd_new(m)<44.5
                PLdd2dd_new(m)=89;
            else
                PLdd2dd_new(m)=40.3+40*log10(Ddd2dd_new(m));
            end
            if Ddd2dd_cell(m)<44.2
                PLdd2dd_cell(m)=38.47+20*log10(Ddd2dd_cell(m));
            else
                PLdd2dd_cell(m)=40.3+40*log10(Ddd2dd_cell(m));
            end
            if SDcell_new(m)<44.5
                SPLcell_new(m)=89;
            else
                SPLcell_new(m)=40.3+40*log10(SDcell_new(m));
            end
            if SDdd_new(m)<44.5
                SPLdd_new(m)=89;
            else
                SPLdd_new(m)=40.3+40*log10(SDdd_new(m));
            end
            if SDdd_eNB(m)<35
                SPLdd_eNB(m)=89.3;
            else
                SPLdd_eNB(m)=35.24+35*log10(SDdd_eNB(m));
            end
            if Ddd_dd(i)<44.2
                SPLdd(i)=38.47+20*log10(Ddd_dd(i));
            else
                SPLdd(i)=40.3+40*log10(Ddd_dd(i));
            end
            TPdd=min(Pmax,P0+a*PLdd_new);
            RPdd=TPdd-PLdd_new-normrnd(0,8,Ncell,1)-N_eNB;
            TPdd2dd=min(Pmax,P0+a*PLdd2dd_new);
            RPdd2dd=TPdd2dd-PLdd2dd_new-normrnd(0,8,Ncell,1)-N_eNB;
            SINRdd=zeros(1,Ncell);
            SINRdd2dd=zeros(1,Ncell);
            Thdd=zeros(1,Ncell);
            Thdd2dd=zeros(1,Ncell);
            Thdd2dd_new=zeros(1,Ncell);
            Tht=zeros(1,Ncell);
            Icell=TPcell-PLdd2dd_cell;
            if newpoint_int(3)==Hall_int(N(m,i),3)||newpoint_int(4)==Hall_int(N(m,i),4)||newpoint_int(5)==Hall_int(N(m,i),5)    % 若没进行多跳传输
                for i=1:Ncell
                    SINRdd(i)=10*log10(10^(RPdd(i)/10)/(10^(Idd(i)/10)+10^(-PN/10)));
                    %                     SINRdd(i)=10*log10(10^(RPdd(i)/10)/(10^(-PN/10)));
                    Thdd(i)=(B*log2(1+10^(SINRdd(i)/10)))/100+10;
%                     Tht(i)=D./Thdd(i);
                end
            elseif Hall_int(N(m,i),3) && newpoint_int(4)~=Hall_int(N(m,i),4)&&newpoint_int(5)~=Hall_int(N(m,i),5)   % 进行了多跳传输，吞吐量则采用木桶效应计算公式
                for i=1:Ncell
                    %                     SINRdd(i)=10*log10(10^(RPdd(i)/10)/(10^(Idd(i)/10)+10^(-PN/10)));
                    SINRdd2dd(i)=10*log10(10^(RPdd2dd(i)/10)/(10^(Icell(i)/10)+10^(-PN/10)));
                    SINRdd(i)=1/2*min( 10*log10(10^(RPdd(i)/10)/(10^(Idd(i)/10)+10^(-PN/10))), 10*log10(10^(RPdd2dd(i)/10)/(10^(Icell(i)/10)+10^(-PN/10))))  ;
                    Thdd(i)=1/2*min(B*log2(1+10^(SINRdd(i)/10)),B*log2(1+10^(SINRdd2dd(i)/10)));
                end
            end
        end
        STPcell=min(Pmax,P0+a*SPLcell_new);
        STPdd=min(Pmax,P0+a*SPLdd_new);
        RPdd=TPdd-PLdd-normrnd(0,8,Ncell,1)-N_UE;
        SRPcell=STPcell-SPLcell_new-normrnd(0,8,Ncell,1)-N_eNB;
        SIdd=STPdd-SPLdd_eNB; 
        SIcell=STPdd-SPLdd_eNB; 
        SSINRcell=zeros(1,Ncell);
        TSINRcell=zeros(1,Ncell);
        for i=1:Ncell
            SSINRcell(i)=10*log10(10^(SRPcell(i)/10)/(10^(SIdd(i)/10)+10^(-PN/10)));
            TSINRcell(i)=10*log10(10^(SRPcell(i)/10)/(10^(-PN/10)));
        end
        SThcell=zeros(1,Ncell);
        TThcell=zeros(1,Ncell);
        for i=1:Ncell
            SThcell(i)=(B*log2(1+10^(SSINRcell(i)/10)))/100+10;
            TThcell(i)=(B*log2(1+10^(TSINRcell(i)/10)))/100+10;
        end
    end
%     G_SINRnew(j,:)=SINRcell;
%     G_SINRddnew(j,:)=SINRdd;
%     G_SINRall(j,:)=[SINRcell,SINRdd];
%     G_SSINRcell(j,:)=SSINRcell;
%     G_TSINRcell(j,:)=TSINRcell;
     
    G_Thnew(j,:)=Thcell;    % SDSP
    G_Thddnew(j,:)=Thdd;    % A2
    G_Thall(j,:)=[Thcell,Thdd];
    SG_Thnew(j,:)=SThcell;  % Random
    TG_Thnew(j,:)=TThcell;  % A1
end

% for i=1:nSum
%    if  G_Thddnew(i,1)==0
%        G_Thddnew(i,:)=G_Thddnew(i-2,:);
%    end
% end

hold on
plot(Pmax_range, mean(throughput_data1, 2), 'r-*', 'LineWidth', 1, 'MarkerIndices',1:1:length(Pmax_range));
h=legend('KDNN算法','Location','northwest')
% h.FontSize=14;
% set(h,'box','off')
xlabel('Pmax/dB','FontSize',12,'fontweight','bold');
ylabel('小区边缘用户和速率/Mbps','FontSize',12,'fontweight','bold');
set(gca,'LineWidth',1.2)
box on
hold off
figure(3)
hold on
range=8:1:28;

M2=hist(TG_Thnew,range);
cdf2= cumsum(M2)/sum(M2);
plot(range,cdf2,'r-*','LineWidth',1);  
h=legend('KDNN算法','Location','southeast');
h.FontSize=14;
set(h,'box','off')
xlabel('小区边缘用户和速率/Mbps')
ylabel('累计分布函数')
set(gca,'LineWidth',1.20)
xlim([8 28])
xticks(8:4:28)
hold off




