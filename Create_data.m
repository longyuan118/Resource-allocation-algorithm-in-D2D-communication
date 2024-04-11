clear,clc,close all;

R=500;     
Ndd=80;  
Ncell=100;  
Lmin=25;  
Lmax=80;  
nSum=1000;
N=zeros(10,3);

for m=1:2:100
%确定蜂窝用户位置
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

%确定D2D用户位置
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
Hddt=[xddt,yddt];       % d2d发送端

Ddd_dd=rand(Ndd,1)*Lmax;
z=rand(Ndd,1)*2*pi;
xddr=zeros(Ndd,1);
yddr=zeros(Ndd,1);
for i=1:Ndd
    xddr(i)=xddt(i)+Ddd_dd(i)*cos(z(i));
    yddr(i)=yddt(i)+Ddd_dd(i)*sin(z(i));
end
Hddr=[xddr,yddr];       % d2d接收端
Hall=[Hcell;Hddt;Hddr];


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

    newpoint=[-300 -400];
    Mdl1=KDTreeSearcher(Hall(1:Ncell+Ndd,:));
    
%     % k=5
%     [n]=knnsearch(Mdl1,newpoint,'k',5);
%     A=tabulate(interest(n));
%     N5_old(m:m+1,:)=A(:,:); 
    % k=7
%     [n]=knnsearch(Mdl1,newpoint,'k',7);
%     A=tabulate(interest(n));
%     N7_old(m:m+1,:)=A(:,:);
    % k=9
%     [n]=knnsearch(Mdl1,newpoint,'k',9);
%     A=tabulate(interest(n));
%     N9_old(m:m+1,:)=A(:,:);
    % k=11
    [n]=knnsearch(Mdl1,newpoint,'k',11);
    A=tabulate(interest(n));
    N11_old(m:m+1,:)=A(:,:);
    % k=13
    [n]=knnsearch(Mdl1,newpoint,'k',13);
    B=tabulate(interest(n));
    N13_old(m:m+1,:)=B(:,:);
    % k=15
    [n]=knnsearch(Mdl1,newpoint,'k',15);
    C=tabulate(interest(n));
    N15_old(m:m+1,:)=C(:,:);
    % k=17
    [n]=knnsearch(Mdl1,newpoint,'k',17);
    D=tabulate(interest(n));
    N17_old(m:m+1,:)=D(:,:);
     % k=19
    [n]=knnsearch(Mdl1,newpoint,'k',19);
    E=tabulate(interest(n));
    N19_old(m:m+1,:)=E(:,:);
    
end
N11=N11_old(randperm(size(N11_old, 1)),:);
N13=N13_old(randperm(size(N13_old, 1)),:);
N15=N15_old(randperm(size(N15_old, 1)),:);
N17=N17_old(randperm(size(N17_old, 1)),:);
N19=N19_old(randperm(size(N19_old, 1)),:);
dlmwrite('data11.txt',N11,'delimiter','\t','newline','pc');
dlmwrite('data13.txt',N13,'delimiter','\t','newline','pc');
dlmwrite('data15.txt',N15,'delimiter','\t','newline','pc');
dlmwrite('data17.txt',N17,'delimiter','\t','newline','pc');
dlmwrite('data19.txt',N19,'delimiter','\t','newline','pc');

N11(:,4)=N13(:,3);
N11(:,5)=N15(:,3);
N11(:,6)=N17(:,3);
N11(:,7)=N19(:,3);
N11(:,[1 2])=N11(:,[2 1]);
dlmwrite('D:\matlab2022\project\data_all.txt',N11,'delimiter','\t','newline','pc');

