clear;clc;close all
% subplot(1,2,1)
% hold on
% grid on
% Knn_Id=1115;
% [Err]=BP_Improve_KNN02(Knn_Id);
% plot(Err,'k-*','MarkerIndices',1:50:1001,'linewidth',1)
% Knn_Id=1317;
% [Err]=BP_Improve_KNN02(Knn_Id);
% plot(Err,'b-^','MarkerIndices',1:50:1001)
% Knn_Id=1519;
% [Err]=BP_Improve_KNN02(Knn_Id);
% plot(Err,'r-o','MarkerIndices',1:50:1001)
% legend('K--11~15','K--13~17','K--15~19')
% xlabel('Number of iterations(Epoch)')
% ylabel('Error value(Err)')
% hold off


% subplot(1,2,2)
hold on
filename=11;
Err=BP_Improve_KNN01(filename);
plot(Err,'k-d','MarkerIndices',1:5:21,'LineWidth',1)
filename=13;
Err=BP_Improve_KNN01(filename);
plot(Err,'b-o','MarkerIndices',1:5:21,'LineWidth',1)
filename=15;
Err=BP_Improve_KNN01(filename);
plot(Err,'r-*','MarkerIndices',1:5:21,'LineWidth',1)
hold off
h=legend('K=11','K=13','K=15')
h.FontSize = 14;
set(h,'box','off')
xlabel('迭代回合数')
ylabel('误差值')
set(gca,'LineWidth',1.2)
xlim([1 21])
xticks(1:5:21)
ylim([0.0085 0.016])
hold off

% filename=17;
% Err=BP_Improve_KNN01(filename);
% plot(Err,'b-')

% filename=19;
% Err=BP_Improve_KNN(filename);
% plot(Err,'b-')


