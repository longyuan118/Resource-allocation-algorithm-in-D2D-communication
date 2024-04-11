clear,clc,close all
hold on

Knn_Id=1115;
[Err]=BP_Improve_KNN02(Knn_Id);
plot(Err,'r-*','MarkerIndices',1:40:201,'LineWidth', 1)
Knn_Id=1317;
[Err]=BP_Improve_KNN02(Knn_Id);
plot(Err,'b-o','MarkerIndices',1:40:201,'LineWidth', 1)
Knn_Id=1519;
[Err]=BP_Improve_KNN02(Knn_Id);
plot(Err,'k-d','MarkerIndices',1:40:201,'LineWidth', 1)
h=legend('K=11~15','K=13~17','K=15~19')
h.FontSize = 14;
set(h,'box','off')
xlabel('迭代回合数','FontSize',12,'fontweight','bold')
ylabel('误差值','FontSize',12,'fontweight','bold')
set(gca,'linewidth',1.2);
xlim([1 201])
xticks(1:40:201)
ylim([0.010 0.016])
hold off
