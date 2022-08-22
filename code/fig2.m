
% generates figure 2 which plots FE coefficients for G=4

cd /Users/aidan/Dropbox/FCR_panel_m1/replication/

addpath(genpath('code/functions'))

%% estimation
[fe, zeta, weights] = FCR(t,y,x,controls,G,1.001,100,250,false);


%% load BM coefficients and plot


load('data/raw/BM_fe_4G.mat')

group1=fe([1:4:25],:);
group2=fe([2:4:26],:);
group3=fe([3:4:27],:);
group4=fe([4:4:28],:);

group1_BM=alpha([1:7],:);
group2_BM=alpha([8:14],:);
group3_BM=alpha([15:21],:);
group4_BM=alpha([22:28],:);


hold on
x=[1970 1975 1980 1985 1990 1995 2000];
plot(x,group2,'-s','Color','#005c8c', 'MarkerFaceColor','#005c8c','LineWidth',1.3)
plot(x,group2_BM,'--o','Color','#a9d8e6','LineWidth',1.3)
plot(x,group3_BM,'--o','Color','#a9d8e6','LineWidth',1.3)
plot(x,group4_BM,'--o','Color','#a9d8e6','LineWidth',1.3)
plot(x,group1_BM,'--o','Color','#a9d8e6','LineWidth',1.3)
plot(x,group1,'-s','Color','#005c8c','MarkerFaceColor','#005c8c','LineWidth',1.3)
plot(x,group3,'-s','Color','#005c8c','MarkerFaceColor','#005c8c','LineWidth',1.3)
plot(x,group4,'-s','Color','#005c8c','MarkerFaceColor','#005c8c','LineWidth',1.3)
xlabel("",'fontsize',1)
ylabel("Group-Time Effects",'fontsize',18)
yticks([-.6 -.4 -.2 0])
grid on
ax = gca;
ax.FontSize = 13;
legend('FCR','GFE','Location','southeast','fontsize',10)
saveas(gcf,"output/fig2",'pdf')
hold off