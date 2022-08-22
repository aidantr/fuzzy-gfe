
% generates figure 1 which plots common coefficients

clear all; clc;

addpath(genpath('/code/functions'))

parpool(250);

%% estimation

%load data
load('data/raw/BM_RHS_panel.mat')
load('data/raw/BM_LHS_panel.mat')

y = BMLHS(:,1);
id = BMLHS(:,2);
controls = BMRHS(:,1:2);
t = BMRHS(:,10);
x = BMRHS(:,1);

% loop over G
for G = 2:15
    [fe{G}, zeta(:,G), weights{G}, ~, SE{G}] = FCR(t,y,x,controls,G,1.001,1000,250,true);
end

save('output/fig1')


%% load BM coefficients and plot figure

%load BM coefficient results
load('data/raw/BM_coeffs.mat')


%vectors of coefficients
lag_dem_coeffs = [baseline_guess(end-1) zeta(1,2:end)];
lag_dem_BM = Vect_par(:,1);
inc_coeffs = [baseline_guess(end) zeta(2,2:end)];
inc_BM = Vect_par(:,2);

for i = 2:15
    fcr_se_dem(i-1) = kmeans_se{:,i}(end-1);
    fcr_se_inc(i-1) = kmeans_se{:,i}(end);
end

FCR_inf_inc = [Vect_inf(1,2) (inc_coeffs(2:end) - 1.96.*fcr_se_inc)];
FCR_inf_dem = [Vect_inf(1,1) (lag_dem_coeffs(2:end) - 1.96.*fcr_se_dem)];
FCR_sup_inc = [Vect_sup(1,2) (inc_coeffs(2:end) + 1.96.*fcr_se_inc)];
FCR_sup_dem = [Vect_sup(1,1) (lag_dem_coeffs(2:end) + 1.96.*fcr_se_dem)];


hold on
x=[1:15];
plot(x,lag_dem_BM,'--o','Color','#a9d8e6','LineWidth',1.5)
plot(x,lag_dem_coeffs,'-s','Color','#005c8c', 'MarkerFaceColor','#005c8c','LineWidth',1.5)
plot(x,Vect_inf_b(:,1),':','Color','#a9d8e6','LineWidth',1.75)
plot(x,Vect_sup_b(:,1),':','Color','#a9d8e6','LineWidth',1.75)
plot(x,FCR_inf_dem,':','Color','#005c8c','LineWidth',1.75)
plot(x,FCR_sup_dem,':','Color','#005c8c','LineWidth',1.75)
axis([1 15 0 1])
xlabel("Groups",'fontsize',18)
ylabel("Coefficient",'fontsize',18)
yticks([0 0.2 0.4 .6 0.8 1])
grid on
ax = gca;
ax.FontSize = 18;
legend('GFE','FCR','Location','northwest','fontsize',10)
saveas(gcf,"output/fig1a",'pdf')
hold off


hold on
x=[1:15];
plot(x,inc_BM,'--o','Color','#a9d8e6','LineWidth',1.5)
plot(x,inc_coeffs,'-s','Color','#005c8c', 'MarkerFaceColor','#005c8c','LineWidth',1.5)
plot(x,Vect_inf_b(:,2),':','Color','#a9d8e6','LineWidth',1.75)
plot(x,Vect_sup_b(:,2),':','Color','#a9d8e6','LineWidth',1.75)
plot(x,FCR_inf_inc,':','Color','#005c8c','LineWidth',1.75)
plot(x,FCR_sup_inc,':','Color','#005c8c','LineWidth',1.75)
axis([1 15 0 .15])
xlabel("Groups",'fontsize',18)
ylabel("Coefficient",'fontsize',18)
yticks([0 0.025 0.05 0.075 0.1 0.125 .15])
grid on
ax = gca;
ax.FontSize = 18;
legend('GFE','FCR','Location','northwest','fontsize',10)
saveas(gcf,"output/fig1b",'pdf')
hold off







