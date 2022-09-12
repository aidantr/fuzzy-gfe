
% generate simulation results shown in table 3

clear;

addpath(genpath('code/functions'))

%options
m = 1.001;
cores = 250;
startvals = 1000;

parpool(cores)

%loop over choice of G
for G = [3 10]
    
    %loop over dataset size
    for Nmax = 1:14

        %loop over simulations
        for rep = 1:10
    
            %set seed
            rng(rep)
    
            %load simulated data
            eval(['load(''data/intermediate/sim_panel_' int2str(G) 'G_' int2str(Nmax) 'N'')'])
            y = YM(:,:,rep);
            N = length(y)/T;
            controls = XM(:,:,rep);
            t = repmat([1:T]',N,1);
            x = XM(:,1,rep);
            true_group = repmat(repelem(ginumb,T,1),Nmax,1);
            clearvars -except dem_bias_r inc_bias_r dem_MSE_r inc_MSE_r misclass_r fe zeta weights SEs SE_r G rep Nmax theta_par y N controls t x true_group
    
            %run FCR
            timer=tic;
            [fe{rep}, zeta(:,rep), weights{rep}] = FCR(t,y,x,controls,G,m,startvals,true,false);
            time_r(rep) = toc(timer);
    
            %bias, MSE, and misclassification
            dem_bias_r(rep) = abs(theta_par(1)-zeta(1,rep));
            inc_bias_r(rep) = abs(theta_par(2)-zeta(2,rep));
            
            dem_MSE_r(rep) = sqrt((theta_par(1)-zeta(1,rep))^2) ;
            inc_MSE_r(rep) = sqrt((theta_par(2)-zeta(2,rep))^2) ;
            
            orderings = perms(1:G);
            [~,max_wgts] = max(weights{rep},[],2);
            for ij=1:size(orderings,1)
                max_wgts_test = orderings(ij,max_wgts)';
                misclass_test(ij) = 1-mean(true_group == max_wgts_test);
            end
    
            misclass_r(rep)=min(misclass_test);
        
        end
    
        %store results for table
        dem_bias(G,Nmax) = mean(dem_bias_r);
        inc_bias(G,Nmax) = mean(inc_bias_r);
        dem_MSE(G,Nmax) = mean(dem_MSE_r);
        inc_MSE(G,Nmax) = mean(inc_MSE_r);
        misclass(G,Nmax) = mean(misclass_r);
        time(G,Nmax) = mean(time_r(2:end));
    
    end
end

save('output/fig3')



% from the replication package for Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)














    
