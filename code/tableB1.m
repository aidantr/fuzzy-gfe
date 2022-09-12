

% generate simulation results shown in table 1

clear;

addpath(genpath('code/functions'))

%options
m = 1.001;
cores = 250;
startvals = 1000;

parpool(cores)

for sval = [1 5 10 100]

    %loop over choice of G
    for G = [3 5 10]
        
        %loop over simulations
        for rep = 1:100
    
            %set seed
            rng(rep)
    
            %load simulated data
            eval(['load(''data/intermediate/sim_panel_' int2str(G) 'G_1N'')'])
            y = YM(:,:,rep);
            N = length(y)/T;
            controls = XM(:,:,rep);
            t = repmat([1:T]',N,1);
            x = XM(:,1,rep);
            true_group = repelem(ginumb,T,1);
            clearvars -except time_r dem_bias_r inc_bias_r dem_MSE_r inc_MSE_r misclass_r fe zeta weights SEs SE_r G rep theta_par y N controls t x true_group
    
            %run FCR
            timer=tic;
            [fe{rep}, zeta(:,rep), weights{rep}] = FCR(t,y,x,controls,G,m,startvals,true,false);
            time_r(rep) = toc(timer);
    
            %bias, MSE, and misclassification
            dem_bias_r(rep) = abs(theta_par(1)-zeta(1,rep));
            inc_bias_r(rep) = abs(theta_par(2)-zeta(2,rep));

            
            orderings = perms(1:G);
            [~,max_wgts] = max(weights{rep},[],2);
            for ij=1:size(orderings,1)
                max_wgts_test = orderings(ij,max_wgts)';
                misclass_test(ij) = 1-mean(true_group == max_wgts_test);
            end
    
            misclass_r(rep)=min(misclass_test);
    
        end
    
        %store results for table
        dem_bias(G,sval) = mean(dem_bias_r);
        inc_bias(G,sval) = mean(inc_bias_r);
        dem_MSE(G,sval) = sqrt(sum(dem_bias_r.^2));
        inc_MSE(G,sval) = sqrt(sum(inc_bias_r.^2));
        misclass(G,sval) = mean(misclass_r);
        time(G,sval) = mean(time_r(2:end));
   
    end 
end

save('output/tableB1')

% from the replication package for Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)








    
