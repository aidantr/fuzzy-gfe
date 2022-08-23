
% generate simulation results shown in table 1


addpath(genpath('code/functions'))

clear;

%options
m = 1.001;
cores = 250;
startvals = 1000;

parpool(cores)

%loop over choice of G
for G = [3 5 10]
    
    %loop over simulations
    for rep = 1:1000

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
        clearvars -except dem_bias_r inc_bias_r dem_MSE_r inc_MSE_r misclass_r fe zeta weights SEs SE_r G rep theta_par y N controls t x true_group

        %run FCR
        timer=tic;
        [fe{rep}, zeta(:,rep), weights{rep}, ~, SEs{rep}] = FCR(t,y,x,controls,G,m,startvals,true,true);
        time_r(rep) = toc(timer);

        %bias, MSE, misclassification, SEs, and coverage
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

        std_errors = SEs{rep};
        SE_r(:,rep) = std_errors(end-1:end);

        reject_dem_r(rep) = abs(theta_par(1)-zeta(1,rep))/std_errors(end-1) < 1.96;
        reject_inc_r(rep) = abs(theta_par(2)-zeta(2,rep))/std_errors(end) < 1.96;
    
    end

    %store results for table
    dem_bias(G) = mean(dem_bias_r);
    inc_bias(G) = mean(inc_bias_r);
    dem_MSE(G) = mean(dem_MSE_r);
    inc_MSE(G) = mean(inc_MSE_r);
    misclass(G) = mean(misclass_r);
    dem_coverage = mean(reject_dem_r);
    inc_coverage = mean(reject_inc_r);
    SE(:,G) = median(SE_r,2);
    time(G) = mean(time_r(2:end));
    
end

save('output/table1')


% from the replication package for "A Fuzzy Clustering Approach to Estimating
% Grouped Fixed-Effects" by Lewis, Melangi, Pilossoph, and Toner-Rodgers
% (2022)
    