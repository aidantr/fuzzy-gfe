
% Generates data used in simulation for table 1

% the code follows the simulation method of "DGP_final.m" in the BM replication package
% takes as input final_data.mat from BM which contains the GFE and
% common coefficient estimates from the application


clear all;
clc;

cd /Users/aidan/Dropbox/FCR_panel_m1/replication/

load('data/raw/final_data');

rng('default');

X0 = ldem_linc;
Y0 = dem;

XX = X0(:,1:2);
YY = Y0;
clear X Y
X = XX;
Y = YY;
T = 7;
K = 2;
N = size(Y,1)/T;
DGP = 1; % IID


for Nmax = [1:14]
    if Nmax == 1
        Nsim = 1000;
    else
        Nsim = 10;
    end
    for G = [3 5 10]
    
        if G==3
            ginumb= [1           3           3           2           3           1 ...
                3           1           1           1           3           3 ...
                1           2           2           2           2           3 ...
                3           3           3           3           2           1 ...
                2           3           3           3           2           3 ...
                1           2           1           3           1           2 ...
                3           3           2           3           3           3 ...
                3           2           3           2           1           3 ...
                3           2           1           2           1           2 ...
                1           3           1           2           1           3 ...
                3           3           3           1           2           1 ...
                3           2           1           2           2           2 ...
                3           3           2           2           2           3 ...
                3           2           3           1           1           2 ...
                1           3           3           1           2           1]';
        end
        
        if G==4
            ginumb=  [   2           4           4           1           4           3 ...
                   3           2           2           3           4           4 ...
                   3           1           1           1           1           4 ...
                   4           4           4           4           1           2 ...
                   1           2           4           4           1           4 ...
                   3           1           2           4           2           1 ...
                   4           4           1           4           4           4 ...
                   4           1           4           1           2           4 ...
                   4           1           3           3           3           1 ...
                   3           4           3           1           3           4 ...
                   4           2           4           3           2           3 ...
                   2           1           3           1           1           1 ...
                   4           4           1           1           1           2 ...
                   4           1           4           3           3           1 ...
                   2           4           4           3           1           3]';
         end
        
        if G==5
            ginumb=      [   1           2           2           4           2           3 ...
                5           1           1           3           2           2 ...
                5           4           4           4           4           2 ...
                2           2           2           2           4           1 ...
                4           1           2           2           4           2 ...
                3           4           1           5           1           5 ...
                2           2           4           2           2           2 ...
                2           4           2           4           1           2 ...
                2           5           3           5           3           4 ...
                3           5           3           5           5           2 ...
                2           1           2           3           1           1 ...
                1           5           3           4           4           5 ...
                2           2           4           4           4           1 ...
                2           4           5           3           3           4 ...
                1           2           2           3           4           3]';
        end
        
        
        if G==10
            ginumb=[ 4           7           7           6           7          10 ...
                9           3           3          10           7           7 ...
                9           6           6           6           6           7 ...
                7           7           7           7           6           8 ...
                6           8           7           7           6           7 ...
                1           6           4           5           8           2 ...
                7           7           6           7           7           7 ...
                7           2           7           2           3           7 ...
                7           2           3           2          10           6 ...
                10           2          10           1           3           7 ...
                7           8           7          10           8           3 ...
                8           2          10           6           6           2 ...
                5           7           6           6           6           4 ...
                7           6           5           3           2           2 ...
                3           7           7           9           6           2
                ]';
        end
        
        
        gi=zeros(N,G);
        
        for g=1:G
            gi(:,g)=(ginumb==g);
        end
        
        
        Ybar_gt=zeros(N*T,1);
        
        Xbar_gt=zeros(N*T,K);
        
        MYbar_gt=zeros(G*T,1);
        MXbar_gt=zeros(G*T,K);
        
        gisum=sum(gi);
        
        for i=1:N
            if gisum(ginumb(i))>1
                for t=1:T
                    Yt=Y(t:T:N*T);
                    Ybar_gt((i-1)*T+t)=mean(Yt(ginumb==ginumb(i)));
                    Xt=X(t:T:N*T,:);
                    Xbar_gt((i-1)*T+t,:)=mean(Xt(ginumb==ginumb(i),:));
                end
            else
                for t=1:T
                    Yt=Y(t:T:N*T);
                    Ybar_gt((i-1)*T+t)=mean(Yt(ginumb==ginumb(i)));
                    Xt=X(t:T:N*T,:);
                    Xbar_gt((i-1)*T+t,:)=Xt(ginumb==ginumb(i),:);
                end
            end
        end
        
        denom=zeros(K,K);
        numer=zeros(K,1);
        
        for i=1:N
            for t=1:T
                denom=denom+(X((i-1)*T+t,:)'-Xbar_gt((i-1)*T+t,:)')*...
                    (X((i-1)*T+t,:)-Xbar_gt((i-1)*T+t,:));
                numer=numer+(X((i-1)*T+t,:)'-Xbar_gt((i-1)*T+t,:)')*...
                    (Y((i-1)*T+t,:)-Ybar_gt((i-1)*T+t,:));
                MYbar_gt((ginumb(i)-1)*T+t)=Ybar_gt((i-1)*T+t);
                MXbar_gt((ginumb(i)-1)*T+t,:)=Xbar_gt((i-1)*T+t,:);
            end
        end
        
        theta_par=denom\numer;
        
        a=MYbar_gt-MXbar_gt*theta_par;
        gitot = kron(gi,eye(T));
        delta_hat = gitot*a;
        
        % group-specific trends
        
        a_r=zeros(T,G);
        for g=1:G
            a_r(:,g)=a((g-1)*T+1:g*T);
        end
        
        
        obj=0;
        ei = zeros(N*T,1);
        for i=1:N
            for t=1:T
                obj=obj+(Y((i-1)*T+t,:)-Ybar_gt((i-1)*T+t,:)-(X((i-1)*T+t,:)-Xbar_gt((i-1)*T+t,:))*theta_par).^2;
                ei((i-1)*T+t) = Y((i-1)*T+t,:)-Ybar_gt((i-1)*T+t,:)-(X((i-1)*T+t,:)-Xbar_gt((i-1)*T+t,:))*theta_par;
            end
        end
        
        sigma2_naive = obj/(N*T-G*T-N-K);
        
        
        for sim=1:Nsim
            for nn = 1:Nmax
                if (DGP == 1)
                    V=randn(N*T,1);
                    err = sqrt(sigma2_naive).*V/std(V);
                    err = reshape(err,T,N)';
                end
                
                if (DGP == 2)
                    eta_var=randn(N,1);
                    e_i_star=kron(eta_var,ones(T,1)).*ei;
                    err = reshape(e_i_star,T,N)';
                end
                
                dm = zeros(N,T+1);
                ym = zeros(N,T);
                Lag_dem = reshape(X(:,1),T,N)';
                Lag_inc = reshape(X(:,2),T,N)';
                Rdelta = reshape(delta_hat,T,N)';
                for i = 1:N
                    dm(i,1) = Lag_dem(i,1);
                    ym(i,:) = Lag_inc(i,:);
                    for t = 2:T+1
                        dm(i,t) = [dm(i,t-1), ym(i,t-1)]*theta_par + Rdelta(i,t-1)+err(i,t-1);
                    end
                end
                
                
                Rdm = reshape(dm(:,2:T+1)',N*T,1);
                Rdm_1 = reshape(dm(:,1:T)',N*T,1);
                
                
                corr(Rdm,Y)
                
                X0 = [Rdm_1 X(:,2)];
                Y0 = Rdm;
                clear X Y
                X = X0;
                Y = Y0;
                
                XM((nn-1)*N*T+1:nn*N*T,:,sim) = X;
                YM((nn-1)*N*T+1:nn*N*T,1,sim) = Y;
            end
        
        end
    
    eval(['save(''data/intermediate/sim_panel_' int2str(G) 'G_' int2str(Nmax) 'N'')'])
    
    end
end
