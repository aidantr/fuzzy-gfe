
% implements the FCR estimator from Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)


function [alpha, zeta, wgts, obj, SE] = FCR(t,y,x,Z,G,m,sims,parallel,std_errors)

    % Arguments:
    
    %  t: time vector
    %  y: dependent variabe
    %  x: heterogeneous variables
    %  Z: controls
    %  G: number of groups
    %  m: fuzzy tuning parameter
    %  sims: number of starting values 
    %  parallel: T/F for whether run in parallel
    %  std_errors: T/F for whether to compute standard errors
    
    % Returns:
    
    % alpha: grouped fixed effect coefficients
    % zeta: coefficients on controls
    % weights: group weights
    % SE: standard errors
    % obj: value of objective function at minimum
    
    
    %fmincon options
    options = optimoptions('fmincon','Display','none','SpecifyObjectiveGradient',false,'StepTolerance',1e-10,'ConstraintTolerance',1e-10,'FunctionTolerance',1e-10,'ObjectiveLimit',sqrt(1e-6./(G)),'MaxIterations',1000,'MaxFunctionEvaluations',1000,'OptimalityTolerance',1e-4);
    
    % time dummies
    T = length(grpstats(t,t));
    timed = dummyvar(t);
    
    % constraints, FE ordering
    A = eye(G-1);
    A = [A zeros(G-1,1)]+[zeros(G-1,1) -A];
    A = [A zeros(G-1,G*(T-1)+size(Z,2))] ;
    b = zeros(G-1,1);
    
    % homogeneous specification
    homogeneous = regstats(y,[Z timed],'linear',{'beta','covb','r'}) ; % baseline homogeneous specification
    baseline_guess  = homogeneous.beta;
   
    % pre-allocate memory
    ResFE       = zeros(G*T,sims);
    ResZeta     = zeros(size(Z,2),sims);
    ResL        = zeros(sims,1) + 1e15;
    ResH        = zeros(G*T+size(Z,2),G*T+size(Z,2),sims);
    ResWgts     = zeros(length(y),G,sims);
    
    %loop over starting values
    if parallel
        parfor i=1:sims
            %starting value
            startval = [unifrnd(-1,0,G*T,1); baseline_guess(2:3)];
    
            % minimize objective wrt a
            [ahold,obj,~,~,~,~, hessian] = fmincon(@(a) objective(a,y,Z,timed(:,1:end),G,m),startval,A,b,[],[],[],[],[],options);
    
            %store output
            coeffs            = reshape(ahold,1,G*T+size(Z,2));
            ResFE(:,i)     = coeffs(1:G*T);
            ResZeta(:,i)   = coeffs(G*T+1:end);
            ResWgts(:,:,i) = weights(ResFE(:,i),ResZeta(:,i),y,timed(:,1:end),Z,G,m) ;
            ResL(i)        = obj;
            ResH(:,:,i)    = hessian;
        end

    else
        for i=1:sims
            %starting value
            startval = [unifrnd(-1,0,G*T,1); baseline_guess(2:3)];
    
            % minimize objective w.r.t. a
            [ahold,obj,~,~,~,~, hessian] = fmincon(@(a) objective(a,y,Z,timed(:,1:end),G,m),startval,A,b,[],[],[],[],[],options);
    
            %store output
            coeffs            = reshape(ahold,1,G*T+size(Z,2));
            ResFE(:,i)     = coeffs(1:G*T);
            ResZeta(:,i)   = coeffs(G*T+1:end);
            ResWgts(:,:,i) = weights(ResFE(:,i),ResZeta(:,i),y,timed(:,1:end),Z,G,m) ;
            ResL(i)        = obj;
            ResH(:,:,i)    = hessian;
        end
    
    end
    
    % choose minimizing coefficients
    [~, index]   = sort(ResL, 'ascend');
    
    alpha             = ResFE(:,index(1));
    zeta              = ResZeta(:,index(1));
    obj               = ResL(index(1));
    wgts           = ResWgts(:,:,index(1)) ;
    
    %compute standard errors
    if std_errors
        H = ResH(:,:,index(1));
        GRAD = gradient([alpha' zeta'],y,Z,timed(:,1:end),G,m,A,b);
        V = GRAD*GRAD';
        N = length(y)/size(timed,2);
        VARCOVAR = (inv(H./(N))*(V./(N))*inv(H./(N)))./N;
        SE = sqrt(diag(VARCOVAR));
    end

end

% from the replication package for Lewis, Melcangi, Pilossoph, and Toner-Rodgers (2022)
