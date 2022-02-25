%%% Compiler of the streamlined code: This example solves the MCF problem
%%% for one instance

%% Initialization
% Loads data and create graphs

run initialization.m            %GC and GD are defined here

%% Alternative for the sample
nn = max(max(fd),max(td));      % Get number of nodes
nl = length(fd) ;               % Get number of links, no repeated links.

    % Def 2: Traffic, and origin destination pairs.
    Tmat = zeros(4,4);
    Tmat(1,4) = 11;
    Tmat(2,3) = 10;
    [s,t,T] = setup_traffic(Tmat);
    
    [nsp,csp] = getsp_s(fd,td,dd,s,t);
    paths = [nsp]; pcosts = [csp];
    
%% Jump simulation
capmatrix = full(GC.adjacency).*GC.distances;

distmatrix = full(GD.adjacency).*GD.distances;

FTCD = [fd,td,cdd,dd];

[kopath,capctr,dij,pathcosts,linklist,kpath] = setuppathproblem_s(paths, pcosts, capmatrix, s, t, T, FTCD); 

    gamma = 9999;               
    inf_cap = 999999; % Higher than any individual traffic flow
    %Solve the LP problem in Cplex
    [sol] = solve_MCF_s(pathcosts, dij ,kopath,gamma,inf_cap);
        
%Handle the solution vector, considering feasibility
[link_flows, link_duals, comm_duals,se_flows, D, F] = sol_handle_s(sol,dij,FTCD,nl,nn,pathcosts,kpath,s,t);    

%% Check if there is flow in superedges
inf_counter = 0; %set up escape condition, tb a search depth variable
while sum(se_flows(:,3)) > 0 || inf_counter < 10
    
    inf_counter = inf_counter + 1;
    % Find more paths for those o-d pairs still not satisfied
    
    [nsp_cg,csp_cg] = getsp_s_rc(s(se_flows(:,3)>0) ...
                                 ,t(se_flows(:,3)>0) ...
                                 ,link_duals);
    
    % Add columns (paths) to the problem
    paths = [nsp;nsp_cg]; pcosts = [csp;csp_cg];
    
    % re-build the problem for cplex
    [kopath,capctr,dij,pathcosts,linklist,kpath] = setuppathproblem_s(paths, pcosts, capmatrix, s, t, T, FTCD); 
    
    % Solve
    [sol] = solve_MCF_s(pathcosts, dij ,kopath,gamma,inf_cap);
        
    %Handle the solution vector, considering feasibility
    [link_flows, link_duals, comm_duals,se_flows, D, F] = sol_handle_s(sol,dij,FTCD,nl,nn,pathcosts,kpath,s,t);
    
end



%% Full case
% [FTCD , distmatrix , capmatrix ] = adj_mats_s(fd, td, cdd, dd); 
% 
% nn = max(max(fd),max(td));   % Get number of nodes
% nl = length(fd) ;           % Get number of links, no repeated links.
% 
%     % Def 2: Traffic, and origin destination pairs.
%     Tmat = rand(nn,nn) * 300;
%     Tmat = Tmat - diag(diag(Tmat));
%     [s,t,T] = setup_traffic(Tmat);
%     
%       
%     % This function yields
%     % The matrix of maximum flows MF
%     % The paths that produce the maximum flow
%     % Costs of the paths found 
%     [MF , paths , pcosts] = maxflowpaths_st(GC,GD,s,t); 
% 
%     % Obtain the shortest paths and costs 
%     [nsp,csp] = getsp_s(fd,td,dd,s,t);
%         
%     % Ok, so we have the paths to be used. Join them all for the initial
%     % set
%     
%     paths = [nsp;paths]; pcosts = [csp;pcosts];
%     
%%
%counter of acceptance rate in MCMC
counter = 0;

%Topology of the network data:
    %Obtain the central location of every link
[XC,YC] = centeroflinks(xlocation,ylocation,from,to);

    %Obtain the SP distance in the original G
    DG =  distances(GD);

%Historical Rainfall load
if exist('HH','var') == 0
    
    HH=xlsread(".\data\data_rainfall.xlsx",'B:B');
    fit = gevfit(HH); k = fit(1); sigma = fit(2); mu = fit(3);

else
end

%Historical Traffic: Here its simulated in the absence of data
% if exist('TT','var') == 0
%     HH=xlsread(".\data\data_rainfall.xlsx",'B:B');
%     [meantt, sdtt] = lognfit(TT);
% else
% end
% meantt=log((mmtt^2)/sqrt(stt + mmtt^2)); 
% sdtt = sqrt(log((stt/(mmtt^2))+1));

% What would the traffic for each link look like. Something increasing for the whole matrix? 

mmtt = 0.0001*10; stt = 0.0001*1000;

%% Block with the MCMC Iterations
%initialize results empty vectors
Fvec = zeros(n,1);
res = zeros(n,1);
Rvec = zeros(n,1);
Tvec = zeros(n,1);
Pf = zeros(n,1);
Pfv = zeros(n,1);
DS={}; 
CDS=zeros(nl,1);
rng(42)


%% Sample simulation
        
        R = gevrnd(k,sigma,mu);                   % Rainfall intensity  ~100mm
    
        H = hzsim(XC , YC, xlocation, ylocation, R, 15);
    
        Tparam = mmtt + exprnd(stt);              % lognrnd(meantt,sdtt); % Small numbers

%% Event produced from Simulation

        Tmat = gentraffics(population, DG,Tparam); % Traffic Event k

        capacity_h = cap_hazard(H , capacity ); % Rainfall Event k
 
        [FTCD , distmatrix , capmatrix ] = adj_mats_s(fd, td, [capacity_h;capacity_h], dd);
        
%% Obtention of the Routing Cost Matrix

        [s,t,T] = setup_traffic(Tmat);
        
        [kopath,capctr,dij,pathcosts,linklist,kpath] = setuppathproblem_s(paths, pcosts, capmatrix ,s,t,T); 

                gamma = 9999;               
                inf_cap = 999999; % Higher than any individual traffic flow
                %Solve the LP problem in Cplex
                [sol] = solve_MCF_s(pathcosts, dij ,kopath,gamma,inf_cap);
        
%Handle the solution vector, considering feasibility
[link_flows, link_duals, comm_duals,se_flows, D, F] = sol_handle_s(sol,dij,linklist,nl,nn,pathcosts,kpath,s,t);
        
%% Enter the solution not found block

if sum(se_flows(:,3))>0
    
    while inf_counter < 11
                
                % Create new columns
                [nsp,csp] = getsp_s_rc(fd,td,dd,s,t,link_duals); %% This should be performed in modified cost space, first trial there are none calulated
                
                paths = [paths;nsp]; 
                pcosts = [pcosts;csp];
                
                
                [kopath,capctr,dij,pathcosts,linklist,kpath] = setuppathproblem_s(paths, pcosts, capmatrix, s, t, T);
               
                [sol] = solve_MCF_s(pathcosts, dij ,kopath);
                
                % Handle the solution vector, considering feasibility
                [link_flows, link_duals, comm_duals,se_flows, D, F] = sol_handle_s(sol,dij,linklist,nl,nn,pathcosts,kpath);
                
                if sum(se_flows(:,3))>0
                inf_counter = inf_counter + 1;
                                
                else
                inf_counter = 11;   % escape condition
                
                end
    end
    
end
            

%% Average Veh. Cost Calculation and storage of results for each mcmc trial

        res(trial) = abs(sum(sum(D-DG.*Tmat))); 
        
        Rvec(trial) = R;
        
        Tvec(trial) = Tparam;
        
        DS{trial} = link_duals(:,4);
        
        Fvec(trial) = F;
        
        %%
        set(gcf, 'Position',  posize);
set(gcf,'color','w');
p3 = plot(GD,'XData',xlocation,'YData',ylocation,'MarkerSize',5);
p3.EdgeCData = abs(sum(DDS,2));
colormap(jet);
pbaspect([2 1 1])
set(gca,'FontSize', fsize)
colorbar