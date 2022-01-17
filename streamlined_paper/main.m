%%% Compiler of the streamlined code: This example solves the MCF problem
%%% for one instance

%% Initialization
% Loads data and create graphs

run initialization.m %GC and GD are defined here

%%
[FTCD , distmatrix , capmatrix ] = adj_mats_s(fd, td, cdd, dd); 

n = max(max(fd),max(td));   % Get number of nodes
nl = length(fd) ;           % Get number of links, no repeated links.
%% Flows matrix section
% Considering there is a certain demand in X (vehicles, tonnes, etc):
% should be same units as the capacity of the links
% Define a traffic matrix with origins s, destinations t and flows T

    % Format:
    % s = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8];  
    % t = [8,7,6,5,4,2,1,8,7,6,4,3,2,1,1];
    % T = [500,100,100,500,20,50,400,200,200,400,500,0,200,300,40]*1;

% Alternatively, the T matrix can be given shaped with origins (rows)
% destinations (columns)
% OD matrices are typical, such as in the case of Argentina (view example)

% Approach 3: Random OD matrix is produced:

rng(123)                            % Set seed
Tmat = rand(n,n) * 100;             % Scale according to capacities
Tmat = Tmat - diag(diag(Tmat));     % Remove diagonal values, as they are self contained in the node

[s,t,T] = setup_traffic(Tmat);      % Convert OD matrix to the s,t,T format described above

%% Pre-caluclations (Guesses)

%    Maximum flow precalculation
% Given the MF flow is relatively simple to calculate
% We make a first guess with the paths and maximum flows obtained.
% Uses the inbuilt Boykov-Kolmogorov algorithm: Computes the maximum flow by constructing two search trees associated with nodes s and t.
[MF , paths , pcosts] = maxflowpaths_st(GC,GD,s,t);

% Another set of good guesses would be in the Shortest paths.
% Currently unused.
%[nsp,csp] = getsp_s(fd,td,dd,s,t);

%% Superedges approach (To cover unfulfilled demand) Not used here.
% Define the superedges for unfulfilled demand
% This is just adding a virtual path = {s , t} for all nodes with pcost = {gamma};
% [nspSE, cspSE] =  getsuperedges_s(10000,s,t) ;


%% Solution of the minimization problem
% Set up the matrices for cplex; 
    %kopath: kommodity - path incidence matrix
    %capctr: list every link and its capacity, and the path it belongs to
    %dij: path - arc incidence matrix
    %pathcosts: the coefficients row for the objective
    %linklist: The sorted link list used in the constraints
    %[kopath,capctr,dij,pathcosts,linklist,kpath] = setuppathproblem_s([paths;nsp], [pcosts;csp], capmatrix ,s,t,T); 
    
[kopath,capctr,dij,pathcosts,linklist,kpath] = setuppathproblem_s([paths], [pcosts], capmatrix ,s,t,T);

%Solve the mutlicommodity flow problem in Cplex
[sol] = solve_MCF_s(pathcosts, dij ,kopath);

%Returns the solution part of cplex.Solution structure
results = sol.Solution; 

%% Plotting
% Results tidied up for plotting
link_res = [linklist , dij(:,1:end-1) * results.x , results.dual(1:nl)];

% Print the flows:
print_net_results(GC, link_res, xlocation, ylocation, Tmat)

%% Plot the paths?
D = zeros(n,n);
for idx = 1:length(pathcosts)
    D(kpath(idx,2), kpath(idx,3)) = D(kpath(idx,2), kpath(idx,3)) + results.x(idx) * pathcosts(idx);
end

%% Reduced costs
[nsp,csp] = getsp_mod_s(fd,td,dd,s,t,results);

%% 
% [kopath,capctr,dij,pathcosts,linklist,kpath] = setuppathproblem_s([paths;nsp], [pcosts;csp], capmatrix ,s,t,T);
% 
% %Solve the LP problem in Cplex
% [sol] = solve_MCF_s(pathcosts, dij ,kopath);
% 
% %Returns the solution part of cplex.Solution structure
% results = sol.Solution; 
% 
% % Results tidied up for plotting
% link_res = [linklist , dij(:,1:end-1) * results.x , results.dual(1:nl)];
% 
% %Print the flows:
% print_net_results(GC, link_res,xlocation,ylocation,Tmat)
% 
% %Obtain the Cost Matrix
% D = zeros(n,n);
% for idx = 1:length(pathcosts)
%     D(kpath(idx,2), kpath(idx,3)) = D(kpath(idx,2), kpath(idx,3)) + results.x(idx) * pathcosts(idx);
% end