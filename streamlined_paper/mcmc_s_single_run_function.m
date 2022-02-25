function [link_flows, link_duals, comm_duals,se_flows, D, F,...
          paths , pcosts] = mcmc_s_single_run_function(FTCD, nl, nn, s, t, T, capmatrix, distmatrix, paths, pcosts, loss_cost, inf_cap, search_depth)

% Considers initialization already ran
[kopath,~,dij,pathcosts,kpath] = setuppathproblem_s(paths, pcosts, capmatrix, s, t, T, FTCD); 

%Solve the LP problem in Cplex
[sol] = solve_MCF_s(pathcosts, dij ,kopath,loss_cost,inf_cap);
        
%Handle the solution vector, considering feasibility
[link_flows, link_duals, comm_duals, se_flows, D, F] = sol_handle_s(sol,dij,FTCD,nl,nn,pathcosts,kpath,s,t);    

%% Check if there is flow in superedges
inf_counter = 0; %set up escape condition, linked to the search depth variable

if sum(se_flows(:,3)) > 0

while sum(se_flows(:,3)) > 0 && inf_counter < search_depth % While either a) there are superedge flows or b) we are still within search depth keep on. 
    
    inf_counter = inf_counter + 1;
    % Find more paths for those o-d pairs still not satisfied
    
    [nsp_cg,csp_cg] = getsp_s_rc(s(se_flows(:,3)>0) ...
                                 ,t(se_flows(:,3)>0) ...
                                 ,link_duals, ...
                                  distmatrix);
    
    % Add columns (paths) to the problem
    paths_update = [paths;nsp_cg]; pcosts_update = [pcosts;csp_cg];
    
    % re-build the problem for cplex
    [kopath,~,dij,pathcosts,kpath] = setuppathproblem_s(paths_update, pcosts_update, capmatrix, s, t, T, FTCD); 
    
    % Solve
    [sol] = solve_MCF_s(pathcosts, dij ,kopath,loss_cost,inf_cap);
    
    se_flows_update = [s',t',sol.Solution.x(end-length(s)+1:end)]; 
    
    if isequal(se_flows_update,se_flows) % If the added paths did not improve the solution, escape the while loop
        
        inf_counter = search_depth;
        
    else
    % If the paths did change the solution, update and keep on searching
    % A better solution could be to find the paths with negative reduced costs. If there is none, then stop.
    % Here I am just running the sp inbuilt algorithm and it needs all
    % positive values.
        
    %Handle the solution vector, considering feasibility
    [link_flows, link_duals, comm_duals,se_flows, D, F] = sol_handle_s(sol,dij,FTCD,nl,nn,pathcosts,kpath,s,t);
        
    paths = paths_update;
    pcosts = pcosts_update;
    
    % F should always be = 0

    end
        
    
end

else
    % optimal solution was found, do nothing    
end