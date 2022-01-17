function [link_flows, link_duals, comm_duals,se_flows, D, F] = sol_handle_s(sol,dij,FTCD,nl,nn,pathcosts,kpath,s,t)

%% Function handling the solution structure provided by cplex
% grab the solution cplex throws and check
% Does .x exist? Store the duals.
% Was the answer feasible?

% F = Indicator of feasibility of solution
% D = HyperDistances matrix
% link_flows

D = zeros(nn,nn);

if strcmp(sol.Solution.statusstring,'optimal')
    
    F = 0;
    
    link_flows = [FTCD , dij(:,1:end-1) * ...
                             sol.Solution.x(1:length(pathcosts)) ];        % Flows are given for paths, so we translate to arcs and append % sol.Solution.x(1:length(pathcosts)) removes the flows on superedges
    
    link_duals = [FTCD , sol.Solution.dual(1:nl)];                         % The duals are given for the capacity constraints, then for the K=T constraint (number of paths)
    
    se_flows = [s',t',sol.Solution.x(end-length(s)+1:end)];                % Flows on the superdeges i.e: non materialized flows
    
    comm_duals = [sol.Solution.dual(nl+length(s)+1:end)];                  % Commodity duals: What one extra unit should minimum pay to be profitable
    
    % Build the hyperdistances matrix
    for idx = 1:length(pathcosts)
        
        D(kpath(idx,2), kpath(idx,3)) = D(kpath(idx,2), kpath(idx,3)) + sol.Solution.x(idx) * pathcosts(idx);
        
    end

else % These cases should only show catch cases of infeasibility, which should not occur now.
    
    if exist('results.x','var')
        
        F = 1;
        
        link_flows = [FTCD , dij(:,1:end-1) * sol.Solution.x ];
        link_duals = [FTCD , sol.Solution.dual(1:nl)];
        comm_duals = [sol.Solution.dual(nl+1:end)];
        
        D(:,:) = Inf;
        
    else
        
        F = 2;
        
        link_flows = [FTCD , zeros(length(FTCD),1) ];
        link_duals = [FTCD , zeros(length(FTCD),1) ];
        comm_duals = [zeros(length(FTCD),1)];
        
        D(:,:) = Inf;        
    end
    
end

end
    
    