function [cplex] = solve_MCF_s(pathcosts, dij ,kopath,gamma,inf_cap)

nk = length(kopath(:,1));

% Build model

   cplex = Cplex('mincostflow');
   cplex.DisplayFunc = [];
   cplex.Model.sense = 'minimize';
   
    % The variables are the amount of flow to send in each path

    % Populate by column: the pathcosts obtained and zeros as lower bounds,
    % the upper bound is the maximum possible flow any single link can
    % handle (given by the highest capacity found in the dij matrix)
    
    % Objective function:
    cplex.addCols(pathcosts, [], zeros(length(pathcosts),1), max(dij(:,end))*ones(length(pathcosts),1));
    
    % add superedges at the end
     cplex.addCols(repmat(gamma,nk,1),... %this is defining the cost of lost demand per unit
         [], zeros(nk,1), ...
         repmat(inf_cap,nk,1));
    
    % Now add dij^p*arc_ij<Cap constraints for each path
    % with each 
    
     % add superedges at the end? We can have irrestricted capacity but
     % need to add the columns
     
     dim = size(dij(:,1:end-1));
     
     dij = [dij(:,1:end-1), zeros(dim(1),nk) , dij(:,end);...
         zeros(nk,dim(2)), eye(nk) , repmat(inf_cap,nk,1)];
   
   for i = 1:length(dij(:,length(dij(1,:))))                                                      % Up to the last column of the dij matrix
                                
     cplex.addRows(0, dij(i,1:(length(dij(1,:)))-1), dij(i,length(dij(1,:))));                    % Add the multipliers 1 for path uses that arc, 0 for no. Upper bound is the capacity.
     
   end
      
    % Add the flow requirements constraints, OD Flow = T
    % So the first and las columns of kopath have the T value
    % And columns are paths
   
    % I will add here the superedges to the kopath matrix
    % That is by adding a final identity matrix block
    % that is for each kommodity, a new final path
    
    kopath = [kopath(:,1:end-1) eye(nk) kopath(:,end)];
    
   for i = 1:nk
    
     cplex.addRows(kopath(i,1), kopath(i,(2:end-1)),kopath(i,end));
     
   end
   
   cplex.solve();
   
end