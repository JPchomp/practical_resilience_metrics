function [nsp,csp] = getsp_s_rc(s,t,link_duals,distmatrix) 

%% this function gets the shortest path, but also handles negative...
% weights which are an issue in matlabs SP. 
% 1) Sum M to all edges
% 2) Find SP 
% 3) Restore original sum as the csp = cost of spath (Cost - N*M) N = num
% of edges = num of nodes - 1

csp={};
nsp={};

m = max(link_duals(:,4)+abs(link_duals(:,5)))*100;               % Big number considering the maximum, and multiplied. Could be the case however that given a VERY negative number this is not enough!

Gtem=digraph(link_duals(:,1),...                                 % from  
             link_duals(:,2),...                                 % to
             link_duals(:,4) - link_duals(:,5));                 % Create graph with positive edges
         
         %% for the simple 2 od case i am not getting dual values!

s = full(s); t = full(t);         
         
for i = 1:length(s)
    
    if s(i) == t(i)
        
    else
    
    [nsp{i,1},~] = Gtem.shortestpath(s(i),t(i));                        % The path is already stored here  
    
    csp_temp = 0;
    
    for k = 1:length(nsp{i,1})-1
        csp_temp = csp_temp + distmatrix( nsp{i,1}(k) , nsp{i,1}(k+1)); % need to calculate like this due to the graph being with duals
    end
    
    csp{i,1} = csp_temp;
    
    end
end
