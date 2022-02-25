function [R , spgf , pcsf ] = maxflowpaths_st(GC,GD,s,t) %Input the graph with capacities and distances and origin destinations

% This function yields
% The matrix of maximum flows MF
% The paths that produce the maximum flow spgf
% Costs of the paths found pcsf

%% We Attempt to get all the paths related to the maximum flow

% Define the size of all the matrices and cells by the num of nodes.

n = length(s); 

R = zeros(1,n) ; mfg = {};

% Define spg, the cell with the maximum flow paths
% Define pcs, the cost of each maximum flow path

spg = cell(n,1); pcs = cell(n,1); 

% Iterate for each origin and destination

for i = 1:n
             
            % Store the max flows and digraphs
            [R(i),mfg{i}] = GC.maxflow(s(i),t(i)); % returns the max flow R and the digraph mfg
 
            %Find all paths within the max flow digraph
            spg{i} = pathbetweennodes(mfg{i}.adjacency , s(i) , t(i));
           
            %% For each path found, find the cost:
            for p = 1 : length(spg{i})
                
                % a) Initialize the cost as zero
                pcs{i}(p,1) = 0;
                
                % For each edge in the path p:
                for k = 1 : length(spg{i}{p})-1
                               
                    
                    
                    
                    % Cycle through each path and sum each links cost
                     pcs{i}(p,1) = pcs{i}(p,1) + ...
                         GD.Edges.Weight(findedge(GD,spg{i}{p}(k),spg{i}{p}(k+1))); %Matlab2020 has a faster way for this % used to be mfg{i} 

% TODO: This probably breaks down if you have more than one edge between nodes. But the Matlab graph definition already had problems with that.
% TODO: Why not use instead the distance adjacency matrix???


                end
            
            
            end

end

%% Considering its use will be for inputting later as a full thing

spgf = {};
pcsf = {};

for i = 1:length(spg)
    spgf = [ spgf ; spg{i} ] ;
    pcsf = [ pcsf ; pcs{i} ] ;
end

%% Obtain digraph of max flows
% Get shortest path
% find path with minimum capacity
% Substract amount from the edges just obtained as SP
% re-run shortest path
% if no result, then end.