function [FTCD , capmatrix, distmatrix ] = adj_mats_s(from, to, capacity, distance)

%% This function currently does not support repeated edges!!

FTCD = [ from to capacity distance ] ;                        % Condensed Data matrix FTCD

DistList=unique(FTCD,'rows');         			              % Master List with Network Properties, filtering repeated links. Only used for the iteration

n=max(max(from),max(to));                                     % Get the number of nodes

capmatrix=sparse(n,n);                                        % Initialize the matrix of capacities
distmatrix=sparse(n,n);                                       % Initialize the matrix of capacities

    for i=1:length(DistList(:,1))

        iter=[DistList(i,1) DistList(i,2)];

        capmatrix(iter(1),iter(2))=DistList(i,3);                % Matrix Containing all capacities between nodes

        distmatrix(iter(1),iter(2))=DistList(i,4);               % Matrix Containing all capacities between nodes (not needed rn) 
    end

end