%% Define Graphs and initial inputs
clear all
clc

%Load Data
run load_data.m

%% 
%Store the initial undirected graph for capacity and distance
UGD = graph(from,to,distance);
UGC = graph(from,to,capacity);


    % Undirect all
    
    fd = [from; to];
    td = [to ; from];
    cdd = [capacity ; capacity];
    dd = [distance ; distance];
    

%Store the directed versions, for OR analysis.

    % The distance graph
    
    GD = digraph(fd,td,dd);
    
    % The capacity graph
    
    GC = digraph(fd,td,cdd);
    
    % Now calling weight for any calls the QoI of each.
    
    