%%% Compiler of the streamlined code: This example solves the MCF problem
%%% for one instance

%% Initialization
% Loads data and create graphs

run initialization.m %GC and GD are defined here
rng(42)

%% Def 1: Main properties of the network and relevant tables
nn = max(max(fd),max(td));   % Get number of nodes
nl = length(fd) ;            % Get number of links, no repeated links.


% Obtain the SP distance in the original G
DG =  distances(GD);

%% Def 2:  BASE Traffic, and origin destination pairs. Ideally come from an OD matrix known

% Option 1: User defined for specific OD
    Tmat = rand(nn,nn) * 500;
    Tmat = Tmat - diag(diag(Tmat));
    [s,t,T_o] = setup_traffic(Tmat);
    
    % This case is not returning duals!
    
    Tmat = sparse(nn,nn); Tmat(1,6) = 5000; Tmat(6,2) = 400;
    [s,t,T_o] = setup_traffic(Tmat);

% Option 2: Historical Traffic fit, ideally for each OD pair
    % if exist('TT','var') == 0
    %     HH=xlsread(".\data\data_traffic.xlsx",'B:B');
    %     [meantt, sdtt] = lognfit(TT);
    % else
    % end
    % meantt=log((mmtt^2)/sqrt(stt + mmtt^2)); 
    % sdtt = sqrt(log((stt/(mmtt^2))+1));

    % option 3    
%     Tmat = sprand(nn,nn,(123)/(1426^2)) * 500;
%     Tmat = Tmat - diag(diag(Tmat));
%     [s,t,T_o] = setup_traffic(Tmat);
    
% Option 3: Estimation with gravitational model
    % Tmat = gentraffics_s(population, DG);       % Traffic Event: al pairs of OD flows
    
% In the absence of data we assume a lognormal with E(x) = ~1, Var(x) = 0.035
mmtt = 0.01; stt = 0.15;

%% Def 3: Base rainfall data, or other historical data applicable to the links fragility function
      
%Historical Rainfall load data
if exist('HH','var') == 0

    HH=xlsread(".\data\data_rainfall.xlsx",'B:B');
    fit = gevfit(HH); k = fit(1); sigma = fit(2); mu = fit(3);

else
end

%Obtain the central location of every link (used to calculate Hij(k))
[XC,YC] = centeroflinks(xlocation,ylocation,from,to);

%% Def 4: Initial pre-processing of paths according to Traffic ODs
paths = []; pcosts=[];

% Initial guess of paths using 
% a) maximum flow 
% [MF , paths , pcosts] = maxflowpaths_st(GC,GD,s,t);  %TODO: Change path isolation algorithm

% b) shortest paths
[nsp,csp] = getsp_s(fd,td,dd,s,t);

% Join into an initial path list
paths = [nsp;paths]; pcosts = [csp;pcosts];

%% Def 5: Simulation parameters and results initialization
% number of simulations to run
N_sim = 100;
n_sim = 1; % initialization of counter

% counter of acceptance rate in MCMC
counter = 0;

% Define key parameters
loss_cost = 9999;  % Cost for the lost demand               
inf_cap = 999999;  % Higher than any individual traffic flow, maximum lost demand that can go over superedges
search_depth = 10; % Persistence of the column generation algorithm

% initialize results empty vectors
Rvec = gevrnd(k,sigma,mu,[N_sim,1]);
Tvec = lognrnd(mmtt,stt,[N_sim,1]);

Fvec = zeros(N_sim,1);
res = zeros(N_sim,1);

Pf = zeros(N_sim,1);
Pfv = zeros(N_sim,1);

EC = {};
DS = sparse(nl,N_sim); 
SE = sparse(length(s),N_sim);
LF = sparse(nl,N_sim); 
CDS=zeros(nl,1);
counter = 0;


%% Begin Simulations
tic
% Small sanity check 
% max(MF)>mean(Tmat)

% display a progress bar
h = waitbar(0,'Please wait...');

% Rainfall event

        R = Rvec(n_sim);                                     % Rainfall intensity  ~100mm
    
        H = hzsim(XC , YC, xlocation, ylocation, R, 15);     % fragility function over all links
        
        capacity_h = cap_hazard( H , capacity );             % reduced capacity of each link       

% Traffic event

        Tparam = Tvec(n_sim);                                % Draw from lognormal
        
        T = T_o * Tparam;                                    % Affect the T vector by Tparam (TODO: Tparam(s,t))   

% Write new network capacities and run MCF_CG

        [FTCD , capmatrix, distmatrix ] = adj_mats_s(fd, td, [capacity_h;capacity_h], dd);

        [link_flows, link_duals, comm_duals,se_flows, D, F, paths , pcosts] = ...
        ...
        mcmc_s_single_run_function(FTCD, nl, nn, s, t, T, capmatrix, distmatrix, paths, pcosts, loss_cost, inf_cap, search_depth);       
        
% Write results     

        res(n_sim) = abs(sum(sum(D-DG.*Tmat*Tparam))); 
        
        EC{n_sim} = D-DG.*Tmat*Tparam;
                
        DS(:,n_sim) = link_duals(:,5);
        
        LF(:,n_sim) = link_flows(:,5);
        
        SE(:,n_sim) = se_flows(:,3);
        
        Fvec(n_sim) = F;

%% 
while n_sim < N_sim
    
    % update progress bar 
    waitbar(n_sim / N_sim)
        
    % update counter
    n_sim = n_sim + 1;

    % Draw a sample of a)traffic 
    R = Rvec(n_sim);     % Rainfall intensity

    % and b) rainfall    
    Tparam = Tvec(n_sim); %

    % likelihood of proposed parameters
    p = log(gevpdf(R,k,sigma,mu)) + log(lognpdf(Tparam, mmtt, stt));

    % likelihood of previous parameters
    p1 = log(gevpdf(Rvec(n_sim-1),k,sigma,mu)) + log(lognpdf(Tvec(n_sim-1), mmtt, stt));
    
    if rand() > exp(p-p1)  % the likelihood is too low, chain stays at previous position
            
        res(n_sim) = res(n_sim-1); Rvec(n_sim) = Rvec(n_sim-1); SE(:,n_sim) = SE(:,n_sim-1); LF(:,n_sim) = LF(:,n_sim-1);
        Tvec(n_sim) = Tvec(n_sim-1); DS(:,n_sim) = DS(:,n_sim-1); Fvec(n_sim) = F; EC{n_sim} = EC{n_sim-1};
        
    else                   % likelihood high, advance with a new event and routing
            
        % Rainfall event
    
        H = hzsim(XC , YC, xlocation, ylocation, R, 15); % fragility function over all links
        
        capacity_h = cap_hazard( H , capacity );         % reduced capacity of each link 
        
        % Traffic event

        Tparam = lognrnd(mmtt,stt);                      % Draw from lognormal
        
        T = T_o * Tparam;                                % Affect the T vector by Tparam (TODO: Tparam(s,t))   
        
        % Write new network capacities and run MCF_CG
        
        [FTCD , capmatrix ] = adj_mats_s(fd, td, [capacity_h; capacity_h], dd); % Here it assumes bidirectionality!
        
        [link_flows, link_duals, comm_duals,se_flows, D, F, paths , pcosts] = ...
        ...
        mcmc_s_single_run_function(FTCD, nl, nn, s, t, T, capmatrix, distmatrix, paths, pcosts, loss_cost, inf_cap, search_depth); 
                
        
        % Write results
        
        res(n_sim) = abs(sum(sum(D-DG.*Tmat*Tparam)));
        
        EC{n_sim} = D-DG.*Tmat*Tparam;
               
        DS(:,n_sim) = link_duals(:,5);
        
        LF(:,n_sim) = link_flows(:,5);
        
        SE(:,n_sim-1) = se_flows(:,3);
        
        Fvec(n_sim) = F;                                 % Should be all zeros

    end
    
        %Empirical Pf and Variance
        counter = counter + (sum(se_flows(:,3)) > 0);
        
        Pf(n_sim) = counter / n_sim;
        % sum(res(1:n_sim) > sum(sum(DG))) /n_sim;
        
        Pfv(n_sim) = var(Pf(1:n_sim));  
        
end
%CDS = CDS./trial; % sample mean of cumulative costs
waitbar(1)
delete(h)
toc

aa = toc;
%% Save results block

% save(".\out\results.mat",'R')


%% Figures Block
color = 'w';
fsize = 10;
posize = [100,100,500,250];

% 3D Histogram of sampled events
figs(1) = figure(1);

set(gcf,'color','w');
X = [Rvec,Tvec];
hist3(X,'CDataMode','auto','FaceColor','interp')
xlabel('IM')
ylabel('Flow Demand')
title('Distribution of F/IM Sampled')
view(135,45)
hold off;

% 2D histogram of expected yearly costs
figs(2) = figure(2);  
set(gcf,'color','w');
nhist(res,'title','Distribution of Extra Km-Vehicle travelled', 'xlabel', 'Extra Vehicle-Km','ylabel','Count')
hold off;

% %histogram for the storm data
% figs(3) = figure(3);
% set(gcf,'color','w');
% bins = 0:5:max(HH);
% h = histfit(HH,20,'gev');
% %h = bar(bins,histc(HH,bins)/length(HH),'histc');hold on;
% h(1).FaceColor = [.9 .9 .9];
% % ygrid = linspace(min(HH),max(HH),100);
% % line(ygrid,gevpdf(ygrid,k,sigma,mu));
% xlabel('Maximum');
% ylabel('Probability Density');
% % xlim([min(HH) max(HH)]);


% Define the link labels correctly
lab = string([]);
for  i = 1 : length(FTCD)
    ll = string([FTCD(i,1),FTCD(i,2)]);
    lab(i) = strcat("(",ll(1),",",ll(2),")");
end

% Post event capacity sample
figs(4) = figure(4);
set(gcf,'color','w');
x = [1:length(capacity)];
vals = [capacity(x), capacity_h(x)];
b = bar(x,vals);
set(gca,'XLim',[0 length(capacity)+1],'XTick',[1:1:length(capacity)]);
xlabel('Link ID');
set(gca,'xticklabel',lab);
ylabel('Capacity');
legend('Original Capacity','Post-Event Capacity',...
       'Location','NE')

DDS = [DS{:,:}];

% Mean "cost averted by unit capacity" bargraph 
figs(5) = figure(5);
set(gcf, 'Position',  posize);
set(gcf,'color','w');
b = bar(abs(sum(DDS,2)));
set(gca,'XLim',[1 24],'XTick',[1:1:24])
xlabel('Link ID');
set(gca,'xticklabel',lab)
set(gca,'xticklabelrotation',35)
ylabel('Accumulated dual costs (veh-km)');
pbaspect([2 1 1])
set(gca,'FontSize', fsize)
...legend('Original Capacity','Post-Event Capacity','Location','NW')
    

set(gcf, 'Position',  posize);
set(gcf,'color','w');
p3 = plot(GD,'XData',xlocation,'YData',ylocation,'MarkerSize',5);
p3.EdgeCData = abs(sum(DDS,2));
colormap(jet);
pbaspect([2 1 1])
set(gca,'FontSize', fsize)
colorbar
% plot(digraph(linklist(1,:),linklist(2,:),abs(sum(DDS,2)))

% Autocorrelation of MCMC
figs(6) = figure(6);
set(gcf,'color','w');
[c,lags] = xcorr(Rvec);
stem(lags,c)
xlabel('Lag');
ylabel('Correlation');
legend('Autocorrelation of MCMC','Location','NW')

%figure7
run plotpf.m;

% % set(gcf,'color','w');
% % ecdf(rres)
% % xlabel('Performance Cost');
% % ylabel('F(x)');
% % legend('P(c<C)','Location','NW')

%figure8: Expected flow lost
loss = mean(SE,2);
heat_tab=table(s',t',full(loss));
heat_tab.Properties.VariableNames{'Var1'} = 'Origin';
heat_tab.Properties.VariableNames{'Var2'} = 'Destination';
heat_tab.Properties.VariableNames{'Var3'} = 'Loss';
heatmap(heat_tab,'Origin','Destination','ColorVariable','Loss')

%figure8: Expected relative flow lost (equity)
loss = mean(SE,2);
heat_tab=table(s',t',full(loss)./T');
heat_tab.Properties.VariableNames{'Var1'} = 'Origin';
heat_tab.Properties.VariableNames{'Var2'} = 'Destination';
heat_tab.Properties.VariableNames{'Var3'} = 'Relative_Loss';
heatmap(heat_tab,'Origin','Destination','ColorVariable','Relative_Loss')