    function [H] = hzsim_s(XC , YC, INT)

%% This function takes in the location of the network nodes 
%  and returns the hazard intensity at each edge/link (H)

% XC = Center of link locations (x)
% YC = Center of link locations (y)

% read a table with a hazard described by points on a grid with a certain intensity.
path=".\data\"+"hazard.xlsx";
data_hz = xlsread(path);

% Keep only points with a hazard intensity
d_hz = data_hz(data_hz(:,3)>0,:);

% Assume a hazard with impact according to  a normal distribution from its
% center.
E_hz = @(x,y) normpdf(sqrt(((x-d_hz(:,1)').^2)+(((y-d_hz(:,2)').^2)))... % use the distance between the Center and HZ as intensity
                      , 0, d_hz(:,4)');                                  % Use the sd of each focus

% Run for each edge center                  
H_base = E_hz(XC,YC);

% Special rowwise operation, for each H_base row we want to multiply by its unique intesnity.
% and then the highest value

dims = size(H_base);
H_max = zeros(dims);

for i  = 1:length(XC)
    H_max(i,:) = H_base(i,:)*d_hz(:,3);
end

H  = INT*max(H_max')';

end