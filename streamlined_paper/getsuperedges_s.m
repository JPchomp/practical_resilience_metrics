function [nsp,csp] = getsuperedges_s(gamma,s,t) 

idx = length(s);

csp=cell(idx,1);
nsp=cell(idx,1);

for k = 1:idx
        
    if s(k) == t(k)
        
    else
    
    nsp{k,1} = [s(k), t(k)];
    csp{k,1} = gamma;
    
    end
    
end
