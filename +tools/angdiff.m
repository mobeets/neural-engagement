function delta = angdiff(a,b)
% returns smallest angle difference between a and b (in degrees)
% answer will be between 0-180. 
    delta = min(360-mod(a-b,360), mod(a-b,360));
end
