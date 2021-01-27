function cnts = thetaCenters(n)
    if nargin < 1
        n = 8;
    end

%     rads = 0:pi/4:2*pi-pi/4;
    rads = linspace(0, 2*pi, n+1);
    rads = rads(1:end-1);
    cnts = rad2deg(rads)';
    
    if all(sqrt((cnts - round(cnts)).^2) < 1e-10)
        cnts = round(cnts);
    end % help round to nearest integer
end
