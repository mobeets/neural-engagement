function clrs_dull = dullColors(clrs, frac)
% function clrs_dull = dullColors(clrs, frac)
%     desaturates colors by the fraction 0 <= frac <= 1 (default: 0.5)
%
    if nargin < 2
        frac = 0.5;
    end
    if numel(frac) > 1 && size(frac,1) == 1
        frac = frac';
    end
    if numel(frac) == 1
        frac = repmat(frac, size(clrs,1), 1);
    end
    clrs_dull = rgb2hsv(clrs);    
    ind = 2; % modify saturation (default)
    isGrayScale = all(clrs_dull(:,1) == 0 & clrs_dull(:,2) == 0);
    if isGrayScale % if all colors are grayscale, modify value
        clrs_dull = clrs + frac*(1-clrs); return;
    end
    clrs_dull(:,ind) = frac.*clrs_dull(:,ind);
    clrs_dull = hsv2rgb(clrs_dull);
end
