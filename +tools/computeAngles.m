function ths = computeAngles(vecs, doFlipV)
    if nargin < 2
        doFlipV = false;
    end
    if ~doFlipV
        vecs(:,2) = -vecs(:,2);
    end
    ths = mod(atan2d(vecs(:,2), vecs(:,1)), 360);
end
