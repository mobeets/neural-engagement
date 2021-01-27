function doFlip = flipEngagementDimsUsingAdjacentDims(engagement_dims)
% given each dim, find the previous and next dim
% if these dims would both be closer were we to flip the current dim,
% then we know to flip it
    
    doFlip = nan(size(engagement_dims,1),1);
    for kk = 1:size(engagement_dims,1)
        % consider the current orientation, and the flipped one
        vdim = engagement_dims(kk,:);
        ypos = vdim;
        yneg = -vdim;

        kprev = kk-1;
        if kprev < 1
            kprev = size(engagement_dims,1);
        end
        knext = kk+1;
        if knext > size(engagement_dims,1)
            knext = 1;
        end

        % compare to prev dim
        yadj = engagement_dims(kprev,:);
        prevSaysFlip = norm(yneg-yadj) < norm(ypos-yadj);

        % compare to next dim
        yadj = engagement_dims(knext,:);
        nextSaysFlip = norm(yneg-yadj) < norm(ypos-yadj);

        doFlip(kk) = prevSaysFlip & nextSaysFlip;
    end
end
