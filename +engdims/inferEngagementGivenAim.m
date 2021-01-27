function [vs, inds] = inferEngagementGivenAim(aims, latents, engagement_info, engagementDimField)
% inputs:
% - aims [N x 1] - list of angles, in deg, between [0,360)
% - latents [N x K] - neural activity
% - engagement_info (struct)
% outputs:
% - vs [N x 1] - inferred engagement for each time step
% 
    if nargin < 4
        engagementDimField = 'engagement_dims';
    end
    vs = nan(numel(aims),1);
    inds = nan(numel(aims),1);
    for t = 1:size(vs,1)
        % interpolate to find aim and engagement dimensions
        angs = tools.angdiff(aims(t), engagement_info.grps_fine);
        [~,ix] = min(angs);
        vs(t) = (latents(t,:) - engagement_info.Ysmu(ix,:)) * ...
            engagement_info.(engagementDimField)(ix,:)';
        inds(t) = find(ix);
    end
    
end
