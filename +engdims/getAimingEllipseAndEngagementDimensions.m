function info = getAimingEllipseAndEngagementDimensions(dtstr, blockMode)
    if nargin < 2
        blockMode = 'pre'; % options: pre, post
    end
    
    % load session
    [B, FactorAnalysisParams] = loadSessionForFindingEngagement(dtstr, blockMode);
    
    % find aiming ellipse and engagement dims, interpolated to grps_fine
    grps_fine = tools.thetaCenters(360*8);
    [Ysmu, Ysmu_anchors] = engdims.getEllipse(B.target, B.latents, ...
        grps_fine);
    [engagement_dims, engagement_dims_anchors, vexps] = ...
        engdims.getEngagementDimensions(B.target, B.latents, ...
        FactorAnalysisParams, grps_fine);
    
    % store results
    clear info;
    info.datestr = dtstr;
    info.grps_fine = grps_fine;
    info.Ysmu = Ysmu;
    info.engagement_dims = engagement_dims;
    info.Ysmu_anchors = Ysmu_anchors;
    info.engagement_dims_anchors = engagement_dims_anchors;
    info.engagement_dim_variance_explained = vexps;
    
end

function [B, FAP] = loadSessionForFindingEngagement(dtstr, blockMode)
    [Bs,D] = io.loadPreprocessedData(dtstr);
    if strcmpi(blockMode, 'pre')
        B = Bs(1);
    elseif strcmpi(blockMode, 'post')
        B = Bs(2);
    else
        error('Invalid blockMode');
    end
    FAP = D.FactorAnalysisParams;
end
