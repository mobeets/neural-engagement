function hs = plotEllipseAndEngagement(Ysmu, engagement_dims_fine, grps_fine, opts)
    if nargin < 4
        opts = struct();
    end
    defopts = struct('viewfcn', [], 'clr', 'k', 'lblnm', 'FA', ...
        'showEllipse', false, 'showTargetAverages', false, ...
        'doFill', false, 'marker', 'o', 'markerSize', 50, ...
        'vmax', 120, 'targetColor', [], ...
        'axisColor', 0.8*ones(1,3), ...
        'ignoreAvgsInHandle', false, ...
        'hideNegativeEngagement', false, 'engagementAmount', 1, 'showTip', true);
    opts = tools.setDefaultOptsWhenNecessary(opts, defopts);
    
    ysc = Ysmu;
    ysc = [ysc; ysc(1,:)];
    if ~isempty(opts.viewfcn)
        vsc = opts.viewfcn(ysc);
    else
        vsc = ysc;
    end
    hs = [];
    if opts.showEllipse
        h = plot3(vsc(:,1), vsc(:,2), vsc(:,3), '.-', 'Color', opts.clr);
        hs = [hs h];
    end    
    if ~isempty(engagement_dims_fine)
        for kk = 1:numel(grps_fine)
            y_pre = ysc(kk,:);
            if size(engagement_dims_fine,1) == 1
                vdim = engagement_dims_fine(1,:);
            else
                vdim = engagement_dims_fine(kk,:);
            end
            pt1 = y_pre + -opts.engagementAmount*vdim;
            if opts.hideNegativeEngagement
                pt1 = y_pre;
            end
            pt2 = y_pre + opts.engagementAmount*vdim;
            pts = [pt1; pt2];
            if ~isempty(opts.viewfcn)
                pts = opts.viewfcn(pts);
            end
            h = plot3(pts(:,1), pts(:,2), pts(:,3), '-', ...
                'Color', opts.axisColor);
            hs = [hs h];
            if opts.showTip
                h = plot3(pts(end,1), pts(end,2), pts(end,3), '.', ...
                    'Color', opts.axisColor, 'MarkerSize', 10);
                hs = [hs h];
            end
        end
    end
    if isempty(opts.targetColor)
        opts.targetColor = plot.targetColors;
    end
    grps = tools.thetaCenters;
    if opts.showTargetAverages
        tpts = vsc(ismember(grps_fine, grps),:);
        if opts.doFill
            h = scatter3(tpts(:,1), tpts(:,2), tpts(:,3), ...
                opts.markerSize, ...
                opts.targetColor, 'filled', opts.marker);
        else
            h = scatter3(tpts(:,1), tpts(:,2), tpts(:,3), ...
                opts.markerSize, ...
                opts.targetColor, opts.marker);
        end
        if ~opts.ignoreAvgsInHandle
            hs = [hs h];
        end
    end

    xlabel([opts.lblnm '_1']);
    ylabel([opts.lblnm '_2']);
    zlabel([opts.lblnm '_3']);
    axis equal; axis tight;
    if strcmpi(opts.lblnm(1), 'v')
        h = plot3(0, 0, 0, 'k+');
        hs = [hs h];
        if ~isnan(opts.vmax)
            xlim(opts.vmax*[-1 1]);
            ylim(xlim);
        end
    end

end
