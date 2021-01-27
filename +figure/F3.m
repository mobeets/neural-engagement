%% This code generates Fig. 3B, and Supplemental Fig. 7

%% Fig. 3B

doSave = false;
saveDir = 'data/figures';
fnm = 'Fig3B';
grps = tools.thetaCenters;
dts = {'20120528'};

lw = 1;
nFinePoints = 32;

for ii = 1:numel(dts)
    dtstr = dts{ii};
    engagement_info = engdims.getAimingEllipseAndEngagementDimensions(dtstr);
    [Bs,D] = io.loadPreprocessedData(dtstr);
    B = Bs(1);

    for mm = 1:3 % toggle to show factors vs. velocities
        if mm == 1
            showFactors = true;
        elseif mm == 2
            showFactors = false;
            showWmpVels = false;
        elseif mm == 3
            showFactors = false;
            showWmpVels = true;
        end
        
        if showFactors
            vmx = nan;
        else
            vmx = 120;
        end
        popts = struct('viewfcn', [], 'clr', 0.8*ones(1,3), 'lblnm', 'FA', ...
            'showTargetAverages', true, 'vmax', vmx, 'showEllipse', true, ...
            'markerSize', 450, 'engagementAmount', 2, 'doFill', true, ...
            'showTip', false);
        if ~showFactors
            if showWmpVels
                popts.viewfcn = @(z) [D.vfcns{2}(z) 0*D.vfcns{2}(z)];
            else
                popts.viewfcn = @(z) [D.vfcns{1}(z) 0*D.vfcns{1}(z)];
            end
        end

        grps_coarse = tools.thetaCenters(nFinePoints);
        grps_fine = engagement_info.grps_fine;
        [~,ixCoarse] = min(pdist2(grps_coarse, grps_fine),[],2);
        Ysmu = engagement_info.Ysmu(ixCoarse,:);
        vdims = engagement_info.engagement_dims(ixCoarse,:);

        plot.init;
        hs = engdims.plotEllipseAndEngagement(Ysmu, vdims, ...
            grps_coarse, popts);
        for kk = 1:numel(hs)
            hs(kk).LineWidth = lw;
        end

        popts.targetColor = 0.7*ones(1,3);
        popts.markerSize = 0.4*popts.markerSize;
        popts.showEllipse = false;
        popts.axisColor = [241 90 41]/255;
        popts.engagementAmount = 3;
        popts.hideNegativeEngagement = true;
        popts.showTip = false;
        hs = engdims.plotEllipseAndEngagement(engagement_info.Ysmu_anchors, ...
            engagement_info.engagement_dims_anchors, grps, popts);
        for kk = 1:numel(hs)
            hs(kk).LineWidth = lw;
        end

        if ~isnan(vmx)
            plot(-vmx*[1 1], vmx*[-1 1], 'k-', 'LineWidth', lw);
            plot(vmx*[1 1], vmx*[-1 1], 'k-', 'LineWidth', lw);
            plot(vmx*[-1 1], vmx*[1 1], 'k-', 'LineWidth', lw);
            plot(vmx*[-1 1], -vmx*[1 1], 'k-', 'LineWidth', lw);
            axis tight; axis equal;
            xlim([-vmx vmx]); ylim(xlim);
            plot(0, 0, '+', 'Color', [0 165 80]/255, ...
                'MarkerSize', 20, 'LineWidth', lw);    
            mnm = io.dtToMnkNm(dtstr);
            text(-45, 0.87*vmx, [mnm(1) dtstr], 'FontSize', 33);
        end
        if showFactors
            zmx = 2;
            plot3([0 zmx], [0 0], [0 0], 'k-', 'LineWidth', lw);
            plot3([0 0], [0 -zmx], [0 0], 'k-', 'LineWidth', lw);
            plot3([0 0], [0 0], [0 zmx], 'k-', 'LineWidth', lw);
            view(60, 24);
        end
        axis off;
        if doSave
            export_fig(gcf, fullfile(saveDir, ...
                [fnm '_' num2str(showFactors) '_' dtstr '.pdf']));
        end
    end
end
