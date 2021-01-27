%% This code generates Fig. 5, Fig. 6B-C, and Supplemental Fig. 9A-B
% (using data from only two example sessions)

%% Fig. 5A-B: show example velocities over time

doSave = false;
saveDir = 'data/figures';
lw = 3;
fnm = 'Fig5A'; showPath = true; showEndpoint = true;
showFirstTrial = true;
showEllipses = false;
msz = 20;
rad = 5;
dtstr = '20120528'; trgs = 180; fnm_suffix = 'helped';
grps = tools.thetaCenters;

[Bs,D] = io.loadPreprocessedData(dtstr);
engagement_info = engdims.getAimingEllipseAndEngagementDimensions(dtstr);

[p_vig, p_base] = engdims.targetsHelpedByEngagement(...
    engagement_info.Ysmu_anchors, engagement_info.engagement_dims_anchors, ...
    D.vfcns{2}, 'progress');
vIsGood = p_vig > p_base;
trgsHelped = grps(vIsGood == 1); trgsHurt = grps(vIsGood == 0);

clrs = {[0 161 228]/255, [237 28 36]/255};
tclrs = plot.targetColors;
vigClr = [241 90 41]/255;

plot.init;

grayClr = [128 130 133]/255;
for kk = 1:numel(trgs)
    trg = trgs(kk);
    if ismember(trg, trgsHelped)
        clr = clrs{1};
    elseif ismember(trg, trgsHurt)
        clr = clrs{2};
    else
        clr = 'k';
    end
    
    trgDir = rad*[cosd(trgs(kk)) sind(trgs(kk))];
    plot([0 trgDir(1)], [0 trgDir(2)], '--', ...
        'Color', clr, 'LineWidth', lw);
end
plot(0, 0, 'k+', 'LineWidth', 3);

ns = [1 1 5 15 50 inf]; ks = [ 1 4 6  10 20];

mus1 = nan(numel(trgs), 2);
mus2 = nan(numel(trgs), 2);
for kk = 1:numel(trgs)
    trg = trgs(kk);
    if ismember(trg, trgsHelped)
        clr = clrs{1};
    elseif ismember(trg, trgsHurt)
        clr = clrs{2};
    else
        clr = 'k';
    end
    
    vdim = engagement_info.engagement_dims_anchors(grps == trg,:);
    vdimv = D.vfcns{2}(vdim);
    vdimv = 20*vdimv/norm(vdimv);
    
    B1 = Bs(1);
    B1 = io.filterTrialsByIdx(B1, B1.target == trg);
    B2 = Bs(2);
    B2 = io.filterTrialsByIdx(B2, B2.target == trg);

    vs1 = grpstats(B1.vel_wmp, B1.trial_index_rel_to_wmp_per_trg);
    vs1 = nanmean(vs1);
    vs2 = grpstats(B2.vel_wmp, B2.trial_index_rel_to_wmp_per_trg);
    [vs2,xss] = tools.smoothInPieces(vs2, ns, ks);
    
    % WMP path
    if showPath
        plot(vs2(:,1), vs2(:,2), '-', ...
            'Color', plot.dullColors(clr), 'LineWidth', lw);
    end    
    
    % engagement axis
    pA = vs1; pB = vs1 + vdimv; dp = pB-pA;
    quiver(pA(1), pA(2), dp(1), dp(2), 0, 'Color', vigClr, ...
        'MaxHeadSize', 0/norm(dp), 'LineWidth', lw);
    
    % intuitive average
    plot(vs1(1), vs1(2), 'o', 'Color', clr, ...
        'MarkerFaceColor', 0.5*ones(1,3), 'MarkerSize', msz, ...
        'LineWidth', lw);
    
    % WMP start and engdims
    if showFirstTrial
        plot(vs2(1,1), vs2(1,2), 'o', ...
            'MarkerSize', msz, 'MarkerFaceColor', 'w', ...
            'Color', clr, 'LineWidth', lw);
    end    
    if showEndpoint
        plot(vs2(end,1), vs2(end,2), 'o', ...
            'MarkerSize', msz, 'Color', 'k', ...
            'MarkerFaceColor', clr, 'LineWidth', lw);
    end
    mus1(kk,:) = vs1;
    mus2(kk,:) = vs2(1,:);
end
if showEllipses
    mus1 = [mus1; mus1(1,:)];
    mus2 = [mus2; mus2(1,:)];
    plot(mus1(:,1), mus1(:,2), 'k-');
    plot(mus2(:,1), mus2(:,2), 'r-');
end

axis equal;
axis tight;
set(gca, 'LineWidth', lw);
set(gca, 'TickDir', 'out');
axis off;
plot.setPrintSize(gcf, struct('width', 5, 'height', 4.5));
if doSave
    export_fig(gcf, fullfile(saveDir, ...
        [fnm '_' dtstr '_' fnm_suffix '.pdf']));
end

%% get engagement at each trial during Int/WMP

grps = tools.thetaCenters;
dts = io.getDates;

Trials = cell(numel(dts), numel(grps), 2);
TrialsAbs = cell(numel(dts), numel(grps), 2);
Engagements = cell(numel(dts), numel(grps), 2);
EngagementsNull = cell(numel(dts), numel(grps), 2);
EngagementsPotent = cell(numel(dts), numel(grps), 2);
ProgInt = cell(numel(dts), numel(grps), 2);
ProgWmp = cell(numel(dts), numel(grps), 2);
PupilSizes = cell(numel(dts), numel(grps), 2);
engagementIncreaseIsGoodForProgress = nan(numel(dts), numel(grps));

for ii = 1:numel(dts)
    dtstr = dts{ii};
    engagement_info = engdims.getAimingEllipseAndEngagementDimensions(dtstr);
    [Bs,D] = io.loadPreprocessedData(dtstr);
    
    % get potent component of engagement dims
    RB = D.decoders(2).fDecoder.RowM2;
    engagement_dims_wmprow = engagement_info.engagement_dims*(RB*RB');
    engagement_info.engagement_dims_wmprow = bsxfun(@times, engagement_dims_wmprow, ...
        1./tools.rowwiseNorm(engagement_dims_wmprow));
    
    % get null component of engagement dims
    NB = D.decoders(2).fDecoder.NulM2;
    engagement_dims_wmpnull = engagement_info.engagement_dims*(NB*NB');
    engagement_info.engagement_dims_wmpnull = bsxfun(@times, engagement_dims_wmpnull, ...
        1./tools.rowwiseNorm(engagement_dims_wmpnull));
    
    for jj = 1:2
        B = Bs(jj);
        B.thetaInferred = tools.computeAngles(B.vel_int, true);
        B.engagement = engdims.inferEngagementGivenAim(B.thetaInferred, ...
            B.latents, engagement_info);
        [B.engagement_potent,inds] = engdims.inferEngagementGivenAim(...
            B.thetaInferred, ...
            B.latents, engagement_info, 'engagement_dims_wmprow');
        B.engagement_null = engdims.inferEngagementGivenAim(...
            B.thetaInferred, ...
            B.latents, engagement_info, 'engagement_dims_wmpnull');
        
        for kk = 1:numel(grps)
            ixg = B.target == grps(kk);
            vs = grpstats(B.engagement(ixg), ...
                B.trial_index_rel_to_wmp_per_trg(ixg));
            vs_pot = grpstats(B.engagement_potent(ixg), ...
                B.trial_index_rel_to_wmp_per_trg(ixg));
            vs_nul = grpstats(B.engagement_null(ixg), ...
                B.trial_index_rel_to_wmp_per_trg(ixg));
            prgs_int = grpstats(B.progress_int(ixg), ...
                B.trial_index_rel_to_wmp_per_trg(ixg));
            prgs_wmp = grpstats(B.progress_wmp(ixg), ...
                B.trial_index_rel_to_wmp_per_trg(ixg));
            if isfield(B, 'pupilSize')
                pups = grpstats(B.pupilSize(ixg,1), ...
                    B.trial_index_rel_to_wmp_per_trg(ixg));
            else
                pups = nan*size(vs);
            end
            
            Trials{ii,kk,jj} = unique(B.trial_index_rel_to_wmp_per_trg(ixg));
            TrialsAbs{ii,kk,jj} = unique(B.trial_index_rel_to_wmp(ixg));
            Engagements{ii,kk,jj} = vs;
            EngagementsPotent{ii,kk,jj} = vs_pot;
            EngagementsNull{ii,kk,jj} = vs_nul;
            ProgInt{ii,kk,jj} = prgs_int;
            ProgWmp{ii,kk,jj} = prgs_wmp;
            PupilSizes{ii,kk,jj} = pups;
        end
    end
    
    % label targets where engagement increase improves behavior metrics
    vdims = engagement_info.engagement_dims_anchors;
    [p_vig, p_base] = engdims.targetsHelpedByEngagement(...
        engagement_info.Ysmu_anchors, vdims, ...
        D.vfcns{2}, 'progress');
    engagementIncreaseIsGoodForProgress(ii,:) = p_vig > p_base;    
end

[trialQuitRanges, trialQuitMasks, trialQuitCounts] = ...
    tools.findTrialsWhereMonkeyQuit(dts);

%% Fig 5-6: plot engagement/progress over trials split by engagement usefulness
% also used for Supplemental Fig. 9A-B

doSave = false;
saveDir = 'data/figures';
trialsBeforeForBaseline = 10;
dtstr = '';
mnkNm = '';
splitByTargetType = true;
N1 = 10; N2 = 75;
grpsToShow = grps;
fnm_prefix = '';
keepQuitTrials = false;
lw = 3;

plotModes = {'engagement', 'engagement_null', 'engagement_potent', ...
    'pupil', 'progress'};
for vv = 1:numel(plotModes)
    plotMode = plotModes{vv};
    plot.init(28);

    xTrials = [-(N1-1):0 1:N2];
    Xs = Trials;
    vIsGood = engagementIncreaseIsGoodForProgress;

    showZero = true;
    if strcmpi(plotMode, 'engagement')
        Ys = Engagements; ylbl = '\Delta neural engagement (a.u.)';
        meanCenter = true;
        doNormalize = true;
        yTickStep = 0.5;
        yoffset = 0.1;
        fnm = 'Fig5C';
    elseif strcmpi(plotMode, 'pupil')
        Ys = PupilSizes; ylbl = 'pupil size (a.u.)';
        meanCenter = true;
        doNormalize = false;
        yTickStep = 0.5;
        yoffset = 0.1;
        fnm = 'Fig5D';
    elseif strcmpi(plotMode, 'progress')
        Ys = ProgWmp;
        ylbl = '\Delta speed to target (mm/s)';
        meanCenter = true;
        doNormalize = false;
        yTickStep = 10;
        yoffset = 1;
        fnm = 'Fig5E';
    elseif strcmpi(plotMode, 'engagement_null')
        Ys = EngagementsNull; ylbl = '\Delta null-engagement (a.u.)';
        meanCenter = true;
        doNormalize = true;
        yTickStep = 0.5;
        yoffset = 0.1;
        fnm = 'Fig6B';
    elseif strcmpi(plotMode, 'engagement_potent')
        Ys = EngagementsPotent; ylbl = '\Delta potent-engagement (a.u.)';
        meanCenter = true;
        doNormalize = true;
        yTickStep = 0.5;
        yoffset = 0.1;
        fnm = 'Fig6C';
    else
        error('Invalid plotMode');
    end

    if ~isempty(dtstr)
        ixd = ismember(dts, dtstr);
        Ys = Ys(ixd,:,:);
        Xs = Xs(ixd,:,:);
        vIsGood = vIsGood(ixd,:);
    elseif ~isempty(mnkNm)
        ixd = io.getMonkeyDateFilter(dts, {mnkNm});
        Ys = Ys(ixd,:,:);
        Xs = Xs(ixd,:,:);
        vIsGood = vIsGood(ixd,:);
    end
    YsMu = nan(numel(xTrials), size(Xs,1), size(Xs,2));
    
    for ii = 1:size(Xs,1)
        for kk = 1:size(Xs,2)
            if ~ismember(grps(kk), grpsToShow)
                continue;
            end
            xs1 = Xs{ii,kk,1};
            xs2 = Xs{ii,kk,2};
            ys1 = Ys{ii,kk,1};
            ys2 = Ys{ii,kk,2};
            if isempty(ys1) && isempty(ys2)
                continue;
            end
            if all(isnan(ys1)) && all(isnan(ys2))
                continue;
            end

            if ~keepQuitTrials
                q1 = trialQuitMasks{ii,kk,1};
                q2 = trialQuitMasks{ii,kk,2};

                ys1(q1) = nan;
                ys2(q2) = nan;
            end
        
            if meanCenter
                if trialsBeforeForBaseline > 0
                    mu = nanmean(ys1((max(1,end-trialsBeforeForBaseline+1)):end));
                else
                    mu = nanmean(ys1);
                end
            else
                mu = 0;
            end
            if doNormalize
                sd = nanstd(ys1);
            else
                sd = 1;
            end
            ys1 = (ys1-mu)/sd;
            ys2 = (ys2-mu)/sd;
            xsc = [xs1; xs2];
            ysc = [ys1; ys2];
            ixc = ismember(xsc, xTrials);
            xsc = xsc(ixc);
            ysc = ysc(ixc);
            ixc = ismember(xTrials, xsc);

            YsMu(ixc,ii,kk) = ysc;
        end
    end
    YsMu = reshape(YsMu, size(YsMu,1), []);    
    
    Ys1 = YsMu(:,vIsGood(:) == 1);
    Ys2 = YsMu(:,vIsGood(:) == 0);    
    Yss = {Ys1, Ys2};
    clrs = {[0 161 228]/255, [237 28 36]/255};
    seClrs = {[225 244 253]/255, [253 233 233]/255};

    [~,mids] = io.dtToMnkNm(dts, false);
    mids = repmat(mids, 1, 8);
    Mnks1 = mids(vIsGood(:) == 1);
    Mnks2 = mids(vIsGood(:) == 0);

    DtInds = repmat(1:numel(dts), 8, 1)';
    Dts1 = DtInds(vIsGood(:) == 1);
    Dts2 = DtInds(vIsGood(:) == 0);
    seClr = [230 231 232]/255;

    if showZero
        plot([-N1 N2], [0 0], '--', 'LineWidth', lw, 'Color', 0.0*ones(1,3));
    end

    for ll = 1:numel(Yss)
        pts = Yss{ll};

        for jj = 1:2
            if jj == 1
                ixc = xTrials <= 0;
            else
                ixc = xTrials > 0;
            end
            xsc = xTrials(ixc)';
            ysc = pts(ixc,:);

            % plot mean and SE
            muc = nanmean(ysc,2);
            sd = nanstd(ysc,[],2);
            ns = sum(~isnan(ysc),2);
            sec = sd./sqrt(ns);
            plot(xsc, muc, '-', ...
                'Color', clrs{ll}, 'LineWidth', 3);
            plot.plotPolygonFromBounds(xsc, muc-sec, muc+sec, seClrs{ll});
        end
    end

    axis tight;
    yl = ylim;
    plot([0 0], [yl(1) yl(2)+yoffset], '--', 'LineWidth', lw, ...
        'Color', 0.8*ones(1,3));
    set(gca, 'LineWidth', lw);
    set(gca, 'TickDir', 'out');
    set(gca, 'XTick', 0:25:max(xlim));
    yticks = 0:yTickStep:max(ylim);
    if min(ylim) < 0
        yticksLeft = 0:yTickStep:-min(ylim);
        yticks = unique([-yticksLeft(end:-1:1) yticks]);
    end
    set(gca, 'YTick', yticks);
    set(gca, 'XColor', 'k', 'YColor', 'k');
    xlabel('# trials rel. to Block 2');
    ylabel(ylbl);
    plot.setPrintSize(gcf, struct('width', 5, 'height', 4.5));
    if doSave
        export_fig(gcf, fullfile(saveDir, [fnm_prefix fnm '.pdf']));
    end
end

%% Fig. 5E: learning rates for T+ and T- targets

doSave = false;
saveDir = 'data/figures';
fnm = 'Fig5F';
binSize = 10;
N2 = 75;
lw = 3;
kSmooth = 8;
pctThresh = 1.0;

% split targets as T+ vs. T-
targetType = engagementIncreaseIsGoodForProgress(:);
clrs = {[0 161 228]/255, [237 28 36]/255};

% smooth progress per target
yss = tools.cellToFixedLengthMat(ProgWmp(:,:,2), N2);
ysc = tools.smooth((1:size(yss,2)), yss, kSmooth, 1, @mean);
assert(size(ysc,1) == size(yss,1));

% scale each target so that performance is relative to max
ysc = ysc - repmat(nanmin(ysc, [], 2), 1, size(ysc,2));
ysc_scaled = ysc./repmat(nanmax(ysc, [], 2), 1, size(ysc,2));

% find trial where we first cross some percentage of max
LearningSpeeds = nan(size(yss,1), 1);
for ii = 1:size(ysc_scaled,1)
    if all(isnan(ysc_scaled(ii,:)))
        continue;
    end
    LearningSpeeds(ii) = find(ysc_scaled(ii,:) >= pctThresh*(1-eps),1);
end
Ys1 = LearningSpeeds(targetType(:) == 1);
Ys2 = LearningSpeeds(targetType(:) == 0);
maxN = max([Ys1; Ys2]);

plot.init(28);
bins = 1:binSize:(maxN+1);
if bins(end) ~= (maxN+1)
    bins = [bins maxN+1];
end

cs1 = histcounts(Ys1, bins); cs1 = cs1/sum(cs1);
cs2 = histcounts(Ys2, bins); cs2 = cs2/sum(cs2);
plot.histogramPretty(100*cs1, bins, clrs{1}, lw);
plot.histogramPretty(100*cs2, bins, clrs{2}, lw);
axis tight;
set(gca, 'XTick', [1 25:25:maxN]);
set(gca, 'YTick', [0:10:max(ylim)]);
xlabel({'# trials to reach', 'peak performance'});
ylabel('% of targets');
h = plot(nanmedian(Ys1), max(ylim)+2, ...
    'v', 'Color', clrs{1}, 'MarkerFaceColor', clrs{1}, 'MarkerSize', 10);
plot(nanmedian(Ys2), h.YData, ...
    'v', 'Color', clrs{2}, 'MarkerFaceColor', clrs{2}, 'MarkerSize', 10);
set(gca, 'LineWidth', lw);

disp(ranksum(Ys1, Ys2));
disp([nanmedian(Ys1) nanmedian(Ys2)]);

plot.setPrintSize(gcf, struct('width', 5, 'height', 4.5));
if doSave
    export_fig(gcf, fullfile(saveDir, [fnm '.pdf']));
end
