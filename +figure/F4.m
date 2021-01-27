%% This code generates Fig. 4C (using data from only two example sessions)

%% gather engagement on first trial(s) per target or per session

dts = io.getDates;

grps = tools.thetaCenters;
Xs = cell(numel(dts), 2, numel(grps));
Ys = cell(numel(dts), 2, numel(grps));
engagementIncreaseIsGoodForProgress = nan(numel(dts), numel(grps));

for ii = 1:numel(dts)
    dtstr = dts{ii};
    engagement_info = engdims.getAimingEllipseAndEngagementDimensions(dtstr);
    [Bs,D] = io.loadPreprocessedData(dtstr);
    
    for jj = 1:2
        B = Bs(jj);
        B.engagement = engdims.inferEngagementGivenAim(B.theta, ...
            B.latents, engagement_info);
        
        % average per trial
        xs = [B.target B.trial_index_rel_to_wmp_per_trg];
        vs = grpstats(B.engagement, xs);
        xsa = unique(xs, 'rows');
        trgs = xsa(:,1);
        trs = xsa(:,2);

        % store engagement separately per target
        for kk = 1:numel(grps)
            ixc = trgs == grps(kk);
            Xs{ii,jj,kk} = trs(ixc);
            Ys{ii,jj,kk} = vs(ixc);
        end
    end
    
    % label targets where engagement increase improves behavior metrics
    vdims = engagement_info.engagement_dims_anchors;
    [p_vig, p_base] = engdims.targetsHelpedByEngagement(...
        engagement_info.Ysmu_anchors, vdims, ...
        D.vfcns{2}, 'progress');
    engagementIncreaseIsGoodForProgress(ii,:) = p_vig > p_base;
end

%% Fig. 4C: histograms of engagement increase on first trial, for T+ and T-

doSave = false;
saveDir = 'data/figures';
maxTrialNumber = 1;
trialsBeforeForBaseline = 10;
showAvgPerTarget = true;
doNormalize = true;
lw = 3;
xTickStep = 3;
binSize = 0.5;
fnm = 'Fig4C';

xlbl = '\Delta neural engagement (a.u.)';
if doNormalize
    ylbl = '% of targets';
    yOffset = 0.02;
else
    ylbl = '# targets';
    yOffset = 2;
end
Xsc = Xs; Ysc = Ys;
curEngagementIncreaseIsGood = engagementIncreaseIsGoodForProgress;

yss = nan(size(Ysc,1), size(Ysc,3));
for ii = 1:size(Ysc,1)
    for kk = 1:size(yss,2)
        ysc = Ysc{ii,2,kk};
        if isempty(ysc)
            continue;
        end
        yss(ii,kk) = nanmean(ysc(1:maxTrialNumber));
        ysc_base = Ysc{ii,1,kk};
        if trialsBeforeForBaseline > 0
            ysc_recent = ysc_base((end-trialsBeforeForBaseline+1):end);
            ysc_base_mean = nanmean(ysc_recent);
        else
            ysc_base_mean = nanmean(ysc_base);
        end
        yss(ii,kk) = (yss(ii,kk) - ysc_base_mean)/nanstd(ysc_base);
    end
end

if ~showAvgPerTarget
    yss = nanmedian(yss,2);
end
yss = yss(:);

ymn = floor(min(yss(:))); ymx = ceil(max(yss(:)));
if ymn < 0
    binsLeftOfZero = [0:binSize:-ymn];
    binsLeftOfZero = -binsLeftOfZero(end:-1:1);
else
    binsLeftOfZero = [];
end
if ymx > 0
    binsRightOfZero = [0:binSize:ymx];
else
    binsRightOfZero = [];
end
bins = unique([binsLeftOfZero binsRightOfZero]);

plot.init(28);

if showAvgPerTarget
    ixHelps = curEngagementIncreaseIsGood(:) == 1;
    ixNoHelps = curEngagementIncreaseIsGood(:) == 0;
    yscs = {yss(ixHelps), yss(ixNoHelps)};
    clrs = {[0 161 228]/255, [237 28 36]/255};
else
    yscs = {yss};
    clrs = {'k'};
end
mus = cell(numel(yscs),1);
for ll = 1:numel(yscs)
    cs = histcounts(yscs{ll}, bins);
    if doNormalize
        cs = cs/sum(cs);
    end
    mus{ll} = nanmedian(yscs{ll});
    [meds, ci_lb, ci_ub] = tools.bootstrapMedian(yscs{ll}, 10000, 0.95);
    [ci_lb median(meds) ci_ub]

    plot.histogramPretty(cs, bins, clrs{ll}, lw);
end
p = ranksum(yscs{1}, yscs{2}); % compare medians
[h,p] = kstest2(yscs{1}, yscs{2}); % compare distributions
p

axis tight;
yl = ylim;
plot([0 0], [yl(1) yl(2)+yOffset], '--', 'LineWidth', lw, ...
    'Color', 0.8*ones(1,3));
for ll = 1:numel(mus)
    plot(mus{ll}, yl(2)+yOffset, 'v', 'MarkerSize', 10, ...
        'LineWidth', 1, 'Color', clrs{ll}, ...
        'MarkerFaceColor', clrs{ll});
end
xticks = 0:xTickStep:max(xlim);
if min(xlim) < 0
    xticksLeft = 0:xTickStep:-min(xlim);
    xticks = unique([-xticksLeft(end:-1:1) xticks]);
end
set(gca, 'XTick', xticks);

set(gca, 'LineWidth', lw);
set(gca, 'TickDir', 'out');
set(gca, 'XColor', 'k', 'YColor', 'k');
if doNormalize
    % label as percent
    set(gca, 'YTickLabel', 100*get(gca, 'YTick'));
end
xlabel(xlbl);
ylabel(ylbl);
plot.setPrintSize(gcf, struct('width', 5, 'height', 4));
if doSave
    export_fig(gcf, fullfile(saveDir, [fnm '.pdf']));
end
