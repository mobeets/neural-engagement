%% This code generates Fig. 2A,B,D,E,H, and Supplemental Fig. 1
% (using data from only two example sessions)

%% Fig. 2A-B: plot example neural activity over trials to same target

doSave = false;
saveDir = 'data/figures';
grps = tools.thetaCenters;

dtstr = '20120528'; trg = 180;

engagement_info = engdims.getAimingEllipseAndEngagementDimensions(dtstr);
[Bs,D] = io.loadPreprocessedData(dtstr);

for vv = 1:2
    if vv == 1
        blkIndToPlot = 1; plotFirstK = inf; msz = 8; outlineFirstTrial = false;
    else
        blkIndToPlot = 2; plotFirstK = 20; msz = 10; outlineFirstTrial = true;
    end
    
    % plot example target's engagement per trial, in factor space
    plot.init;

    if ~isnan(plotFirstK) && ~isinf(plotFirstK)
        clrs = gray(plotFirstK+10);
        clrs = clrs(5:(end-5),:);
    end
    meanColor = 0.8*ones(1,3);
    axisColor = [241 90 41]/255;

    for jj = 1:2
        B = Bs(jj);
        B = io.filterTrialsByIdx(B, B.target == trg);
        B.engagement = engdims.inferEngagementGivenAim(B.theta, B.latents, engagement_info);

        % average per trial
        vs = grpstats(B.engagement, B.trial_index_rel_to_wmp_per_trg);
        zs = grpstats(B.latents, B.trial_index_rel_to_wmp_per_trg);
        trgs = grpstats(B.target, B.trial_index_rel_to_wmp_per_trg);
        trs = unique(B.trial_index_rel_to_wmp_per_trg);

        % find mean and engagement dim
        if jj == 1
            coeff = pca(zs); vdim = var(zs(:,1))*coeff(:,1);
            mu = nanmean(zs);
        end

        % only plot first K trials
        if ~isnan(plotFirstK) && ~isinf(plotFirstK)
            vs = vs(1:plotFirstK);
            zs = zs(1:plotFirstK,:);
            trgs = trgs(1:plotFirstK);
            trs = trs(1:plotFirstK);
        else
            clrs = 0.8*ones(numel(trs), 3);
        end

        if jj ~= blkIndToPlot
            continue;
        end    

        % plot scatter
        for kk = 1:size(zs,1)
            if outlineFirstTrial && kk == 1
                plot3(zs(kk,1), zs(kk,2), zs(kk,3), 'o', ...
                    'Color', 'k', 'MarkerFaceColor', clrs(kk,:), ...
                    'MarkerSize', msz);
            else
                plot3(zs(kk,1), zs(kk,2), zs(kk,3), 'o', ...
                    'Color', clrs(kk,:), 'MarkerFaceColor', clrs(kk,:), ...
                    'MarkerSize', msz);
            end
        end

        % plot mean and engagement axis
        plot3(mu(1) + vdim(1)*[-1 1], ...
            mu(2) + vdim(2)*[-1 1], mu(3) + vdim(3)*[-1 1], ...
            '-', 'Color', axisColor, 'LineWidth', 3);
        scatter3(mu(1), mu(2), mu(3), 200, meanColor, 'o', 'filled');
        scatter3(mu(1), mu(2), mu(3), 200, 'k', 'o', 'LineWidth', 2);

    end
    if ~isnan(plotFirstK) && ~isinf(plotFirstK)
        colormap(clrs);
    end
    axis equal;
    xlabel('z_1'); ylabel('z_2'); zlabel('z_3');
    axis off;

    if doSave
        fnm = [dtstr '_Blk' num2str(blkIndToPlot) '_' ...
            num2str(trg) '_trs1-' num2str(plotFirstK)];
        fnm = fullfile(saveDir, [fnm '.pdf']);
        export_fig(gcf, fnm);
    end
end

%% gather engagement for all sessions, broken up by pauses and block changes

dts = io.getDates;
grps = tools.thetaCenters;

trialGapThreshMins = 1.5;
onlyKeepMaxPause = true;
Tms = cell(numel(dts), 2, 3);
Trgs = cell(numel(dts), 2, 3);
TrialsAbs = cell(numel(dts), 2, 3); % absolute
TrialsRel = cell(numel(dts), 2, 3); % relative to event change
InCors = cell(numel(dts), 2, 3);
Ys = cell(numel(dts), 2, 3);
for ii = 1:numel(dts)
    dtstr = dts{ii};
    engagement_info = engdims.getAimingEllipseAndEngagementDimensions(dtstr);
    [Bs,D] = io.loadPreprocessedData(dtstr);
    for jj = 1:2
        B = Bs(jj);
        B.engagement = engdims.inferEngagementGivenAim(B.theta, ...
            B.latents, engagement_info);
        
        % average per trial
        vs = grpstats(B.engagement, B.trial_index);
        trgs = grpstats(B.target, B.trial_index);
        incors = grpstats(~B.isCorrect, B.trial_index);
        trs = unique(B.trial_index);
        tms = tools.trialTimesToSeconds(B.trialTime)/60; % minutes

        % average per target during Int control
        if jj == 1
            mus = nan(numel(grps),1);
            sds = nan(numel(grps),1);
            for kk = 1:numel(grps)
                vsc = vs(trgs == grps(kk));
                mus(kk) = nanmean(vsc);
                sds(kk) = nanstd(vsc);
            end
            minTrialTime = min(tms);
        end    
        for kk = 1:numel(grps)
            ix = trgs == grps(kk);
            vs(ix) = (vs(ix) - mus(kk))/sds(kk);
        end
        tms = tms - minTrialTime;

        % plot gaps in time separately (only keep max pause)
        trialGaps = diff(tms);
        trialGapInds = find(trialGaps > trialGapThreshMins);
        if numel(trialGapInds) > 2 && onlyKeepMaxPause
            % keep end, and longest middle
            [~,maxind] = max(trialGaps(trialGapInds(1:end-1)));
            trialGapInds = trialGapInds([maxind end]);
        end
        if numel(tms) ~= numel(vs)
            warning(dtstr, 'Trial timing error.');
            continue;
        end
        inds = [1; trialGapInds+1; numel(tms)+1];
        for kk = 1:(numel(inds)-1)
            cinds = inds(kk):(inds(kk+1)-1);
            tmsc = tms(cinds); vsc = vs(cinds);
            trsc = trs(cinds); trgsc = trgs(cinds); incorsc = incors(cinds);
            Tms{ii,jj,kk} = tmsc;
            Ys{ii,jj,kk} = vsc;
            Trgs{ii,jj,kk} = trgsc;
            TrialsAbs{ii,jj,kk} = trsc;
            % 1-indexed for this event:
            TrialsRel{ii,jj,kk} = trsc - min(trsc) + 1;
            InCors{ii,jj,kk} = incorsc;
        end
    end
end

%% Fig. 2D (and Supplemental Fig. 1): plot example session

doSave = false;
saveDir = 'data/figures';
fnm = 'Fig2D';
tickTimeMins = 5;
tickTrialCount = 200;
lw = 3;
xAxisShowsMinutes = false;
skipIncorrects = true;
kSmooth = 10;
yMax = 2;
vigClr = [241 90 41]/255;
FontSize = 28;
dtstr = '20120528';

xs = Tms(strcmpi(dts, dtstr),:,:);
trs = TrialsAbs(strcmpi(dts, dtstr),:,:);
ys = Ys(strcmpi(dts, dtstr),:,:);
incor = InCors(strcmpi(dts, dtstr),:,:);

plot.init(FontSize);
timeRefs = [];
for jj = 1:2
    for kk = 1:size(xs,3)
        xsc = xs{1,jj,kk};
        trsc = trs{1,jj,kk};
        timeRefs = [timeRefs; xsc trsc];
        ysc = ys{1,jj,kk};
        incorc = incor{1,jj,kk};
        if skipIncorrects
            ysc(incorc == 1) = nan;
        end
        if jj == 2 && kk == 1
            pertStart = min(xsc);
        end
        if kSmooth > 0
            [ysc,xsc] = tools.smooth(xsc,ysc,kSmooth);
        end
        plot(xsc, ysc, 'LineWidth', lw, 'Color', vigClr);
    end
end

% find closest trial for each clock time interval
timeTicks = 0:tickTimeMins:max(xlim);
set(gca, 'XTick', timeTicks);
if xAxisShowsMinutes
    xlabel('clock time (min)');
else
    trialTicks = 0:tickTrialCount:max(timeRefs(:,2));
    ds = pdist2(trialTicks', timeRefs(:,2));
    [~,ix] = min(ds,[],2);
    timeTicks = timeRefs(ix,1);
    set(gca, 'XTick', timeTicks);
    set(gca, 'XTickLabel', trialTicks);
    xlabel('# trials');
end

axis tight;
xl = xlim; xlim([xl(1)-0.3 xl(2)]);
ylim([min(ylim) yMax]);
plot(pertStart*[1 1], ylim, 'k-', 'LineWidth', lw);
plot(xlim, [0 0], '--', 'LineWidth', lw, 'Color', 'k');
ylabel('neural engagement (a.u.)');
set(gca, 'LineWidth', lw);
set(gca, 'TickDir', 'out');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'YTick', unique(round(get(gca, 'YTick'))));
plot.setPrintSize(gcf, struct('width', 12, 'height', 5));

if doSave
    if xAxisShowsMinutes
        export_fig(gcf, fullfile(saveDir, [fnm '_' dtstr '_mins.pdf']));
    else
        export_fig(gcf, fullfile(saveDir, [fnm '_' dtstr '.pdf']));
    end
end

%% Fig. 2E: plot average across sessions

doSave = false;
saveDir = 'data/figures';
fnm = 'Fig2E';
maxTrialNumber = 50;
xTickStep = 25;
yRange = [-0.57 1.1];
splitTypes = {'start', 'pause', 'pert'};
% splitTypes = {'combined'};
lw = 3;
kSmooth = 3;
trialsBeforeForBaseline = 0;
keepIncorrects = true;

meanClr = [240 90 41]/255;
seClr = [230 231 232]/255;

for ll = 1:numel(splitTypes)
    splitType = splitTypes{ll};
    
    maxn = maxTrialNumber;
    yss = nan(size(Ys,1), maxTrialNumber);

    for ii = 1:size(Ys,1)        
        ysc_base = 0;
        if strcmpi(splitType, 'start')
            ysc = Ys{ii,1,1};
            if ~keepIncorrects
                incors = InCors{ii,1,1};
                ysc(incors == 1) = nan;
            end
        elseif strcmpi(splitType, 'pause')
            ysc = Ys{ii,1,2};
            if ~keepIncorrects
                incors = InCors{ii,1,2};
                ysc(incors == 1) = nan;
            end
            if trialsBeforeForBaseline > 0
                ysc_base = Ys{ii,1,1};
                if isempty(ysc_base)
                    ysc_base = 0;
                else
                    ysc_base = ysc_base((end-trialsBeforeForBaseline+1):end);
                    ysc_base = nanmedian(ysc_base);
                end                
            end
        elseif strcmpi(splitType, 'pert')
            ysc = Ys{ii,2,1};
            if ~keepIncorrects
                incors = InCors{ii,2,1};
                ysc(incors == 1) = nan;
            end
            if trialsBeforeForBaseline > 0
                ysc_base = Ys{ii,1,3};
                if isempty(ysc_base)
                    ysc_base = Ys{ii,1,2};
                end
                if isempty(ysc_base)
                    ysc_base = 0;
                else
                    ysc_base = ysc_base((end-trialsBeforeForBaseline+1):end);
                    ysc_base = nanmedian(ysc_base);
                end 
            end
        elseif strcmpi(splitType, 'combined')
            ysc_a = Ys{ii,1,1};
            ysc_b = Ys{ii,1,2};
            ysc_c = Ys{ii,2,1};
            n_min = min([numel(ysc_a) numel(ysc_b) numel(ysc_c)]);
            ysc = [ysc_a(1:n_min) ysc_b(1:n_min) ysc_c(1:n_min)];
            ysc = nanmean(ysc,2);
            if ~keepIncorrects
                error('Cannot used ''combined'' and ignore incorrects');
            end
            if trialsBeforeForBaseline > 0
                error('Cannot used ''combined'' with trialsBeforeForBaseline > 0');
            end
        else
            error('Invalid splitType');
        end

        ysc = ysc - ysc_base;
        if kSmooth > 0
            ysc = tools.smooth(ysc, kSmooth);
        end
        if numel(ysc) > maxTrialNumber
            ysc = ysc(1:maxTrialNumber);
        end
        yss(ii,1:numel(ysc)) = ysc;
    end

    plot.init(28);
    yss = yss(:,1:maxTrialNumber);
    mu = nanmean(yss);
    sd = nanstd(yss);
    n = sum(~isnan(yss));
    se = sd./sqrt(n);
    h = plot(mu, 'LineWidth', lw, 'Color', meanClr);
    plot.plotPolygonFromBounds(1:numel(mu), mu-se, mu+se, seClr);

    xlabel('# trials');
    ylabel('neural engagement (a.u.)');
    set(gca, 'TickDir', 'out');
    set(gca, 'YTick', unique(round(get(gca, 'YTick'))));
    set(gca, 'XTick', sort([1 0:xTickStep:max(xlim)]));
    set(gca, 'LineWidth', lw);
    axis tight;
    xlim([1 max(xlim)]);
    if ~isempty(yRange)
        ylim(yRange);
    end
    plot(xlim, [0 0], '--', 'LineWidth', lw, 'Color', 'k');
    set(gca, 'XColor', 'k', 'YColor', 'k');
    if ll > 1
        set(gca, 'YColor', 'w');
    end
    plot.setPrintSize(gcf, struct('width', 5, 'height', 5));

    if doSave
        export_fig(gcf, fullfile(saveDir, [fnm '_' splitType '.pdf']));
    end
end

%% Fig. 2G: Pearson's correlation between pupil size and neural engagement
% also Supplemental Fig. 3D-E

bind = 2;
doSave = false;
saveDir = 'data/figures';
fnm = 'Fig2G';
minTrialCount = 200;
kSmooth = 30;
minTime = 7; maxTime = 20;

dts = io.getDates;
grps = tools.thetaCenters;
if ~exist('trialQuitRanges', 'var')
    [trialQuitRanges, trialQuitMasks, trialQuitCounts] = ...
        tools.findTrialsWhereMonkeyQuit(dts);
end

Pts = cell(numel(dts), 2, 2);
for ii = 1:numel(dts)
    dtstr = dts{ii};
    engagement_info = engdims.getAimingEllipseAndEngagementDimensions(dtstr);
    [Bs,D] = io.loadPreprocessedData(dtstr, true);
    
    for jj = 1:2
        B = Bs(jj);
        if ~isfield(B, 'pupilSize')
            continue;
        end
        ixt = (B.time >= minTime) & (B.time <= maxTime);
        
        % find neural engagement
        B.thetaInferred = tools.computeAngles(B.vel_int, true);
        B.engagement = engdims.inferEngagementGivenAim(B.thetaInferred, ...
            B.latents, engagement_info);
        
        % z-score engagement per target
        V = grpstats(B.engagement(ixt), B.trial_index(ixt));
        trgs = grpstats(B.target(ixt), B.trial_index(ixt));
        for kk = 1:numel(grps)
            ixg = trgs == grps(kk);
            V(ixg) = (V(ixg) - nanmean(V(ixg)))/nanstd(V(ixg));
        end
        
        % find pupil size
        ixt = (B.time >= minTime) & (B.time <= maxTime);
        P = grpstats(B.pupilSize(ixt,1), B.trial_index(ixt), @nanmedian);
        assert(numel(P) == numel(V), ...
            'pupil size and engagement have different trials');
        
        % skip quit trials
        assert(size(trialQuitRanges,1) == numel(dts));
        trsQuit = trialQuitRanges{ii,jj};
        for ll = 1:size(trsQuit,1)
            inds = trsQuit(ll,1):trsQuit(ll,2);
            P(inds) = nan;
            V(inds) = nan;
        end
        
        % skip any sessions without enough trials
        if sum(~isnan(P)) < minTrialCount
            continue;
        end
        
        Pts{ii,jj,1} = P;
        Pts{ii,jj,2} = V;
    end
end

rhos = nan(size(Pts,1), 2);
for jj = 1:2
    for ii = 1:size(Pts,1)
        P = Pts{ii,jj,1};
        V = Pts{ii,jj,2};
        if isempty(P) || isempty(V) || all(isnan(P))
            continue;
        end

        % smoothing
        if kSmooth > 1
            P = tools.smooth((1:numel(P))', P, kSmooth);
            V = tools.smooth((1:numel(V))', V, kSmooth);
        end
        
        % find correlation
        rhos(ii,jj) = corr(P, V, 'rows', 'complete');
    end
end
cpts = rhos(:,bind);

plot.init(20);
rng('default');
rng(127); % for reproducing random x placement
xs = nan(size(cpts));
bins = prctile(cpts, 0:10:100);
for ii = 1:(numel(bins)-1)
    ix = cpts >= bins(ii) & cpts <= bins(ii+1);
    cx = sum(ix)*rand(sum(ix),1)/1.2;
    xs(ix) = cx - mean(cx);
end

mnks = io.monkeyNames;
mkrs = {'o', 'v', 'x'};
for mm = 1:numel(mnks)
    ixm = io.getMonkeyDateFilter(dts, mnks(mm))';
    plot(xs(ixm), cpts(ixm), mkrs{mm}, ...
        'LineWidth', 1.5, 'Color', 0.5*ones(1,3), ...
        'MarkerFaceColor', 'w', 'MarkerSize', 5);
end

[meds, ci_lb, ci_ub] = tools.bootstrapMedian(cpts, 10000, 0.95);
mu = median(meds);
prcs = [ci_lb, ci_ub];
plot(0*[1 1], [prcs(1) prcs(2)], 'k-', 'LineWidth', 3);
plot(0, mu, 'ko', 'LineWidth', 3, ...
    'MarkerFaceColor', 'w', 'MarkerSize', 10);    

axis tight;
xl = xlim;
xlim([xl(1)-0.2 xl(2)+0.2]);
ylim([-max(abs(ylim)) max(abs(ylim))]);
set(gca, 'XTick', []);
set(gca, 'YTick', -1:0.5:1);
set(gca, 'TickDir', 'out');
set(gca, 'LineWidth', 2);
set(gca, 'XColor', 'w', 'YColor', 'k');
ylim([-1 1]);

plot.setPrintSize(gcf, struct('width', 1.0, 'height', 3));
if doSave
    export_fig(gcf, fullfile(saveDir, ...
        [fnm '_Block' num2str(bind) '.pdf']));
end

%% Fig. 2H: find var. explained by engagement

doSave = false;
saveDir = 'data/figures';
fnm = 'Fig2H';

dts = io.getDates;
grps = tools.thetaCenters;
pct_var_exp_per_trg = nan(numel(dts), numel(grps), 2);
pct_var_exp_across_trg = nan(numel(dts), 2, 2);

for ii = 1:numel(dts)
    dtstr = dts{ii};
    engagement_info = engdims.getAimingEllipseAndEngagementDimensions(dtstr);
    [Bs,D] = io.loadPreprocessedData(dtstr, true);
    B = Bs(1);
    
    ixt = (B.time >= 3) & (B.time <= 8); % pre-feedback window
    Zs = grpstats(B.latents(ixt,:), B.trial_index(ixt));
    trgs = grpstats(B.target(ixt), B.trial_index(ixt));
    
    vs = nan(numel(trgs), 1);
    mus = nan(numel(grps), size(Zs,2));
    for kk = 1:numel(grps)
        ixc = trgs == grps(kk);
        vdim = engagement_info.engagement_dims_anchors(kk,:);
        zs = Zs(ixc,:);
        mus(kk,:) = nanmean(zs);
        vs(ixc) = bsxfun(@minus, zs, mus(kk,:)) * vdim';
        pct_var_exp_per_trg(ii,kk) = var(vs(ixc))/sum(var(zs));
    end
    pct_var_exp_across_trg(ii) = var(vs)/sum(var(Zs));
end

%
pct_vars = {pct_var_exp_across_trg, pct_var_exp_per_trg};
nms = {'Total variance', 'Variance per target'};

plot.init;
for ll = 1:numel(pct_vars)
    pct_var = pct_vars{ll};

    % plot mean percent var explained, per monkey:
    disp('-------');
    disp(['% variance explained by engagement, ' nms{ll}]);
    psc = 100*pct_var;
    psc = psc(~isnan(psc));
    mu = nanmean(psc);
    sd = nanstd(psc);
    se = sd/sqrt(numel(psc));
    disp(['All: ' sprintf('%0.2f +/- %0.2f', [mu sd])]);

    xs = ll + (rand(numel(psc),1)-0.5)/2;
    ys = psc;
    plot(xs, ys, 'o', 'LineWidth', 1, 'Color', 0.5*ones(1,3), ...
        'MarkerFaceColor', 'w', 'MarkerSize', 5);
    
    prcs = prctile(psc, [25 75]);
    mu = nanmedian(psc);
    plot(ll*[1 1], [prcs(1) prcs(2)], 'k-', 'LineWidth', 2);
    plot(ll, mu, 'ko', 'LineWidth', 2, ...
        'MarkerFaceColor', 'w', 'MarkerSize', 10);  
end

set(gca, 'XTick', 1:numel(pct_vars));
set(gca, 'XTickLabel', nms);
set(gca, 'XTickLabelRotation', 45);

axis tight;
xl = xlim;
xlim([xl(1)-1 xl(2)+1]);
ylim([0 100]);
% set(gca, 'XTick', []);
set(gca, 'TickDir', 'out');
set(gca, 'XColor', 'k', 'YColor', 'k');
ylabel({'% variance explained by', 'neural engagement axes'});

plot.setPrintSize(gcf, struct('width', 2.5, 'height', 4));

if doSave
    export_fig(gcf, fullfile(saveDir, [fnm '.pdf']));
end
