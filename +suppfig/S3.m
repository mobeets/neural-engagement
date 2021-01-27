%% This code generates Supplemental Fig. 3A-C
% (using provided example sessions)

doSave = false;
saveDir = 'data/figures';
fnm = 'SuppFig3';
FontSize = 28;
lw = 2;
dts = io.getDates;
minTime = 1; maxTime = 20;
dtstr = '20120528';

engagement_info = engdims.getAimingEllipseAndEngagementDimensions(dtstr);
[Bs,D] = io.loadPreprocessedData(dtstr, true);

[p_vig, p_base] = engdims.targetsHelpedByEngagement(...
    engagement_info.Ysmu_anchors, ...
    engagement_info.engagement_dims_anchors, D.vfcns{2}, 'progress');
targetType = p_vig > p_base;

B = Bs(2);
B.thetaInferred = tools.computeAngles(B.vel_int, true);
B.engagement = engdims.inferEngagementGivenAim(B.thetaInferred, ...
    B.latents, engagement_info);

% skip quit trials
assert(size(trialQuitRanges,1) == numel(dts));
trsQuit = trialQuitRanges{strcmpi(dts, dtstr),B.block_index};
for ll = 1:size(trsQuit,1)
    inds = trsQuit(ll,1):trsQuit(ll,2);
    ixtr = ismember(B.trial_index_rel_to_wmp, inds);
    B.engagement(ixtr) = nan;
    B.pupilSize(ixtr,:) = nan;
end

% find high vs. low engagement trials (per target)
ixt = (B.time >= 7) & (B.time <= 20);
values = grpstats(B.engagement(ixt), B.trial_index(ixt));
atrs = unique(B.trial_index(ixt));
trgs = grpstats(B.target(ixt), B.trial_index(ixt));
lowTrials = [];
highTrials = [];
for kk = 1:numel(grps)
    ixg = trgs == grps(kk);
    med = nanmedian(values(ixg));
    lowTrials = [lowTrials; atrs(ixg & values <= med)];
    highTrials = [highTrials; atrs(ixg & values > med)];
end

% find pupil size on high and low engagement trials
P1 = {}; P2 = {};
trs = [lowTrials; highTrials];
pupilMean = nanmean(B.pupilSize(:,1));
pupilStd = nanstd(B.pupilSize(:,1));
for ii = 1:numel(trs)
    tms = B.time(B.trial_index == trs(ii));
    vs = B.pupilSize(B.trial_index == trs(ii),1);
    ix = (tms >= minTime) & (tms <= maxTime);
    
    vs = vs(ix);
    tms = tms(ix);
    ps = nan(maxTime-minTime+1,1);
    ps(ismember(minTime:maxTime, tms)) = vs;
    
    if ismember(trs(ii), highTrials)
        P1 = [P1 ps];
    elseif ismember(trs(ii), lowTrials)
        P2 = [P2 ps];
    end
end
P1 = (cell2mat(P1)' - pupilMean)/pupilStd;
P2 = (cell2mat(P2)' - pupilMean)/pupilStd;
Ps = {P1, P2};

plot.init(FontSize);

clrs = {[165, 3, 252]/255, [210, 172, 230]/255};

for ll = 1:numel(Ps)
    mu = nanmean(Ps{ll},1);
    se = nanstd(Ps{ll},[],1)/sqrt(size(Ps{ll},1));

    xs = ((minTime:maxTime)-1)*45;
    h = plot(xs, mu, '-', 'LineWidth', lw, 'Color', clrs{ll});
    plot.plotPolygonFromBounds(xs, mu-se, mu+se, plot.dullColors(h.Color, 0.2));
end
axis tight;
% mark freeze period:
plot(45*6*[1 1], ylim, '--', 'LineWidth', lw, 'Color', 'k');
xlabel('time rel. to target onset (ms)');
ylabel('pupil size (a.u.)');
set(gca, 'LineWidth', 2);
if max(ylim) > 1
    set(gca, 'YTick', [-1 0 1]);
else
    set(gca, 'YTick', [-0.5 0 0.5]);
    if min(ylim) > -0.5
        ylim([-0.5 max(ylim)]);
    end
end
if ~doSave
    title(B.datestr);
end

plot.setPrintSize(gcf, struct('width', 5, 'height', 4));
if doSave
    export_fig(gcf, fullfile(saveDir, [fnm '_' dtstr '.pdf']));
end
