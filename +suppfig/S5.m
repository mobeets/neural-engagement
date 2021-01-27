%% This code generates Supplemental Fig. 5
% (using data from only two example sessions)

% requirements: https://github.com/mobeets/fa
% addpath('~/code/fa');

%% fit FA and cosine tuning on trials with high vs. low engagement

dts = io.getDates;
grps = tools.thetaCenters;
Stats = [];
for ii = 1:numel(dts)
    dtstr = dts{ii}
    [Bs,D] = io.loadPreprocessedData(dtstr);
    engagement_info = engdims.getAimingEllipseAndEngagementDimensions(dtstr);
    
    % get engagement
    B = Bs(1);
    B.thetaInferred = tools.computeAngles(B.vel_int, true);
    B.engagement = engdims.inferEngagementGivenAim(B.thetaInferred, ...
        B.latents, engagement_info);

    % find high vs. low engagement trials (per target)
    V = grpstats(B.engagement, B.trial_index);
    Trgs = grpstats(B.target, B.trial_index);
    Trs = unique(B.trial_index);
    lowTrials = [];
    highTrials = [];
    for kk = 1:numel(grps)
        ixg = Trgs == grps(kk);
        prcs = prctile(V(ixg), [25 50 75]);
        lowTrials = [lowTrials; Trs(ixg & V <= prcs(2))];
        highTrials = [highTrials; Trs(ixg & V > prcs(2))];
    end

    % for z-scoring spike counts per unit
    ymu = D.FactorAnalysisParams.spikeCountMean;
    ysd = D.FactorAnalysisParams.spikeCountStd;
    
    % fit FA to trial splits to compare amount of private variance
    ndims = 10;
    Y = B.spikes;
    Y = (Y - repmat(ymu, size(Y,1), 1))./repmat(ysd, size(Y,1), 1);
    ixLow = ismember(B.trial_index, lowTrials);
    fa_low = tools.fitFA(Y(ixLow,:), ndims, false, 'showPlots', false);
    ixHigh = ismember(B.trial_index, highTrials);
    fa_high = tools.fitFA(Y(ixHigh,:), ndims, false, 'showPlots', false);
    fa_low = fa_low{1};
    fa_high = fa_high{1};
    
    % fit cosine tuning to each unit on high vs. low trials
    Y = grpstats(B.spikes, B.trial_index);
    Y = (Y - repmat(ymu, size(Y,1), 1))./repmat(ysd, size(Y,1), 1);
    ixLow = ismember(Trs, lowTrials);
    mdls_low = tools.fitCosineTuning(deg2rad(Trgs(ixLow)), Y(ixLow,:));
    m_low = ([mdls_low.r_max] - [mdls_low.r_0])/2;
    b_low = ([mdls_low.r_max] + [mdls_low.r_0])/2;
    ixHigh = ismember(Trs, highTrials);
    mdls_high = tools.fitCosineTuning(deg2rad(Trgs(ixHigh)), Y(ixHigh,:));
    m_high = ([mdls_high.r_max] - [mdls_high.r_0])/2;
    b_high = ([mdls_high.r_max] + [mdls_high.r_0])/2;
    
    % save
    clear cdata;
    cdata.datestr = B.datestr;
    cdata.block_index = B.block_index;
    cdata.fa_ndims = ndims;
    cdata.fa_low = fa_low;
    cdata.fa_high = fa_high;
    cdata.modulationDepth_low = m_low;
    cdata.modulationDepth_high = m_high;
    cdata.baselineRate_low = b_low;
    cdata.baselineRate_high = b_high;
    Stats = [Stats; cdata];
end

%% plot Supplemental Fig. 5A-B

doSave = false;
saveDir = 'data/figures';
fnm_base = 'SuppFig5A';
fnm_mod = 'SuppFig5B';
lw = 3;
msz = 5;

b_low = [Stats.baselineRate_low];
b_high = [Stats.baselineRate_high];
m_low = [Stats.modulationDepth_low];
m_high = [Stats.modulationDepth_high];
clr = 0.4*ones(1,3);

% compare baselines
plot.init(28);
scatter(b_low, b_high, msz, clr, 'filled');
vmn = floor(prctile([b_low b_high], 0.5));
vmx = ceil(prctile([b_low b_high], 99.5));
axis equal; xlim([vmn vmx]); ylim([vmn vmx]);
plot(xlim, ylim, 'k--', 'LineWidth', lw);
set(gca, 'LineWidth', lw);
set(gca, 'XTick', [0 vmx]);
set(gca, 'YTick', [0 vmx]);
xlabel({'baseline rate (a.u.),', 'low engagement trials'});
ylabel({'baseline rate (a.u.),', 'high engagement trials'});
plot.setPrintSize(gcf, struct('width', 5, 'height', 4));
if doSave
    export_fig(gcf, fullfile(saveDir, [fnm_base '.pdf']));
end

% compare modulation depths
plot.init(28);
scatter(m_low, m_high, msz, clr, 'filled')
vmx = ceil(prctile([m_low m_high], 99.5));
axis equal; xlim([0 vmx]); ylim([0 vmx]);
plot(xlim, ylim, 'k--', 'LineWidth', lw);
set(gca, 'LineWidth', lw);
set(gca, 'XTick', [0 vmx]);
set(gca, 'YTick', [0 vmx]);
xlabel({'modulation depth (a.u.),', 'low engagement trials'});
ylabel({'modulation depth (a.u.),', 'high engagement trials'});
plot.setPrintSize(gcf, struct('width', 5, 'height', 4));
if doSave
    export_fig(gcf, fullfile(saveDir, [fnm_mod '.pdf']));
end

[p1,h1,ss1] = signtest(b_low, b_high);
[p2,h2,ss2] = signtest(m_low, m_high);
[p1 p2]

%% collect FA results across sessions

neuralPts = [];
for ii = 1:numel(Stats)
    cdata = Stats(ii);
    fa_low = cdata.fa_low;
    fa_high = cdata.fa_high;
    
    Ph_low = fa_low.FactorAnalysisParams.Ph;
    Ph_high = fa_high.FactorAnalysisParams.Ph;
    L_low = fa_low.FactorAnalysisParams.L;
    L_high = fa_high.FactorAnalysisParams.L;
    Shared_low = diag(L_low*L_low');
    Shared_high = diag(L_high*L_high');
    SOT_low = Shared_low./(Shared_low + Ph_low);
    SOT_high = Shared_high./(Shared_high + Ph_high);
    
    neuralPts = [neuralPts; ii*ones(numel(Ph_low),1) ...
        Ph_low Ph_high ...
        Shared_low Shared_high ...
        SOT_low SOT_high];
end
sessionAvg = grpstats(neuralPts, neuralPts(:,1));

%% plot Supplemental Fig. 5C-D

doSave = false;
saveDir = 'data/figures';
fnm = 'SuppFig5C-D';
cPts = neuralPts;
FontSize = 28;
lw = 3;
msz = 5;
clr = 0.4*ones(1,3);

plot.init(FontSize);
nrows = 1; ncols = 2; c = 1;
plot.subplot(nrows, ncols, c, FontSize); c = c+1;
scatter(cPts(:,2), cPts(:,3), msz, clr, 'filled');
axis tight; vmx = max([xlim ylim]); axis equal;
vmn = min(prctile(cPts(:,2:3), 1)); vmx = max(prctile(cPts(:,2:3), 99));
xlim([vmn vmx]); ylim(xlim);
plot(xlim, ylim, 'k--', 'LineWidth', lw);
set(gca, 'LineWidth', lw);
set(gca, 'XTick', [0.5 1]);
set(gca, 'YTick', get(gca, 'XTick'));
xlabel({'private variance (a.u.),', 'low engagement trials'});
ylabel({'private variance (a.u.),', 'high engagement trials'});

plot.subplot(nrows, ncols, c, FontSize); c = c+1;
scatter(100*cPts(:,4), 100*cPts(:,5), msz, clr, 'filled');
axis equal; axis tight;
% xlim([0 100]); ylim([0 100]);
vmn = 100*min(prctile(cPts(:,4:5), 1)); vmx = 100*max(prctile(cPts(:,4:5), 99));
xlim([vmn vmx]); ylim(xlim);
plot(xlim, ylim, 'k--', 'LineWidth', lw);
set(gca, 'LineWidth', lw);
set(gca, 'XTick', [0 50 100]);
set(gca, 'YTick', get(gca, 'XTick'));
xlabel({'shared variance (a.u.),', 'low engagement trials'});
ylabel({'shared variance (a.u.),', 'high engagement trials'});

[p1,h1,ss1] = signtest(cPts(:,2), cPts(:,3));
[p2,h2,ss2] = signtest(cPts(:,4), cPts(:,5));
[p1 p2]

plot.setPrintSize(gcf, struct('width', 10, 'height', 4));
if doSave
    export_fig(gcf, fullfile(saveDir, [fnm '.pdf']));
end
