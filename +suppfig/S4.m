%% This code generates Supplemental Fig. 4
% (using data from only two example sessions)

%% find mean # of units that increased given increase in vigor dim

dts = io.getDates;

grps = tools.thetaCenters;
eng_neuron_impacts = nan(numel(dts), numel(grps));
eng_neuron_impacts_null = nan(numel(dts), numel(grps));

for ii = 1:numel(dts)
    dtstr = dts{ii};
    engagement_info = engdims.getAimingEllipseAndEngagementDimensions(dtstr);
    [~,D] = io.loadPreprocessedData(dtstr);
    L = (engagement_info.engagement_dims_anchors / ...
        D.FactorAnalysisParams.spikeRot) * D.FactorAnalysisParams.L';
    pcts = mean(sign(L) == 1,2);
    pcts(pcts < 0.5) = 1 - pcts(pcts < 0.5);
    eng_neuron_impacts(ii,:) = pcts;
    
    % null distribution
    dims_null = rand(size(engagement_info.engagement_dims_anchors))-0.5;
    dims_null = bsxfun(@times, dims_null, 1./tools.rowwiseNorm(dims_null));
    L = (dims_null/D.FactorAnalysisParams.spikeRot) * D.FactorAnalysisParams.L';
    pcts = mean(sign(L) == 1,2);
    pcts(pcts < 0.5) = 1 - pcts(pcts < 0.5);
    eng_neuron_impacts_null(ii,:) = pcts;
end

[meds, ci_lb, ci_ub] = tools.bootstrapMedian(eng_neuron_impacts_null(:), 10000);
100*[ci_lb median(meds) ci_ub] % 59.6785   60.8946   62.5000

[meds, ci_lb, ci_ub] = tools.bootstrapMedian(eng_neuron_impacts(:), 10000);
100*[ci_lb median(meds) ci_ub] % 97.6471   97.6744   97.7273

%% plot distribution of the above

doSave = false;
saveDir = 'data/figures';
fnm = 'SuppFig4';
yOffset = 5;
lw = 3;

pts_null = eng_neuron_impacts_null(:);
pts = eng_neuron_impacts(:);

plot.init(28);
ymn = 0; ymx = 1;
bins = linspace(0,100,31);
yscs = {100*pts_null, 100*pts};
clrs = {0.7*ones(1,3), 'k'};

[p,h] = ranksum(yscs{1}, yscs{2});
p

mus = cell(numel(yscs),1);
for ll = 1:numel(yscs)
    cs = histcounts(yscs{ll}, bins);
    cs = 100*cs/sum(cs);
    plot.histogramPretty(cs, bins, clrs{ll}, lw);
    mus{ll} = nanmedian(yscs{ll});
end

axis tight;
xlim([min(bins) max(bins)]);
yl = ylim;
for ll = 1:numel(mus)
    plot(mus{ll}, yl(2)+yOffset, 'v', 'MarkerSize', 10, ...
        'LineWidth', 1, 'Color', clrs{ll}, ...
        'MarkerFaceColor', clrs{ll});
end

set(gca, 'LineWidth', lw);
set(gca, 'TickDir', 'out');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'YTick', 0:20:max(ylim));

xlabel({'% of units with same neural', 'engagement coefficient sign'});
ylabel('% of targets');
plot.setPrintSize(gcf, struct('width', 5, 'height', 4));
if doSave
    export_fig(gcf, fullfile(saveDir, [fnm '.pdf']));
end
