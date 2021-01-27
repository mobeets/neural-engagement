%% This code generates Supplemental Fig. 10
% (using data from only two example sessions)

%% compare angle between vigor dims estimated during Int vs. end of WMP

dts = io.getDates;
grps = tools.thetaCenters;
first_prin_angles = nan(numel(dts), numel(grps));

for ii = 1:numel(dts)
    dtstr = dts{ii};
    eng_info_pre = engdims.getAimingEllipseAndEngagementDimensions(dtstr, 'pre');
    eng_info_post = engdims.getAimingEllipseAndEngagementDimensions(dtstr, 'post');
    
    for kk = 1:numel(grps)
        ytrg = eng_info_post.Ysmu_anchors(kk,:);
        yrefs = eng_info_pre.Ysmu;
        [~,ixd] = min(pdist2(ytrg, yrefs));
        
        vtrg = eng_info_post.engagement_dims_anchors(kk,:);
        vref = eng_info_pre.engagement_dims(ixd,:);
        first_prin_angles(ii,kk) = rad2deg(subspace(vtrg', vref'));
    end
end

%% find average principal angle between two random vectors

% find principal angles between random n-d vectors
K = 10; % # of dimensions
nboots = 50000;
rndvecs = randn(2,nboots,K);
rndangs = nan(nboots,1);
for ll = 1:nboots
    a1 = zeros(1,K); a1(1) = 1;
    a2 = squeeze(rndvecs(2,ll,:))';
    a2 = a2/norm(a2);
    rndangs(ll,:) = rad2deg(tools.prinangle(a1', a2'));
end
chanceAngs = tools.computeAngles([cosd(rndangs) ...
    sind(rndangs)], true);
chanceAng = tools.computeAngles(nanmean([cosd(rndangs) ...
    sind(rndangs)],1),true);

%% plot first principal angle of dims during Int vs. end of WMP

doSave = false;
fnm = 'SuppFig10';
saveDir = 'data/figures';
lw = 3;
dataClr = 'k';
grayClr = 0.7*ones(1,3);

plot.init(28);
bins = linspace(0,90,21);
for ll = 1:2
    if ll == 2
        ysc = first_prin_angles(:);
        clr = dataClr;
    else
        ysc = chanceAngs;
        clr = grayClr;
    end
    cs = histcounts(ysc,bins);
    cs = 100*cs/sum(cs);

    plot.histogramPretty(cs, bins, clr, lw);
end
ys = first_prin_angles(:);
mu = tools.computeAngles(nanmedian([cosd(ys) sind(ys)]), true);
sprintf('Change in %% variance (mean): %0.1f', mu);

[p,h] = ranksum(ys, chanceAngs);
p

axis tight;
xlim([0 90]);

yl = ylim;
yoffset = 2;
plot(mu, yl(2)+yoffset, 'v', 'MarkerSize', 10, ...
    'LineWidth', 1, 'Color', dataClr, 'MarkerFaceColor', dataClr);

muChance = tools.computeAngles(nanmedian(...
    [cosd(chanceAngs) sind(chanceAngs)]), true);
plot(muChance, yl(2)+yoffset, 'v', 'MarkerSize', 10, ...
    'LineWidth', 1, 'Color', grayClr, 'MarkerFaceColor', grayClr);

xlabel({'angle (degrees) between', 'neural engagement axes', ...
    'before/after learning'});
ylabel('% of targets');
set(gca, 'LineWidth', lw);
set(gca, 'XTick', [0 45 90]);
set(gca, 'TickDir', 'out');
set(gca, 'XColor', 'k', 'YColor', 'k');
plot.setPrintSize(gcf, struct('width', 5, 'height', 4.5));
if doSave
    export_fig(gcf, fullfile(saveDir, [fnm '.pdf']));
end
