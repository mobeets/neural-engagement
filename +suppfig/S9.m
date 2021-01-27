%% This code generates Supplemental Fig. 9C-D
% (using data from only two example sessions)

%% plot progress traces (smoothed, and normalized to max) for T1 and T2
% required inputs: ProgWmp, vigorIncreaseIsGoodForProgress

pctThreshes = 0.75:0.05:1.0;
kSmooths = 4:2:16;
N2 = 75;
Ysc = ProgWmp(:,:,2);

% take first N trials per target, smooth, and find where each crosses threshold
yss = tools.cellToFixedLengthMat(Ysc, N2);
LearningSpeeds = nan(size(yss,1), numel(pctThreshes), numel(kSmooths));
for kk = 1:numel(kSmooths)
    % smooth all targets by the same amount
    ysc = tools.smooth((1:N2), yss, kSmooths(kk), 1, @mean);
    assert(size(ysc,1) == size(yss,1));
    % scale each target so that performance is relative to max
    ysc = ysc - repmat(nanmin(ysc, [], 2), 1, size(ysc,2));
    ysc_scaled = ysc./repmat(nanmax(ysc, [], 2), 1, size(ysc,2));

    for ll = 1:numel(pctThreshes)
        for ii = 1:size(ysc_scaled,1)
            if all(isnan(ysc_scaled(ii,:)))
                continue;
            end
            % find trial where we first cross some percentage of max
            LearningSpeeds(ii,ll,kk) = find(ysc_scaled(ii,:) >= pctThreshes(ll)*(1-eps),1);
        end
    end
end

%% plot how median of histograms depend on each parameter choice

doSave = false;
saveDir = 'data/figures';
fnm = 'SuppFig9';
lw = 2;
FontSize = 14;
targetType = engagementIncreaseIsGoodForProgress(:);
clrs = {[0 161 228]/255, [237 28 36]/255};
pctThreshToShow = 1.0;
kSmoothToShow = 8;

plot.init(FontSize);
paramNames = {'# trials smoothing', 'threshold (% of max)'};
mnks = io.monkeyNames; mnks = mnks([1 3 2]);
nrows = numel(paramNames); ncols = numel(mnks); c = 1;

for kk = 1:numel(paramNames)
    for mm = 1:numel(mnks)
        ixm = repmat(io.getMonkeyDateFilter(dts, mnks(mm))', 1, numel(grps));
        ixm = ixm(:);
        plot.subplot(nrows,ncols,c,FontSize); c = c+1;
        
        for jj = 1:2
            if kk == 1
                xs = kSmooths;
                ys = squeeze(LearningSpeeds(ixm & targetType == -jj+2, ...
                    pctThreshes == pctThreshToShow, :));
            else
                xs = 100*pctThreshes;
                ys = squeeze(LearningSpeeds(ixm & targetType == -jj+2, ...
                    :, kSmooths == kSmoothToShow));
            end            
            plot(xs, nanmedian(ys,1), '.-', 'LineWidth', lw, 'Color', clrs{jj});
        end
        set(gca, 'LineWidth', lw);
        axis tight;
        xlabel(paramNames{kk});
        ylabel({'# trials to', 'peak performance'});
        title(['Monkey ' mnks{mm}(1)]);
    end
end

plot.setPrintSize(gcf, struct('width', 8, 'height', 4.5));
if doSave
    export_fig(gcf, fullfile(saveDir, [fnm '.pdf']));
end
