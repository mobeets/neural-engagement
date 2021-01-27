%% This code generates Supplemental Fig. 2
% (using data from only one example BCI session)

%% load avg hand speeds per trial from all sessions

dtsBci = io.getDates;
SpdsBciPerTrial = cell(numel(dtsBci), 2);
SpdsBciWithinTrial = cell(numel(dtsBci), 2);
for ll = 1:numel(dtsBci)
    dtstr = dtsBci{ll};
    [Bs,D] = io.loadPreprocessedData(dtstr);
    for jj = 1:2
        B = Bs(jj);
        
        ixt = B.time >= 7 & B.time <= 20;
        SpdsBciPerTrial{ll,jj} = grpstats(B.handSpeed(ixt), B.trial_index(ixt));
        
        ixt = B.time <= 40;
        SpdsBciWithinTrial{ll,jj} = grpstats(B.handSpeed(ixt), B.time(ixt));
    end
end
mnkNmsBci = io.monkeyNames;
mnkNmsBci = mnkNmsBci(1:2);

%% plot hand speeds within a trial

saveDir = 'data/figures';
doSave = false;
fnm = 'SuppFig2A';
isBciData = true;
clr = 0.2*ones(1,3);

for mm = 1:numel(mnkNmsBci)    
    ixm = io.getMonkeyDateFilter(dtsBci, mnkNmsBci(mm));
    if sum(ixm) == 0
        continue;
    end
    spds = SpdsBciWithinTrial(ixm,2);
    minTime = 7*45;
    maxTime = 1000;
    
    plot.init(20);
    for jj = 1:size(spds,2)
        ys = spds(:,jj);
        N = floor(prctile(cellfun(@numel, ys), 5));
        ys = ys(cellfun(@numel, ys) >= N);
        ys = cell2mat(cellfun(@(y) y(1:N), ys, 'uni', 0)')';
        
        mus = prctile(ys, 50, 1);
        lbs = prctile(ys, 25, 1);
        ubs = prctile(ys, 75, 1);
        
        xs = 1:N;
        xs = xs*45;
        ixt = xs >= minTime;
        h = plot(xs(ixt) - minTime, mus(ixt), '-', 'Color', clr, 'LineWidth', 2);
        plot.plotPolygonFromBounds(xs(ixt) - minTime, ...
            lbs(ixt), ubs(ixt), plot.dullColors(h.Color));
    end
    xlim([0 maxTime]);
    ylim([0 0.5]);
    xlabel('time rel. to cursor release (ms)');
    ylabel('hand speed (m/s)');
    title(['Monkey ' mnkNmsBci{mm}(1) ' (BCI)']);
    set(gca, 'YTick', [0 0.25 0.5]);
    set(gca, 'LineWidth', 2);
    plot.setPrintSize(gcf, struct('width', 3, 'height', 2.5));
    if doSave
        export_fig(gcf, fullfile(saveDir, [fnm '_' mnkNmsBci{mm} '.pdf']));
    end
end

%% plot hand speeds over trials

doSave = false;
saveDir = 'data/figures';
fnm = 'SuppFig2D';
clr = 0.2*ones(1,3);

for mm = 1:numel(mnkNmsBci)
    ixm = io.getMonkeyDateFilter(dtsBci, mnkNmsBci(mm));
    if sum(ixm) == 0
        continue;
    end
    spds = SpdsBciPerTrial(ixm,:);
    plot.init(20);
    for jj = 1:size(spds,2)
        ys = spds(:,jj);
        N = floor(prctile(cellfun(@numel, ys), 5));
        ys = ys(cellfun(@numel, ys) >= N);
        ys = cell2mat(cellfun(@(y) y(1:N), ys, 'uni', 0)')';
        mus = prctile(ys, 50, 1);
        lbs = prctile(ys, 25, 1);
        ubs = prctile(ys, 75, 1);
        
        if jj == 1
            xs = -N:-1;
        elseif jj == 2
            xs = 1:N;
        end
        h = plot(xs, mus, '-', 'Color', clr, 'LineWidth', 2);
        plot.plotPolygonFromBounds(xs, lbs, ubs, ...
            plot.dullColors(h.Color, 0.8));
    end
    axis tight;    
    ylabel('hand speed (m/s)');
    xlabel('# trials, rel. to Block 2');
    xlim([-50 50]);
    title(['Monkey ' mnkNmsBci{mm}(1) ' (BCI)']);
    
    plot([0 0], ylim, 'k--', 'LineWidth', 2);
    set(gca, 'LineWidth', 2);
    if mm == 1
        set(gca, 'YTick', [0 0.0015 0.003]);
    else
        set(gca, 'YTick', [0 0.025 0.05]);
    end
    plot.setPrintSize(gcf, struct('width', 4.5, 'height', 2.25));
    if doSave
        export_fig(gcf, fullfile(saveDir, [fnm '_' mnkNmsBci{mm} '.pdf']));
    end
end
