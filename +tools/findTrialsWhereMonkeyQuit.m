function [trialQuitRanges, trialQuitMasks, trialQuitCounts] = ...
    findTrialsWhereMonkeyQuit(dts, nIncorrectInARow)
    if nargin < 2
        nIncorrectInARow = 5;
    end
    grps = tools.thetaCenters;
    
    trialQuitCounts = nan(numel(dts), 2);
    trialQuitRanges = cell(numel(dts), 2);
    trialQuitMasks = cell(numel(dts), numel(grps), 2);    

    for ii = 1:numel(dts)
        Bs = io.loadPreprocessedData(dts{ii});
        for jj = 1:2
            B = Bs(jj);
            iscors = grpstats(B.isCorrect, B.trial_index);
            trgs = grpstats(B.target, B.trial_index);
            trs = grpstats(B.trial_index, B.trial_index);
            if jj == 1
                trs = trs - max(trs);
            elseif jj == 2
                trs = trs - min(trs) + 1;
            end
            
            % get trials where we have nIncorrectInARow incorrects in a row
            cs = tools.smooth((1:numel(iscors))', iscors, ...
                nIncorrectInARow, 1, @(ys) ~any(ys));
            trialQuitCounts(ii,jj) = sum(diff(cs) == 1);

            % find trial ranges of each quit window
            trialStarts = find(diff(cs) == 1)+1;
            trialEnds = nan(numel(trialStarts),1);
            for ll = 1:numel(trialStarts)
                nextCor = find(iscors((trialStarts(ll)+1):end) == true, 1, 'first');
                if isempty(nextCor)
                    trialEnds(ll) = numel(iscors);
                else
                    trialEnds(ll) = trialStarts(ll) + nextCor-1;
                end
            end
            trialRange = [trialStarts trialEnds];
            trialQuitRanges{ii,jj} = trialRange;

            % make mask per target
            for kk = 1:numel(grps)
                ctrs = trs(trgs == grps(kk));
                mask = false(size(ctrs));
                for ll = 1:size(trialRange,1)
                    mask(ismember(ctrs, trialRange(ll,1):trialRange(ll,2))) = true;
                end
                trialQuitMasks{ii,kk,jj} = mask;
            end
        end
    end
    
    % count ignored trials relative to total number of trials
    cs = arrayfun(@(ii) sum(cellfun(@sum, trialQuitMasks(ii,:,2))), ...
        1:size(trialQuitMasks,1));
    ns = arrayfun(@(ii) sum(cellfun(@numel, trialQuitMasks(ii,:,2))), ...
        1:size(trialQuitMasks,1));
    pctIgnored = 100*sum(cs)/sum(ns);
    disp(['Monkey quit on ' sprintf('%0.2f', pctIgnored) ...
        '% of all trials during Block 2.']);
end
