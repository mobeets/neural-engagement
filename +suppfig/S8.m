%% This code generates Supplemental Fig. 8
% (using data from only two example sessions)
% data required: ProgInt, ProgWmp (from figure.F5and6)

% ProgCur = ProgInt(:,:,1); fnm = 'SuppFig8A';
ProgCur = ProgWmp(:,:,1); fnm = 'SuppFig8B';

baselineProgress = cellfun(@nanmean, ProgCur);
targetType = nan(size(baselineProgress));
mnks = io.monkeyNames;
for mm = 1:numel(mnks)
    ixm = io.getMonkeyDateFilter(dts, mnks(mm));
    cpts = baselineProgress(ixm,:);
    targetType(ixm,:) = cpts > nanmedian(cpts(:));
end
targetType = targetType(:);
clrs = {[0.2 0.5 0.2], 0.6*ones(1,3)};

% Now, run the cell block in figure.F5and6 corresponding to Fig. 5E
%   after commenting out the lines defining targetType and clrs
