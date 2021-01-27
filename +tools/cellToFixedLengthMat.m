function yss = cellToFixedLengthMat(Ys, N, takeLastN)
% convert ys, a Tx1 cell array, each of which is a variable length array,
% into a TxN matrix (with nans-filled where necessary)

    if nargin < 3
        takeLastN = false;
    end
    yss = nan(numel(Ys), N);
    for ii = 1:numel(Ys)
        ysc = Ys{ii};
        Nc = min([N numel(ysc)]);
        if takeLastN
            yss(ii,1:Nc) = ysc(end-Nc+1:end);
        else
            yss(ii,1:Nc) = ysc(1:Nc);
        end
    end
end