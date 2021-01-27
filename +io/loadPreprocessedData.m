function [Bs,D] = loadPreprocessedData(dtstr, includeFreezePeriod)
    if nargin < 2
        includeFreezePeriod = false;
    end
    d = load(fullfile('data', 'preprocessed', [dtstr '.mat']));
    Bs = d.Bs;
    D = d.D;
    if ~includeFreezePeriod
        for jj = 1:numel(Bs)
            Bs(jj) = io.filterTrialsByIdx(Bs(jj), Bs(jj).time >= 7);
        end
    end
    
    vfcn_int = tools.makeVelFcn(D.decoders(1).fDecoder, false);
    vfcn_wmp = tools.makeVelFcn(D.decoders(2).fDecoder, false);
    D.vfcns = {vfcn_int, vfcn_wmp};
end
