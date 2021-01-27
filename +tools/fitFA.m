function params = fitFA(Y, ndims, doZScore, varargin)
    if nargin < 2
        ndims = 1:20;
    end
    if nargin < 3
        doZScore = true;
    end
    ix = ~any(isnan(Y),2);
    Y = Y(ix,:);
    if doZScore
        mu = mean(Y,1);
        sdev = std(Y,[],1);
        Y = (Y - repmat(mu, size(Y,1), 1))./repmat(sdev, size(Y,1), 1);
    end
    if ~exist('crossvalidate_fa', 'file')
        error(['Missing FA source code. ' ...
            'Download from: https://github.com/mobeets/fa']);
    end
    dim = crossvalidate_fa(Y', 'zDimList', ndims, varargin{:});
    params = {dim.estParams};
    sumPE = [dim.sumPE];
    sumLL = [dim.sumLL];
    for ii = 1:numel(params)
        params{ii}.FactorAnalysisParams.L = params{ii}.L;
        params{ii}.FactorAnalysisParams.Ph = params{ii}.Ph;
        params{ii}.FactorAnalysisParams.d = params{ii}.d;
        if doZScore
            params{ii}.NormalizeSpikes.mean = mu;
            params{ii}.NormalizeSpikes.std = sdev;
        end
        params{ii}.sumPE = sumPE(ii);
        params{ii}.sumLL = sumLL(ii);
        params{ii} = rmfield(params{ii}, 'L');
        params{ii} = rmfield(params{ii}, 'Ph');
        params{ii} = rmfield(params{ii}, 'd');
    end
    
end
