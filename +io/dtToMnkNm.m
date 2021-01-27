function [mnkNms, mnkIds] = dtToMnkNm(dts, splitJeffy)
    if nargin < 2
        splitJeffy = true;
    end
    if ~iscell(dts)
        assert(ischar(dts), 'must provide cell or char');
        dts = {dts};
        providedCell = false;
    else
        providedCell = true;
    end
    mnkNms = cell(numel(dts),1);
    mnkIds = nan(numel(dts),1);
    for ii = 1:numel(dts)
        dtstr = dts{ii};
        [mnkNms{ii}, mnkIds(ii)] = dtToMnkNm_inner(dtstr, splitJeffy);
    end
    if ~providedCell
        assert(numel(dts) == 1, 'internal error');
        mnkNms = mnkNms{1};
        mnkIds = mnkIds(1);
    end
end

function [mnkNm, mnkId] = dtToMnkNm_inner(dtstr, splitJeffy)
    if nargin < 2
        splitJeffy = false;
    end
    if strcmpi(dtstr(1:4), '2012')
        mnkNm = 'Jeffy';
        mnkId = 1;
        mm = str2double(dtstr(6));
        if splitJeffy && mm >= 5
            mnkNm = 'Jeffy-late';
            mnkId = 2;
        elseif splitJeffy
            mnkNm = 'Jeffy-early';
            mnkId = 1;
        end
    elseif strcmpi(dtstr(1:4), '2013')
        mnkNm = 'Lincoln';
        mnkId = 2;
        if splitJeffy
            mnkId = 3;
        end
    else
        mnkNm = 'Nelson';
        mnkId = 3;
        if splitJeffy
            mnkId = 4;
        end
    end
end
