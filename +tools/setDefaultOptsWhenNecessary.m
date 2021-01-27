function opts = setDefaultOptsWhenNecessary(opts, defopts, warnIfNoDefault)
    if nargin < 3
        % if opts has a field that defopts does not, give a warning
        warnIfNoDefault = false;
    end
    fns = fieldnames(defopts);
    for ii = 1:numel(fns)
        if ~isfield(opts, fns{ii})
            opts.(fns{ii}) = defopts.(fns{ii});
        end
    end
    if ~warnIfNoDefault
        return;
    end
    extra_fns = setdiff(fieldnames(opts), fns);
    if numel(extra_fns) > 0
        warning(['The following fields in opts had no defaults: ''' ...
            strjoin(extra_fns, ''', ''') '''. Are these fields valid?']);
    end
end
