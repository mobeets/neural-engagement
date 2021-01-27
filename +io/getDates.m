function dts = getDates()
% function dts = getDates()
% 
    files = dir(fullfile('data', 'preprocessed', '201*.mat'));
    dts = cellfun(@(d) d(1:end-4), {files.name}, 'uni', 0);
end
