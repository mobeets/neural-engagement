function secs = trialTimesToSeconds(tms)
    if isa(tms, 'double')
        if all(tms < 24*60*60)
            % already been converted to seconds, so do nothing
            secs = tms;
            return;
        end
        tms = num2str(tms);
    end
    hhmmss = tms(:,end-5:end);
    hh = str2num(hhmmss(:,1:2));
    mm = str2num(hhmmss(:,3:4));
    ss = str2num(hhmmss(:,5:6));
    secs = 60*60*hh + 60*mm + ss;
end
