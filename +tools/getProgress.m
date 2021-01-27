function [prog, angErr, progDist] = getProgress(y, curPos, trgPos, vfcn, v)
    if nargin < 4
        vfcn = [];
    end
    if nargin < 5
        v = vfcn(y);
    end
    if isempty(curPos)
        curPos = zeros(size(trgPos));
    end
    
    goal = trgPos - curPos;
    nrm = sqrt(sum(goal.*goal,2));
    goal = bsxfun(@times, goal, 1./nrm);
    [nt1, nd1] = size(v);
    [nt2, nd2] = size(goal);
    assert((nt1 == nt2) & (nd1 == nd2));
    prog = sum(v.*goal, 2);
    if nargout > 1
        % compute angular error (in degrees)
        angErr = tools.getAngularError(v, goal);
    end
    if nargout > 2
        % compute progress in terms of decreased distance to target
        origDists = sqrt(sum((trgPos - curPos).^2,2));
        deltaPos = (45/1000)*v;
        newDists = sqrt(sum((trgPos - (curPos + deltaPos)).^2,2));
        progDist = origDists - newDists;
    end
end
