function vfcn = makeVelFcn(dec, useSteadyState)
    
    if useSteadyState
        % i.e., suppose v_{t-1} = v_t
        decfcn = @(d, y) bsxfun(@plus, ...
            ((eye(size(d.M1)) - d.M1)\d.M2)*y', ...
            (eye(size(d.M1)) - d.M1)\d.M0)';
    else
        % i.e., v_t = B*z_t + c
        decfcn = @(d, y) (d.M2*y' + repmat(d.M0, 1, size(y,1)))';
    end
    vfcn = @(y) decfcn(dec, y);
end
