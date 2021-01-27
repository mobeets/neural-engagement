function nrms = rowwiseNorm(Xs)
    nrms = sqrt(sum(Xs.^2,2));
end
