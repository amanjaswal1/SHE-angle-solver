function X = halton_points(idx, d)
%HALTON_POINTS  Deterministic, evenly spread points in the unit cube (0,1)^d.
%
%   X = halton_points(idx, d)
%
%   idx  vector of point indices (1, 2, 3, ...);  d  dimension (<= 15).
%   X    numel(idx) x d. Row k uses the radical inverse of idx(k) in the
%        first d prime bases (the Halton sequence).
%
%   Used instead of RAND for multi-start searches, so MATLAB and Octave -
%   whose random-number generators differ - explore exactly the same
%   starting points and find exactly the same solutions.

    primes_list = [2 3 5 7 11 13 17 19 23 29 31 37 41 43 47];
    if d > numel(primes_list)
        error('halton_points:d', 'At most %d dimensions.', numel(primes_list));
    end
    idx = idx(:);
    X = zeros(numel(idx), d);
    for j = 1:d
        b = primes_list(j);
        for k = 1:numel(idx)
            n = idx(k); f = 1; r = 0;
            while n > 0
                f = f / b;
                r = r + f * mod(n, b);
                n = floor(n / b);
            end
            X(k, j) = r;
        end
    end
end
