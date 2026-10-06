function test_halton()
% Halton sequence: radical inverses in bases 2 and 3 (known values), all
% points strictly inside (0, 1).
    X = halton_points(1:4, 2);
    expect = [1/2 1/3; 1/4 2/3; 3/4 1/9; 1/8 4/9];
    assert(max(abs(X(:) - expect(:))) < 1e-15, 'wrong Halton values');
    Y = halton_points(1:500, 6);
    assert(all(Y(:) > 0 & Y(:) < 1), 'points must lie inside the unit cube');
end
