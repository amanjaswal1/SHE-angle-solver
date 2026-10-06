function test_newton_quadratic()
% Near the solution, Newton with the exact Jacobian converges quadratically:
% the observed order log(e_k+1/e_k)/log(e_k/e_k-1) should approach 2.
    H = [5 7 11 13 17];
    [~, info] = she_solve(0.6, H, [10 20 30 40 50 60] * pi / 180);
    e = info.history;
    e = e(e > 1e-14);                              % drop values at round-off
    p = log(e(3:end) ./ e(2:end-1)) ./ log(e(2:end-1) ./ e(1:end-2));
    assert(max(p) > 1.8, sprintf('best observed order %.2f, expected ~2', max(p)));
end
