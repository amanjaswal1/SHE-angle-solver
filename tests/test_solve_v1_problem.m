function test_solve_v1_problem()
% The v1 problem (6 angles, eliminate 5,7,11,13,17) at m = 0.6 from a fixed
% guess: residual at machine precision and a feasible, ordered pattern.
    H = [5 7 11 13 17];
    alpha0 = [10 20 30 40 50 60] * pi / 180;      % v1 README initial guess
    [alpha, info] = she_solve(0.6, H, alpha0);
    assert(info.converged, 'did not converge');
    assert(info.feasible, 'solution not ordered inside (0, pi/2)');
    b = she_harmonics(alpha, [1 H]);
    assert(abs(b(1) - 0.6) < 1e-12, 'fundamental not at target');
    assert(max(abs(b(2:end))) < 1e-12, 'harmonics not eliminated');
end
