function test_sweep_continuity()
% Every traced branch must be feasible at every point and change smoothly
% with m (no jumps between branches).
    rng(11);
    br = she_sweep(0.30:0.02:0.70, [5 7], 20, 5);
    assert(~isempty(br), 'no branches found');
    for i = 1:numel(br)
        a = br{i}.alpha;
        assert(all(a(:, 1) > 0) && all(a(:, end) < pi/2), 'angles outside (0, pi/2)');
        assert(all(all(diff(a, 1, 2) > 0)), 'angles not increasing');
        if size(a, 1) > 1
            assert(max(max(abs(diff(a, 1, 1)))) <= 0.15, 'branch jumped');
        end
        for k = 1:numel(br{i}.m)
            assert(norm(she_residual(a(k, :), br{i}.m(k), [5 7])) < 1e-10, 'residual too large');
        end
    end
end
