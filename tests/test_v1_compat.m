function test_v1_compat()
% SHE_Solver (v2 wrapper) must reproduce the published v1 residual exactly.
    rng(7);
    for trial = 1:20
        alpha = sort(rand(1, 6)) * pi / 2;
        I_dc = 50 + 300 * rand;
        I1 = I_dc * rand;
        d = SHE_Solver(alpha, I_dc, I1) - she_v1_reference(alpha, I_dc, I1);
        assert(max(abs(d)) < 1e-9 * I_dc, sprintf('trial %d: max diff %.3g', trial, max(abs(d))));
    end
end
