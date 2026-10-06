function test_jacobian_fd()
% Analytic Jacobian must match central finite differences.
    H = [5 7 11 13 17];
    alpha = [5 12 21 33 47 61] * pi / 180;
    J = she_jacobian(alpha, H);
    h = 1e-6;
    J_fd = zeros(size(J));
    for k = 1:numel(alpha)
        e = zeros(size(alpha)); e(k) = h;
        J_fd(:, k) = (she_residual(alpha + e, 0.5, H) - she_residual(alpha - e, 0.5, H)) / (2*h);
    end
    rel = max(abs(J(:) - J_fd(:))) / max(abs(J(:)));
    assert(rel < 1e-8, sprintf('Jacobian relative error %.3g', rel));
end
