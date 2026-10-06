function J = she_jacobian(alpha, harmonics)
%SHE_JACOBIAN  Analytic Jacobian of SHE_RESIDUAL with respect to alpha.
%
%   J = she_jacobian(alpha, harmonics)
%
%   d b_n / d alpha_k = -(4/pi) * (-1)^(k-1) * sin(n * alpha_k)
%
%   Rows follow [1, harmonics], columns follow alpha.

    alpha = alpha(:).';
    n = [1; harmonics(:)];
    sgn = (-1) .^ (0:numel(alpha) - 1);
    J = -(4 / pi) * sin(n * alpha) .* sgn;    % implicit expansion over rows
end
