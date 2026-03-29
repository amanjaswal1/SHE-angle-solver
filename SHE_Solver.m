% Description: Nonlinear system of equations for Selective Harmonic
%              Elimination (SHE) angle computation in a 6-pulse CSI.
%              Solves for 6 switching angles that eliminate harmonics
%              5, 7, 11, 13, 17 and control the fundamental.

function F = SHE_nonlinear_system_with_jacobian(alpha, I_dc, I_1_target)

% --- Compute harmonic amplitudes ---
a_1  = compute_a_n(1,  alpha, I_dc);
a_5  = compute_a_n(5,  alpha, I_dc);
a_7  = compute_a_n(7,  alpha, I_dc);
a_11 = compute_a_n(11, alpha, I_dc);
a_13 = compute_a_n(13, alpha, I_dc);
a_17 = compute_a_n(17, alpha, I_dc);

% --- SHE constraint equations (F = 0 at solution) ---
% 6 equations for 6 unknowns (alpha_1 ... alpha_6)
F = zeros(6, 1);
F(1) = a_1  - I_1_target;   % Fundamental = target
F(2) = a_5;                  % Eliminate 5th harmonic
F(3) = a_7;                  % Eliminate 7th harmonic
F(4) = a_11;                 % Eliminate 11th harmonic
F(5) = a_13;                 % Eliminate 13th harmonic
F(6) = a_17;                 % Eliminate 17th harmonic

end

% ---------------------------------------------------------------
% Helper: Fourier coefficient for nth harmonic (half-wave symmetry)
% ---------------------------------------------------------------
function a_n = compute_a_n(n, alpha, I_dc)

% Standard CSI Fourier series (symmetric 6-pulse pattern):
% a_n = (4*I_dc / (pi*n)) * sum_{k=1}^{6} (-1)^(k-1) * cos(n*alpha_k)

term = 0;
for k = 1:length(alpha)
    term = term + (-1)^(k-1) * cos(n * alpha(k));
end

a_n = (4 * I_dc / (pi * n)) * term;

end
