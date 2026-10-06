function F = SHE_Solver(alpha, I_dc, I_1_target)
%SHE_SOLVER  v1 residual function, kept for backward compatibility.
%
%   F = SHE_Solver(alpha, I_dc, I_1_target)
%
%   6 switching angles; fundamental set to I_1_target (absolute units) and
%   harmonics 5, 7, 11, 13, 17 eliminated, for a DC link of I_dc. Identical
%   output to v1 (checked in tests/test_v1_compat.m).
%
%   v1 note: the published v1 README called this function by its internal
%   name SHE_nonlinear_system_with_jacobian, which MATLAB cannot find (it
%   uses the file name). Call SHE_Solver, or better use the v2 functions:
%       alpha = she_solve(I_1_target / I_dc, [5 7 11 13 17], alpha0);
%
%   See also SHE_SOLVE, SHE_RESIDUAL.

    F = I_dc * she_residual(alpha, I_1_target / I_dc, [5 7 11 13 17]);
end
