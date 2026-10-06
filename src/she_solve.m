function [alpha, info] = she_solve(m, harmonics, alpha0, opts)
%SHE_SOLVE  Solve the SHE equations with a damped Newton method.
%
%   [alpha, info] = she_solve(m, harmonics, alpha0)
%   [alpha, info] = she_solve(m, harmonics, alpha0, opts)
%
%   m          Target fundamental, per unit of the DC level (0 < m < 4/pi).
%   harmonics  Odd harmonic orders to eliminate, e.g. [5 7 11 13 17].
%   alpha0     Initial guess, N = numel(harmonics) + 1 angles [rad].
%   opts       Optional struct: tol (1e-12), max_iter (50).
%
%   alpha      Solution [rad], each angle mapped into [0, pi].
%   info       Struct: converged, feasible, iterations, residual_norm,
%              history (residual norm at every iteration).
%
%   Method: Newton steps d = -J \ F with the analytic Jacobian, a
%   backtracking line search on ||F|| for robustness far from the solution,
%   and a Levenberg-Marquardt step when J is close to singular. Converges
%   quadratically near a solution (see results/figures). No toolbox needed;
%   MATLAB's fsolve is NOT required.
%
%   feasible = converged AND 0 < alpha_1 < ... < alpha_N < pi/2.

    if nargin < 4, opts = struct(); end
    tol      = get_opt(opts, 'tol', 1e-12);
    max_iter = get_opt(opts, 'max_iter', 50);

    alpha = alpha0(:).';
    F = she_residual(alpha, m, harmonics);
    nF = norm(F);
    history = nF;
    it = 0;
    while nF > tol && it < max_iter
        it = it + 1;
        J = she_jacobian(alpha, harmonics);
        if rcond(J) > 1e-12
            d = -(J \ F).';
        else                                   % Levenberg-Marquardt fallback
            mu = 1e-6 * max(1, norm(J, 'fro'))^2;
            d = -((J.' * J + mu * eye(numel(alpha))) \ (J.' * F)).';
        end
        % Backtracking line search (Armijo condition on ||F||)
        step = 1;
        while step > 1e-10
            trial = alpha + step * d;
            F_trial = she_residual(trial, m, harmonics);
            if norm(F_trial) <= (1 - 1e-4 * step) * nF
                break
            end
            step = step / 2;
        end
        alpha = trial;
        F = F_trial;
        nF = norm(F);
        history(end + 1) = nF; %#ok<AGROW>
        if step <= 1e-10
            break                              % stalled: no descent direction
        end
    end

    % cos(n*alpha) is unchanged by alpha -> alpha + 2*pi and alpha -> -alpha,
    % so map every angle into [0, pi] without changing the residual.
    alpha = mod(alpha, 2*pi);
    alpha(alpha > pi) = 2*pi - alpha(alpha > pi);
    nF = norm(she_residual(alpha, m, harmonics));
    converged = nF <= max(tol, 1e-10);
    % A solution of the equations is only a usable switching pattern if the
    % angles are strictly increasing inside the first quarter period.
    feasible = converged && alpha(1) > 0 && alpha(end) < pi/2 && all(diff(alpha) > 1e-9);

    info = struct('converged', converged, 'feasible', feasible, ...
                  'iterations', it, 'residual_norm', nF, 'history', history);
end

function v = get_opt(opts, name, default)
    if isfield(opts, name), v = opts.(name); else, v = default; end
end
