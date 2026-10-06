function branches = she_sweep(m_grid, harmonics, n_starts, seed_every)
%SHE_SWEEP  Map every SHE solution branch over a range of modulation index.
%
%   branches = she_sweep(m_grid, harmonics)
%   branches = she_sweep(m_grid, harmonics, n_starts, seed_every)
%
%   m_grid      Increasing vector of fundamental targets (per unit).
%   harmonics   Harmonic orders to eliminate.
%   n_starts    Initial guesses per seed point (default 60).
%   seed_every  Run the multi-start search at every k-th grid point
%               (default 5).
%
%   branches    Cell array of structs with fields
%                 m      1 x K   modulation index along the branch
%                 alpha  K x N   feasible angles [rad] at each m
%
%   Two-stage method:
%     1. Multi-start: sorted guesses in (0, pi/2) from a Halton sequence
%        (deterministic, so MATLAB and Octave find the same solutions),
%        solved with SHE_SOLVE, discover solutions at the seed points.
%     2. Continuation: each new solution is traced forwards and backwards
%        along m_grid, using a secant predictor (linear extrapolation of the
%        last two points) as the Newton starting guess. A branch ends where
%        Newton fails, the angles leave the feasible region, or the step
%        would jump to a different branch.
%   Multi-start alone rarely lands on solutions near the edge of the
%   feasible range; continuation follows the branch there smoothly.
%   Multi-start can never prove that every branch was found.

    if nargin < 3 || isempty(n_starts),   n_starts = 60;  end
    if nargin < 4 || isempty(seed_every), seed_every = 5; end
    N = numel(harmonics) + 1;
    branches = {};

    n_seed = 0;
    for j = 1:seed_every:numel(m_grid)
        m = m_grid(j);
        starts = sort(halton_points(n_seed * n_starts + (1:n_starts), N), 2) * pi / 2;
        n_seed = n_seed + 1;
        for s = 1:n_starts
            [a, info] = she_solve(m, harmonics, starts(s, :));
            if ~info.feasible || on_known_branch(branches, m, a)
                continue
            end
            [m_f, a_f] = trace(m_grid, j, a, +1, harmonics);
            [m_b, a_b] = trace(m_grid, j, a, -1, harmonics);
            br.m     = [fliplr(m_b), m, m_f];
            br.alpha = [flipud(a_b); a; a_f];
            branches{end + 1} = br; %#ok<AGROW>
        end
    end
end

function [m_out, a_out] = trace(m_grid, j, a, direction, harmonics)
% Natural-parameter continuation from grid index j in one direction.
    m_out = [];
    a_out = zeros(0, numel(a));
    prev = a;
    prev2 = [];
    k = j + direction;
    while k >= 1 && k <= numel(m_grid)
        if isempty(prev2)
            guess = prev;
        else
            guess = 2 * prev - prev2;          % secant predictor
        end
        [a_new, info] = she_solve(m_grid(k), harmonics, guess);
        if ~info.feasible || max(abs(a_new - prev)) > 0.15
            break                              % failed, or jumped branches
        end
        m_out(end + 1) = m_grid(k); %#ok<AGROW>
        a_out(end + 1, :) = a_new; %#ok<AGROW>
        prev2 = prev;
        prev = a_new;
        k = k + direction;
    end
end

function tf = on_known_branch(branches, m, a)
% True if solution a at modulation index m already lies on a traced branch.
    tf = false;
    for i = 1:numel(branches)
        idx = find(abs(branches{i}.m - m) < 1e-12, 1);
        if ~isempty(idx) && max(abs(branches{i}.alpha(idx, :) - a)) < 1e-6
            tf = true;
            return
        end
    end
end
