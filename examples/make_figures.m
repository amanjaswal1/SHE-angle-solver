% MAKE_FIGURES  Regenerate every figure and number quoted in README.md.
%
%   run('examples/make_figures.m')   ->  results/figures/*.png + printed numbers
%   MATLAB R2019b+ or GNU Octave 8+, no toolboxes. All starting guesses come
%   from a Halton sequence (no random numbers), so MATLAB and Octave produce
%   identical results. ~2 min in Octave.

clear; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));
out_dir = fullfile(here, '..', 'results', 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
heavy_check = false;   % true: repeat the branch search with 3x more starting points (~5 min)

H = [5 7 11 13 17];                  % the v1 problem: 6 angles
m_grid = 0.02:0.01:1.27;             % up to the square-wave limit 4/pi
angle_colors = [0.000 0.447 0.741; 0.850 0.325 0.098; 0.929 0.694 0.125;
                0.494 0.184 0.556; 0.466 0.674 0.188; 0.301 0.745 0.933];

%% Solution map: every branch over m, and its line-to-line THD -------------
t0 = tic;
br = she_sweep(m_grid, H, 100, 3);   % 100 Halton starts at every 3rd grid point
t_sweep = toc(t0);

% THD along each branch, and the best implementable branch at m = 0.6:
% lowest line THD among patterns whose narrowest pulse is >= 2 degrees.
thd = cell(size(br));
m_star = 0.6; min_pulse_deg = 2;
best = [Inf, 0, 0];                 % [thd, branch, index]
for i = 1:numel(br)
    thd{i} = zeros(size(br{i}.m));
    for k = 1:numel(br{i}.m)
        thd{i}(k) = she_line_thd(br{i}.alpha(k, :));
    end
    k6 = find(abs(br{i}.m - m_star) < 1e-9, 1);
    if ~isempty(k6) && thd{i}(k6) < best(1) && ...
            she_min_pulse(br{i}.alpha(k6, :)) * 180/pi >= min_pulse_deg
        best = [thd{i}(k6), i, k6];
    end
end
m_max = max(cellfun(@(b) max(b.m), br));

fig = figure('Visible', 'off', 'Color', 'w');
subplot(3, 1, [1 2]); hold on;
for i = 1:numel(br)
    for k = 1:6
        plot(br{i}.m, br{i}.alpha(:, k) * 180/pi, '-', 'Color', angle_colors(k, :), 'LineWidth', 1.6);
    end
end
for k = 1:6
    h(k) = plot(NaN, NaN, '-', 'Color', angle_colors(k, :), 'LineWidth', 2.5); %#ok<SAGROW>
end
plot([m_star m_star], [0 90], 'k:', 'LineWidth', 1);
grid on; box on; xlim([0 1.3]); ylim([0 90]);
ylabel('switching angle [deg]');
legend(h, {'\alpha_1', '\alpha_2', '\alpha_3', '\alpha_4', '\alpha_5', '\alpha_6'}, ...
       'Location', 'northwest', 'NumColumns', 3);
title(sprintf('Solution branches found: %d, feasible up to m = %.2f (square wave: 1.27)', numel(br), m_max));

subplot(3, 1, 3); hold on;
branch_colors = [0.2 0.2 0.2; 0.000 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188; ...
                 0.494 0.184 0.556; 0.635 0.078 0.184];
for i = 1:numel(br)
    c = branch_colors(1 + mod(i - 1, size(branch_colors, 1)), :);
    plot(br{i}.m, 100 * thd{i}, '-', 'Color', c, 'LineWidth', 1.8);
end
plot(m_star, 100 * best(1), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
grid on; box on; xlim([0 1.3]); ylim([0 150]);
xlabel('modulation index m (fundamental, per unit of DC level)');
ylabel('line THD [%]');
title('Same harmonics eliminated, different quality: rank branches by line-to-line THD');
save_png(fig, fullfile(out_dir, 'fig1_solution_map.png'), 1000, 900);

fprintf('\nSolution map (eliminate %s):\n', mat2str(H));
fprintf('  %d branches traced over m = %.2f..%.2f in %.1f s\n', numel(br), m_grid(1), m_grid(end), t_sweep);
for i = 1:numel(br)
    fprintf('  branch %d: m = %.2f .. %.2f\n', i, min(br{i}.m), max(br{i}.m));
end
fprintf('  highest feasible m found: %.2f\n', m_max);
if heavy_check
    br2 = she_sweep(m_grid, H, 200, 2);
    fprintf('  heavier search (200 starts at every 2nd point): %d branches, highest m %.2f\n', ...
            numel(br2), max(cellfun(@(b) max(b.m), br2)));
end
fprintf('  at m = %.1f (branch: line THD, narrowest pulse):\n', m_star);
for i = 1:numel(br)
    k6 = find(abs(br{i}.m - m_star) < 1e-9, 1);
    if ~isempty(k6)
        fprintf('    branch %d: %.1f %%, %.2f deg\n', i, 100*thd{i}(k6), ...
                she_min_pulse(br{i}.alpha(k6, :)) * 180/pi);
    end
end
fprintf('  chosen: branch %d (lowest THD with narrowest pulse >= %g deg)\n', best(2), min_pulse_deg);

%% Waveform and spectrum of the lowest-THD solution at m = 0.6 --------------
alpha = br{best(2)}.alpha(best(3), :);
fprintf('\nChosen pattern at m = %.1f (degrees): %s\n', m_star, mat2str(alpha * 180/pi, 5));
[theta, v] = she_waveform(alpha, 3 * 2^14);
n = 1:49;
b = she_harmonics(alpha, n).';           % row, like b_fft
X = fft(v) / numel(v);
b_fft = -2 * imag(X(n + 1));

fig = figure('Visible', 'off', 'Color', 'w');
subplot(2, 1, 1);
plot(theta * 180/pi, v, 'Color', palette('blue'), 'LineWidth', 1.6); hold on;
for k = 1:numel(alpha)
    plot(alpha(k) * 180/pi * [1 1], [-1.2 1.2], ':', 'Color', palette('grey'));
end
grid on; xlim([0 360]); ylim([-1.25 1.25]);
set(gca, 'XTick', 0:45:360);
xlabel('angle [deg]'); ylabel('level [p.u.]');
title(sprintf('Chosen pattern at m = %.1f (dotted: the 6 solved angles)', m_star));

subplot(2, 1, 2); hold on;
odd = mod(n, 2) == 1;
trip = odd & mod(n, 3) == 0;
elim = ismember(n, H);
keep = odd & ~trip & ~elim;
stem_bar(n(keep), abs(b(keep)), palette('blue'));
stem_bar(n(trip), abs(b(trip)), [0.70 0.70 0.70]);
stem_bar(n(elim), max(abs(b(elim)), 4e-3), palette('red'));   % drawn as a stub
plot(n(odd), abs(b_fft(odd)), 'k.', 'MarkerSize', 10);
hl = [plot(NaN, NaN, '-', 'Color', palette('blue'), 'LineWidth', 5), ...
      plot(NaN, NaN, '-', 'Color', [0.70 0.70 0.70], 'LineWidth', 5), ...
      plot(NaN, NaN, '-', 'Color', palette('red'), 'LineWidth', 5), ...
      plot(NaN, NaN, 'k.', 'MarkerSize', 10)];
grid on; box on; xlim([0 50]); ylim([0 1.1 * max(abs(b))]);
set(gca, 'XTick', [1 5 7 11 13 17 19 23 25 29 31 35 37 41 43 47 49]);
xlabel('harmonic order n'); ylabel('|b_n| [p.u.]');
legend(hl, {'kept harmonics', 'triplen (cancel line-to-line)', 'eliminated (exactly 0)', ...
        'FFT of the waveform'}, 'Location', 'northeast');
title('Spectrum: closed-form coefficients (bars) agree with an FFT of the waveform (dots)');
save_png(fig, fullfile(out_dir, 'fig2_waveform_spectrum.png'), 1000, 820);

fprintf('  |b_n| of eliminated harmonics: %s\n', mat2str(abs(b(elim)), 3));
fprintf('  first uneliminated non-triplen harmonic: n = 19, |b_19| = %.4f\n', abs(b(19)));
fprintf('  max |closed form - FFT| over odd n <= 49: %.1e\n', max(abs(b(odd) - b_fft(odd))));

%% Newton convergence ------------------------------------------------------
fig = figure('Visible', 'off', 'Color', 'w'); hold on;
n_runs = 0; its = [];
starts = sort(halton_points(1:200, 6), 2) * pi / 2;   % deterministic starting guesses
for s = 1:200
    [~, info] = she_solve(m_star, H, starts(s, :));
    if info.feasible
        its(end + 1) = info.iterations; %#ok<SAGROW>
        if n_runs < 8
            semilogy(0:numel(info.history) - 1, info.history, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4);
            n_runs = n_runs + 1;
        end
    end
end
set(gca, 'YScale', 'log');
grid on; box on; ylim([1e-17 10]);
xlabel('Newton iteration'); ylabel('||F(\alpha)||');
set(gca, 'YTick', 10 .^ (-16:4:0));
title('Damped Newton, analytic Jacobian: quadratic convergence near a root');
save_png(fig, fullfile(out_dir, 'fig3_newton_convergence.png'), 900, 520);

[~, info] = she_solve(m_star, H, [10 20 30 40 50 60] * pi / 180);
e = info.history(info.history > 1e-14 & info.history < 1e-1);   % local phase
p = log(e(3:end) ./ e(2:end-1)) ./ log(e(2:end-1) ./ e(1:end-2));
fprintf('\nNewton (from the v1 README guess): %d iterations, final ||F|| = %.1e\n', ...
        info.iterations, info.residual_norm);
fprintf('  residual history: %s\n', mat2str(info.history, 2));
fprintf('  observed order once ||F|| < 0.1: %s\n', mat2str(p, 3));
fprintf('Starting guesses at m = %.1f: %d of 200 reached a feasible solution, median %g iterations\n\n', ...
        m_star, numel(its), median(its));
