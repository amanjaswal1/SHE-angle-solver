% DEMO_V1_PROBLEM  Solve the original v1 problem with the v2 solver.
%
%   Current-source inverter, DC-link current I_dc = 200 A. Choose 6 switching
%   angles so the fundamental of the PWM line current is 120 A (m = 0.6) and
%   harmonics 5, 7, 11, 13, 17 are eliminated.

clear;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));

I_dc = 200;                           % DC-link current [A]
I_1_target = 120;                     % desired fundamental amplitude [A]
H = [5 7 11 13 17];                   % harmonics to eliminate
alpha0 = [10 20 30 40 50 60] * pi/180;   % initial guess (as in v1)

[alpha, info] = she_solve(I_1_target / I_dc, H, alpha0);

fprintf('converged: %d, feasible: %d, iterations: %d, ||F|| = %.1e\n', ...
        info.converged, info.feasible, info.iterations, info.residual_norm);
fprintf('angles [deg]: %s\n', sprintf('%.4f  ', alpha * 180/pi));
fprintf('narrowest pulse: %.2f deg\n', she_min_pulse(alpha) * 180/pi);
amps = I_dc * she_harmonics(alpha, [1 H]);
fprintf('I_1 = %.6f A\n', amps(1));
fprintf('I_%d = %.1e A\n', [H; amps(2:end).']);
fprintf('line-to-line THD: %.1f %%\n', 100 * she_line_thd(alpha));
