function thd = she_line_thd(alpha, n_points)
%SHE_LINE_THD  THD of the line-to-line waveform produced by a SHE pattern.
%
%   thd = she_line_thd(alpha)
%   thd = she_line_thd(alpha, n_points)
%
%   In a balanced three-phase, three-wire system the triplen harmonics
%   (3, 9, 15, ...) cancel between phases, so the quality of a pattern is
%   judged on the line-to-line quantity v_ab(t) = v_a(t) - v_a(t - 2*pi/3).
%
%   THD = sqrt(RMS^2 - V1^2/2) / (V1/sqrt(2)),  V1 = sqrt(3) * b_1.
%
%   The RMS comes from the sampled waveform (n_points must be divisible by
%   3 so the 120-degree shift is exact; default 3*2^14), the fundamental
%   from the closed-form coefficient. Returned as a fraction (0.1 = 10 %).

    if nargin < 2, n_points = 3 * 2^14; end
    if mod(n_points, 3) ~= 0
        error('she_line_thd:n_points', 'n_points must be divisible by 3.');
    end
    [~, va] = she_waveform(alpha, n_points);
    vab = va - circshift(va, [0, n_points / 3]);     % v_a(t) - v_a(t - 120 deg)
    V1 = sqrt(3) * abs(she_harmonics(alpha, 1));
    thd = sqrt(max(mean(vab .^ 2) - V1^2 / 2, 0)) / (V1 / sqrt(2));
end
