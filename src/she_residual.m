function F = she_residual(alpha, m, harmonics)
%SHE_RESIDUAL  SHE equations F(alpha) = 0.
%
%   F = she_residual(alpha, m, harmonics)
%
%   alpha      N switching angles [rad].
%   m          Target fundamental amplitude, per unit of the DC level
%              (the full square wave has 4/pi = 1.273).
%   harmonics  The N-1 odd harmonic orders to eliminate.
%
%   F(1)   = b_1(alpha) - m          (fundamental control)
%   F(i+1) = b_h(i)(alpha)           (harmonic h(i) set to zero)

    if numel(alpha) ~= numel(harmonics) + 1
        error('she_residual:size', ...
              'Need N = numel(harmonics) + 1 angles (got %d angles, %d harmonics).', ...
              numel(alpha), numel(harmonics));
    end
    F = she_harmonics(alpha, [1, harmonics(:).']);
    F(1) = F(1) - m;
end
