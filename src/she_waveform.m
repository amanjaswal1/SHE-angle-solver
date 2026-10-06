function [theta, v] = she_waveform(alpha, n_points)
%SHE_WAVEFORM  Sample one fundamental period of the SHE switching pattern.
%
%   [theta, v] = she_waveform(alpha)
%   [theta, v] = she_waveform(alpha, n_points)
%
%   alpha     Switching angles [rad], increasing, inside (0, pi/2).
%   n_points  Samples per period (default 2^14).
%   theta     Angle grid [rad], 0 <= theta < 2*pi.
%   v         Waveform in per unit of the DC level: values -1, 0 or +1.
%
%   Built directly from the definition (first quarter: 0, +1, 0, ... toggling
%   at each alpha_k; then half- and quarter-wave symmetry), independently of
%   the Fourier formula in SHE_HARMONICS, so an FFT of v can be used to
%   check that formula.

    if nargin < 2, n_points = 2^14; end
    theta = (0:n_points - 1) * (2*pi / n_points);

    % Fold every angle into the first quarter period using the symmetries
    th = mod(theta, pi);                     % half-wave: v(t + pi) = -v(t)
    q  = min(th, pi - th);                   % quarter-wave: v(pi - t) = v(t)
    level = mod(sum(q(:) >= alpha(:).', 2), 2).';   % 0,1,0,1,... after each alpha
    sgn = ones(size(theta));
    sgn(mod(theta, 2*pi) >= pi) = -1;
    v = sgn .* level;
end
