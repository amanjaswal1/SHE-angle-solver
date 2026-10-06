function b = she_harmonics(alpha, n)
%SHE_HARMONICS  Fourier amplitudes of a quarter-wave symmetric 3-level wave.
%
%   b = she_harmonics(alpha, n)
%
%   alpha  Switching angles [rad], 0 < alpha_1 < ... < alpha_N < pi/2.
%   n      Harmonic orders (odd integers), e.g. [1 5 7 11 13 17].
%   b      Amplitude of each harmonic, per unit of the DC-link level
%          (multiply by I_dc or V_dc for absolute units).
%
%   Waveform over the first quarter period: 0 until alpha_1, +1 until
%   alpha_2, 0 until alpha_3, ... (alternating), extended with half- and
%   quarter-wave symmetry. Its Fourier sine coefficients are
%
%       b_n = 4/(n*pi) * sum_k (-1)^(k-1) * cos(n * alpha_k)
%
%   and all even harmonics and cosine terms are zero.
%
%   See also SHE_RESIDUAL, SHE_JACOBIAN, SHE_WAVEFORM.

    alpha = alpha(:).';                       % row
    n = n(:);                                 % column
    sgn = (-1) .^ (0:numel(alpha) - 1);       % +1, -1, +1, ...
    b = (4 ./ (n * pi)) .* (cos(n * alpha) * sgn.');
end
