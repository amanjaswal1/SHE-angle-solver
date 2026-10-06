function w = she_min_pulse(alpha)
%SHE_MIN_PULSE  Narrowest pulse or notch in a SHE pattern [rad].
%
%   w = she_min_pulse(alpha)
%
%   Real switches need a minimum on/off time, so a pattern whose narrowest
%   pulse is too short cannot be implemented even though it solves the SHE
%   equations exactly. With quarter-wave symmetry the interval widths are
%       2*alpha_1,  alpha_2 - alpha_1, ...,  alpha_N - alpha_N-1,  2*(pi/2 - alpha_N)
%   (the first and last intervals continue across 0 and pi/2 by symmetry).

    alpha = alpha(:).';
    w = min([2 * alpha(1), diff(alpha), 2 * (pi/2 - alpha(end))]);
end
