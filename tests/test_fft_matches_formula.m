function test_fft_matches_formula()
% The closed-form Fourier coefficients must agree with an FFT of the
% switching waveform built directly from its definition (independent check).
    alpha = [7 14 26 31 44 58] * pi / 180;
    n_pts = 2^16;
    [~, v] = she_waveform(alpha, n_pts);
    X = fft(v) / n_pts;
    n = 1:2:49;
    b_fft = -2 * imag(X(n + 1));
    b_formula = she_harmonics(alpha, n).';
    err = max(abs(b_fft - b_formula));
    assert(err < 3e-4, sprintf('FFT vs formula: max error %.3g', err));   % ~ 2*pi/n_pts
    assert(max(abs(X(3:2:51))) < 1e-12, 'even harmonics (X(n+1), n even) should be zero');
end
