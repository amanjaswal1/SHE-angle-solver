function test_line_thd()
% Six-step check: a single pulse from 0 to pi/2 in the first quarter (one
% angle alpha = 0+) is a square wave; its line-to-line waveform is the
% quasi-square six-step wave whose THD is known: sqrt(pi^2/9 - 1) = 31.08 %.
    thd = she_line_thd(1e-9);
    assert(abs(thd - sqrt(pi^2/9 - 1)) < 1e-3, sprintf('six-step THD %.5f', thd));
end
