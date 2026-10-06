function test_input_checks()
% Wrong number of angles must raise a clear error.
    threw = false;
    try, she_residual([0.1 0.2], 0.5, [5 7]); catch, threw = true; end
    assert(threw, 'size mismatch should error');
end
