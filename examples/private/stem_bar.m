function stem_bar(x, y, c)
%STEM_BAR  Thin vertical bars (portable alternative to bar/stem styling).
    for i = 1:numel(x)
        plot([x(i) x(i)], [0 y(i)], '-', 'Color', c, 'LineWidth', 5);
    end
end
