function fig = init(FontSize, FontName)
    if nargin < 1
        FontSize = 18;
    end
    if nargin < 2
        FontName = 'Helvetica';
    end
    fig = figure;
    set(gcf, 'color', 'w');
    set(gcf, 'PaperPositionMode', 'auto'); % save the fig how it looks on screen
    plot.subplot(nan, nan, nan, FontSize, FontName);
end
