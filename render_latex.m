function fig = render_latex(str, fs, h)
% offscreen figure
fig = figure('Menubar','none', 'Color','white', ...
    'Units','inches', 'Position',[100 100 8 h]);
axis off
str = strrep(str, '\left(\begin{array}', '\left[\begin{array}');
str = strrep(str, '\end{array}\right)', '\end{array}\right]');

text(0.5, 0.5, ['$$' str '$$'], 'Interpreter','latex', 'FontSize',fs, ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle')
% snapnow
