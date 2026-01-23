function makeAxesDataSquare(ax, padFrac)
    if nargin < 2 || isempty(padFrac), padFrac = 0.02; end
    drawnow;

    xl = xlim(ax);
    yl = ylim(ax);

    dx = diff(xl); dy = diff(yl);
    if dx <= 0, dx = 1; xl = xl + [-0.5 0.5]*dx; end
    if dy <= 0, dy = 1; yl = yl + [-0.5 0.5]*dy; end

    xl = xl + [-1 1]*padFrac*dx;
    yl = yl + [-1 1]*padFrac*dy;

    spanX = diff(xl);
    spanY = diff(yl);
    span  = max(spanX, spanY);

    cx = mean(xl);
    cy = mean(yl);

    xlim(ax, [cx - span/2, cx + span/2]);
    ylim(ax, [cy - span/2, cy + span/2]);

    pbaspect(ax,[1 1 1]);
    ax.LooseInset = [0 0 0 0];
end