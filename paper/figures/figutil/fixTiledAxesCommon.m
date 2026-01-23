function fixTiledAxesCommon(ax)
    set(ax, 'Visible','on');
    set(ax, 'XTick',[],'YTick',[], 'XTickLabel',[],'YTickLabel',[]);
    pbaspect(ax,[1 1 1]);
    box(ax,'on');
    ax.LooseInset = [0 0 0 0];
end