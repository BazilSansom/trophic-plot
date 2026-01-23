function h = addPanelLabel(ax, str, varargin)
%ADDPANELLABEL Add a small in-axes label without affecting layout.
%
% h = addPanelLabel(ax, str, 'Location',"nw", 'FontSize',10, 'Background',true)

p = inputParser;
p.addParameter('Location',"nw");      % "nw","ne","sw","se"
p.addParameter('FontSize',10);
p.addParameter('Background',true);
p.addParameter('Margin',1);
p.parse(varargin{:});
o = p.Results;

switch lower(string(o.Location))
    case "nw", x=0.02; y=0.98; ha='left';  va='top';
    case "ne", x=0.98; y=0.98; ha='right'; va='top';
    case "sw", x=0.02; y=0.02; ha='left';  va='bottom';
    case "se", x=0.98; y=0.02; ha='right'; va='bottom';
    otherwise, x=0.02; y=0.98; ha='left';  va='top';
end

args = {'Units','normalized','HorizontalAlignment',ha,'VerticalAlignment',va,...
        'FontSize',o.FontSize,'Interpreter','none'};

if o.Background
    args = [args, {'BackgroundColor','w','Margin',o.Margin}];
end

h = text(ax, x, y, str, args{:});
end
