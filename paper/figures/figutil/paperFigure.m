function fig = paperFigure(style, varargin)
%PAPERFIGURE  Create a figure with consistent publication sizing (centimeters).
%
%   fig = paperFigure()          % defaults to 'double'
%   fig = paperFigure('double')  % standard 1x2 comparison figure footprint
%   fig = paperFigure('single')  % single-column footprint
%   fig = paperFigure('square')  % square footprint (e.g. 2x2 panels)
%
% Optional name-value:
%   'Position' : [x y w h] in centimeters (overrides preset)
%   'Color'    : figure color (default 'w')
%
% Notes:
% - exportgraphics respects the on-screen figure size well.
% - Keep all sizing changes centralized here for reproducibility.

if nargin < 1 || isempty(style)
    style = 'double';
end

ip = inputParser;
ip.addParameter('Position', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==4));
ip.addParameter('Color', 'w');
ip.parse(varargin{:});
opt = ip.Results;

style = lower(string(style));

% ---- presets (cm) ----
% Tune these ONCE and everything downstream stays consistent.
switch style
    case "double"
        % Your default for 1x2 comparison figures (two square-ish panels)
        % Adjust if needed later in ONE place.
        pos = [2 2 18 9];   % [x y width height] cm

    case "single"
        % Single-column style
        pos = [2 2 9 7];    % cm

    case "square"
        % Square-ish canvas (e.g. 2x2 panels)
        pos = [2 2 14 14];  % cm

    otherwise
        error('paperFigure:BadStyle', ...
            'Unknown style "%s". Use ''double'',''single'',''square''.', style);
end

if ~isempty(opt.Position)
    pos = opt.Position;
end

% ---- create figure ----
fig = figure('Color', opt.Color);
set(fig, ...
    'Units', 'centimeters', ...
    'Position', pos, ...
    'PaperUnits', 'centimeters', ...
    'PaperPositionMode', 'auto', ...
    'InvertHardcopy', 'off');

% Vector-friendly default
try
    set(fig, 'Renderer', 'painters');
catch
    % ignore if renderer set fails in some contexts
end

end
