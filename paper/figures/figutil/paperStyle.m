function S = paperStyle(varargin)
%PAPERSTYLE  Centralised style defaults for paper figures.

% Defaults
S = struct();

% Typography
S.FontName        = "";      % e.g. "Times New Roman" (leave "" to use MATLAB default)
S.TitleFontSize   = 10;
S.TitleFontWeight = "normal";
S.TitleInterpreter= "none";

% PlotTFL / tflPlot aesthetics we may want to standardise
S.ShowLabels      = true;

% If we want global node/label sizing consistency across all figures,
% we can pass these through to plotTFL / tflPlot via plotopts:
S.NodeSize        = 50;
S.LabelFontSize   = 9;

% Axes cosmetics
S.PadFrac         = 0.02;
S.TileSpacing     = "compact";
S.Padding         = "compact";

% Optional overrides via name/value
if ~isempty(varargin)
    if mod(numel(varargin),2) ~= 0
        error('paperStyle:BadArgs','Name/value pairs required.');
    end
    for k = 1:2:numel(varargin)
        S.(varargin{k}) = varargin{k+1};
    end
end
end
