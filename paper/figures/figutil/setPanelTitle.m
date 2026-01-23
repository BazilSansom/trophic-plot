function h = setPanelTitle(ax, txt, S)
%SETPANELTITLE  Consistent panel titles (single source of truth).
%
%   h = setPanelTitle(ax, txt)
%   h = setPanelTitle(ax, txt, S)   % style struct from paperStyle()
%
% Uses paperStyle() defaults if S not provided.

if nargin < 3 || isempty(S)
    S = paperStyle();
end

if nargin < 2
    txt = "";
end

% Defaults / fallbacks
fs = [];
fw = 'normal';
interp = 'none';

if isfield(S,'TitleFontSize') && ~isempty(S.TitleFontSize)
    fs = S.TitleFontSize;
end
if isfield(S,'TitleFontWeight') && ~isempty(S.TitleFontWeight)
    fw = char(S.TitleFontWeight);
end
if isfield(S,'TitleInterpreter') && ~isempty(S.TitleInterpreter)
    interp = char(S.TitleInterpreter);
end

% Apply
if isempty(fs)
    h = title(ax, txt, 'FontWeight', fw, 'Interpreter', interp);
else
    h = title(ax, txt, 'FontSize', fs, 'FontWeight', fw, 'Interpreter', interp);
end

end
