function add_toolbox_path()
%ADD_TOOLBOX_PATH Add toolbox folder (and subfolders) to MATLAB path.
repoRoot = fileparts(fileparts(mfilename('fullpath')));  % .../paper -> repo root
toolboxRoot = fullfile(repoRoot, "toolbox");
addpath(genpath(toolboxRoot));
end
