% tools/run_matlab_tests.m
% Robust test runner (works across MATLAB versions and from any working dir)

thisFile = mfilename("fullpath");
toolsDir = fileparts(thisFile);
repoRoot = fileparts(toolsDir);  % .../trophic-plot

% Add toolbox code to path (new structure)
toolboxDir = fullfile(repoRoot, "toolbox");
assert(isfolder(toolboxDir), "Toolbox folder not found: %s", toolboxDir);
addpath(genpath(toolboxDir));

% Discover and run tests
testsFolder = fullfile(repoRoot, "tests");
assert(isfolder(testsFolder), "Tests folder not found: %s", testsFolder);

suite = matlab.unittest.TestSuite.fromFolder(testsFolder, "IncludingSubfolders", true);
assert(~isempty(suite), "No tests discovered under: %s", testsFolder);

results = run(suite);
disp(table(results));

% Assert success (older MATLAB compatibility)
if exist("assertSuccess","file") == 2
    assertSuccess(results);
else
    assert(all([results.Passed]), "Some tests failed.");
end
