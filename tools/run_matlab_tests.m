% tools/run_matlab_tests.m
% Robust test runner (works across MATLAB versions)

% repoRoot = parent folder of this file's folder (.../trophic-plot/tools -> .../trophic-plot)
repoRoot = fileparts(fileparts(mfilename("fullpath")));

addpath(genpath(fullfile(repoRoot, "src")));

testsFolder = fullfile(repoRoot, "tests");
assert(isfolder(testsFolder), "Tests folder not found: %s", testsFolder);

suite = matlab.unittest.TestSuite.fromFolder(testsFolder, "IncludingSubfolders", true);
assert(~isempty(suite), "No tests discovered under: %s", testsFolder);

results = run(suite);
disp(table(results));
assertSuccess(results);
