function runAllTests()
%RUNALLTESTS Entry point for local runs and CI.
%
% Delegates to tools/run_matlab_tests.m (single source of truth).

    % Find repo root from this file location
    repoRoot = fileparts(mfilename("fullpath"));

    runner = fullfile(repoRoot, "tools", "run_matlab_tests.m");
    assert(isfile(runner), "Test runner not found: %s", runner);

    % Call the runner
    run(runner);
end
