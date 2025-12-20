function tests = testTrophicLevels
%TFL_TEST_01_TROPHIC_LEVELS  Unit tests for trophic_levels.m
%
% Uses MATLAB's function-based unit test style so it plays nicely with
% runtests and CI.

    tests = functiontests(localfunctions);
end


%% ---- Test 1: simple 3-node feed-forward chain ----
function testSimpleChain(testCase)

    W = sparse([0 1 0; ...
                0 0 1; ...
                0 0 0]);

    [h, info] = trophic_levels(W, 'ComputeCoherence', true);

    % One component, no isolates
    verifyEqual(testCase, info.nComp, 1);
    verifyEqual(testCase, numel(h), 3);

    % Monotone increasing levels (up to gauge)
    verifyGreaterThanOrEqual(testCase, h(2) - h(1), 0);
    verifyGreaterThanOrEqual(testCase, h(3) - h(2), 0);

    % Some layering: coherence in (0,1], incoherence in [0,1)
    verifyGreaterThanOrEqual(testCase, info.C_global, 0);
    verifyLessThanOrEqual(testCase, info.C_global, 1);
end


%% ---- Test 2: single bidirectional pair (symmetric 2-node) ----
function testSymmetricPair(testCase)

    W = sparse([0 1; ...
                1 0]);

    h = [];
    info = struct();

    testCase.verifyWarning(@runAndCapture, 'trophic:SymmetricAdjacency');

    % Connected, 2 nodes
    verifyEqual(testCase, info.nComp, 1);
    verifyEqual(testCase, numel(h), 2);

    % Symmetry => both nodes same level (up to rounding)
    verifyEqual(testCase, h(1), h(2), 'AbsTol', 1e-10);

    % Many edges frustrated => incoherence > 0
    verifyGreaterThan(testCase, info.F0_global, 0);
    verifyLessThanOrEqual(testCase, info.F0_global, 1);

    function runAndCapture()
        [h, info] = trophic_levels(W, 'ComputeCoherence', true);
    end
end



%% ---- Test 3: multi-component + isolates ----
function testMultiComponentAndIsolates(testCase)

    % Component A: 1 -> 2 -> 3
    % Component B: 4 -> 5
    % Isolate:     6
    W = sparse(6,6);
    W(1,2) = 1;
    W(2,3) = 1;
    W(4,5) = 1;

    [h, info] = trophic_levels(W, 'ComputeCoherence', true);

    % Three weak components: {1,2,3}, {4,5}, {6}
    verifyEqual(testCase, info.nComp, 3);

    % Isolate height is 0 by convention
    verifyEqual(testCase, h(6), 0);

    % Component-level incoherence:
    %   - Component with edges should have 0 <= F0 <= 1
    %   - Edge-free component (the isolate) has F0 = 0 by convention
    verifyEqual(testCase, numel(info.F0_comp), 3);


    % Edge-free component (the isolate): coherence metrics undefined => NaN
    cIso = info.compIdx(6);
    verifyTrue(testCase, isnan(info.F0_comp(cIso)));
    verifyTrue(testCase, isnan(info.C_comp(cIso)));

end


%% ---- Test 4: h0 options on a single component ----
function testH0OptionsSingleComponent(testCase)

    % Simple 3-node chain
    W = sparse([0 1 0; ...
                0 0 1; ...
                0 0 0]);

    [h_min, info_min] = trophic_levels(W, 'h0', 'min', 'ComputeCoherence', false);
    [h_wm,  info_wm]  = trophic_levels(W, 'h0', 'wm',  'ComputeCoherence', false);
    [h_sm,  info_sm]  = trophic_levels(W, 'h0', 'sm',  'ComputeCoherence', false);

    % Single component in each case
    verifyEqual(testCase, info_min.nComp, 1);
    verifyEqual(testCase, info_wm.nComp, 1);
    verifyEqual(testCase, info_sm.nComp, 1);

    % After centring by simple mean, all choices differ only by numerical noise
    d_min = h_min - mean(h_min);
    d_wm  = h_wm  - mean(h_wm);
    d_sm  = h_sm  - mean(h_sm);

    verifyLessThan(testCase, norm(d_min - d_wm), 1e-8);
    verifyLessThan(testCase, norm(d_min - d_sm), 1e-8);
end


%% ---- Test 5: custom random targets (tau != 1) ----
function testRandomTargets(testCase)

    rng(123);  % reproducible

    % Simple 3-node chain
    W = sparse([0 1 0; ...
                0 0 1; ...
                0 0 0]);

    n = size(W,1);
    Tau = zeros(n);

    [iE, jE] = find(W);
    Tau(sub2ind(size(Tau), iE, jE)) = 0.5 + 1.5*rand(numel(iE),1);  % U(0.5,2.0)

    [h, info] = trophic_levels(W, 'targets', Tau, 'ComputeCoherence', true); %#ok<ASGLU>

    % F0 is a non-negative weighted mean of squared deviations, but not
    % necessarily <= 1 for arbitrary targets. Just enforce basic sanity.
    verifyTrue(testCase, isfinite(info.F0_global));
    verifyGreaterThanOrEqual(testCase, info.F0_global, 0);
end
