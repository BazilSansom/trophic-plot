function [h, F0, eta, info] = incoherence(W, varargin)
%INCOHERENCE  Compute trophic levels and incoherence statistics.
%
%   [h, F0, eta] = incoherence(W)
%   [h, F0, eta, info] = incoherence(W, Name, Value, ...)
%
% Thin compatibility wrapper around TROPHIC_LEVELS.
%
% Inputs
% ------
% W : square nonnegative adjacency matrix
%
% Name-value options
% ------------------
% Passed through to trophic_levels, including for example:
%   'h0'                : 'min' (default), 'wm', or 'sm'
%   'symtest'           : 1 or 0
%   'targets'           : scalar 1 or target matrix
%   'ComputeCoherence'  : true or false
%
% Outputs
% -------
% h    : trophic levels
% F0   : global trophic incoherence
% eta  : sqrt(F0/(1-F0)); Inf when F0=1, 0 when F0=0
% info : full info struct returned by trophic_levels
%
% Notes
% -----
% This function is retained mainly for convenience/backward compatibility.
% For richer output, prefer calling trophic_levels directly.

    [h, info] = trophic_levels(W, varargin{:});

    if isfield(info, 'F0_global') && ~isempty(info.F0_global)
        F0 = info.F0_global;
    else
        % Ensure incoherence is computed if caller disabled it
        [h, info] = trophic_levels(W, varargin{:}, 'ComputeCoherence', true);
        F0 = info.F0_global;
    end

    if F0 <= 0
        eta = 0;
    elseif F0 >= 1
        eta = Inf;
    else
        eta = sqrt(F0 / (1 - F0));
    end
end