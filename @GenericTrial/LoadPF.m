function trial = LoadPF(trial, varargin)
% Loads place fields
    trial = trial.Load({{'PF', varargin{:}}});
end