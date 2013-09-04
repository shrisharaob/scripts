function trial = LoadCR(trial)
% trial = LoadCR(trial)
% Loads clu and res @ sample rate
    trial = trial.Load({{'CluRes'}});
end