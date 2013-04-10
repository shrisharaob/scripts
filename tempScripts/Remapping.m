function Remapping(filebase, varargin)
% Remapping(trial)
% 
[arena, roi] = DefaultArgs(varargin, {{'bigSquare', 'midSquare'}, {'CA3'});
gt = GenericTrial();
gt.filebase = filebase;
trialNames = [];
for kArena = 1 : length(arena)
trialNames = [trialNames; cell2mat(Beh(strcmp(Beh(:,5) , arena{kArena} &  strcmp(Beh(:,2),filebase), 4)];
end

nTrials = length(trialNames);  
for ktr = 1 : nTrials
   gt.trialName = trialName{kTr};
   gt = gt.LoadPF; 
   if ~isempty(gt.pfObject)
clus = cell2mat(gt.GetRegionClu(roi);
