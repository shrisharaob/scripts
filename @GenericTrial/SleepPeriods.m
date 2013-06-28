function sleepPeriods = SleepPeriods(gt, varargin)
% sleepPeriods = SleepPeriods(gt)
% returns sleep before and after trial

   [outFs] = DefaultArgs(varargin, {0});
   sts = gt.LoadStatePeriods('SWS', 0);
   preTrialPeriods = sts(sts(:, 2) < gt.trialPeriods(1), :);
   postTrialPeriods = sts(sts(:, 1) > gt.trialPeriods(2), :);
   if outFs ~= 0 
       preTrialPeriods = ConvertFs(preTrialPeriods, gt.lfpSampleRate, outFs);
       postTrialPeriods = ConvertFs(postTrialPeriods, gt.lfpSampleRate, outFs);
   end
   sleepPeriods = {preTrialPeriods, postTrialPeriods};
end


