function out = CorrSignf(timeStamps, varargin)

[IS_DISPLAY, nShuffle, shuffleType] = DefaultArgs(varargin, {0, 200, 'circular'}