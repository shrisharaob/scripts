function out = ButterBandPass(x, fs, fRange, varargin)
% ButterBandPass(x, fs, fRange, varargin)

[order] = DefaultArgs(varargin, {6});

[b a ]=butter(order,fRange(1) / (fs/2), 'high');
tempx=filtfilt(b,a,x);
[b a ]=butter(order,fRange(2) / (fs/2), 'low');

out=filtfilt(b,a,tempx);

end