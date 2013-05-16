function out=ButterBandPass(x,fs,f1,f2,varargin)
%x,fs,f1,f2,order(default=6)

if nargin>4 && ~isempty(varargin(1)) 
    order=varargin(1);else order=6;

[b a ]=butter(order,f1/(fs/2),'high');
tempx=filtfilt(b,a,x);
[b a ]=butter(order,f2/(fs/2),'low');
out=filtfilt(b,a,tempx);

end