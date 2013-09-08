% cirmean - calculates mean phase and radius of a given
%           circular statistics
%
% function [theta, r]  = circmean(data,returnmode);
% (data has to be in the unit of radians)
% if returnmode = 1
function [theta, r]  = circmean(data, returnmode)
if nargin<2 |isempty(returnmode)
    returnmode =0;
end
if returnmode & nargout>1
    error('this mode is returning both theta and r as two column vector');
end
   
% convert the data into complex represenation

% if length(data_rad)==0
% 	theta = 0;
% 	r = 0;
% 	return;
% end
% 
% data = data_rad(find(~isnan(data_rad)));
data = data(find(~isnan(data)));
datai = exp(data * 1i);

datas = mean(datai);
r = abs(datas)./repmat(size(data,1),1,size(data,2));

if returnmode 
    theta = [angle(datas) abs(datas)];   
else
    theta = angle(datas);
end
