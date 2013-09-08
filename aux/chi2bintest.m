function [p, Q]= chi2bintest(x, n)
	
% Usage: [p, Q]= chi2bintest(x, n)
% 
% The chi-squared binary test. 
% 
% Given a number of samples with a binary outcome ("is or is not something"), 
% this function tests the hypothesis that the samples are independent. 
% If Q > chi2(p, nu), the hypothesis is rejected. 
% 
% If x is a vector, each item is the result of a sample.
% n must be a scalar (same total for all samples) or a vector of the same length as x.
% 
% If x is a matrix, each row is the result of a sample. Each column is an independent 
% sample series. n must be a scalar (same total for all samples) or 
% a vector of the same length as the number of samples in x (same totals for all sample 
% series) or a matrix of the same size as x. 
% 
% If you find any errors, please let me know: .
% 
% ARGUMENTS:
% x     Relative frequencies of the outcome, if 0<=x<=1.
%       Absolut number count, if integer x>=0.
% n     Number of tries.
% p     The prob ability value, calculated from Q.
% Q     The resulting Q-value.
% 
% EXAMPLE 1
% In region A, 324 of 556 cows were red, whereas in region B 98 of 260 were red.
% [p, Q]= chi2bintest([324; 98], [556; 260])
% p=
%    4.2073e-08
% Q=
%    30.0515
% With an error risk of about 4e-08, we can claim that the samples are independent.
% 
% EXAMPLE 2
% In three regions, reindeers were cheched for four various deseases.
% p= chi2bintest([28,39,8,4782; 17,7,5,903; 21,14,6,2322], [15996; 10127; 12476])
% p =
%    0.9869    0.0007    0.9973         0
% It seems that the first and third deseases are region independent, whereas the others 
% have regional dependencies (are independent to each other).
% 
% SEE ALSO:   chi2test, soon to be uploaded to 
%             http://www.mathworks.com/matlabcentral/fileexchange/
% 
% HISTORY:    v.1.0, first working version, 2007-08-30.
% 
% COPYRIGHT:  (c) 2007 Peder Axensten. Use at own risk.

% KEYWORDS:   chi-squared test, chi-squared, chi2, binary, test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Check the arguments.
	if(nargin ~= 2),			error('Two and only two arguments required!');						end
	if(ndims(x) > 2),			error( 'First argument (x) must be a vector or a 2d matrix!');		end
	if(ndims(n) > 2),			error('Second argument (n) must be a vector or a 2d matrix!');		end
	if(size(x, 1) == 1),		x=	x';																end
	if(size(n, 1) == 1),		n=	n';																end
	
	if(size(n,1) == 1),			n=	repmat(n, size(x, 1), 1);										end
	if(size(n,2) == 1),			n=	repmat(n, 1, size(x, 2));										end
	if(any(size(x) ~= size(n))),error('Size of second argument (n) doesn''t match size of first argument (x)!');	end
	
	if(any(~isreal(x))),		error('All values of first argument (x) must be real values!');		end
	if(any(~isreal(n))),		error('All values of second argument (n) must be real values!');	end
	
	if(any(n<1 | (n ~= floor(n))))
		error('All values of second argument (n) must be positive integers');
	end
	
	% Calculate Q = sum( (a-np*)^2/(np*(1-p*)) )
	s=		size(x, 1);
	ep=		repmat((sum(x)./sum(n)), s, 1);
	Q=		sum((x-n .* ep).^2./(n.*ep.*(1 - ep)));
	
	% Calculate cdf of chi-squared to Q. Degrees of freedom, v, is s-1.
	p=		1 - gammainc(Q/2, (s-1)/2);
end
