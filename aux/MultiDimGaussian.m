function N = MultiDimGaussian(x, mu, sigma)
    % History:
    %   Shrisha Rao - Created  

% dimension = length(mu)

dimension = length(mu);


N = (1 / 2 * pi) ^ (0.5 * dimension) * sqrt(1 / det(sigma)) * exp(- 0.5 * (x - mu)' * inv(sigma) * (x - mu));
end

% mu = [0, 0]';
% sigma = 10* eye(2);
% x = -10:.1:10;
% y = -10:.1:10;
% [X, Y] = meshgrid(x,y);
% X = X(:);
% Y = Y(:);
% 
% for i = 1: length(X)
%     y(i) = MultiDimGaussian([X(i);Y(i)],mu,sigma);
% end
% 
% y =reshape(y,length(x),[]);
% 
% mesh(y)