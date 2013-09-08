function samples = SampleDistri(x, n)
% samples = SampleDistri(x, n)
% get samples from the emprical distribution of x, small sample size effects ??

    support = unique(x);
    m = length(support);
    counts = histc(x, support);
    p = counts ./ sum(counts);
    uni = rand(1,n);
    cumprob = [0; cumsum(p)];
    samples=zeros(1,n);
    for j=1:m
        ind=find((uni>cumprob(j)) & (uni<=cumprob(j+1)));
        samples(ind)=support(j);
    end
end