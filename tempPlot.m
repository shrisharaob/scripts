

 figure;
 hold on;
for kDir = 1 :2
    switch kDir
        case 1
            a2b = ccgSegsA2B.a2b;
            c = 'b';
        otherwise
            a2b = ccgSegsA2B.b2a;
            c = 'r';
    end
    nPairs = length(a2b);
    mA2BOffset = nan(nPairs, 1);
    mA2BfirstPeak = nan(nPairs, 1);
    comDist = nan(nPairs, 1);
    CHECK1 = 0;
    CHECK2 = 0;
    for mPair = 1 : nPairs
        
        if ~cellfun(@isempty, a2b(mPair))
            if ~isempty(a2b{mPair}.offset)
                mA2BOffset(mPair) = a2b{mPair}.offset;
                CHECK1 = 1;
            end
            if ~isempty(a2b{mPair}.firstPeak)
                mA2BfirstPeak(mPair) = a2b{mPair}.firstPeak;
                CHECK2 = 1;
            end
            if CHECK1 && CHECK2, comDist(mPair) = pfPars.comDist(mPair); end
        end
    end
    
   
    if kDir == 1, plot(comDist, mA2BOffset,'*', 'Color', c);
    else plot(comDist, -mA2BOffset,'*', 'Color', c); end
    grid on;
%     hold on;
%     plot(comDist, mA2BfirstPeak,'*', 'color', c)
end
legend({'a2b','b2a'})
