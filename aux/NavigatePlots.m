function NavigatePlots(plotFunc, loopPar, varargin)
% NavigatePlots(plotFunc, lopPar) 
% looppar is the last arg of plotFunc 
% [IF_REPORTFIG, IF_WAITFORBUTTONPRESS,IF_GCF, filebase, fileTag ]  -                         
% [ 0, 1, 0, '', []]
%
[IF_REPORTFIG, IF_WAITFORBUTTONPRESS,IF_GCF, filebase, fileTag ] = DefaultArgs(varargin, {0, 1, 0,[],[]});
 fieldNames = fieldnames(plotFunc);
 
    plotFuncHandle = plotFunc.(fieldNames{1});

   IS_FULL_CYCLE = 0; 
   persistent IS_DONE;
   IS_DONE = 0;
   nDisplayed= 0;
    kLoopIndx = 1;
    nItems = length(loopPar);
    if ~IF_GCF, h = figure; else h = gcf; hold on; end
    set(h,'CloseRequestFcn',@my_closefcn);
          [~, nClmns] = size(loopPar);
          if nClmns > 2  
              loopPar = loopPar(:);
          end
    while(~IS_DONE)
        
        nDisplayed = nDisplayed + 1;
        if nDisplayed ~= 1 && ~IF_GCF, clf; end
        if nDisplayed ~= 1, clf; end
        if iscell(plotFunc.args)
            feval(plotFuncHandle, plotFunc.args{:} , loopPar(kLoopIndx, :));
        elseif  isempty(plotFunc.args)
            feval(plotFuncHandle, loopPar(kLoopIndx, :));
        else
            feval(plotFuncHandle, plotFunc.args , loopPar(kLoopIndx, :));
        end
        if ishandle(h)
            if IF_WAITFORBUTTONPRESS
                try
                    waitforbuttonpress;
                    navKey = double(get(gcf,'CurrentCharacter'));
                    if navKey == 28 | navKey == 31
                        dir = 'down';
                    elseif navKey == 30 | navKey == 29
                        dir = 'up';
                    elseif navKey == double('x')
                        dir = 0;
                    else
                        dir = -1;
                    end
                   
                catch err
                    IS_DONE = 1;
                end
            else
%                 pause(1);
                disp(num2str(nDisplayed));
                if kLoopIndx == nItems  
                    dir = 0; % exit
                else
                    dir = 'up'; 
                    
                end
            end
            
                                
            switch dir
                case 'up'
                    kLoopIndx = mod(kLoopIndx + 1, nItems);
                    if mod(kLoopIndx, nItems) == 0
                        kLoopIndx = 1;
                    end
                case 'down'
                    kLoopIndx = mod(kLoopIndx - 1, nItems);
                    if mod(kLoopIndx, nItems) == 0
                        kLoopIndx = nItems;
                    end
                case 0
                    IS_DONE = 1;
                otherwise
                    kLoopIndx = 1;
            end
            
            
        else return;
        end

        if IF_REPORTFIG && ~IS_FULL_CYCLE
            if isempty(filebase)
                if isa(plotFunc.args{1}, 'GenericPF')
                   filebase = plotFunc.args{1}.filebase;
                else
                   filebase = input('enter filebase :');
                end
            end
            if isempty(fileTag)
                fileName = [filebase '.' func2str(plotFuncHandle)];
            else
                fileName = [filebase '.' func2str(plotFuncHandle), '.' fileTag '.'];
            end
            if nDisplayed == nItems
                IS_FULL_CYCLE = 1;
            end
            reportfig(gcf,fileName , 0, [num2str(nDisplayed)], [],0);
        end           
    end
end
    
    function my_closefcn(src, evnt)
    delete(gcf);
        IS_DONE = 1;
    end
        