function PlotNeuronQuality(filebase)
    IS_DONE = 1;
    load(['../analysis/' filebase '.NeuronQuality.mat']);
    
    
    nElectrodes = length(nq);
    fields = fieldnames(nq) ;
    nFields = length(fieldnames(nq) );

    h = figure;
    kField = 1;
    while(IS_DONE)
%     for kField = 1 : nFields
        for kElectrode = 1 : nElectrodes
            subplot(4, 3, kElectrode);
            
            plot(nq(kElectrode).(genvarname(fields{kField}))','*-');
            title(['electrode# ' num2str(kElectrode)]);
            if strcmp(genvarname(fields{kField}), 'AvSpk')
                xlabel('time');
            else
                xlabel('cluster#')
                nClusInElectrode = length(nq(kElectrode).Clus);
                if nClusInElectrode
                    cluId = nq(kElectrode).Clus;
                    set(gca,'Xtick', cluId );
                    set(gca,'XTickLabel',mat2cell(nq(kElectrode).Clus, ones(1, nClusInElectrode), 1));
                end
            end
        end
        
        set(h,'Name',[ '   ' fields{kField}]);
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
        
        switch dir
            case 'up'
                kField = mod(kField + 1, nFields);
                if mod(kField, nFields) == 0
                    kField = 1;
                end
            case 'down'
                kField = mod(kField - 1, nFields);
                if mod(kField, nFields) == 0
                    kField = nFields;
                end
            case 0
                IS_DONE = 0;
            otherwise
                kField = 1;
        end
        
       
        
                 
       
    end
