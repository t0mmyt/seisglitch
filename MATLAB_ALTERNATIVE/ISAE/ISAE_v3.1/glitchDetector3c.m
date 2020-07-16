function [glLoc3c, refGlLoc3c] = glitchDetector3c(data, times, lGl, corMin3c, chNumber, nShapes, gllisttmp, fs, glDiff)
    
    delayMax = 2/fs/86400;
    glitch = cell(chNumber,1);
    lenGl = cell(chNumber,1);
    lData = cell(chNumber,1);
    for ch=1:chNumber
        glitch{ch} = lGl{ch}(:,1);
        lenGl{ch} = length(glitch{ch});
        lData{ch} = length(data{ch});
    end

    %% compute correlation coefficients
    corcoef = cell(chNumber,1);
    glitchNorm = cell(chNumber,1);
    for ch=1:chNumber
        glitchNorm{ch}=sqrt(sum(glitch{ch}.*glitch{ch}));
    end
    for ch=1:chNumber
        [corfun,lags]=xcorr(data{ch},glitch{ch},'none');
        i0=find(lags==0);
        x2=data{ch}.*data{ch};
        xnorm2=conv(x2,ones(lenGl{ch},1),'full');
        dataNorm=sqrt(xnorm2(lenGl{ch}:end));
        corcoef{ch}=corfun(i0:end)./(glitchNorm{ch}*dataNorm);
    end

    
    glLoc3c = cell(chNumber, 1);
    %% find peaks of correlation that are greater than corMin3c and linked to other glitches
    for ch=1:chNumber
        pks = islocalmax(abs(corcoef{ch}), 'FlatSelection', 'all');
        pks = find(pks==1);
        indcor = find(abs(corcoef{ch})>corMin3c);
        indcor = intersect(indcor, pks);
        indcor = indcor + glDiff(ch);
    
        otherCh = [1:ch-1 ch+1:chNumber];
        otherTimes = vertcat(gllisttmp(otherCh).gltimes);
        otherTimes = sort(otherTimes);
        timeCor = times{ch}(indcor);
        [timeCor, idxsorttimeCor] = sort(timeCor);
        indcor = indcor(idxsorttimeCor);
        indcor2 = zeros(length(indcor),1);
        idxIndcor2 = 0;
        idxAfter = 1;
        timeAfter = otherTimes(idxAfter);
        for i = 1:length(timeCor)
            while timeAfter < timeCor(i) && idxAfter < length(otherTimes)
                idxAfter = idxAfter+1;
                timeAfter = otherTimes(idxAfter);
            end
            idxBefore = max(idxAfter - 1, 1);
            timeBefore = otherTimes(idxBefore);
            delay = min([abs(timeBefore-timeCor(i)), abs(timeAfter-timeCor(i))]);
            if delay <= delayMax
                idxIndcor2 = idxIndcor2 + 1;
                indcor2(idxIndcor2) = indcor(i);
            end
        end
        indcor = indcor2(1:idxIndcor2);
        glLoc3c{ch} = sort(indcor-glDiff(ch));
    end
        

    %% refine glitch locations
    pbCor=cell(chNumber,1);
    nPos = cell(chNumber, 1);
    refGlLoc3c = cell(chNumber,1);
    refCorGl = cell(chNumber,1);
    for ch=1:chNumber
        pbCor{ch}=0;
        nPos{ch} = length(glLoc3c{ch});
        refGlLoc3c{ch} = zeros(nPos{ch},1);
        nGl = 0;
        idxBestGl = zeros(nPos{ch},1);
        refCor = zeros(length(glLoc3c{ch}),1);
        x2=data{ch}.*data{ch};
        xnorm2=conv(x2,ones(lenGl{ch},1),'full');
        dataNorm=sqrt(xnorm2(lenGl{ch}:end));
        for detectedGlitch=(glLoc3c{ch})'
            nGl = nGl + 1;
            corCoef = cell(nShapes{ch},1);
            lags = cell(nShapes{ch},1);
            corMax = zeros(nShapes{ch},1);
            indMax = zeros(nShapes{ch},1);
            for i=1:nShapes{ch}
                [corCoef{i}, lags{i}] = xcorr(...
                    data{ch}(detectedGlitch-10:detectedGlitch+lenGl{ch}+9),...
                    lGl{ch}(:,i));
                i0=find(lags{i}==0);
                corCoef{i} = corCoef{i}(i0:end);
                glitchNorm{ch}=sqrt(lGl{ch}(:,i)'*lGl{ch}(:,i));
                corCoef{i} = (corCoef{i}/glitchNorm{ch})./...
                    dataNorm(detectedGlitch-10:detectedGlitch+lenGl{ch}+9);
                [corMax(i), indMax(i)] = max(abs(corCoef{i}));
            end
            % best position for oversampled glitches
            [refCor(nGl),idxBestGl(nGl)] = max(corMax);
            glLoc3c{ch}(nGl) = glLoc3c{ch}(nGl) + indMax(idxBestGl(nGl))-11;
        end
        refCorGl{ch}=refCor;
        refGlLoc3c{ch}=idxBestGl;
    end

end