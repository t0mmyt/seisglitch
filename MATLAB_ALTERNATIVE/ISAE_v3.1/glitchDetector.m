function [glLoc, refGlLoc] = glitchDetector(data, lGl, corMin, chNumber, nShapes)

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


    %% search for best correlation coefficients
    glLoc = cell(chNumber, 1);
    listCorGlitch = cell(chNumber, 1);
    nGlSel = cell(chNumber,1);

    % first search in data
    for ch=1:chNumber    
        indcor=find(abs(corcoef{ch})>corMin);

        nindcor=length(indcor);
        if nindcor > 0
            glitchBeginning=max([lenGl{ch}+1 , indcor(1)]);
            glitchEnd=max([lenGl{ch}+2 , indcor(1)+1]);
        else
            glitchBeginning=lenGl{ch}+1 ;
            glitchEnd=lenGl{ch}+2 ;
        end

        j=0;

        % initialisation of the lists
        corGlitch = zeros(nindcor,1);
        glitchPos = zeros(nindcor,1);
        glitchShift = zeros(nindcor,1);
        glitchSign = zeros(nindcor,1);
        glitchShape = zeros(nindcor, lenGl{ch});
        glitchAmp = zeros(nindcor,1);

        for corNumber=1:nindcor
            if (corNumber<nindcor)&&(indcor(corNumber+1)<lData{ch}-2*lenGl{ch})&& ...
                    (indcor(corNumber)>lenGl{ch})
                if (indcor(corNumber+1)-indcor(corNumber))==1
                    glitchEnd=indcor(corNumber+1);
                else
                    [cor,imax]=max(abs(corcoef{ch}(glitchBeginning:glitchEnd)));
                    glitchPosTemp=glitchBeginning+imax-1;
                    j=j+1;
                    corGlitch(j)=cor;
                    glitchPos(j)=glitchPosTemp;
                    glitchShift(j)=mean(data{ch}(glitchPosTemp:glitchPosTemp+3));
                    glitchSign(j)=sign(corcoef{ch}(glitchPosTemp));
                    glitchShape(j,:)=glitchSign(j)*data{ch}(glitchPosTemp:glitchPosTemp+lenGl{ch}-1);
                    glitchAmp(j)=max(glitchSign(j)*(data{ch}(glitchPosTemp:glitchPosTemp+lenGl{ch}-1)-glitchShift(j)));
                    corcoef{ch}(glitchBeginning-5:glitchEnd+5)=0.0;
                    glitchBeginning=indcor(corNumber+1);
                    glitchEnd=indcor(corNumber+1);
                end
            end
        end

        nGlSel{ch} = j; % number of glitches selected

        % cut the lists
        corGlitch = corGlitch(1:nGlSel{ch},:);
        glitchPos = glitchPos(1:nGlSel{ch},:);
        glitchShift = glitchShift(1:nGlSel{ch},:);
        glitchSign = glitchSign(1:nGlSel{ch},:);
        glitchShape = glitchShape(1:nGlSel{ch},:);
        glitchAmp = glitchAmp(1:nGlSel{ch},:);

        listCorGlitch{ch} = corGlitch;
        glLoc{ch} = glitchPos;
    end


    %% refine glitch locations
    pbCor=cell(chNumber,1);
    nPos = cell(chNumber, 1);
    refGlLoc = cell(chNumber,1);
    refCorGl = cell(chNumber,1);
    for ch=1:chNumber
        pbCor{ch}=0;
        nPos{ch} = length(glLoc{ch});
        refGlLoc{ch} = zeros(nPos{ch},1);
        nGl = 0;
        idxBestGl = zeros(nPos{ch},1);
        refCor = zeros(length(glLoc{ch}),1);
        x2=data{ch}.*data{ch};
        xnorm2=conv(x2,ones(lenGl{ch},1),'full');
        dataNorm=sqrt(xnorm2(lenGl{ch}:end));	
        for detectedGlitch=(glLoc{ch})'
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
            glLoc{ch}(nGl) = glLoc{ch}(nGl) + indMax(idxBestGl(nGl))-11;
        end
        refCorGl{ch}=refCor;
        refGlLoc{ch}=idxBestGl;
    end
    

end