function [dgdata, glList] = deglitchAlgo(mseed, xml, corMin, do3c, nofigure, dgdata1c, glList1c)
% this function is the one that actually performs the glitch detection and
% removal; it should be modified as little as possible; it is much better
% to modify the outer functions that call deglitchAlgo.


    %% deglitch parameters
    
    fcutlow = 0.001; % low pass removal frequency
    
    if iscell(mseed)
        chNum = length(mseed);
    else
        chNum = 1;
        tmp = mseed;
        mseed = cell(1,1);
        mseed{1} = tmp;
        clear tmp
    end
    
    fsGlitch = 100; % sampling frequency of the zpk of the glitch


    %% load metadata
    zpk = cell(chNum,1);
    locCode = cell(chNum,1);
    chName = cell(chNum,1);
    for ch=1:chNum
        for i=1:length(mseed{ch})
            if mseed{ch}(i).ChannelFullName ~= mseed{ch}(1).ChannelFullName
                error("In the mseed input cell, each row of the cell must contain only one SEIS channel")
            end
        end
        locCode{ch} = mseed{ch}(1).ChannelFullName(end-5:end-4);
        chName{ch} = mseed{ch}(1).ChannelFullName(end-2:end);
        zpk{ch} = rdXMLmetadata(locCode{ch}, chName{ch}, xml);
    end

    
    %% check all channels are similar sensors
    for ch = 1:chNum
        if ~prod(chName{ch}(1:2) == chName{1}(1:2)) || ~prod(locCode{ch} == locCode{1})
            error("not similar channels... do separate deglitching")
        end
    end

    
    %% divide fcutlow by 4 to improve detection of SP longer glitches iff it is an SP channel
    SPlist = ["65.EH";"65.SH";"66.SH";"67.SH";"68.SH";"65.MH";"66.MH";"67.MH";"65.LH";"70.EH";"70.SH";"71.SH";"72.SH";"73.SH";"70.MH";"71.MH";"72.MH";"73.LH";"70.VH";"71.VH";"72.VH";"73.VH";"70.UH";"71.UH";"72.UH";"73.UH";"70.RH";"71.RH";"72.RH";"73.RH"];

    testSP = false;
    for i=1:length(SPlist)
        if locCode{1} + "." + chName{1}(1:2) == SPlist(i)
            testSP = true;
        end
    end
    if testSP
        fcutlow = fcutlow/4;
    end


    %% load data
    times = cell(chNum,1);
    data = cell(chNum,1);
    fs = cell(chNum,1);
    head = cell(chNum,1);
    lData = zeros(chNum,1);
    if do3c
        for ch=1:chNum
            [times{ch}, ~, fs{ch}] = dataLoader(mseed{ch});
            data{ch} = dgdata1c{ch};
            lData(ch) = length(data{ch});    
        end        
    else
        for ch=1:chNum
            [times{ch}, data{ch}, fs{ch}] = dataLoader(mseed{ch});
            lData(ch) = length(data{ch});    
        end
    end
 
    for ch=1:chNum
        if fs{ch}~=fs{1}
            error("The channels have different frequencies...")
        end
    end
    fs = fs{1};
    fcutlow = fcutlow*fs/2;

    if min(lData)>100*fs
        %% generate glitch
        lGl = cell(chNum,1);
        nShapes = cell(chNum,1);
        glDiff = zeros(chNum,1);
        for ch=1:chNum
            [lGl{ch}, nShapes{ch}, glDiff(ch)] = glitchGenerator(zpk{ch}, fsGlitch, fs);
        end

        %% process data
        processedData = cell(chNum,1);
        filteredPart = cell(chNum, 1);
        for ch=1:chNum
            [processedData{ch},filteredPart{ch}] = dataProcessor(data{ch}, fs, fcutlow);
        end
        if ~nofigure
            lgdNames = [];
            for ch =1:chNum
                lgdNames = [lgdNames; mseed{ch}(1).ChannelFullName(end)];
            end
            figure()
            ax = zeros(1,2);
            ax(1) = subplot(2,1,1);
            for ch =1:chNum
                plot(datetime(times{ch}, 'ConvertFrom', 'datenum'), data{ch})
                hold on
            end
            xlabel("UTC")
            ylabel("amplitude in counts")
            legend(lgdNames)
            if ~do3c
                title(["loaded data from channels " + mseed{1}(1).ChannelFullName(1:end-1) + "? in 1c deglitch"], 'FontSize', 10)
            else
                title(["loaded data from channels " + mseed{1}(1).ChannelFullName(1:end-1) + "? in 3c deglitch"], 'FontSize', 10)
            end
            ax(2) = subplot(2,1,2);
            for ch = 1:chNum
                plot(datetime(times{ch}, 'ConvertFrom', 'datenum'), processedData{ch})
                hold on
            end
            xlabel("UTC")
            ylabel("amplitude in counts")
            legend(lgdNames)
            title(["processed data from channel " + mseed{ch}(1).ChannelFullName(1:end-1) + "?"], 'FontSize', 10)
            linkaxes(ax, 'x')
        end


        %% loop to remove glitches
        if do3c
            loopmax = 1;
        else
            loopmax = 2;
        end
        glTimes = cell(chNum,loopmax);
        dataGlitches = cell(chNum,loopmax);
        glLocKept=cell(chNum,1);
        refGlLocKept=cell(chNum,1);
        glampKept = cell(chNum,1);
        loopData = processedData;

        for loop = 1:loopmax
            % detect the glitches by cross-correlation
            if do3c
                [glLoc, refGlLoc] = glitchDetector3c(loopData, times, lGl, corMin, chNum, nShapes, glList1c, fs, glDiff);
            else
                [glLoc, refGlLoc] = glitchDetector(loopData, lGl, corMin, chNum, nShapes);
            end
            % remove the glitches
            deglitchedData = cell(chNum, 1);
            for ch=1:chNum
                glTimes{ch,loop} = times{ch}(glLoc{ch});   
                [deglitchedData{ch}, dataGlitches{ch,loop}, glamp] = glitchRemover(loopData{ch},...
                    glLoc{ch}, refGlLoc{ch}, lGl{ch});
                glampKept{ch} = [glampKept{ch}; glamp];
            end
            for ch = 1:chNum
                glLoc{ch} = glLoc{ch} + glDiff(ch);
            end
            for ch = 1:chNum
                glLocKept{ch}=[glLocKept{ch} ; glLoc{ch}];
                refGlLocKept{ch}=[refGlLocKept{ch} ; refGlLoc{ch}];
            end
            loopData = deglitchedData;

        end
        
        if ~nofigure
            figure()
            ax = zeros(chNum,1);
            for ch = 1:chNum
                ax(ch) = subplot(chNum, 1, ch);
                plot(datetime(times{ch}, 'ConvertFrom', 'datenum'), processedData{ch}, 'k')
                hold on
                plot(datetime(times{ch}, 'ConvertFrom', 'datenum'), deglitchedData{ch}, 'b')
                for loop = 1:loopmax
                    if loop == 1
                        plot(datetime(times{ch}, 'ConvertFrom', 'datenum'), dataGlitches{ch,loop}, 'r', 'HandleVisibility', 'on')
                    else
                        plot(datetime(times{ch}, 'ConvertFrom', 'datenum'), dataGlitches{ch,loop}, 'r', 'HandleVisibility', 'off')
                    end
                end
                legend('raw data', 'deglitched data', 'glitches')
                ylabel('counts')
                xlabel('UTC')
                if ~do3c
                    title(["processed data from channel " + mseed{ch}(1).ChannelFullName + " in 1c deglitch"], 'FontSize', 10)
                else
                    title(["processed data from channel " + mseed{ch}(1).ChannelFullName + " in 3c deglitch"], 'FontSize', 10)
                end
                    
            end
            linkaxes(ax, 'x')
        end


        %% unfilter the data
        dgdata = cell(chNum, 1);
        glList(chNum) = struct();
        for ch=1:chNum
            dgdata{ch} = deglitchedData{ch} + filteredPart{ch};
            gltimes=times{ch}(glLocKept{ch});
            glList(ch).glLoc = glLocKept{ch};
            glList(ch).gltimes = gltimes;
            glList(ch).glamp = glampKept{ch};
        end
        
    end 
    
end
