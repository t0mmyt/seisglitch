% turns the 'glitch list' structure from deglitch into a text file

%% for 1c glitches
tic()
saveDir = "./saveDir/";
listCh = dir(saveDir+"*.mat");

glitchList = fopen([saveDir, 'glitchList1c.txt'], 'w');

fprintf(glitchList,'%17s   %21s   %18s\n','channel','YYYY-MM-ddTHH:mm:ss.S', 'amplitude (counts)');

for chnum = 1:length(listCh)
    
    filename = listCh(chnum).name;
    if ~contains(filename, "3c") % to select only the glitches found by a 1c deglitching - depends on how the .mat files are named
        savedMat = load(saveDir+filename);
        for i = 1:length(savedMat.glTimesUTC)
            if savedMat.glAmpCounts(i)>=0
                fprintf(glitchList, '%17s   %21s   %1s%-.9e\n', savedMat.channel(i,:), ...
                    datetime(savedMat.glTimesUTC(i), 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd''T''HH:mm:ss.S'), ...
                    '+', abs(savedMat.glAmpCounts(i)));
            else
                fprintf(glitchList, '%17s   %21s   %1s%-.9e\n', savedMat.channel(i,:), ...
                    datetime(savedMat.glTimesUTC(i), 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd''T''HH:mm:ss.S'), ...
                    '-', abs(savedMat.glAmpCounts(i)));
            end

        end
    end

end

fclose(glitchList);
toc()

%% for 3c deglitching
tic()
saveDir = "./saveDir/";
listCh = dir(saveDir+"*.mat");

glitchList = fopen([saveDir, 'glitchList3c.txt'], 'w');

fprintf(glitchList,'%17s   %21s   %18s\n','channel','YYYY-MM-ddTHH:mm:ss.S', 'amplitude (counts)');

for chnum = 1:length(listCh)
    
    filename = listCh(chnum).name;
    if contains(filename, "3c") % to select only the glitches found by a 3c deglitching - depends on how the .mat files are named
        savedMat = load(saveDir+filename);
        for i = 1:length(savedMat.glTimesUTC)
            if savedMat.glAmpCounts(i)>=0
                fprintf(glitchList, '%17s   %21s   %1s%-.9e\n', savedMat.channel(i,:), ...
                    datetime(savedMat.glTimesUTC(i), 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd''T''HH:mm:ss.S'), ...
                    '+', abs(savedMat.glAmpCounts(i)));
            else
                fprintf(glitchList, '%17s   %21s   %1s%-.9e\n', savedMat.channel(i,:), ...
                    datetime(savedMat.glTimesUTC(i), 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd''T''HH:mm:ss.S'), ...
                    '-', abs(savedMat.glAmpCounts(i)));
            end

        end
    end

end

fclose(glitchList);
toc()

          
          
          
          
          
          
