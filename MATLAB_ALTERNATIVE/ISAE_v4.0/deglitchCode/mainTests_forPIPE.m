% Code to test the deglitchSEIS algorithm. Put the data in a 'data'
% directory in the same directory as the deglitchSEIS code, or update this
% 'mainTests.m' script to point to the correct directory.


% close all;
% clear all;
% % 
% %CAS SANS DECTIC
% metadata = '/Users/m.drilleau/Documents/InSight/PIPE/Automat_DATA_detigli_3c_160620/DATA/METADATA/XB.ELYSE.StationXML.xml'
% datadir = '/Users/m.drilleau/Documents/InSight/PIPE/Automat_DATA_detigli_3c_160620/DATA//RAW_DATA/SOL170/DETICNOISE/XB.ELYSE.02.BH';
% outputdir='/Users/m.drilleau/Documents/InSight/PIPE/Automat_DATA_detigli_3c_160620/DATA//RAW_DATA/SOL170/DEGLITCH_3c/'
% Compo_dispo='WVU'
% 
% CAS AVEC DETIC 
% metadata = '/Users/m.drilleau/Documents/InSight/PIPE/Automat_DATA_detigli_3c_160620/DATA//METADATA/XB.ELYSE.StationXML.xml'
% datadir = '/Users/m.drilleau/Documents/InSight/PIPE/Automat_DATA_detigli_3c_160620/DATA/RAW_DATA/SOL212/DETICNOISE/XB.ELYSE.02.BH';
% outputdir='/Users/m.drilleau/Documents/InSight/PIPE/Automat_DATA_detigli_3c_160620/DATA//RAW_DATA/SOL212/DEGLITCH_3c/'
% Compo_dispo='WVU'

%% load parameters

Compo_dispo=sort(Compo_dispo);
chNum = length(Compo_dispo);

disp('Deglitch the following files:')
list = ls([datadir '*']);
len_list = length(list);
inc = (len_list - chNum)/chNum;

nbSOL = outputdir(end-18:end-13);

i1 = 1;
i2 = inc;
for ch=1:chNum;
    filenames{ch} = list(i1:i2);  
    disp(filenames{ch})
    OutPutName{ch} = [outputdir, datadir(end-13:end), Compo_dispo(ch), '/', datadir(end-13:end), Compo_dispo(ch), '.', nbSOL];
    saveName{ch} = [outputdir, datadir(end-13:end), Compo_dispo(ch), '/saveDir/'];
    % to save the glitch list; needs a directory saveDir
    if ~isfolder(saveName{ch})
       mkdir(saveName{ch})
    end
    i1 = i1 + 1 + inc;
    i2 = i1 + inc-1;
end
disp(' ')

do3c = 1; % 1 if you want a more extensive detection using the 3 components at the same time
nofigure = 1; %0 if you want some figures
doSave = 1; %to not save the mseed

corMin = 0.88; %0.9; % threshold for the cross-correlation; empirically chosen but can be modified
corMin3c = 0.5; %0.5; % threshold for 3c cross-correlation; empirically chosen but can be modified

xml = xmlread(metadata);


mseed = cell(chNum,1);
for ch = 1:chNum
    mseed{ch} = rdmseed(filenames{ch});
end

% %% check load
% figure();
% for ch = 1:chNum
%     plot(datetime(cat(1,mseed{ch}(:).t), 'ConvertFrom', 'datenum'), cat(1,mseed{ch}(:).d))
%     hold on
% end
% xlabel("UTC")
% ylabel("counts")
% legend("02.BHU","02.BHV","02.BHW")
% % legend("02.BHU","02.BHV")
% title("load check")


%% main code
tic()
if ~do3c
    [dgdata, glList] = deglitchSEIS(mseed, xml, corMin, do3c, [], nofigure);
else
    [dgdata, glList, glList1c] = deglitchSEIS(mseed, xml, corMin, do3c, corMin3c, nofigure);
end
toc()


% %% check results
% figure();
% colors = ["b", "r", "k"];
% % colors = ["b", "r"];
% for ch = 1:chNum
%     plot(datetime(cat(1,mseed{ch}(:).t), 'ConvertFrom', 'datenum'), cat(1,mseed{ch}(:).d), 'Color', colors(ch), 'LineStyle', '-')
%     hold on
%     plot(datetime(cat(1,mseed{ch}(:).t), 'ConvertFrom', 'datenum'), dgdata{ch}, 'Color', colors(ch), 'LineStyle', '--')
% end
% xlabel("UTC")
% ylabel("counts")
% legend("02.BHU raw","02.BHU deglitched","02.BHV raw","02.BHV deglitched", "02.BHW raw", "02.BHW deglitched")
% % legend("02.BHU raw","02.BHU deglitched","02.BHV raw","02.BHV deglitched")
% title("results check")


%% save results in .mseed format
if doSave
    for ch = 1:chNum
        mkmseed(OutPutName{ch}, dgdata{ch}, cat(1,mseed{ch}(:).t), mseed{ch}(1).SampleRate,'onefile')
    end
end

%% provide glitch list
if iscell(mseed)
    chNum = length(mseed);
else
    chNum = 1;
    tmp = mseed;
    mseed = cell(1,1);
    mseed{1} = tmp;
    clear tmp
end
if do3c
    glStruct = struct('channel', [], 'glTimesUTC', [], 'gldelay', [], 'glAmpCounts', [], 'y0', [], 'slope', []);
    for ch = 1:chNum
        
        glitchList = fopen([saveName{ch}, 'glitchList.txt'], 'w');
        fprintf(glitchList,'%-15s  %23s  %9s  %18s  %-15s   %10s\n','channel','YYYY-MM-ddTHH:mm:ss.SSS', 'delay (s)', 'amplitude (counts)', 'y0 (counts)', 'slope (counts/s)');
        
        for i = 1:length(glList1c(ch).gltimes)
            glStruct.channel = [glStruct.channel; mseed{ch}(1).ChannelFullName];
        end
        glStruct.glTimesUTC = [glStruct.glTimesUTC; glList1c(ch).gltimes];
        glStruct.glAmpCounts = [glStruct.glAmpCounts; glList1c(ch).glamp];
        glStruct.gldelay = [glStruct.gldelay; glList1c(ch).gldelay];
        glStruct.y0 = [glStruct.y0; glList1c(ch).y0];
        glStruct.slope = [glStruct.slope; glList1c(ch).slope];
        
        for i = 1:length(glStruct.glTimesUTC);
            fprintf(glitchList, ['%15s  %23s  %-.3e  %+.9e', blanks(3), ' %+.9e   %+.9e\n'], glStruct.channel(i,:), ...
                    datetime(glStruct.glTimesUTC(i), 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd''T''HH:mm:ss.S'), ...
                    glStruct.gldelay(i,:), glStruct.glAmpCounts(i), glStruct.y0(i),glStruct.slope(i) );                            
        end   
        
        glStruct = struct('channel', [], 'glTimesUTC', [], 'gldelay', [], 'glAmpCounts', [], 'y0', [], 'slope', []);
        
    end
    glStruct3c = struct('channel', [], 'glTimesUTC', [], 'gldelay', [], 'glAmpCounts', [], 'y0', [], 'slope', []);
    if ~isempty(glList)
        for ch = 1:chNum
            
            glitchList = fopen([saveName{ch}, 'glitchList3c.txt'], 'w');
            fprintf(glitchList,'%-15s  %23s  %9s  %18s  %-15s   %10s\n','channel','YYYY-MM-ddTHH:mm:ss.SSS', 'delay (s)', 'amplitude (counts)', 'y0 (counts)', 'slope (counts/s)');
            
            for i = 1:length(glList(ch).gltimes)
                glStruct3c.channel = [glStruct3c.channel; mseed{ch}(1).ChannelFullName];
            end
            glStruct3c.glTimesUTC = [glStruct3c.glTimesUTC; glList(ch).gltimes];
            glStruct3c.glAmpCounts = [glStruct3c.glAmpCounts; glList(ch).glamp];
            glStruct3c.gldelay = [glStruct3c.gldelay; glList(ch).gldelay];
            glStruct3c.y0 = [glStruct3c.y0; glList(ch).y0];
            glStruct3c.slope = [glStruct3c.slope; glList(ch).slope];
            
            for i = 1:length(glStruct3c.glTimesUTC);
                fprintf(glitchList, ['%15s  %23s  %-.3e  %+.9e', blanks(3), ' %+.9e   %+.9e\n'], glStruct3c.channel(i,:), ...
                    datetime(glStruct3c.glTimesUTC(i), 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd''T''HH:mm:ss.S'), ...
                    glStruct3c.gldelay(i,:), glStruct3c.glAmpCounts(i), glStruct3c.y0(i),glStruct3c.slope(i) );
           end
           glStruct3c = struct('channel', [], 'glTimesUTC', [], 'gldelay', [], 'glAmpCounts', [], 'y0', [], 'slope', []);
        end
    end
%     save([saveName, '.mat'], '-struct', 'glStruct')
%     saveName3c = [saveName '3c.mat'];
%     save(saveName3c, '-struct', 'glStruct3c')
else
%     glStruct = struct('channel', [], 'glTimesUTC', [], 'glAmpCounts', []);
    glStruct = struct('channel', [], 'glTimesUTC', [], 'gldelay', [], 'glAmpCounts', [], 'y0', [], 'slope', []);
    for ch = 1:chNum
        for i = 1:length(glList(ch).gltimes)
            glStruct.channel = [glStruct.channel; mseed{ch}(1).ChannelFullName];
        end
        glStruct.glTimesUTC = [glStruct.glTimesUTC; glList(ch).gltimes];
        glStruct.glAmpCounts = [glStruct.glAmpCounts; glList(ch).glamp];
        glStruct.gldelay = [glStruct.gldelay; glList(ch).gldelay];
        glStruct.y0 = [glStruct.y0; glList(ch).y0];
        glStruct.slope = [glStruct.slope; glList(ch).slope];
        
    end
    save([saveName, '.mat'], '-struct', 'glStruct')
end
