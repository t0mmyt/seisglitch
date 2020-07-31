% Code to test the deglitchSEIS algorithm. Put the data in a 'data'
% directory in the same directory as the deglitchSEIS code, or update this
% 'mainTests.m' script to point to the correct directory.


%% load parameters

directory = './data/';
metadata = [directory, 'XB.ELYSE.metadata.xml'];
xml = xmlread(metadata);

filenames{1} = [directory, 'XB.ELYSE.02.BHU.mseed'];
filenames{2} = [directory, 'XB.ELYSE.02.BHV.mseed'];
filenames{3} = [directory, 'XB.ELYSE.02.BHW.mseed'];

mseed = cell(3,1);
for ch = 1:3
    mseed{ch} = rdmseed(filenames{ch});
end

do3c = 1; % 1 if you want a more extensive detection using the 3 components at the same time
nofigure = 0; %0 if you want some figures

corMin = 0.9; % threshold for the cross-correlation; empirically chosen but can be modified
corMin3c = 0.5; % threshold for 3c cross-correlation; empirically chosen but can be modified

saveName = '/saveDir/saveGlStruct'; % to save the glitch list; needs a directory saveDir


%% check load
% figure();
% for ch = 1:3
%     plot(datetime(cat(1,mseed{ch}(:).t), 'ConvertFrom', 'datenum'), cat(1,mseed{ch}(:).d))
%     hold on
% end
% xlabel("UTC")
% ylabel("counts")
% legend("02.BHU","02.BHV","02.BHW")
% title("load check")


%% main code
tic()
if ~do3c
    [dgdata, glList] = deglitchSEIS(mseed, xml, corMin, do3c, [], nofigure);
else
    [dgdata, glList, glList1c] = deglitchSEIS(mseed, xml, corMin, do3c, corMin3c, nofigure);
end
toc()


%% check results
% figure();
% colors = ["b", "r", "k"];
% for ch = 1:3
%     plot(datetime(cat(1,mseed{ch}(:).t), 'ConvertFrom', 'datenum'), cat(1,mseed{ch}(:).d), 'Color', colors(ch), 'LineStyle', '-')
%     hold on
%     plot(datetime(cat(1,mseed{ch}(:).t), 'ConvertFrom', 'datenum'), dgdata{ch}, 'Color', colors(ch), 'LineStyle', '--')
% end
% xlabel("UTC")
% ylabel("counts")
% legend("02.BHU raw","02.BHU deglitched","02.BHV raw","02.BHV deglitched", "02.BHW raw", "02.BHW deglitched")
% title("results check")


%% save results in .mseed format
for ch = 1:3
    mkmseed(replace(mseed{ch}(1).ChannelFullName, ':', '.'), dgdata{ch}, cat(1,mseed{ch}(:).t), mseed{ch}(1).SampleRate)
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
    glStruct = struct('channel', [], 'glTimesUTC', [], 'glAmpCounts', []);
    for ch = 1:chNum
        for i = 1:length(glList1c(ch).gltimes)
            glStruct.channel = [glStruct.channel; mseed{ch}(1).ChannelFullName];
        end
        glStruct.glTimesUTC = [glStruct.glTimesUTC; glList1c(ch).gltimes];
        glStruct.glAmpCounts = [glStruct.glAmpCounts; glList1c(ch).glamp];
    end
    glStruct3c = struct('channel', [], 'glTimesUTC', [], 'glAmpCounts', []);
    for ch = 1:chNum
        for i = 1:length(glList(ch).gltimes)
            glStruct3c.channel = [glStruct3c.channel; mseed{ch}(1).ChannelFullName];
        end
        glStruct3c.glTimesUTC = [glStruct3c.glTimesUTC; glList(ch).gltimes];
        glStruct3c.glAmpCounts = [glStruct3c.glAmpCounts; glList(ch).glamp];
    end
    save([saveName, '.mat'], '-struct', 'glStruct')
    saveName3c = [saveName '3c.mat'];
    save(saveName3c, '-struct', 'glStruct3c')
else
    glStruct = struct('channel', [], 'glTimesUTC', [], 'glAmpCounts', []);
    for ch = 1:chNum
        for i = 1:length(glList(ch).gltimes)
            glStruct.channel = [glStruct.channel; mseed{ch}(1).ChannelFullName];
        end
        glStruct.glTimesUTC = [glStruct.glTimesUTC; glList(ch).gltimes];
        glStruct.glAmpCounts = [glStruct.glAmpCounts; glList(ch).glamp];
    end
    save([saveName, '.mat'], '-struct', 'glStruct')
end
