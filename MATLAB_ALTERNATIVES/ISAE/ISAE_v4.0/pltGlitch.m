
% plot the data at the time of the glitches in the glitch list previously
% found


%% data loading

directory = '/home/deos/b.pinot/Documents/insight/dataRepository/02BHoneWeek/';

filenames{1} = [directory, 'XB.ELYSE.02.BHU.mseed'];
filenames{2} = [directory, 'XB.ELYSE.02.BHV.mseed'];
filenames{3} = [directory, 'XB.ELYSE.02.BHW.mseed'];


mseed = cell(3,1);
for ch = 1:3
    mseed{ch} = rdmseed(filenames{ch});
end
fs = mseed{1}(1).SampleRate;

times = cell(3,1);
data = cell(3,1);

for ch = 1:3
    times{ch} = cat(1,mseed{ch}(:).t);
    [times{ch},I] = unique(times{ch},'sorted');
    times{ch} = times{ch}(I);
    data{ch} = cat(1,mseed{ch}(:).d);
    data{ch} = data{ch}(I);
end


%% glitch reading

saveFile = '/home/deos/b.pinot/Documents/insight/glitches/saveDirTest/glitchList1c.txt';
[tSt, delay, amp, y0, slope] = rdList(saveFile);

for ch = 1:3
    tSt{ch} = datenum(tSt{ch});
end


%% glitch plot

tBefore = 10/86400;
tAfter = 45/86400;
sensors = 'UVW';

for ch = 1:3
    
    figure()
        
    for i = 1:length(tSt{ch})  
        idxB = find(times{ch} >= tSt{ch}(i) - tBefore);
        idxA = find(times{ch} <= tSt{ch}(i) + tAfter);
        j0 = find(times{ch} >= tSt{ch}(i), 1, 'first');
        z0 = y0{ch}(i);
        idx = intersect(idxA, idxB);
        ytrend = zeros(length(idx),1);
        for j = 1:length(idx)
            ytrend(j) = slope{ch}(i)*(idx(j)-j0)/fs + z0;
        end
        timeVect = (times{ch}(idx)-tSt{ch}(i))*86400 - delay{ch}(i);
        plot(timeVect, (data{ch}(idx)-ytrend)/amp{ch}(i), 'Color', '#C0C0C0')
        hold on
    end
    
    xlabel("time (s)")
    ylabel("normalised")
    title("Glitches found on 02.BH" + sensors(ch) + " for one week, june-july 2019")
    ylim([-2, 4])
    fig = gcf();
    saveas(fig, ['glitchesOneWeekSum_02.BH', sensors(ch), '.png'])
    
end
    
    