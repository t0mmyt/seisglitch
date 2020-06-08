function [processedData, filteredPart] = dataProcessor(data, fs, fcutlow)
%% Filtering by subtracting low pass

flmax=fs/2;
wlh=fcutlow/flmax;
norderl=2;
[bb,ab] = butter(norderl,wlh,'low');

filteredPart = filtfilt(bb,ab,data(:));
processedData=data-filteredPart;

end