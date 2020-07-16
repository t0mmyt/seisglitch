function [processedData, filteredPart] = dataProcessor(data, fs, fcutlow, filtertype)

if filtertype == "cut_low"
    %% Filtering by subtracting low pass

    wlh=fcutlow/fs;
    norderl=2;
    [bb,ab] = butter(norderl,wlh,'low');

    filteredPart = filtfilt(bb,ab,data(:));
    processedData=data-filteredPart;
    
elseif filtertype == "cut_high"
     %% Filtering by subtracting high pass
    if fcutlow < fs
        wlh=fcutlow/fs;
        norderl=2;
        [bb,ab] = butter(norderl,wlh,'high');
        filteredPart = filtfilt(bb,ab,data(:));
    else
        filteredPart = zeros(length(data), 1);
    end
    processedData=data-filteredPart;   

end