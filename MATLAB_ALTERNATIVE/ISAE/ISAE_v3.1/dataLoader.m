function [times, data, fs] = dataLoader(mseed)

data = cat(1,mseed(:).d);
times = cat(1,mseed(:).t);
fs = mseed(1).SampleRate;

end