function [tSt, delay, amp, y0, slope] = rdList(saveFile)
% read the glitch list (the text list from 'glitchList2txt.m')

glitchList = fopen(saveFile, 'r');
tline = fgetl(glitchList);
disp(tline)
tline = fgetl(glitchList);
tSt = cell(3,1);
amp = cell(3,1);
delay = cell(3,1);
y0 = cell(3,1);
slope = cell(3,1);
while ischar(tline)
    C = textscan(tline, ['%15s  %23s  %.3f  %.9f', blanks(3), ' %.9f   %.9f\n']);
    if prod(C{1}{1}=='XB:ELYSE:02:BHU')
        ch = 1;
    elseif prod(C{1}{1}=='XB:ELYSE:02:BHV')
        ch = 2;
    elseif prod(C{1}{1}=='XB:ELYSE:02:BHW')
        ch = 3;
    end
    tSt{ch} = [tSt{ch}; datetime(C{2}{1}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS')];
    delay{ch} = [delay{ch}; C{3}];
    amp{ch} = [amp{ch}; C{4}];
    y0{ch} = [y0{ch}; C{5}];
    slope{ch} = [slope{ch}; C{6}];
    tline = fgetl(glitchList);
end

for ch = 1:3
    disp(['number of glitches for channel ', int2str(ch), ':'])
    disp(length(tSt{ch}))
end