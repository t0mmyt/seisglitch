function [dgdata, glList, glList1c, dgdata1c] = deglitchSEIS(mseed, xml, corMin, do3c, corMin3c, nofigure)
% Main code for deglitching SEIS data
%
% This algorithm detects glitches with a cross-correlation and then removes
% them from the data by a least-squares fit. It either works in a 1c
% process(for 1 component process) or a 3c process. 1c is the faster one, 
% 3c uses in addition the results from a 1c to compare each component 
% together and find more glitches that way (doing a new cross-correlation
% with a lower threshold). When selecting 3c, it automatically does the
% necessary previous 1c.
%
% The inputs handled by this function are VBB VEL channels and SP channels.
%
% OUTPUTS: (note that the 2 first outputs are for 1c deglitching if do3c is
% not selected, and for the 3c deglitching otherwise; the two following
% outputs are only used for 3c deglitching, to give the results of the
% necessary 1c deglitching that is also conducted)
% - dgdata: deglitched data; it is a Matlab cell. Each row of the cell
% corresponds to a channel (e.g. dgdata{1} : BHU, dgdata{2} : BHV...). In
% each row, there are only the deglitched data time series. Use the 'mseed'
% input and the 'mkmseed.m' code to generate the .mseed corresponding to
% the deglitched data. In case of 3c, this includes all the successive
% deglitching that has been done in 1c and then in 3c.
% - glList: list of detected glitches. It is a (n by 1) structure, n being
% the number of channels, containing fields gltimes (glitch time in UTC
% datenum), glamp (glitch amplitude in counts), glLoc (idx of the glitch in
% the data). It is advised to save it as a .mat (see 'mainTests.m'
% provided), and then turn it into a text file (see 'glitchList2txt.m').
% For 3c deglitching, it only provides the ADDITIONAL glitches that are
% provided from the 3c; in other words, the sum of all glitches is
% gList1c + glList3c (no overlap).
%
% INPUTS: 
% - mseed, a cell of mseeds already read by rdmseed; or only one
% mseed read by rdmseed. Each cell row is a different component, at similar
% times as the other components (especially for the 3c deglitching).
% - xml, a metadata file already opened as a java DOM node
% - corMin, the cross-correlation coefficient threshold to detect as a
% glitch; empirically 0.9 is a good choice.
% - do3c: 1 = yes, 0 = no, only 1c.
% - corMin3c, the cross-correlation coefficient threshold for the
% deglitch3c process; should be lower than the 1c.
% - nofigure: 1 = yes, no figure display; 0 = display figures
%
% contact/author: baptiste.pinot@isae-supaero.fr

%    arguments
%        mseed
%        xml
%        corMin
%        do3c (1,1) logical
%        corMin3c
%        nofigure (1,1) logical
%    end
    
    oneChannel = 0;
    
    if iscell(mseed)
        if length(mseed) == 1
            oneChannel = 1;
        end
    end    
    
    if ~do3c
        % do not do 3c simultaneous deglitching
        % on each channel (separately), detects glitches that have a higher
        % correlation coefficient than corMin and then removes them
        [dgdata, glList] = deglitchAlgo(mseed, xml, corMin, 0, nofigure);
        dgdata1c = [];
        glList1c = [];
        
    else
        % do 3c simultaneous deglitching
        
        % first find glitches with 1c process
        % same as 1c process
        [dgdata1c, glList1c] = deglitchAlgo(mseed, xml, corMin,0, nofigure);
        
        goOn = 0;
        for ch = 1:length(glList1c)
            if ~isempty(glList1c(ch).glLoc)
                goOn = 1;
            end
        end
        if goOn && ~oneChannel
            % then add the glitches found with 3c process
            % To make it work, one therefore always needs corMin3c < corMin
            [dgdata, glList] = deglitchAlgo(mseed, xml, corMin3c, do3c, nofigure, dgdata1c, glList1c);
        elseif ~goOn
            disp("Unable to perform deglitch3c: no glitches detected in first pass")
            dgdata = dgdata1c;
            glList = [];
        elseif oneChannel
            disp("Unable to perform deglitch3c: only one input component")
            dgdata = dgdata1c;
            glList = [];
        end
        
    end
    
end