function zpk = rdXMLmetadata(locCode, chName, xml)
% returns a zpk structure from a dataless in file metadata
% the zpk structure contains fiels 'locCode', 'chName', 'startDate',
% 'endDate', and 'zeros' and 'poles' if applicable, as well as 'gain' and
% 'volts2counts' that give the stage 1 gain and stage 2 (digitalisation) gain,
% respectively.
%
% locCode is the code number while chName is the name of the channel.
% For example, 03.BHU has locCode = "03" and chName = "BHU".
% It is important to put double quotes around both inputs.
% Note that these inputs can be multiple, in which case you have to write
% them line by line and they must have the same number of lines;
% for example, ["00"; "00"] as locCode and ["BHU"; "BHV"] is valid.



    %% opens selected channels (00.VMU, O3.VKI etc)
    % returns the list of item number in station with the correct channels
    
    chNum = length(locCode(:,1));
    
    child = xml.getFirstChild();
    % child.getNodeName()

    listCh = cell(chNum,1);
    for ch = 1:chNum
        listCh{ch} = [];
        children = child.getChildNodes();
        for i = 0:children.getLength()-1
            if string(children.item(i).getNodeName()) == "Network"
                network = children.item(i).getChildNodes();
                for k = 0:network.getLength()-1
                    if string(network.item(k).getNodeName()) == "Station"
                        station = network.item(k).getChildNodes();
                        for j = 0:station.getLength()-1
                            if string(station.item(j).getNodeName()) == "Channel"
                                if string(station.item(j).getAttribute("code")) == chName(ch,:) && string(station.item(j).getAttribute("locationCode")) == locCode(ch,:)
                                    listCh{ch}(end+1) = j;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    %% zpk structure generation
    totLength = 0;
    for ch = 1:chNum
        totLength = totLength + length(listCh{ch});
    end
    zpk(totLength) = struct();


    %% gets end time and start time of each trace
    
    cntrZpk = 0;
    for ch = 1:chNum
        startDate = NaT(length(listCh{ch}), 1);
        endDate = NaT(length(listCh{ch}), 1);
        for j = 1:length(listCh{ch})
            cntrZpk = cntrZpk + 1;
            zpk(cntrZpk).locCode = string(station.item(listCh{ch}(j)).getAttribute("locationCode"));
            zpk(cntrZpk).channel = string(station.item(listCh{ch}(j)).getAttribute("code"));
            attr = station.item(listCh{ch}(j)).getAttributes();
            for i=0:attr.getLength()-1
                if string(attr.item(i).getNodeName()) == "startDate"
                    content = char(attr.item(i).getTextContent());
                    if length(content) == 19
                        content = [content,'.000'];
                    end
                    startDate(j) = datetime(content, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
                    zpk(cntrZpk).startDate = startDate(j);
                end
                if string(attr.item(i).getNodeName()) == "endDate"
                    content = char(attr.item(i).getTextContent());
                    if length(content) == 19
                        content = [content,'.000'];
                    end
                    endDate(j) = datetime(content, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
                    zpk(cntrZpk).endDate = endDate(j);
                end
            end
        end
    end
    

    %% finds instrument response
    cntrZpk = 0;
    for ch = 1:chNum
        for j = 1:length(listCh{ch})
            transferFunction = false;
            channel = station.item(listCh{ch}(j)).getChildNodes();
            for i=0:channel.getLength()-1
                if string(channel.item(i).getNodeName()) == "Response"
                    response = channel.item(i).getChildNodes();
                    clear stageTF stageDigit
                    for k =0:response.getLength()-1
                        if string(response.item(k).getNodeName()) == "Stage"
                            attr = response.item(k).getAttributes();
                            for attrNum = 0:attr.getLength()-1
                                if string(attr.item(attrNum).getNodeName()) == "number"
                                    if attr.item(attrNum).getTextContent == '1'
                                        stageTF = response.item(k).getChildNodes(); % stage transfer function
                                    end
                                    if attr.item(attrNum).getTextContent == '2'
                                        stageDigit = response.item(k).getChildNodes(); % stage digitalisation
                                    end
                                end
                            end
                        end
                    end
                    
                    for l=0:stageTF.getLength()-1
                        if string(stageTF.item(l).getNodeName()) == "PolesZeros"
                            polesZeros = stageTF.item(l);
                            cntrZpk = cntrZpk + 1;
                            zpk(cntrZpk).poles = [];
                            zpk(cntrZpk).zeros = [];
                            for ll = 0:polesZeros.getLength()-1
                                if string(polesZeros.item(ll).getNodeName()) == "NormalizationFactor"
                                    zpk(cntrZpk).gain = str2double(polesZeros.item(ll).getTextContent());
                                end
                            end
                            transferFunction = true;
                        end
                    end
                    
                    % find stage gain
                    if transferFunction
                        for l=0:stageTF.getLength()-1
                            if string(stageTF.item(l).getNodeName()) == "StageGain"
                                stageGain = stageTF.item(l).getChildNodes();
                                for ll = 0:stageGain.getLength()-1
                                    if string(stageGain.item(ll).getNodeName()) == "Value"
                                        zpk(cntrZpk).gain = zpk(cntrZpk).gain*str2double(stageGain.item(ll).getTextContent());
                                    end
                                end
                            end
                        end    
                    end
                    
                    % find poles and zeros and update total gain
                    if transferFunction
                        for l = 0:polesZeros.getLength()-1
        %                     disp(polesZeros.item(l).getNodeName())
                            if polesZeros.item(l).getLength() > 1 
                                if string(polesZeros.item(l).getNodeName) == "Pole"
                                    pole = 0;
                                    for pz = 0:polesZeros.item(l).getLength()-1
                                        if string(polesZeros.item(l).item(pz).getNodeName()) == "Imaginary"
                                            pole = pole + str2num(polesZeros.item(l).item(pz).getTextContent())*1i;
                                        else
                                            pole = pole + str2num(polesZeros.item(l).item(pz).getTextContent());
                                        end                                    
        %                                 disp(polesZeros.item(l).item(pz).getTextContent())
                                    end
                                    zpk(cntrZpk).poles = [zpk(cntrZpk).poles; pole];
                                elseif string(polesZeros.item(l).getNodeName) == "Zero"
                                    zero = 0;
                                    for pz = 0:polesZeros.item(l).getLength()-1
                                        if string(polesZeros.item(l).item(pz).getNodeName()) == "Imaginary"
                                            zero = zero + str2num(polesZeros.item(l).item(pz).getTextContent())*1i;
                                        else
                                            zero = zero + str2num(polesZeros.item(l).item(pz).getTextContent());
                                        end                                    
        %                                 disp(polesZeros.item(l).item(pz).getTextContent())
                                    end
                                    zpk(cntrZpk).zeros = [zpk(cntrZpk).zeros; zero];                            
                                end
                            else
        %                         disp(polesZeros.item(l).getTextContent())
                            end
        %                     disp("-----")
                        end
                    end
                        
                    % calculate stage 2 gain (i.e. volts to counts)
                    if transferFunction
                        for l=0:stageDigit.getLength()-1
                            if string(stageDigit.item(l).getNodeName()) == "StageGain"
                                stageGain = stageDigit.item(l).getChildNodes();
                                for ll = 0:stageGain.getLength()-1
                                    if string(stageGain.item(ll).getNodeName()) == "Value"
                                        zpk(cntrZpk).volts2counts= str2double(stageGain.item(ll).getTextContent());
                                    end
                                end
                            end
                        end    
                    end
                end
            end
        end
    end
    
 
end
