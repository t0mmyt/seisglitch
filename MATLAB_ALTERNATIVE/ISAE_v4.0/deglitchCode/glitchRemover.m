function [deglitchedData, dataGlitches, glInfo] = glitchRemover(data, glLoc, refGlLoc, lGl, glDiff)

% data = processedData{ch};
% glLoc = glLoc{ch};
% refGlLoc = refGlLoc{ch};
% lGl= lGl{ch};

lenGl = length(lGl(:,1));
nShapes = length(lGl(1,:));
glNumber = length(glLoc);

% inverse problem
A = cell(nShapes,1);
AA = cell(nShapes,1);
for j=1:nShapes
    % A{j} = (glitch, 0 to 1, 1 to 0)
    A{j} = [lGl(:,j) ([0:lenGl-1]/(lenGl-1))' (1-[0:lenGl-1]/(lenGl-1))'];
    AA{j} = pinv(A{j});
end

deglitchedData = data;
dataGlitches = zeros(length(data), 1);
glAmp = zeros(glNumber,1);
y0 = zeros(glNumber,1);
slope = zeros(glNumber,1);
model = cell(glNumber,nShapes);

for nGl=1:glNumber % for each selected glitch
    dataToInvert = data(glLoc(nGl):glLoc(nGl)+lenGl-1);
    j=refGlLoc(nGl);
    model{nGl,j} = (AA{j}*dataToInvert)'; % A CHANGER : MIS UN POINT .
    glAmp(nGl) = model{nGl,j}(1);
    slope(nGl) = (model{nGl,j}(2)-model{nGl,j}(3))/lenGl;
    y0(nGl) = model{nGl,j}(3) + slope(nGl)*glDiff;
    deglitchedData(glLoc(nGl):glLoc(nGl)+lenGl-1) = ...
        data(glLoc(nGl):glLoc(nGl)+lenGl-1) ...
        - model{nGl,j}(1).*lGl(:,j);
    dataGlitches(glLoc(nGl):glLoc(nGl)+lenGl-1)= ...
        dataGlitches(glLoc(nGl):glLoc(nGl)+lenGl-1)+ ...
        (model{nGl,j}(1).*lGl(:,j)+ ...
        (model{nGl,j}(2)-model{nGl,j}(3)).*([0:lenGl-1]/(lenGl-1))' + ...
        model{nGl,j}(3).*ones(lenGl,1));
end

glInfo = struct('glAmp', glAmp, 'y0', y0, 'slope', slope);

end














