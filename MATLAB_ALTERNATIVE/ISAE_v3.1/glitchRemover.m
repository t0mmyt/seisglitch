function [deglitchedData, dataGlitches, glAmp] = glitchRemover(data, glLoc, refGlLoc, lGl)

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
    A{j} = [lGl(:,j) ([0:lenGl-1]/(lenGl-1))' (1-[0:lenGl-1]/(lenGl-1))'];
    AA{j} = pinv(A{j});
end

deglitchedData = data;
dataGlitches = zeros(length(data), 1);
glAmp = zeros(glNumber,1);
model = cell(glNumber,nShapes);

for nGl=1:glNumber % for each selected glitch
    dataToInvert = data(glLoc(nGl):glLoc(nGl)+lenGl-1);
    j=refGlLoc(nGl);
    model{nGl,j} = (AA{j}*dataToInvert)'; % A CHANGER : MIS UN POINT .
    glAmp(nGl) = model{nGl,j}(1);
    deglitchedData(glLoc(nGl):glLoc(nGl)+lenGl-1) = ...
        data(glLoc(nGl):glLoc(nGl)+lenGl-1) ...
        - model{nGl,j}(1).*lGl(:,j);
    dataGlitches(glLoc(nGl):glLoc(nGl)+lenGl-1)= ...
        dataGlitches(glLoc(nGl):glLoc(nGl)+lenGl-1)+ ...
        (model{nGl,j}(1).*lGl(:,j)+ ...
        (model{nGl,j}(2)-model{nGl,j}(3)).*([0:lenGl-1]/(lenGl-1))' + ...
        model{nGl,j}(3).*ones(lenGl,1));
end

end














