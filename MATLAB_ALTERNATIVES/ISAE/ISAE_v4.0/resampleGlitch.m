function newGlitch = resampleGlitch(glitch,freqDiv,PFO,gain)

lenGl=length(glitch);
lenPFO=length(PFO);

newSample=1;
for sample=1:freqDiv:lenGl
    if sample < lenPFO
        newGlitch(newSample)=0;
    else
        newGlitch(newSample)=0;
        for i=1:lenPFO
            newGlitch(newSample)=newGlitch(newSample)+PFO(i)*glitch(sample-lenPFO+i);
        end
    end
    newGlitch(newSample)=newGlitch(newSample)/gain;
    newSample=newSample+1;
end

end
