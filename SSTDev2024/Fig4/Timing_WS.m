%Open spikenums in Fallall.mat 
stims=[];
    fstruct2 = dir('*DataFile_*.mat');
    for jj=1:length(fstruct2)
       load(fstruct2(jj).name,'stimTimes');
       stimTimes2=stimTimes+(jj-1)*60;
       stims=[stims stimTimes2];
    end
    
stims=stims*30;

imagesc(spikenums)
hold on;

for i=1:length(stims)
    xline(stims(i),'red')
end