Ag=7;  %Age of the Animal
Dep='Sup';  %either 'Sup' or 'Deep'

BasePath='F:\Calcium Imaging\Modulation_of_SST_int\Analysis\Evoked\p7';

a=Age==Ag;
tf1 = strcmp('Sup',Dep);
if tf1==1
b=Depth<=300;
else
b=Depth>=200;
end

cases=a&b;
cases=find(cases==1);
C={};
C2={};


for kk=1:length(cases)
    a=fileList(cases(kk)).folder;
    cd(a)
    %%kk check running cases
    load('Fallall.mat', 'spike_rates');
    spike_rates(isnan(spike_rates))=0;
    
   % for i=1:size(spike_rates,1)
   % spike_rates(i,:)=normalize(spike_rates(i,:),'center','mean');
   % end
    
    fstruct = dir('*allClusters*.mat');
    load(fstruct.name,'idd2');
    idd=zeros(1,size(spike_rates,1));
    idd(idd2)=1;

    load('Fall.mat', 'Fneu', 'F','iscell');
    isc=find(iscell(:,1)==1);
    Fneu=Fneu(isc,:);
    F=F(isc,:);
    
    Fnew=zeros(size(F));
    for i=1:length(Fneu(:,1))
        Fnew(i,:)=normalize((F(i,:)-0.7*Fneu(i,:))/median(Fneu(i,:)));
    end

    stims=[];
    fstruct2 = dir('*DataFile_*.mat');
    for jj=1:length(fstruct2)
       load(fstruct2(jj).name,'stimTimes');
       stimTimes2=stimTimes+(jj-1)*60;
       stims=[stims stimTimes2];
    end
    
stims=stims*30;

if size(F,2)==9000|size(F,2)==7200|size(F,2)==5400|size(F,2)==9400|size(F,2)==9428|size(F,2)==9576|size(F,2)==9708|size(F,2)==7640
    
    
    
    B1=[];
    B2=[];
    for jj=1:size(Fneu,1)
        A1=[];
        A2=[];
        if idd(jj)==0
            for j=2:length(stims)-1
                try
                    ab=spike_rates(jj,stims(j)-15:stims(j)+15);% change for plotting
                    %ab=find(ab==max(ab));
                    A1(end+1,:)=ab;
                catch
                end
            end
            
            tmp=nanmean(A1);
            tmpMean=mean(tmp(1:15));
            tmp2=tmp(16:end);
            maxRespWS=max(tmp2);
            
            if maxRespWS>tmpMean+(1.984*std(tmp(1:15))) %One-tail significance 0.05-0.02 df=1000
               B1(end+1,:)=find(tmp2==max(tmp2));
            end
        else
            for j=2:length(stims)-1
                try
                    ab=spike_rates(jj,stims(j)-15:stims(j)+15);% change for plotting
                    %ab=find(ab==max(ab));
                    A2(end+1,:)=ab;
                catch
                end
            end
            tmp=nanmean(A2);
            tmpMean=mean(tmp(1:15));
            tmp2=tmp(16:end);
            maxRespWS=max(tmp2);
             if maxRespWS>tmpMean+(1.984*std(tmp(1:15))); %One-tail significance 0.05-0.02 df=1000
               B2(end+1,:)=find(tmp2==max(tmp2));
             end
        end
    end

C{kk}=B1';
C2{kk}=B2';
end
end

cd(BasePath);

CC=cell2mat(C); %Variable givin the latency to max cell response after ws (all cells)
CC=CC'; 

CC2=cell2mat(C2);%Variable givin the latency to max cell response after ws (Red cells)
CC2=CC2';



