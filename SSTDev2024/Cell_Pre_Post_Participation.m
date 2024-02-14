%% Calculate the fraction of cells recruited after WS (45 frames after WS (Line 63))
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
    load('Fallall.mat', 'spikenums','spike_rates');
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
    
    pre=[];
    post=[];
    for j=2:length(stims)-1
        try
            ab=spikenums(:,stims(j)-45:stims(j)+45);% change for plotting
            a_1=sum(ab(:,1:45),2);
            a_2=sum(ab(:,47:end),2);
            a_1(a_1>0)=1;
            a_2(a_2>0)=1;
            as1=sum(a_1)/size(ab,1)*100;
            as2=sum(a_2)/size(ab,1)*100;
            pre(end+1)=as1;
            post(end+1)=as2;
        catch
        end
    end
     pre=mean(pre);
     post=mean(post);
end

C{kk}=pre;
C2{kk}=post;
end

cd(BasePath);

CC=cell2mat(C); %Variable giving the % of cells that are active before WS (Average of all Stims in a movie)
CC=CC'; 

CC2=cell2mat(C2);%Variable giving the % of cells that are active after WS (Average of all Stims in a movie)
CC2=CC2';


