clear
%load Fallall for FOV1
[FileName,PathName,FilterIndex]=uigetfile('*.mat');
addpath(PathName);
load(FileName,'spikenums');
load(FileName,'Coor');
load(FileName,'spike_rates');

Coor1=Coor;
spikenums1=spikenums;
spike_rates1=spike_rates;

clear spikenums Coor spike_rates

%load Fallall for FOV2
[FileName,PathName,FilterIndex]=uigetfile('*.mat');
addpath(PathName);
load(FileName,'spikenums');
load(FileName,'Coor');
load(FileName,'spike_rates');

Coor2=Coor;
spikenums2=spikenums;
spike_rates2=spike_rates;

Coor_C=horzcat(Coor1,Coor2);
spikenums_C=vertcat(spikenums1,spikenums2);
spike_rates_C=vertcat(spike_rates1,spike_rates2);

%Get location of cells
xcen=[];
ycen=[];
zcen=[];
for i=1:length(Coor_C)
    try
    cnttemp=Coor_C{i};
    xcen(i)=mean(cnttemp(2,:));
    ycen(i)=mean(cnttemp(1,:));
    catch
    xcen(i)=250;
    ycen(i)=250;
    end
    
    if i<=length(Coor1)
       zcen(i)=1; 
    else
        zcen(i)=250;
    end
    
end

%pixel size
xcen=xcen*1.2695;
ycen=ycen*1.2695;

Coordinates=vertcat(xcen,ycen,zcen)';
Distance_pairs= squareform(pdist(Coordinates)); %matrix with the distance between cell pairs

Correlation_pairs=zeros(length(xcen),length(xcen));

for i=1:length(xcen)
    for j=i:length(xcen)
    X = spike_rates_C(i,:);
    Y = spike_rates_C(j,:);
    d = nancorr(X',Y');
    Correlation_pairs(i,j)=d;
    end
end

Correlation_pairs = Correlation_pairs - diag(diag(Correlation_pairs));
Distance_pairs = Distance_pairs - diag(diag(Distance_pairs));

mask = triu(true(size(Distance_pairs)),1);
Correlation_vector = Correlation_pairs(mask);
Distance_vector=Distance_pairs(mask);

figure;scatter(Distance_vector,Correlation_vector)

SS = FileName(1:end-4);
SSS='Scatter_correlation_distance_General.fig';
St=strcat(SS,SSS);
savefig(St)

Superficial_correlation=Correlation_pairs(1:length(Coor1),1:length(Coor1));
Superficial_distance=Distance_pairs(1:length(Coor1),1:length(Coor1));

mask = triu(true(size(Superficial_distance)),1);
Correlation_S = Superficial_correlation(mask);
Distance_S=Superficial_distance(mask);

figure;scatter(Distance_S,Correlation_S)

SS = FileName(1:end-4);
SSS='Scatter_correlation_distance_Superficial.fig';
St=strcat(SS,SSS);
savefig(St)



Deep_correlation=Correlation_pairs(length(Coor1)+1:end,length(Coor1)+1:end);
Deep_distance=Distance_pairs(length(Coor1)+1:end,length(Coor1)+1:end);

mask = triu(true(size(Deep_distance)),1);
Correlation_D = Deep_correlation(mask);
Distance_D=Deep_distance(mask);

figure;scatter(Distance_D,Correlation_D)

SS = FileName(1:end-4);
SSS='Scatter_correlation_distance_Deep.fig';
St=strcat(SS,SSS);

savefig(St)


Between_correlation=Correlation_pairs(1:length(Coor1),length(Coor1)+1:end);
Between_distance=Distance_pairs(1:length(Coor1),length(Coor1)+1:end);


%show spatial correlation and cell clusters

B=Correlation_pairs;
B=triu(B)+triu(B,1)';


Race = B;
[NCell,NRace] = size(Race);
[IDX2,sCl,M,S] = kmeansopt(Race,100,'var');
NCl = max(IDX2);
[~,x2] = sort(IDX2);
MSort = M(x2,x2);

figure;
imagesc(MSort)

SS = FileName(1:end-4);
SSS='Cell_Clusters.fig';
St=strcat(SS,SSS);

savefig(St)

figure;
scatter3(xcen,ycen,zcen,[],IDX2,'filled')
xlim([0 650])
ylim([0 650])
zlim([0 250])


SS = FileName(1:end-4);
SSS='Spatial_clusters.fig';
St=strcat(SS,SSS);

savefig(St)

%% Spatial organization of Large and Small events

% make surrogate distribution
spikenums22=spikenums_C(:,1:end);
for j=1:500
thresh=zeros(length(spike_rates_C(:,1)),length(spikenums22(1,:)));
parfor i=1:length(spike_rates_C(:,1))
thresh(i,:)=circshift(spikenums22(i,:), randi([1, length(spikenums22(1,:))]));
end
summ=sum(thresh,1);
summ2{j}=summ;
end
tress = cell2mat(summ2);
thresh = prctile(tress,99);

for i=1:length(spikenums_C)
    y(i)=sum(spikenums_C(:,i));
end


y = medfilt1(y,10);
[pks,locs] = findpeaks(y,'MinPeakDistance',9,'MinPeakHeight',thresh);

figure;
imagesc(spikenums_C)
hold on;
plot(y,'color','r')

findpeaks(y,'MinPeakDistance',9,'MinPeakHeight',thresh);
SS = FileName(1:end-4);
SSS='Events_Planes_combined.fig';
St=strcat(SS,SSS);
savefig(St)

for i=1:length(spikenums_C)
    y(i)=sum(spikenums_C(:,i));
end

y = medfilt1(y,10);


up=find(y>thresh);

A=up;
A(end+1)=2;   % adds new endpoint to very end of A so code picks up end of last group of consecutive values
I_1=find(diff(A)~=1);  % finds where sequences of consecutive numbers end
[m,n]=size(I_1);   % finds dimensions of I_1 i.e. how many sequences of consecutive numbers you have
startpoint=1;    % sets start index at the first value in your array
seq=cell(1,n);  % had to preallocate because without, it only saved last iteration of the for loop below
               % used n because this array is a row vector
for i=1:n
    End_Idx=I_1(i);   %set end index
    seq{i}=A(startpoint:End_Idx);  %finds sequences of consecutive numbers and assigns to cell array
    startpoint=End_Idx+1;   %update start index for the next consecutive sequence
end

for i=1:length(seq)
    len_seq(i)=length(seq{i});
end


participation=zeros(length(seq),length(spikenums_C(:,1)));
for i=1:length(seq)
    aaa=[];
    bbb=[];
    aaa=seq{i};
    bbb=spikenums_C(:,aaa);
    if length(aaa)>1
    participation(i,:)=sum(bbb')/length(aaa);
    else
    participation(i,:)=bbb;
    end
end

participation=participation';


maxpeakcells=zeros(length(seq),1);
for i=1:length(seq)
    aaa=[];
    bbb=[];
    ccc=[];
    aaa=seq{i};
    bbb=spikenums_C(:,aaa);
    ccc=sum(bbb');
    ccc=ccc';
    ccc=find(ccc>0);
    ccc=length(ccc)/length(spikenums_C(:,1));
    maxpeakcells(i)=ccc;
end

%find Large and Small events
large3=0.8;
long3=find(maxpeakcells>=large3);
short3=find(maxpeakcells<large3);

%Binarize matrix with spikes
participation2=participation;
participation2(participation2>0)=1;

L_events=participation2(:,long3);
S_events=participation2(:,short3);

%% Cluster different event types (Large and Small)
%Large events
%start cluster algo, run time should probably be improved
Race = S_events;
[NCell,NRace] = size(Race);
[IDX2,sCl,M,S] = kmeansopt(Race,100,'var');
% M = CovarM(Race);
% IDX2 = kmedoids(M,NCl);
NCl = max(IDX2);

[~,x2] = sort(IDX2);
MSort = M(x2,x2);

%Race clusters
R = cell(0);
CellScore = zeros(NCell,NCl);
CellScoreN = zeros(NCell,NCl);
for i = 1:NCl
    R{i} = find(IDX2==i);
    CellScore(:,i) = sum(Race(:,R{i}),2);
    CellScoreN(:,i) = CellScore(:,i)/length(R{i});
end
%Assign cells to cluster with which it most likely spikes
[~,CellCl] = max(CellScoreN,[],2);
%Remove cells with less than 2 spikes in a given cluster
CellCl(max(CellScore,[],2)<2) = 0;
[X1,x1] = sort(CellCl);

figure
subplot(1,2,1)
imagesc(MSort)
colormap jet
axis image

ylabel('RACE #')

subplot(1,2,2)
imagesc(Race(x1,x2),[-1 1.2])
axis image
xlabel('RACE #')
ylabel('Cell #')
SS = FileName(1:end-4);
SSS='Clusters_S_events.fig';
St=strcat(SS,SSS);

savefig(St)


%% Save Clusters
SS = FileName(1:end-4);
SSS='Clusters_S_events.mat';
St=strcat(SS,SSS);
save(St,'IDX2')

%% Remove cluster non-statistically significant

sClrnd = zeros(1,1000);
for i = 1:1000
    sClrnd(i) = kmeansoptrnd(Race,10,NCl);
end
% TODO: replace max by percentile (done)
NClOK = sum(sCl>prctile(sClrnd,95));
sClOK = sCl(1:NClOK)';

%save('NClustersOK.mat','NClOK')

RaceOK = Race(:,IDX2<=NClOK);
NRaceOK = size(RaceOK,2);

%% Statistical definition of cell assemblies

NShuf = 5000;
%Count number of participation to each cluster
CellP = zeros(NCell,NCl); CellR = zeros(NCell,NCl);
for i = 1:NCl
    CellP(:,i) = sum(Race(:,IDX2 == i),2);
    CellR(:,i) = CellP(:,i)/sum(IDX2 == i);
end

%Test for statistical significance
CellCl = zeros(NCl,NCell); %Binary matrix of cell associated to clusters
for j = 1:NCell
    %Random distribution among Clusters
    RClr = zeros(NCl,NShuf);
    Nrnd = sum(Race(j,:));
    parfor l = 1:NShuf
        Random = randperm(NRace);
        Random = Random(1:Nrnd);
        Racer = zeros(1,NRace);
        Racer(Random) = 1;
        for i = 1:NCl
            RClr(i,l) = sum(Racer(:,IDX2 == i),2);
        end
    end
    RClr = sort(RClr,2);
    %         ThMin = mean(Random) - 2*std(Random);
    %Proba above 95th percentile
    ThMax = RClr(:,round(NShuf*(1-0.05/NCl))); 
    for i = 1:NCl
        CellCl(i,j) = double(CellP(j,i)>ThMax(i));% - double(RCl(:,j)<ThMin);
    end
end

A0 = find(sum(CellCl) == 0); %Cells not in any cluster
A1 = find(sum(CellCl) == 1); %Cells in one cluster
A2 = find(sum(CellCl) >= 2); %Cells in several clusters

%Keep cluster where they participate the most
for i = A2
    [~,idx] = max(CellR(i,:));
    CellCl(:,i) = 0;
    CellCl(idx,i) = 1;
end
C0 = cell(0);
k = 0;
for i = 1:NCl
    if length(find(CellCl(i,:)))>5
        k = k+1;
        C0{k} = find(CellCl(i,:));
    end
end

SS = FileName(1:end-4);
SSS='CellList_S_events.mat';
St=strcat(SS,SSS);
save(St,'C0');

%Participation rate to its own cluster
SS = FileName(1:end-4);
SSS='ClusPart_S_events.mat';
St=strcat(SS,SSS);
CellRCl = max(CellR([A1 A2],:),[],2);
save(St,'CellRCl')
%% Assign RACE to groups of cells

NCl = length(C0);
[NCell,NRace] = size(Race);
%Cell count in each cluster
RCl = zeros(NCl,NRace);
PCl = zeros(NCl,NRace);
for i = 1:NCl
    RCl(i,:) = sum(Race(C0{i},:));
end

RCln = zeros(NCl,NRace);
for j = 1:NRace
    %Random distribution among Clusters
    RClr = zeros(NCl,NShuf);
    Nrnd = sum(Race(:,j));
    parfor l = 1:NShuf
        Random = randperm(NCell);
        Random = Random(1:Nrnd);
        Racer = zeros(NCell,1);
        Racer(Random) = 1;
        for i = 1:NCl
            RClr(i,l) = sum(Racer(C0{i}));
        end
    end
    %         ThMin = mean(Random) - 2*std(Random);
    RClr = sort(RClr,2);
    %         ThMin = mean(Random) - 2*std(Random);
    %Proba above 95th percentile
    ThMax = RClr(:,round(NShuf*(1-0.05/NCl)));
    for i = 1:NCl
        PCl(i,j) = double(RCl(i,j)>ThMax(i));% - double(RCl(:,j)<ThMin);
    end
    %Normalize (probability)
    RCln(:,j) = RCl(:,j)/sum(Race(:,j));
end

SS = FileName(1:end-4);
SSS='CellCluster_S_events.mat';
St=strcat(SS,SSS);

save(St,'PCl');

%% Show Sorted Rasterplot

%load('RaceCl2')
%load('CellClList2.mat')
%load('Race.mat')
[NCell,NRace] = size(Race);
NCl = length(C0);

%Recreate CellCl (equivalent of RCl for the cells)
CellCl = zeros(NCl,NCell);
for i = 1:NCl
    CellCl(i,C0{i}) = 1;
end

NCl = length(C0);

Cl0 = find(sum(PCl,1) == 0);
Cl1 = find(sum(PCl,1) == 1);
Cl2 = find(sum(PCl,1) == 2);
Cl3 = find(sum(PCl,1) == 3);
Cl4 = find(sum(PCl,1) == 4);

Bin = 2.^(0:NCl-1);

%Sort Cl1
[~,x01] = sort(Bin*PCl(:,Cl1));
Cl1 = Cl1(x01);

%Sort Cl2
[~,x02] = sort(Bin*PCl(:,Cl2));
Cl2 = Cl2(x02);

%Sort Cl3
[~,x03] = sort(Bin*PCl(:,Cl3));
Cl3 = Cl3(x03);

RList = [Cl0 Cl1 Cl2 Cl3 Cl4];
%x1 from DetectRace

[X1,x1] = sort(Bin*CellCl);

figure
imagesc(Race(x1,RList))
colormap hot
axis image

SS = FileName(1:end-4);
SSS='CellCluster_S_events.fig';
St=strcat(SS,SSS);
savefig(St)


%% Show 3D correlation
Clus=zeros(1,length(Coor_C));
Clus(CellCl(1,:)==1)=1;
Clus(CellCl(2,:)==1)=2;

a=find(Clus>0);

figure;
scatter3(xcen(a),ycen(a),zcen(a),[],Clus(a),'filled')

SS = FileName(1:end-4);
SSS='Spatial_clusters_S_events.fig';
St=strcat(SS,SSS);

savefig(St)

SS = FileName(1:end-4);
SSS='volumetric_all.mat';
St=strcat(SS,SSS);

save(St,'-v7.3');

