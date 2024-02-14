clear
%load volumetric_Fallall
[FileName3,PathName3,FilterIndex3]=uigetfile('*.mat');
load(FileName3,'spikenums1','spikenums2');

spikenums_C=vertcat(spikenums1,spikenums2);

%% make surrogate distribution
for j=1:500
thresh=zeros(length(spikenums_C(:,1)),length(spikenums_C(1,:)));
parfor i=1:length(spikenums_C(:,1))
thresh(i,:)=circshift(spikenums_C(i,:), randi([1, length(spikenums_C(1,:))]));
end
summ=sum(thresh,1);
summ2{j}=summ;
end
tress = cell2mat(summ2);
thresh = prctile(tress,99);

y=sum(spikenums_C,1);

y = medfilt1(y,10);
[pks,locs] = findpeaks(y,'MinPeakDistance',10,'MinPeakHeight',thresh);

figure;
imagesc(spikenums_C)
hold on;
plot(y,'color','r')

findpeaks(y,'MinPeakDistance',10,'MinPeakHeight',thresh);
SS = FileName3(1:end-4);
SSS='Events_Combined.fig';
St=strcat(SS,SSS);
savefig(St)

%% postprocessing
lb=10; %frames before
hb=10; %frames after vector of concatenated vector should be center at 0
if (locs(end)+hb)>length(spikenums_C(1,:))
    locs=locs(1:end-1);
    pks=pks(1:end-1);
end

if (locs(1)-lb)<1
    locs=locs(2:end);
    pks=pks(2:end);
end

att=size(spikenums_C,1);

t1=NaN(length(locs),att);
t2=zeros(length(locs),att);


for i=1:length(pks)
    for j=1:att
        axxx=spikenums_C(j,locs(i)-lb:locs(i)+hb);
        ax=min(find(axxx==1));
        if ax>=1
            t1(i,j)=ax;
        end
    end
end

tttt=round(length(-lb:hb)/2);
t11=t1-tttt;
t12=t11;
t12(~isnan(t11))=1;
t12(isnan(t12))=0;

%% ab and ab2 variables
t1=NaN(length(locs),length(spikenums_C(:,1)));
t2=zeros(length(locs),length(spikenums_C(:,1)));


for i=1:length(pks)
    for j=1:length(spikenums_C(:,1))
        axxx=spikenums_C(j,locs(i)-10:locs(i)+10);
        ax=min(find(axxx==1));
        if ax>=1
            t1(i,j)=ax;
        end
    end
end

t11=t1-10;


for i=1:length(t11(1,:))
    ab(i)=((length(t11(:,i))-sum(isnan(t11(:,i))))/length(locs))*100;
    ab2(i)=nanmean(t11(:,i));
end

clearvars t1 t2 tttt axxx ax lb hb summ summ2 tress

%%
mn1=[];
mn2=[];
mmn1=[];
mmn2=[];
figure;
hold on;
for i = 1:att
    if i<= size(spikenums1,1)
        scatter(ab2(i),ab(i),'r') %%this are deep planes
        mn1(end+1)=ab2(i);
        mmn1(end+1)=ab(i);
    end
    if i> size(spikenums1,1)
        scatter(ab2(i),ab(i),'b')
        mn2(end+1)=ab2(i);
        mmn2(end+1)=ab(i);
    end
end

pd = fitdist(mn2','Normal');
x_pdf = [-20:0.5:20];
y = pdf(pd,x_pdf);

pd2 = fitdist(mn1','Normal');
x_pdf2 = [-20:0.5:20];
y2 = pdf(pd2,x_pdf2);


hold on
scale = max(mmn2);
scale2=max(mmn1);
h=area((x_pdf),(y.*scale))
h(1).FaceColor = [0 0 1];
alpha(h(1),0.2)

h=area((x_pdf2),(y2.*scale2));
h(1).FaceColor = [1 0 0];
alpha(h(1),0.2)


SS = 'red';
SSS='OverlapScatter_superficial_deep.fig';
St=strcat(SS,SSS);

savefig(St);


mn1=mn1(~isnan(mn1));
mn2=mn2(~isnan(mn2));
[h,p,ci,stats]=ttest2(mn1,mn2)
ordercells=mean(mn2)-mean(mn1)

SS = FileName3(1:end-4);
SSS='stats_C.mat';
St=strcat(SS,SSS);
save(St,'h','p','ci','stats','ordercells');


