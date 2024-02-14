%% Calculate the distribution of cell reponse in LSE (or High synchronous events) 
%Compares distributions between RED cells and other cells at each age

Ag=8;  %Age of the Animal
Dep='Sup';  %either 'Sup' or 'Deep'

BasePath='F:\Calcium Imaging\Modulation_of_SST_int\Analysis\Spont_Mix_ages_Distribution_analysis_p_eight';

a=Age==Ag;
tf1 = strcmp('Sup',Dep);
if tf1==1
b=Depth<=500;
else
b=Depth>=250;
end

cases=a&b;
cases=find(cases==1);

red_l={};
nonred_l={};
partred_l={};
partnonred_l={};

for kk=1:length(cases)
    
  a=fileList(cases(kk)).folder;
    cd(a)
try        
load('Fallall.mat', 'ab', 'ab2','spikenums','thresh');
fstruct = dir('*allCLusters*.mat');


load(fstruct.name,'idd2');
idd=zeros(1,length(ab));
idd(idd2)=1;

for i=1:length(spikenums)
    y(i)=sum(spikenums(:,i));
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

nn=1:length(spikenums(:,1));
novl=setdiff(nn,idd2);

maxpeakcells=zeros(length(seq),1);
for i=1:length(seq)
    aaa=[];
    bbb=[];
    ccc=[];
    aaa=seq{i};
    bbb=spikenums(:,aaa);
    ccc=sum(bbb');
    ccc=ccc';
    ccc=find(ccc>0);
    ccc=length(ccc)/length(spikenums(:,1));
    maxpeakcells(i)=ccc;
end

participation=zeros(length(ab2),length(seq));
for i=1:length(seq)
    aaa=[];
    bbb=[];
    ccc=[];
    aaa=seq{i};
    bbb=spikenums(:,aaa);
    ax=[];
    for j=1:length(ab2)
        axx=min(find(bbb(j,:)==1));
        if axx>=1
            ax(j)=axx;
        else
            ax(j)=NaN;
        end
    end
    ax=ax';
    ax=ax-(len_seq(i)/2);
    ax=ax.*(1/30);
    participation(:,i)=ax;
end



large2=0.8;
long2=find(maxpeakcells>=large2);

firing_red=participation(idd2,long2);
firing_nonred=participation(novl,long2);


red_l{kk}=reshape(firing_red,[],1);
nonred_l{kk}=reshape(firing_nonred,[],1);


participation=zeros(length(seq),length(spikenums(:,1)));
for i=1:length(seq)
    aaa=[];
    bbb=[];
    aaa=seq{i};
    bbb=spikenums(:,aaa);
    if length(aaa)>1
    participation(i,:)=sum(bbb')/length(aaa);
    else
    participation(i,:)=bbb;
    end
end

participation=participation';

firing_red=participation(idd2,long2);
firing_nonred=participation(novl,long2);


partred_l{kk}=reshape(firing_red,[],1);
partnonred_l{kk}=reshape(firing_nonred,[],1);


catch
end
end

cd(BasePath);


mn1=cell2mat(red_l');
mn1 = mn1(~isnan(mn1))';
mn2=cell2mat(nonred_l');
mn2 = mn2(~isnan(mn2))';
mmn1=cell2mat(partred_l');
mmn1 = mmn1(~isnan(mmn1))';
mmn2=cell2mat(partnonred_l');
mmn2 = mmn2(~isnan(mmn2))';

pd = fitdist(mn2','Normal');
x_pdf = [-30:0.5:10];
y = pdf(pd,x_pdf);

pd2 = fitdist(mn1','Normal');
x_pdf2 = [-30:0.5:10];
y2 = pdf(pd2,x_pdf2);


hold on
scale = max(mmn2);
scale2=max(mmn1);
h=area((x_pdf),(y.*scale))
h(1).FaceColor = [1 0 0];
alpha(h(1),0.2)

h=area((x_pdf2),(y2.*scale2));
h(1).FaceColor = [0 0 1];
alpha(h(1),0.2)

St=strcat('P',num2str(Ag),'_Depth_',Dep,'High');
savefig(St);


mn1=mn1(~isnan(mn1));
mn2=mn2(~isnan(mn2));
[h,p,ci,stats]=ttest2(mn1,mn2)
ordercells=mean(mn2)-mean(mn1)

St2=strcat(St,'_stats.mat');
save(St2,'h','p','ci','stats','ordercells');
