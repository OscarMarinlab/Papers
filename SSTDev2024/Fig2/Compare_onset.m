Ag=7;  %Age of the Animal
Dep='Dep';  %either 'Sup' or 'Deep'

BasePath='E:\Calcium Imaging\Modulation_of_SST_int\Analysis\Evoked';

a=Age==Ag;
tf1 = strcmp('Sup',Dep);
if tf1==1
b=Depth<=200;
else
b=Depth>=200;
end

cases=a&b;
cases=find(cases==1);

for kk=1:length(cases)
    
  a=fileList(cases(kk)).folder;
    cd(a)
try        
load('Fallall.mat', 'ab', 'ab2');
fstruct = dir('*allCLusters*.mat');


load(fstruct.name,'idd2');
idd=zeros(1,length(ab));
idd(idd2)=1;

mn1=[];
mn2=[];
mmn1=[];
mmn2=[];
for i = 1:length(ab2)
    if idd(i)==1
        mn1(end+1)=ab2(i);
        mmn1(end+1)=ab(i);
    end
    if idd(i)==0
        mn2(end+1)=ab2(i);
        mmn2(end+1)=ab(i);
    end
end

mn11{kk}=mn1;
mn22{kk}=mn2;
mmn11{kk}=mmn1;
mmn22{kk}=mmn2;
catch
end
end

cd(BasePath);


mn1=cell2mat(mn11);
mn2=cell2mat(mn22);
mmn1=cell2mat(mmn11);
mmn2=cell2mat(mmn22);

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
h(1).FaceColor = [0 0 1];
alpha(h(1),0.2)

h=area((x_pdf2),(y2.*scale2));
h(1).FaceColor = [1 0 0];
alpha(h(1),0.2)

St=strcat('P',num2str(Ag),'_Depth_',Dep);
savefig(St);


mn1=mn1(~isnan(mn1));
mn2=mn2(~isnan(mn2));
[h,p,ci,stats]=ttest2(mn1,mn2)
ordercells=mean(mn2)-mean(mn1)

St2=strcat(St,'_stats.mat');
save(St2,'h','p','ci','stats','ordercells');