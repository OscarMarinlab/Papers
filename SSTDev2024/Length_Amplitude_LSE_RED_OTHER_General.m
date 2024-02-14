Ag=10;  %Age of the Animal
Dep='Dep';  %either 'Sup' or 'Deep'

BasePath='F:\Calcium Imaging\Modulation_of_SST_int\Analysis\Spontaneous';
a=Age==Ag;
tf1 = strcmp('Sup',Dep);
if tf1==1
    b=Depth<=200;
else
    b=Depth>200;
end

cases=a&b;
cases=find(cases==1);


for kk=1:length(cases)
    a=fileList(cases(kk)).folder;
    cd(a)
    
    load('Fallall.mat', 'F','Fneu','iscell');
    Fnew=zeros(size(F));
    for i=1:length(iscell(:,1))
        Fnew(i,:)=(F(i,:)-0.7*Fneu(i,:))/median(Fneu(i,:));
    end
    
    a1=iscell(:,1)==1;
    a1=find(iscell(:,1)==1);
    Fnew=Fnew(a1,:);
    
    try
        fstruct = dir('*Fallallparticipation*.mat');
        load(fstruct.name);
        
        length_Hevents=len_seq(long2);
        length_Levents=len_seq(short2);
        
        participation2=participation;
        participation2(participation>0)=1;
        participation2=sum(participation2);
        
        participation_Hevents=participation2(long2);
        participation_Levents=participation2(short2);
        
        Ampl_Hevents=[];
        Ampl_Levents=[];
        seq2=seq(long2);
        for i=1:size(Fnew,1)
            for j=1:length(seq2)
                Ampl_Hevents(i,j)=max(Fnew(i,seq2{j}));
            end
        end
        
        
        seq2=seq(short2);
        for i=1:size(Fnew,1)
            for j=1:length(seq2)
                Ampl_Levents(i,j)=max(Fnew(i,seq2{j}));
            end
        end
        
        if long2>0
            Ampl_Hevents=mean(Ampl_Hevents,2);
            Ampl_Levents=mean(Ampl_Levents,2);
            
            nn=1:size(Fnew,1);
            novl=setdiff(nn,idd2);
            
            
            Ampl_HeventsSST=Ampl_Hevents(idd2);
            Ampl_HeventsOther=Ampl_Hevents(novl);
            
            Ampl_LeventsSST=Ampl_Levents(idd2);
            Ampl_LeventsOther=Ampl_Levents(novl);
            
            length_Hevents_sum{kk}=length_Hevents;
            length_Levents_sum{kk}= length_Levents;
            participation_Hevents_sum{kk}= participation_Hevents;
            participation_Levents_sum{kk}=participation_Levents;
            Ampl_Hevents_sum{kk}=Ampl_Hevents;
            Ampl_Levents_sum{kk}=Ampl_Levents;
            Ampl_HeventsSST_sum{kk}=Ampl_HeventsSST;
            Ampl_HeventsOther_sum{kk}=Ampl_HeventsOther;
            Ampl_LeventsSST_sum{kk}=Ampl_LeventsSST;
            Ampl_LeventsOther_sum{kk}=Ampl_LeventsOther;
        else
        end
    catch
    end
end

cd(BasePath);

 length_Hevents_sum=cell2mat(length_Hevents_sum);
 length_Levents_sum=cell2mat(length_Levents_sum);
 participation_Hevents_sum=cell2mat(participation_Hevents_sum);
 participation_Levents_sum=cell2mat(participation_Levents_sum);
 Ampl_Hevents_sum=cell2mat(Ampl_Hevents_sum');
 Ampl_Levents_sum=cell2mat(Ampl_Levents_sum');
 Ampl_HeventsSST_sum=cell2mat(Ampl_HeventsSST_sum');
 Ampl_HeventsOther_sum=cell2mat( Ampl_HeventsOther_sum');
 Ampl_LeventsSST_sum=cell2mat(Ampl_LeventsSST_sum');
 Ampl_LeventsOther_sum=cell2mat(Ampl_LeventsOther_sum');
 
 clearvars -except length_Hevents_sum length_Levents_sum participation_Hevents_sum participation_Levents_sum Ampl_Hevents_sum Ampl_Levents_sum Ampl_HeventsSST_sum Ampl_HeventsOther_sum Ampl_LeventsSST_sum  Ampl_LeventsOther_sum