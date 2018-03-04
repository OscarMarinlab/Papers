
close all, clear all; 

i = dir ('*pA.png'); 
%%
Genotype = 'Het'; 
data = zeros([],5); 
data = mat2dataset(data, 'VarNames', {'MeanIPSC', 'Row', 'cellid', 'Genotype', 'Layer'}); 
%%
for k = 1:numel(i)
    FileName = i(k).name; 
    %%
    Layer = FileName(1:2); 
    cellid = strrep(FileName(4:end),'_pA.png',''); 
    M = csvread(strcat(cellid, '_pA.csv')); 
    C_avg = mean(M);
    Th = C_avg/mean(C_avg); 
    %%
    % set threshold, column average has to be greater than 15% mean IPSC
    B = Th >0.15; 
    SubM = M(:, B); % subset matrix
    %%
    Row_avg = mean(SubM, 2); 
    Row_avg(:,2) = (1:6);
    Row_avg = mat2dataset(Row_avg, 'VarNames', {'MeanIPSC', 'Row'}); 
    Row_avg.cellid = repmat(cellid, 6, 1);
    Row_avg.Genotype = repmat(Genotype, 6, 1);
    Row_avg.Layer = repmat(Layer, 6, 1);
   
    %%
    data = vertcat(data, Row_avg);
    
end
export(data,'File', strcat(Genotype,'15psCRACM_avg.csv'),'Delimiter',',')

