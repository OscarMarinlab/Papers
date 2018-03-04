close all; 
clear all; 
%%

i = dir('*.tif'); %lists all stk stacks in current directory
px = 1.1033; %pixel to um conversion - check on apotome
ch1 = 1; % reporter channel


Brainno = '130516_B2';
Genotype = 'vgat-cre';
Age = 'E14.5';
A= struct('Brainno', {Brainno}, 'Genotype', {Genotype}, 'Age', {Age}, 'ch', {ch1}, 'px', {px});
outputfile = strcat (Brainno, '_', Genotype, '.csv'); 

%%

for k = 1:numel(i) % loops through all images in the directory
    
    filename = i(k).name;
    z=0;
    z1 =0; 
    while z1 ==0; 
        vgat((2*k-1:2*k), :) = vgat_e16 (filename, A); 
        questions = [{'Redo same image? Yes = 1, No = 2'}];
        default_answers = [{'2'}];
        answer = inputdlg(questions,'Redo',1,default_answers);
        close all; 
        a = str2double(answer(1)); 
        if a ==2
            z1 = 1;
        else
            z1 = 0;
        end
        
    end
end

export (vgat, 'File', outputfile,'Delimiter',','); 

