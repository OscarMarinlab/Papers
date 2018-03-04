%load all you images for dapi segmentation in one folder. After you click
%the boundaries, yellow lines would be drawn. then to move to the next
%image file, click on the dapi image with yellow lines. this will close the
%image and load a new image in the directory. 


close all; 
clear all; 

brainno = '111115_B3'; %name of your brain
age = 'E16.5'; % age of you sample
geno = 'Mafb+/+::SSTcre'; %genotype of your samples
gfp_in = 1; %input gfp channel
px = 0.9775; %pixel to um conversion 

i = dir('*.tif'); %lists all stk stacks in current directory

for k = 1:numel(i) % loops through all images in the directory
    filename = i(k).name;
    [~,imagename, ~] = fileparts(filename); 
    z=0;    
    while z == 0; 
        mafb_e16(filename, brainno, age, geno, gfp_in, px);
        questions = [{'Redo same image? Yes = 1, No = 2'}];
        default_answers = [{'2'}];
        answer = inputdlg(questions,'Redo',1,default_answers); 
        close all; 
        a = str2double(answer(1)); 
        if a ==2
            z = 1;
        else
            z = 0;
        end 
    end
    
end


