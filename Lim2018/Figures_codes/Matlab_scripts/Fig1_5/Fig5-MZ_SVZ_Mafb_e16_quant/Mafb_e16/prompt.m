z=0; 
while z == 0; 
    questions = [{'Redo GFP segment? Yes = 1, No = 2'}];
    default_answers = [{'2'}];
    answer = inputdlg(questions,'GFP',1,default_answers);
    a = str2double(answer(1));
    if a == 2
        z=1; 
    else
        z=0;
    end
end 
