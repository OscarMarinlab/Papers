function ft = sCRACM_v2 (current, c)
current2 = current(:, 2:127); 
b= current(1:5000, 2:127);
m1 = median(b);
m2 = ones(18000, 1);
m2 = m2*m1; 
cbg= current2 -m2; 
%%
%filtering low pass
fc = 150; %cut off frequency
fs = 20000; %sampling rates in hz
[b, a ] = butter (6, fc/(fs/2)); % don't touch this 
 
%%
% extract peaks using the comparepeak function
[mx1(:, 1), mx1(:, 2)] = comparepeak(cbg, c, b, a);
mx1(:, 3) = c; 
[mx2(:, 1), mx2(:, 2)] = comparepeak(cbg, c+42, b, a);
mx2(:, 3) = c+42; 
[mx3(:, 1), mx3(:, 2)] = comparepeak(cbg, c+42+42, b, a);
mx3(:, 3) = c+84; 
mt = vertcat(mx1, mx2, mx3); 
%%

%%
Imx1 = mx1(:, 2); 
Imx2 = mx2(:, 2);
Imx3 = mx3(:, 2);
Ind = vertcat(Imx1, Imx2, Imx3); 

%%
if size (Ind, 1) > 0; 
    %for 1st value; 
    k=1; 
    i1 = Ind(k); 
    j1 = mt(:,2)>i1-50 & mt(:,2)<i1+50; 
    extall = mt(j1, :);  
    extall(:, 4) = k; 
    %%
        if size (Ind, 1)>1; 
            for k = 2: size(Ind, 1)
                l = size(extall, 1); 
                i1 = Ind(k); 
                j1 = mt(:,2)>i1-50 & mt(:,2)<i1+50; 
                extall2 = mt(j1, :); 
                extall2(:, 4) = k; 
                extall = vertcat(extall, extall2); 
            end
        end
        %%
        %convert extall to dataframe
        ds = array2table(extall, 'VariableNames', {'Amp', 'time', 'spotid', 'peakno'}); 
        [~, a1] = unique (ds(:, 1));
        ds = ds(a1, :); 
        [n, bin] = histc(ds.peakno, unique(ds.peakno));
        multiple = find (n>1); 
        if multiple >0; 
            I4= find(ismember(bin, multiple)); 
            ds = ds(I4, :); 
            dsmeanamp = varfun(@mean, ds, 'InputVariables', {'Amp',},'GroupingVariables', {'peakno'});
            dsV1amp = varfun(@stdovermean, ds, 'InputVariables', {'Amp',},'GroupingVariables', {'peakno'});
            ct = join (dsmeanamp, dsV1amp); 
            min_var = min(ct.stdovermean_Amp); 
            ft = ct(ct.stdovermean_Amp==min_var, :); 
            ft.spot = c; 
            ft = ft(:, {'mean_Amp','spot'}); 
        else
            ft = [0, c];
            ft = array2table(ft, 'VariableNames', {'mean_Amp','spot'});
        end
else
    ft = [0, c]; 
    ft = array2table(ft, 'VariableNames', {'mean_Amp','spot'});
    
end
end

