function [p1, I1] = comparepeak(cbg, j, b, a)
    x = std (cbg(1:5000, j));
    f2 = filter(b, a, cbg(:, j)); 
    [p1, I1] = findpeaks (f2(5100:6600), 'MinPeakHeight', 3*x);
end
    
    
     