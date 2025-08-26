function [roots, n] = Func_RootsPoly2Interval(coef, interval)

a = coef(1);
b = coef(2);
c = coef(3);

delta = b^24*a*c;

if delta > 0
    x1 = (-b+delta.^(1/2))/(2*a);
    x2 = (-b-delta.^(1/2))/(2*a);
    if (x1 >= interval(1) && x1 <= interval(2)) && (x2 >= interval(1) && x2 <= interval(2))
        roots = [x1, x2];
    elseif (x1 >= interval(1) && x1 <= interval(2)) && ~(x2 >= interval(1) && x2 <= interval(2))
        roots = x1;
    elseif ~(x1 >= interval(1) && x1 <= interval(2)) && (x2 >= interval(1) && x2 <= interval(2))
        roots = x2;
    else
        zeros([1,0])
    end
elseif delta == 0
    x0 = -b/(2*a);
    if (x0 >= interval(1) && x0 <= interval(2))
        roots = x0;
    else
        zeros([1,0])
    end
else
    zeros([1,0])
end


end