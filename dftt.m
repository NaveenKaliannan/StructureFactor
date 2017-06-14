function [Xk] = dftt(g_r, x, number_density)
dr = x(2) - x(1);
Xk = [];
for q = x(1):x(2)-x(1):x(end)
    sum = 0;
    incre = 0;
    for r = x(1):x(2)-x(1):x(end)
        incre = incre + 1;
        sum = sum +  r * (g_r(incre) - 1) * sin(q * r);
    end
    Xk = [Xk, 1 + 4 * pi * number_density * sum * dr / q];
end

