function [t0] = Func_root_parmLine(A, B, X0)

arguments (Input)
    A (:,1) {mustBeNumeric}
    B (:,1) {mustBeNumeric}
    X0 (:,1) {mustBeNumeric}
end

arguments (Output)
    t0
end

% a*t+c = x0
% b*t+d = y0

a = A(1);
b = A(2);
c = B(1);
d = B(2);

x0 = X0(1);
y0 = X0(2);

eps = 1e-5;

if ~Func_equal_eps(a, 0, eps)
    t0 = (x0 - c)/a;
elseif ~Func_equal_eps(b, 0, eps)
    t0 = (y0 - d)/b;
else
    t0 = Inf;
end


end