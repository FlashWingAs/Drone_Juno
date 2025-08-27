function x0 = Func_intersect_parametrized(A1,B1, A2, B2)
%FUNC_INTERSECT_PARAMETRIZED Returns the interesection coordinates between 2 straight lines
%   Detailed explanation goes here
arguments (Input)
    A1 (:,1) {mustBeNumeric}
    B1 (:,1) {mustBeNumeric}
    A2 (:,1) {mustBeNumeric}
    B2 (:,1) {mustBeNumeric}
end

arguments (Output)
    x0
end

% A1*t + B1 == A2*t + B2

% if A1 ~= A2
%     t0 = (A1-A2)\(B2-B1);
%     x0 = A1*t0+B1;
% else
%     x0 = Inf([size(B1)]);
% end

% a1*t+c1 = a2*t+c2
% b1*t+d1 = b2*t+d2

a1 = A1(1);
b1 = A1(2);
c1 = B1(1);
d1 = B1(2);
a2 = A2(1);
b2 = A2(2);
c2 = B2(1);
d2 = B2(2);

eps = 1e-5;

if ~Func_equal_eps(a1, a2, eps)
    t0 = (c2-c1)/(a1-a2);
    x0 = A1*t0+B1;
elseif ~Func_equal_eps(b1, b2, eps)
    t0 = (d2-d1)/(b1-b2);
    x0 = A2*t0+B2;
else
    x0 = [Inf; Inf];
end