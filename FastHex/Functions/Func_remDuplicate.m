function Y = Func_remDuplicate(X)
%FUNC_REM Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    X (:,:) {mustBeNumeric};
end

arguments (Output)
    Y;
end

[n,m] = size(X);
eps = 1e-5;
indicesToRemove = false(1, size(X, 2));

for i = 1:m-1
    for j = i+1:m
        k = 1;
        flag = true;
        while flag
            if ~Func_equal_eps(X(k,i), X(k,j), eps)
                flag = false;
            else
                if k == n
                    indicesToRemove(j) = true;
                    flag = false;
                else
                    k = k+1;
                end
            end
        end

    end
end

Y = X(:, ~indicesToRemove);

end