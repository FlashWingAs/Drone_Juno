function comb = Func_combinaisons(vec,rep)
if rep < 2
    comb = vec;
else
    n = numel(vec);
    c = n^rep;
    comb = zeros([rep, c]);
    for i = 1:rep
        for k = 1:n^(i-1)
            for j = 1:n
                start = (j-1+(k-1)*n)*n^(rep-i)+1;
                stop = (j+(k-1)*n)*n^(rep-i);
                % disp("-----------------")
                % disp("i="+num2str(i))
                % disp("k="+num2str(k));
                % disp("j="+num2str(j));
                % disp("start="+num2str(start));
                % disp("stop="+num2str(stop));
                % disp("-----------------")
                comb(i,start:stop) = vec(j);
            end
        end
    end
end

end