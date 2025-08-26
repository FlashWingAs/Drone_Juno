function Progress_Log(done,text, LEN)

if ~done
    fprintf(text)
    n = size(text, 2);
    if n + 4 >= LEN
        LEN = n + 1;
    end
    for k = 1 : LEN - n
        fprintf('.')
    end
    fprintf('\n')
else
    fprintf('\b\b\b\b\bDONE\n')
end