function bool = Func_equal_eps(var1,var2, eps)

bool = var1 >= var2 - eps && var1 <= var2 + eps;

end