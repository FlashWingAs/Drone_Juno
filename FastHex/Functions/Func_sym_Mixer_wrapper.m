function Mixer = Func_sym_Mixer_wrapper(alpha, dead_zone, Mixer_flat)

if alpha >= dead_zone
    Mixer = AutoFunc_sym_Mixer(alpha);
else
    Mixer = Mixer_flat;
end

end