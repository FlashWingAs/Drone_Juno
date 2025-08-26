function [Mixer, Mixer_flat, Mixer_augmented] = Func_sym_Mixer_wrapper(alpha)

if alpha <= pi/4
    [Mixer, Mixer_flat, ~] = AutoFunc_sym_Mixer(alpha);
else
    [Mixer, ~, Mixer_flat] = AutoFunc_sym_Mixer(alpha);
end
Mixer_augmented = [Mixer; zeros([2,8])];
end