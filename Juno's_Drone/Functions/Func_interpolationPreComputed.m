function [kFmot] = Func_interpolationPreComputed(alphaValue, yawValue, alphaMesh, yawMesh, kFmot_memory, deadzone_alpha)
%FUNC_INTERPOLATIONPRECOMPUTED Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    alphaValue (1,1) {mustBeNumeric}
    yawValue   (1,1) {mustBeNumeric}
    alphaMesh  (:,:) {mustBeNumeric}
    yawMesh    (:,:) {mustBeNumeric}
    kFmot_memory (:,:,2,8) {mustBeNumeric}
    deadzone_alpha (1,1) {mustBeNumeric}
end

arguments (Output)
    kFmot
end

if alphaValue < 0 + deadzone_alpha
    alphaValue = 0 + deadzone_alpha;
elseif alphaValue > pi/2 - deadzone_alpha
    alphaValue = pi/2 - deadzone_alpha;
end

yawValue = wrapToPi(yawValue);


meshSize = size(yawMesh);
kFmot = zeros(2,8);

[~, iAlpha] = min(abs(alphaMesh-alphaValue));

for i = 1:2
    for j = 1:8
        kFmot(i, j) = interp1(yawMesh, reshape(kFmot_memory(:, iAlpha, i, j), meshSize), yawValue, "nearest");
    end
end

end