function onTheRight = Func_RightOrLeft(Vec1,Vec2)

arguments (Input)
    Vec1 (:,1) {mustBeNumeric}
    Vec2 (:,1) {mustBeNumeric}
end

arguments (Output)
    onTheRight
end

R = [0,1;-1,0];

onTheRight = dot(R*Vec1/norm(Vec1), Vec2/norm(Vec2)) > 0;

end