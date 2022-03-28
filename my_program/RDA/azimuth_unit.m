function [x, y] = azimuth_unit(Z,y_a,x_r)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
max_val = max(max(abs(Z)));
[a_max,r_max]=find(abs(Z)==max_val);
q=[];
for i=1:size(Z,1)
    if abs(Z(i,r_max))>=(2^(-0.5))*max_val
        q=[q,y_a(i)];
    end
end
z=[];
for i=1:size(Z,2)
    if abs(Z(a_max,i))>=(2^(-0.5))*max_val
        z=[z,x_r(i)];
    end
x = max(q)-min(q);
y = max(z)-min(z);
end

