function [x, y] = azimuth_unit(Z,y_a,x_r)
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
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

