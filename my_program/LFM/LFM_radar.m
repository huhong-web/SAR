function [T,B,Rmin,Rmax,R,RCS] = LFM_radar(~,~)
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
if nargin==0
    T = 10e-6;
    B = 30e6;
    Rmin = 10000;Rmax=15000;
    R = [11000,13000,13002];
    RCS = [1 1 1];
end

