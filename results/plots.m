clear all; close all; clc;

%% Eclipse 10x10
p_10x10 = load('eclipse/10x10-pressure.dat')
for x=1:10
    for y=1:10
        p_diag = p_10x10(x,y);
    end
end