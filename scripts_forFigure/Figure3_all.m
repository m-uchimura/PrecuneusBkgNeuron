%%
%this code draw PSTH related figures (Figure3 and Figure 13)
%%
clear 
close all
addpath(genpath('../toolboxes'))
%%
Figure3_psth_raw  %Figure 3A,B and Figure 13E-G
Figure3_gaze      %Figure 3C
Figure3_histogram %Figure 3D,E,F

%%
if 0
cd figs
figure(1)
print('fig3_mean.ai', '-dpdf', '-painters','-bestfit')
figure(2)
print('fig13_psth_mean.ai', '-dpdf', '-painters','-bestfit')
cd ../
end