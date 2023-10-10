clear, close all
addpath(genpath(pwd));

% ------------------Select model ------------------1
% Case='Sphere'; %333 Sphere
% Case='Hyperboloid'; %333 Spindle
Case='Cone'; %333 Cone

% -----------------------------------------------------
load([Case,'.mat']) 
viewPoint=[-25 18];
figure 
Plot(Ntub,viewPoint,'texxt','ftexxt','Snapology',4,1,'on',[])
axis off
axis equal
figure
Plot(Ntub2D,viewPoint,'texxt','ftexxt','Snapology',4,1,'on',[])
