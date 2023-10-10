clear,close all
% Create some data
materials.name= {'Cork','Balsa Wood','Aerogels','Metal Foams','Polymer Foams','this work', 'Rubber',...
    'this work_steel','this work_al'};
materials.densityRange = [200 300; 100 200; 1 15; 200 800; 10 100;  45 420; 0.8e3 1.5e3;...
    [45 420].*(7850/1200);  [45 420].*(2700/1200); ];  % in kg/m^3
materials.bulk_modulusRange= [0.02 0.05; 0.01 0.05; 0.01 0.1; 0.1 1; 0.01 0.1; 0.25e-3 10e-3;...
    1e-3 2000e-3; ...
    [0.25e-3 10e-3].*(210/2.7);  [0.25e-3 10e-3].*(70/2.7);];  % in Gpa

facecolor{1}=[255 240 245]./256;
facecolor{2}=[255 193 193]./256;
facecolor{3}=[250 240 230]./256;
% facecolor{4}=[240 248 255]./256;
facecolor{4}=[230, 240, 255]./256;
facecolor{5}=[238 224 229]./256;
facecolor{6}=[238 200 229]./256;
% facecolor{7}=[230, 255, 238]./256;

 facecolor{7}=[214, 245, 214]./256;
facecolor{8}=[238 200 229]./256;
facecolor{9}=[238 200 229]./256;
% Create a figure
figure
% Plot a rectangle
% 'Position' vector has the form [x y width height]
for i=1:length(materials.name)-2
    Position=[materials.densityRange(i,1)  materials.bulk_modulusRange(i,1),...
        materials.densityRange(i,2)-materials.densityRange(i,1),...
        materials.bulk_modulusRange(i,2)-materials.bulk_modulusRange(i,1)];
    
    rectangle('Position',Position,'EdgeColor',facecolor{i},...
        'LineWidth',1,'Curvature',[0.5,1],'FaceColor',facecolor{i});
%     text( Position(1)+ Position(3)/2, Position(2)+ Position(4)/6, ...
%         materials.name{i}, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
%     
    hold on
end

% Create a scatter plot
% scatter(density, bulk_modulus, 'filled');
Bulk(1).property = linspace(1e6*1e-9, 1e7*1e-9, 10);
Bulk(1).name='pentamode';
Bulk(1).Density=linspace(500, 1200, 10);

Bulk(2).property = linspace(0.1*143, 0.45*143, 10);
Bulk(2).name='Plate lattice';
Bulk(2).Density=linspace(0.05*7890, 0.4*7890, 10);

Bulk(3).property = linspace(0.05*5, 0.35*5, 10);
Bulk(3).name='pentamode 2';
Bulk(3).Density=linspace(0.05*1190, 0.6*1190, 10);

Bulk(4).property =[0.996e-3, 0.36e-3];
Bulk(4).name='sphere2sphere';
Bulk(4).Density= [ 62.5 ,320];


Bulk(5).property =[1.846e-3, 8.23e-3] ;
Bulk(5).name='sphere2cube';
Bulk(5).Density=[187, 62.5];


markerColor{1}=[43, 189, 87]./256;
markerColor{2}=[43, 118, 189]./256;
markerColor{3}=[42, 12, 128]./256;

markerColor{4}=[196, 16, 94]./256;
markerColor{5}=[196, 16, 94]./256;
% To plot it on a logarithmic scale, use semilogx or semilogy or loglog
% figure
for i=1:length(Bulk)
scatter(Bulk(i).Density,Bulk(i).property,'filled','MarkerEdgeColor',[101, 126, 150]./256, 'MarkerFaceColor',markerColor{i},'MarkerFaceAlpha',0.3)
hold on
end
% Set the scales to logarithmic
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([ 0.5e0  1e4])
% Annotate the points
grid on;

% Add titles and labels
% title('Ashby Chart: Bulk Modulus vs Density');
xlabel('Density ( kg/m^3)');
ylabel('Bulk Modulus (GPa)');
fig = figure(1);
% Define the figure size [left, bottom, width, height]
% figureSize = [646.6,267.4,497.6,238.4];  %  setting before 9/15
figureSize = [646.6,267.4,597.6,438.4];  %  setting after 9/15
% Set the figure size
set(fig, 'Position', figureSize);
% Get current axes handle
ax = gca;
% Set font size of tick labels
ax.FontSize = 13;  % Adjust the font size as needed