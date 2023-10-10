clear,close all
% Create some data
materials.name= {'Cork','Balsa Wood','Aerogels','Metal Foams','Polymer Foams','this work','Rubber',...
    'this work_steel','this work_al'};
materials.densityRange = [200 300; 100 200; 1 15; 200 800; 10 100;  42 430; 0.8e3 3e3;...
    [42 430].*(7850/1200);  [42 430].*(2700/1200);  ];  % in kg/m^3
materials.Shear_modulusRange= [0.005 0.020;   0.05 0.15;   0.1e-3 0.01;  0.01  100; ...
    0.001 0.100; 0.0164e-3 0.49e-3; 0.1e-3 10e-3;...
    [0.0164e-3 0.49e-3].*(210/2.7) ; [0.0164e-3 0.49e-3].*(70/2.7);];  % in Gpa

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
    Position=[materials.densityRange(i,1)  materials.Shear_modulusRange(i,1),...
        materials.densityRange(i,2)-materials.densityRange(i,1),...
        materials.Shear_modulusRange(i,2)-materials.Shear_modulusRange(i,1)];
    
    rectangle('Position',Position,'EdgeColor',facecolor{i},...
        'LineWidth',1,'Curvature',[0.5,1],'FaceColor',facecolor{i});
%     text( Position(1)+ Position(3)/2, Position(2)+ Position(4)/6, ...
%         materials.name{i}, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
%     
    hold on
end
% Create a scatter plot
% scatter(density, bulk_modulus, 'filled');
pentamodeShear = linspace(1e3*1e-9, 1e5*1e-9, 10);
pentamodeDensity = linspace(500, 1200, 10);

Shear(1).property = linspace(1e3*1e-9, 1e5*1e-9, 10);
Shear(1).name='pentamode';
Shear(1).Density=linspace(500, 1200, 10);

Shear(2).property = linspace(0.1*77, 0.45*77, 10);
Shear(2).name='Plate lattice';
Shear(2).Density=linspace(0.05*7890, 0.4*7890, 10);

Shear(3).property = linspace(0.15, 0.45, 10);
Shear(3).name='gear';
Shear(3).Density=linspace(0.4*8900, 0.5*8900, 10);

Shear(4).property =[0.0675e-3, 0.24e-3];
Shear(4).name='sphere2Sphere';
Shear(4).Density=[ 62.5 ,320];

Shear(5).property =[0.0164e-3, 0.48e-3] ;
Shear(5).name='sphere2cube';
Shear(5).Density=[187, 62.5];


% To plot it on a logarithmic scale, use semilogx or semilogy or loglog
% figure
markerColor{1}=[43, 189, 87]./256;
markerColor{2}=[43, 118, 189]./256;
markerColor{3}=[42, 42, 128]./256;
markerColor{4}=[196, 16, 94]./256;
markerColor{5}=[196, 16, 94]./256;
for i=1:length(Shear)
scatter(Shear(i).Density,Shear(i).property,'filled','MarkerEdgeColor',[101, 126, 150]./256, 'MarkerFaceColor',markerColor{i},'MarkerFaceAlpha',0.3)
hold on
end
% Set the scales to logarithmic
yticks([ 1e-5 1e-3 1e-1 1e1])
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([ 0.5e0  1e4])
% Annotate the points
grid on;

% Add titles and labels
% title('Ashby Chart: Bulk Modulus vs Density');
xlabel('Density (kg/m^3)');
ylabel('Shear Modulus (GPa)');

fig = figure(1);
% Define the figure size [left, bottom, width, height]
figureSize = [646.6,267.4,597.6,438.4];  %  setting after 9/15
% Set the figure size
set(fig, 'Position', figureSize);
% Get current axes handle
ax = gca;
% Set font size of tick labels
ax.FontSize = 13;  % Adjust the font size as needed