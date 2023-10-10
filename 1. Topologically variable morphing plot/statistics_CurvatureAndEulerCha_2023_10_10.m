% Curvature chart
clear, close all
% addpath('Statistics of Gassuan Curvature\All 2D morphing');
% auto get curvatures of specified STL file

%% Get a list of all STL files in the current folder
filename = 'allcurvature.mat';

if exist(filename, 'file') == 2
    disp('File exists')
    load([filename])
else
    disp('File does not exist')
    files = dir('*.stl');
    getderivatives=0;
    % Get a list of all STL files in the current folder
    files = dir('*.stl');
    
    % Group all STL files according to their names
    % Create a Map to store the file names for each set of letters
    fileGroups = containers.Map();
    % Loop over all files
    for i = 1:length(files)
        % Get the file name
        filename = files(i).name;
        % Remove the '.stl' extension
        filename_no_ext = erase(filename, '.stl');
        % Extract letters from the file name
        letters = regexp(filename_no_ext, '[a-zA-Z]', 'match');
        % Concatenate the letters into a single string
        letters_str = strjoin(letters, '');
        % If this set of letters is not yet in the Map, add it
        if ~isKey(fileGroups, letters_str)
            fileGroups(letters_str) = {filename};
        else
            % If this set of letters is already in the Map, append the file name
            fileGroups(letters_str) = [fileGroups(letters_str) {filename}];
        end
    end
    % Get all the keys from the Map
    % Get all the keys from the Map
    letters = keys(fileGroups);
    %% Calculate Curvatures
    % Loop over the keys
    for i = 1:length(letters)
        % Get the current key
        current_letter = letters{i};
        % Get the files for the current key
        filenames = fileGroups(current_letter);
        % Display the key and the corresponding files
        fprintf('Letters: %s\n', current_letter);
        % Loop over the files
        for j = 1:length(filenames)
            FV=[]; PrincipalCurvature=[];
            STLfile=stlread(filenames{j});
            FV.vertices=STLfile.Points;
            FV.faces=STLfile.ConnectivityList;
            % Now you can do something with the vertices and faces...
            [PrincipalCurvatures,PrincipalDir1,PrincipalDir2,FaceCMatrix,VertexCMatrix,Cmagnitude]= GetCurvatures( FV ,getderivatives);
            [Pc1, Pc2]=filterXY(PrincipalCurvatures(1,:),PrincipalCurvatures(2,:));
            GausianCurvature=Pc1.*Pc2;
            % get Gaussian Curvature
            GausianCurvature_fliter =filter(GausianCurvature);
            rangeGausianCurvature=[min(GausianCurvature_fliter) max(GausianCurvature_fliter)];
            % get Mean Curvature
            meanCurvature=0.5*(Pc1+Pc2);
            meanCurvature_fliter =filter(meanCurvature);
            rangemeanCurvature=[min(meanCurvature_fliter) max(meanCurvature_fliter)];
            ref(i).GausianCurvature_fliter{j}=GausianCurvature_fliter; %
            ref(i).meanCurvature_fliter{j}=meanCurvature_fliter; %
            ref(i).PrincipalCurvatures{j}=[Pc1; Pc2;]; %
            ref(i).GausianMeanCurvatures{j}=[GausianCurvature; meanCurvature;]; %
            ref(i).curvatureRange(j,:)=[rangeGausianCurvature , rangemeanCurvature]; % K1 K2 H1 H2
            ref(i).Name{j}=filenames{j};
        end
    end
    %     save(filename,'ref')
end

%% plot range of curvatures in x=1 plane
% Create a new figure
fig=figure;
view([44,27]);
k=0;
 for i= 1: length(ref) 
        for j=1:length(ref(i).GausianMeanCurvatures)
        x=[]; y=[]; log10x=[]; log10y=[];
        %     for j=1:length( ref(i).PrincipalCurvatures)
        % Define the rectangle parameters
        x = ref(i).GausianMeanCurvatures{j}(1,:); % x-coordinates
        y=ref(i).GausianMeanCurvatures{j}(2,:); % y-coordinates
        %                 x = ref(i).PrincipalCurvatures{j}(1,:); % x-coordinates
        %                 y =ref(i).PrincipalCurvatures{j}(2,:); % y-coordinates
        x=abs(x);  y=abs(y);
        filename_no_ext = erase(ref(i).Name{j}, '.stl');
        name =strjoin(regexp(filename_no_ext, '[a-zA-Z]', 'match'), '');
        if strcmp(name,'kai')
            planeIdx=2;
            facecolor=findColor2(j,planeIdx,5);
        elseif strcmp(name,'Xo')
            planeIdx=0;
            facecolor=findColor2(j,planeIdx,7);
        else
            planeIdx=1;
            facecolor=findColor2(k,planeIdx,5);
        end
        z = planeIdx* ones(1,length(x)); % z-coordinates
        %         facecolor=findColor2(i,14);
        
        P=[x; y]';
        %         loglog(x,y,'.')
        for kl=1:length(x)
            log10x(kl)=x(kl)/abs(x(kl))*log10(abs(x(kl)));
            log10y(kl)=y(kl)/abs(y(kl))*log10(abs(y(kl)));
        end
        [log10x,log10y]=excludeInfNaN(log10x,log10y);
        scatter3(x, y, z,'MarkerFaceColor',facecolor ,'MarkerEdgeColor',facecolor,...
            'MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.3,'SizeData',10);
        set(gca, 'XScale', 'log', 'YScale', 'log');
        K = boundary(log10x', log10y',0);
        hold on
        set(gca, 'XScale', 'log', 'YScale', 'log');
        fill3(x(K)', y(K)', z(K)', facecolor, 'FaceAlpha', 0.2, 'EdgeColor',facecolor, 'EdgeAlpha',0.3,'LineWidth',1);
        hold on
        %         if   planeIdx==1
        %         k=k+1;
        %         end
    end
    if   planeIdx==1
        k=k+1;
    end
end
% Set the view angle
% coordinates of the rectangle
% Patch XYZ planes
xlimValuelow=1e-9;   ylimValuelow=1e-9;
xlimValuehigh=1e2;   ylimValuehigh=1e2;
% xlimValuehigh=1e0;   ylimValuehigh=1e0;
% x = [xlimValuelow  xlimValuehigh  xlimValuehigh  xlimValuelow];
xticks([1e-9 1e-5 1e-1])
% y = [ylimValuelow  ylimValuelow   ylimValuehigh  ylimValuehigh];
yticks([1e-9 1e-5 1e-1])
% create the rectangle
% trans=0.06;
% patch('XData',x,'YData',y,'ZData',[1 1 1 1],'FaceColor',[162  62 125]./256, 'EdgeColor',[162  62 125]./256,'LineWidth',2,'faceAlpha',trans,'edgeAlpha',trans/10);
% % coordinates of the rectangle
% patch('XData',x,'YData',y,'ZData',[2 2 2 2],'FaceColor',[255,218,185]./256, 'EdgeColor',[255,218,185]./256,'LineWidth',2,'faceAlpha',trans*2,'edgeAlpha',trans/10);
% % coordinates of the rectangle
% patch('XData',x,'YData',y,'ZData',[0 0 0 0],'FaceColor',[62 162 102]./256, 'EdgeColor',[62 162 102]./256,'LineWidth',2,'faceAlpha',trans,'edgeAlpha',trans/10);
set(fig, 'Position', [400, 300, 470, 450]);
% view([0,90]);
view([44,23]);
grid on;
hold on
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([xlimValuelow   xlimValuehigh])
ylim([ylimValuelow   ylimValuehigh])
zlim([0 2])
zticks([0 1 2])
% Label the axes
xlabel('$|H| (mm^{-1})$', 'Interpreter', 'latex','FontSize', 14);
ylabel('$|K| (mm^{-2})$', 'Interpreter', 'latex','FontSize', 14);
zlabel('\chi','FontSize', 14);
ax = gca;  % get current axes
ax.XAxis.FontSize = 9;  % Set x-axis tick label font size
ax.YAxis.FontSize = 9;  % Set y-axis tick label font size
ax.ZAxis.FontSize = 9;  % Set y-axis tick label font size
% axis equal

% ----------------------------------------------------------------------------------
% ------------------------------ Sub- Function ----------------------------------
function [x,y]=excludeInfNaN(x,y)
% create an array with finite and non-finite elements
% create a logical array that is true for each finite element of A
finiteIdx1 = find(isfinite(x)==1);
finiteIdx2 = find(isfinite(y)==1);
% use the logical array to index into A and create a new array B that excludes non-finite elements
x=x(unique([finiteIdx1 finiteIdx2]));
y=y(unique([finiteIdx1 finiteIdx2]));
end

function facecolor=findColor2(stlIdx,planeIdx,numUnitPoly)
switch  planeIdx
    case 2
        switch mod(stlIdx,numUnitPoly)
            case 0
                facecolor=[255,99,71]./256;   % original
            case 1
                facecolor=[255,127,80]./256;    % original
            case 2
                facecolor=[205,92,92]./256;
            case 3
                facecolor=[240,128,128]./256;
            case 4
                facecolor=[233,150,122]./256;
            case 5
                facecolor=[250,128,114]./256;
        end
    case 1   %  E.C.= 1 plane
        switch mod(stlIdx,numUnitPoly)
            case 0  % baiyun
                facecolor=[0,128,0]./256;  % green
            case 1  % yihui
                facecolor=[128,0,128]./256;   % original   purple
            case 2 % elastom
                facecolor=[244,164,96]./256;    % original    % orange
            case 3  % 4D
                facecolor=[0,0,128]./256;  % blue
            case 4  % maha
                facecolor=[0,128,128]./256;  % c
            case 5
                facecolor=[100,149,237]./256;
            case 6
                facecolor=[192,192,192]./256;
            case 7
                facecolor=[128,128,128]./256;
        end
    case 0  %  E.C.= 1 plane
        switch mod(stlIdx,numUnitPoly)
            case 0
                facecolor=[255,99,71]./256;   % original  % orange
            case 1
                facecolor=[255,127,80]./256;    % original    % gREEN
            case 2
                facecolor=[205,92,92]./256;   %
            case 3
                facecolor=[240,128,128]./256;
            case 4
                facecolor=[233,150,122]./256;
            case 5
                facecolor=[250,128,114]./256;
        end
end
end



function facecolor=findColor(i,numUnitPoly)
switch mod(i,numUnitPoly)
    case 0
        facecolor=[255 240 245]./256;   % original
    case 1
        %                   facecolor=[250 240 230]./256;
        facecolor=[255 193 193]./256;    % original
    case 2
        facecolor=[250 240 230]./256;
    case 3
        facecolor=[240 248 255]./256;
    case 4
        facecolor=[238 224 229]./256;
    case 5
        facecolor=[238 200 229]./256;
    case 6
        facecolor=[245 200 210]./256;
    case 7
        facecolor=[245 220 210]./256;
    case 8
        facecolor=[255 193 193]./256;
    case 9
        facecolor=[250 240 230]./256;
    case 10
        facecolor=[240 248 255]./256;
    case 11
        facecolor=[238 200 229]./256;
    case 12
        facecolor=[238 224 229]./256;
    case 13
        facecolor=[245 220 210]./256;
    case 14
        facecolor=[245 200 210]./256;
end
end

function data_filtered =filter(data)
mu = mean(data);
sigma = std(data);

% 计算离群值的阈值
threshold_high = mu + 3*sigma;
threshold_low = mu - 3*sigma;

% 找到在这个范围内的数据
data_filtered = data(data > threshold_low & data < threshold_high);
end

function [x_filtered, y_filtered]=filterXY(x, y)  %  this uses interpolation method

xQ1 = quantile(x, 0.25);   yQ1 = quantile(y, 0.25);
xQ3 = quantile(x, 0.75);   yQ3 = quantile(y, 0.75);

% Calculate the IQR
xIQR = xQ3 - xQ1;   yIQR = yQ3 - yQ1;

% Determine the threshold for outliers
xlowerBound = xQ1 - 1.5*xIQR;     ylowerBound = yQ1 - 1.5*yIQR;
xupperBound = xQ3 + 1.5*xIQR;     yupperBound = yQ3 + 1.5*yIQR;

% Identify the outliers
x_outliers =find( (x < xlowerBound) | (x > xupperBound));
y_outliers =find( (y < ylowerBound) | (y > yupperBound));
outliers=unique([x_outliers,y_outliers]);
% Remove the outliers
x_filtered=x;
x_filtered(outliers)=[];
y_filtered=y;
y_filtered(outliers)=[];

end

% function [x_filtered, y_filtered]=filterXY(x, y)  %  this uses interpolation method
% x_mu = mean(x);      y_mu = mean(y);
% x_sigma = std(x);     y_sigma = std(y);
%
% x_threshold_high = x_mu + 3*x_sigma;    y_threshold_high = y_mu + 3*y_sigma;
% x_threshold_low = x_mu - 3*x_sigma;       y_threshold_low = y_mu - 3*y_sigma;
%
% % x_filtered = x(x > x_threshold_low & x < x_threshold_high);
% % y_filtered = y(y > y_threshold_low & y < y_threshold_high);
% x_outliers_idx = find(x< x_threshold_low & x > x_threshold_high);
% y_outliers_idx = find(y< y_threshold_low & y > y_threshold_high);
%
% % Remove the outliers from your data
% x_filtered =x;   y_filtered =y;
% x_filtered(unique([x_outliers_idx  y_outliers_idx])) = [];
% y_filtered(unique([x_outliers_idx  y_outliers_idx])) = [];
%
% % Now you can interpolate your data. For example, let's use linear interpolation:
% % xq = min(x):0.1:max(x);  % query points
% % yq = interp1(x, y, xq, 'linear');
% end