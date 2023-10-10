clear, close all, clc
addpath(genpath(pwd));
% ------------------------- Main Parameters of shapeing STL --------------------------
% --------- Written By Chao Song from Westlake University, China-------------------
Thickness=1; % Thickness of the hinges (mm)    !!!!!It CANNOT Be too large!!!!!

w1=0.5; % Half of the width of the hinges and plates (mm)
w2=0.5; % Width of the surfaces covered on hinges and plates (mm)

% These two are Kai's parameters, not using in this version
% widthScaleFactor = 0.8;
% plateThickness = 0.06; % (Fixed) thickness of plate.

% ---------------------------------------------------------------------------------------------
Case='shell323 Spindle objunitWith2D';
% Case = '2  2  3foldableSphereSymmetryOptimized';
% Case='10   4   1OutSphereAssem_05_31_NoOverLapping';
% Case='10   4   1OutSphereAssem_Cube_05_31_NoOverLapping';
load([Case,'.mat'])

Ntub.face=Ntub.face;
Ntub.node=Ntub.node*50; % Initial scale was too small, multiply 20 times
Ntub.node=Ntub.node+max(max(abs(Ntub.node)))+5; % '.STL' files are not using negative values as coordinate data


%%  ---------------------------------------Visualize-----------------------------------
viewPoint=[-23 30];
% figure(1)
% Plot(Ntub,viewPoint,'textx','ftexxt','Snapology',4,1,'on',[])
% axis on
% axis equal

figure(2)
hold on
plot3(Ntub.node(:,1),Ntub.node(:,2),Ntub.node(:,3),'*r')
axis on
axis equal

for i=1:size(Ntub.edge,1)
    plot3([Ntub.node(Ntub.edge(i,1),1),Ntub.node(Ntub.edge(i,2),1)],[Ntub.node(Ntub.edge(i,1),2),Ntub.node(Ntub.edge(i,2),2)],[Ntub.node(Ntub.edge(i,1),3),Ntub.node(Ntub.edge(i,2),3)],'k');
end

for i=1:size(Ntub.face,2)
    patch([Ntub.node(Ntub.face{1,i}(1),1),Ntub.node(Ntub.face{1,i}(2),1),Ntub.node(Ntub.face{1,i}(3),1)],...
        [Ntub.node(Ntub.face{1,i}(1),2),Ntub.node(Ntub.face{1,i}(2),2),Ntub.node(Ntub.face{1,i}(3),2)],...
        [Ntub.node(Ntub.face{1,i}(1),3),Ntub.node(Ntub.face{1,i}(2),3),Ntub.node(Ntub.face{1,i}(3),3)],'b');
end

%%  face smoothing for avoiding sharp angles in single face
UnitFacePairIdx=[1 2; 3 4; 5 6; 7 8; 9 10; 11  12; 13 14; 15 16; 17 18; 19 20;];
% 

% [Ntub,Ntub2D]=SmoothTriangularPanel(Ntub,Ntub2D,assemPoly);

figure(1)
Plot(Ntub,viewPoint,'textx','ftexxt','Snapology',4,1,'on',[])
axis on
axis equal


figure(3)
hold on
plot3(Ntub.node(:,1),Ntub.node(:,2),Ntub.node(:,3),'*r')
axis on
axis equal

for i=1:size(Ntub.edge,1)
    plot3([Ntub.node(Ntub.edge(i,1),1),Ntub.node(Ntub.edge(i,2),1)],[Ntub.node(Ntub.edge(i,1),2),Ntub.node(Ntub.edge(i,2),2)],[Ntub.node(Ntub.edge(i,1),3),Ntub.node(Ntub.edge(i,2),3)],'k');
end

for i=1:size(Ntub.face,2)
    patch([Ntub.node(Ntub.face{1,i}(1),1),Ntub.node(Ntub.face{1,i}(2),1),Ntub.node(Ntub.face{1,i}(3),1)],...
        [Ntub.node(Ntub.face{1,i}(1),2),Ntub.node(Ntub.face{1,i}(2),2),Ntub.node(Ntub.face{1,i}(3),2)],...
        [Ntub.node(Ntub.face{1,i}(1),3),Ntub.node(Ntub.face{1,i}(2),3),Ntub.node(Ntub.face{1,i}(3),3)],'b');
end


%% Find edges as hinges between faces.
hingeEx2Face = [];  % Hinges connecting two faces.
hingeEx3Face = [];  % Hinges connecting three faces.
for i = 1 : 1 : size(Ntub.face, 2)
    for j = i + 1 : 1 : size(Ntub.face, 2)
        h1 = intersect(Ntub.face{i}, Ntub.face{j}); % Find the nodes shared by each 2 faces.
        if length(h1) == 2
                    e11=find(Ntub.face{i}==h1(1));
                    e12=find(Ntub.face{i}==h1(2));
                    e21=find(Ntub.face{j}==h1(1));
                    e22=find(Ntub.face{j}==h1(2));
                    hingeEx2Face = [hingeEx2Face; h1, i, j, e11, e12, e21, e22];
                    % Each row of hingeEx takes the form of [n1, n2, i, j, e11, e12, e21, e22],
                    % where n1, n2 are the nodes of the hinge, and i, j are the indices of the two faces connected,
                    % e11, e12 indicate the position of two nodes in Ntub.face{i},
                    % e21, e22 indicate the position of two nodes in Ntub.face{j}.
        end
    end
end



Not2Plate=[];
Hinge3Plate=[];
for i=1:size(hingeEx2Face,1)  % Find these hinges shared by three plates, exclued them from hingeEx2Face
    Node_temp1=hingeEx2Face(i,1:2);
    for j=i+1:size(hingeEx2Face,1)
        if Node_temp1 == hingeEx2Face(j,1:2)
            for k=j+1:size(hingeEx2Face,1)
                if Node_temp1 == hingeEx2Face(k,1:2)
                    Hinge3Plate=[Hinge3Plate; i, j, k];
                    Not2Plate=[Not2Plate, i, j, k];
                end
            end
        end
    end
end

for i=1:size(Hinge3Plate)
    hingeEx3Face=[hingeEx3Face;hingeEx2Face(Hinge3Plate(i,1),1:2),hingeEx2Face(Hinge3Plate(i,1),3:4),hingeEx2Face(Hinge3Plate(i,2),4),hingeEx2Face(Hinge3Plate(i,1),5:8),hingeEx2Face(Hinge3Plate(i,2),7:8)];
    % Each row of hingeEx takes the form of [n1, n2, i, j, k, e11, e12, e21, e22, e31, e32],
    % where n1, n2 are the nodes of the hinge, and i, j, k are the indices of the three faces connected,
    % e11, e12 indicate the position of two nodes in Ntub.face{i},
    % e21, e22 indicate the position of two nodes in Ntub.face{j}.
    % e31, e32 indicate the position of two nodes in Ntub.face{k}.
end





hingeEx2Face(Not2Plate,:)=[];



%% Generate plates from given faces.


NormalVectorPlate=zeros(length(Ntub.face),3); % Find the normal vectors of plates.
for i=1:length(Ntub.face)
    NormalVectorPlate(i,:)=local_find_normal(Ntub.node(Ntub.face{1,i}(1),:),Ntub.node(Ntub.face{1,i}(2),:),Ntub.node(Ntub.face{1,i}(3),:),0);
end
NormalVectorPlateW1=NormalVectorPlate.*w1; % Make the vector length as the same as the 'w1'.

NormalVectorPlateW2=NormalVectorPlate.*(w2+w1);

MidPoint=zeros(length(Ntub.face),3); % Find the midpoint of plates.
for i=1:length(Ntub.face)
    MidPoint(i,:)=(Ntub.node(Ntub.face{1,i}(1),:)+Ntub.node(Ntub.face{1,i}(2),:)+Ntub.node(Ntub.face{1,i}(3),:))./3;
end

%Test if the NormalVectorPlate and MidPoint are correct
for i=1:length(Ntub.face)
    plot3([MidPoint(i,1),MidPoint(i,1)+NormalVectorPlateW1(i,1)], ...
        [MidPoint(i,2),MidPoint(i,2)+NormalVectorPlateW1(i,2)], ...
        [MidPoint(i,3),MidPoint(i,3)+NormalVectorPlateW1(i,3)] ...
        ,'g')
end

% % Test if hingeEx present correctly
% for i=1:size(hingeEx,1)
%     plot3([Ntub.node(Ntub.face{hingeEx(i,4)}(hingeEx(i,7)),1),Ntub.node(Ntub.face{hingeEx(i,4)}(hingeEx(i,8)),1)],...
%         [Ntub.node(Ntub.face{hingeEx(i,4)}(hingeEx(i,7)),2),Ntub.node(Ntub.face{hingeEx(i,4)}(hingeEx(i,8)),2)],...
%         [Ntub.node(Ntub.face{hingeEx(i,4)}(hingeEx(i,7)),3),Ntub.node(Ntub.face{hingeEx(i,4)}(hingeEx(i,8)),3)],'r')
% end


FaceOffset=zeros(3,3,length(Ntub.face)); % to Offset triangulars
for i=1:length(Ntub.face)
    [FaceOffset(1,:,i),FaceOffset(2,:,i),FaceOffset(3,:,i)] = OffsetTriangular(Ntub.node(Ntub.face{1,i}(1),:),Ntub.node(Ntub.face{1,i}(2),:),Ntub.node(Ntub.face{1,i}(3),:),Thickness);
    plot3([FaceOffset(:,1,i);FaceOffset(1,1,i)],...
        [FaceOffset(:,2,i);FaceOffset(1,2,i)],...
        [FaceOffset(:,3,i);FaceOffset(1,3,i)],'k')
end

figure(4)
hold on
axis on
axis equal

PlateNodes=zeros(6,3,length(Ntub.face)); % generate plates from triangulars
for i=1:length(Ntub.face)
    [PlateNodes(1,:,i),PlateNodes(2,:,i),PlateNodes(3,:,i),PlateNodes(4,:,i),PlateNodes(5,:,i),PlateNodes(6,:,i)] = Triangular2Plates(FaceOffset(1,:,i),FaceOffset(2,:,i),FaceOffset(3,:,i),NormalVectorPlateW1(i,:),-NormalVectorPlateW1(i,:));
    % PlateNodes is a 3D matrix (NodeNumber,NodeCoordinates,FaceNumber);
    % FaceNumber is the number of the face;
    % NodeNumber is the node number in this face, from 1 to 6 (each face has two extruded faces, with 6 nodes in total);
    % NodeCoordinates are the x, y, z coordinate of this node.
    
    plot3([PlateNodes(1:3,1,i);PlateNodes(1,1,i)],...
    [PlateNodes(1:3,2,i);PlateNodes(1,2,i)],...
    [PlateNodes(1:3,3,i);PlateNodes(1,3,i)],'k')

    plot3([PlateNodes(4:6,1,i);PlateNodes(4,1,i)],...
    [PlateNodes(4:6,2,i);PlateNodes(4,2,i)],...
    [PlateNodes(4:6,3,i);PlateNodes(4,3,i)],'k')

    plot3([PlateNodes(1,1,i);PlateNodes(4,1,i)],...
    [PlateNodes(1,2,i);PlateNodes(4,2,i)],...
    [PlateNodes(1,3,i);PlateNodes(4,3,i)],'k')

    plot3([PlateNodes(2,1,i);PlateNodes(5,1,i)],...
    [PlateNodes(2,2,i);PlateNodes(5,2,i)],...
    [PlateNodes(2,3,i);PlateNodes(5,3,i)],'k')

    plot3([PlateNodes(3,1,i);PlateNodes(6,1,i)],...
    [PlateNodes(3,2,i);PlateNodes(6,2,i)],...
    [PlateNodes(3,3,i);PlateNodes(6,3,i)],'k')

    patch([PlateNodes(1,1,i);PlateNodes(4,1,i);PlateNodes(5,1,i);PlateNodes(2,1,i)],...
    [PlateNodes(1,2,i);PlateNodes(4,2,i);PlateNodes(5,2,i);PlateNodes(2,2,i)],...
    [PlateNodes(1,3,i);PlateNodes(4,3,i);PlateNodes(5,3,i);PlateNodes(2,3,i)],'b')

    patch([PlateNodes(1,1,i);PlateNodes(4,1,i);PlateNodes(6,1,i);PlateNodes(3,1,i)],...
    [PlateNodes(1,2,i);PlateNodes(4,2,i);PlateNodes(6,2,i);PlateNodes(3,2,i)],...
    [PlateNodes(1,3,i);PlateNodes(4,3,i);PlateNodes(6,3,i);PlateNodes(3,3,i)],'b')

    patch([PlateNodes(2,1,i);PlateNodes(5,1,i);PlateNodes(6,1,i);PlateNodes(3,1,i)],...
    [PlateNodes(2,2,i);PlateNodes(5,2,i);PlateNodes(6,2,i);PlateNodes(3,2,i)],...
    [PlateNodes(2,3,i);PlateNodes(5,3,i);PlateNodes(6,3,i);PlateNodes(3,3,i)],'b')

    patch(PlateNodes(4:6,1,i),...
    PlateNodes(4:6,2,i),...
    PlateNodes(4:6,3,i),'b')

    patch(PlateNodes(1:3,1,i),...
    PlateNodes(1:3,2,i),...
    PlateNodes(1:3,3,i),'b')
end












%% Generate surface of plates.
SurfacePlatesNodes1=zeros(6,3,length(Ntub.face));
SurfacePlatesNodes2=zeros(6,3,length(Ntub.face));
% SurfacePlatesNodes1 and SurfacePlatesNodes2 are two surfaces on two directions
% They are both 3D matrix (NodeNumber,NodeCoordinates,FaceNumber);
% FaceNumber is the number of the face;
% NodeNumber is the node number in this face, from 1 to 6 (each face has two extruded faces, with 6 nodes in total);
% NodeCoordinates are the x, y, z coordinate of this node.

for i=1:length(Ntub.face)
    [SurfacePlatesNodes1(1,:,i),SurfacePlatesNodes1(2,:,i),SurfacePlatesNodes1(3,:,i),SurfacePlatesNodes1(4,:,i),SurfacePlatesNodes1(5,:,i),SurfacePlatesNodes1(6,:,i)] = Triangular2Plates(FaceOffset(1,:,i),FaceOffset(2,:,i),FaceOffset(3,:,i),NormalVectorPlateW2(i,:),NormalVectorPlateW1(i,:));
    [SurfacePlatesNodes2(1,:,i),SurfacePlatesNodes2(2,:,i),SurfacePlatesNodes2(3,:,i),SurfacePlatesNodes2(4,:,i),SurfacePlatesNodes2(5,:,i),SurfacePlatesNodes2(6,:,i)] = Triangular2Plates(FaceOffset(1,:,i),FaceOffset(2,:,i),FaceOffset(3,:,i),-NormalVectorPlateW1(i,:),-NormalVectorPlateW2(i,:));
    % SurfacePlatesNodes1(i,j,k)
    % k shows which face it is in Ntub.face
    % i is the number of node in that face, (1,:,:), (2,:,:), (3,:,:) are the OUTERS and (4,:,:), (5,:,:), (6,:,:) are the INNERS
    % j is the x, y, z coordinates of the node (i,:,k)

    % SurfacePlatesNodes2(i,j,k)
    % k shows which face it is in Ntub.face
    % i is the number of node in that face, (1,:,:), (2,:,:), (3,:,:) are the INNERS and (4,:,:), (5,:,:), (6,:,:) are the OUTERS
    % !!! (^THAT'S DIFFERENT) !!!
    % j is the x, y, z coordinates of the node (i,:,k)   


    
    patch([SurfacePlatesNodes1(1,1,i);SurfacePlatesNodes1(4,1,i);SurfacePlatesNodes1(5,1,i);SurfacePlatesNodes1(2,1,i)],...
    [SurfacePlatesNodes1(1,2,i);SurfacePlatesNodes1(4,2,i);SurfacePlatesNodes1(5,2,i);SurfacePlatesNodes1(2,2,i)],...
    [SurfacePlatesNodes1(1,3,i);SurfacePlatesNodes1(4,3,i);SurfacePlatesNodes1(5,3,i);SurfacePlatesNodes1(2,3,i)],'g')

    patch([SurfacePlatesNodes1(1,1,i);SurfacePlatesNodes1(4,1,i);SurfacePlatesNodes1(6,1,i);SurfacePlatesNodes1(3,1,i)],...
    [SurfacePlatesNodes1(1,2,i);SurfacePlatesNodes1(4,2,i);SurfacePlatesNodes1(6,2,i);SurfacePlatesNodes1(3,2,i)],...
    [SurfacePlatesNodes1(1,3,i);SurfacePlatesNodes1(4,3,i);SurfacePlatesNodes1(6,3,i);SurfacePlatesNodes1(3,3,i)],'g')

    patch([SurfacePlatesNodes1(2,1,i);SurfacePlatesNodes1(5,1,i);SurfacePlatesNodes1(6,1,i);SurfacePlatesNodes1(3,1,i)],...
    [SurfacePlatesNodes1(2,2,i);SurfacePlatesNodes1(5,2,i);SurfacePlatesNodes1(6,2,i);SurfacePlatesNodes1(3,2,i)],...
    [SurfacePlatesNodes1(2,3,i);SurfacePlatesNodes1(5,3,i);SurfacePlatesNodes1(6,3,i);SurfacePlatesNodes1(3,3,i)],'g')

    patch(SurfacePlatesNodes1(4:6,1,i),...
    SurfacePlatesNodes1(4:6,2,i),...
    SurfacePlatesNodes1(4:6,3,i),'g')

    patch(SurfacePlatesNodes1(1:3,1,i),...
    SurfacePlatesNodes1(1:3,2,i),...
    SurfacePlatesNodes1(1:3,3,i),'g')


    patch([SurfacePlatesNodes2(1,1,i);SurfacePlatesNodes2(4,1,i);SurfacePlatesNodes2(5,1,i);SurfacePlatesNodes2(2,1,i)],...
    [SurfacePlatesNodes2(1,2,i);SurfacePlatesNodes2(4,2,i);SurfacePlatesNodes2(5,2,i);SurfacePlatesNodes2(2,2,i)],...
    [SurfacePlatesNodes2(1,3,i);SurfacePlatesNodes2(4,3,i);SurfacePlatesNodes2(5,3,i);SurfacePlatesNodes2(2,3,i)],'g')

    patch([SurfacePlatesNodes2(1,1,i);SurfacePlatesNodes2(4,1,i);SurfacePlatesNodes2(6,1,i);SurfacePlatesNodes2(3,1,i)],...
    [SurfacePlatesNodes2(1,2,i);SurfacePlatesNodes2(4,2,i);SurfacePlatesNodes2(6,2,i);SurfacePlatesNodes2(3,2,i)],...
    [SurfacePlatesNodes2(1,3,i);SurfacePlatesNodes2(4,3,i);SurfacePlatesNodes2(6,3,i);SurfacePlatesNodes2(3,3,i)],'g')

    patch([SurfacePlatesNodes2(2,1,i);SurfacePlatesNodes2(5,1,i);SurfacePlatesNodes2(6,1,i);SurfacePlatesNodes2(3,1,i)],...
    [SurfacePlatesNodes2(2,2,i);SurfacePlatesNodes2(5,2,i);SurfacePlatesNodes2(6,2,i);SurfacePlatesNodes2(3,2,i)],...
    [SurfacePlatesNodes2(2,3,i);SurfacePlatesNodes2(5,3,i);SurfacePlatesNodes2(6,3,i);SurfacePlatesNodes2(3,3,i)],'g')

    patch(SurfacePlatesNodes2(4:6,1,i),...
    SurfacePlatesNodes2(4:6,2,i),...
    SurfacePlatesNodes2(4:6,3,i),'g')

    patch(SurfacePlatesNodes2(1:3,1,i),...
    SurfacePlatesNodes2(1:3,2,i),...
    SurfacePlatesNodes2(1:3,3,i),'g')



end




%% Generate hinges (connecting two plates) from the plates

HingeNodesTwo=zeros(8,3,size(hingeEx2Face,1));
% HingeNodesTwo is a 3D matrix (NodeNumber,NodeCoordinates,FaceNumber);
% FaceNumber is the number of the hinges connected with three hinges;
% NodeNumber is the node number in this hinge, from 1 to 8 (each hinge has 8 nodes in total, 4 nodes for each plate connected);
% NodeCoordinates are the x, y, z coordinate of this node.

for i=1:size(hingeEx2Face,1) % Find these hinges' nodes in plates
    HingeNodesTwo(1,:,i)=PlateNodes(hingeEx2Face(i,5),:,hingeEx2Face(i,3));
    HingeNodesTwo(2,:,i)=PlateNodes(hingeEx2Face(i,6),:,hingeEx2Face(i,3));
    HingeNodesTwo(3,:,i)=PlateNodes(hingeEx2Face(i,7),:,hingeEx2Face(i,4));
    HingeNodesTwo(4,:,i)=PlateNodes(hingeEx2Face(i,8),:,hingeEx2Face(i,4));
    HingeNodesTwo(5,:,i)=PlateNodes(hingeEx2Face(i,5)+3,:,hingeEx2Face(i,3));
    HingeNodesTwo(6,:,i)=PlateNodes(hingeEx2Face(i,6)+3,:,hingeEx2Face(i,3));
    HingeNodesTwo(7,:,i)=PlateNodes(hingeEx2Face(i,7)+3,:,hingeEx2Face(i,4));
    HingeNodesTwo(8,:,i)=PlateNodes(hingeEx2Face(i,8)+3,:,hingeEx2Face(i,4));

    % One Plate: (1,:,i), (5,:,i)     --     (2,:,i), (6,:,i) 
    % Another Plate: (3,:,i), (7,:,i)     --     (4,:,i), (8,:,i) 

    for j=1:3
    plot3([HingeNodesTwo(j,1,i),HingeNodesTwo(j+1,1,i)], ...
        [HingeNodesTwo(j,2,i),HingeNodesTwo(j+1,2,i)], ...
        [HingeNodesTwo(j,3,i),HingeNodesTwo(j+1,3,i)],'r')
    plot3([HingeNodesTwo(j+4,1,i),HingeNodesTwo(j+5,1,i)], ...
        [HingeNodesTwo(j+4,2,i),HingeNodesTwo(j+5,2,i)], ...
        [HingeNodesTwo(j+4,3,i),HingeNodesTwo(j+5,3,i)],'r')
    end
    plot3([HingeNodesTwo(4,1,i),HingeNodesTwo(1,1,i)], ...
        [HingeNodesTwo(4,2,i),HingeNodesTwo(1,2,i)], ...
        [HingeNodesTwo(4,3,i),HingeNodesTwo(1,3,i)],'r')
    plot3([HingeNodesTwo(8,1,i),HingeNodesTwo(5,1,i)], ...
        [HingeNodesTwo(8,2,i),HingeNodesTwo(5,2,i)], ...
        [HingeNodesTwo(8,3,i),HingeNodesTwo(5,3,i)],'r')
end


%% Generate surface of hinges (connecting two plates)
SurfaceHingeNodesTwo1=zeros(8,3,size(hingeEx2Face,1));
for i=1:size(hingeEx2Face,1) % Find these hinges' nodes in plates
    SurfaceHingeNodesTwo1(1,:,i)=SurfacePlatesNodes1(hingeEx2Face(i,5),:,hingeEx2Face(i,3));
    SurfaceHingeNodesTwo1(2,:,i)=SurfacePlatesNodes1(hingeEx2Face(i,6),:,hingeEx2Face(i,3));
    SurfaceHingeNodesTwo1(3,:,i)=SurfacePlatesNodes1(hingeEx2Face(i,7),:,hingeEx2Face(i,4));
    SurfaceHingeNodesTwo1(4,:,i)=SurfacePlatesNodes1(hingeEx2Face(i,8),:,hingeEx2Face(i,4));
    SurfaceHingeNodesTwo1(5,:,i)=SurfacePlatesNodes1(hingeEx2Face(i,5)+3,:,hingeEx2Face(i,3));
    SurfaceHingeNodesTwo1(6,:,i)=SurfacePlatesNodes1(hingeEx2Face(i,6)+3,:,hingeEx2Face(i,3));
    SurfaceHingeNodesTwo1(7,:,i)=SurfacePlatesNodes1(hingeEx2Face(i,7)+3,:,hingeEx2Face(i,4));
    SurfaceHingeNodesTwo1(8,:,i)=SurfacePlatesNodes1(hingeEx2Face(i,8)+3,:,hingeEx2Face(i,4));

    for j=1:3
    plot3([SurfaceHingeNodesTwo1(j,1,i),SurfaceHingeNodesTwo1(j+1,1,i)], ...
        [SurfaceHingeNodesTwo1(j,2,i),SurfaceHingeNodesTwo1(j+1,2,i)], ...
        [SurfaceHingeNodesTwo1(j,3,i),SurfaceHingeNodesTwo1(j+1,3,i)],'r')
    plot3([SurfaceHingeNodesTwo1(j+4,1,i),SurfaceHingeNodesTwo1(j+5,1,i)], ...
        [SurfaceHingeNodesTwo1(j+4,2,i),SurfaceHingeNodesTwo1(j+5,2,i)], ...
        [SurfaceHingeNodesTwo1(j+4,3,i),SurfaceHingeNodesTwo1(j+5,3,i)],'r')
    end
    plot3([SurfaceHingeNodesTwo1(4,1,i),SurfaceHingeNodesTwo1(1,1,i)], ...
        [SurfaceHingeNodesTwo1(4,2,i),SurfaceHingeNodesTwo1(1,2,i)], ...
        [SurfaceHingeNodesTwo1(4,3,i),SurfaceHingeNodesTwo1(1,3,i)],'r')
    plot3([SurfaceHingeNodesTwo1(8,1,i),SurfaceHingeNodesTwo1(5,1,i)], ...
        [SurfaceHingeNodesTwo1(8,2,i),SurfaceHingeNodesTwo1(5,2,i)], ...
        [SurfaceHingeNodesTwo1(8,3,i),SurfaceHingeNodesTwo1(5,3,i)],'r')

end

SurfaceHingeNodesTwo2=zeros(8,3,size(hingeEx2Face,1));
for i=1:size(hingeEx2Face,1) % Find these hinges' nodes in plates
    SurfaceHingeNodesTwo2(1,:,i)=SurfacePlatesNodes2(hingeEx2Face(i,5),:,hingeEx2Face(i,3));
    SurfaceHingeNodesTwo2(2,:,i)=SurfacePlatesNodes2(hingeEx2Face(i,6),:,hingeEx2Face(i,3));
    SurfaceHingeNodesTwo2(3,:,i)=SurfacePlatesNodes2(hingeEx2Face(i,7),:,hingeEx2Face(i,4));
    SurfaceHingeNodesTwo2(4,:,i)=SurfacePlatesNodes2(hingeEx2Face(i,8),:,hingeEx2Face(i,4));
    SurfaceHingeNodesTwo2(5,:,i)=SurfacePlatesNodes2(hingeEx2Face(i,5)+3,:,hingeEx2Face(i,3));
    SurfaceHingeNodesTwo2(6,:,i)=SurfacePlatesNodes2(hingeEx2Face(i,6)+3,:,hingeEx2Face(i,3));
    SurfaceHingeNodesTwo2(7,:,i)=SurfacePlatesNodes2(hingeEx2Face(i,7)+3,:,hingeEx2Face(i,4));
    SurfaceHingeNodesTwo2(8,:,i)=SurfacePlatesNodes2(hingeEx2Face(i,8)+3,:,hingeEx2Face(i,4));

    for j=1:3
    plot3([SurfaceHingeNodesTwo2(j,1,i),SurfaceHingeNodesTwo2(j+1,1,i)], ...
        [SurfaceHingeNodesTwo2(j,2,i),SurfaceHingeNodesTwo2(j+1,2,i)], ...
        [SurfaceHingeNodesTwo2(j,3,i),SurfaceHingeNodesTwo2(j+1,3,i)],'r')
    plot3([SurfaceHingeNodesTwo2(j+4,1,i),SurfaceHingeNodesTwo2(j+5,1,i)], ...
        [SurfaceHingeNodesTwo2(j+4,2,i),SurfaceHingeNodesTwo2(j+5,2,i)], ...
        [SurfaceHingeNodesTwo2(j+4,3,i),SurfaceHingeNodesTwo2(j+5,3,i)],'r')
    end
    plot3([SurfaceHingeNodesTwo2(4,1,i),SurfaceHingeNodesTwo2(1,1,i)], ...
        [SurfaceHingeNodesTwo2(4,2,i),SurfaceHingeNodesTwo2(1,2,i)], ...
        [SurfaceHingeNodesTwo2(4,3,i),SurfaceHingeNodesTwo2(1,3,i)],'r')
    plot3([SurfaceHingeNodesTwo2(8,1,i),SurfaceHingeNodesTwo2(5,1,i)], ...
        [SurfaceHingeNodesTwo2(8,2,i),SurfaceHingeNodesTwo2(5,2,i)], ...
        [SurfaceHingeNodesTwo2(8,3,i),SurfaceHingeNodesTwo2(5,3,i)],'r')
end






%% Generate hinges (connecting three plates) from the plates


HingeNodesThree=zeros(12,3,size(hingeEx3Face,1));
% HingeNodesTwo is a 3D matrix (NodeNumber,NodeCoordinates,FaceNumber);
% FaceNumber is the number of the hinge connected with three hinges;
% NodeNumber is the node number in this hinge, from 1 to 12 (each hinge has 12 nodes in total, 4 nodes for each plate connected);
% NodeCoordinates are the x, y, z coordinate of this node.

for i=1:size(hingeEx3Face,1) % Find these hinges' nodes in plates
    HingeNodesThree(1,:,i)=PlateNodes(hingeEx3Face(i,6),:,hingeEx3Face(i,3)); 
    HingeNodesThree(2,:,i)=PlateNodes(hingeEx3Face(i,7),:,hingeEx3Face(i,3));
    HingeNodesThree(3,:,i)=PlateNodes(hingeEx3Face(i,8),:,hingeEx3Face(i,4)); 
    HingeNodesThree(4,:,i)=PlateNodes(hingeEx3Face(i,9),:,hingeEx3Face(i,4));
    HingeNodesThree(5,:,i)=PlateNodes(hingeEx3Face(i,10),:,hingeEx3Face(i,5)); 
    HingeNodesThree(6,:,i)=PlateNodes(hingeEx3Face(i,11),:,hingeEx3Face(i,5));

    HingeNodesThree(7,:,i)=PlateNodes(hingeEx3Face(i,6)+3,:,hingeEx3Face(i,3));
    HingeNodesThree(8,:,i)=PlateNodes(hingeEx3Face(i,7)+3,:,hingeEx3Face(i,3));
    HingeNodesThree(9,:,i)=PlateNodes(hingeEx3Face(i,8)+3,:,hingeEx3Face(i,4));
    HingeNodesThree(10,:,i)=PlateNodes(hingeEx3Face(i,9)+3,:,hingeEx3Face(i,4));
    HingeNodesThree(11,:,i)=PlateNodes(hingeEx3Face(i,10)+3,:,hingeEx3Face(i,5));
    HingeNodesThree(12,:,i)=PlateNodes(hingeEx3Face(i,11)+3,:,hingeEx3Face(i,5));

    % plate #1 contains (1,:,i), (7,:,i) on one side and (2,:,i), (8,:,i) on another side
    % plate #2 contains (3,:,i), (9,:,i) on one side and (4,:,i), (10,:,i) on another side
    % plate #3 contains (5,:,i), (11,:,i) on one side and (6,:,i), (12,:,i) on another side



%     HingeNodesThree(15,:,i)=


    for j=1:12
    plot3(HingeNodesThree(j,1,i), ...
        HingeNodesThree(j,2,i), ...
        HingeNodesThree(j,3,i),'ko')
    end
end

figure
hold on
axis equal
box on




for i=1:length(Ntub.face)
    [SurfacePlatesNodes1(1,:,i),SurfacePlatesNodes1(2,:,i),SurfacePlatesNodes1(3,:,i),SurfacePlatesNodes1(4,:,i),SurfacePlatesNodes1(5,:,i),SurfacePlatesNodes1(6,:,i)] = Triangular2Plates(FaceOffset(1,:,i),FaceOffset(2,:,i),FaceOffset(3,:,i),NormalVectorPlateW2(i,:),NormalVectorPlateW1(i,:));
    [SurfacePlatesNodes2(1,:,i),SurfacePlatesNodes2(2,:,i),SurfacePlatesNodes2(3,:,i),SurfacePlatesNodes2(4,:,i),SurfacePlatesNodes2(5,:,i),SurfacePlatesNodes2(6,:,i)] = Triangular2Plates(FaceOffset(1,:,i),FaceOffset(2,:,i),FaceOffset(3,:,i),-NormalVectorPlateW1(i,:),-NormalVectorPlateW2(i,:));
    
    
    patch([SurfacePlatesNodes1(1,1,i);SurfacePlatesNodes1(4,1,i);SurfacePlatesNodes1(5,1,i);SurfacePlatesNodes1(2,1,i)],...
    [SurfacePlatesNodes1(1,2,i);SurfacePlatesNodes1(4,2,i);SurfacePlatesNodes1(5,2,i);SurfacePlatesNodes1(2,2,i)],...
    [SurfacePlatesNodes1(1,3,i);SurfacePlatesNodes1(4,3,i);SurfacePlatesNodes1(5,3,i);SurfacePlatesNodes1(2,3,i)],'g','FaceAlpha',.1)

    patch([SurfacePlatesNodes1(1,1,i);SurfacePlatesNodes1(4,1,i);SurfacePlatesNodes1(6,1,i);SurfacePlatesNodes1(3,1,i)],...
    [SurfacePlatesNodes1(1,2,i);SurfacePlatesNodes1(4,2,i);SurfacePlatesNodes1(6,2,i);SurfacePlatesNodes1(3,2,i)],...
    [SurfacePlatesNodes1(1,3,i);SurfacePlatesNodes1(4,3,i);SurfacePlatesNodes1(6,3,i);SurfacePlatesNodes1(3,3,i)],'g','FaceAlpha',.1)

    patch([SurfacePlatesNodes1(2,1,i);SurfacePlatesNodes1(5,1,i);SurfacePlatesNodes1(6,1,i);SurfacePlatesNodes1(3,1,i)],...
    [SurfacePlatesNodes1(2,2,i);SurfacePlatesNodes1(5,2,i);SurfacePlatesNodes1(6,2,i);SurfacePlatesNodes1(3,2,i)],...
    [SurfacePlatesNodes1(2,3,i);SurfacePlatesNodes1(5,3,i);SurfacePlatesNodes1(6,3,i);SurfacePlatesNodes1(3,3,i)],'g','FaceAlpha',.1)

    patch(SurfacePlatesNodes1(4:6,1,i),...
    SurfacePlatesNodes1(4:6,2,i),...
    SurfacePlatesNodes1(4:6,3,i),'g','FaceAlpha',.1)

    patch(SurfacePlatesNodes1(1:3,1,i),...
    SurfacePlatesNodes1(1:3,2,i),...
    SurfacePlatesNodes1(1:3,3,i),'g','FaceAlpha',.1)


    patch([SurfacePlatesNodes2(1,1,i);SurfacePlatesNodes2(4,1,i);SurfacePlatesNodes2(5,1,i);SurfacePlatesNodes2(2,1,i)],...
    [SurfacePlatesNodes2(1,2,i);SurfacePlatesNodes2(4,2,i);SurfacePlatesNodes2(5,2,i);SurfacePlatesNodes2(2,2,i)],...
    [SurfacePlatesNodes2(1,3,i);SurfacePlatesNodes2(4,3,i);SurfacePlatesNodes2(5,3,i);SurfacePlatesNodes2(2,3,i)],'g','FaceAlpha',.1)

    patch([SurfacePlatesNodes2(1,1,i);SurfacePlatesNodes2(4,1,i);SurfacePlatesNodes2(6,1,i);SurfacePlatesNodes2(3,1,i)],...
    [SurfacePlatesNodes2(1,2,i);SurfacePlatesNodes2(4,2,i);SurfacePlatesNodes2(6,2,i);SurfacePlatesNodes2(3,2,i)],...
    [SurfacePlatesNodes2(1,3,i);SurfacePlatesNodes2(4,3,i);SurfacePlatesNodes2(6,3,i);SurfacePlatesNodes2(3,3,i)],'g','FaceAlpha',.1)

    patch([SurfacePlatesNodes2(2,1,i);SurfacePlatesNodes2(5,1,i);SurfacePlatesNodes2(6,1,i);SurfacePlatesNodes2(3,1,i)],...
    [SurfacePlatesNodes2(2,2,i);SurfacePlatesNodes2(5,2,i);SurfacePlatesNodes2(6,2,i);SurfacePlatesNodes2(3,2,i)],...
    [SurfacePlatesNodes2(2,3,i);SurfacePlatesNodes2(5,3,i);SurfacePlatesNodes2(6,3,i);SurfacePlatesNodes2(3,3,i)],'g','FaceAlpha',.1)

    patch(SurfacePlatesNodes2(4:6,1,i),...
    SurfacePlatesNodes2(4:6,2,i),...
    SurfacePlatesNodes2(4:6,3,i),'g','FaceAlpha',.1)

    patch(SurfacePlatesNodes2(1:3,1,i),...
    SurfacePlatesNodes2(1:3,2,i),...
    SurfacePlatesNodes2(1:3,3,i),'g','FaceAlpha',.1)
end



for i=1:size(hingeEx3Face,1)

    for j=[1 3 5 7 9 11]
    plot3(HingeNodesThree(j,1,i), ...
        HingeNodesThree(j,2,i), ...
        HingeNodesThree(j,3,i),'r*')
    end

    for j=[2 4 6 8 10 12]
    plot3(HingeNodesThree(j,1,i), ...
        HingeNodesThree(j,2,i), ...
        HingeNodesThree(j,3,i),'b*')
    end

end



%  对三个plate共享的hinge上的点进行排序,
%  第一步是找到顺序矩阵（HingeNodesThreeOrder），
%  第二步是用顺序矩阵对（HingeNodesThree）进行排序，得到（HingeNodesThreeSorted）。

HingeNodesThreeHalf=HingeNodesThree([1 7 3 9 5 11],:,:);  % 1和7（第一个plate），3和9（第二个plate），5和11（第三个plate）

HingeNodesThreeOrder=zeros(size(HingeNodesThreeHalf,3),6);  %  第一步是找到顺序矩阵（HingeNodesThreeOrder）。
for i=1:size(HingeNodesThreeOrder,1)
% for i=2
    HingeNodesThreeOrder(i,1:2)=[1,2];
    temp1=min(norm(HingeNodesThreeHalf(1,:,i)-HingeNodesThreeHalf(3,:,i)),norm(HingeNodesThreeHalf(1,:,i)-HingeNodesThreeHalf(4,:,i)));
    temp2=min(norm(HingeNodesThreeHalf(2,:,i)-HingeNodesThreeHalf(3,:,i)),norm(HingeNodesThreeHalf(2,:,i)-HingeNodesThreeHalf(4,:,i)));
    if temp1 > temp2
        temp3=norm(HingeNodesThreeHalf(2,:,i)-HingeNodesThreeHalf(3,:,i));
        temp4=norm(HingeNodesThreeHalf(2,:,i)-HingeNodesThreeHalf(4,:,i));
        if temp3 > temp4
            HingeNodesThreeOrder(i,3:4)=[4,3];
        else
            HingeNodesThreeOrder(i,3:4)=[3,4];
        end
        temp5=norm(HingeNodesThreeHalf(1,:,i)-HingeNodesThreeHalf(5,:,i));
        temp6=norm(HingeNodesThreeHalf(1,:,i)-HingeNodesThreeHalf(6,:,i));
        if temp5 > temp6
            HingeNodesThreeOrder(i,5:6)=[5,6];
        else
            HingeNodesThreeOrder(i,5:6)=[6,5];
        end
    else
        temp3=norm(HingeNodesThreeHalf(1,:,i)-HingeNodesThreeHalf(3,:,i));
        temp4=norm(HingeNodesThreeHalf(1,:,i)-HingeNodesThreeHalf(4,:,i));
        if temp3 > temp4
            HingeNodesThreeOrder(i,5:6)=[3,4];
        else
            HingeNodesThreeOrder(i,5:6)=[4,3];
        end
        temp5=norm(HingeNodesThreeHalf(2,:,i)-HingeNodesThreeHalf(5,:,i));
        temp6=norm(HingeNodesThreeHalf(2,:,i)-HingeNodesThreeHalf(6,:,i));
        if temp5 > temp6
            HingeNodesThreeOrder(i,3:4)=[6,5];
        else
            HingeNodesThreeOrder(i,3:4)=[5,6];
        end
    end
end


HingeNodesThreeSorted = HingeNodesThree;

for i=1:size(HingeNodesThree,3) %  第二步是用顺序矩阵对（HingeNodesThree）进行排序，得到（HingeNodesThreeSorted）。
    temp1=HingeNodesThree([1 7 3 9 5 11],:,i);
    temp2=HingeNodesThree([2 8 4 10 6 12],:,i);
    HingeNodesThreeSorted([1 3 5 7 9 11],:,i)=temp1(HingeNodesThreeOrder(i,:),:);
    HingeNodesThreeSorted([2 4 6 8 10 12],:,i)=temp2(HingeNodesThreeOrder(i,:),:);
end


for i=1:size(HingeNodesThreeSorted,3) % 在HingeNodesThreeSorted添加两个面的中点
    HingeNodesThreeSorted(13,:,i)=(HingeNodesThreeSorted(1,:,i)+ ...   % (13,:,i) is the middle point of node (1,:,i), (3,:,i), (5,:,i), (7,:,i), (9,:,i), (11,:,i)
        HingeNodesThreeSorted(3,:,i)+ ...
        HingeNodesThreeSorted(5,:,i)+ ...
        HingeNodesThreeSorted(7,:,i)+ ...
        HingeNodesThreeSorted(9,:,i)+ ...
        HingeNodesThreeSorted(11,:,i))./6;

    HingeNodesThreeSorted(14,:,i)=(HingeNodesThreeSorted(2,:,i)+ ...   % (14,:,i) is the middle point of node (2,:,i), (4,:,i), (6,:,i), (8,:,i), (10,:,i), (12,:,i)
        HingeNodesThreeSorted(4,:,i)+ ...
        HingeNodesThreeSorted(6,:,i)+ ...
        HingeNodesThreeSorted(8,:,i)+ ...
        HingeNodesThreeSorted(10,:,i)+ ...
        HingeNodesThreeSorted(12,:,i))./6;

    % 通过HingeNodesThreeSorted中点，每个面插入三个插值点。
    HingeNodesThreeSorted(15,:,i)=(HingeNodesThreeSorted(13,:,i)+(HingeNodesThreeSorted(7,:,i)+HingeNodesThreeSorted(9,:,i))./2)./2;
    HingeNodesThreeSorted(16,:,i)=(HingeNodesThreeSorted(14,:,i)+(HingeNodesThreeSorted(8,:,i)+HingeNodesThreeSorted(10,:,i))./2)./2;

    HingeNodesThreeSorted(17,:,i)=(HingeNodesThreeSorted(13,:,i)+(HingeNodesThreeSorted(1,:,i)+HingeNodesThreeSorted(11,:,i))./2)./2;
    HingeNodesThreeSorted(18,:,i)=(HingeNodesThreeSorted(14,:,i)+(HingeNodesThreeSorted(2,:,i)+HingeNodesThreeSorted(12,:,i))./2)./2;

    HingeNodesThreeSorted(19,:,i)=(HingeNodesThreeSorted(13,:,i)+(HingeNodesThreeSorted(3,:,i)+HingeNodesThreeSorted(5,:,i))./2)./2;
    HingeNodesThreeSorted(20,:,i)=(HingeNodesThreeSorted(14,:,i)+(HingeNodesThreeSorted(4,:,i)+HingeNodesThreeSorted(6,:,i))./2)./2;
end
% Using above code, we turn the hinge in order 1-3-19-5-7-15-9-11-17-1 (逆时针), 中点是 13
% 另一侧， 2-4-20-6-8-16-10-12-18-2 (逆时针), 中点是 14



%% Generate the surface of hinges (connecting three plates) from the plates 


SurfaceHingeNodesThree=zeros(12,3,size(hingeEx3Face,1));

for i=1:size(hingeEx3Face,1) % Find these hinges' nodes in plates
    SurfaceHingeNodesThree(1,:,i)=SurfacePlatesNodes1(hingeEx3Face(i,6),:,hingeEx3Face(i,3)); 
    SurfaceHingeNodesThree(2,:,i)=SurfacePlatesNodes1(hingeEx3Face(i,7),:,hingeEx3Face(i,3));
    SurfaceHingeNodesThree(3,:,i)=SurfacePlatesNodes1(hingeEx3Face(i,8),:,hingeEx3Face(i,4)); 
    SurfaceHingeNodesThree(4,:,i)=SurfacePlatesNodes1(hingeEx3Face(i,9),:,hingeEx3Face(i,4));
    SurfaceHingeNodesThree(5,:,i)=SurfacePlatesNodes1(hingeEx3Face(i,10),:,hingeEx3Face(i,5)); 
    SurfaceHingeNodesThree(6,:,i)=SurfacePlatesNodes1(hingeEx3Face(i,11),:,hingeEx3Face(i,5));

    SurfaceHingeNodesThree(7,:,i)=SurfacePlatesNodes2(hingeEx3Face(i,6)+3,:,hingeEx3Face(i,3));
    SurfaceHingeNodesThree(8,:,i)=SurfacePlatesNodes2(hingeEx3Face(i,7)+3,:,hingeEx3Face(i,3));
    SurfaceHingeNodesThree(9,:,i)=SurfacePlatesNodes2(hingeEx3Face(i,8)+3,:,hingeEx3Face(i,4));
    SurfaceHingeNodesThree(10,:,i)=SurfacePlatesNodes2(hingeEx3Face(i,9)+3,:,hingeEx3Face(i,4));
    SurfaceHingeNodesThree(11,:,i)=SurfacePlatesNodes2(hingeEx3Face(i,10)+3,:,hingeEx3Face(i,5));
    SurfaceHingeNodesThree(12,:,i)=SurfacePlatesNodes2(hingeEx3Face(i,11)+3,:,hingeEx3Face(i,5));
end




SurfaceHingeNodesThreeSorted = SurfaceHingeNodesThree;

for i=1:size(SurfaceHingeNodesThree,3) %  用顺序矩阵对（SurfaceHingeNodesThree）进行排序，得到（SurfaceHingeNodesThreeSorted）。
    temp1=SurfaceHingeNodesThree([1 7 3 9 5 11],:,i);
    temp2=SurfaceHingeNodesThree([2 8 4 10 6 12],:,i);
    SurfaceHingeNodesThreeSorted([1 3 5 7 9 11],:,i)=temp1(HingeNodesThreeOrder(i,:),:);
    SurfaceHingeNodesThreeSorted([2 4 6 8 10 12],:,i)=temp2(HingeNodesThreeOrder(i,:),:);
end


for i=1:size(SurfaceHingeNodesThreeSorted,3) % 在SurfaceHingeNodesThreeSorted添加两个面的中点
    SurfaceHingeNodesThreeSorted(13,:,i)=(SurfaceHingeNodesThreeSorted(1,:,i)+ ...   % (13,:,i) is the middle point of node (1,:,i), (3,:,i), (5,:,i), (7,:,i), (9,:,i), (11,:,i)
        SurfaceHingeNodesThreeSorted(3,:,i)+ ...
        SurfaceHingeNodesThreeSorted(5,:,i)+ ...
        SurfaceHingeNodesThreeSorted(7,:,i)+ ...
        SurfaceHingeNodesThreeSorted(9,:,i)+ ...
        SurfaceHingeNodesThreeSorted(11,:,i))./6;

    SurfaceHingeNodesThreeSorted(14,:,i)=(SurfaceHingeNodesThreeSorted(2,:,i)+ ...   % (14,:,i) is the middle point of node (2,:,i), (4,:,i), (6,:,i), (8,:,i), (10,:,i), (12,:,i)
        SurfaceHingeNodesThreeSorted(4,:,i)+ ...
        SurfaceHingeNodesThreeSorted(6,:,i)+ ...
        SurfaceHingeNodesThreeSorted(8,:,i)+ ...
        SurfaceHingeNodesThreeSorted(10,:,i)+ ...
        SurfaceHingeNodesThreeSorted(12,:,i))./6;

    % 通过SurfaceHingeNodesThreeSorted中点，每个面插入三个插值点。
    SurfaceHingeNodesThreeSorted(15,:,i)=(SurfaceHingeNodesThreeSorted(13,:,i)+(SurfaceHingeNodesThreeSorted(7,:,i)+SurfaceHingeNodesThreeSorted(9,:,i))./2)./2;
    SurfaceHingeNodesThreeSorted(16,:,i)=(SurfaceHingeNodesThreeSorted(14,:,i)+(SurfaceHingeNodesThreeSorted(8,:,i)+SurfaceHingeNodesThreeSorted(10,:,i))./2)./2;

    SurfaceHingeNodesThreeSorted(17,:,i)=(SurfaceHingeNodesThreeSorted(13,:,i)+(SurfaceHingeNodesThreeSorted(1,:,i)+SurfaceHingeNodesThreeSorted(11,:,i))./2)./2;
    SurfaceHingeNodesThreeSorted(18,:,i)=(SurfaceHingeNodesThreeSorted(14,:,i)+(SurfaceHingeNodesThreeSorted(2,:,i)+SurfaceHingeNodesThreeSorted(12,:,i))./2)./2;

    SurfaceHingeNodesThreeSorted(19,:,i)=(SurfaceHingeNodesThreeSorted(13,:,i)+(SurfaceHingeNodesThreeSorted(3,:,i)+SurfaceHingeNodesThreeSorted(5,:,i))./2)./2;
    SurfaceHingeNodesThreeSorted(20,:,i)=(SurfaceHingeNodesThreeSorted(14,:,i)+(SurfaceHingeNodesThreeSorted(4,:,i)+SurfaceHingeNodesThreeSorted(6,:,i))./2)./2;

    SurfaceHingeNodesThreeSorted(21:40,:,i)=HingeNodesThreeSorted(1:20,:,i);

end

% Using above code, we turn the hinge in order 1-3-19-5-7-15-9-11-17-1 (逆时针，外圈), 中点是 13；21-23-39-25-27-35-29-31-37-21 (逆时针，内圈)
% 另一侧， 2-4-20-6-8-16-10-12-18-2 (逆时针，外圈), 中点是 14；22-24-40-26-28-36-30-32-38-22 (逆时针，内圈)

















%% Plots
















%% Write STL file


% Step #1, Write the plate part, using [PlateNodes]
fid=fopen('Plate.stl','w'); % STL file of the plate
fprintf(fid,'solid\n');

STLOrderPlateMatrix=...    % Each solid has 8 triangulars 
    [1	3	2;
    4	5	6;
    1	2	4;
    2	5	4;
    2	3	6;
    2	6	5;
    1	6	3;
    1	4	6];
for j=1:size(PlateNodes,3)
    Node_temp=PlateNodes(:,:,j);
    for i=1:size(STLOrderPlateMatrix,1)
        e=STLOrderPlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end

fprintf(fid,'endsolid');
fclose(fid);



% Step #2, Write the plate surface part, using [SurfacePlatesNodes1] and [SurfacePlatesNodes2],


fid=fopen('PlateSurface.stl','w'); % STL file of the plate surface
fprintf(fid,'solid\n');

STLOrderSurfacePlateMatrix=...    % Each solid has 8 triangulars 
    [1	3	2;
    4	5	6;
    1	2	4;
    2	5	4;
    2	3	6;
    2	6	5;
    1	6	3;
    1	4	6];
for j=1:size(SurfacePlatesNodes1,3)
    Node_temp=SurfacePlatesNodes1(:,:,j);
    for i=1:size(STLOrderSurfacePlateMatrix,1)
        e=STLOrderSurfacePlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end

for j=1:size(SurfacePlatesNodes2,3)
    Node_temp=SurfacePlatesNodes2(:,:,j);
    for i=1:size(STLOrderSurfacePlateMatrix,1)
        e=STLOrderSurfacePlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end

fprintf(fid,'endsolid');
fclose(fid);



% Step #3, Write the hinges, including hinges (connecting two plates) and hinges (connecting three plates)

fid=fopen('Hinge.stl','w'); % STL file of the hinge connecting two surfaces
fprintf(fid,'solid\n');
STLOrderHingeTwoPlateMatrix=...    % Each solid has 6 faces and each face has 2 triangulars, totally 12 triangulars
    [6	2	1;
    1	5	6;
    5	1	3;
    3	7	5;
    7	3	4;
    4	8	7;
    8	4	2;
    2	6	8;
    3	1	2;
    2	4	3;
    5	7	8;
    8	6	5];

for j=1:size(HingeNodesTwo,3) % the hinge connecting two surfaces
    Node_temp=HingeNodesTwo(:,:,j);
    for i=1:size(STLOrderHingeTwoPlateMatrix,1)
        e=STLOrderHingeTwoPlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end

STLOrderHingeThreePlateMatrix=...  % top and bottom faces have 7+7 triangulars, 侧面 have 18 triangulars
    [3	1	17;
    17	19	3;
    15	17	11;
    11	9	15;
    19	15	7;
    7	5	19;
    19	17	15;
    4	18	2;
    18	4	20;
    16	12	18;
    12	16	10;
    20	8	16;
    8	20	6;
    20	16	18;
    1	3	4;
    3	19	20;
    19	5	6;
    5	7	8;
    7	15	16;
    15	9	10;
    9	11	12;
    11	17	18;
    17	1	2;
    4	2	1;
    20	4	3;
    6	20	19;
    8	6	5;
    16	8	7;
    10	16	15;
    12	10	9;
    18	12	11;
    2	18	17];

for j=1:size(HingeNodesThreeSorted,3) % the hinge connecting three surfaces
    Node_temp=HingeNodesThreeSorted(:,:,j);
    for i=1:size(STLOrderHingeThreePlateMatrix,1)
        e=STLOrderHingeThreePlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end

fprintf(fid,'endsolid');
fclose(fid);


% Step #4, Write the hinges surface, including hinges surface (connecting two plates) and hinges surface (connecting three plates)


fid=fopen('HingeSurface.stl','w'); % STL file of the hinge connecting two surfaces
fprintf(fid,'solid\n');
STLOrderSurfaceHingeTwoPlateMatrix=...    % Each solid has 6 faces and each face has 2 triangulars, totally 12 triangulars
    [6	2	1;
    1	5	6;
    5	1	3;
    3	7	5;
    7	3	4;
    4	8	7;
    8	4	2;
    2	6	8;
    3	1	2;
    2	4	3;
    5	7	8;
    8	6	5];

for j=1:size(SurfaceHingeNodesTwo1,3) % the hinge connecting two surfaces
    Node_temp=SurfaceHingeNodesTwo1(:,:,j);
    for i=1:size(STLOrderSurfaceHingeTwoPlateMatrix,1)
        e=STLOrderSurfaceHingeTwoPlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end

for j=1:size(SurfaceHingeNodesTwo2,3) % the hinge connecting two surfaces
    Node_temp=SurfaceHingeNodesTwo2(:,:,j);
    for i=1:size(STLOrderSurfaceHingeTwoPlateMatrix,1)
        e=STLOrderSurfaceHingeTwoPlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end




STLOrderSurfaceHingeThreePlateMatrix=... 
[23	39	3;
39	19	3;
19	39	25;
25	5	19;
24	4	40;
40	4	20;
20	26	40;
26	20	6;
3	19	20;
19	5	6;
5	25	26;
25	39	40;
39	23	24;
23	3	4;
20	4	3;
6	20	19;
26	6	5;
40	26	25;
24	40	39;
4	24	23;
31	37	11;
37	17	11;
17	37	21;
21	1	17;
32	12	38;
38	12	18;
18	22	38;
22	18	2;
11	17	18;
17	1	2;
1	21	22;
21	37	38;
37	31	32;
31	11	12;
18	12	11;
2	18	17;
22	2	1;
38	22	21;
32	38	37;
12	32	31;
27	35	7;
35	15	7;
15	35	29;
29	9	15;
28	8	36;
36	8	16;
16	30	36;
30	16	10;
7	15	16;
15	9	10;
9	29	30;
29	35	36;
35	27	28;
27	7	8;
16	8	7;
10	16	15;
30	10	9;
36	30	29;
28	36	35;
8	28	27];


for j=1:size(SurfaceHingeNodesThreeSorted,3) % the hinge connecting two surfaces
    Node_temp=SurfaceHingeNodesThreeSorted(:,:,j);
    for i=1:size(STLOrderSurfaceHingeThreePlateMatrix,1)
        e=STLOrderSurfaceHingeThreePlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end

fprintf(fid,'endsolid');
fclose(fid);

















%% Write STL into only two files

% Hard materials, only the plate

fid=fopen('HardMaterials.stl','w'); % STL file of the plate
fprintf(fid,'solid\n');

STLOrderPlateMatrix=...    % Each solid has 8 triangulars 
    [1	3	2;
    4	5	6;
    1	2	4;
    2	5	4;
    2	3	6;
    2	6	5;
    1	6	3;
    1	4	6];
for j=1:size(PlateNodes,3)
    Node_temp=PlateNodes(:,:,j);
    for i=1:size(STLOrderPlateMatrix,1)
        e=STLOrderPlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end

fprintf(fid,'endsolid');
fclose(fid);


% soft materials, the plate surface, and hinge, and hinge surface


fid=fopen('SoftMaterials.stl','w'); % STL file of the plate
fprintf(fid,'solid\n');
  
STLOrderSurfacePlateMatrix=...    % Each solid has 8 triangulars 
    [1	3	2;
    4	5	6;
    1	2	4;
    2	5	4;
    2	3	6;
    2	6	5;
    1	6	3;
    1	4	6];
for j=1:size(SurfacePlatesNodes1,3)
    Node_temp=SurfacePlatesNodes1(:,:,j);
    for i=1:size(STLOrderSurfacePlateMatrix,1)
        e=STLOrderSurfacePlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end

for j=1:size(SurfacePlatesNodes2,3)
    Node_temp=SurfacePlatesNodes2(:,:,j);
    for i=1:size(STLOrderSurfacePlateMatrix,1)
        e=STLOrderSurfacePlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end

STLOrderHingeTwoPlateMatrix=...    % Each solid has 6 faces and each face has 2 triangulars, totally 12 triangulars
    [6	2	1;
    1	5	6;
    5	1	3;
    3	7	5;
    7	3	4;
    4	8	7;
    8	4	2;
    2	6	8;
    3	1	2;
    2	4	3;
    5	7	8;
    8	6	5];

for j=1:size(HingeNodesTwo,3) % the hinge connecting two surfaces
    Node_temp=HingeNodesTwo(:,:,j);
    for i=1:size(STLOrderHingeTwoPlateMatrix,1)
        e=STLOrderHingeTwoPlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end

STLOrderHingeThreePlateMatrix=...  % top and bottom faces have 7+7 triangulars, 侧面 have 18 triangulars
    [3	1	17;
    17	19	3;
    15	17	11;
    11	9	15;
    19	15	7;
    7	5	19;
    19	17	15;
    4	18	2;
    18	4	20;
    16	12	18;
    12	16	10;
    20	8	16;
    8	20	6;
    20	16	18;
    1	3	4;
    3	19	20;
    19	5	6;
    5	7	8;
    7	15	16;
    15	9	10;
    9	11	12;
    11	17	18;
    17	1	2;
    4	2	1;
    20	4	3;
    6	20	19;
    8	6	5;
    16	8	7;
    10	16	15;
    12	10	9;
    18	12	11;
    2	18	17];

for j=1:size(HingeNodesThreeSorted,3) % the hinge connecting three surfaces
    Node_temp=HingeNodesThreeSorted(:,:,j);
    for i=1:size(STLOrderHingeThreePlateMatrix,1)
        e=STLOrderHingeThreePlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end


STLOrderSurfaceHingeTwoPlateMatrix=...    % Each solid has 6 faces and each face has 2 triangulars, totally 12 triangulars
    [6	2	1;
    1	5	6;
    5	1	3;
    3	7	5;
    7	3	4;
    4	8	7;
    8	4	2;
    2	6	8;
    3	1	2;
    2	4	3;
    5	7	8;
    8	6	5];

for j=1:size(SurfaceHingeNodesTwo1,3) % the hinge connecting two surfaces
    Node_temp=SurfaceHingeNodesTwo1(:,:,j);
    for i=1:size(STLOrderSurfaceHingeTwoPlateMatrix,1)
        e=STLOrderSurfaceHingeTwoPlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end

for j=1:size(SurfaceHingeNodesTwo2,3) % the hinge connecting two surfaces
    Node_temp=SurfaceHingeNodesTwo2(:,:,j);
    for i=1:size(STLOrderSurfaceHingeTwoPlateMatrix,1)
        e=STLOrderSurfaceHingeTwoPlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end




STLOrderSurfaceHingeThreePlateMatrix=... 
[23	39	3;
39	19	3;
19	39	25;
25	5	19;
24	4	40;
40	4	20;
20	26	40;
26	20	6;
3	19	20;
19	5	6;
5	25	26;
25	39	40;
39	23	24;
23	3	4;
20	4	3;
6	20	19;
26	6	5;
40	26	25;
24	40	39;
4	24	23;
31	37	11;
37	17	11;
17	37	21;
21	1	17;
32	12	38;
38	12	18;
18	22	38;
22	18	2;
11	17	18;
17	1	2;
1	21	22;
21	37	38;
37	31	32;
31	11	12;
18	12	11;
2	18	17;
22	2	1;
38	22	21;
32	38	37;
12	32	31;
27	35	7;
35	15	7;
15	35	29;
29	9	15;
28	8	36;
36	8	16;
16	30	36;
30	16	10;
7	15	16;
15	9	10;
9	29	30;
29	35	36;
35	27	28;
27	7	8;
16	8	7;
10	16	15;
30	10	9;
36	30	29;
28	36	35;
8	28	27];


for j=1:size(SurfaceHingeNodesThreeSorted,3) % the hinge connecting two surfaces
    Node_temp=SurfaceHingeNodesThreeSorted(:,:,j);
    for i=1:size(STLOrderSurfaceHingeThreePlateMatrix,1)
        e=STLOrderSurfaceHingeThreePlateMatrix(i,:);
        NV=local_find_normal(Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),0); % Normal vector of the triangular
        local_write_facet(fid,Node_temp(e(1),:),Node_temp(e(2),:),Node_temp(e(3),:),NV)
    end
end

fprintf(fid,'endsolid');
fclose(fid);





%% Subfunctions


function [Ntub,Ntub2D]=SmoothTriangularPanel(Ntub,Ntub2D,assemPoly)
for ifa=1:length(Ntub.face)
vertices =Ntub.node(Ntub.face{ifa},:);  % Example triangle with vertices at (0,0,0), (1,0,0), and (0,1,0)
resultSharpAngles(ifa) = judgeSharpAngles(vertices(1,:), vertices(2,:), vertices(3,:));
end
kj=1;
for jj=1:length(assemPoly)
    for as=1:length(assemPoly(jj).polyfaceIdx)/2
    assFacePair(kj,:)=[assemPoly(jj).polyfaceIdx(2*as-1)  assemPoly(jj).polyfaceIdx(2*as)];
    kj=kj+1;
    end
end
assFacePair=unique(sort(assFacePair,2),'rows');
sharpFaceIdx=find(resultSharpAngles==1);
uniRow=[];
for ic=1:  length(sharpFaceIdx)
[irr0,icc]=find(assFacePair==sharpFaceIdx(ic));
uniRow=[uniRow, irr0];
end
uniRow=unique(uniRow);
for irr1=1:length(uniRow)
    irr=uniRow(irr1);
mutualHinge=intersect(Ntub.face{assFacePair(irr,1)},Ntub.face{assFacePair(irr,2)});
face1Rep=[Ntub.face{assFacePair(irr,1)};...
    Ntub.face{assFacePair(irr,1)}(2)  Ntub.face{assFacePair(irr,1)}(3) Ntub.face{assFacePair(irr,1)}(1);...
    Ntub.face{assFacePair(irr,1)}(3)  Ntub.face{assFacePair(irr,1)}(1) Ntub.face{assFacePair(irr,1)}(2);];

face2Rep=[Ntub.face{assFacePair(irr,2)};...
    Ntub.face{assFacePair(irr,2)}(2)  Ntub.face{assFacePair(irr,2)}(3) Ntub.face{assFacePair(irr,2)}(1);...
    Ntub.face{assFacePair(irr,2)}(3)  Ntub.face{assFacePair(irr,2)}(1) Ntub.face{assFacePair(irr,2)}(2);];

face1=Ntub.face{assFacePair(irr,1)};    face2=Ntub.face{assFacePair(irr,2)};

remPface1=setxor(face1,mutualHinge);
remPface2=setxor(face2,mutualHinge);

Ntub.face{assFacePair(irr,1)}=[remPface1, face2Rep((face2Rep(:,2)==remPface2),1), remPface2];
Ntub.face{assFacePair(irr,2)}=[remPface2, face1Rep((face1Rep(:,2)==remPface1),1), remPface1];
Ntub2D.face{assFacePair(irr,1)}=[remPface1, face2Rep((face2Rep(:,2)==remPface2),1), remPface2];
Ntub2D.face{assFacePair(irr,2)}=[remPface2, face2Rep((face1Rep(:,2)==remPface1),1), remPface1];
end
%  rebuild edges;
Ntub.edge=[];  
Ntub.edge=unique(sort(buildEdgeFromFace(Ntub.face),2),'rows');
Ntub2D.edge=[]; 
Ntub2D.edge=unique(sort(buildEdgeFromFace(Ntub2D.face),2),'rows');
end



function result = judgeSharpAngles(A, B, C)

    % Calculate the lengths of the sides of the triangle
    AB = norm(A - B);
    BC = norm(B - C);
    CA = norm(C - A);

    % Calculate the angles of the triangle using the cosine rule
    angleA = acos((BC^2 + CA^2 - AB^2) / (2 * BC * CA));
    angleB = acos((CA^2 + AB^2 - BC^2) / (2 * CA * AB));
    angleC = acos((AB^2 + BC^2 - CA^2) / (2 * AB * BC));

    % Store the angles in a vector
    angles = [angleA, angleB, angleC];

    % Sort the angles
    sortedAngles = sort(angles);

    % Check if the sum of the two smallest angles is less than the half of the largest angle
    if sum(sortedAngles(1:2)) < 0.7*sortedAngles(3)
        result = true;
    else
        result = false;
    end

end

 function edge=buildEdgeFromFace(face)
 for i=1:length(face)
     edge(3*i-2,:)= [ face{i}(1)   face{i}(2) ]; 
     edge(3*i-1,:)= [ face{i}(2)   face{i}(3) ]; 
     edge(3*i,:)= [ face{i}(3)   face{i}(1) ]; 
 end
 end


 %to generate the normal direction vector n, from the vertex of triangular
function n = local_find_normal(p1,p2,p3,opt)
if ~opt
v2 = p2-p1;
v1 = p3-p1;
v3 = cross(v1,v2);
n = v3 ./ sqrt(sum(v3.*v3));
else
    n=opt;
end
end



% to offset these triangulars by 't'
function [f1,f2,f3] = OffsetTriangular(a1,a2,a3,t)
a4=(a1+a2+a3)./3;

b1=(a2+a3)./2;
b2=(a3+a1)./2;
b3=(a1+a2)./2;


c1=(a4-b1)./norm(a4-b1).*t;
c2=(a4-b2)./norm(a4-b2).*t;
c3=(a4-b3)./norm(a4-b3).*t;

d1=b1+c1;
d2=b2+c2;
d3=b3+c3;


e1=(d2+d3)./2;
e2=(d3+d1)./2;
e3=(d1+d2)./2;


f1=e1+(e1-d1);
f2=e2+(e2-d2);
f3=e3+(e3-d3);

end


function [b1,b2,b3,c1,c2,c3] = Triangular2Plates(a1,a2,a3,normalvector1,normalvector2)
b1=a1+normalvector1;
b2=a2+normalvector1;
b3=a3+normalvector1;

c1=a1+normalvector2;
c2=a2+normalvector2;
c3=a3+normalvector2;
end




%to write a single small triangular to the file
function local_write_facet(fid,p1,p2,p3,opt)
    n = local_find_normal(p1,p2,p3,opt);

        fprintf(fid,'facet normal %.7E %.7E %.7E\r\n', n(1),n(2),n(3) );
        fprintf(fid,'outer loop\r\n');
        fprintf(fid,'vertex %.7E %.7E %.7E\r\n', p1);
        fprintf(fid,'vertex %.7E %.7E %.7E\r\n', p2);
        fprintf(fid,'vertex %.7E %.7E %.7E\r\n', p3);
        fprintf(fid,'endloop\r\n');
        fprintf(fid,'endfacet\r\n');

end



















%% 生成 STLOrderHingeThreePlateMatrix 用的
% % order=[1 3 19 5 7 15 9 11 17 1];
% % 
% % triangular1=zeros(length(order)-1,3);
% % for i=1:length(order)-1
% %     triangular1(i,:)=[order(i),order(i+1),order(i+1)+1];
% % end
% % 
% % 
% % order=[1 3 19 5 7 15 9 11 17 1]+1;
% % 
% % triangular2=zeros(length(order)-1,3);
% % for i=1:length(order)-1
% %     triangular2(i,:)=[order(i+1),order(i),order(i)-1];
% % end
