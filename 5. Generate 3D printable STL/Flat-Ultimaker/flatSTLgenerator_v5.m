
clear, close all
Case='Shape_Lockable Sphere';
% Case='singleUnitRegularunitWith2D'; 
load([Case,'.mat'])
% Ntub2D.face=Ntub2D.face';
%% If only use single unit in a assembly
% % singlePolyIdx=[1];
% % polyNodeIdx=assemPoly(singlePolyIdx).nodeInMergeIdx;
% % polyfaceIdx=assemPoly(singlePolyIdx).polyfaceIdx;
% % % Ntub.node=Ntub.node(polyNodeIdx,:);
% % Ntub.face=Ntub.face(polyfaceIdx);
% % [Ntub.node,Ntub.face]=reNumberPartialFaceSnapology(Ntub.node,Ntub.face');
% % Ntub.edge=buildEdgeFromFace(Ntub.face);
% % 
% % Ntub2D.face=Ntub2D.face(polyfaceIdx);
% % [Ntub2D.node,Ntub2D.face]=reNumberPartialFaceSnapology(Ntub2D.node,Ntub2D.face');
% % Ntub2D.edge=buildEdgeFromFace(Ntub2D.face);
% % 
% Ntub.face=Ntub.face';
% Ntub2D.face=Ntub2D.face';
%% assign basic geometry
thickPLA = 1.5;   % mm; 
thickTPU = 0.5;   % mm
thickIso = 0.5;  % mm ; thickness of TPU isolation layer
d_re = 1;   % retract distance of PLA plate from the edge
d_vh = 0.6; % width of vertical hinge

scaleFactor =1;

%% scale the original structure to printable size

Ntub.node = scaleFactor*Ntub.node;
Ntub2D.node = scaleFactor*Ntub2D.node;


%% identify the printing layer of each global face
UnitLayer = [1;2;2;1;3;4;4;3];    % UnitLayer(i)=j ----  unit-i is in j-layer

%  UnitLayer = [1;2;2;1;3;4;4;3;5;5];  % 2  2  3 Full units
% LayerUnits{1}=[1 3 5 7 9]; LayerUnits{2}=[2 4 6 8]; LayerUnits{3}=[10 12 14 16 18];
% LayerUnits{4}=[11 13 15 17];   LayerUnits{5}=[19 20 21 22 23]; 
%  UnitLayer =[];
% for i=1:length(LayerUnits)
%      UnitLayer(LayerUnits{i})=i*ones(length(LayerUnits{i}),1);
% end 

upper = [7, 8, 9, 10, 11, 12, 13, 14, 19, 20];   % 
FaceLayer = zeros( size(Ntub2D.face, 2),1);   % FaceLayer(i)=j ---- global face-i is in the j-th layer

for i = 1:length(UnitLayer)
    for j = 1:length(assemPoly(i).polyfaceIdx)
        gfaceIdx = assemPoly(i).polyfaceIdx(j);    % global index of face in (j-th face of i-th unit)
        if ismember(j,upper)                       % check if it is in the upper layer of each unit
            FaceLayer(gfaceIdx) = UnitLayer(i)+1;  % upper face of its unit - global layer=unit layer+1
        else
            FaceLayer(gfaceIdx) = UnitLayer(i);    % lower face of its unit - global layer=unit layer
        end
        
    end
end

%%

upper = [7, 8, 9, 10, 11, 12, 13, 14, 19, 20];

widthScaleFactor = 0.85;
plateThickness = 0.1; % (Fixed) thickness of plate.
upwardScaleFactor = 0.5; % Upper half of the thickness of hinge (as ratio to half plateThickness).
downwardScaleFactor = 0.5; % Lower half of the thickness of hinge (as ratio to half plateThickness).

z1 = 0.5 * plateThickness; % Lower reference plane.
z2 = 2 * plateThickness; % Upper reference plane. 0.5-plateThickness gap between upper and lower plates.
%%  Visualize 
viewPoint=[-165 25];
figure(1)
Ntub.node=Ntub.node-sum(Ntub.node(:,2),1)/size(Ntub.node,1);
Ntub.node=Ntub.node-[ 0 0  min(Ntub.node(:,3))];
Plot(Ntub,viewPoint,'texxt','ftexxt','Snapology',4,1,'on',[])
% axis on
% axis equal
figure(2)
% hold on

Ntub2D.node=Ntub2D.node-sum(Ntub2D.node(:,2),1)/size(Ntub2D.node,1);
% Ntub2D.node=Ntub2D.node-
Plot(Ntub2D,viewPoint,'texxt','ftxext','Snapology',4,1,'on',[])
% axis on
axis equal
figure(3)
Plot(assemPoly,[-23 30],'texxt','ftext','Snapology',4,1,'on',[])
axis on
axis equal
figure(4)
Plot(assem,[-23 30],'texxt','ftext','Snapology',4,0.7,'on',[])
axis on
axis equal
%% Find edges as hinges between faces.
hingeEx = [];
for i = 1 : 1 : size(Ntub2D.face, 2)
    for j = i + 1 : 1 : size(Ntub2D.face, 2)
        h = intersect(Ntub2D.face{i}, Ntub2D.face{j}); % Find the nodes shared by each 2 faces.
        if length(h) == 2 % If a pair of faces share 2 nodes, they are connected by a hinge.
            hingeEx = [hingeEx; h, i, j, 0]; % Each row of hingeEx takes the form of [n1, n2, i, j, type], where n1, n2 are the nodes of the hinge, and i, j are the indices of the two faces connected.
        end
    end
end


%% classify hinges
for i = 1:size(hingeEx,1)
    if hingeEx(i,5)~=0 
        continue
    end
    k = 0;
    indx = [];
    for j = i+1:size(hingeEx,1)
        if hingeEx(i,1)==hingeEx(j,1) & hingeEx(i,2)==hingeEx(j,2)   % check if this hinge has been used for several times
            k = k+1;
            indx(k) = j;
        end
    end
    
    if k==0 % this hinge is connecting 2 faces; it might be flat(1) or vertical(2) hinge     
        faceID1 = hingeEx(i,3);
        faceID2 = hingeEx(i,4);        
        if FaceLayer(faceID1)==FaceLayer(faceID2)
            hingeEx(i,5) = 1;       % hingeEx(:,5) refers to the type of hinge; 1-flat;2-vertical;3-T shape
        else
            hingeEx(i,5) = 2;
        end
    else     % this hinge is T-shape-hinge(3) connecting 3 or more faces
        hingeEx(indx,5) = 3;
        hingeEx(i,5) = 3;
    end
end

%% Find spare edges (not shared between faces).
spEdge = [];
for i = 1 : 1 : size(Ntub2D.edge, 1)
    sp = Ntub2D.edge(i, :); % For each edge of the structure.
    if ~ismember(sp, hingeEx(:, 1 : 2), 'rows') % If the edge is not a hinge, it is a spare edge.
        j = 1;
        h = intersect(sp, Ntub2D.face{j});
        while length(h) ~= 2
            j = j + 1;
            h = intersect(sp, Ntub2D.face{j});
        end % Find the index of the face that contains the spare edge.
        spEdge = [spEdge; sp, j]; % Each row of spEdge takes the form of [n1, n2, j], where n1, n2 are the nodes of the spare edge and j is the index of the face that contains the edge.
    end
end


%% generate the PLA plate 
platePLAInfo = cell(size(Ntub2D.face, 2), 1);

for i = 1:size(Ntub2D.face, 2)
    PLAnode = zeros(6,3);
    % calculate (x,y) coordinates
    [spIn, spLoc] = ismember(i, spEdge(:, 3));
    for j = 1:3
        nodeIDc = Ntub2D.face{i}(  mod(j +3-1,3)+1  ) ;    % node index of central point
        nodeIDs1 = Ntub2D.face{i}(  mod(j+1 +3-1,3)+1  ); % node index of side point 1
        nodeIDs2 = Ntub2D.face{i}(  mod(j-1 +3-1,3)+1  ); % node index of side point 2
        v1 = (Ntub2D.node(nodeIDs1,1:2)-Ntub2D.node(nodeIDc,1:2))/norm(Ntub2D.node(nodeIDs1,1:2)-Ntub2D.node(nodeIDc,1:2));
        v2 = (Ntub2D.node(nodeIDs2,1:2)-Ntub2D.node(nodeIDc,1:2))/norm(Ntub2D.node(nodeIDs2,1:2)-Ntub2D.node(nodeIDc,1:2));
        sin_theta = norm(cross([v1,0],[v2,0]));
        
        if spLoc==0   % this face doesn't contain spare edge, so that all edges of triangle need to be retracted
            Pi = Ntub2D.node(nodeIDc,1:2) + d_re/sin_theta*(v1+v2);
        else    % this face has spare edge, so spare edge doesn'e need to be retracted
            sp = spEdge(spLoc, 1 : 2);    % nodes of spare edge
            if ismember(nodeIDc, sp)   % this node is located in the spare edge
                nodeIDs = setdiff(sp,nodeIDc);  % the other node of the spare edge
                if nodeIDs==nodeIDs1
                    Pi = Ntub2D.node(nodeIDc,1:2) + d_re/sin_theta*(v1);  % along v1 direction
                else
                    Pi = Ntub2D.node(nodeIDc,1:2) + d_re/sin_theta*(v2);  % along v2 direction
                end
            else   % this node is not located in the spare edge
                Pi = Ntub2D.node(nodeIDc,1:2) + d_re/sin_theta*(v1+v2);
            end
            
        end
        
        PLAnode(j,1:2) = Pi;
        PLAnode(j+3,1:2) = Pi;
        NodeGIndex(j) = nodeIDc;  % global index of this node   
  
    end
    
    % calculate z-coordinate    
    FaceHeight(i,1) = thickPLA/2 + ( FaceLayer(i)-1 )*(thickPLA + thickIso);    % z-coordinate of the central plane of each layer
    FaceHeight(i,2) = FaceHeight(i,1) - thickPLA/2;                             % z-coordinate of the lower surface of each layer
    FaceHeight(i,3) = FaceHeight(i,1) + thickPLA/2;                             % z-coordinate of the upper surface of each layer
   
    PLAnode(1:3,3) = FaceHeight(i,2);       % lower surface
    PLAnode(4:6,3) = FaceHeight(i,3);       % upper surface
    
    platePLAInfo{i} = struct('NodeGIndex', NodeGIndex, 'PLAnode', PLAnode);
      
end




%% generate the flat TPU plate and vertical hinge 
plateTPUInfo = cell(size(Ntub2D.face, 2), 1);
HingeVertInfo = cell(size(Ntub2D.face, 2), 1);
verHinge = hingeEx(  find(hingeEx(:,5)==2)  ,:);

for i = 1:size(Ntub2D.face, 2)
    TPUnode = zeros(6,3);
    HingeVert = zeros(4,3);
    [k1,k2] = find(verHinge(:, 3:4)==i);   % check if any vertical hinge connect i-face
    
    FaceHeight(i,1) = thickPLA/2 + ( FaceLayer(i)-1 )*(thickPLA + thickIso);    % z-coordinate of the central plane of each layer
    FaceHeight(i,2) = FaceHeight(i,1) - thickTPU/2;                             % z-coordinate of the lower surface of each layer
    FaceHeight(i,3) = FaceHeight(i,1) + thickTPU/2;                             % z-coordinate of the upper surface of each layer
    
    if isempty(k1)  % this face doesn't contain vertical hinge, three nodes of triangular TPU plate share the same coordinates with original faces
        TPUnode(1:3,:) = Ntub2D.node(Ntub2D.face{i},:);
        TPUnode(4:6,:) = Ntub2D.node(Ntub2D.face{i},:);
        NodeGIndex = Ntub2D.face{i};
    else    % this face contains vertical hinge
        P1Indx = verHinge(k1,1);                            % two nodes that lie in the vertical hinge
        P2Indx = verHinge(k1,2);
        P3Indx = setdiff(Ntub2D.face{i},[P1Indx,P2Indx]);   % the node that doesn't lie in the vertical hinge
        
        for j = 1:3
            nodeIDc = Ntub2D.face{i}(  mod(j +3-1,3)+1  ) ;    % node index of central point
            if nodeIDc==P3Indx
                TPUnode(j,:) = Ntub2D.node(nodeIDc,:);
%                 TPUnode(j,:) = Ntub2D.node(Ntub2D.face{i}(j),:);
            else     % this node is located in the vertical hinge
                k = find(verHinge(k1,1:2)==nodeIDc);
                nodeIDs1 = Ntub2D.face{i}(  mod(j+1 +3-1,3)+1  ); % node index of side point 1
                nodeIDs2 = Ntub2D.face{i}(  mod(j-1 +3-1,3)+1  ); % node index of side point 2
                v1 = (Ntub2D.node(nodeIDs1,1:2)-Ntub2D.node(nodeIDc,1:2))/norm(Ntub2D.node(nodeIDs1,1:2)-Ntub2D.node(nodeIDc,1:2));
                v2 = (Ntub2D.node(nodeIDs2,1:2)-Ntub2D.node(nodeIDc,1:2))/norm(Ntub2D.node(nodeIDs2,1:2)-Ntub2D.node(nodeIDc,1:2));
                sin_theta = norm(cross([v1,0],[v2,0]));

                % calculate the (x,y) coordinate
                if nodeIDs1==P3Indx
                    TPUnode(j,1:2) = Ntub2D.node(nodeIDc,1:2) - d_re/sin_theta*(v1);  % along v1 direction
                    HingeVert(k,1:2) =  Ntub2D.node(nodeIDc,1:2) - d_re/sin_theta*(v1);
                    HingeVert(k+2,1:2) =  Ntub2D.node(nodeIDc,1:2) - (d_vh/2)/sin_theta*(v1);
                else
                    TPUnode(j,1:2) = Ntub2D.node(nodeIDc,1:2) - d_re/sin_theta*(v2);  % along v2 direction
                    HingeVert(k,1:2) =  Ntub2D.node(nodeIDc,1:2) - d_re/sin_theta*(v2);
                    HingeVert(k+2,1:2) =  Ntub2D.node(nodeIDc,1:2) - (d_vh/2)/sin_theta*(v2);
                end
                
                % calculate the z-coordinate of vertical hinge
                
                face1 = verHinge(k1,k2+2);   % current face
                face2 = verHinge(k1,setdiff([1,2],k2)); 
                
                if FaceLayer(face1)<FaceLayer(face2)    % this face is the lower surface of the vertical hinge
                    HingeVert(k,3) = FaceHeight(i,3);
                    HingeVert(k+2,3) = FaceHeight(i,3);
                else         % this face is the upper surface of the vertical hinge
                    HingeVert(k,3) = FaceHeight(i,2);
                    HingeVert(k+2,3) = FaceHeight(i,2);
                end
            end
            
            TPUnode(j+3,:) = TPUnode(j,:);
            NodeGIndex(j) = nodeIDc;  % global index of this node
            
            %             figure(5)
            %             plot(Ntub2D.node(nodeIDc,1),Ntub2D.node(nodeIDc,2),'ko');
%             hold on
%             axis equal
%             plot(TPUnode(j,1),TPUnode(j,2),'ro');
        end

    end
    % calculate z-coordinate of TPU plate
    TPUnode(1:3,3) = FaceHeight(i,2);       % lower surface
    TPUnode(4:6,3) = FaceHeight(i,3);       % upper surface
    
    plateTPUInfo{i} = struct('NodeGIndex', NodeGIndex, 'TPUnode', TPUnode);
    HingeVertInfo{i} = struct('HingeVert', HingeVert);
end




%% generate T shape hinge (type-3)
T_Hinge = hingeEx(  find(hingeEx(:,5)==3)  ,:);
[TH_stored II] = sortrows(T_Hinge);
T_Hinge = T_Hinge(II,:);

HingeTInfo = cell(size(T_Hinge,1), 1);
for i = 1:3:size(T_Hinge,1)
   % check which face is the tip of T
   if FaceLayer(T_Hinge(i,3))==FaceLayer(T_Hinge(i,4))
       FaceTip = T_Hinge(i+1,3);
   else 
       FaceTip = T_Hinge(i,4);
   end
   
   % extend the TPU plate of tip face at T-hinge
   P1Indx = T_Hinge(i,1);                            % two nodes that lie in the vertical hinge
   P2Indx = T_Hinge(i,2);
   P3Indx = setdiff(Ntub2D.face{FaceTip},[P1Indx,P2Indx]);   % the rest of node that doesn't lie in the vertical hinge of tip face
  
   % extend the first point
   nodeIDc = P1Indx;
   nodeIDs1 = P2Indx; % node index of side point 1
   nodeIDs2 = P3Indx; % node index of side point 2
   v1 = (Ntub2D.node(nodeIDs1,1:2)-Ntub2D.node(nodeIDc,1:2))/norm(Ntub2D.node(nodeIDs1,1:2)-Ntub2D.node(nodeIDc,1:2));
   v2 = (Ntub2D.node(nodeIDs2,1:2)-Ntub2D.node(nodeIDc,1:2))/norm(Ntub2D.node(nodeIDs2,1:2)-Ntub2D.node(nodeIDc,1:2));
   sin_theta = norm(cross([v1,0],[v2,0]));
   T_node(1,1:2) = Ntub2D.node(nodeIDc,1:2) - (d_vh/2)/sin_theta*(v2);
   T_node(2,1:2) = Ntub2D.node(nodeIDc,1:2) + (d_vh/2)/sin_theta*(v2);
   
   % extend the second point
   nodeIDc = P2Indx;
   nodeIDs1 = P1Indx; % node index of side point 1
   nodeIDs2 = P3Indx; % node index of side point 2
   v1 = (Ntub2D.node(nodeIDs1,1:2)-Ntub2D.node(nodeIDc,1:2))/norm(Ntub2D.node(nodeIDs1,1:2)-Ntub2D.node(nodeIDc,1:2));
   v2 = (Ntub2D.node(nodeIDs2,1:2)-Ntub2D.node(nodeIDc,1:2))/norm(Ntub2D.node(nodeIDs2,1:2)-Ntub2D.node(nodeIDc,1:2));
   sin_theta = norm(cross([v1,0],[v2,0]));
   T_node(3,1:2) = Ntub2D.node(nodeIDc,1:2) - (d_vh/2)/sin_theta*(v2);
   T_node(4,1:2) = Ntub2D.node(nodeIDc,1:2) + (d_vh/2)/sin_theta*(v2);
   
   T_node(5:8,1:2) =  T_node(1:4,1:2);
   
   
   % identify the z-coordinate
   FaceHeight(i,1) = thickPLA/2 + ( FaceLayer(FaceTip)-1 )*(thickPLA + thickIso);    % z-coordinate of the central plane of tip layer
   FaceHeight(i,2) = FaceHeight(i,1) - thickTPU/2;                             % z-coordinate of the lower surface of tip layer
   FaceHeight(i,3) = FaceHeight(i,1) + thickTPU/2;                             % z-coordinate of the upper surface of tip layer
   if FaceLayer(FaceTip)>FaceLayer(T_Hinge(i,3))  % reverse-T shape
        T_node(5:8,3) = FaceHeight(i,2);   % upper surface of T hinge
        T_node(1:4,3) = T_node(5:8,3)-(thickPLA + thickIso);  % lower surface of T hinge
   else
       T_node(1:4,3) = FaceHeight(i,3);   % upper surface of T hinge
       T_node(5:8,3) = T_node(1:4,3)+(thickPLA + thickIso);  % lower surface of T hinge
   end
   
  HingeTInfo{i} = struct('T_node', T_node);
    
end


%% plot the generated PLA plate  
figure(6)
fac_tri = [1 2 3;4 5 6];
fac_vertface = [1 2 5 4;2 5 6 3;1 4 6 3;];
PLA_vertice=[];
PLA_faceUnit={[1 2 3;],[4 5 6;], [1 2 5 4;], [2 5 6 3;], [1 4 6 3;], };
k=1;
for i = 1:size(Ntub2D.face, 2)
% for i = 1:3
    vert = platePLAInfo{i}.PLAnode;
    patch('Faces', fac_tri, 'Vertices', vert, 'FaceColor', 'white', 'FaceAlpha', 1)
    patch('Faces', fac_vertface, 'Vertices', vert, 'FaceColor', 'white', 'FaceAlpha', 1)
    %  Write vertice and faces into .obj file
    PLA_vertice=[PLA_vertice; vert;];
    for j=1:length(PLA_faceUnit)
    PLA_face{k}=PLA_faceUnit{j}+(i-1)*size(platePLAInfo{i}.PLAnode,1);
    k=k+1;
    end
end
Write2OBJ('PLA', PLA_vertice, PLA_face')
view(3)
axis equal


%% plot the generated TPU plate
figure(6)
fac_tri = [1 2 3;4 5 6];
fac_vertface = [1 2 5 4;2 5 6 3;1 4 6 3];

TPU_vertice=[];
TPU_faceUnit={[1 2 3;],[4 5 6;], [1 2 5 4;], [2 5 6 3;], [1 4 6 3;], };
k=1;
for i = 1:size(Ntub2D.face, 2)
% for i = 1:2
    vert = plateTPUInfo{i}.TPUnode;
    patch('Faces', fac_tri, 'Vertices', vert, 'FaceColor', 'b', 'FaceAlpha', 1)
    patch('Faces', fac_vertface, 'Vertices', vert, 'FaceColor', 'b', 'FaceAlpha', 1)
    TPU_vertice=[TPU_vertice; vert;];
    for j=1:length(TPU_faceUnit)
    TPU_faceT1{k}=TPU_faceUnit{j}+(i-1)*size( plateTPUInfo{i}.TPUnode,1);
    k=k+1;
    end
end
view(3)
axis equal

%% plot the vertical hinge (type-2)
figure(6)
fac = [1 2 4 3;5 6 8 7;1 2 6 5;2 4 8 6;4 8 7 3;1 3 7 5];
sizeTPUver=size(TPU_vertice,1);
% k=1;
for i = 1:size(verHinge,1)
    face1 = verHinge(i,3);
    face2 = verHinge(i,4);
    vert = [HingeVertInfo{face1}.HingeVert;HingeVertInfo{face2}.HingeVert];
%     % for write obj
    TPU_vertice=[TPU_vertice; vert;];
    patch('Faces', fac, 'Vertices', vert, 'FaceColor', 'b', 'FaceAlpha', 1)
    for j=1:size(fac,1)
    TPU_faceT1{k}=fac(j,:)+(i-1)*size(vert,1)+sizeTPUver;
    k=k+1;
    end
end

%% plot the T hinge (type-3)
figure(6)
fac = [1 2 4 3;5 6 8 7;1 2 6 5;2 4 8 6;4 8 7 3;1 3 7 5];
sizeTPUver=size(TPU_vertice,1);
% k=1;
ii=1;
for i = 1:3:size(T_Hinge,1) 
    vert =  HingeTInfo{i}.T_node;
    TPU_vertice=[TPU_vertice; vert;];
    patch('Faces', fac, 'Vertices', vert, 'FaceColor', 'b', 'FaceAlpha', 1);
    for j=1:size(fac,1)
        TPU_faceT1{k}=fac(j,:)+(ii-1)*size(vert,1)+sizeTPUver;
        k=k+1;
    end
     ii=ii+1;
end
Write2OBJ('TPU', TPU_vertice, TPU_faceT1')


% save TPU and PLA together
TPUandPLA_vertice=[TPU_vertice; PLA_vertice;];
sizeVerTpu=size(TPU_vertice,1);
for  i=1:length(PLA_face)
    PLA_face2{i}=PLA_face{i}+sizeVerTpu;
end
TPUandPLA_face=[TPU_faceT1';  PLA_face2';];
Write2OBJ('TPUandPLA', TPUandPLA_vertice, TPUandPLA_face)

% --------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------- Sub-Function -----------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------
function [newNode,newFace]=reNumberPartialFaceSnapology(node,face)
faceNode=[];
for i=1:size(face,1)
    faceNode=[faceNode;face{i},];
end
newNodeIdx=unique(faceNode);
newFace=face; newNode=node(newNodeIdx,:);
% re-number the face
for i=1:size(face,1)
    for j=1:length(face{i})
        newFace{i}(j)=find(newNodeIdx==face{i}(j));
    end
end
% for i=1:size(edge,1)
%     for j=1:length(edge(i,:))
%         newEdge(i,j)=find(newNodeIdx==edge(i,j));
%     end
% end
end
 function edge=buildEdgeFromFace(face)
 for i=1:length(face)
     edge(3*i-2,:)= [ face{i}(1)   face{i}(2) ]; 
     edge(3*i-1,:)= [ face{i}(2)   face{i}(3) ]; 
     edge(3*i,:)= [ face{i}(3)   face{i}(1) ]; 
 end
 edge=unique(sort(edge,2),'rows');
  end