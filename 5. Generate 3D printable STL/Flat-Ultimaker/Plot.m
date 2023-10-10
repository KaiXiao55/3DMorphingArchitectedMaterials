function  Plot(Polyhedron,viewPoint,Text,fText,object,numUnitPoly,trans,lightOn,transVer)
if isempty(transVer)==0
    for i=1:length(Polyhedron)
            Polyhedron(i).node=Polyhedron(i).node+transVer;
    end
end
% Loop from i-th polyhedron
for i=1:length(Polyhedron)
    % Select color for every polyhedron
    if strcmp(object,'polyhedra')==1
        switch mod(i,numUnitPoly)
            case 0
                %                 facecolor=[250 240 230]./256;
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
    if strcmp(object,'interPolyhedra')==1
        facecolor=[199,21,133]/256;
    elseif strcmp(object,'Snapology')==1
        excolor= [147,112,219]/256 ;   % Original
%           incolor=excolor;
        incolor= [255 193 193]./256;    % Original
        %         incolor=excolor;
        Lightpoint=2*[cosd(viewPoint(2))*sind(viewPoint(1)), -cosd(viewPoint(2))*cosd(viewPoint(1)), sind(viewPoint(2))];
    elseif  strcmp(object,'FlatState')==1
        facecolor=[255 240 245]./256;
    end
    if strcmp(Text,'text')==1
        for kk=1:size(Polyhedron(i).node,1)  % text the nodes
            text(Polyhedron(i).node(kk,1),Polyhedron(i).node(kk,2),Polyhedron(i).node(kk,3),num2str(kk),'fontsize',12,'color','k');
        end
        bodyCenter=sum(Polyhedron(i).node)/size(Polyhedron(i).node,1);
        text(bodyCenter(1),bodyCenter(2),bodyCenter(3),num2str(i),'fontsize',16,'color','r');
    end
%     Polyhedron(i).face
%     length(Polyhedron(i).face)
    for j=1:length(Polyhedron(i).face)
        if isempty(Polyhedron(i).face{j})==0
            if strcmp(fText,'ftext')==1
                indexX=sum(Polyhedron(i).node(:,1))/size(Polyhedron(i).node,1);
                indexY=sum(Polyhedron(i).node(:,2))/size(Polyhedron(i).node,1);
                indexZ=sum(Polyhedron(i).node(:,3))/size(Polyhedron(i).node,1);
                text(indexX, indexY, indexZ,num2str(i),'fontsize',14,'Color','b');  % indicate j-th polyhedron this is
                faceCenter=sum(Polyhedron(i).node([Polyhedron(i).face{j}(1:end)],:))/length(Polyhedron(i).face{j});
                text( faceCenter(1), faceCenter(2), faceCenter(3),num2str(j),'fontsize',12,'Color','m');

            end
            if strcmp(object,'polyhedra')==1 || strcmp(object,'FlatState')
                    transp=0.1; %  edgecolor='w';  
                    edgecolor=[192,192,192]/256;  
                    lineWidth=0.01;   %  transp=0.14;
                    edgealpha=0.1;
                patch('Faces',Polyhedron(i).face{j},'Vertices',Polyhedron(i).node,'FaceColor',facecolor,...
                    'facealpha',trans,'EdgeColor',edgecolor,'LineWidth',lineWidth,...
                    'EdgeAlpha',edgealpha); % EdgeColor [192,192,192]/256
            elseif strcmp(object,'Snapology')==1
                nodes=Polyhedron(i).node;  faces=Polyhedron(i).face{j};
%                 j
                faNormal(j,:)=cross((nodes(faces(3),:)-nodes(faces(2),:)),(nodes(faces(1),:)-nodes(faces(2),:)))...
                    /norm(cross((nodes(faces(3),:)-nodes(faces(2),:)),(nodes(faces(1),:)-nodes(faces(2),:))));
                if dot(faNormal(j,:),Lightpoint)>0
                    patch('Faces',Polyhedron(i).face{j},'Vertices',Polyhedron(i).node,'FaceColor',incolor,'facealpha',trans,'LineWidth',0.01,'EdgeColor','k','edgealpha',1); %,,'edgealpha',0.1
                elseif dot(faNormal(j,:),Lightpoint)<=0
                    patch('Faces',Polyhedron(i).face{j},'Vertices',Polyhedron(i).node,'FaceColor',excolor,'facealpha',trans,'LineWidth',0.01,'EdgeColor','k','edgealpha',1); %,'k','edgealpha',0.1
                end
            end
        end
    end
end
% End loop
view(viewPoint(1),viewPoint(2))
if strcmp(lightOn,'on')==1
axis equal
axis off
grid off

if strcmp(object,'Snapology')==1
    light('color','w','style','infinite','position',Lightpoint)
    lighting flat
end
material dull
set(gcf,'color','w')
end
end