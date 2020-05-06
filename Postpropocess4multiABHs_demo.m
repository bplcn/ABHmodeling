% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%     This code generate plate with an acoustic black hole(ABH), the element 
% type can be solid(C3D8), continous shell(SC8R) and common thick shell(S4).
% The damp rubber is modeled directly.
%                                               Bairyou Ben 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% v1.01: the damp can be added on the both sides and film on the damp can 
% also be included in the model.
%                                                        2018/11/27
% 6.May.2020: add the function of adding multi-ABHs.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% input parameters:
% xBound: length of the plate,
% yBound: wide of the plate
% R1: outer radius of the ABH
% R2: inner radius of the ABH
% Rd: radius of the damp
% gapX= xlocation of ABH
% gapY= ylocation of ABH
% lengElem1: inner element length 
% lengElem2; outer element length
% h: thickness of the plate
% dh: thickness of the damp
% numLayer: number of plates layer 
% dampLayer: number of damps layer 
% elemType: type of element
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
clear all;
close all;
clc;
tic;
%% Parameter
% size of plate
xBound=300;
yBound=200;

% size of ABH
R1=70;  % outer radius of ABH
R2=10;  % inner radius of ABH
Rd=35;  % radius of damper

gapX=71; % spacing between ABH in X
gapY=71; % spacing between ABH in Y
Xcenter = [gapX;-gapX;-gapX;gapX];
Ycenter = [gapY;gapY;-gapY;-gapY];

E_Modulu=200; % Elastic of damper

% size of elements and mesh
lengElem1=4;  % inner element size 
lengElem2=1;  % outer element size 
h=5.2;  % height/thickness of the plate
dh=2; % thickness of damper
numLayer=2; % layer of plate elements

% element options
elemType='Shell';   % the element type: ('SC'/'Solid'/'Shell')
if strcmp(elemType,'Shell')
    distributionType='Node';    % 'Node'/'Element'
end
dampSwitch =true;   % damp switch controling if damp is added, true or false 
if dampSwitch
    dampLayer=2;% layer of damper elements
    dampLoc = 'Down'; % 'Up' or 'Down'
    filmSwitch = false;  % true or false 
    if filmSwitch
        filmLayer = 1;  % definition of layers of film elements
        filmThickness = 0.1; % thickness of film
    end
end
tol = lengElem1/10.0;

if strcmp(elemType,'SC')
    numLayer=1;
    dampLayer=1;
    if filmSwitch
        filmLayer = 1;
    end
elseif strcmp(elemType,'Shell')
    numLayer=0;
    dampLayer=0;   
end

if strcmp(elemType,'Shell')&&dampSwitch
    if filmSwitch
       error('The function of adding film on the damp for shell element has not been developed in this version!') 
    end
end


%% Python Script Writer
%
% basic parameters 
% 
cd Model
fID=fopen('Parameter.py','w');
fprintf(fID,'%s\n',['xBound=',num2str(xBound)]);
fprintf(fID,'%s\n',['yBound=',num2str(yBound)]);
fprintf(fID,'%s\n',['R1=',num2str(R1)]);
fprintf(fID,'%s\n',['R2=',num2str(R2)]);
fprintf(fID,'%s\n',['Rd=',num2str(Rd)]);
fprintf(fID,'%s\n',['lengElem1=',num2str(lengElem1)]);
fprintf(fID,'%s\n',['lengElem2=',num2str(lengElem2)]);

fprintf(fID,'%s\n',['gapX=',num2str(gapX)]);
fprintf(fID,'%s\n',['gapY=',num2str(gapY)]);

fclose(fID);
system(['abaqus cae noGUI=DrawCircle.py']);

%% Read .inp
%
% open file
%
fID=fopen(['2Dmesh.inp'],'r');
StrInfor=textscan(fID, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter',',');
fclose(fID);
Str1Cell=StrInfor{1};
clear StrInfor
nodelineLoc=find(strcmp(Str1Cell,'*Node'));     % colume of the start of nodes
elemlineLoc=find(strcmp(Str1Cell,'*Element'));  % colume of the start of elements

%
% get node Information
%
fID=fopen(['2Dmesh.inp'],'r');
StrInfor=textscan(fID, '%f %f %f ','Delimiter',',','HeaderLines',nodelineLoc);
fclose(fID);
nodeInfor=[StrInfor{1},StrInfor{2},StrInfor{3}];
numNode=size(nodeInfor,1);  % number of total nodes
numNodeInput=numNode;       % initial number
clear StrInfor

%
% get element Information
%
fID=fopen(['2Dmesh.inp'],'r');
StrInfor=textscan(fID, '%f %f %f %f %f','Delimiter',',','HeaderLines',elemlineLoc);
fclose(fID);
elemInfor=[StrInfor{1},StrInfor{2},StrInfor{3},StrInfor{4},StrInfor{5}];
numElem=size(elemInfor,1);  % number of total elements
numElemInput=numElem;       % initial number
clear StrInfor

%% Modify mesh
%
% get the coordinate
%
TopZCoord=h*ones(numNode,1);
Dampelement = [];
for kABH = 1:4 
    xcentertemp = Xcenter(kABH);
    ycentertemp = Ycenter(kABH);
    
    nodeRadius=sqrt((nodeInfor(:,2)-xcentertemp).^2+(nodeInfor(:,3)-ycentertemp).^2);   % radical coordination of nodes
    elemCenLocation=1/4*[(nodeInfor(elemInfor(:,2),2)+nodeInfor(elemInfor(:,3),2)...
        +nodeInfor(elemInfor(:,4),2)+nodeInfor(elemInfor(:,5),2)), ...
        (nodeInfor(elemInfor(:,2),3)+nodeInfor(elemInfor(:,3),3)+...
        nodeInfor(elemInfor(:,4),3)+nodeInfor(elemInfor(:,5),3))];  % radical coordination of elements

    %
    % modify thickness of the plate element
    %
    elemRadius=sqrt((elemCenLocation(:,1)-xcentertemp).^2+(elemCenLocation(:,2)-ycentertemp).^2); 
    [BHnode,~]=find(nodeRadius<=R1);
    [Bottonnode,~]=find(nodeRadius<=R2);
    if dampSwitch
        [Dampelementtemp,~]=find(elemRadius<=Rd);
    end
    Dampelement = [Dampelement;Dampelementtemp];
    % >>>>>>>>> modify the thickness (according to the given function )>>>>>>>>
    
    p1=(h-0.2)/((R1-R2)^2);
    TopZCoord(BHnode)=p1*(nodeRadius(BHnode)-R2).^2+0.2;
    TopZCoord(Bottonnode)=0.2;
end
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if strcmp(elemType,'Shell')
    % shell element S4R
    if strcmp(distributionType,'Node')
        nodeInfor3D=nodeInfor;
        elemInfor3D=elemInfor;
    elseif strcmp(distributionType,'Element')
        nodeInfor3D=nodeInfor;
        elemInfor3D=elemInfor;
        TopZElem=h*ones(numElem,1);
        % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        p1=(h-0.2)/((R1-R2)^2);
        TopZElem(elemRadius<=R1)=p1*(elemRadius(elemRadius<=R1)-10).^2+0.2;
        TopZElem(elemRadius<=R2)=0.2;
        % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    end
    numelemtemp=size(elemInfor,1);
    if dampSwitch
        dampElemInfor=elemInfor(Dampelement,:);
        dampElemInfor(:,1)=numelemtemp+dampElemInfor(:,1);
    end
    
else
    error('the function of continous shell & solid element has been removed from this version.')
    
    % continous shell or solid element
    nodeInfor3D=zeros((numLayer+1)*numNode,4);
    nodeInfor3D(:,1)=transpose(1:(numLayer+1)*numNode);
    for klayer=0:numLayer
        nodeInfor3D((klayer*numNode+1):(klayer*numNode+numNode),2:end)=...
            [nodeInfor(:,2:3),TopZCoord*klayer/numLayer];
    end
    
    elemInfor3D=zeros(numLayer*numElem,9);
    elemInfor3D(:,1)=transpose(1:numLayer*numElem);
    for klayer=1:numLayer
        elemInfor3D(((klayer-1)*numElem+1):(klayer*numElem),2:end)=[elemInfor(:,2:end)+(klayer-1)*numNode,elemInfor(:,2:end)+klayer*numNode];
    end
    
    if dampSwitch
        % obtain nodes on the damp
        dampNode=unique([elemInfor(Dampelement,2);elemInfor(Dampelement,3);...
            elemInfor(Dampelement,4);elemInfor(Dampelement,5)]);
        if strcmp(dampLoc,'Up')
            % create new nodes for the damp
            for kdlayer=1:dampLayer
                nodeInfor3D=[nodeInfor3D;...
                    [dampNode+numLayer*numNode+kdlayer*numNode,...
                    nodeInfor(dampNode,2:end),...
                    TopZCoord(dampNode)+kdlayer/dampLayer*dh]];
            end
            % create damp element
            dampElemInfor=[];
            for kdlayer=1:dampLayer
                dampElemInfor=[dampElemInfor;...
                    [numLayer*numElem+(kdlayer-1)*numElem+Dampelement,...
                    elemInfor(Dampelement,2:end)+numLayer*numNode+(kdlayer-1)*numNode,...
                    elemInfor(Dampelement,2:end)+numLayer*numNode+kdlayer*numNode]];
            end
            
            if filmSwitch
                 filmElemInfor=[];
                 for kdlayer=1:filmLayer
                     nodeInfor3D=[nodeInfor3D;...
                         [dampNode+dampLayer*numNode+numLayer*numNode+kdlayer*numNode,...
                         nodeInfor(dampNode,2:end),...
                         TopZCoord(dampNode)+dh*ones(size(dampNode,1),1)+kdlayer/filmLayer*filmThickness]];
                 end
                 for kdlayer=1:filmLayer
                     filmElemInfor=[filmElemInfor;...
                         [numLayer*numElem+dampLayer*numElem+(kdlayer-1)*numElem+Dampelement,...
                         elemInfor(Dampelement,2:end)+numLayer*numNode+dampLayer*numNode+(kdlayer-1)*numNode,...
                         elemInfor(Dampelement,2:end)+numLayer*numNode+dampLayer*numNode+kdlayer*numNode]];
                 end
            end
            
        elseif strcmp(dampLoc,'Down')
            for kdlayer=1:dampLayer
                nodeInfor3D=[nodeInfor3D;...
                    [dampNode+numLayer*numNode+kdlayer*numNode,...
                    nodeInfor(dampNode,2:end),...
                    zeros(size(dampNode,1),1)-kdlayer/dampLayer*dh]];
            end
            % create damp element
            dampElemInfor=[numLayer*numElem+Dampelement,...
                elemInfor(Dampelement,2:end)+numLayer*numNode+numNode,...
                elemInfor(Dampelement,2:end)];
            
            for kdlayer=2:dampLayer
                dampElemInfor=[dampElemInfor;...
                    [numLayer*numElem+(kdlayer-1)*numElem+Dampelement,...
                    elemInfor(Dampelement,2:end)+numLayer*numNode+(kdlayer)*numNode,...
                    elemInfor(Dampelement,2:end)+numLayer*numNode+(kdlayer-1)*numNode]];
            end 
            if filmSwitch
                for kdlayer=1:filmLayer
                     nodeInfor3D=[nodeInfor3D;...
                         [dampNode+dampLayer*numNode+numLayer*numNode+kdlayer*numNode,...
                         nodeInfor(dampNode,2:end),...
                         -dh*ones(size(dampNode,1),1)-kdlayer/filmLayer*filmThickness]];
                end
                filmElemInfor = [];
%                 filmElemInfor=[numLayer*numElem+dampLayer*numElem+Dampelement,...
%                     elemInfor(Dampelement,2:end)+numLayer*numNode+numNode*dampLayer+numNode,...
%                     elemInfor(Dampelement,2:end)+numNode*dampLayer+numNode*dampLayer];
                for kdlayer=1:filmLayer
                    filmElemInfor=[filmElemInfor;...
                        [numLayer*numElem+dampLayer*numElem+(kdlayer-1)*numElem+Dampelement,...
                        elemInfor(Dampelement,2:end)+numLayer*numNode+dampLayer*numNode+kdlayer*numNode,...
                        elemInfor(Dampelement,2:end)+numLayer*numNode+dampLayer*numNode+(kdlayer-1)*numNode,...
                        ]];
                end
            end
        end
    end
end

%% New .inp generation
if dampSwitch
    if strcmp(elemType,'Shell')
        if strcmp(distributionType,'Node')
            fileName=['3Dmesh_',elemType,'_nodeDistribution.inp'];
        elseif strcmp(distributionType,'Element')
            fileName=['3Dmesh_',elemType,'_elementDistribution.inp'];
        end
    else
         fileName=['3Dmesh_',elemType,'.inp'];
    end
else
    if strcmp(elemType,'Shell')
        if strcmp(distributionType,'Node')
            fileName=['3Dmesh_',elemType,'_nodeDistribution_noDamp.inp'];
        elseif strcmp(distributionType,'Element')
            fileName=['3Dmesh_',elemType,'_elementDistribution_noDamp.inp'];
        end
    else
         fileName=['3Dmesh_',elemType,'_noDamp.inp'];
    end
end

%
% heading
%
fID=fopen(fileName,'w');
fprintf(fID,'%s\n',['*Heading']);
abc=datestr(now);
fprintf(fID,'%s\n',['ABAQUS job created by Bairyou Ben and on ' abc(1:11) ' at ' abc(13:end)] );

%
% part
%
fprintf(fID,'%s\n','**');
fprintf(fID,'%s\n','** Part');
fprintf(fID,'%s\n','**');
fprintf(fID,'%s\n','*Part, name=Part-1');

% nodes
fprintf(fID,'%s\n',['*Node']);
fclose(fID);
dlmwrite(fileName,nodeInfor3D,'-append','precision',10);

% element
fID=fopen(fileName,'a');
if strcmp(elemType,'SC')
    fprintf(fID,'%s\n',['*Element,type=SC8R,elset=Plate']);
elseif strcmp(elemType,'Solid')
    fprintf(fID,'%s\n',['*Element,type=C3D8,elset=Plate']);
elseif strcmp(elemType,'Shell')
    fprintf(fID,'%s\n',['*Element,type=S4,elset=Plate']);
end
fclose(fID);
dlmwrite(fileName,elemInfor3D,'-append','precision',10);
fID=fopen(fileName,'a');
if dampSwitch
    if strcmp(elemType,'SC')
        fprintf(fID,'%s\n',['*Element,type=SC8R,elset=Damp']);
    elseif strcmp(elemType,'Solid')
        fprintf(fID,'%s\n',['*Element,type=C3D8,elset=Damp']);
    elseif strcmp(elemType,'Shell')
        fprintf(fID,'%s\n',['*Element,type=S4,elset=Damp']);
    end
    fclose(fID);
    dlmwrite(fileName,dampElemInfor,'-append','precision',10);
    fID=fopen(fileName,'a');
    
    if filmSwitch
       if strcmp(elemType,'SC')
           fprintf(fID,'%s\n',['*Element,type=SC8R,elset=film']);
        elseif strcmp(elemType,'Solid')
            fprintf(fID,'%s\n',['*Element,type=C3D8,elset=film']);
       end
       fclose(fID);
       dlmwrite(fileName,filmElemInfor,'-append','precision',10);
       fID=fopen(fileName,'a');
    end
end

% boundary nodes set
LeftBoundary=find(abs(nodeInfor3D(:,2)+xBound)<tol);
RightBoundary=find(abs(nodeInfor3D(:,2)-xBound)<tol);
UpBoundary=find(abs(nodeInfor3D(:,3)+yBound)<tol);
DownBoundary=find(abs(nodeInfor3D(:,3)-yBound)<tol);
BoundaryNodes=[LeftBoundary;RightBoundary;UpBoundary;DownBoundary];
BoundaryNodes=unique(BoundaryNodes);
BoundaryNodes=nodeInfor3D(BoundaryNodes,1);
fprintf(fID,'%s\n','*Nset, nset=Boundary');
for bnnum=1:size(BoundaryNodes,1)
    if mod(bnnum,8)==0
        fprintf(fID,'%s\n',num2str(BoundaryNodes(bnnum,1)));
    else
        fprintf(fID,'%s,',num2str(BoundaryNodes(bnnum,1)));
    end
    if bnnum==size(BoundaryNodes,1) && mod(bnnum,8)~=0
        fprintf(fID,'\n');
    end
end

% section
if strcmp(elemType,'Shell') % shell
    if strcmp(distributionType,'Node')
        fprintf(fID,'%s\n',['*Shell Section,elset=Plate, material=Shell, offset=SNEG, NODAL THICKNESS']);
        fprintf(fID,'%s\n',[',']);
        fprintf(fID,'%s\n',['*NODAL THICKNESS']);
        for kn_number = 1:numNode
            fprintf(fID,'%s\n',[num2str(nodeInfor3D(kn_number,1)),',',num2str(TopZCoord(kn_number))]);
        end
    elseif strcmp(distributionType,'Element')
        fprintf(fID,'%s\n',['*Shell Section,elset=Plate, material=Shell, offset=SNEG, SHELL THICKNESS=Distribution-1']);
        fprintf(fID,'%s\n',[',']);
        fprintf(fID,'%s\n',['*DISTRIBUTION, LOCATION=ELEMENT, Name=Distribution-1,  Table=DISTRIBUTION-1_Table']);
        fprintf(fID,'%s\n',[', 1.0']);
        for kn_number = 1:numElem
            fprintf(fID,'%s\n',[num2str(elemInfor(kn_number,1)),',',num2str(TopZElem(kn_number))]);
        end
        
    end
    if dampSwitch
        fprintf(fID,'%s\n',['*Shell Section,elset=Damp, offset=SPOS, material=Rubber']);
        fprintf(fID,'%s\n',[num2str(dh)]);
    end
    
else % continous shell & solid
    if strcmp(elemType,'Solid')
        fprintf(fID,'%s\n',['*Solid Section,elset=Plate, material=Shell']);
        fprintf(fID,'%s\n',[',']);
    elseif strcmp(elemType,'SC')
        fprintf(fID,'%s\n',['*Shell Section,elset=Plate, material=Shell']);
        fprintf(fID,'%s\n',['1.0']);
    end
    if dampSwitch
        if strcmp(elemType,'Solid')
            fprintf(fID,'%s\n',['*Solid Section,elset=Damp, material=Rubber']);
            fprintf(fID,'%s\n',['']);
        elseif strcmp(elemType,'SC')
            fprintf(fID,'%s\n',['*Shell Section,elset=Damp, material=Rubber']);
            fprintf(fID,'%s\n',['1.0']);
        end
        if filmSwitch
            if strcmp(elemType,'Solid')
                fprintf(fID,'%s\n',['*Solid Section,elset=film, material=Alum']);
                fprintf(fID,'%s\n',['']);
            elseif strcmp(elemType,'SC')
                fprintf(fID,'%s\n',['*Shell Section,elset=film, material=Alum']);
                fprintf(fID,'%s\n',['1.0']);
            end
        end
    end
end
fprintf(fID,'%s\n','*End Part');
fprintf(fID,'%s\n','**');

%
% assembly
%
fprintf(fID,'%s\n','**');
fprintf(fID,'%s\n','** ASSEMBLY');
fprintf(fID,'%s\n','**');
fprintf(fID,'%s\n','*Assembly, name=Assembly');
fprintf(fID,'%s\n','** ');
fprintf(fID,'%s\n','*Instance, name=Part-1-1, part=Part-1');
fprintf(fID,'%s\n','*End Instance');
fprintf(fID,'%s\n','**');
fprintf(fID,'%s\n','*End Assembly');

if strcmp(elemType,'Shell')&&strcmp(distributionType,'Element')
    fprintf(fID,'%s\n',['*Distribution Table, name=DISTRIBUTION-1_Table']);
    fprintf(fID,'%s\n',['length']);
end

%
% boundary
%
fprintf(fID,'%s\n','**');
fprintf(fID,'%s\n','** BOUNDARY CONDITIONS');
fprintf(fID,'%s\n','**');
fprintf(fID,'%s\n','*Boundary');
if strcmp(elemType,'Shell') % shell
    fprintf(fID,'%s\n','Part-1-1.Boundary,1,6,0');
else
    fprintf(fID,'%s\n','Part-1-1.Boundary,1,3,0');
end
%
% material
%
fprintf(fID,'%s\n',['*Material, name=Shell']);
fprintf(fID,'%s\n',['*Elastic']);
fprintf(fID,'%s\n',['71000,0.33']);
fprintf(fID,'%s\n',['*Density']);
fprintf(fID,'%s\n',['0.00000000282']);
fprintf(fID,'%s\n',['*Damping, structural=0.002']);

fprintf(fID,'%s\n',['*Material, name=Rubber']);
fprintf(fID,'%s\n',['*Elastic']);
fprintf(fID,'%s\n',[num2str(E_Modulu),',0.45']);
fprintf(fID,'%s\n',['*Density']);
fprintf(fID,'%s\n',['0.00000000185']);
fprintf(fID,'%s\n',['*Damping, structural=0.3']);

fprintf(fID,'%s\n',['*Material, name=Alum']);
fprintf(fID,'%s\n',['*Elastic']);
fprintf(fID,'%s\n',['71000,0.33']);
fprintf(fID,'%s\n',['*Density']);
fprintf(fID,'%s\n',['0.00000000282']);
fprintf(fID,'%s\n',['*Damping, structural=0.002']);
%
% step
%
fprintf(fID,'%s\n',['** ----------------------------------------------------------------']);
fprintf(fID,'%s\n',['** ']);
fprintf(fID,'%s\n',['** STEP: Step-1']);
fprintf(fID,'%s\n',['** ']);
fprintf(fID,'%s\n',['*Step, name=Step-1, nlgeom=NO, perturbation']);
% fprintf(fID,'%s\n',['*Frequency, eigensolver=LANCZOS, normalization=displacement']);
fprintf(fID,'%s\n',['*Frequency, eigensolver=Lanczos, sim, acoustic coupling=projection, normalization=mass']);
fprintf(fID,'%s\n',[', , 8950., , , ']);
fprintf(fID,'%s\n',['** ']);
fprintf(fID,'%s\n',['** OUTPUT REQUESTS']);
fprintf(fID,'%s\n',['** ']);
fprintf(fID,'%s\n',['** ']);
fprintf(fID,'%s\n',['**  FIELD OUTPUT: F-Output-1']);
fprintf(fID,'%s\n',['** ']);
fprintf(fID,'%s\n',['*Output, field, variable=PRESELECT']);
fprintf(fID,'%s\n',['*End Step']);

fprintf(fID,'%s\n',['** ----------------------------------------------------------------']);
fprintf(fID,'%s\n',['** ']);
fprintf(fID,'%s\n',['** STEP: Step-2']);
fprintf(fID,'%s\n',['** ']);
fprintf(fID,'%s\n',['*Step, name=Step-2, nlgeom=NO, perturbation, unsymm=YES']);
fprintf(fID,'%s\n',['*Complex Frequency, friction damping=NO']);
fprintf(fID,'%s\n',[', , 8950., ']);
fprintf(fID,'%s\n',['** ']);
fprintf(fID,'%s\n',['** OUTPUT REQUESTS']);
fprintf(fID,'%s\n',['** ']);
fprintf(fID,'%s\n',['** ']);
fprintf(fID,'%s\n',['**  FIELD OUTPUT: F-Output-2']);
fprintf(fID,'%s\n',['** ']);
fprintf(fID,'%s\n',['*Output, field, variable=PRESELECT']);
fprintf(fID,'%s\n',['*End Step']);

fclose(fID);
% %system(['abaqus cae noGUI=Mesh_SolidC3D20.py']);%SolidC3D20 by WXD added
% system(['echo y|cmd/c abaqus job=',fileName(1:end-4),' cpus=2 int']);

cd ..

% 
% %
% % reuslt reader
% %u
% system(['abaqus viewer noGUI=modalreader.py']);
% modalInfor=load('dampInfor.txt');
% % plot(modalInfor(:,1),modalInfor(:,2));
% sum(modalInfor(:,2).^2)
toc;