%% Code to solve 2D trusses systems with stiffness method
%Author: Rolando Valdez Guzmán
%Aka: Tutoingeniero
%Youtube channel: https://www.youtube.com/channel/UCU1pdvVscOdtLpRQBp-TbWg
%Versión: 1.0
%Updated: 22/jul/2020

%References: "A First Course in the Finite Element Method" by Daryl. L.
%Logan

%% Variables:
%numelements = A scalar that defines the number of trusses in the system.

%E =  Elastici modulus of each truss. If every truss has the same E, use
%indices to define the same value n times in a vector (e.g: E(1:numelements) = X)
%If some or every truss has different values for E, write each one in a
%vector (e.g: (E = [E1 E2 E3 ...])

%area = Area of the cross section of each truss. Just like with E, if all
%the trusses have the same area, use indices (e.g: area(1:numelements) = X)
%If some or every truss has different values for area, write each one in a
%vector (e.g: (area = [area1 area2 area3 ...])

%nodes = The [X,Y] coordinates of each node in the units that are defined in
%your problem. I recommend placing the origin on the first node.
%e.g [x1 y1, x2 y2, ...]

%NodesUnion = The indices of the nodes that make each truss. Every truss is
%between two nodes.

%Displacements = Boundary conditions for each XY componentof each node 
%(0 if the node component is fixed or 1 if it can move). e.g: [dx1 dy1 dx2 dy2 ...]

%Forces = Force vector of the system (0 if there isn't a force acting on a
%certain node or an X value if there is indeed a force. Depending on the
%direction of each force, you may have to use a positive or negative sign.
%e.g : [Fx1 Fy1 Fx2 Fy2 ...]

%L = Length of each truss.

clc; clear ;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%

%Solve a 2D truss system with the stiffness method. You can use any of the
%three examples down below in this section, all taken straight out of the
%book I mention up here, or you can modify these examples a little to solve
%your own system, you just need to follow the same structure.

%NOTE: This code ONLY works if ALL the nodal displacements are unknown, the
%algorithm won't solve correctly systems that have one or more known nodal
%displacements.

%Ejemplo 3.5 del libro

numelements=3;                            
E(1:numelements)=30*10^6;                  
area(1:numelements)=2;                     
nodes=[0 0 ; 0 120; 120 120; 120 0];        
NodesUnion=[1 2;1 3;1 4];                   
Displacements=[1 1 0 0 0 0 0 0];         
Forces=[0 -10000 0 0 0 0 0 0 ];           

%Problem 3.22

% numelements=3;
% E(1:numelements)=10*10^6;
% area(1:numelements)=1;
% nodes=[0 0 ; -100 100*tand(60); -100  0 ;-100 -100*tand(30)];
% NodesUnion=[1 2;1 3;1 4];
% Displacements=[1 1 0 0 0 0 0 0 ];
% Forces=[1000 1000 0 0 0 0 0 0];

%Problem 3.29

% numelements=3;
% E(1:numelements)=210*10^9;
% area(1:numelements)=0.0004;
% nodes=[0 0 ; 0 2 ; -3 0 ; -5*sind(30) -5*sind(60)];
% NodesUnion=[1 2 ; 1 3 ; 1 4];
% Displacements=[1 1 0 0 0 0 0 0];
% Forces=[0 -40000 0 0 0 0 0 0];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Algorithm~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Construction of every stiffness matrix of each truss

L=zeros(1,numelements);              %Initialize vectors to be filled
angle=zeros(1,numelements);         
for i=1:numelements
    index=NodesUnion(i,:);
    b=nodes(index(2),:);                               %Pick Y2 y Y1
    a=nodes(index(1),:);                               %Pick X2 y X1
    L(i)=norm(b-a); %Calculate length using distance between two points
    
    %Determine the angle of the truss
    dx=nodes(index(2),1) - nodes(index(1),1);
    dy=nodes(index(2),2) - nodes(index(1),2);   
    %Determine the quadrant in which the truss is
    angle(i)=atand(dy/dx);                                 %Quadrant 1
    if sign(dy)==1 && sign(dx)==-1                          %Quadrant 2
        angle(i)=180+angle(i);
    elseif dy==0 && sign(dx)==-1 %Over the horizontal but on quadrant 2
        angle(i)=180;
    elseif sign(dy)==-1 && sign(dx)==-1                     %Quadrant 3
        angle(i)=180+angle(i);
    elseif sign(dy)==-1 && sign(dx)==1                      %Quadrant 4
        angle(i)=360+angle(i);
    end  
end
L
angle
k=(E.*area)./L                           %AE/L constant for each truss

%Assembly of each truss stiffness matrix inside the global stiffness matrix

A=zeros(4,4,numelements);     %Initialize the stiffness matrices
for i=1:numelements
    %Calculate the stiffness matrix of each truss first
    A(:,:,i)=k(i)*transformadalineal(angle(i));  %TRANSFORMADALINEAL is down below 
    
    %Make A array to cells but the 4x4 matrices of A will be divided in 4
    %submatrices of 2x2
    j=NodesUnion(i,:);                      
    B(:,:,i)=mat2cell(A(:,:,i),[2 2],[2 2]);
    
    %Assign each submatrix in the correct indices of the global stiffness
    %matrix
    C(j(1),j(1),i)=B(1,1,i);
    C(j(1),j(2),i)=B(1,2,i);
    C(j(2),j(1),i)=B(2,1,i);
    C(j(2),j(2),i)=B(2,2,i);
end
A

S=2*size(nodes,1);                         %Initiliaze the global matrix
m=cell(S/2,S/2);
for i=1:size(nodes,1)
    for j=1:size(nodes,1)
        %On each lap we pick up the elements with the same index (i,j),
        %then overlap them and sum them.
        clear x
        x(:,:,:)=cell2mat(reshape(C(i,j,:),1,[],numelements));
        m(i,j)={sum(x,3)};
        
        %If this (i,j) index has no submatrix assigned, assing a 2x2
        %submatrix of zeros
        if size(m{i,j})==[0 0]  
            m(i,j)={zeros(2,2)};
        end
    end
end
MG=cell2mat(m)          %Make the global matrix a numeric array

%Reduce the global matrix
v=find(Displacements==0);                   
MGR=MG;
MGR(v,:)=0;     %Find the rows and columns that have full zeros
MGR(:,v)=0;
indicefil=zeros(1,S);
indicecol=zeros(1,S);
for i=1:S
    if MGR(i,:)==0
        indicefil(i)=i;
    end
    if MGR(:,i)==0
        indicecol(i)=i;
    end
end
MGR(indicefil~=0,:)=[];    %Erase the rows and columns full of zeros
MGR(:,indicecol~=0)=[]
Forces(indicefil~=0)=[]  %Erase the rows of the force vector too.

%Solve the system
d=MGR\Forces';         
dfinal=zeros(S,1);
k=1;
for i=1:length(Displacements) %Insert the new found displacements inside the original displacements vector
    if Displacements(i)==0
        dfinal(i,1)=0;
    else
        dfinal(i,1)=d(k);
        k=k+1;
    end
end

%Coefficient vector [-C -S C S] to calculate the stresses
Ve={@(x) -cosd(x) @(x) -sind(x) @(x) cosd(x) @(x) sind(x)};

stressesC=zeros(numelements,4);
for i=1:numelements
    for j=1:4
        stressesC(i,j)=feval(Ve{1,j},angle(i)); %Evaluate Ve for each angle of each truss
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Results~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

d2=mat2cell(dfinal,2*ones(1,size(nodes,1)),1); %Divide dfinal in 2x1 submatrices
Stresses=zeros(1,numelements);
Flocal=zeros(numelements,4);
j=1;
for i=1:numelements     %Calculate the stresses of each truss
    index=NodesUnion(i,:);
    Stresses(i)=(E(i)./L(i))*stressesC(i,:)*[d2{index(1,1)};d2{index(1,2)}];
    Flocal(i,:)=A(:,:,i)*[d2{index(1,1)};d2{index(1,2)}]; 
    j=j+2;
end

dfinal                              %Print results
Stresses
Reacciones=MG*dfinal               %Reactions of each node
Flocal                             %Local reactions

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Graph~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

d3=reshape(dfinal,[2,size(nodes,1)])';
NodosDesp=nodes + d3;

for j=1:numelements
    index=NodesUnion(j,:);
    line([nodes(index(1),1) nodes(index(2),1)],...
        [nodes(index(1),2) nodes(index(2),2)],...
        'LineWidth',1.5,'Color','k');
    hold on
    line([NodosDesp(index(1),1) NodosDesp(index(2),1)],...
        [NodosDesp(index(1),2) NodosDesp(index(2),2)],...
        'LineWidth',1,'Color','b');
end

for i=1:size(nodes,1)
    plot(nodes(i,1),nodes(i,2),'ro','MarkerSize',6,'MarkerFaceColor','r');
    hold on
end
axis([min(nodes(:,1))-10 max(nodes(:,2))+10 min(nodes(:,1))-10 max(nodes(:,2)+10)])
grid on
zoom on
NodosDespx=NodosDesp(:,1);
NodosDespy=NodosDesp(:,2);
Nocerox=find(d3(:,1)~=0);
Noceroy=find(d3(:,2)~=0);
plot(NodosDespx(Nocerox),NodosDespy(Noceroy),'og','MarkerSize',6,'MarkerFaceColor','g');

%~~~~~~~~~~~~~~~~~~~~~~~~~transformadalineal function~~~~~~~~~~~~~~~~~~~~~~%

function [M]=transformadalineal(theta)
m={@(x) (cosd(x))^2 @(x) sind(x)*cosd(x) @(x) -(cosd(x))^2 @(x) -sind(x)*cosd(x);...
    @(x) sind(x)*cosd(x) @(x) (sind(x))^2 @(x) -sind(x)*cosd(x) @(x) -(sind(x))^2;...
    @(x) -(cosd(x))^2  @(x) -sind(x)*cosd(x) @(x) (cosd(x))^2 @(x) sind(x)*cosd(x);...
    @(x) -sind(x)*cosd(x) @(x) -(sind(x))^2 @(x) sind(x)*cosd(x) @(x) (sind(x))^2};

M=zeros(4,4);
for i=1:4
    for j=1:4
        M(i,j)=feval(m{i,j},theta);
    end
end
end