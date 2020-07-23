%% Code to solve 3D trusses systems with stiffness method
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

%nodes = The [X,Y,Z] coordinates of each node in the units that are defined in
%your problem. I recommend placing the origin on the first node.
%e.g [x1 y1 z1, x2 y2 z2, ...]

%NodesUnion = The indices of the nodes that make each truss. Every truss is
%between two nodes.

%Displacements = Boundary conditions for each XYZ componentof each node 
%(0 if the node component is fixed or 1 if it can move). 
%e.g: [dx1 dy1 dz1 dx2 dy2 dz2...]

%Forces = Force vector of the system (0 if there isn't a force acting on a
%certain node or an X value if there is indeed a force. Depending on the
%direction of each force, you may have to use a positive or negative sign.
%e.g : [Fx1 Fy1 Fz1 Fx2 Fy2 Fz2 ...]

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


%Example 3.8 of the book

numelements=3;
E(1:numelements)=1.2*10^6;
area=[0.302 0.729 0.187];
nodes=[72 0 0 ; 0 36 0 ; 0 36 72 ; 0 0 -48];
NodesUnion=[1 2;1 3;1 4];
Displacements=[1 0 1 0 0 0 0 0 0 0 0 0];
Forces=[0 0 -1000 0 0 0 0 0 0 0 0 0];

%Example 3.9 of the book

% numelements=3;
% E(1:numelements)=210*10^9;
% area(1:numelements)=10*10^-4;
% nodes=[12 -3 -4 ; 0 0 0 ; 12 -3 -7 ; 14 6 0];
% NodesUnion=[1 2 ; 1 3 ; 1 4];
% Displacements=[1 1 1 0 0 0 0 0 0 0 0 0];
% Forces=[20000 0 0 0 0 0 0 0 0 0 0 0];

%Problem 3.40 of the book

% numelements=4;
% E(1:numelements)=210*10^9;
% area(1:numelements)=10*10^-4;
% nodes=[4 4 3 ; 0 4 0 ; 0 4 6 ; 4 0 3 ; 8 -1 1];
% NodesUnion=[1 2 ;1 3 ; 1 4 ; 1 5];
% Displacements=[1 1 1 0 0 0 0 0 0 0 0 0 0 0 0];
% Forces=[0 -10000 0 0 0 0 0 0 0 0 0 0 0 0 0];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Algorithm~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Construction of every stiffness matrix of each truss

L = zeros(1,numelements);   %Initialize vectors to be filled
Cx = zeros(1,numelements);
Cy = zeros(1,numelements);
Cz = zeros(1,numelements);
LAMDA = zeros(6,6);

for i = 1:numelements
    index = NodesUnion(i,:);                                           %Define point 1 and 2 by indices
    P1 = [nodes(index(1),1) nodes(index(1),2) nodes(index(1),3)];    %XYZ coordinates of point 1
    P2 = [nodes(index(2),1) nodes(index(2),2) nodes(index(2),3)];    %XYZ coordinates of point 2
    L(i) = norm(P1-P2);                             %Calculate length using distance between two points
    
    %Calculate direction cosines 
    Cx(i) = (P2(1) - P1(1))/ L(i);
    Cy(i) = (P2(2) - P1(2))/ L(i);
    Cz(i) = (P2(3) - P1(3))/ L(i);
    
    %Calculate lambda for each truss and assembly the local stiffness
    %matrix
    lamda = [Cx(i)^2 Cx(i)*Cy(i) Cx(i)*Cz(i) ; Cy(i)*Cx(i) Cy(i)^2 Cy(i)*Cz(i) ;...
           Cz(i)*Cx(i) Cz(i)*Cy(i) Cz(i)^2];
    LAMDA(:,:,i) = [lamda -lamda ; -lamda lamda];
end
k = (E.*area)./L;                       %AE/L constant of each truss
A = zeros(6,6);     

%Assembly of each truss stiffness matrix inside the global stiffness matrix

for i = 1:numelements
    %Calculate the stiffness matrix of each truss first
    A(:,:,i) = k(i)*LAMDA(:,:,i);
    
    %Make A array to cells but the 6x6 matrices of A will be divided in 4
    %submatrices of 3x3 
    j = NodesUnion(i,:);                    
    B(:,:,i) = mat2cell(A(:,:,i),[3 3],[3 3]);
    
    %Assign each submatrix in the correct indices of the global stiffness
    %matrix
    C(j(1),j(1),i) = B(1,1,i);
    C(j(1),j(2),i) = B(1,2,i);
    C(j(2),j(1),i) = B(2,1,i);
    C(j(2),j(2),i) = B(2,2,i);
end
A

S = 3*size(nodes,1);                          %Initiliaze the global matrix
m = cell(S/3,S/3);
for i = 1:size(nodes,1)
    for j = 1:size(nodes,1)
        %On each lap we pick up the elements with the same index (i,j),
        %then overlap them and sum them.
        clear x
        x(:,:,:) = cell2mat(reshape(C(i,j,:),1,[],numelements));
        m(i,j) = {sum(x,3)};
        
        %If this (i,j) index has no submatrix assigned, assing a 2x2
        %submatrix of zeros
        if size(m{i,j}) == [0 0]
            m(i,j) = {zeros(3,3)};
        end
        
    end
end
MG = cell2mat(m)          %Make the global matrix a numeric array

%Reduce the global matrix
v = find(Displacements==0);                   
MGR = MG;
MGR(v,:) = 0;     %Find the rows and columns that have full zeros
MGR(:,v) = 0;
indicefil = zeros(1,S);
indicecol = zeros(1,S);
for i = 1:S
    if MGR(i,:) == 0
        indicefil(i) = i;
    end
    if MGR(:,i) == 0
        indicecol(i) = i;
    end
end
MGR(indicefil~=0,:) = [];    %Erase the rows and columns full of zeros
MGR(:,indicecol~=0) = []
Forces(indicefil~=0) = [];  %Erase the rows of the force vector too.


%Solve the system
d = MGR\Forces';         
dfinal = zeros(S,1);
k = 1;
for i = 1:length(Displacements) %Insert the new found displacements inside the original displacements vector
    if Displacements(i) == 0
        dfinal(i,1) = 0;
    else
        dfinal(i,1) = d(k);
        k = k+1;
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Results~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

d2 = mat2cell(dfinal,3*ones(1,size(nodes,1)),1); %Divide dfinal in 3x1 submatrices
Stresses = zeros(1,numelements);
Flocal = zeros(numelements,6);
j = 1;
for i = 1:numelements  %Calculate the stresses of each truss
    index = NodesUnion(i,:);
    Stresses(i) = (E(i)./L(i)) * [-Cx(i) -Cy(i) -Cz(i) Cx(i) Cy(i) Cz(i)] * [d2{index(1,1)} ; d2{index(1,2)}];
    Flocal(i,:) = A(:,:,i)*[d2{index(1,1)} ; d2{index(1,2)}];
    j = j + 2;
end

dfinal                              %Print results
Stresses
Reacciones=MG*dfinal               %Reactions of each node
Flocal                             %Local reactions

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Graph~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

d3 = reshape(dfinal,[3,size(nodes,1)])';
NodosDesp = nodes + d3;
d4 = reshape(Displacements,[3,size(nodes,1)])';

for j = 1:numelements
    index = NodesUnion(j,:);
    line([nodes(index(1),1) nodes(index(2),1)],...
        [nodes(index(1),2) nodes(index(2),2)],...
        [nodes(index(1),3) nodes(index(2),3)],...
        'LineWidth',1.5,'Color','k');
    hold on
    line([NodosDesp(index(1),1) NodosDesp(index(2),1)],...
        [NodosDesp(index(1),2) NodosDesp(index(2),2)],...
        [NodosDesp(index(1),3) NodosDesp(index(2),3)],...
        'LineWidth',1,'Color','b');
end

for i = 1:size(nodes,1)
    plot3(nodes(i,1),nodes(i,2),nodes(i,3),'ro','MarkerSize',6,'MarkerFaceColor','r');
    hold on
end
grid on
NodosDespx = NodosDesp(:,1);
NodosDespy = NodosDesp(:,2);
NodosDespz = NodosDesp(:,3);
[fil,col] = find(d4 ~= 0);
plot3(NodosDespx(fil),NodosDespy(fil),NodosDespz(fil),'og','MarkerSize',6,'MarkerFaceColor','g');
view(45,45);