%% Code to solve 1D trusses systems with stiffness method
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
%indexes to define the same value n times in a vector (e.g: E(1:numelements) = X)
%If some or every truss has different values for E, write each one in a
%vector (e.g: (E = [E1 E2 E3 ...])

%area = Area of the cross section of each truss. Just like with E, if all
%the trusses have the same area, use indexes (e.g: area(1:numelements) = X)
%If some or every truss has different values for area, write each one in a
%vector (e.g: (area = [area1 area2 area3 ...])

%nodes = The X coordinates of each node in the units that are defined in
%your problem. I recommend placing the origin on the first node.

%NodesUnion = The indexes of the nodes that make each truss. Every truss is
%between two nodes.

%Displacements = Boundary conditions for each node (0 if the node is fixed
%or 1 if it can move)

%Forces = Force vector of the system (0 if there isn't a force acting on a
%certain node or an X value if there is indeed a force. Depending on the
%direction of each force, you may have to use a positive or negative sign.

%L = Length of each truss.

clc; clear ;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%

%Solve a 1D truss system with the stiffness method. You can use any of the
%three examples down below in this section, all taken straight out of the
%book I mention up here, or you can modify these examples a little to solve
%your own system, you just need to follow the same structure.

%NOTE: This code ONLY works if ALL the nodal displacements are unknown, the
%algorithm won't solve correctly systems that have one or more known nodal
%displacements.


%Example 3.1 of the book 

numelements=3;                     %Number of elements
E=[30*10^6 30*10^6 15*10^6];        %Elastic modulus of each truss
area=[1 1 2];                       %Cross section areas of each truss
nodes=[0 ; 30 ; 60 ; 90];           %Coordinates of each node
NodesUnion=[1 2 ; 2 3 ; 3 4];       %Nodes that make each truss
Displacements=[0 1 1 0];          %Nodal displacements
Forces=[0 3000 0 0];               %Forces acting on each node
L=[30 30 30];                       %Length of each truss

%Problem 3.1 of the book

% numelements=3;                      %Number of elements
% E(1:numelements)=10*10^6;        %Elastic modulus of each truss
% area(1:numelements)=1;              %Cross section areas of each truss
% nodes=[0 ; 10 ; 20 ; 30];           %Coordinates of each node
% NodesUnion=[1 2 ; 2 3 ; 3 4];       %Nodes that make each truss
% Displacements=[0 1 1 0];          %Nodal displacements
% Forces=[0 0 1000 0];               %Forces acting on each node
% L=[10 10 10];                       %Length of each truss

%Problem 3.6 of the book

% numelements=3;                     %Number of elements
% E=[30*10^6 10*10^6 10*10^6];        %Elastic modulus of each truss
% area(1:numelements)=2;               %Cross section areas of each truss
% nodes=[0 ; 50 ; 80 ; 80];           %Coordenadas de cada nodo
% NodesUnion=[1 2 ; 2 3 ; 2 4];       %Nodes that make each truss
% Displacements=[0 1 0 0];          %Nodal displacements
% Forces=[0 8000 0 0];               %Forces acting on each node
% L=[50 30 30];                      %Length of each truss

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Algorithm~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Construction of every stiffness matrix of each truss
A = zeros(2,2,numelements);
for i = 1:numelements
    A(:,:,i) = ((E(i).*area(i))./L(i))*[1 -1 ; -1 1];   
end
A

%Assembly of each truss stiffness matrix inside the global stiffness matrix
S=size(nodes,1); 
MG=zeros(S,S);                 %Initiliaze the global stiffnes matrix empty                                    
for i=1:numelements
    m=NodesUnion(i,1);
    n=NodesUnion(i,2);
    MG([m n],[m n],i)=A(:,:,i);      %Put each stiffness matrix inside a 3D array 
end
MG=sum(MG,3)                        %Overlap the 3D array and sum the three matrices

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
d = MGR\Forces'         
k=1;
dfinal=zeros(S,1);
for i=1:length(Displacements) %Insert the new found displacements inside the original displacements vector
    if Displacements(i)==0
        dfinal(i,1)=0;
    else
        dfinal(i,1)=d(k);
        k=k+1;
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Results~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

Stresses=zeros(1,numelements);
for i=1:numelements     %Calculate the stresses of each truss
    indice=NodesUnion(i,:);
    Stresses(i)=(E(i)./L(i))*[-1 1]*[dfinal(indice(1,1));dfinal(indice(1,2))];
end

dfinal                              %Print results
Stresses
Reactions=MG*dfinal               %Calculate the reactions of all trusses.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Graphics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

line([nodes(1) nodes(end)],[0 0],'Color','k','LineWidth',2);
hold on
plot(nodes',zeros(1,length(nodes)),'or','MarkerSize',6,'MarkerFaceColor','r');
grid on
NodosDesp=nodes'+dfinal';
Nocero=find(dfinal~=0);
hold on
plot(NodosDesp(Nocero),zeros(1,length(Nocero)),'og','MarkerSize',6,'MarkerFaceColor','g');
axis([nodes(1)-10 nodes(end)+10 -1 1]);
legend({'Truss';'Fixed node' ; 'Free node'});