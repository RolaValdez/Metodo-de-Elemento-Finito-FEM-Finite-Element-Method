%% Código para resolver sistemas de elementos barra en 3D
%Autor: Rolando Valdez Guzmán
%Alias: Tutoingeniero
%Canal de Youtube: https://www.youtube.com/channel/UCU1pdvVscOdtLpRQBp-TbWg
%Versión: 2.0
%Actualizado: 15/sep/2022

%Referencias: "A First Course in the Finite Element Method" por Daryl. L.
%Logan

%% ~~~~~~~~~~~~~~INSTRUCCIONES DE USO! LEER DETALLADAMENTE~~~~~~~~~~~~~~~~
% Variables:
%numelementos = Es un escalar que define el número de barras que tiene el
%sistema.

%E =  Módulo de elasticidad de las barras. si todos tienen el mismo usa
%índices en E para definir el mismo valor n veces en un vector (por ejemplo: E(1:numelementos) = X)
%Si algunas o todas la barras tienen diferentes valores para E, escribe
%cada uno en un vector (E = [E1 E2 E3 ...])

%area = Area de la sección transversal de cada barra. Al igual que con E,
%si todas las barras tienen la misma área usa índices (por ejemplo: area(1:numelementos) = X)
%Si algunas o todas la barras tienen diferentes valores de area, escribe
%cada uno en un vector (area = [area1 area2 area3 ...])

%nodos = Son las coordenadas [X,Y,Z] de cada nodo en las unidades que defina tu
%problema. Recomiendo siempre poner el origen en el primer nodo. Ejemplo:
%[x1 y1 z1; x2 y2 z1; ...]

%UnionNodos = Son los índices de los nodos que conforman a cada barra. Cada
%barra se conforma por una línea que va de un nodo a otro.

%Desplazamientos = Condiciones de frontera para cada componente XY de cada nodo 
%(0 si el nodo está empotrado y 1 si el nodo puede moverse). Por ejemplo:
%[dx1 dy1 dz1 dx2 dy2 dz2 dx3 dy3 dz3...]

%Fuerzas = Vector de fuerzas del sistema (0 si no hay una fuerza actuando
%en dicho nodo y un valor cualquiera si hay una fuerza. Dependiendo del
%sentido de cada fuerza, se deberá usar un signo negativo o positivo).
%[Fx1 Fy1 Fz1 Fx2 Fy2 Fz2 Fx3 Fy3 Fz3...]

clc; clear;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%

%Resuelve un sistema de barras en 3D con el método de elemento finito.
%Puedes usar cualquiera de los tres ejemplos que dejo aquí en esta sección,
%todos del libro que dejo arriba o basarte en los ejemplos para resolver tu
%propio sistema, sólo debes seguir la misma estructura.

%NOTA: Este código sólo funciona si se desconocen TODOS los desplazamientos
%nodales, el algoritmo NO resolverá sistemas con desplazamientos conocidos
%correctamente.

%Ejemplo 3.8 del libro

% numelementos=3;
% E(1:numelementos)=1.2*10^6;
% area=[0.302 0.729 0.187];
% nodos=[72 0 0 ; 0 36 0 ; 0 36 72 ; 0 0 -48];
% UnionNodos=[1 2;1 3;1 4];
% Desplazamientos=[1 0 1 0 0 0 0 0 0 0 0 0];
% Fuerzas=[0 0 -1000 0 0 0 0 0 0 0 0 0];

%Ejemplo 3.9 del libro

% numelementos=10;
% E(1:numelementos)=158000;
% area=[20 20 80 20 80 60 80 20 80 20 ];
% nodos=[250,0,0 ; 250,400,0 ; 0 0 0 ; 0 400 0; 0 0 250; 0 400 250; 0 0 600; 0 400 600];
% UnionNodos=[1 5;2 6;3 5; 3 6; 4 6 ; 5 6 ; 5 7 ; 6 7 ; 6 8 ; 7 8 ];
% Desplazamientos=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 1 1];
% Fuerzas=[0 0 0    0 0 0    0 0 0    0 0 0    0 0 0   10000 5000 0   4000 -5000 -10000      0 -10000 -10000 ];

%Problema 3.40 del libro

% numelementos=4;
% E(1:numelementos)=210*10^9;
% area(1:numelementos)=10*10^-4;
% nodos=[4 4 3 ; 0 4 0 ; 0 4 6 ; 4 0 3 ; 8 -1 1];
% UnionNodos=[1 2 ;1 3 ; 1 4 ; 1 5];
% Desplazamientos=[1 1 1 0 0 0 0 0 0 0 0 0 0 0 0];
% Fuerzas=[0 -10000 0 0 0 0 0 0 0 0 0 0 0 0 0];

% numelementos = 18;
% E(1:numelementos) = 210*10^9;
% area(1:numelementos)= 0.1;
% nodos=[-3 2 0; 3 2 0; 3 -2 0; -3 -2 0; -1.5 1 7.2; 1.5 1 7.2; 1.5 -1 7.2; -1.5 -1 7.2];
% UnionNodos = [1 2 ; 2 3 ; 3 4; 1 4; 1 3; 5 6; 6 7; 7 8; 5 8; 6 8; 1 5; 2 6; 3 7;...
%               4 8; 2 5; 3 6; 4 7; 1 8];
% Desplazamientos=[0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1];
% Fuerzas=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 -80000 60000 0 0 60000 0 0 0 0 -80000];

numelementos = 25;
E(1:numelementos) = 10*10^7;
area(1:numelementos) = 4;
nodos=[-37.5 0 200; 37.5 0 200; -37.5 37.5 100; 37.5 37.5 100; 37.5 -37.5 100; ...
       -37.5 -37.5 100; -100 100 0; 100 100 0; 100 -100 0; -100 -100 0];
UnionNodos = [1 2 ; 3 6; 5 6; 4 5; 3 4; 3 7; 6 10; 5 9; 4 8; 1 4; 1 5; 2 6; 2 3;...
              10 3; 6 7; 9 4; 5 8; 7 4; 3 8; 10 5; 6 9; 2 5; 2 4; 1 3; 1 6];
Desplazamientos=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0];
Fuerzas=[1000 10000 -5000 0 10000 -5000 500 0 0 0 0 0 0 0 0 500 0 0 0 0 0, ...
         0 0 0 0 0 0 0 0 0];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Algoritmo~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Construcción de las matrices de rigidez de cada barra.
L = zeros(1,numelementos);
Cx = zeros(1,numelementos); Cy = zeros(1,numelementos);
Cz = zeros(1,numelementos); LAMDA = zeros(6,6);

for i = 1:numelementos
    indice = UnionNodos(i,:);                                   
    P1 = [nodos(indice(1),1) nodos(indice(1),2) nodos(indice(1),3)];
    P2 = [nodos(indice(2),1) nodos(indice(2),2) nodos(indice(2),3)];
    L(i) = norm(P1-P2);
    Cx(i) = (P2(1) - P1(1))/ L(i);
    Cy(i) = (P2(2) - P1(2))/ L(i);
    Cz(i) = (P2(3) - P1(3))/ L(i);
    lamda = [Cx(i)^2 Cx(i)*Cy(i) Cx(i)*Cz(i) ; Cy(i)*Cx(i) Cy(i)^2 Cy(i)*Cz(i) ;...
           Cz(i)*Cx(i) Cz(i)*Cy(i) Cz(i)^2];
    LAMDA(:,:,i) = [lamda -lamda ; -lamda lamda];
end
k = (E.*area)./L;                       %Constante AE/L para cada elemento
A = zeros(6,6);     %Inicializamos las matrices de rigidez de cada elemento

%Ensamble de la matriz global de rigidez
for i = 1:numelementos
    %Calculo de la matriz de rigidez de cada elemento (AE/L)*k_local
    A(:,:,i) = k(i)*LAMDA(:,:,i);
    
    %Convertimos el arreglo A en celdas pero las matrices de 6x6 de A las
    %dividimos en 4 paquetes de 3x3
    j = UnionNodos(i,:);                    %Par de índice de cada elemento
    B(:,:,i) = mat2cell(A(:,:,i),[3 3],[3 3]);
    
    %Asignamos cada paquete en los índices correspondientes de la matriz
    %global de rigidez
    C(j(1),j(1),i) = B(1,1,i);
    C(j(1),j(2),i) = B(1,2,i);
    C(j(2),j(1),i) = B(2,1,i);
    C(j(2),j(2),i) = B(2,2,i);
end

S = 3*size(nodos,1); m = cell(S/3,S/3);
for i = 1:size(nodos,1)
    for j = 1:size(nodos,1)
        %En cada vuelta recogemos todos los elementos con el mismo índice
        %(i,j), los superponemos y sumamos entre sí.
        clear x
        x(:,:,:) = cell2mat(reshape(C(i,j,:),1,[],numelementos));
        m(i,j) = {sum(x,3)};
        
        %Si en ese índice (i,j) no se asignó ningú paquete, metemos un paquete de ceros de 2x2
        if size(m{i,j}) == [0 0]
            m(i,j) = {zeros(3,3)};
        end
        
    end
end
MG = cell2mat(m);          %Convertimos la matriz global en un arreglo numérico

%Reducir la matriz global
v = find(Desplazamientos == 0);                
MGR = MG; MGR(v,:) = 0; MGR(:,v) = 0;
indicefil = zeros(1,S); indicecol = zeros(1,S);
for i = 1:S
    if MGR(i,:) == 0
        indicefil(i) = i;
    end
    if MGR(:,i) == 0
        indicecol(i) = i;
    end
end
MGR(indicefil~=0,:) = [];    %Eliminar filas y columnas de ceros para tener la matriz global reducida
MGR(:,indicecol~=0) = [];
Fuerzas(indicefil~=0) = [];  %Eliminar filas y columnas de ceros de las fuerzas

%Calcular los desplazamientos nodales en ambas direcciones
d = MGR\Fuerzas'; dfinal = zeros(S,1);
k = 1;
for i = 1:length(Desplazamientos) %Insertar los desplazamientos calculados en el vector de desplazamientos original
    if Desplazamientos(i) == 0
        dfinal(i,1) = 0;
    else
        dfinal(i,1) = d(k);
        k = k+1;
    end
    if k > length(d)
        break
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Resultados~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

d2 = mat2cell(dfinal,3*ones(1,size(nodos,1)),1); %Dividimos dfinal en paquetes de 3x1
Esfuerzos = zeros(1,numelementos); Flocal = zeros(numelementos,6);
j = 1;
for i = 1:numelementos
    indice = UnionNodos(i,:);
    Esfuerzos(i) = (E(i)./L(i)) * [-Cx(i) -Cy(i) -Cz(i) Cx(i) Cy(i) Cz(i)] * [d2{indice(1,1)} ; d2{indice(1,2)}];
%     Flocal(i,:) = A(:,:,i)*[d2{indice(1,1)} ; d2{indice(1,2)}];
    j = j + 2;
end
Reacciones = MG*dfinal;       

for i = 1:length(nodos)
    nstr{i} = ['Nodo ',num2str(i)];
end
for i = 1:numelementos
    bstr{i} = ['Elemento ', num2str(i)];
end

d3 = reshape(dfinal,[3,size(nodos,1)])';
Reac = reshape(Reacciones,[3,size(nodos,1)])';
T =  array2table([d3(:,1), d3(:,2), d3(:,3), Reac(:,1), Reac(:,2), Reac(:,3)],...
    'VariableNames',{'Desplazamientos nodales (X)'; 'Desplazamientos nodales (Y)'; ...
    'Desplazamientos nodales (Z)'; 'Reacciones nodales (X)'; 'Reacciones nodales (Y)';...
    'Reacciones nodales (Z)'}, 'RowNames',nstr)
T2 = array2table(Esfuerzos','VariableNames',{'Esfuerzos por barra'},...
    'RowNames',bstr)


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Gráfica~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

FAmp = 500;
NodosDesp = nodos + d3*FAmp;
for j = 1:numelementos
    indice = UnionNodos(j,:);
    ArmaduraOG.Barra(j) = line([nodos(indice(1),1) nodos(indice(2),1)],...
        [nodos(indice(1),2) nodos(indice(2),2)],...
        [nodos(indice(1),3) nodos(indice(2),3)],...
        'LineWidth',1.5,'Color','k');
    hold on
    ArmaduraN.Barra(j) = line([NodosDesp(indice(1),1) NodosDesp(indice(2),1)],...
        [NodosDesp(indice(1),2) NodosDesp(indice(2),2)],...
        [NodosDesp(indice(1),3) NodosDesp(indice(2),3)],...
        'LineWidth',1,'Color','b');
end

for i = 1:size(nodos,1)
    NodosStruct.N(i) = line(nodos(i,1),nodos(i,2),nodos(i,3),'Color','r',...
        'Marker','o','MarkerSize',6,'MarkerFaceColor','r');
    hold on
end

NodosDespx = NodosDesp(:,1); NodosDespy = NodosDesp(:,2); NodosDespz = NodosDesp(:,3);
[fil,col] = find(d3 ~= 0);
NodosDStruct.N = plot3(NodosDespx(fil),NodosDespy(fil),NodosDespz(fil),'og','MarkerSize',6,'MarkerFaceColor','g');

title(['Sistema resuelto. Amplificación x',num2str(FAmp)]); grid on
legend([ArmaduraOG.Barra(1), ArmaduraN.Barra(1), NodosStruct.N(1), NodosDStruct.N(1)], ...
    {'Armadura original'; 'Armadura deformada'; 'Nodos originales'; 'Nodos desplazados'});
axis equal
view(45,45)
