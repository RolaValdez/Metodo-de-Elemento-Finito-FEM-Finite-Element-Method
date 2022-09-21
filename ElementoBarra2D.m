%% Código para resolver sistemas de elementos barra en 2D
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

%nodos = Son las coordenadas [X,Y] de cada nodo en las unidades que defina tu
%problema. Recomiendo siempre poner el origen en el primer nodo. Ejemplo:
%[x1 y1 ; x2 y2 ; ...]

%UnionNodos = Son los índices de los nodos que conforman a cada barra. Cada
%barra se conforma por una línea que va de un nodo a otro.

%Desplazamientos = Condiciones de frontera para cada componente XY de cada nodo 
%(0 si el nodo está empotrado y 1 si el nodo puede moverse). Por ejemplo:
%[dx1 dy1 dx2 dy2 dx3 dy3 ...] = [0 0 1 1 0 0 ...]

%Fuerzas = Vector de fuerzas del sistema (0 si no hay una fuerza actuando
%en dicho nodo y un valor cualquiera si hay una fuerza. Dependiendo del
%sentido de cada fuerza, se deberá usar un signo negativo o positivo).
%[Fx1 Fy1 Fx2 Fy2 Fx3 Fy3 ...] = [0 0 1000 0 -1000 -5000 ...]

clc; clear;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%

%Resuelve un sistema de barras en 2D con el método de elemento finito.
%Puedes usar cualquiera de los tres ejemplos que dejo aquí en esta sección,
%todos del libro que dejo arriba o basarte en los ejemplos para resolver tu
%propio sistema, sólo debes seguir la misma estructura.

%NOTA: Este código sólo funciona si se desconocen TODOS los desplazamientos
%nodales, el algoritmo NO resolverá sistemas con desplazamientos conocidos
%correctamente.

%Ejemplo 3.5 del libro

% numelementos=3;                             %Número de elementos
% E(1:numelementos)=30*10^6;                  %Módulo de elasticidad por elemento
% area(1:numelementos)=2;                     %Areas de cada elemento
% nodos=[0 0 ; 0 120; 120 120; 120 0];        %Coordenadas de cada nodo
% UnionNodos=[1 2;1 3;1 4];                   %Nodos que unen a cada elemento
% Desplazamientos=[1 1 0 0 0 0 0 0];          %Desplazamientos por nodo
% Fuerzas=[0 -10000 0 0 0 0 0 0 ];            %Fuerzas por nodo

%Problema 3.22

% numelementos=3;
% E(1:numelementos)=10*10^6;
% area(1:numelementos)=1;
% nodos=[0 0 ; -100 100*tand(60); -100  0 ;-100 -100*tand(30)];
% UnionNodos=[1 2;1 3;1 4];
% Desplazamientos=[1 1 0 0 0 0 0 0 ];
% Fuerzas=[1000 1000 0 0 0 0 0 0];

%Problema 3.24

% numelementos = 6;
% E(1:numelementos) = 10*10^6;
% area(1:numelementos) = 1;
% nodos=[0 0 ; 240 0; 240 180; 0 180];
% UnionNodos=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
% Desplazamientos=[0 0 1 1 1 1 0 0];
% Fuerzas=[0 0 0 1000 0 1000 0 0];

%Problema 3.65

numelementos = 11;
E(1:numelementos) = 10*10^6;
area(1:numelementos) = 4;
nodos=[0 0 ; 0 108; 36 36; 36 72; 36 108; 72 72; 108 108];
UnionNodos=[1 2; 2 3; 2 4; 2 5; 1 3; 3 4; 4 5; 4 6; 5 6; 6 7; 5 7];
Desplazamientos=[0 1 0 0 1 1 1 1 1 1 1 1 1 1];
Fuerzas=[0 0 0 0 0 0 0 0 0 0 0 0 0 -40000];

% numelementos = 13;
% E(1:numelementos) = 1*10^6;
% area(1:numelementos)=0.1;
% nodos=[0 0 ; 5 5 ; 5 0 ; 10 5; 10 0; 15 5; 15 0; 20 0];
% UnionNodos=[1 2 ; 1 3 ; 2 3; 3 5; 3 4; 2 4; 4 5; 5 7; 4 7; 4 6; 6 7; 7 8; 6 8];
% Desplazamientos=[0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0];
% Fuerzas=[0 0 0 0 0 -1000 0 0 0 -1000 0 0 0 -1000 0 0];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Algoritmo~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Construcción de las matrices de rigidez de cada barra.
L = zeros(1,numelementos); grados = zeros(1,numelementos); 
for i = 1:numelementos
    indice = UnionNodos(i,:);
    b = nodos(indice(2),:);                               %Recogemos Y2 y Y1
    a = nodos(indice(1),:);                               %Recogemos X2 y X1
    L(i) = norm(b-a); %Calculamos longitud usando la distancia entre dos puntos
    
    %Determinamos el ángulo del elemento con respecto a la horizontal
    dx = nodos(indice(2),1) - nodos(indice(1),1);
    dy = nodos(indice(2),2) - nodos(indice(1),2);   
    %Determinamos en qué cuadrante se encuentra el elemento
    grados(i) = atand(dy/dx);                                 %Cuadrante 1
    if sign(dy) == 1 && sign(dx) == -1                          %Cuadrante 2
        grados(i) = 180 + grados(i);
    elseif dy == 0 && sign(dx) == -1 %Sobre la horizontal pero en el cuadrante 2
        grados(i) = 180;
    elseif sign(dy) == -1 && sign(dx) == -1                     %Cuadrante 3
        grados(i) = 180 + grados(i);
    elseif sign(dy) == -1 && sign(dx) == 1                      %Cuadrante 4
        grados(i) = 360 + grados(i);
    end  
end
k = (E.*area)./L;                        %Constante AE/L para cada elemento

%Ensamble de la matriz global de rigidez.
A = zeros(4,4,numelementos);     %Inicializamos las matrices de rigidez de cada elemento
for i = 1:numelementos
    %Calculo de la matriz de rigidez de cada elemento (AE/L)*k_local
    A(:,:,i) = k(i)*transformadalineal(grados(i));  %TRANSFORMADALINEAL se encuentra hasta abajo  
    
    %Convertimos el arreglo A en celdas pero las matrices de 4x4 de A las
    %dividimos en 4 paquetes de 2x2
    j = UnionNodos(i,:);                      %Par de índice de cada elemento
    B(:,:,i) = mat2cell(A(:,:,i),[2 2],[2 2]);
    
    %Asignamos cada paquete en los índices correspondientes de la matriz
    %global de rigidez
    C(j(1),j(1),i) = B(1,1,i);
    C(j(1),j(2),i) = B(1,2,i);
    C(j(2),j(1),i) = B(2,1,i);
    C(j(2),j(2),i) = B(2,2,i);
end

%Dimensiones de la matriz global
S = 2*size(nodos,1); m = cell(S/2,S/2);
for i=1:size(nodos,1)
    for j=1:size(nodos,1)
        %En cada vuelta recogemos todos los elementos con el mismo índice
        %(i,j), los superponemos y sumamos entre sí.
        clear x
        x(:,:,:) = cell2mat(reshape(C(i,j,:),1,[],numelementos));
        m(i,j) = {sum(x,3)};
        
        %Si en ese índice (i,j) no se asignó ningú paquete, metemos un paquete de ceros de 2x2
        if size(m{i,j}) == [0 0]  
            m(i,j) = {zeros(2,2)};
        end
    end
end
MG = cell2mat(m);         %Convertimos la matriz global en un arreglo numérico

%Reducir la matriz global
v = find(Desplazamientos == 0); MGR = MG;
MGR(v,:) = 0; MGR(:,v) = 0;
indicefil = zeros(1,S); indicecol = zeros(1,S);
for i = 1:S
    if MGR(i,:) == 0
        indicefil(i) = i;
    end
    if MGR(:,i) == 0
        indicecol(i) = i;
    end
end
MGR(indicefil~=0,:)=[];    %Eliminar filas y columnas de ceros para tener la matriz global reducida
MGR(:,indicecol~=0) = [];
Fuerzas(indicefil~=0)=[]; %Eliminar filas y columnas de ceros de las fuerzas

%Calcular los desplazamientos nodales en ambas direcciones
d = MGR\Fuerzas';
dfinal = zeros(S,1); k=1;
for i = 1:length(Desplazamientos) %Insertar los desplazamientos calculados en el vector de desplazamientos original
    if Desplazamientos(i) == 0
        dfinal(i,1) = 0;
    else
        dfinal(i,1) = d(k);
        k = k + 1;
    end
end

%Vector de coeficientes [-C -S C S] para calcular los esfuerzos
Ve = {@(x) -cosd(x) @(x) -sind(x) @(x) cosd(x) @(x) sind(x)};

esfuerzosC = zeros(numelementos,4);
for i = 1:numelementos
    for j = 1:4
        esfuerzosC(i,j) = feval(Ve{1,j},grados(i)); %Evaluar Ve para cada grado de cada elemento
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Resultados~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

d2 = mat2cell(dfinal,2*ones(1,size(nodos,1)),1); %Dividimos dfinal en paquetes de 2x1
Esfuerzos = zeros(1,numelementos); Flocal = zeros(numelementos,4);
j=1;
for i = 1:numelementos     %Calcular los esfuerzos por elemento
    indice = UnionNodos(i,:);
    Esfuerzos(i) = (E(i)./L(i))*esfuerzosC(i,:)*[d2{indice(1,1)};d2{indice(1,2)}];
%     Flocal(i,:) = A(:,:,i)*[d2{indice(1,1)};d2{indice(1,2)}]; 
    j = j+2;
end
Reacciones = MG*dfinal;       

for i = 1:length(nodos)
    nstr{i} = ['Nodo ',num2str(i)];
end
for i = 1:numelementos
    bstr{i} = ['Elemento ', num2str(i)];
end

d3 = reshape(dfinal,[2,size(nodos,1)])';
Reac = reshape(Reacciones,[2,size(nodos,1)])';
T =  array2table([d3(:,1), d3(:,2), Reac(:,1), Reac(:,2)],'VariableNames',...
    {'Desplazamientos nodales (X)'; 'Desplazamientos nodales (Y)'; ...
    'Reacciones nodales (X)'; 'Reacciones nodales (Y)'}, 'RowNames',nstr)
T2 = array2table(Esfuerzos','VariableNames',{'Esfuerzos por barra'},...
    'RowNames',bstr)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Gráfica~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
FAmp = 3;
NodosDesp = nodos + d3*FAmp;
for j = 1:numelementos
    indice = UnionNodos(j,:);
    ArmaduraOG.Barra(j) = line([nodos(indice(1),1) nodos(indice(2),1)],...
        [nodos(indice(1),2) nodos(indice(2),2)],...
        'LineWidth',1.5,'Color','k');
    hold on
    ArmaduraN.Barra(j) = line([NodosDesp(indice(1),1) NodosDesp(indice(2),1)],...
        [NodosDesp(indice(1),2) NodosDesp(indice(2),2)],...
        'LineWidth',1,'Color','b');
end

for i = 1:size(nodos,1)
    NodosStruct.N(i) = line(nodos(i,1),nodos(i,2),'Color','r',...
        'Marker','o','MarkerSize',6,'MarkerFaceColor','r');
    hold on
end

NodosDespx = NodosDesp(:,1); NodosDespy = NodosDesp(:,2);
Nocerox = find(d3(:,1)~=0); Noceroy = find(d3(:,2)~=0);
if length(Nocerox) < length(Noceroy)
    for i = 1:length(Noceroy)
        NodosDStruct.N(i) = line(NodosDespx(Noceroy(i)),NodosDespy(Noceroy(i)),...
            'Color','g','Marker','o','MarkerSize',6,'MarkerFaceColor','g');
    end
else
    for i = 1:length(Noceroy)
        NodosDStruct.N(i) = line(NodosDespx(Nocerox(i)),NodosDespy(Nocerox(i)),...
            'Color','g','Marker','o','MarkerSize',6,'MarkerFaceColor','g');
    end
end

title(['Sistema resuelto. Amplificación x',num2str(FAmp)]); grid on
legend([ArmaduraOG.Barra(1), ArmaduraN.Barra(1), NodosStruct.N(1), NodosDStruct.N(1)], ...
    {'Armadura original'; 'Armadura deformada'; 'Nodos originales'; 'Nodos desplazados'});
lim = max(nodos,[],'all');
axis([(0 -lim/10), (lim + lim/10), (-lim/2 -lim/10), (lim/2 + lim/10)]);

%~~~~~~~~~~~~~~~~~~~~~~~~~función transformadalineal~~~~~~~~~~~~~~~~~~~~~~%

function [M] = transformadalineal(theta)
m={@(x) (cosd(x))^2 @(x) sind(x)*cosd(x) @(x) -(cosd(x))^2 @(x) -sind(x)*cosd(x);...
    @(x) sind(x)*cosd(x) @(x) (sind(x))^2 @(x) -sind(x)*cosd(x) @(x) -(sind(x))^2;...
    @(x) -(cosd(x))^2  @(x) -sind(x)*cosd(x) @(x) (cosd(x))^2 @(x) sind(x)*cosd(x);...
    @(x) -sind(x)*cosd(x) @(x) -(sind(x))^2 @(x) sind(x)*cosd(x) @(x) (sind(x))^2};

M = zeros(4,4);
for i = 1:4
    for j = 1:4
        M(i,j) = feval(m{i,j},theta);
    end
end
end
end
