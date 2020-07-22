%% Código para resolver sistemas de elementos tipo barra en 2D
%Autor: Rolando Valdez Guzmán
%Alias: Tutoingeniero
%Canal de Youtube: https://www.youtube.com/channel/UCU1pdvVscOdtLpRQBp-TbWg
%Versión: 1.0
%Actualizado: 22/jul/2020

%Referencias: "A First Course in the Finite Element Method" por Daryl. L.
%Logan

%% Variables:
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
%[dx1 dy1 dx2 dy2 dx3 dy3 ...]

%Fuerzas = Vector de fuerzas del sistema (0 si no hay una fuerza actuando
%en dicho nodo y un valor cualquiera si hay una fuerza. Dependiendo del
%sentido de cada fuerza, se deberá usar un signo negativo o positivo).
%[Fx1 Fy1 Fx2 Fy2 Fx3 Fy3 ...]

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

numelementos=3;                             %Número de elementos
E(1:numelementos)=30*10^6;                  %Módulo de elasticidad por elemento
area(1:numelementos)=2;                     %Areas de cada elemento
nodos=[0 0 ; 0 120; 120 120; 120 0];        %Coordenadas de cada nodo
UnionNodos=[1 2;1 3;1 4];                   %Nodos que unen a cada elemento
Desplazamientos=[1 1 0 0 0 0 0 0];          %Desplazamientos por nodo
Fuerzas=[0 -10000 0 0 0 0 0 0 ];            %Fuerzas por nodo

%Problema 3.22

% numelementos=3;
% E(1:numelementos)=10*10^6;
% area(1:numelementos)=1;
% nodos=[0 0 ; -100 100*tand(60); -100  0 ;-100 -100*tand(30)];
% UnionNodos=[1 2;1 3;1 4];
% Desplazamientos=[1 1 0 0 0 0 0 0 ];
% Fuerzas=[1000 1000 0 0 0 0 0 0];

%Problema 3.29

% numelementos=3;
% E(1:numelementos)=210*10^9;
% area(1:numelementos)=0.0004;
% nodos=[0 0 ; 0 2 ; -3 0 ; -5*sind(30) -5*sind(60)];
% UnionNodos=[1 2 ; 1 3 ; 1 4];
% Desplazamientos=[1 1 0 0 0 0 0 0];
% Fuerzas=[0 -40000 0 0 0 0 0 0];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Algoritmo~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Construcción de las matrices de rigidez de cada barra.

L=zeros(1,numelementos);              %Inicializamos el vector de longitudes
grados=zeros(1,numelementos);         %Inicializamos el vector de inclinaciones
for i=1:numelementos
    indice=UnionNodos(i,:);
    b=nodos(indice(2),:);                               %Recogemos Y2 y Y1
    a=nodos(indice(1),:);                               %Recogemos X2 y X1
    L(i)=norm(b-a); %Calculamos longitud usando la distancia entre dos puntos
    
    %Determinamos el ángulo del elemento con respecto a la horizontal
    dx=nodos(indice(2),1) - nodos(indice(1),1);
    dy=nodos(indice(2),2) - nodos(indice(1),2);   
    %Determinamos en qué cuadrante se encuentra el elemento
    grados(i)=atand(dy/dx);                                 %Cuadrante 1
    if sign(dy)==1 && sign(dx)==-1                          %Cuadrante 2
        grados(i)=180+grados(i);
    elseif dy==0 && sign(dx)==-1 %Sobre la horizontal pero en el cuadrante 2
        grados(i)=180;
    elseif sign(dy)==-1 && sign(dx)==-1                     %Cuadrante 3
        grados(i)=180+grados(i);
    elseif sign(dy)==-1 && sign(dx)==1                      %Cuadrante 4
        grados(i)=360+grados(i);
    end  
end
L
grados
k=(E.*area)./L                           %Constante AE/L para cada elemento

%Ensamble de la matriz global de rigidez.

A=zeros(4,4,numelementos);     %Inicializamos las matrices de rigidez de cada elemento
for i=1:numelementos
    %Calculo de la matriz de rigidez de cada elemento (AE/L)*k_local
    
    A(:,:,i)=k(i)*transformadalineal(grados(i));  %NECESITAS LA FUNCION TRANSFORMADALINEAL!  
    
    %Convertimos el arreglo A en celdas pero las matrices de 4x4 de A las
    %dividimos en 4 paquetes de 2x2
    j=UnionNodos(i,:);                      %Par de índice de cada elemento
    B(:,:,i)=mat2cell(A(:,:,i),[2 2],[2 2]);
    
    %Asignamos cada paquete en los índices correspondientes de la matriz
    %global de rigidez
    C(j(1),j(1),i)=B(1,1,i);
    C(j(1),j(2),i)=B(1,2,i);
    C(j(2),j(1),i)=B(2,1,i);
    C(j(2),j(2),i)=B(2,2,i);
end
A

S=2*size(nodos,1);                          %Dimensiones de la matriz global
m=cell(S/2,S/2);
for i=1:size(nodos,1)
    for j=1:size(nodos,1)
        
        %En cada vuelta recogemos todos los elementos con el mismo índice
        %(i,j), los superponemos y sumamos entre sí.
        clear x
        x(:,:,:)=cell2mat(reshape(C(i,j,:),1,[],numelementos));
        m(i,j)={sum(x,3)};
        
        %Si en ese índice (i,j) no se asignó ningú paquete, metemos un paquete de ceros de 2x2
        if size(m{i,j})==[0 0]  
            m(i,j)={zeros(2,2)};
        end
    end
end
MG=cell2mat(m)          %Convertimos la matriz global en un arreglo numérico

%Reducir la matriz global
v=find(Desplazamientos==0);                   
MGR=MG;
MGR(v,:)=0;     
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
MGR(indicefil~=0,:)=[];    %Eliminar filas y columnas de ceros para tener la matriz global reducida
MGR(:,indicecol~=0)=[]
Fuerzas(indicefil~=0)=[]  %Eliminar filas y columnas de ceros de las fuerzas

%Calcular los desplazamientos nodales en ambas direcciones
d=MGR\Fuerzas';         
dfinal=zeros(S,1);
k=1;
for i=1:length(Desplazamientos) %Insertar los desplazamientos calculados en el vector de desplazamientos original
    if Desplazamientos(i)==0
        dfinal(i,1)=0;
    else
        dfinal(i,1)=d(k);
        k=k+1;
    end
end

%Vector de coeficientes [-C -S C S] para calcular los esfuerzos
Ve={@(x) -cosd(x) @(x) -sind(x) @(x) cosd(x) @(x) sind(x)};

esfuerzosC=zeros(numelementos,4);
for i=1:numelementos
    for j=1:4
        esfuerzosC(i,j)=feval(Ve{1,j},grados(i)); %Evaluar Ve para cada grado de cada elemento
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Resultados~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

d2=mat2cell(dfinal,2*ones(1,size(nodos,1)),1); %Dividimos dfinal en paquetes de 2x1
Esfuerzos=zeros(1,numelementos);
Flocal=zeros(numelementos,4);
j=1;
for i=1:numelementos     %Calcular los esfuerzos por elemento
    indice=UnionNodos(i,:);
    Esfuerzos(i)=(E(i)./L(i))*esfuerzosC(i,:)*[d2{indice(1,1)};d2{indice(1,2)}];
    Flocal(i,:)=A(:,:,i)*[d2{indice(1,1)};d2{indice(1,2)}]; 
    j=j+2;
end

dfinal                              %Imprimir resultados
Esfuerzos
Reacciones=MG*dfinal               %Reacciones globales
Flocal                    %Reacciones locales por elemento

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Gráfica~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

d3=reshape(dfinal,[2,size(nodos,1)])';
NodosDesp=nodos + d3;

for j=1:numelementos
    indice=UnionNodos(j,:);
    line([nodos(indice(1),1) nodos(indice(2),1)],...
        [nodos(indice(1),2) nodos(indice(2),2)],...
        'LineWidth',1.5,'Color','k');
    hold on
    line([NodosDesp(indice(1),1) NodosDesp(indice(2),1)],...
        [NodosDesp(indice(1),2) NodosDesp(indice(2),2)],...
        'LineWidth',1,'Color','b');
end

for i=1:size(nodos,1)
    plot(nodos(i,1),nodos(i,2),'ro','MarkerSize',6,'MarkerFaceColor','r');
    hold on
end
axis([min(nodos(:,1))-10 max(nodos(:,2))+10 min(nodos(:,1))-10 max(nodos(:,2)+10)])
grid on
zoom on
NodosDespx=NodosDesp(:,1);
NodosDespy=NodosDesp(:,2);
Nocerox=find(d3(:,1)~=0);
Noceroy=find(d3(:,2)~=0);
plot(NodosDespx(Nocerox),NodosDespy(Noceroy),'og','MarkerSize',6,'MarkerFaceColor','g');

