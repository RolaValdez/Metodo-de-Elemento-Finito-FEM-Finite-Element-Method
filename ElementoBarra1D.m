%% Código para resolver sistemas de elementos barra en 1D
%Autor: Rolando Valdez Guzmán
%Alias: Tutoingeniero
%Canal de Youtube: https://www.youtube.com/channel/UCU1pdvVscOdtLpRQBp-TbWg
%Versión: 2.0
%Actualizado: 22/jul/2022

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

%kr = Vector con las constantes de resorte de cada elemento. Si un elemento
%es un resorte, se agrega el valor de la constante en el índice
%correspondiente. Si no hay ningún resorte, este es un vector de ceros.

%nodos = Son las coordenadas X de cada nodo en las unidades que defina tu
%problema. Recomiendo siempre poner el origen en el primer nodo.

%UnionNodos = Son los índices de los nodos que conforman a cada barra. Cada
%barra es unida por dos nodos, por lo que si la barra 1 empieza en 0 (nodo 1) y
%termina en 3 (nodo 2), la barra se define con los índices [1 2], y así
%para cada barra.

%Desplazamientos = Condiciones de frontera para cada nodo (0 si el nodo
%está empotrado y 1 si el nodo puede moverse).

%Fuerzas = Vector de fuerzas del sistema (0 si no hay una fuerza actuando
%en dicho nodo y un valor cualquiera si hay una fuerza. Dependiendo del
%sentido de cada fuerza, se deberá usar un signo negativo o positivo).

clc; clear ;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%

%Resuelve un sistema de barras en 1D con el método de elemento finito.
%Puedes usar cualquiera de los tres ejemplos que dejo aquí en esta sección,
%todos del libro que dejo arriba o basarte en los ejemplos para resolver tu
%propio sistema, sólo debes seguir la misma estructura.

%NOTA: Este código sólo funciona si se desconocen TODOS los desplazamientos
%nodales, el algoritmo NO resolverá sistemas con desplazamientos conocidos
%correctamente.

%Ejemplo 3.1 del libro 

numelementos=3;                     %Número de elementos
E=[30*10^6 30*10^6 15*10^6];        %Módulo de elasticidad por elemento
area=[1 1 2];                       %Areas de cada elemento
kr = [0, 0, 0];                  %Constante de resorte de cada elemento
nodos=[0 ; 30 ; 60 ; 90];           %Coordenadas de cada nodo
UnionNodos=[1 2 ; 2 3 ; 3 4];       %Nodos que unen a cada elemento
Desplazamientos=[0 1 1 0];          %Desplazamientos por nodo
Fuerzas=[0 3000 0 0];               %Fuerzas por nodo

%Problema 3.1 del libro

% numelementos=3;                     %Número de elementos
% E(1:numelementos)=10*10^6;        %Módulo de elasticidad por elemento
% area(1:numelementos)=1;                       %Areas de cada elemento
% kr = [0, 0, 0];                  %Constante de resorte de cada elemento
% nodos=[0 ; 10 ; 20 ; 30];           %Coordenadas de cada nodo
% UnionNodos=[1 2 ; 2 3 ; 3 4];       %Nodos que unen a cada elemento
% Desplazamientos=[0 1 1 0];          %Desplazamientos por nodo
% Fuerzas=[0 0 1000 0];               %Fuerzas por nodo

% numelementos=3;                     %Número de elementos
% E(1:numelementos)= 30*10^6;        %Módulo de elasticidad por elemento
% area(1:numelementos)=4;                       %Areas de cada elemento
% kr = [0, 0, 0];                  %Constante de resorte de cada elemento
% nodos=[0 ; 30 ; 60 ; 90];           %Coordenadas de cada nodo
% UnionNodos=[1 2 ; 2 3 ; 3 4];       %Nodos que unen a cada elemento
% Desplazamientos=[0 1 1 0];          %Desplazamientos por nodo
% Fuerzas=[0 4000 -10000 0];               %Fuerzas por nodo

%Ejemplo con resorte

% numelementos = 3;                     %Número de elementos
% E = [15*10^6, NaN, 15*10^6];        %Módulo de elasticidad por elemento
% area = [3,NaN,3];
% kr = [0, 5000, 0];                       %Areas de cada elemento
% nodos=[0 ; 30 ; 30 ; 60];           %Coordenadas de cada nodo
% UnionNodos=[1 2 ; 2 3 ; 3 4];       %Nodos que unen a cada elemento
% Desplazamientos=[0 1 1 0];          %Desplazamientos por nodo
% Fuerzas=[0 0 10000 0];               %Fuerzas por nodo

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Algoritmo~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Construcción de las matrices de rigidez por cada elemento
A = zeros(2,2,numelementos); L = diff(nodos);
for i = 1:numelementos
    if kr(i) ~= 0
        A(:,:,i) = kr(i)*[1 -1 ; -1 1];  %%Matrices de rigidez por elemento
    else
        A(:,:,i) = ((E(i).*area(i))./L(i))*[1 -1 ; -1 1];  %%Matrices de rigidez por elemento
    end
end

%Ensamblar las matrices de rigidez de cada elemento en la matriz global
S = size(nodos,1); MG = zeros(S,S);                       %Matriz global vacía
for i=1:numelementos
    m = UnionNodos(i,1);
    n = UnionNodos(i,2);
    MG([m n],[m n],i) = A(:,:,i);      %Metemos las matrices de rigidez en tres matrices globales vacías 
end
MG = sum(MG,3);                       %Superponemos todas las matrices globales para tener una sola

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
MGR(indicefil ~= 0,:) = [];    %Eliminar filas y columnas de ceros para tener la matriz global reducida
MGR(:,indicecol ~= 0) = [];
Fuerzas(indicefil ~= 0) = [];  %Eliminar filas y columnas de ceros de las fuerzas

%Calcular los desplazamientos por nodo
d = MGR\Fuerzas'; k = 1; d2 = zeros(S,1);
for i = 1:length(Desplazamientos) %Insertar los desplazamientos calculados en el vector de desplazamientos original
    if Desplazamientos(i) == 0
        d2(i,1) = 0;
    else
        d2(i,1) = d(k);
        k = k + 1;
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Resultados~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

Esfuerzos = zeros(1,numelementos);
for i = 1:numelementos     %Calcular los esfuerzos por elemento
    indice = UnionNodos(i,:);
    Esfuerzos(i) = (E(i)./L(i))*[-1 1]*[d2(indice(1,1));d2(indice(1,2))];
end

Reacciones = MG*d2;              %Calcular reacciones por elemento
for i = 1:length(nodos)
    nstr{i} = ['Nodo ',num2str(i)];
end
for i = 1:numelementos
    bstr{i} = ['Elemento ', num2str(i)];
end
T =  array2table([d2,Reacciones],'VariableNames',...
    {'Desplazamientos nodales'; 'Reacciones nodales'}, 'RowNames',nstr)
T2 = array2table(Esfuerzos','VariableNames',{'Esfuerzos por barra'},...
    'RowNames',bstr)


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Gráfica~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

FAmp = 1000;
Nocero = find(d2 ~= 0); NodosDesp = nodos' + d2'; NodosDesp2 = nodos' + d2'*FAmp;

line([nodos(1) nodos(end)],[0 0],'Color','k','LineWidth',2); hold on
plot(nodos',zeros(1,length(nodos)),'or','MarkerSize',6,'MarkerFaceColor','r');
hold on
plot(NodosDesp(Nocero),zeros(1,length(Nocero)),'og','MarkerSize',6,'MarkerFaceColor','g');
hold on
plot(NodosDesp2(Nocero),zeros(1,length(Nocero)),'ob','MarkerSize',6,'MarkerFaceColor','b');
title(['Sistema resuelto. Amplificación x',num2str(FAmp)]); grid on
legend({'Barra'; 'Nodos restringidos'; 'Nodos móviles'; 'Nodos desplazados'});
axis([nodos(1)-10, nodos(end)+10, -1, 1]);
