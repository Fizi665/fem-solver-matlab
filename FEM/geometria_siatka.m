clear all
close all
clc

% tworzenie geometrii
wierzcholki=[0,-1; 6,-1; 6,3; 5,3; 5,1; 1,1; 1,3; 0,3];
figure
plot(wierzcholki(:,1),wierzcholki(:,2),'-')
grid on
hold on
xlabel('x [cm]')
ylabel('y [cm]')

% tworzenie siatki
rozmiar=1;

% generowanie punktow na bazie maksymalnych rozmiarow 
[xGrid, yGrid] = meshgrid(min(wierzcholki(:,1)):rozmiar:max(wierzcholki(:,1)), min(wierzcholki(:,2)):rozmiar:max(wierzcholki(:,2)));
% wizualizacja punktów 
plot(xGrid, yGrid, 'r.'); 


% identyfikacja ktore punkty sa wewnatrz obiektu
in = inpolygon(xGrid, yGrid, wierzcholki(:,1), wierzcholki(:,2));
xInside = xGrid(in);
yInside = yGrid(in);

% wizualizacja punktów
plot(xInside, yInside, 'bo'); 

% znajdowanie par ktorych dystans od siebie = rozmiar
distances = pdist2([xInside, yInside], [xInside, yInside]);
distances = distances.*(tril(ones(length(distances))*inf)+1); % Eliminacja duplikatow przez nadanie wartosci inf macierzy trojkatnej
[rows, cols] = find(distances == rozmiar);

% tworzenie punktow miedzy parami
newX = [];
newY = [];
for i = 1:length(rows)
    x1 = xInside(rows(i));
    y1 = yInside(rows(i));
    x2 = xInside(cols(i));
    y2 = yInside(cols(i));
    
    newX = [newX; (x1 + x2) / 2];
    newY = [newY; (y1 + y2) / 2];
end

% wizualizacja punktów
plot(newX, newY, 'g*'); 


% zamiana cm na m (nie dziala jak da się od razu m)
wierzcholki = wierzcholki*0.01;
rozmiar=rozmiar*0.01;
newX = newX*0.01;
newY = newY*0.01;
xInside = xInside*0.01;
yInside = yInside*0.01;
figure
plot(wierzcholki(:,1),wierzcholki(:,2),'-',xInside, yInside, 'bo',newX, newY, 'g*')
grid on
%% tworzenie elementow
P = [xInside, yInside];
M = [newX, newY];

nodes = [P; M];

nP = size(P,1);   % liczba węzłów narożnych
findNode = @(pt) find( abs(nodes(:,1)-pt(1))<1e-8 & ...
                       abs(nodes(:,2)-pt(2))<1e-8 );
ELEMENTS = [];

dx = rozmiar;
dy = rozmiar;

for i = 1:nP
    x = P(i,1);
    y = P(i,2);

    % narożniki
    p1 = [x,     y];
    p2 = [x+dx,  y];
    p3 = [x+dx,  y+dy];
    p4 = [x,     y+dy];

    % sprawdzenie czy wszystkie narożniki istnieją
    if all( ismembertol([p1;p2;p3;p4], P, 1e-8, 'ByRows', true) )

        % środki boków
        m12 = (p1+p2)/2;
        m23 = (p2+p3)/2;
        m34 = (p3+p4)/2;
        m41 = (p4+p1)/2;

        % ID węzłów
        el = [
            findNode(p1)
            findNode(p2)
            findNode(p3)
            findNode(p4)
            findNode(m12)
            findNode(m23)
            findNode(m34)
            findNode(m41)
        ];

        ELEMENTS = [ELEMENTS; el'];
    end
end


figure; hold on; grid on
plot(nodes(:,1), nodes(:,2), 'ko')

for e = 1:size(ELEMENTS,1)
    el = ELEMENTS(e,[1:4 1]); % tylko narożniki
    plot(nodes(el,1), nodes(el,2),'b-','LineWidth',1)
end

for i=1:size(nodes,1)
    text(nodes(i,1)+rozmiar*0.05, nodes(i,2)+rozmiar*0.08, num2str(i))
end

%% inicjalizacja macierzy globalnych
load fem_symbolic

nN = size(nodes,1);

Kk = zeros(nN);
Kp = zeros(nN);
Kalfa = zeros(nN);
Rq = zeros(nN,1);
Rbeta = zeros(nN,1);

%% warunki objetosciowe
kx = 55;
ky = 55;
Q = 800;
p = 0;

elements_Q = [1 5 7 9 11 13];

for e = 1:size(ELEMENTS,1)
    el_nodes = ELEMENTS(e,:);
    
    % przewodzenie
    Ke = kk_fun(kx,ky);
    Kk(el_nodes,el_nodes) = Kk(el_nodes,el_nodes) + Ke;

    % pojemność (p=0 → nic nie wnosi)
    Kp(el_nodes,el_nodes) = Kp(el_nodes,el_nodes) + kp_fun(p);

    % źródło objętościowe
    if ismember(e,elements_Q)
        Rq(el_nodes) = Rq(el_nodes) + rq_fun(Q)';
    end
end

%% warunki brzegowe
% przenikanie
h = 85;
Tinf = 22;
beta = h*Tinf;

edge_alfa = [2 26];
isEdge = @(edge,el) all(ismember(edge,el));

% konwekcja
edges_beta = [
    2 5
    5 10
    10 8
    8 22
    22 24
    24 29
    29 26
];

%% agregacja brzegow
for e = 1:size(ELEMENTS,1)
    el = ELEMENTS(e,:);
    
    edges = {
        el([1 2 3])
        el([3 4 5])
        el([5 6 7])
        el([7 8 1])
    };

    for b = 1:4
        ed = edges{b};

        % alfa
        if isEdge(edge_alfa,ed)
            Kalfa(el,el) = Kalfa(el,el) + kalfa_fun{b}(h);
        end

        % beta
        for k=1:size(edges_beta,1)
            if isEdge(edges_beta(k,:),ed)
                Rbeta(el) = Rbeta(el) + rbeta_fun{b}(beta)';
            end
        end
    end
end

%% wynik
U = (Kk + Kp + Kalfa) \ (Rq + Rbeta);

%% wizualizacja
figure; hold on; axis equal
title('Rozkład temperatury + siatka FEM')
xlabel('x [cm]')
ylabel('y [cm]')

% kolorowanie węzłów
scatter(nodes(:,1), nodes(:,2), 60, U, 'filled')

% rysowanie elementów
for e = 1:size(ELEMENTS,1)
    el = ELEMENTS(e,[1:4 1]); % tylko narożniki
    plot(nodes(el,1), nodes(el,2), 'k-', 'LineWidth', 0.8)
end

colorbar
colormap(turbo)
grid on

