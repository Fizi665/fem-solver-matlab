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

nN = size(nodes,1);        % liczba węzłów globalnych
nE = size(ELEMENTS,1);     % liczba elementów

K = zeros(nN);             % globalna macierz sztywności
R = zeros(nN,1);           % globalny wektor obciążeń

%% parametry globalne
kx = 55;
ky = 55;
p  = 0;

Q     = 800;
h     = 85;
Tinf  = 22;

%% definicje warunków
elements_with_Q = [1 5 7 9 11 13];

edge_alfa  = [2 26];
edges_beta = [2 5; 5 10; 10 8; 8 22; 22 24; 24 29; 29 26];

tol = 1e-10;

segments_alfa = [];
segments_beta = [];

% alfa
segments_alfa = expandEdgeLine(nodes, edge_alfa(1), edge_alfa(2), tol);

% beta (wiele linii)
for k = 1:size(edges_beta,1)
    seg = expandEdgeLine(nodes, edges_beta(k,1), edges_beta(k,2), tol);
    segments_beta = [segments_beta; seg];
end


%% pętla po elementach
for e = 1:nE

    nodes_e = ELEMENTS(e,:);      % globalne ID węzłów elementu
    Xe = nodes(nodes_e,1);
    Ye = nodes(nodes_e,2);

    %% ---------- OBJĘTOŚĆ ----------
    if ismember(e, elements_with_Q)
        q_e = Q;
    else
        q_e = 0;
    end

    Kk_e = Kk_sym(kx,ky);
    Kp_e = Kp_sym(p);
    Rq_e = Rq_sym(q_e);

    %% ---------- BRZEGI ----------
    Kalfa_e = zeros(8);
    Rbeta_e = zeros(1,8);

    % lokalne boki: [1-2], [2-3], [3-4], [4-1]
    edgeLoc = [1 2; 2 3; 3 4; 4 1];

for ed = 1:4

    g1 = nodes_e(edgeLoc(ed,1));
    g2 = nodes_e(edgeLoc(ed,2));

    alfa_e = 0;
    beta_e = 0;

    % --- alfa ---
    for k = 1:size(segments_alfa,1)
        if (g1 == segments_alfa(k,1) && g2 == segments_alfa(k,2)) || ...
           (g2 == segments_alfa(k,1) && g1 == segments_alfa(k,2))
            alfa_e = h;
            break
        end
    end

    % --- beta ---
    for k = 1:size(segments_beta,1)
        if (g1 == segments_beta(k,1) && g2 == segments_beta(k,2)) || ...
           (g2 == segments_beta(k,1) && g1 == segments_beta(k,2))
            beta_e = h*Tinf;
            break
        end
    end

    Kalfa_e = Kalfa_e + Kalfa_edge{ed}(alfa_e);
    Rbeta_e = Rbeta_e + Rbeta_edge{ed}(beta_e);
end

    %% ---------- MACIERZE ELEMENTU ----------
    Ke = Kk_e + Kp_e + Kalfa_e;
    Re = Rq_e + Rbeta_e;

    %% ---------- AGREGACJA ----------
    for i = 1:8
        I = nodes_e(i);
        R(I) = R(I) + Re(i);
        for j = 1:8
            J = nodes_e(j);
            K(I,J) = K(I,J) + Ke(i,j);
        end
    end
end

%% ---------- ROZWIĄZANIE ----------
u = K \ R;

%% ---------- WIZUALIZACJA ----------
figure; hold on; axis equal
title('Rozkład temperatury + siatka FEM')
xlabel('x [cm]')
ylabel('y [cm]')

% kolorowanie węzłów
scatter(nodes(:,1), nodes(:,2), 60, u, 'filled')

% rysowanie elementów
for e = 1:size(ELEMENTS,1)
    el = ELEMENTS(e,[1:4 1]); % tylko narożniki
    plot(nodes(el,1), nodes(el,2), 'k-', 'LineWidth', 0.8)
end

colorbar
colormap(turbo)
grid on
%% funkcje pomocnicze
function segments = expandEdgeLine(nodes, nStart, nEnd, tol)

    P0 = nodes(nStart,:);
    P1 = nodes(nEnd,:);
    v  = P1 - P0;
    L  = norm(v);

    segments = [];

    for i = 1:size(nodes,1)
        Pi = nodes(i,:) - P0;

        % współliniowość (iloczyn wektorowy w 2D)
        if abs(det([v; Pi])) < tol
            t = dot(Pi,v)/dot(v,v);
            if t > -tol && t < 1+tol
                segments = [segments; i t];
            end
        end
    end

    % sortowanie wzdłuż linii
    segments = sortrows(segments,2);
    segments = segments(:,1);

    % tworzenie par [i i+1]
    segments = [segments(1:end-1) segments(2:end)];
end
