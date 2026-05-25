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
rozmiar=0.25; %działa bez zarzutow dla 1,0.5,0.25

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


% zamiana cm na m (nie dziala jak da się od razu m, trzeba by dodać tolerancję)
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
findNode = @(pt) find( abs(nodes(:,1)-pt(1))<1e-10 & ...
                       abs(nodes(:,2)-pt(2))<1e-10 );
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
Tinf  = 22+273;

%% definicje warunków
%nowe
BC_poly = [
  0.00 0.00
  0.00 0.03
  0.01 0.03
  0.01 0.01
  0.05 0.01
  0.05 0.03
  0.06 0.03
  0.06 0.00
];

tol = 1e-8;
%stare
% elements_with_Q = [1 5 7 9 11 13];
% 
% edge_alfa  = [2 26];
% edges_beta = [2 5; 5 10; 10 8; 8 22; 22 24; 24 29; 29 26];
% 
% tol = 1e-10;
% 
% % --- ALFA ---
% alfa_nodes = nodes_on_line(edge_alfa(1), edge_alfa(2), nodes, P, tol);
% alfa_segments = [alfa_nodes(1:end-1)' alfa_nodes(2:end)'];
% 
% % --- BETA ---
% beta_segments = [];
% for k = 1:size(edges_beta,1)
%     bn = nodes_on_line(edges_beta(k,1), edges_beta(k,2), nodes, P, tol);
%     beta_segments = [beta_segments;
%                      bn(1:end-1)' bn(2:end)'];
% end
% alfa_segments=beta_segments;


%% pętla po elementach
for e = 1:nE

    nodes_e = ELEMENTS(e,:);      % globalne ID węzłów elementu
    Xe = nodes(nodes_e,1);
    Ye = nodes(nodes_e,2);

    %% ---------- OBJĘTOŚĆ ----------
    %nowe
if min(Ye) <-rozmiar*0.5
    q_e = Q;
else
    q_e = 0;
end
    %było
    % if ismember(e, elements_with_Q)
    %     q_e = Q;
    % else
    %     q_e = 0;
    % end


    Kk_e = Kk_sym(kx,ky);
    Kp_e = Kp_sym(p);
    Rq_e = Rq_sym(q_e);

    %% ---------- BRZEGI ----------
    Kalfa_e = zeros(8);
    Rbeta_e = zeros(8,1);

    % lokalne boki: [1-3], [3-5], [5-7], [7-1]
    edgeLoc = [1 2; 2 3; 3 4; 4 1];

for ed = 1:4

    g1 = nodes_e(edgeLoc(ed,1));
    g2 = nodes_e(edgeLoc(ed,2));

%nowe
p1 = nodes(g1,:);
p2 = nodes(g2,:);

if edge_fully_on_poly(p1,p2,BC_poly,tol)
    alfa_e = h;
    beta_e = h*Tinf;
else
    alfa_e = 0;
    beta_e = 0;
end

%stare
% %% --- alfa ---
% alfa_e = 0;
% for k = 1:size(alfa_segments,1)
%     if (g1 == alfa_segments(k,1) && g2 == alfa_segments(k,2)) || ...
%        (g2 == alfa_segments(k,1) && g1 == alfa_segments(k,2))
%         alfa_e = h;
%         break
%     end
% end
% 
% %% --- beta ---
% beta_e = 0;
% for k = 1:size(beta_segments,1)
%     if (g1 == beta_segments(k,1) && g2 == beta_segments(k,2)) || ...
%        (g2 == beta_segments(k,1) && g1 == beta_segments(k,2))
%         beta_e = h*Tinf;
%         break
%     end
% end


% nowe
if alfa_e ~= 0
    Kalfa_e = Kalfa_e + Kalfa_edge{ed}(alfa_e);
end

if beta_e ~= 0
    Rbeta_e = Rbeta_e + Rbeta_edge{ed}(beta_e)';
end
%było
%Kalfa_e = Kalfa_e + Kalfa_edge{ed}(alfa_e);
%Rbeta_e = Rbeta_e + Rbeta_edge{ed}(beta_e)'
% %;
end

    %% ---------- MACIERZE ELEMENTU ----------
    Ke = Kk_e + Kp_e + Kalfa_e;
    %Ke = Kk_e + Kalfa_e;
    Re = Rq_e' + Rbeta_e;

%do pokazania gdzie jest konwekcja
% for i=1:size(BC_poly,1)-1
%     plot(BC_poly(i:i+1,1),BC_poly(i:i+1,2),'r','LineWidth',3)
% end

% zamiana kolejnosci pkt w elemencie
temp=nodes_e;
nodes_e(:,2)=temp(:,5);
nodes_e(:,3)=temp(:,2);
nodes_e(:,4)=temp(:,6);
nodes_e(:,5)=temp(:,3);
nodes_e(:,6)=temp(:,7);
nodes_e(:,7)=temp(:,4);
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
u = K\R;

%% ---------- WIZUALIZACJA ----------
% figure; hold on; axis equal
% title('Rozkład temperatury + siatka FEM')
% xlabel('x [m]')
% ylabel('y [m]')
% 
% % kolorowanie węzłów
% scatter(nodes(:,1), nodes(:,2), 60, u-273, 'filled')
% 
% % rysowanie elementów
% for e = 1:size(ELEMENTS,1)
%     el = ELEMENTS(e,[1:4 1]); % tylko narożniki
%     plot(nodes(el,1), nodes(el,2), 'k-', 'LineWidth', 0.8)
% end
% 
% colorbar
% % colormap(turbo)
% grid on

%%%% lepsza wizualizacja
ksi = linspace(-1,1,20);
eta = linspace(-1,1,20);
[KS, ET] = meshgrid(ksi, eta);
N1 = 0.25*(1-KS).*(1-ET);
N2 = 0.25*(1+KS).*(1-ET);
N3 = 0.25*(1+KS).*(1+ET);
N4 = 0.25*(1-KS).*(1+ET);
figure; hold on; axis equal
title('Rozkład temperatury (FEM)')
xlabel('x [m]')
ylabel('y [m]')

for e = 1:size(ELEMENTS,1)
    el = ELEMENTS(e,1:4);

    % współrzędne węzłów elementu
    Xe = nodes(el,1);
    Ye = nodes(el,2);
    Te = u(el) - 273;

    % interpolacja temperatury
    Tloc = N1*Te(1) + N2*Te(2) + N3*Te(3) + N4*Te(4);

    % interpolacja geometrii
    Xloc = N1*Xe(1) + N2*Xe(2) + N3*Xe(3) + N4*Xe(4);
    Yloc = N1*Ye(1) + N2*Ye(2) + N3*Ye(3) + N4*Ye(4);

    surf(Xloc, Yloc, zeros(size(Tloc)), Tloc, ...
        'EdgeColor','none');
end

view(2)
colorbar
colormap(turbo)
% %%% wizualizacja 3
% figure; hold on; axis equal
% title('Rozkład temperatury')
% xlabel('x [m]')
% ylabel('y [m]')
% 
% patch('Faces', ELEMENTS(:,1:4), ...
%       'Vertices', nodes, ...
%       'FaceVertexCData', u-273, ...
%       'FaceColor','interp', ...
%       'EdgeColor','k');
% 
% colorbar
% colormap(turbo)


%% funkcje pomocnicze
%nowe
function onBC = edge_fully_on_poly(p1, p2, poly, tol)
% TRUE tylko jeśli CAŁY bok elementu leży na JEDNYM odcinku BC

onBC = false;

for i = 1:size(poly,1)-1
    q1 = poly(i,:);
    q2 = poly(i+1,:);

    % oba punkty boku muszą leżeć na tym samym odcinku BC
    if point_on_segment(p1,q1,q2,tol) && ...
       point_on_segment(p2,q1,q2,tol)
        onBC = true;
        return
    end
end
end

function tf = point_on_segment(p, a, b, tol)
% punkt p leży na odcinku a-b

v1 = p - a;
v2 = b - a;

% współliniowość
if norm(cross([v1 0],[v2 0])) > tol
    tf = false;
    return
end

% zawarcie w przedziale
t = dot(v1,v2) / dot(v2,v2);
tf = (t > -tol) && (t < 1+tol);
end


%stare
% function line_nodes = nodes_on_line(n1, n2, nodes, P, tol)
%     % n1, n2 – węzły końcowe linii
%     % nodes – wszystkie węzły
%     % P – tylko węzły narożne
%     % tol – tolerancja
% 
%     p1 = nodes(n1,:);
%     p2 = nodes(n2,:);
% 
%     dir = p2 - p1;
%     L = norm(dir);
%     dir = dir / L;
% 
%     line_nodes = [];
% 
%     for i = 1:size(P,1)
%         pi = P(i,:);
%         v = pi - p1;
% 
%         % czy punkt leży na prostej
%         if norm(cross([v 0],[dir 0])) < tol
%             % czy jest między n1 i n2
%             t = dot(v,dir);
%             if t > -tol && t < L+tol
%                 line_nodes(end+1) = i; %#ok<AGROW>
%             end
%         end
%     end
% 
%     % sortowanie wzdłuż linii
%     coords = nodes(line_nodes,:);
%     proj = (coords - p1)*dir';
%     [~,idx] = sort(proj);
%     line_nodes = line_nodes(idx);
% end
