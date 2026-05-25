clear all
close all
clc
%szukany max q dla ktorego max(u)-35=0 (<0)
%% tu zaczynamy petle szukajaca q
Q=[];
f_u=[];
Q(1)=800;
f_u(1)=-12.96;
Q(2)=800000;
f_u(2)=25;
wi=3;
while f_u(wi-1) < -0.01 || f_u(wi-1)>0
Q(wi)=Q(wi-1)-((f_u(wi-1)*(Q(wi-1)-Q(wi-2)))/(f_u(wi-1)-f_u(wi-2)));

%% podstawy modelu
wierzcholki=[0,-1; 6,-1; 6,3; 5,3; 5,1; 1,1; 1,3; 0,3];
rozmiar=1; 
[xGrid, yGrid] = meshgrid(min(wierzcholki(:,1)):rozmiar:max(wierzcholki(:,1)), min(wierzcholki(:,2)):rozmiar:max(wierzcholki(:,2)));
in = inpolygon(xGrid, yGrid, wierzcholki(:,1), wierzcholki(:,2));
xInside = xGrid(in);
yInside = yGrid(in);
distances = pdist2([xInside, yInside], [xInside, yInside]);
distances = distances.*(tril(ones(length(distances))*inf)+1); % Eliminacja duplikatow przez nadanie wartosci inf macierzy trojkatnej
[rows, cols] = find(distances == rozmiar);
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

wierzcholki = wierzcholki*0.01;
rozmiar=rozmiar*0.01;
newX = newX*0.01;
newY = newY*0.01;
xInside = xInside*0.01;
yInside = yInside*0.01;
P = [xInside, yInside];
M = [newX, newY];
nodes = [P; M];
nP = size(P,1);
findNode = @(pt) find( abs(nodes(:,1)-pt(1))<1e-10 & ...
                       abs(nodes(:,2)-pt(2))<1e-10 );
ELEMENTS = [];
dx = rozmiar;
dy = rozmiar;
for i = 1:nP
    x = P(i,1);
    y = P(i,2);
    p1 = [x,     y];
    p2 = [x+dx,  y];
    p3 = [x+dx,  y+dy];
    p4 = [x,     y+dy];
    if all( ismembertol([p1;p2;p3;p4], P, 1e-8, 'ByRows', true) )
        m12 = (p1+p2)/2;
        m23 = (p2+p3)/2;
        m34 = (p3+p4)/2;
        m41 = (p4+p1)/2;
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
load fem_symbolic
nN = size(nodes,1);        
nE = size(ELEMENTS,1);     
K = zeros(nN);
R = zeros(nN,1);
kx = 55;
ky = 55;
p  = 0;
h     = 85;
Tinf  = 22+273.15;
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


% pętla po elementach
for e = 1:nE

    nodes_e = ELEMENTS(e,:);
    Xe = nodes(nodes_e,1);
    Ye = nodes(nodes_e,2);
    if min(Ye) <-rozmiar*0.5 
        q_e = Q(wi);
    else
        q_e = 0;
    end
    Kk_e = Kk_sym(kx,ky);
    Kp_e = Kp_sym(p);
    Rq_e = Rq_sym(q_e);
    Kalfa_e = zeros(8);
    Rbeta_e = zeros(8,1);
    edgeLoc = [1 2; 2 3; 3 4; 4 1];
    for ed = 1:4
        g1 = nodes_e(edgeLoc(ed,1));
        g2 = nodes_e(edgeLoc(ed,2));
        p1 = nodes(g1,:);
        p2 = nodes(g2,:);
        if edge_fully_on_poly(p1,p2,BC_poly,tol)
            alfa_e = h;
            beta_e = h*Tinf;
        else
            alfa_e = 0;
            beta_e = 0;
        end
        if alfa_e ~= 0
            Kalfa_e = Kalfa_e + Kalfa_edge{ed}(alfa_e);
        end
        if beta_e ~= 0
            Rbeta_e = Rbeta_e + Rbeta_edge{ed}(beta_e)';
        end 
    end
    Ke = Kk_e + Kp_e + Kalfa_e;
    Re = Rq_e' + Rbeta_e;
    temp=nodes_e;
    nodes_e(:,2)=temp(:,5);
    nodes_e(:,3)=temp(:,2);
    nodes_e(:,4)=temp(:,6);
    nodes_e(:,5)=temp(:,3);
    nodes_e(:,6)=temp(:,7);
    nodes_e(:,7)=temp(:,4);
    for i = 1:8
        I = nodes_e(i);
        R(I) = R(I) + Re(i);
        for j = 1:8
             J = nodes_e(j);
             K(I,J) = K(I,J) + Ke(i,j);
        end
    end
end
u = K\R;

f_u(wi)=max(u)-35-273.15
wi=wi+1;
end

odp=Q(wi-1)

%% funkcje pomocnicze
function onBC = edge_fully_on_poly(p1, p2, poly, tol)
onBC = false;
for i = 1:size(poly,1)-1
    q1 = poly(i,:);
    q2 = poly(i+1,:);
    if point_on_segment(p1,q1,q2,tol) && ...
       point_on_segment(p2,q1,q2,tol)
        onBC = true;
        return
    end
end
end
function tf = point_on_segment(p, a, b, tol)
v1 = p - a;
v2 = b - a;
if norm(cross([v1 0],[v2 0])) > tol
    tf = false;
    return
end
t = dot(v1,v2) / dot(v2,v2);
tf = (t > -tol) && (t < 1+tol);
end