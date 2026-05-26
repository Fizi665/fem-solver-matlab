% =========================================================================
% fem_heat_solver.m
% =========================================================================
% Author: Wiktor Komorowski
%
% Description:
% 2D steady-state heat transfer simulation using the Finite Element Method.
%
% Features:
% - Structured mesh generation
% - 8-node quadrilateral finite elements
% - Internal heat generation
% - Convective boundary conditions
% - Temperature field visualization
%
% =========================================================================

clear
close all
clc

%% Geometry definition

vertices = [
    0,-1
    6,-1
    6,3
    5,3
    5,1
    1,1
    1,3
    0,3
];

%% Geometry visualization

figure

plot(vertices(:,1), vertices(:,2), '-', 'LineWidth', 3)

grid on
hold on

xlabel('x [cm]')
ylabel('y [cm]')

%% Mesh generation

elementSize = 0.25;

[xGrid, yGrid] = meshgrid( ...
    min(vertices(:,1)):elementSize:max(vertices(:,1)), ...
    min(vertices(:,2)):elementSize:max(vertices(:,2)));

plot(xGrid, yGrid, 'r.', 'MarkerSize', 15)

%% Identify points inside the geometry

insideDomain = inpolygon( ...
    xGrid, yGrid, ...
    vertices(:,1), vertices(:,2));

xInside = xGrid(insideDomain);
yInside = yGrid(insideDomain);

plot(xInside, yInside, 'bo', 'LineWidth', 3)

%% Generate midside nodes

distances = pdist2( ...
    [xInside, yInside], ...
    [xInside, yInside]);

distances = distances .* ...
    (tril(ones(length(distances))*inf) + 1);

[rows, cols] = find(distances == elementSize);

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

plot(newX, newY, 'g*', 'LineWidth', 3)

%% Convert dimensions from cm to m

vertices = vertices * 0.01;
elementSize = elementSize * 0.01;

newX = newX * 0.01;
newY = newY * 0.01;

xInside = xInside * 0.01;
yInside = yInside * 0.01;

%% Create node arrays

cornerNodes = [xInside, yInside];
midNodes = [newX, newY];

nodes = [cornerNodes; midNodes];

numberOfCornerNodes = size(cornerNodes,1);

findNode = @(pt) find( ...
    abs(nodes(:,1)-pt(1)) < 1e-10 & ...
    abs(nodes(:,2)-pt(2)) < 1e-10);

%% Element generation

ELEMENTS = [];

dx = elementSize;
dy = elementSize;

for i = 1:numberOfCornerNodes

    x = cornerNodes(i,1);
    y = cornerNodes(i,2);

    p1 = [x,     y];
    p2 = [x+dx,  y];
    p3 = [x+dx,  y+dy];
    p4 = [x,     y+dy];

    if all(ismembertol([p1;p2;p3;p4], ...
            cornerNodes, ...
            1e-8, ...
            'ByRows', true))

        m12 = (p1+p2)/2;
        m23 = (p2+p3)/2;
        m34 = (p3+p4)/2;
        m41 = (p4+p1)/2;

        element = [
            findNode(p1)
            findNode(p2)
            findNode(p3)
            findNode(p4)
            findNode(m12)
            findNode(m23)
            findNode(m34)
            findNode(m41)
        ];

        ELEMENTS = [ELEMENTS; element'];
    end
end

%% Mesh visualization

figure
hold on
grid on

plot(nodes(:,1), nodes(:,2), 'ko', 'LineWidth', 3)

xlabel('x [m]')
ylabel('y [m]')

for e = 1:size(ELEMENTS,1)

    element = ELEMENTS(e,[1:4 1]);

    plot(nodes(element,1), ...
         nodes(element,2), ...
         'b-', 'LineWidth', 3)
end

%% Load symbolic FEM functions

load fem_symbolic

%% Global matrices initialization

numberOfNodes = size(nodes,1);
numberOfElements = size(ELEMENTS,1);

K = zeros(numberOfNodes);
R = zeros(numberOfNodes,1);

%% Material and thermal parameters

kx = 55;
ky = 55;

p = 0;

Q = 800;

h = 85;
Tinf = 22 + 273.15;

%% Convection boundary polygon

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

tolerance = 1e-8;

%% FEM assembly loop

for e = 1:numberOfElements

    elementNodes = ELEMENTS(e,:);

    Xe = nodes(elementNodes,1);
    Ye = nodes(elementNodes,2);

    %% Internal heat source

    if min(Ye) < -elementSize * 0.5
        q_e = Q;
    else
        q_e = 0;
    end

    %% Volume matrices

    Kk_e = Kk_sym(kx,ky);
    Kp_e = Kp_sym(p);

    Rq_e = Rq_sym(q_e);

    %% Boundary matrices

    Kalfa_e = zeros(8);
    Rbeta_e = zeros(8,1);

    edgeConnectivity = [
        1 2
        2 3
        3 4
        4 1
    ];

    for edge = 1:4

        g1 = elementNodes(edgeConnectivity(edge,1));
        g2 = elementNodes(edgeConnectivity(edge,2));

        p1 = nodes(g1,:);
        p2 = nodes(g2,:);

        if edge_fully_on_poly(p1,p2,BC_poly,tolerance)

            alpha_e = h;
            beta_e = h * Tinf;

        else

            alpha_e = 0;
            beta_e = 0;
        end

        if alpha_e ~= 0
            Kalfa_e = Kalfa_e + ...
                Kalfa_edge{edge}(alpha_e);
        end

        if beta_e ~= 0
            Rbeta_e = Rbeta_e + ...
                Rbeta_edge{edge}(beta_e)';
        end
    end

    %% Element matrices

    Ke = Kk_e + Kp_e + Kalfa_e;
    Re = Rq_e' + Rbeta_e;

    %% Node reordering

    temp = elementNodes;

    elementNodes(:,2) = temp(:,5);
    elementNodes(:,3) = temp(:,2);
    elementNodes(:,4) = temp(:,6);
    elementNodes(:,5) = temp(:,3);
    elementNodes(:,6) = temp(:,7);
    elementNodes(:,7) = temp(:,4);

    %% Global assembly

    for i = 1:8

        I = elementNodes(i);

        R(I) = R(I) + Re(i);

        for j = 1:8

            J = elementNodes(j);

            K(I,J) = K(I,J) + Ke(i,j);
        end
    end
end

%% Solve FEM system

u = K \ R;

%% Temperature limits

u_min = min(u) - 273.15
u_max = max(u) - 273.15

%% Temperature field visualization

ksi = linspace(-1,1,20);
eta = linspace(-1,1,20);

[KS, ET] = meshgrid(ksi, eta);

%% Bilinear interpolation functions

N1 = 0.25 * (1 - KS) .* (1 - ET);
N2 = 0.25 * (1 + KS) .* (1 - ET);
N3 = 0.25 * (1 + KS) .* (1 + ET);
N4 = 0.25 * (1 - KS) .* (1 + ET);

%% Temperature distribution plot

figure

hold on
axis equal

title('Temperature Distribution [deg C]')

xlabel('x [m]')
ylabel('y [m]')

for e = 1:size(ELEMENTS,1)

    element = ELEMENTS(e,1:4);

    %% Element node coordinates

    Xe = nodes(element,1);
    Ye = nodes(element,2);

    %% Element temperatures

    Te = u(element) - 273.15;

    %% Temperature interpolation

    Tloc = ...
        N1 * Te(1) + ...
        N2 * Te(2) + ...
        N3 * Te(3) + ...
        N4 * Te(4);

    %% Geometry interpolation

    Xloc = ...
        N1 * Xe(1) + ...
        N2 * Xe(2) + ...
        N3 * Xe(3) + ...
        N4 * Xe(4);

    Yloc = ...
        N1 * Ye(1) + ...
        N2 * Ye(2) + ...
        N3 * Ye(3) + ...
        N4 * Ye(4);

    %% Surface plot

    surf( ...
        Xloc, ...
        Yloc, ...
        zeros(size(Tloc)), ...
        Tloc, ...
        'EdgeColor', 'none');
end

view(2)

colorbar
colormap(turbo)

xlabel('x [m]')
ylabel('y [m]')

%% =========================================================================
% Helper functions
% =========================================================================

function onBoundary = edge_fully_on_poly( ...
    p1, p2, polygon, tolerance)

% Returns true if the entire edge lies on
% a single polygon boundary segment

onBoundary = false;

for i = 1:size(polygon,1)-1

    q1 = polygon(i,:);
    q2 = polygon(i+1,:);

    if point_on_segment(p1,q1,q2,tolerance) && ...
       point_on_segment(p2,q1,q2,tolerance)

        onBoundary = true;
        return
    end
end
end

%% -------------------------------------------------------------------------

function isOnSegment = point_on_segment( ...
    p, a, b, tolerance)

% Returns true if point p lies on segment a-b

v1 = p - a;
v2 = b - a;

%% Collinearity check

if norm(cross([v1 0],[v2 0])) > tolerance

    isOnSegment = false;
    return
end

%% Segment inclusion check

t = dot(v1,v2) / dot(v2,v2);

isOnSegment = ...
    (t > -tolerance) && ...
    (t < 1 + tolerance);
end