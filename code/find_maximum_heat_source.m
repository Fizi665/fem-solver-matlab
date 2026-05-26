% =========================================================================
% find_maximum_heat_source.m
% =========================================================================
% Author: Wiktor Komorowski
%
% Description:
% Determines the maximum internal heat generation value using
% the secant method for a 2D steady-state FEM heat transfer model.
%
% Target condition:
%   max(T) = 35 degC
%
% Features:
% - Structured mesh generation
% - 8-node quadrilateral finite elements
% - Internal heat generation
% - Convective boundary conditions
% - Secant method root finding
%
% Required file:
%   fem_symbolic.mat
%
% =========================================================================

clear
close all
clc

%% Initial values for the secant method

Q = [];
f_u = [];

Q(1) = 800;
f_u(1) = -12.96;

Q(2) = 800000;
f_u(2) = 25;

iteration = 3;

%% Secant method loop

while f_u(iteration-1) < -0.01 || f_u(iteration-1) > 0
    Q(iteration)=Q(iteration-1)-((f_u(iteration-1)*(Q(iteration-1)-Q(iteration-2)))/(f_u(iteration-1)-f_u(iteration-2)));

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

    %% Mesh generation

    elementSize = 0.25;

    [xGrid, yGrid] = meshgrid( ...
        min(vertices(:,1)):elementSize:max(vertices(:,1)), ...
        min(vertices(:,2)):elementSize:max(vertices(:,2)));

    insideDomain = inpolygon( ...
        xGrid, yGrid, ...
        vertices(:,1), vertices(:,2));

    xInside = xGrid(insideDomain);
    yInside = yGrid(insideDomain);

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

    %% Load symbolic FEM functions

    load fem_symbolic

    %% Global matrix initialization

    numberOfNodes = size(nodes,1);
    numberOfElements = size(ELEMENTS,1);

    K = zeros(numberOfNodes);
    R = zeros(numberOfNodes,1);

    %% Material and thermal parameters

    kx = 55;
    ky = 55;

    p = 0;

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

            q_e = Q(iteration);

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

            if edge_fully_on_poly( ...
                    p1, p2, BC_poly, tolerance)

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

    %% Objective function evaluation

    f_u(iteration) = max(u) - 35 - 273.15;

    iteration = iteration + 1;
end

%% Final result

maximumHeatSource = Q(iteration-1)

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