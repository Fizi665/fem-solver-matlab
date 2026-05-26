% =========================================================================
% generate_symbolic_element.m
% =========================================================================
% Author: Wiktor Komorowski
%
% Description:
% Symbolic derivation of an 8-node finite element for 2D steady-state
% heat transfer analysis.
%
% The script generates:
% - Conductivity matrix
% - Reaction matrix
% - Heat source vector
% - Convection boundary matrices
%
% Symbolic expressions are converted into MATLAB functions and saved
% into the file:
%   fem_symbolic.mat
%
% =========================================================================

clear
close all
clc

%% Element size definition

elementSize = 0.25 * 0.01;

%% Symbolic variable definitions

syms u1 u2 u3 u4 u5 u6 u7 u8
syms a b s t cc
syms c0 c1 c2 c3 c4 c5 c6 c7
syms kx ky Tinf h Q alfa beta q p

u = [u1; u2; u3; u4; u5; u6; u7; u8];

%% Shape function matrix

A = [ ...
    ones(8,1), ...
    [-a;0;a;a;a;0;-a;-a], ...
    [-b;-b;-b;0;b;b;b;0], ...
    [a*b;0;-a*b;0;a*b;0;-a*b;0], ...
    [a^2;0;a^2;a^2;a^2;0;a^2;a^2], ...
    [b^2;b^2;b^2;0;b^2;b^2;b^2;0], ...
    [-b*a^2;0;-b*a^2;0;b*a^2;0;b*a^2;0], ...
    [-a*b^2;0;a*b^2;0;a*b^2;0;-a*b^2;0] ...
];

N = [1 s t s*t s^2 t^2 s^2*t t^2*s] * inv(A);

%% Derivative matrix

B = [
    diff(N,s)
    diff(N,t)
];

%% Material property matrix

C = [
    kx 0
    0 ky
];

%% Conductivity matrix

fk(t,s) = B' * C * B;

kk = int( ...
    int(fk, t, [-b,b]), ...
    s, [-a,a]);

%% Reaction matrix

fp(t,s) = -p * N' * N;

kp = int( ...
    int(fp, t, [-b,b]), ...
    s, [-a,a]);

%% Heat source vector

fq(t,s) = q * N;

rq = int( ...
    int(fq, t, [-b,b]), ...
    s, [-a,a]);

%% Boundary 1-3

Nc13 = subs(N,[s,t],[cc,-b]);

falfa13 = -alfa * Nc13.' * Nc13;
kalfa13 = simplify(int(falfa13, cc, [-a,a]));

fbeta13 = beta * Nc13;
rbeta13 = simplify(int(fbeta13, cc, [-a,a]));

%% Boundary 3-5

Nc35 = subs(N,[s,t],[a,cc]);

falfa35 = -alfa * Nc35' * Nc35;
kalfa35 = int(falfa35, cc, [-b,b]);

fbeta35 = beta * Nc35;
rbeta35 = int(fbeta35, cc, [-b,b]);

%% Boundary 5-7

Nc57 = subs(N,[s,t],[-cc,b]);

falfa57 = -alfa * Nc57' * Nc57;
kalfa57 = int(falfa57, cc, [-a,a]);

fbeta57 = beta * Nc57;
rbeta57 = int(fbeta57, cc, [-a,a]);

%% Boundary 7-1

Nc71 = subs(N,[s,t],[-a,-cc]);

falfa71 = -alfa * Nc71' * Nc71;
kalfa71 = int(falfa71, cc, [-b,b]);

fbeta71 = beta * Nc71;
rbeta71 = int(fbeta71, cc, [-b,b]);

%% Convert symbolic expressions into MATLAB functions

a_val = elementSize / 2;
b_val = elementSize / 2;

Kk_sym = matlabFunction( ...
    subs(kk,{a,b},{a_val,b_val}), ...
    'Vars',[kx ky]);

Kp_sym = matlabFunction( ...
    subs(kp,{a,b},{a_val,b_val}), ...
    'Vars',p);

Rq_sym = matlabFunction( ...
    subs(rq,{a,b},{a_val,b_val}), ...
    'Vars',q);

%% Boundary convection matrices

Kalfa_edge = {
    matlabFunction(subs(-kalfa13,{a,b},{a_val,b_val}),'Vars',alfa)
    matlabFunction(subs(-kalfa35,{a,b},{a_val,b_val}),'Vars',alfa)
    matlabFunction(subs(-kalfa57,{a,b},{a_val,b_val}),'Vars',alfa)
    matlabFunction(subs(-kalfa71,{a,b},{a_val,b_val}),'Vars',alfa)
};

Rbeta_edge = {
    matlabFunction(subs(rbeta13,{a,b},{a_val,b_val}),'Vars',beta)
    matlabFunction(subs(rbeta35,{a,b},{a_val,b_val}),'Vars',beta)
    matlabFunction(subs(rbeta57,{a,b},{a_val,b_val}),'Vars',beta)
    matlabFunction(subs(rbeta71,{a,b},{a_val,b_val}),'Vars',beta)
};

%% Save symbolic FEM functions

save fem_symbolic Kk_sym Kp_sym Rq_sym Kalfa_edge Rbeta_edge