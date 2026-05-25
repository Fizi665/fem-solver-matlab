clear all
close all
clc

rozmiar =0.25 * 0.01;
% wartrości liczbowe wymagane do nadania okreśonym elementom, brzegom

%kx = 55;
%ky = 55; 
%Tinf = 22;
%h=85;
%Q=800;
%alfa = h;
%beta = h*Tinf;
%q=Q
%p=0

%% Obliczanie symboliczne elementu
syms u u1 u2 u3 u4 u5 u6 u7 u8 a b s t c0 c1 c2 c3 c4 c5 c6 c7 cc kx ky Tinf h Q alfa beta q p
u=[u1 ; u2 ; u3 ; u4; u5;u6; u7; u8];
c=transpose( [c0 c1 c2 c3 c4 c5 c6 c7] );
A=[ ones(8,1) [-a;0;a;a;a;0;-a;-a] [-b;-b;-b;0;b;b;b;0] [a*b;0;-a*b;0; a*b;0;-a*b;0] [a^2;0; a^2; a^2; a^2; 0; a^2; a^2] [b^2; b^2; b^2; 0; b^2; b^2; b^2; 0] [-b*a^2; 0; -b*a^2;0; b*a^2; 0; b*a^2; 0] [-a*b^2; 0; a*b^2;0; a*b^2; 0; -a*b^2; 0]];
N=[1 s t s*t s^2 t^2 s^2*t t^2*s]*inv(A);
% F=factor(N(1))
% c=A^(-1)*u

B=[diff(N,s); diff(N,t)];

C=[kx,0;0,ky];

fk(t,s) = B'*C*B;
kk = int(int(fk,t,[-b,b]),s,[-a,a]); 

fp(t,s) = -p*N'*N;
kp = int(int(fp,t,[-b,b]),s,[-a,a]); 

fq(t,s) = q*N;
rq = int(int(fq,t,[-b,b]),s,[-a,a]);

%brzeg 1-3
Nc13=subs(N,[s,t],[cc,-b]);
falfa13 =-alfa*Nc13.'*Nc13;
kalfa13 = simplify(int( falfa13,cc,[-a,a]));
fbeta13= beta*Nc13; %*b ??l
rbeta13 =simplify( int( fbeta13,cc,[-a,a]));

%brzeg 3-5
Nc35=subs(N,[s,t],[a,cc]);
falfa35=-alfa*Nc35'*Nc35;
kalfa35 = int( falfa35,cc,[-b,b]);
fbeta35= beta*Nc35;
rbeta35 = int( fbeta35,cc,[-b,b]);

%brzeg 5-7
Nc57=subs(N,[s,t],[-cc,b]);
falfa57=-alfa*Nc57'*Nc57;
kalfa57 = int( falfa57,cc,[-a,a]);
fbeta57= beta*Nc57;
rbeta57 = int( fbeta57,cc,[-a,a]);

%brzeg 7-1
Nc71=subs(N,[s,t],[-a,-cc]);
falfa71=-alfa*Nc71'*Nc71;
kalfa71 = int( falfa71,cc,[-b,b]);
fbeta71= beta*Nc71;
rbeta71 = int( fbeta71,cc,[-b,b]);

%u=(kk+kp+kalfa)^(-1)*(rq+rbeta)
%% przejscie z symboliki na funkcje
a_val = rozmiar/2;
b_val = rozmiar/2;
Kk_sym = matlabFunction(subs(kk,{a,b},{a_val,b_val}), 'Vars',[kx ky]);
Kp_sym = matlabFunction(subs(kp,{a,b},{a_val,b_val}), 'Vars',p);
Rq_sym = matlabFunction(subs(rq,{a,b},{a_val,b_val}), 'Vars',q);

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

save fem_symbolic Kk_sym Kp_sym Rq_sym Kalfa_edge Rbeta_edge
