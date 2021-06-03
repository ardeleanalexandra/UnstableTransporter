clear all
clc
A=[0 1 0 0 0 0;
    0 -0.4733 -29.1464 0.4733 0 0;
    0 0 0 1 0 0;
    0 0.1401 29.1651 -0.1401 0 0;
    0 0 0 0 0 1;
    0 0 0 0 0 -2.6349];
B=[0 0;
    109.3664 109.3664;
    0 0;
    -32.3752 -32.3752;
    0 0;
    -0.5055 0.5055];

C=eye(6);
D=zeros(6,2);
sys=ss(A,B,C,D);
Hp=tf(sys);

G11=Hp(1,1);
G12=Hp(2,1);
G13=Hp(3,1);
G14=Hp(4,1);
G15=Hp(5,1);
G16=Hp(6,1);
G21=Hp(1,2);
G22=Hp(2,2);
G23=Hp(3,2);
G24=Hp(4,2);
G25=Hp(5,2);
G26=Hp(6,2);

G_11=minreal(zpk(G11))
G_12=minreal(zpk(G13))
G_21=minreal(zpk(G21))
G_22=minreal(zpk(G23))

[num11,den11]=tfdata(G_11,'V');
[num12,den12]=tfdata(G_12,'V');
[num21,den21]=tfdata(G_21,'V');
[num22,den22]=tfdata(G_22,'V');

G=[G_11 G_12;
   G_21 G_22];

G_inv=inv(G)

s=tf('s');
H11=minreal(zpk(4.7003e14/(s+5.552)/(s+0.3322)))
H22=minreal(zpk(1.5878e15/(s+4.532)/(s+5.552)))

[num_h11,den_h11]=tfdata(H11,'v');
[num_h22,den_h22]=tfdata(H22,'v');

H=[H11 0;
    0 H22];
%decuplor
D=minreal(zpk(G_inv*H))

D11=D(1,1)
D22=D(2,2)
D12=D(1,2)
D21=D(2,1)

[num_d12,den_d12]=tfdata(D12,'V');
[num_d21,den_d21]=tfdata(D21,'V');
[num_d11,den_d11]=tfdata(D11,'V');
[num_d22,den_d22]=tfdata(D22,'V');

%partile fixe pt controller
G1=minreal(zpk(H11))
G2=minreal(zpk(H22))
%% Calcul controller 1
sigma=0.1;
tr=2;

% =>
zeta=(-log(sigma))/(sqrt(pi^2+(log(sigma))^2));
wn=4/(zeta*tr);
cv=wn/(2*zeta);
wb=wn*sqrt(1-2*zeta^2+sqrt(2-4*zeta^2+4*zeta^4));

G01=tf(wn^2,[1 2*zeta*wn wn^2])
Gc1=minreal(zpk(1/G1*G01/(1-G01)))


[num_gc1,den_gc1]=tfdata(Gc1,'V')

Gs1=series(Gc1,H11)
Gc01=feedback(Gs1,1)
step(Gc01)
%% Calcul controller 2
sigma=0.05;
tr=1;

% =>
zeta=(-log(sigma))/(sqrt(pi^2+(log(sigma))^2));
wn=4/(zeta*tr);
cv=wn/(2*zeta);
wb=wn*sqrt(1-2*zeta^2+sqrt(2-4*zeta^2+4*zeta^4));

G02=tf(wn^2,[1 2*zeta*wn wn^2])

Gc2=minreal(zpk(1/G2*G02/(1-G02)))

[num_gc2,den_gc2]=tfdata(Gc2,'V');

Gs2=series(Gc2,H22)
Gc02=feedback(Gs2,1)
step(Gc02)