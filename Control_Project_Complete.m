 clear;close all;clc

%% Constants

UA = 1200;
RHO = 1200;
CP = .42;
RHOJ = 1000;
CPJ = 4.18;
K1 = 6e5;
K2 = 2e6;
E1 = 5e4;
E2 = 5.5e4;
CAO = 7;
T0 = 300;
F = 4;
FJ = 1;
TJ0 = 240;
DH1 = -2e4;
DH2 = -3e4;
V = 6;
VJ = 5;
R = 8.314;

%% Steady-State Guesses
% 
% CAS = 0.0226301208683143;
% CBS = 0.0150450773785251;
% CCS = 6.96232480175316;
% TS = 753.722140349192;
% TJS = 354.584864018407;

CAS = 6.99677782598635;
CBS = 0.00322159294563358;
CCS = 5.81073225911324e-07;
TS = 281.113901760693;
TJS = 249.17038700984;

%% Vars for Matrix

Av = F/V; %Fo/V
Bv = K1*exp(-E1/(R*TS));
Cv = K2*exp(-E2/(R*TS));
Dv = E1/(R*TS^2);
Ev = E2/(R*TS^2);
Fv = UA/(RHO*CP*V);
Gv = UA/(RHOJ*CPJ*VJ);
Hv = FJ/VJ;
Iv = DH1/(RHO*CP);
Jv = DH2/(RHO*CP);

%%  PART 1 Matrix vectors, Matrix and eigen vector
M1 = [(-Av-Bv),0,0,(-Bv*CAS*Dv),0];
M2 = [(Bv),(-Av-Cv),0,(Bv*Dv*CAS-Cv*Ev*CBS),0];
M3 = [0,Cv,-Av,(Cv*Ev*CBS),0];
M4 = [(-Bv*Iv),(-Cv*Jv),0,(-Fv-Av-Bv*Dv*CAS*Iv-Cv*Ev*CBS*Jv),Fv];
M5 = [0,0,0,Gv,(-Hv-Gv)];

M = [M1;M2;M3;M4;M5];

J = eig(M);

%% Part 3 to find Kc vals and test potential issue values. TEST matrix will have 14 columns
% and 7 rows each column pair corresponds to a Kc value and each row, one
% of the test equations. odd columns are KC < val even are KC > Val. So I
% have to find a column of all positive values.
a = M1(1,1);
b = M1(1,4);
c = M2(1,1);
d = M2(1,2);
e = M2(1,4);
i = M4(1,1);
j = M4(1,2);
k = M4(1,4);
l = M4(1,5);
m = Av;
n = M5(1,5);
o = M5(1,4);
p = Hv;

A = d + a + n + k;
B = n*d+a*n+a*d-b*i-j*e+k*d+k*a+k*n-l*o;
C = (a*n*d-b*i*n-b*i*d+j*b*c-j*e*n-j*e*a+k*n*d+k*a*n+k*a*d-l*o*a-l*o*d);
D = -b*i*n*d + j*b*c*n - j*e*a*n + k*a*n*d - l*o*a*d;

Kc1 = -B/(l*p)
Kc2 = -C/(l*p*a+l*p*d) 
Kc3 = -D/(l*p*a*d)

x0 = [-1e10 1e10];
syms x1
fun4 = @(x) ((-A*(B+x*l*p)-(-x*l*p*a-x*l*p*d-C))/-A);
fun5 = (((-A*(B+x1*l*p)-(-x1*l*p*a-x1*l*p*d-C))/-A).*(-x1*l*p*a-x1*l*p*d-C)-(-A)*(D+x1*l*p*a*d))./((-A*(B+x1*l*p)-(-x1*l*p*a-x1*l*p*d-C))/-A) == 0;
fun6 = @(x) (D+x*l*p*a*d);
Kc4 = fzero(fun4,x0)
Kc55 = solve(fun5,x1);
Kc5 = double(Kc55)'
Kc6 = fzero(fun6,x0)

KCT = [Kc1-.1,Kc1+.1;Kc2-.1,Kc2+.1;Kc3-.1,Kc3+.1;Kc4-.1,Kc4+.1;Kc5(1,1)-.001,Kc5(1,1)+.001;Kc5(1,2)-.01,Kc5(1,2)+.01;Kc6-.1,Kc6+.1];
KCTT = [Kc1;Kc2;Kc3;Kc4;Kc5(1);Kc5(2);Kc6];

for jj = 1:7
    for ii = 1:2
        TEST(1,ii+2*(jj-1)) = KCT(jj,ii)*l*p +B;
        TEST(2,ii+2*(jj-1)) = KCT(jj,ii)*(-l*p*a-l*p*d)-C;
        TEST(3,ii+2*(jj-1)) = KCT(jj,ii)*(l*p*a*d)+D;
        TEST(4,ii+2*(jj-1)) = ((-A*(B+KCT(jj,ii)*l*p)-(-KCT(jj,ii)*l*p*a-KCT(jj,ii)*l*p*d-C))/-A);
        TEST(5,ii+2*(jj-1)) = (((-A*(B+KCT(jj,ii)*l*p)-(-KCT(jj,ii)*l*p*a-KCT(jj,ii)*l*p*d-C))/-A)*(-KCT(jj,ii)*l*p*a-KCT(jj,ii)*l*p*d-C)-(-A)*(D+KCT(jj,ii)*l*p*a*d))/((-A*(B+KCT(jj,ii)*l*p)-(-KCT(jj,ii)*l*p*a-KCT(jj,ii)*l*p*d-C))/-A);
        TEST(6,ii+2*(jj-1)) = (((-A*(B+KCT(jj,ii)*l*p)-(-KCT(jj,ii)*l*p*a-KCT(jj,ii)*l*p*d-C))/-A)*(-KCT(jj,ii)*l*p*a-KCT(jj,ii)*l*p*d-C)-(-A)*(D+KCT(jj,ii)*l*p*a*d))/((-A*(B+KCT(jj,ii)*l*p)-(-KCT(jj,ii)*l*p*a-KCT(jj,ii)*l*p*d-C))/-A);
        TEST(7,ii+2*(jj-1)) = KCT(jj,ii)*(l*p*a*d)+D;
    end
end


%% Part 4/5/6 design the PI controller aspect by making ubar Kc*e(s)*(1+1/taui*s
time = 0;
h = 1e-5;
finaltime = 50;

INTSTEP = floor(finaltime/h +.5);

YMATL = zeros(INTSTEP,1);
ZMATL = zeros(INTSTEP,1);
CAMATL = zeros(INTSTEP,1);
CBMATL = zeros(INTSTEP,1);
TIMEMAT = zeros(INTSTEP,1);
TJMATL = zeros(INTSTEP,1);
UMATL = zeros(INTSTEP,1);

TMATNL = zeros(INTSTEP,1);
ZMATNL = zeros(INTSTEP,1);
CAMATNL = zeros(INTSTEP,1);
CBMATNL = zeros(INTSTEP,1);
TJMATNL = zeros(INTSTEP,1);
UMATNL = zeros(INTSTEP,1);

%arbitrarily changed steady state values using random addition

TMATNL(1,1) = TS;
CAMATNL(1,1) = CAS;
CBMATNL(1,1) = CBS;
TJMATNL(1,1) = TJS;


%guessed Kc value between range and guess tauI value to find stability.
KCL = 6.35;
TAUIL = 4.8;
KCNL = 6.35;
TAUINL = 4.8;
SPC = 1;
YSPNL = TS + SPC;
YSPL = SPC;
DISTUSC = 0;
DISTNL = T0+DISTUSC;
DISTL = DISTUSC;
UMATL(1,1) = KCL*(YSPL-YMATL(1,1));
UMATNL(1,1) = KCL*(YSPNL-TMATNL(1,1))+TJ0;
for ii = 1:(INTSTEP-1)
    
   
    %LINEARIZED SYSTEM
    UMATL(ii+1,1) = KCL*(YSPL-YMATL(ii,1)) + (KCL/TAUIL)*ZMATL(ii,1);
    CAMATL(ii+1,1) = CAMATL(ii,1) + h*(a*CAMATL(ii,1) + b*YMATL(ii,1));
    CBMATL(ii+1,1) = CBMATL(ii,1) + h*(c*CAMATL(ii,1) + d*CBMATL(ii,1) + e*YMATL(ii,1));
    YMATL(ii+1,1) = YMATL(ii,1) + h*(i*CAMATL(ii,1) + j*CBMATL(ii,1) + k*YMATL(ii,1) + l*TJMATL(ii,1) + m*DISTL);
    TJMATL(ii+1,1) = TJMATL(ii,1) + h*(n*TJMATL(ii,1) + o*YMATL(ii,1) + p*UMATL(ii,1));
    ZMATL(ii+1,1) = ZMATL(ii,1) + h*(YSPL-YMATL(ii,1));
    
    %nonlinear system
    UMATNL(ii+1,1) = KCNL*(YSPNL-TMATNL(ii,1)) + (KCNL/TAUINL)*ZMATNL(ii,1) + TJ0;
    CAMATNL(ii+1,1) = CAMATNL(ii,1) + h*(Av*(CAO-CAMATNL(ii,1)) - CAMATNL(ii,1)*K1*exp(-E1/(R*TMATNL(ii,1))));
    CBMATNL(ii+1,1) = CBMATNL(ii,1) + h*(CBMATNL(ii,1)*-Av + CAMATNL(ii,1)*K1*exp(-E1/(R*TMATNL(ii,1))) -CBMATNL(ii,1)*K2*exp(-E2/(R*TMATNL(ii,1))));
    TMATNL(ii+1,1) = TMATNL(ii,1)+ h*(Fv*(TJMATNL(ii,1) - TMATNL(ii,1)) + Av*(DISTNL-TMATNL(ii,1)) - CAMATNL(ii,1)*K1*exp(-E1/(R*TMATNL(ii,1)))*Iv - CBMATNL(ii,1)*K2*exp(-E2/(R*TMATNL(ii,1)))*Jv);
    TJMATNL(ii+1,1) = TJMATNL(ii,1) + h*(Hv*(UMATNL(ii,1)-TJMATNL(ii,1)) +  Gv*(TMATNL(ii,1)-TJMATNL(ii,1)));
    ZMATNL(ii+1,1) = ZMATNL(ii,1) + h*(YSPNL-TMATNL(ii,1));
    
    TIMEMAT(ii+1,1) = TIMEMAT(ii,1) + h;
    
end

PLOTVEC = ones(1,length(TIMEMAT));
plot(TIMEMAT,YMATL+TS)
hold on
plot(TIMEMAT,PLOTVEC*YSPNL,'red')
plot(TIMEMAT,TMATNL)
title('Linear Output')
xlabel('Time (h)')
ylabel('Temperature of Reactor (K)')
legend('Temperature','Set Point','Location','Southeast')
hold off

% subplot(1,2,1)
% PLOTVEC = ones(1,length(TIMEMAT));
% plot(TIMEMAT,YMATL+TS)
% hold on
% plot(TIMEMAT,PLOTVEC*YSPNL,'red')
% %title('Linear Output')
% xlabel('Time (h)')
% ylabel('Temperature of Reactor (K)')
% legend('Temperature','Set Point','Location','Southeast')
% hold off
% 
% subplot(1,2,2)
% plot(TIMEMAT,TMATNL)
% hold on
% plot(TIMEMAT,PLOTVEC*YSPNL,'red')
% %title('Non-linear')
% xlabel('Time (h)')
% ylabel('Temperature of Reactor (K)')
% legend('Temperature','Set Point','Location','Southeast')
% 
% figure(2)
% subplot(1,2,1)
% plot(TIMEMAT,UMATL+TJ0)
% %title('Linear Input')
% xlabel('Time (h)')
% ylabel('Temperature of Jacket Inlet (K)')
% 
% subplot(1,2,2)
% plot(TIMEMAT,UMATNL)
% %title('Non-linear Input')
% xlabel('Time (h)')
% ylabel('Temperature of Jacket Inlet (K)')
% 
% 
