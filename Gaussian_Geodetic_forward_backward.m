function Gaussian_Geodetic_forward_backward()
%coded on April 16, 2020 by jjwang
clc;clear;
nselect=1;
switch nselect
    case 1
    double B1=39.6224;
    double L1=112.9263;
    double A1=61.2926;
    double S=100.4335;
    case 2
B1=30.29582042;
L1=120.0540218;
A1=247.2750428;
S=28230.93486;%meter
end
n=1;
[B2,L2,A2]=Gaussian_forward(B1,L1,A1,S,n);
[A12,A21,S12]=Gaussian_backward(B1,L1,B2,L2,n);
%-------------------------------------------------------------------------%
B1dms=degree2degreeminutesecond(B1);
L1dms=degree2degreeminutesecond(L1);
A1dms=degree2degreeminutesecond(A1);
%
B2dms=degree2degreeminutesecond(B2);
L2dms=degree2degreeminutesecond(L2);
A2dms=degree2degreeminutesecond(A2);
%
A12dms=degree2degreeminutesecond(A12);
A21dms=degree2degreeminutesecond(A21);
%
N=size(B1dms,1);
for i=1:N
input_D=[sprintf('B1=%-d°%-d''%-f''',B1dms(i,1),B1dms(i,2),B1dms(i,3) ) '   '...
 sprintf('L1=%-d°%-d''%-f''',L1dms(i,1),L1dms(i,2),L1dms(i,3) ) '   '...
 sprintf('A1=%-d°%-d''%-f''',A1dms(i,1),A1dms(i,2),A1dms(i,3) ) '   ' ...
 sprintf('S=%13.6f m',S(i))];
forward_D=[sprintf('B2=%-d°%-d''%-f''',B2dms(i,1),B2dms(i,2),B2dms(i,3) ) '   '...
 sprintf('L2=%-d°%-d''%-f''',L2dms(i,1),L2dms(i,2),L2dms(i,3) ) '   '...
 sprintf('A2=%-d°%-d''%-f''',A2dms(i,1),A2dms(i,2),A2dms(i,3) ) ];
backward_D=[sprintf('A12=%-d°%-d''%-f''',A12dms(i,1),A12dms(i,2),A12dms(i,3) ) '   '...
 sprintf('A21=%-d°%-d''%-f''',A21dms(i,1),A21dms(i,2),A21dms(i,3) ) '   '...
 sprintf('S=%13.6f m',S12(i))];
 input_D
 forward_D
 backward_D
 disp('------------------------------------------------------------------');
end
% %along the parallel circle
% B1=30;
% L1=10;
% B2=30;
% L2=40;
% n=1;
% [A1,A2,S]=Gaussian_backward(B1,L1,B2,L2,n)
% disp('------------------------------------------------------------------');
% %along the meridian circle
% B1=30;
% L1=10;
% B2=50;
% L2=10;
% n=1;
% [A1,A2,S]=Gaussian_backward(B1,L1,B2,L2,n)

B1=47+46/60+52.647/3600;
L1=35+49/60+36.330/3600;
B2=48+4/60+09.6384/3600;
L2=36+14/60+45.0505/3600;
n=1;
[A1,A2,S]=Gaussian_backward(B1,L1,B2,L2,n);
num2str(degrees2dms(A1))
num2str(degrees2dms(A2))
num2str(S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B2,L2,A2]=Gaussian_forward(B1,L1,A12,S,n)
switch n
    case 1 
        a=6378245.0;
        b=6356863.0188;
    case 2
        a=6378140.0;
        b=6356775.2882;
    case 3
        a=6378137.0;
        b=6356752.3142;
end
e=sqrt(a^2-b^2)/a;
e1=e/sqrt(1-e^2);
c=a*sqrt(1+e1^2);
NN=length(B1);
B11=deg2rad(B1);
L11=deg2rad(L1);
A11=deg2rad(A12);
B2=[];
L2=[];
A2=[];
for i=1:NN
    B1=B11(i);
    L1=L11(i);
    A1=A11(i);
    A1=A1+(A1<0)*2*pi;
    B20=1.0e+9;
    L20=1.0e+9;
    A20=1.0e+9;
    epsilon=1.0e-10;
    %------------------initial values of dB, dL and dA--------------------%
    B=B1;
    A=A1;
    t=tan(B);
    eta=e1*cos(B);
    V=sqrt(1+eta^2);
    N=c/V;
    dB0=V^2/N*S*cos(A);
    dL0=sin(A)/(N*cos(B))*S;
    dA0=S/N*sin(A)*t;
    %---------------------------------------------------------------------%
    B21=B1+dB0;
    L21=L1+dL0;
    A21=A1+dA0;
    %
    while abs(B21-B20)>epsilon||abs(L21-L20)>epsilon||abs(A21-A20)>epsilon
        B20=B21;
        L20=L21;
        A20=A21;
        B=1/2*(B1+B20);
        A=1/2*(A1+A20);
        t=tan(B);
        eta=e1*cos(B);
        V=sqrt(1+eta^2);
        N=c/V;
        %F10=B2-B1
        F10=  V^2/N*S*cos(A) + S^3/(24*N^3)*cos(A)*( sin(A)^2*V^2*(2+ 2*eta^2 + 3*t^2) -3*eta^2*cos(A)^2*V^2*(  1 + eta^2  - t^2 + 4*eta^2*t^2 )  ) ;
        %F20=L2-L1
        F20=sin(A)/(N*cos(B))*S+S^3/(24*N^3*cos(B))*sin(A)*( t^2*sin(A)^2 - cos(A)^2*(1+eta^2-9*eta^2*t^2) );
        %F30=A21-A12
        F30=S/N*sin(A)*t+S^3/(24*N^3)*t*sin(A)*(     cos(A)^2*(2+7*eta^2+5*eta^4+9*eta^2*t^2 ) +sin(A)^2*( 2+2*eta^2+t^2 )    );
        %
        B21=F10+B1;
        L21=F20+L1;
        A21=F30+A1;
    end
    B2=[B2;B21];
    L2=[L2;L21];
    A2=[A2;A21+pi*(A1<=pi)-pi*(A1>pi)];
end
B2=rad2deg(B2);
L2=rad2deg(L2);
A2=rad2deg(A2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A1,A2,S]=Gaussian_backward(B1,L1,B2,L2,n)
nmethod=2;
switch n
    case 1 
        a=6378245.0;
        b=6356863.0188;
    case 2
        a=6378140.0;
        b=6356775.2882;
    case 3
        a=6378137.0;
        b=6356752.3142;
end
e=sqrt(a^2-b^2)/a;
e1=e/sqrt(1-e^2);
c=a*sqrt(1+e1^2);
B11=deg2rad(B1);
L11=deg2rad(L1);
B22=deg2rad(B2);
L22=deg2rad(L2);
NN=length(B1);
S=[];
A12=[];
A21=[];
epsilon=1.0e-15;
for i=1:NN
            B=1/2*(B11(i)+B22(i));
            t=tan(B);
            eta=e1*cos(B);
            V=sqrt(1+eta^2);
            N=c/V;
            q1=1/(24*N^2)*t^2;
            q2=1/(24*N^2)*( 1 + eta^2 - 9*eta^2*t^2 );
            p1=1/(24*N^2)*( 2 + 3*t^2 + 2*eta^2 );
            p2=1/(24*N^2)*3*eta^2*( 1 + eta^2 - t^2 + 4*eta^2*t^2 );
            deltaB=B22(i)-B11(i);
            deltaL=L22(i)-L11(i);
     switch nmethod
            case 1
               SsinA=deltaL*(N*cos(B));%dL0=sin(A)/(N*cos(B))*S, the first-order term (see 4-221)
               ScosA=deltaB*(N/V^2);%dB0=V^2/N*S*cos(A), also the first-order term (see 4-221)
               SsinA_new=deltaL*N*cos(B)-SsinA^3*q1+SsinA*ScosA^2*q2;
               ScosA_new=deltaB*N/V^2-ScosA*SsinA^2*p1+ScosA^3*p2;
            case 2 %iteration
               SsinA=1.0e+9;
               ScosA=1.0e+0;
               SsinA_new=deltaL*(N*cos(B));%dL0=sin(A)/(N*cos(B))*S, the first-order term (see 4-221)
               ScosA_new=deltaB*(N/V^2);%dB0=V^2/N*S*cos(A), also the first-order term (see 4-221)
            while abs(SsinA-SsinA_new)>epsilon||abs(ScosA-ScosA_new)>epsilon
               SsinA=SsinA_new;
               ScosA=ScosA_new;
               SsinA_new=deltaL*N*cos(B)-SsinA^3*q1+SsinA*ScosA^2*q2;
               ScosA_new=deltaB*N/V^2-ScosA*SsinA^2*p1+ScosA^3*p2; 
            end  
      end
            A=atan2(SsinA_new, ScosA_new);%Am
            if abs(A)<1.0e-8 %be careful! when A=0, the formula S*sin(Am)/sin(Am) breaks down.
                SS=ScosA_new/cos(A);
            else
                SS=SsinA_new/sin(A);
            end
            S=[S;SS];
            %dA=A21-A12
            dA=S/N*sin(A)*t+S^3/(24*N^3)*t*sin(A)*(  cos(A)^2*(2+7*eta^2+5*eta^4+9*eta^2*t^2 ) +sin(A)^2*( 2+2*eta^2+t^2 )  );
            A1=A-1/2*dA;
            A2=A+1/2*dA+pi;
            %
            A2=A2+(A2<0)*2*pi;
            A2=A2-(A2>2*pi)*2*pi;
            A1=A1+(A1<0)*2*pi;
            A1=A1-(A1>2*pi)*2*pi;
            A12=[A12;A1];
            A21=[A21;A2];
            
end
A1=rad2deg(A12);
A2=rad2deg(A21);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=degree2degreeminutesecond(x)
x0=x;
n=length(x0);
for i=1:n
    x=x0(i);
    nflag=0;
if x<0
    x=-x;
    nflag=1;
end
y(i,1)=floor(x);
y(i,2)=floor( (x-y(i,1))*60 );
y(i,3)=( x - y(i,1) - y(i,2)/60 ).*3600;
if nflag==1
    y(i,1)=-y(i,1);
end
end