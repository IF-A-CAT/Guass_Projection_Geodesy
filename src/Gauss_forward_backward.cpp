#include"Gauss_forward_backward.h"
#include<cmath>
void Gauss_forward_calculate(double Data2[3],double B1,double L1,double S,double A12)
{
    const double a=6378137;
    const double b=6356752.3142;
    const double f=1/298.257223563;
    const double e=sqrt((a*a-b*b)/(a*a));
    const double e1=e/sqrt(1-e*e);
    const double c=a*sqrt(1+e1*e1);
    double Bm,Am,V,N,F10,F20,F30,t;
    double B20=100,L20=100,A20=100;
    double B21=B1,L21=L1,A21=(A12+M_PI)>2*M_PI?(A12-M_PI):(A12+M_PI);
    t=tan(B1);
    V=sqrt(1+e1*e1*cos(B1)*cos(B1));
    N=c/V;
    B21=B1+V*V/N*S*cos(A12);
    L21=L1+sin(A12)/(N*cos(B1))*S;
    A21=A12+S/N*sin(A12)*t;
    while((abs(B21-B20)>1e-10||abs(L21-L20)>1e-10||abs(A21-A20))>1e-10)
    {
        B20=B21;
        L20=L21;
        A20=A21;
        Bm=(B1+B20)/2;
        Am=(A12+A20)/2;
        V=sqrt(1+e1*e1*cos(Bm)*cos(Bm));
        N=c/V;
        t=tan(Bm);
        F10=V*V/N*S*cos(Am)+S*S*S/(24*N*N*N)*cos(Am)*(sin(Am)*sin(Am)*V*V*(2+2*e1*e1*cos(Bm)*cos(Bm)+3*t*t
        )-3*e1*e1*cos(Bm)*cos(Bm)*cos(Am)*cos(Am)*V*V*(1+e1*e1*cos(Bm)*cos(Bm)-t*t+4*e1*e1*cos(Bm)*cos(Bm)*t*t));
        //  V^2/N*S*cos(A) + S^3/(24*N^3)*cos(A)*( sin(A)^2*V^2*(2+ 2*eta^2 + 3*t^2) -3*eta^2*cos(A)^2*V^2*(  1 + eta^2  - t^2 + 4*eta^2*t^2 )  ) ;
        F20=sin(Am)/(N*cos(Bm))*S+S*S*S/(24*N*N*N*cos(Bm))*sin(Am)*( t*t*sin(Am)*sin(Am)-
        cos(Am)*cos(Am)*(1+e1*e1*cos(Bm)*cos(Bm)-9*e1*e1*cos(Bm)*cos(Bm)*t*t));
        F30=S/N*sin(Am)*t+S*S*S/(24*N*N*N)*t*sin(Am)*(cos(Am)*cos(Am)*(2+7*e1*e1*cos(Bm)*cos(Bm)+
        5*e1*e1*cos(Bm)*cos(Bm)*e1*e1*cos(Bm)*cos(Bm)+9*e1*e1*cos(Bm)*cos(Bm)*t*t)+sin(Am)*sin(Am)*(2+2*e1*e1*cos(Bm)*cos(Bm)+t*t));
        B21=F10+B1;
        L21=F20+L1;
        A21=F30+A12;

    }
    Data2[0]=B21;
    Data2[1]=L21;
    Data2[2]=A21;
}

void Gauss_backward_calculate(double Data[3],double B1,double B2,double L1,double L2)
{
    const double a=6378137;
    const double b=6356752.314;
    const double f=1/298.257223563;
    const double e=sqrt((a*a-b*b)/(a*a));
    const double e1=e/sqrt(1-e*e);
    const double c=a*sqrt(1+e1*e1);
    double Bm,t,V,N,eta,q1,q2,p1,p2,deltaB,deltaL,A,SS,
    SsinA,ScosA,SsinA_new,ScosA_new,A1,A2,dA;
    Bm=(B1+B2)/2;
    t=tan(Bm);
    eta=e1*cos(Bm);
    V=sqrt(1+eta*eta);
    N=c/V;
    q1=1/(24*N*N)*t*t;
    q2=1/(24*N*N)*(1+eta*eta-9*eta*eta*t*t);
    p1=1/(24*N*N)*(2+2*eta*eta+3*t*t);
    p2=1/(24*N*N)*3*eta*eta*(1+eta*eta-t*t+4*eta*eta*t*t);
    deltaB=B2-B1;
    deltaL=L2-L1;
    SsinA=100;ScosA=100;
    SsinA_new=deltaL*(N*cos(Bm));
    ScosA_new=deltaB*(N/(V*V));
    while(abs(SsinA-SsinA_new)+abs(ScosA-ScosA_new)>1e-15)
    {
        SsinA=SsinA_new;
        ScosA=ScosA_new;
        SsinA_new=deltaL*N*cos(Bm)-SsinA*SsinA*SsinA*q1+SsinA*ScosA*ScosA*q2;
        ScosA_new=deltaB*N/(V*V)-ScosA*SsinA*SsinA*p1+ScosA*ScosA*ScosA*p2;
    }
    A=atan2(SsinA_new,ScosA_new);
    if(abs(A)<1e-10)
        SS=ScosA_new/cos(A);
    else
        SS=SsinA_new/sin(A);
    Data[0]=SS;
    dA=SS/N*sin(A)*t+SS*SS*SS/(24*N*N*N)*t*sin(A)*(cos(A)*cos(A)*(2+7*eta*eta+5*eta*eta*eta*eta+9*eta*eta*t*t)
    +sin(A)*sin(A)*(2+2*eta*eta+t*t));
    A1=A-dA/2;
    A2=A+dA/2+M_PI;
    A1=A1<0?A1+2*M_PI:A1;
    A1=A1>2*M_PI?A1-2*M_PI:A1;
    A2=A2<0?A2+2*M_PI:A2;
    A2=A2>2*M_PI?A2-2*M_PI:A2;
    Data[1]=A1;
    Data[2]=A2;
}