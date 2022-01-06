#include"Gauss_projection_fbward.h"
#include<cmath>
void Guass_f_projection(double B,double L,double *X,double *Y)
{
    // const double a=6378137;
    // const double b=6356752.3142;
    // const double f=1/298.257223563;
    // const double e=sqrt((a*a-b*b)/(a*a));
    // const double e1=e/sqrt(1-e*e);
    // const double c=a*sqrt(1+e1*e1);
    double Brad=B/180*M_PI;
    double Lrad=L/180*M_PI;
    double N,x,y,L0,l;
    double a_coef[5];
    L0=int(L/6)*6+3;
    l=Lrad-L0/180*M_PI;
    // double m_coef[5];
    // V=sqrt(1+e1*e1*cos(Brad)*cos(Brad));
    N=6399596.652-(21565.045-(108.996-0.603*cos(Brad)*cos(Brad))*cos(Brad)*cos(Brad))*cos(Brad)*cos(Brad);
    // X_C=111133.005*B-16038.528*sin(2*Brad)+16.833*sin(4*Brad)-0.022*sin(6*Brad);
    // m_coef[0]=a*(1-e*e);
    // m_coef[1]=3/2*e*e*m_coef[0];
    // m_coef[2]=5/4*e*e*m_coef[1];
    // m_coef[3]=7/6*e*e*m_coef[2];
    // m_coef[4]=9/8*e*e*m_coef[3];
    a_coef[0]=32144.5189-(135.3646-(0.7034-0.0041*cos(Brad)*cos(Brad))*cos(Brad)*cos(Brad))*cos(Brad)*cos(Brad);
    a_coef[1]=(0.25+0.00253*cos(Brad)*cos(Brad))*cos(Brad)*cos(Brad)-0.04167;
    a_coef[2]=(0.167*cos(Brad)*cos(Brad)-0.083)-0.04167;
    a_coef[3]=(0.3333333+0.001123*cos(Brad)*cos(Brad))*cos(Brad)*cos(Brad)-0.1666667;;
    a_coef[4]=0.00878-(0.1702-0.20382*cos(Brad)*cos(Brad))*cos(Brad)*cos(Brad);
    // X_C=a_coef[0]*Brad-a_coef[1]/2*sin(2*Brad)+a_coef[2]/4*sin(4*Brad)-a_coef[3]/6*sin(6*Brad)+a_coef[4]/8*sin(8*Brad);
    x=6367452.1328*Brad-(a_coef[0]-(0.5+(a_coef[1]+a_coef[2]*l*l)*l*l)*l*l*N)*cos(Brad)*sin(Brad);
    y=(1+(a_coef[3]+a_coef[4]*l*l)*l*l)*l*N*cos(Brad);
    *X=x;*Y=y;
}

void Gauss_b_projection(double X,double Y,double *B,double *L)
{
    double Beta,Bf;
    double N,Z,b_coef[4];
    Beta=X/6367452.133;
    Bf=Beta+(50228976+(293697+(2383+22*cos(Beta)*cos(Beta))*cos(Beta)*cos(Beta))*cos(Beta)*cos(Beta))*(1e-10)*sin(Beta)*cos(Beta);
    N=6399596.652-(21565.045-(108.996-0.603*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf);
    Z=Y/N/cos(Bf);
    b_coef[0]=(0.5+0.00336975*cos(Bf)*cos(Bf))*sin(Bf)*cos(Bf);
    b_coef[1]=0.333333-(0.1666667-0.001123*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf);
    b_coef[2]=0.25+(0.161612+0.005617*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf);
    b_coef[3]=0.2-(0.16667-0.00878*cos(Bf)*cos(Bf))*cos(Bf)*cos(Bf);
    *B=Bf-(1-(b_coef[2]-0.147*Z*Z)*Z*Z)*Z*Z*b_coef[0];
    *L=(1-(b_coef[1]-b_coef[3]*Z*Z)*Z*Z)*Z;
}