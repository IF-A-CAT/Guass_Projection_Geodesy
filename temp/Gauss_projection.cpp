#include<iostream>
#include<iomanip>
#include<cmath>
#include"Gauss_projection_fbward.h"
using namespace std;

int main()
{
    double B=40,L=120.5;
    int n=int(L/6)+1;
    double X,Y;
    Guass_f_projection(B,L,&X,&Y);
    cout<<"-----------------正算结果（单位：m）---------------"<<endl;
    cout<<"X: "<<setprecision(8)<<X<<endl;
    cout<<"Y: "<<setprecision(8)<<Y<<endl;
    Gauss_b_projection(X,Y,&B,&L);
    cout<<"-----------------反算结果（单位：°）---------------"<<endl;
    cout<<"B: "<<setprecision(8)<<B*180/M_PI<<endl;
    cout<<"l: "<<setprecision(8)<<L*180/M_PI<<endl;
    cout<<"L: "<<setprecision(8)<<n*6-3+L*180/M_PI<<endl;
    system("pause");
    return 0;
}