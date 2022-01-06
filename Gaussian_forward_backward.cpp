#include<iostream>
#include"Gauss_forward_backward.h"
#include<iomanip>
using namespace std;
#include<cmath>
int main()
{
    const double a=6378137;
    const double b=6356752.3142;
    const double f=1/298.257223563;
    double B1=47.4652647/180*M_PI;
    double L1=35.493633/180*M_PI;
    double A12=44.1213664/180*M_PI;
    double S=44797.2826;
    double L11,B11,L22,B22;
    double Data2[3],Data1[3];
    Gauss_forward_calculate(Data1,B1,L1,S,A12);
    B11=B1;L11=L1;B22=Data1[0];L22=Data1[1];
    Data1[2]=Data1[2]>M_PI?Data1[2]-M_PI:Data1[2]+M_PI;
    cout<<"-----------------正算结果(单位：°)------------------"<<endl;
    cout<<"B2: "<<setprecision(8)<<Data1[0]*180/M_PI<<endl;
    cout<<"L2: "<<setprecision(8)<<Data1[1]*180/M_PI<<endl;
    cout<<"A21: "<<setprecision(8)<<Data1[2]*180/M_PI<<endl;
    Gauss_backward_calculate(Data2,B11,B22,L11,L22);
    cout<<"----------------反算结果（单位：°）-----------------"<<endl;
    cout<<"S: "<<setprecision(8)<<Data2[0]<<endl;
    cout<<"A12: "<<setprecision(8)<< Data2[1]*180/M_PI<<endl;
    cout<<"A21: "<<setprecision(8)<<Data2[2]*180/M_PI<<endl;
    // double S1=Data2[2];
    system("pause");
    return 0;

}