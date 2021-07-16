#include <iostream>
#include <math.h>
using namespace std;
#define N 10
#define M 2
#define MAX_EQNS (N-2)
double M_[MAX_EQNS][MAX_EQNS],
M_0[MAX_EQNS][MAX_EQNS],
MWM[MAX_EQNS][MAX_EQNS],
x_m[MAX_EQNS],
b_m[MAX_EQNS],
w[MAX_EQNS],
mwb[MAX_EQNS],
LU[MAX_EQNS][MAX_EQNS],
Inv[MAX_EQNS][MAX_EQNS];
int ps[MAX_EQNS];

double T[N];
double q[M];
double dx = 0.5;
double dFdT[N];
double dRdq[N][M];
double dRdT[N][N];
double dFdTfull[N];

void LU_decompose(int size)
{
    int i,j,k,pivotindex;
    double scales[MAX_EQNS];
    double normrow,pivot,size1,biggest,mult;

    for (i=0;i<size;i++) //заполнение начальными данными
    {
        ps[i]=i;//маппинг изначального порядка на переставленный.
        normrow=0;//максимум в итой строке

        for (j=0;j<size;j++) {
            LU[i][j]=M_[i][j];
            if (normrow<fabs(LU[i][j]))
                normrow=fabs(LU[i][j]);
        }
        if (normrow!=0)
            scales[i]=1.0/normrow;//для общих множителей
        else {
            scales[i]=0.0;
            //     err_code(DIV_ZERO);
        }
    }
    //метод гаусса с частичным упорядочиванием

    for (k=0;k<size-1;k++)
    {
        biggest=0;
        for (i=k; i<size;i++) {
            size1=fabs(LU[ps[i]][k])*scales[ps[i]];
            if (biggest<size1) {
                biggest=size1;
                pivotindex=i;
            }
        }

        if (biggest==0) {
            pivotindex=0;
        }

        if (pivotindex!=k) {
            j=ps[k];
            ps[k]=ps[pivotindex];
            ps[pivotindex]=j;
        }

        pivot=LU[ps[k]][k];

        for (i=k+1;i<size;i++) {
            mult=LU[ps[i]][k]/pivot;
            LU[ps[i]][k]=mult;

            if (mult!=0.0) {
                for (j=k+1; j<size;j++)
                    LU[ps[i]][j]-=mult*LU[ps[k]][j];
            }
        }
    }
}

void LU_solve(int num)
{
    int i,j;
    double dot;

    for (i=0;i<num;i++) {
        dot=0;
        for (j=0;j<i;j++)
            dot+=LU[ps[i]][j]*x_m[j];

        x_m[i]=b_m[ps[i]]-dot;
    }

    for (i=num-1; i>=0;i--) {
        dot=0.0;
        for (j=i+1;j<num;j++)
            dot+=LU[ps[i]][j]*x_m[j];

        x_m[i]=(x_m[i]-dot)/LU[ps[i]][i];
    }
}


double evaldFdTj(int j)
{
    if(j==N-2)
        return 2 * (T[N-2] - 15);
    else if(j==N-7)
        return 2 * (T[N-7] - 15);
    else
        return 0;
}
/*double evalRi(int i)
{
    return T[i] - 0.5 * (q[0] * exp(-(i*dx) * (i*dx)) * dx * dx - 0.5 * (q[1] * exp(-((i-5)*dx) * ((i-5)*dx)) * dx * dx  + T[i-1] + T[i+1]);
}*/
double evaldRidqj(int i, int j)
{
    if (j==0)
        return  -0.5 * exp(-(i*dx) * (i*dx)) * dx * dx ;
    if (j==1)
        return  -0.5 * exp(-((i-5)*dx) * ((i-5)*dx)) * dx * dx ;
}
double evaldRidTj(int i, int j)
{
    if(i==j)
        return  1.0;
    else if (fabs(i-j) == 1)
        return -0.5;
    else
        return 0.0;
}
void solveJac()
{
    for(int itn = 0 ; itn<10000; itn++)
    {
        for (int i = 1; i < N-1; i++) {
            T[i] = 0.5 * (q[0] * exp(-(i*dx) * (i*dx)) * dx * dx + q[1] * exp(-((i-5)*dx) * ((i-5)*dx)) * dx * dx  + T[i-1] + T[i+1]);
        }
    }
}
int main()
{
    for (int i = 0; i < N; i++) {
        T[i] = 0;
    }
    for(int itn = 0; itn <= 10000; itn++)
    {
        solveJac();
        if(itn%100 == 0)
        {
            cout << "      T____________________";
            cout << endl;
            for (int i = 0; i < N; i++) {
                cout << "  " << T[i];
            }
            cout << endl;
        }

        for (int i = 0; i < N; i++) {
            dFdT[i] = evaldFdTj(i);
            for (int j = 0; j < N; j++) {
                dRdT[i][j] = evaldRidTj(i,j);
            }
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                dRdq[i][j] = evaldRidqj(i,j);
            }
        }
       /* cout << "      dFdT____________________";
        cout << endl;
        for (int i = 0; i < N; i++) {
            cout << "      " << dFdT[i];
        }
        cout << endl;*/


        /*cout << "      dRdq____________________";
        cout << endl;
        for (int i = 0; i < N; i++) {
            cout << "      " << dRdq[i];
        }
        cout << endl;


        cout << "      dRdT____________________";
        cout << endl;

        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) {
                cout << "      " << dRdT[i][j];
            }
            cout << endl;
        }
        cout << endl;*/
        for (int idx = 0; idx < M; idx++) {
            for (int j = 0; j < N-2; j++) {
                for (int i = 0; i < N-2; i++) {
                    M_[i][j] = dRdT[i+1][j+1];
                }
                b_m[j] = -dRdq[j+1][idx];
            }
            LU_decompose(N-2);
            LU_solve(N-2);
            for (int j = 0; j < N-2; j++) {
                q[idx] -= 0.1 * x_m[j] * dFdT[j+1];
            }
        }


        /* cout << "      LU____________________";
        cout << endl;

        for (int j = 0; j < N-2; j++) {
            for (int i = 0; i < N-2; i++) {
                cout << "      " << LU[i][j];
            }
            cout << endl;
        }
        cout << endl;


        cout << "      x_m____________________";
        cout << endl;*/

        /*for (int j = 0; j < N-2; j++) {
            cout << "      " << x_m[j];
        }
        cout << endl;*/

    }
    return 0;
}
















