#include <stdio.h>
#include <math.h>

//init matrix in dynamic memory
void initMatrix(double** * A, const int& size)
{

    int j;

    //Matrix init
    (*A) = new double*[size];
    for (j = 0; j < size; ++j)
    {
        (*A)[j] = new double[size];
    }

}

//delete matrix from dynamic memory
void deleteMatrix(double** * A, const int& size)
{

    int j;

    //Delete matrix
    for (j = 0; j < size; ++j)
    {
        delete[] (*A)[j];
    }

    delete[] (*A);
}

//init array in dynamic memory
void initVec(double** A,const int& size)
{
    (*A) = new double[size];
}

//delete array from dynamic memory
void deleteVec(double** A)
{
    delete[] (*A);
}

// do the PLU
int matrixPLU( double** * A, int** P,const int& size)
{
    int k1,i1,j1;
    double temp;

    //veszem a -1 szeresét a P-nek(ami az f);
    /*for(int k0 = 0; k0 < size; k0++)
    {
      (*P)[k0] = -1 * (*P)[k0];
    }
    */

    //P init
    for (i1 = 0; i1 < size; ++i1)
    {
        (*P)[i1] = i1;
    }

    for (k1 = 0; k1 < size - 1; ++k1)
    {
        int maxRow = k1;

        //which is the biggest value in column
        for (i1 = k1+1; i1 < size; ++i1)
            if(fabs((*A)[i1][k1]) > fabs((*A)[maxRow][k1]))
                maxRow = i1;

        //swipe rows
        if(maxRow != k1)
        {
            for (i1 = 0; i1 < size; ++i1)
            {
                temp = (*A)[k1][i1];
                (*A)[k1][i1] = (*A)[maxRow][i1];
                (*A)[maxRow][i1] = temp;
            }

            //save swaps
            int intTemp = (*P)[k1];
            (*P)[k1] = (*P)[maxRow];
            (*P)[maxRow] = intTemp;
        }

        //singular?
        if(fabs((*A)[k1][k1]) < 1e-15)
            return 1;

        //do the LU
        for (i1 = k1 + 1; i1 < size; ++i1)
        {
            (*A)[i1][k1] /= (*A)[k1][k1];

            for (j1 = k1 + 1; j1 < size; ++j1)
            {
                (*A)[i1][j1] -= (*A)[i1][k1] * (*A)[k1][j1];
            }
        }

    }

    //singular?
    if(fabs((*A)[size-1][size-1]) < 1e-15)
        return 1;

    return 0;
}

//simple matrix vector multiply
void matrixVecMultiply(double** T, double** Y, double** * M, const int& size)
{
    int i,j;

    for (i = 0; i < size; ++i)
    {
        (*T)[i] = 0;
        for (j = 0; j < size; ++j)
        {
            (*T)[i] += (*M)[i][j] * (*Y)[j];
        }
    }
}

//inner product
double innerProduct(double** X, double** Y,const int& size)
{
    double re = 0.0;
    int i;

    for (i = 0; i < size; ++i)
    {
        re += (*X)[i] * (*Y)[i];
    }

    return re;
}

double normaV(double* A, int size)
{
    double max = fabs(A[0]);
    for(int i = 1; i < size; i++)
    {

        if(fabs(A[i]) > max)
        {
            max = fabs(A[i]);
        }

    }

    return max;
}

void f(int s, double* temp, double* A, int size)
{
    switch(s)
    {
    case 1:
        temp[0]= (-1 * pow(A[0], 2)) + A[2] + 3;
        temp[1]= (-1 * A[0]) + 2 * pow(A[1], 2) - pow(A[2], 2) - 3;
        temp[2]= A[1] - 3 * pow(A[2], 2) + 2;
        break;
    case 2:
        temp[0]= 2 * pow(A[0], 2) - A[1] - 1;
        temp[1]= -1 * A[0] + 2 * pow(A[1],2) - 1;
        break;
    case 3:
        temp[0]= -4 * A[0] + cos(2 * A[0] - A[1]) - 3;
        temp[1]= sin(A[0]) - 3  * A[1] - 2;
        break;
    case 4:
        temp[0]= A[0] *pow(A[1],2) -4* A[0] * A[1] + 4 * A[0] - 1;
        temp[1]= exp(A[0]-1) - A[1] + 1;
        break;
    }
}

void doJacobiMatrix(int s, double** * jMatrix, double* fgv)
{
    switch(s)
    {
    case 1:
        (*jMatrix)[0][0] = -2 * fgv[0];
        (*jMatrix)[0][1] = 0;
        (*jMatrix)[0][2] = 1;
        (*jMatrix)[1][0] = -1;
        (*jMatrix)[1][1] = 4 * fgv[1];
        (*jMatrix)[1][2] = -2 * fgv[2];
        (*jMatrix)[2][0] = 0;
        (*jMatrix)[2][1] = 1;
        (*jMatrix)[2][2] = -6 * fgv[2];
        break;
    case 2:
        (*jMatrix)[0][0] = 4 * fgv[0];
        (*jMatrix)[0][1] = -1;
        (*jMatrix)[1][0] = -1;
        (*jMatrix)[1][1] = 4 * fgv[1];
        break;
    case 3:
        (*jMatrix)[0][0] = /*  2*sin(fgv[1]-2*fgv[0])-4;*/ -4 - 2*sin(2 * fgv[0] - fgv[1]);
        (*jMatrix)[0][1] = /*-1 * sin(fgv[1]-2 * fgv[1]); */sin(2*fgv[0]-fgv[1]);
        (*jMatrix)[1][0] = cos(fgv[0]);
        (*jMatrix)[1][1] = -3;
        break;
    case 4:
        (*jMatrix)[0][0] = pow(fgv[1],2)-4*fgv[1]+4;
        (*jMatrix)[0][1] = 2 * fgv[0] * fgv[1] - 4 * fgv[0];
        (*jMatrix)[1][0] = exp(fgv[0]-1);
        (*jMatrix)[1][1] = -1;
        break;
    }
}

int main()
{
    int N;

    // numbers of process
    scanf("%d", &N);

    int i;
    for (i = 0; i < N; ++i)
    {
        int maxit;
        double epszilon;
        int sorszam;


        scanf("%d", &sorszam);
        scanf("%d", &maxit);
        scanf("%lf", &epszilon);

        int vecSize = 2; //vector size

        if(sorszam == 1)
            vecSize = 3;

        //vector init
        double* xKezdo;
        initVec(&xKezdo, vecSize);

        //read element of vector
        for (int j = 0; j < vecSize; ++j)
        {
            scanf("%lf", &xKezdo[j]);
        }


//        for(int j = 0; j < vecSize; j++)
//        {
//            printf("%lf ,", xKezdo[j]);
//        }



        double t = 1.0;

        double* fx;
        initVec(&fx, vecSize);

        f(sorszam, fx, xKezdo, vecSize);

        double normaF0 = normaV(fx, vecSize);
        double normaFK = normaF0;

        bool kilep = false;

        int j;
        for(j = 1; j < maxit; j++)
        {
            if(j != 1) // ez is utólag lett belerakva mármint a feltétel
            {
                f(sorszam, fx, xKezdo, vecSize);
            }

            double** jacobiMatrix;
            initMatrix(&jacobiMatrix, vecSize);

            doJacobiMatrix(sorszam, &jacobiMatrix, xKezdo);

//            for(int m = 0; m < vecSize; m++)
//            {
//                for(int n = 0; n < vecSize; n++)
//                {
//                    printf("%lf ,", jacobiMatrix[m][n]);
//                }
//                printf("\n");
//            }
//            printf("---");

            int* P;
            P = new int[vecSize];

            int sing = matrixPLU(&jacobiMatrix, &P, vecSize);
            if(sing)
            {
                printf("szingularis ");
                for (int m = 0; m < vecSize; ++m)
                {
                    printf("%.8lf ", xKezdo[m]);
                }
                printf("\n");
                deleteMatrix(&jacobiMatrix, vecSize);
                delete[] P;
                break;
            }

            double* temp;
            initVec(&temp, vecSize);

            for (int m = 0; m < vecSize; ++m)
            {
                temp[m] = -1 * fx[P[m]]; //ha nem jó figyelni  a -1 *  re
            }

            for (int m = 0; m < vecSize; ++m)
            {
                double tempD = 0.0;
                for (int n = 0; n < m; ++n)
                {
                    tempD += jacobiMatrix[m][n] * temp[n];
                }

                temp[m] = temp[m] - tempD;

            }

            for (int m = vecSize - 1; m >= 0; --m)
            {
                double tempD = 0.0;
                for (int n = m + 1; n < vecSize; ++n)
                {
                    tempD += jacobiMatrix[m][n] * temp[n];
                }

                temp[m] = temp[m] - tempD;
                temp[m] /= jacobiMatrix[m][m];
            }

//            normaFK = normaV(xKezdo, vecSize);
            bool siker = false;
            int l;
            for (l = 0; l < 8; ++l)
            {


                double* y;
                initVec(&y, vecSize);

                for (int m = 0; m < vecSize; ++m)
                {
                    y[m] = xKezdo[m] + (t * temp[m]);
                }

                double* fxT;
                initVec(&fxT, vecSize);

                f(sorszam, fxT, y, vecSize);

                double cond1 = normaV(fxT, vecSize);
                //normaFK = normaV(xKezdo, vecSize);
                if (cond1 <= normaFK)
                {
                    for (int m = 0; m < vecSize; ++m)
                    {
                        xKezdo[m] = y[m];
                    }

                    normaFK = cond1;
                    siker =true;
                    delete[] fxT;
                    delete[] y;
                    break;
                }


                t /= static_cast<double>(2);

                if(t <= pow(10, -3))
                {
                    printf("sikertelen ");
                    for (int m = 0; m < vecSize; ++m)
                    {
                        printf("%.8lf ", y[m]);
                    }
                    printf("\n");
                    delete[] fxT;
                    delete[] y;
                    kilep = true;
                    break;
                }



                delete[] fxT;
                delete[] y;
            }
            if(l == 8 && siker==true)
            {
                printf("sikertelen ");
                for (int m = 0; m < vecSize; ++m)
                {
                    printf("%.8lf ", xKezdo[m]);
                }
                printf("\n");
                break;

            }

            if(l == 0)
            {
                double tempT = t * 1.5;
                if(tempT < 1)
                {
                    t = tempT;
                }
                else
                {
                    t = 1.0;
                }
            }

            if(kilep)
            {
                delete[] P;
                delete[] temp;
                deleteMatrix(&jacobiMatrix, vecSize);
                break;
            }


            if(normaFK <= epszilon * (1+normaF0)) //ha sok wa akkor kiszedni az =
            {
                printf("siker ");
                for (int m = 0; m < vecSize; ++m)
                {
                    printf("%.8lf ", xKezdo[m]);
                }
                printf("%.8lf %d", normaFK, j);
                printf("\n");
                delete[] P;
                delete[] temp;
                deleteMatrix(&jacobiMatrix, vecSize);
                break;
            }



            delete[] P;
            delete[] temp;
            deleteMatrix(&jacobiMatrix, vecSize);
        }

        if(j == maxit)
        {
            printf("maxit");
            printf("\n");
        }


        delete[] fx;
        delete[] xKezdo;

    }


    return 0;
}

