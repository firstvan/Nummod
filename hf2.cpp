#include <stdio.h>
#include <math.h>


void initMatrix(double *** A, const int & size){

    int j;

    //Matrix init
    *A = new double*[size];
    for (j = 0; j < size; ++j)
    {
        (*A)[j] = new double[size];
    }

}

void deleteMatrix(double *** A, const int &size){

    int j;

    //Delete matrix
    for (j = 0; j < size; ++j)
    {
        delete[] (*A)[j];
    }

    delete[] (*A);
}

void initVec(double ** A,const int& size){
    *A = new double[size];
}

void deleteVec(double ** A){
    delete[] *A;
}

int matrixPLU(double *** A, int ** P,const int& size){
    int k,i,j;
    double temp;

    for (i = 0; i < size; ++i)
    {
        (*P)[i] = i;
    }

    for (k = 0; k < size - 1; ++k)
    {
        int maxRow = k;

        for (i = k+1; i < size; ++i)
            if(fabs((*A)[i][k]) > fabs((*A)[maxRow][k]))
                maxRow = i;

        if(maxRow != k){
            for (i = 0; i < size; ++i)
            {
                temp = (*A)[k][i];
                (*A)[k][i] = (*A)[maxRow][i];
                (*A)[maxRow][i] = temp;
            }

            int tempI = (*P)[k];
            (*P)[k] = (*P)[maxRow];
            (*P)[maxRow] = tempI;
        }

        if(fabs((*A)[k][k]) < 1e-15)
            return 1;

        for (i = k + 1; i < size; ++i)
        {
            (*A)[i][k] /= (*A)[k][k];

            for (j = k + 1; j < size; ++j)
            {
                (*A)[i][j] -= (*A)[i][k] * (*A)[k][j];
            }
        }

    }

    if(fabs((*A)[size-1][size-1]) < 0e-15)
        return 1;

    return 0;
}

void matrixVecMultiply(double ** T, double ** Y, double *** M, const int& size){
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

double innerProduct(double ** X, double ** Y,const int& size){
    double re = 0;
    int i;

    for (i = 0; i < size; ++i)
    {
        re += (*X)[i] * (*Y)[i];
    }

    return re;
}

int main(){
    int N;

    scanf("%d", &N);

    int i;
    for (i = 0; i < N; ++i)
    {
        int j, k;              //iteration
        int nMatrixM;       //Matrix size
        scanf("%d", &nMatrixM);

        //Matrix init
        double ** matrix;
        initMatrix(&matrix, nMatrixM);


        //add to matrix elements
        for (j = 0; j < nMatrixM; ++j)
        {
            for (k = 0; k < nMatrixM; ++k)
            {
                scanf("%lf", &matrix[j][k]);
            }
        }

        int m;              //taskNumber
        scanf("%d", &m);

        for (j = 0; j < m; ++j)
        {
            double c;        //shifting
            scanf("%lf", &c);

            int maxit;       //maximal of iteration
            scanf("%d", &maxit);

            double epsilon;
            scanf("%lf", &epsilon);

            double * y0;
            initVec(&y0, nMatrixM);

            for (k = 0; k < nMatrixM; ++k)
            {
                scanf("%lf", &y0[k]);
            }

            int * p;
            p = new int[nMatrixM];

            int singular = matrixPLU(&matrix, &p, nMatrixM);

            if(singular){
                printf("%.8lf\n", c);
                delete[] p;
                delete[] y0;
                continue;
            }

            for (k = 0; k < nMatrixM; ++k)
            {
                matrix[k][k] -= c;
            }

            double norma = 0.0;

            for (k = 0; k < nMatrixM; ++k)
            {
                norma += y0[k] * y0[k];
            }

            if (fabs(norma) > 1e-15)
            {
                norma = sqrt(norma);
            }


            if (fabs(norma) < 1e-15)
            {
                printf("kezdovektor\n");
                delete[] p;
                delete[] y0;
                continue;
            }

            double * y;
            initVec(&y, nMatrixM);

            for (k = 0; k < nMatrixM; ++k)
            {
                y[k] = y0[k] / norma;
            }

            double lambda, lambda1;
            double * tempVec;
            initVec(&tempVec, nMatrixM);
            matrixVecMultiply(&tempVec, &y, &matrix, nMatrixM);

            lambda = innerProduct(&tempVec, &y, nMatrixM);


            int l;
            for (l = 0; l < maxit; ++l)
            {
                int m,n;
                for (m = 0; m < nMatrixM; ++m)
                {
                    double tempD = y[m];
                    y[m] = y[p[m]];
                    y[p[m]] = tempD;
                }

                for (m = 0; m < nMatrixM; ++m)
                {
                    double tempD = 0;
                    for (n = 0; n < m; ++n)
                    {
                        tempD += matrix[i][j] * y[j];
                    }

                    y[j] -= tempD;
                }

                for (m = nMatrixM - 1; m >= 0; --m)
                {
                    double tempD = 0;
                    for (n = m + 1; n < nMatrixM; ++n)
                    {
                        tempD += matrix[m][n] * y0[n];
                    }

                    y0[m] = y[m] - tempD;
                    if(fabs(y0[m]) > 1e-15)
                        y0[m] /= matrix[m][m];
                }

                norma = innerProduct(&y0, &y0,nMatrixM);

                norma = sqrt(norma);

                for (m = 0; m < nMatrixM; ++m)
                {
                    y[m] = y0[m] / norma;
                }

                matrixVecMultiply(&tempVec, &y, &matrix, nMatrixM);

                lambda1 = innerProduct(&tempVec, &y, nMatrixM);

                if (fabs(lambda1 - lambda) <= epsilon * (1 + lambda1))
                {
                    break;
                }

                lambda = lambda1;

            }

            if (l == maxit)
            {
                printf("maxit\n");
            }
            else
            {
                matrixVecMultiply(&tempVec, &y, &matrix, nMatrixM);

                for (k = 0; k < nMatrixM; ++k)
                {
                    tempVec[k] -= lambda1 * tempVec[k];
                }

                if (innerProduct(&tempVec, &tempVec, nMatrixM) < epsilon)
                {
                    printf("siker %lf ", lambda1);
                    for (m = 0; m < nMatrixM; ++m)
                    {
                       printf("%lf ", tempVec[m]); 
                    }

                    printf("%.8lf ", innerProduct(&tempVec, &tempVec, nMatrixM));
                    printf("%d\n", l);
                }
                else
                {
                    printf("Sikertelen\n");
                }
            }



            delete[] tempVec;
            delete[] y;
            delete[] p;
            delete[] y0;
        }


        for (i = 0; i < nMatrixM; ++i)
        {
            delete[] matrix[i];
        }

        delete[] matrix;
    }

    return 0;
}
