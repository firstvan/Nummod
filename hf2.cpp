#include <stdio.h>
#include <math.h>


void initMatrix(double *** A, const int & size){

    int j;

    //Matrix init
    (*A) = new double*[size];
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
    (*A) = new double[size];
}

void deleteVec(double ** A){
    delete[] (*A);
}

int matrixPLU(double *** A, int ** P,const int& size){
    int k1,i1,j1;
    double temp;

    for (i1 = 0; i1 < size; ++i1)
    {
       (*P)[i1] = i1;
    }

    for (k1 = 0; k1 < size - 1; ++k1)
    {
        int maxRow = k1;

        for (i1 = k1+1; i1 < size; ++i1)
            if(fabs((*A)[i1][k1]) > fabs((*A)[maxRow][k1]))
                maxRow = i1;

        if(maxRow != k1){
            for (i1 = 0; i1 < size; ++i1)
            {
                temp = (*A)[k1][i1];
                (*A)[k1][i1] = (*A)[maxRow][i1];
                (*A)[maxRow][i1] = temp;
            }

            int tempI = (*P)[k1];
            (*P)[k1] = (*P)[maxRow];
            (*P)[maxRow] = tempI;
        }

        if(fabs((*A)[k1][k1]) < 1e-15)
            return 1;

        for (i1 = k1 + 1; i1 < size; ++i1)
        {
            (*A)[i1][k1] /= (*A)[k1][k1];

            for (j1 = k1 + 1; j1 < size; ++j1)
            {
                (*A)[i1][j1] -= (*A)[i1][k1]* (*A)[k1][j1];
            }
        }

    }

    if(fabs((*A)[size-1][size-1]) < 1e-15)
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
    double re = 0.0;
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

        int m1;              //taskNumber
        scanf("%d", &m1);

        for (j = 0; j < m1; ++j)
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


            for (k = 0; k < nMatrixM; ++k)
            {
                matrix[k][k] -= c;
            }

            int * p0 = new int[nMatrixM];


            for (k = 0; k < nMatrixM; ++k)
            {
                p0[k] = 0;
            }

            int singular = matrixPLU(&matrix, &p0, nMatrixM);

            if(singular){
                printf("%.8lf\n", c);
                delete[] p0;
                delete[] y0;
                continue;
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
                delete[] p0;
                delete[] y0;
                continue;
            }

            double * y;
            initVec(&y, nMatrixM);

            for (k = 0; k < nMatrixM; ++k)
            {
                y[k] = y0[k] / norma;
            }

            double lambda = 0.0, lambda1 = 0.0;
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
                    y0[m] = y[p0[m]];
                }

                for (m = 0; m < nMatrixM; ++m)
                {
                    for (n = 0; n < m; ++n)
                    {
                        y0[m] -= matrix[m][n] * y0[n];
                    }

                }

                for (m = nMatrixM - 1; m >= 0; --m)
                {
                    double tempD = 0;
                    for (n = m + 1; n < nMatrixM; ++n)
                    {
                        tempD += matrix[m][n] * y0[n];
                    }

                    y0[m] = y0[m] - tempD;
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
                
                double lambdaDisc = lambda1-lambda;
                double cond = epsilon * (1+fabs(lambda1));
//                printf("%.8lf, %.8lf\n", fabs(lambdaDisc), cond);
                if (fabs(lambdaDisc) <= cond) 
                {
                    break;
                }

                lambda = lambda1;

            }
            printf("%d, %d \n", l, maxit);
            if (l == maxit)
            {
                printf("maxit\n");
            }
            else
            {
                matrixVecMultiply(&tempVec, &y, &matrix, nMatrixM);

                for (k = 0; k < nMatrixM; ++k)
                {
                    tempVec[k] -= lambda1;
                }

                if (innerProduct(&tempVec, &tempVec, nMatrixM) < epsilon)
                {
                    printf("siker %lf ", lambda1);
                    int m;
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
            delete[] p0;
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
