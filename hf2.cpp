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
int matrixPLU(double** * A, int** P,const int& size)
{
    int k1,i1,j1;
    double temp;

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

int main()
{
    int N;

    // matrix numbers
    scanf("%d", &N);

    int i;
    for (i = 0; i < N; ++i)
    {
        int j, k;              //iteration
        int nMatrixM;       //Matrix size
        scanf("%d", &nMatrixM);

        //Matrix init
        double** matrix1;
        initMatrix(&matrix1, nMatrixM);

        //add to matrix elements
        for (j = 0; j < nMatrixM; ++j)
        {
            for (k = 0; k < nMatrixM; ++k)
            {
                scanf("%lf", &matrix1[j][k]);
            }
        }

        int m1;              //taskNumber
        scanf("%d", &m1);

        for (j = 0; j < m1; ++j)
        {
            //matrix for plu, and it is important to not change the base matrix

            printf( "%d .feladat\n", j);
            double** matrix;
            initMatrix(&matrix, nMatrixM);
            int miter;
            for (k = 0; k < nMatrixM; ++k)
            {
                for (miter = 0; miter < nMatrixM; ++miter)
                {
                    matrix[k][miter] = matrix1[k][miter];
                }
            }

            double c;        //shifting
            scanf("%lf", &c);

            int maxit;       //maximal of iteration
            scanf("%d", &maxit);

            double epsilon;  //convergence criteria
            scanf("%lf", &epsilon);

            //y vector
            double* y0;
            initVec(&y0, nMatrixM);

            for (k = 0; k < nMatrixM; ++k)
            {
                scanf("%lf", &y0[k]);
            }

            //sub shift from matrix
            for (k = 0; k < nMatrixM; ++k)
            {
                matrix[k][k] -= c;
            }

            //p vector to PLU
            int* p0 = new int[nMatrixM];

            for (k = 0; k < nMatrixM; ++k)
            {
                p0[k] = 0;
            }

            int singular = matrixPLU(&matrix, &p0, nMatrixM);

//            for (k = 0; k < nMatrixM; ++k)
//            {
//                for (int l = 0; l < nMatrixM; ++l)
//                {
//                    printf("%lf ", matrix[k][l]);
//                }
//                printf("\n");
//            }

            //if the matrix singular, "singular" is 1
            if(singular)
            {
                printf("%.8lf \n", c);
                delete[] p0;
                delete[] y0;
                deleteMatrix(&matrix, nMatrixM);
                continue;
            }

            //norma
            double norma = innerProduct(&y0, &y0, nMatrixM);

            if (fabs(norma) > 1e-15)
            {
                norma = sqrt(norma);
            }

            // norma is null y0 was 0
            if (fabs(norma) < 1e-15)
            {
                printf("kezdovektor\n");
                delete[] p0;
                delete[] y0;
                deleteMatrix(&matrix, nMatrixM);
                continue;
            }

            //temp array
            double* y;
            initVec(&y, nMatrixM);

            // norma
            for (k = 0; k < nMatrixM; ++k)
            {
                y[k] = y0[k] / norma;
            }

            double lambda = 0.0;

            double* tempVec;
            initVec(&tempVec, nMatrixM);

            matrixVecMultiply(&tempVec, &y, &matrix1, nMatrixM);
            printf("Ay vektor: ");
            for (k = 0; k < nMatrixM; ++k)
            {
                printf("%lf, ", tempVec[k]);
            }

            lambda = innerProduct(&tempVec, &y, nMatrixM);
            printf("lambda régi: %lf \n", lambda);

            int l;
            for (l = 0; l < maxit; ++l)
            {
                int m,n;

                //swipe y0 rows with plu
                for (m = 0; m < nMatrixM; ++m)
                {
                    y0[m] = y[p0[m]];
                }

                printf("permutált y vekor");
                for (m = 0; m < nMatrixM; ++m)
                {
                    printf("%lf ", y0[m]);
                }
                printf("\n");
                //Lz = y0
                for (m = 0; m < nMatrixM; ++m)
                {
                    double tempD = 0.0;
                    for (n = 0; n < m; ++n)
                    {
                        tempD += matrix[m][n] * y0[n];
                    }

                    y0[m] = y0[m] - tempD;

                }

                printf("z vekor");
                for (m = 0; m < nMatrixM; ++m)
                {
                    printf("%lf ", y0[m]);
                }
                printf("\n");
                //Uy0 =y
                for (m = nMatrixM - 1; m >= 0; --m)
                {
                    double tempD = 0.0;
                    for (n = m + 1; n < nMatrixM; ++n)
                    {
                        tempD += matrix[m][n] * y0[n];
                    }

                    y0[m] = y0[m] - tempD;
                    if(fabs(y0[m]) > 1e-15)
                        y0[m] /= matrix[m][m];
                }

                printf("x vekor\n");
                for (m = 0; m < nMatrixM; ++m)
                {
                    printf("%lf ", y0[m]);
                }
                printf("\n");

                norma = innerProduct(&y0, &y0,nMatrixM);
                norma = sqrt(norma);

                printf("y normálás után:\n");
                for (m = 0; m < nMatrixM; ++m)
                {
                    y[m] = y0[m] / norma;
                    printf("%lf ", y[m]);
                }
                printf("\n");

                matrixVecMultiply(&tempVec, &y, &matrix1, nMatrixM);
                double lambda1 = innerProduct(&tempVec, &y, nMatrixM);

                printf("uj lambda: %lf \n", lambda1);
                double lambdaDisc = lambda1 - lambda;
                printf(" lambda kül: %lf \n", lambdaDisc);
                double cond = epsilon * (1+fabs(lambda1));
                printf(" küszöb: %lf \n", cond);

                if (fabs(lambdaDisc) <= cond)
                {
                    lambda = lambda1;
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

                for (k = 0; k < nMatrixM; ++k)
                {
                    tempVec[k] -= lambda * y[k];
                }

                double inner = innerProduct(&tempVec, &tempVec, nMatrixM);

                if (inner <= epsilon)
                {
                    printf("siker %.8lf ", lambda);
                    int m;
                    for (m = 0; m < nMatrixM; ++m)
                    {
                        printf("%.8lf ", y[m]);
                    }

                    printf("%.8lf ", inner);
                    printf("%d\n", l);
                }
                else
                {
                    printf("sikertelen %.8lf ", lambda);
                    int m;
                    for (m = 0; m < nMatrixM; ++m)
                    {
                        printf("%.8lf ", y[m]);
                    }

                    printf("%.8lf ", inner);
                    printf("%d\n", l);
                }
            }



            delete[] tempVec;
            delete[] y;
            delete[] p0;
            delete[] y0;
            deleteMatrix(&matrix, nMatrixM);
        }

        deleteMatrix(&matrix1, nMatrixM);
    }
    return 0;
}
