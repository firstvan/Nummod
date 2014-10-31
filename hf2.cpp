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

void deleteVec(double ** A,const int& size){
    delete[] *A;
}

int matrixPLU(double *** A, double ** P,const int& size){
    int k,i,j;
    double temp;

    for (i = 0; i < size; ++i)
    {
        (*P)[i] = i;
    }

    for (k = 0; i < size - 1; ++i)
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

    return 0;
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

            double * y;
            initVec(&y, nMatrixM);

            for (k = 0; k < nMatrixM; ++k)
            {
                scanf("%lf", &y[j]);
            }

            double * p;
            initVec(&p, nMatrixM);

            int singular = matrixPLU(&matrix, &p, nMatrixM);

            if(!singular){
                printf("singular matrix, c=%.8lf", c);
                deleteVec(&p, nMatrixM);
                deleteVec(&y, nMatrixM);
                break;
            }

            for (k = 0; k < nMatrixM; ++k)
            {
                matrix[k][k] -= c;
            }

            deleteVec(&p, nMatrixM);
            deleteVec(&y, nMatrixM);
        }

        deleteMatrix(&matrix, nMatrixM);
    }

    return 0;
}
