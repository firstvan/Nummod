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

double normaV(double * A, int size)
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

int main()
{
    int N;

    // numbers of process
    scanf("%d", &N);

    int i;
    for (i = 0; i < N; ++i)
    {
        int vecSize; //vector size
        scanf("%d", &vecSize);

        //vector init
        double* xKezdo;
        initMatrix(&xKezdo, vecSize);

        //read element of vector
            for (int j = 0; j < vecSize; ++j)
            {
                scanf("%lf", &xKezdo[j]);
            }
            
            for(int j = 0; j < vecSize; j++)
	    {
	      printf("%lf ,", xKezdo[j]);
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
