#include <stdio.h>
#include <math.h>


int main(int argc, char ** argv){
    int taskNumber;
    int valtozas = 0;
    scanf("%d", &taskNumber);

    int taskProcess;
    for(taskProcess = 0; taskProcess < taskNumber; taskProcess++){
        int sing = 0;
        double temp;
        int dim;
        
        scanf("%d", &dim);

        double x[dim];
        int i,j;
        double matrix[dim][dim];

        for(i = 0; i < dim; i++)
            for(j = 0; j < dim; j++)
                scanf("%lf", &matrix[i][j]);

        double b[dim];
        double p[dim];
        
        for(i = 0; i < dim; i++){
            p[i] = i;
        }

        for(i = 0; i < dim; i++)
            scanf("%lf", &b[i]);

        int k;

        for(k = 0; k < dim-1; k++){
            int maxRow = k;
            for(i = k+1; i < dim; i++){
                if(fabs(matrix[i][k]) > fabs(matrix[maxRow][k]))
                    maxRow = i;
            }

            

                if(maxRow != k){
                        for(i = 0; i < dim; i++){
                            temp = matrix[k][i];
                            matrix[k][i] = matrix[maxRow][i];
                            matrix[maxRow][i] = temp;
                        }
                    valtozas++;
                    
                    p[k] = maxRow;
                    p[maxRow] = k;

                }

                if(fabs(matrix[k][k]) < 1e-15)
                    sing = 1;
            
                p[k] = maxRow;
                p[maxRow] = k;
                    
                for(i = k+1; i < dim - 1; i++){
                
                    matrix[i][k] = matrix[i][k]/matrix[k][k];

                    for(j = k+1; j < dim - 1; j++){
                        matrix[i][j] = matrix[i][j] - (matrix[i][k] * matrix[k][j]);     
                    }
                
          }
        }
          if(fabs(matrix[dim-1][dim-1]) < 1e-15)
                sing = 1;
          else{
                double n_b[dim];

                for(i = 0; i < dim; i++){
                    int intTemp = p[i];
                     n_b[i] = b[intTemp];
                }
                    
                temp = 0;       
                    
                //b vektor már nem kell így felhasználom
                for(i = 0; i < dim; i++){
                    for(j = 0; j < i; j++)
                        temp += (matrix[i][j] * b[j]);    
                       
                    b[i] = n_b[i] - temp;    
                }
                
                temp = 0;

                for(i = dim - 1; i >= 0; i--){
                    for(j = i+1; j < dim; j++)
                        temp += (matrix[i][j] * x[j]);
                
                    x[i] = (b[i] - temp)/ matrix[i][j];
                } 

                    
         }

         if(sing == 1)
             printf("singularis");
         else
         {
            double det = 1;
            
            for(i = 0; i < dim; i++)
                det *= matrix[i][i];

            if(valtozas % 2 == 1)
                det *= -1;
            
            printf("%lf ",det);

            for(i = 0; i < dim; i++)
                printf("%lf ", x[i]);
             
         }
          
        
     printf("\n");   
    }
}
