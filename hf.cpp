#include <stdio.h>
#include <math.h>


int main(int argc, char ** argv){
    int taskNumber;
        scanf("%d", &taskNumber);

    int taskProcess;
    for(taskProcess = 0; taskProcess < taskNumber; taskProcess++){
        int csere = 0;
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
        for(i = 0; i < dim; i++)
            scanf("%lf", &b[i]);
        
        double p[dim];
        
        for(i = 0; i < dim; i++){
            p[i] = i;
        }

        int k;

        for(k = 0; k < dim - 1; k++){
            
            int maxRow = k;
            for(i = k+1; i < dim; i++){
                if(fabs(matrix[i][k]) > fabs(matrix[maxRow][k])){
                   maxRow = i;
                }
            }     //jó       

            if(maxRow != k){
                for(i = 0; i < dim; i++){
                    temp = matrix[k][i];
                    matrix[k][i] = matrix[maxRow][i];
                    matrix[maxRow][i] = temp;
                }
                
                csere++;
                    
                p[k] = maxRow;
                p[maxRow] = k;

            }

            if(fabs(matrix[k][k]) < 1e-15)
                sing = 1;
            
               
            for(i = k+1; i < dim; i++){
                matrix[i][k] = matrix[i][k]/matrix[k][k];
                
                for(j = k+1; j < dim; j++){
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
            	temp = 0;  
                for(j = 0; j < i; j++)
                    temp += (matrix[i][j] * b[j]);    
                      
                b[i] = n_b[i] - temp;
            }

            for(i = dim - 1; i >= 0; i--){
            	temp = 0;
                for(j = i+1; j < dim; j++)
 	                temp += (matrix[i][j] * x[j]);
                

            	x[i] = (b[i] - temp);
            	if(fabs(x[i]) > 1e-15)
            		x[i] = x[i] / matrix[i][i];
            }                     
         }

         if(sing == 1)
             printf("singularis");
         else
         {
            double det = 1;
            
            for(i = 0; i < dim; i++)
                det *= matrix[i][i];

            if(csere % 2 == 1)
            	det *= -1;
            
            printf("%lf ",det);

            for(i = 0; i < dim; i++)
                printf("%lf ", x[i]);
             
         }

     printf("\n");   
    }
}
