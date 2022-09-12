#include "../../SCFT.h"

void removeColumn(Eigen:: MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

void SCFT::mixDensityAndersonRelax() {
     /* Column to update on this iteration */
     int col = anderIter % anderM; // column to update
     double temp;

     /* Loop through each component to fill in matrices*/
     for(int i = 0; i < numComponents+1; i++) {
         /* Fill matrices */
         for(int k = 0; k < M; k++) {
              X_And(i*M + k,col)  = X[i*M+k];
              G_And(i*M + k,col)  = F[i*M+k];
         }
     }

     if(anderIter > 500 && anderIter < 10*minIters ) {
         int colRemove = (anderIter+1) % anderM;

         /* Initialize matrices specific to this run */
         Eigen::VectorXd lambdas;
         Eigen::VectorXd g_k     = G_And.col(col);
         Eigen::VectorXd x_k     = X_And.col(col);
         Eigen::MatrixXd dG_temp  = G_And;
         Eigen::MatrixXd dX_temp  = X_And;

         /* Make difference matrices */
         for(int i = anderIter + anderM; i > anderIter; i--) {
            dG_temp.col(i % anderM) = G_And.col(i % anderM) - G_And.col((i-1) % anderM);
            dX_temp.col(i % anderM) = X_And.col(i % anderM) - X_And.col((i-1) % anderM);
         }

         /* Remove the "zero" column */
         removeColumn(dG_temp,colRemove);
         removeColumn(dX_temp,colRemove);

         /* Solve the least square problem to obtain the coefficients */
         lambdas = dG_temp.fullPivHouseholderQr().solve(g_k);

         /* Update the density */
         for(int i = 0; i < numComponents+1; i++) {
             for(int k = 0; k < M; k++) {
                   if(i < numComponents) {
                       temp = (1-mixF)*(x_k(i*M+k) - lambdas.dot(dX_temp.row(i*M+k)) ) + mixF*(x_k(i*M+k) + g_k(i*M+k) - lambdas.dot( dX_temp.row(i*M+k) + dG_temp.row(i*M+k)));
                       if(temp > 1) temp = 1;
                       if(temp < 0) temp = X[i*M+k]*0.9;
                       X[i*M+k] = temp;
                   }
                   else {
                       temp = (1-0.25*mixF)*(x_k(i*M+k) - lambdas.dot(dX_temp.row(i*M+k)) ) + 0.25*mixF*(x_k(i*M+k) + g_k(i*M+k) - lambdas.dot( dX_temp.row(i*M+k) + dG_temp.row(i*M+k)));
                       if(fabs(temp) < 1E-20) temp = 0;
                       X[i*M+k]  = temp;
                   }
             }
         }
      }
      else {
          mixDensityPicardRelax();
      }

     /* Increase iteration parameter */
     anderIter++;
     if(anderIter > 20*minIters) anderIter = 0;
}

void SCFT::mixDensityAndersonRelaxCont() {
     /* Column to update on this iteration */
     int col = anderIter % anderM; // column to update
     double temp;

     /* Loop through each component to fill in matrices*/
     for(int i = 0; i < M*(numComponents+1); i++) {
          X_And_Cont(i,col)  = X[i];
          G_And_Cont(i,col)  = F[i];
     }
     X_And_Cont(M*(numComponents+1),col)  = X_Cont;
     G_And_Cont(M*(numComponents+1),col)  = F_Cont;


     if(anderIter > 500 && anderIter < 19*minIters) {
         int colRemove = (anderIter+1) % anderM;

         /* Initialize matrices specific to this run */
         Eigen::VectorXd lambdas;
         Eigen::VectorXd g_k      = G_And_Cont.col(col);
         Eigen::VectorXd x_k      = X_And_Cont.col(col);
         Eigen::MatrixXd dG_temp  = G_And_Cont;
         Eigen::MatrixXd dX_temp  = X_And_Cont;

         /* Make difference matrices */
         for(int i = anderIter + anderM; i > anderIter; i--) {
            dG_temp.col(i % anderM) = G_And_Cont.col(i % anderM) - G_And_Cont.col((i-1) % anderM);
            dX_temp.col(i % anderM) = X_And_Cont.col(i % anderM) - X_And_Cont.col((i-1) % anderM);
         }

         /* Remove the "zero" column */
         removeColumn(dG_temp,colRemove);
         removeColumn(dX_temp,colRemove);

         /* Solve the least square problem to obtain the coefficients */
         lambdas = dG_temp.fullPivHouseholderQr().solve(g_k);

         /* Update the density */
         for(int i = 0; i < numComponents+1; i++) {
             for(int k = 0; k < M; k++) {
                   if(i < numComponents) {
                       temp = (1-mixF)*(x_k(i*M+k) - lambdas.dot(dX_temp.row(i*M+k)) ) + mixF*(x_k(i*M+k) + g_k(i*M+k) - lambdas.dot( dX_temp.row(i*M+k) + dG_temp.row(i*M+k)));
                       if(temp > 1) temp = 1;
                       if(temp < 0) temp = X[i*M+k]*0.9;
                       X[i*M+k] = temp;
                   }
                   else {
                       temp = (1-0.25*mixF)*(x_k(i*M+k) - lambdas.dot(dX_temp.row(i*M+k)) ) + 0.25*mixF*(x_k(i*M+k) + g_k(i*M+k) - lambdas.dot( dX_temp.row(i*M+k) + dG_temp.row(i*M+k)));
                       if(fabs(temp) < 1E-20) temp = 0;
                       X[i*M+k]  = temp;
                   }
             }
         }

         int i = M*(numComponents+1);
         X_Cont  = (1-mixF)*(x_k(i) - lambdas.dot(dX_temp.row(i)) ) + mixF*(x_k(i) + g_k(i) - lambdas.dot( dX_temp.row(i) + dG_temp.row(i)));
      }
      else {
           mixDensityPicardRelax();
           X_Cont = (1-mixF)*X_Cont+ mixF*(X_Cont + F_Cont);
      }

     /* Increase iteration parameter */
     anderIter++;
     if(anderIter > 20*minIters) anderIter = 0;
}
