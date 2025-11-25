#include "mainfx.h"
#include <cmath>
#include <vector>
#include <iostream>


void PredUfill(double M, double N, std::vector<std::vector<double>> &MatU, std::vector<std::vector<double>> &MatV, std::vector<std::vector<double>> &PredU
    , std::vector<std::vector<double>> &RU, double rho, double dx, double dy, double dt, double mu){

    double un, ue, us, uw;
    double uN, uE, uS, uW, uP;
    double vs,vn;
    double q = dy/dx, dV = dx*dy;

    for(int i = 1; i < M+1; i++){
        for(int j = 1; j < N+2; j++){
            un = avg(MatU[i][j],MatU[i-1][j])*(i!=1);
            ue = avg(MatU[i][j],MatU[i][j+1]);
            us = avg(MatU[i][j],MatU[i+1][j])*(i!=M);
            uw = avg(MatU[i][j],MatU[i][j-1]);
            vn = avg(MatV[i][j-1],MatV[i][j])*(j!=1&&j!=N+1)*(i!=1);
            vs = avg(MatV[i+1][j-1],MatV[i+1][j])*(j!=1&&j!=N+1)*(i!=M);
    
            uN = MatU[i-1][j];
            uE = MatU[i][j+1];
            uS = MatU[i+1][j];
            uW = MatU[i][j-1];
            uP = MatU[i][j];
            
            
            PredU[i][j] = MatU[i][j] - RU[i][j]; //aixi fem servir el RU calculat al timestep previ

            RU[i][j] = dt*(diffx(uN,uE,uS,uW,uP,mu,q)-convx(rho,dx,un,ue,us,uw,vn,vs,q))/(2*rho*dV);
            PredU[i][j] += 3*RU[i][j];
        }
    }
}


void PredVfill(double M, double N, std::vector<std::vector<double>> &MatU, std::vector<std::vector<double>> &MatV, std::vector<std::vector<double>> &PredV
    , std::vector<std::vector<double>> &RV, std::vector<std::vector<double>> &Temp, std::vector<std::vector<double>> &nuTemp 
    , double dt, double beta, double g, double rho, double dx, double dy, double mu){
       
    double vn, ve, vs, vw;
    double vN, vE, vS, vW, vP;
    double ue,uw;
    double dV = dx*dy, q = dy/dx;

    for(int i = 1; i < M+2; i++){
        for(int j = 1; j < N+1; j++){
            vn = avg(MatV[i][j],MatV[i-1][j]);
            ve = avg(MatV[i][j],MatV[i][j+1])*(j!=N);
            vs = avg(MatV[i][j],MatV[i+1][j]);
            vw = avg(MatV[i][j],MatV[i][j-1])*(j!=1);
            ue = avg(MatU[i-1][j+1],MatU[i][j+1])*(j!=N)*(i!=1&&i!=M+1);
            uw = avg(MatU[i-1][j],MatU[i][j])*(j!=1)*(i!=1&&i!=M+1);
    
            vN = MatV[i-1][j];
            vE = MatV[i][j+1];
            vS = MatV[i+1][j];
            vW = MatV[i][j-1];
            vP = MatV[i][j];
    
            //A predictor V tenim el terme de boussinesq en positiu, pero multiplicat pel vector g, que es negatiu
            PredV[i][j] = MatV[i][j] - RV[i][j] //de nou aixi fem servir el RV del timestep previ
                                     - 0.5*beta*(avg(Temp[i-1][j],Temp[i][j])-0.5)*g*dt;// el - es pq el vector g es -g                    

            RV[i][j] = dt*(diffy(vN,vE,vS,vW,vP,mu,q)-convy(rho,dx,vn,ve,vs,vw,ue,uw,q))/(2*rho*dV);
            PredV[i][j] += 3*RV[i][j] + 1.5*beta*(avg(nuTemp[i-1][j],nuTemp[i][j])-0.5)*g*dt;
            //no se si tenia el signe del boussinesq malament...
        }
    }
}


void ConjugateGradient(std::vector<double> &PressureVec, std::vector<double> &R, std::vector<double> &ct, double M, double N, double epsilon, double maxiter
    , std::vector<double> &P, std::vector<double> &CoefN, std::vector<double> &CoefE, std::vector<double> &CoefS
    , std::vector<double> &CoefW, std::vector<double> &CoefP,int k){
    //to be clear the matrix must be positive definite and symmetric
    //the specific one I'm working on for this problem is symmetric
    //el que si que no es positive-definite, justament pel que estava fent
    //per comprovar que el multiplicador de matrius anava, i es que
    //A·(1....1) = 0 aixi que (1...1)^T A (1...1) = 0 
    //En fi tornem al gauss seidel suposo tot i que justament tampoc tinc que la diagonal
    //domini

    matVecMult2DGrid(PressureVec,CoefN,CoefE,CoefS,CoefW,CoefP,M,N);
    R = CombLin(PressureVec,matVecMult2DGrid(PressureVec,CoefN,CoefE,CoefS,CoefW,CoefP,M,N),1,-1);
    int iter = 0;
    P = R;
    double b = 0, a;
    while(maxAbsValue(R)>epsilon && iter < maxiter) {
        
        P = CombLin(R,P,1,b);//aixo surt al final del loop a la wikipedia pero ho començo ja aqui amb b = 0;

        a = dotP(R,R)/dotP(P,matVecMult2DGrid(P,CoefN,CoefE,CoefS,CoefW,CoefP,M,N));
        b = dotP(R,R); //Poso ja el R_n al calcular ara R_{n+1}
        PressureVec = CombLin(PressureVec,P,1,a);
        R = CombLin(R,matVecMult2DGrid(P,CoefN,CoefE,CoefS,CoefW,CoefP,M,N),1,-a);
        b = dotP(R,R)/b;          //fent aixi el cocient evito possibles
                                    //errors de fer 1/a(n) primer 
                                    //he de posar aqui la b, no com la P, perque
                                    //fent aixo ja deixa de ser 0 la b i em fotria el 
                                    //primer loop

                                    
        iter++;
    }
    if(k%20 == 0){
        std::cout << "Timestep " << k << " iter:" << iter << " itermax:" << maxiter << std::endl;
    }
}    



void ErrorPrint(std::vector<std::vector<double>> Pressure, double M, double N, double rho, double dx, double dt, std::vector<std::vector<double>> &PredU
    , std::vector<std::vector<double>> &PredV, double &PressureError, double q, int k, double &StationaryError){


        double C,K;

        //Aqui l'error de tot aixo
        for(int i = 2; i < M+1 ; i++){
            for(int j = 2; j < N+1; j++){
                C = (((i!=1)*Pressure[i-1][j] + (i!=M)*Pressure[i+1][j])/q 
                + ((j!=N)*Pressure[i][j+1] + (j!=1)*Pressure[i][j-1])*q);
                K = rho*dx*(PredV[i][j] + PredU[i][j+1]*q - PredV[i+1][j] - PredU[i][j]*q)/dt;

                if(PressureError < fabs(Ppcoef(N,i,j,q,M)*Pressure[i][j] - C + K)){
                    PressureError = fabs(Ppcoef(N,i,j,q,M)*Pressure[i][j] - C + K);
                }
            }
        }
        if(k%20==0){
            std::cout << "Pressure Error in timestep " << k << " is: " << PressureError << std::endl;
            std::cout << "Stationary Error in timestep " << k << " is: " << StationaryError << std::endl;
        }
    }



void TempChange(std::vector<std::vector<double>> &Temp, double N, double M, std::vector<std::vector<double>> &MatV, std::vector<std::vector<double>> &MatU
    , std::vector<std::vector<double>> &nuTemp, double dt, double lambda, double dy, double dx, double rho, double cp
    , double &StationaryError){

    double Tn,Te,Ts,Tw,vn,ue,vs,uw;
    double q = dy/dx;
    double Delta;

    Temp = nuTemp;

    for(int i = 1; i < M+1; i++){
        for(int j = 1; j < N+1; j++){
            Tn = avg(Temp[i-1][j],Temp[i][j]); 
            Te = avg(Temp[i][j+1],Temp[i][j]);
            Ts = avg(Temp[i+1][j],Temp[i][j]); 
            Tw = avg(Temp[i][j-1],Temp[i][j]);          
            vn = MatV[i][j];
            ue = MatU[i][j+1];
            vs = MatV[i+1][j];
            uw = MatU[i][j];

            Delta = dt*(lambda*(Temp[i-1][j]*(1+(i==1))/q
                            + Temp[i][j+1]*(1+(j==N))*q
                            + Temp[i+1][j]*(1+(i==N))/q
                            + Temp[i][j-1]*(1+(j==1))*q 
                            - Ppcoefreborn(N,i,j,q,M)*Temp[i][j])/(rho*cp*dx) 
                            - (Tn*vn - Ts*vs + (Te*uw - Tw*uw)*q))/dy;

            nuTemp[i][j] = Temp[i][j] + Delta;

            if(StationaryError < fabs(Delta)){
                StationaryError = fabs(Delta);
            }  
        }
    }

    //ara parets adiabatiques adalt i abaix
    //ho poso al nuTemp
    //posem aixo aqui amb el nutemp perque ara el Temp = nuTemp esta
    //abans del loop, per conservar Temp^n i Temp^n+1 dels timesteps.
    for(int i = 1; i < N+1;i++){
        nuTemp[0][i] = nuTemp[1][i];  
        nuTemp[M+1][i] = nuTemp[M][i]; 
    }

}

void NewVelocity(std::vector<std::vector<double>> &MatV, std::vector<std::vector<double>> &MatU, std::vector<std::vector<double>> &PredU
    , std::vector<std::vector<double>> &PredV, double &Vmax, double &vmax, double &Umax, double &umax, int k, double dt, int M
    , int N, double rho, std::vector<std::vector<double>> &Pressure, double L, int TSteps, double &StationaryError
    , std::vector<std::vector<double>>& MatPreU, std::vector<std::vector<double>>& MatPreV){

    double dx = L/(N-1), dy = L/(M-1);  
    double I,J;  

    for(int i = 2; i < M+1; i++){
        for(int j = 1; j < N+1; j++){
            MatV[i][j] = PredV[i][j] - dt*(Pressure[i-1][j]-Pressure[i][j])/(rho*dy);
            if(StationaryError < fabs(MatV[i][j]-MatPreV[i][j])){
                StationaryError = fabs(MatV[i][j]-MatPreV[i][j]);
            }
            if(k==TSteps && fabs(MatV[i][j])>Vmax && i == 1+floor(M/2)){
                Vmax = fabs(MatV[i][j]);
                J = (double)(j-1)*dx; 
            }
            if(fabs(MatV[i][j])>vmax){
                vmax = MatV[i][j];
                //J = (double)(i-0.5)*dx; 
                if(fabs(vmax*dt/dy) > 0.35){
                    std::cout << "One of the convergence criteria fails, v*dt/dy too large" << std::endl;
                    std::cout << "vmax: " << vmax << " k: " << k << std::endl;
                    return;
                }
            }
        }
    } 
    for(int i = 1; i < M+1; i++){
        for(int j = 2; j < N+1; j++){
            MatU[i][j] = PredU[i][j] - dt*(Pressure[i][j]-Pressure[i][j-1])/(rho*dx);
            if(StationaryError < fabs(MatU[i][j]-MatPreU[i][j])){
                StationaryError = fabs(MatU[i][j]-MatPreU[i][j]);
            }
            if(k==TSteps && fabs(MatU[i][j])>Umax && j == 1+floor(N/2)){
                Umax = fabs(MatU[i][j]);
                I = L-(double)(i-1)*dy;
            }
            if(fabs(MatU[i][j])>umax){
                umax = MatU[i][j];
                if(fabs(umax*dt/dx) > 0.35){
                    std::cout << "One of the convergence criteria fails, v*dt/dx too large" << std::endl;
                    std::cout << "umax: " << umax << " k: " << k << std::endl;
                    return;
                }
            }
        }
    }

    if(k%20 == 0){
        double sum_sq = 0;
        for (int i = 1; i < M+1; i++) {
            for (int j = 1; j < N+1; j++) {
                double div = (MatU[i][j+1] - MatU[i][j]) * dy
                            + (MatV[i][j] - MatV[i+1][j]) * dx;//Ara si, confiem en el Pep
                sum_sq += div * div;
            }
        }
        double rms_div = sqrt(sum_sq / (N*M));
        std::cout << " RMS divergence: " << rms_div << std::endl; 
    }
  
    if(k%200 == 0){
        std::cout << "k: " << k << " ,vmax : " << vmax << std::endl;

    }
}

void GaussSeidel(int M, int N, double dy, double dx, double rho, double dt 
, std::vector<std::vector<double>> &Pressure, std::vector<std::vector<double>>& PredU
, std::vector<std::vector<double>>& PredV, int maxiter, double Omega, double epsilon
, double Error, bool B, std::vector<double>& CoefP, std::vector<double>& CoefN
, std::vector<double>& CoefE, std::vector<double>& CoefS, std::vector<double>& CoefW,int k){
//Working on this shit still, time to see if my conjugate convergence works tho
    double P, q = dy/dx;
    int iter = 0;
    double C,K;

    do{
        B = 1;
        for(int i = 1; i < M+1; i++){
            for(int j = 1; j < N+1; j++){
                if (i == 1 && j == 1){
                    continue; //skipping the fixed point since we want to keep it at 10.
                }
                C = (CoefN[idx(i-1,j-1,N)]*Pressure[i-1][j] + CoefS[idx(i-1,j-1,N)]*Pressure[i+1][j] 
                + CoefE[idx(i-1,j-1,N)]*Pressure[i][j+1] + CoefW[idx(i-1,j-1,N)]*Pressure[i][j-1]);
                K = rho*dx*(PredV[i][j] + PredU[i][j+1]*q 
                - PredV[i+1][j] - PredU[i][j]*q)/dt;  //com sigui la merda aquesta que no se com estava malament (segueix malament)
                
                P = (-C + K)/CoefP[idx(i-1,j-1,N)]; //He multiplicat adalt per -1, perque crec que el
                                                    //CoefP esta en negatiu ja pq aixi la suma de coefs
                                                    //a la matriu a cada file donava 0 
                B = B && (fabs(Pressure[i][j]-P)<epsilon);
                Pressure[i][j] = (1-Omega)*Pressure[i][j]+Omega*P;
           
            }                                                                                   
        }
        //now to fill in the boundary condition dP/dn = 0.
        //Due to these rows and columns being identical to their adjacent ones
        //they are already checked in the |Pressure[i][j]-P|<epsilon
        for(int i = 1; i < N+1 ; i++){
            Pressure[0][i] = Pressure[1][i];
            Pressure[M+1][i] = Pressure[M][i];
        }
        for(int i = 1; i < M+1; i++){
            Pressure[i][0] = Pressure[i][1];
            Pressure[i][N+1] = Pressure[i][N];
        }
        iter++;
    }while(!B && iter < maxiter); 

    if(k%20 == 0){
        std::cout << "Timestep " << k << " iter:" << iter << " itermax:" << maxiter << std::endl;
    }
}


