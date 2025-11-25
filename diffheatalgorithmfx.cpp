#include <iostream>
#include <fstream>
#include <math.h> // To use pow
#include <vector> // To use vectors
#include <chrono>
#include <iomanip>
#include <string> 

#include "coefs_stokes.h"
#include "operations_vectors.h"
#include "mainfx.h"
#include "matrix_print.h"
using namespace std;


int main(void){

    double L=1.0, rho = 1, mu = 0.71;
    double dt, t, TSteps;
    int iter=0, maxiter=30000;
    double pun, pue, pus, puw;
    double puN, puE, puS, puW, puP;
    double pvn, pve, pvs, pvw;
    double pvN, pvE, pvS, pvW, pvP;
    double epsilon = 0.000001,StationaryEpsilon = 0.00001;
    double C, K, Press;
    double Tn,Te,Ts,Tw;
    double Omega = 1.9;
    double cp=1,g,beta=1,Thot=1,Tcold=0,Tinf=0.5,lambda=1,Ra;
    double PressureError = 0, StationaryError = 0;
    double I,J;
    double vmax,umax,Umax,Vmax;
    int M,N;
    //Aqui hi ha unes variables adimensaionals que son les tipiques de x/L, T' = T-Tcold/(Thot-Tcold), uL/alpha
    //Llavors, com que per mi la alpha es 1 per l'eleccio de mu = 0.71, en lloc de cp = 0.71, L = 1, i Thot = 1,
    //Tcold = 0, totes les variables d'aqui ja estan en aqusta base, i haurien de quadrar amb les taules donades
    //i tal.
    
    




    
    std::cout << "Enter the value of the Rayleigh number: " << std::endl;
    cin >> Ra;
    if(Ra!=1000 && Ra!=10000 && Ra!= 100000 && Ra!= 1000000){
        return 1;
    }
    g = Ra*0.71;
    std::cout << "With the given values, the Prandtl number is Pr = " << mu*cp/lambda << std::endl;
    std::cout << "The gravitational acceleration in this planet is g = " << g << std::endl;

    t = 10;
    std::cout << "Duration of the simulation of the system in seconds unless convergence is achieved earlier is: " << t << " seconds." << std::endl;
    

    std::cout << "Enter the number of horizontal nodes, N (must be an integer): ";
    cin >> N;
    std::cout << "Enter the number of vertical nodes, M (must be an integer): ";
    cin >> M; 
    double dx = L/(N-1); //millor al final no introdueixo un dy...
    double dy = L/(M-1); //ara que pringat
    double q = dy/dx;
    double dV = dx*dy;

    
    dt = 0.001;
    if(dt > 0.35*dx){
        dt = 0.20*dx;
    }
    if(dt > 0.35*dy){
        dt = 0.20*dy;
    }
    if(dt > 0.2*rho*dx*dx/mu){
        dt = 0.12*rho*dx*dx/mu;
    }
    if(dt > 0.2*rho*dy*dy/mu){
        dt = 0.12*rho*dy*dy/mu;
    } 
    //CFL (Courant–Friedrichs–Lewy) condition – convective stability and
    //Diffusion (viscous) stability criterion
    //no se segur si es dx*dy ara, pero em sembla probable, es veu que es min(dx,dy)
    std::cout << "convective stability and CFL conditions C_conv = dt/dx, C_diff = mu*dt/(rho*dx*dx), and their dy equivalents yield dt = " << dt << std::endl;
       
    TSteps = floor(t/dt);
    std::cout << "Time steps to be carried out unless stationary state is reached before: " << TSteps << std::endl;



    using namespace std::chrono;
    // Start global timer
    auto start = high_resolution_clock::now();


    std::vector<double> CoefN(M*N,0);
    std::vector<double> CoefE(M*N,0);
    std::vector<double> CoefS(M*N,0);
    std::vector<double> CoefW(M*N,0);
    std::vector<double> CoefP(M*N,0);
    getCoefs(CoefN,CoefE,CoefS,CoefW,CoefP,M,N,q);
    //print = matVecMult2DGrid(test,CoefN,CoefE,CoefS,CoefW,CoefP,M,N);
    //Sembla que aquests dos van be
    


    //We will talk about the matrix for the horizontal velocity, MatU, to explain the idea.
    //Later on, we will treat the nodes that touch the edges of the material separately, as we know
    //that velocity is 0 there, at the same time, we need a way to treat the other edge nodes that
    //aren't placed at the edges of the material, as they have a single virtual neighbour outside the
    //material, with again, velocity 0. So the extra 2 nodes will be this virtual velocity = 0.
    std::vector<std::vector<double>> MatU(M+2, std::vector<double>(N+3, 0));   //M+2 Rows, N+3 Columns
    std::vector<std::vector<double>> MatV(M+3, std::vector<double>(N+2, 0));

    //previous timestep velocities
    std::vector<std::vector<double>> RU(M+2, std::vector<double>(N+3, 0));   //M+2 Rows, N+3 Columns
    std::vector<std::vector<double>> RV(M+3, std::vector<double>(N+2, 0));   
    std::vector<std::vector<double>> MatPreU(M+2, std::vector<double>(N+3, 0));   //M+2 Rows, N+3 Columns
    std::vector<std::vector<double>> MatPreV(M+3, std::vector<double>(N+2, 0));   
    

    //predictor velocity matrices
    std::vector<std::vector<double>> PredU(M+2, std::vector<double>(N+3, 0));    //M+2 rows, N+3 columns
    std::vector<std::vector<double>> PredV(M+3, std::vector<double>(N+2, 0)); 

    //Pressure matrix and copy matrix for Gauss-Seidel
    //two more rows and columns for the extra node equal to 0 for the boundary nodes because of dP/dn = 0.
    std::vector<std::vector<double>> Pressure(M+2, std::vector<double>(N+2, 1));
    //temperature matrices
    //Due to the similarities between a particular equation in the energy equation we will also employ
    //two meshes for the temperature.
    std::vector<std::vector<double>> Temp(M+2, std::vector<double>(N+2,0));//M+2 rows 
    std::vector<std::vector<double>> nuTemp(M+2, std::vector<double>(N+2,0));//because I'm doing temperature explicitly

    for(int i = 1; i < M+2; i++){
        Temp[i][0] = 1; 
        Temp[i][N+1] = 0;
        nuTemp[i][0] = 1;
        nuTemp[i][N+1] = 0;
    }

    //For the recursive algorithm
    // std::vector<double> R(M*N, 0);  // Size M*N vector
    // std::vector<double> P(M*N, 0);  // Size M*N vector 2 
    // double a,b;
    // std::vector<double> PressureVec(M*N,0);
    // std::vector<double> ct(M*N,0);
    // std::vector<double> dummy(M*N,0);


    //sincerament els noms son mes o menys adaptats de les lletres usades a la wikipedia 
    //per aquest algoritme

    auto t1 = high_resolution_clock::now();

    std::cout << "q: " << q << " ,1/q: " << 1/q << " ,dy: " << dy << " ,dx: " << dx << std::endl;


    for(int k = 1; k < TSteps + 1; k++){
        

        //Step 1//
        //Compute the predictor velocity matrices to use afterwards//
        //Aqui canvia el primer step al afegir temperatura pel boussinesq term
        //i tambe interessa conservar que div(v^p) = dt/rho * nabla^2 P^n+1, en plan mantenir el step 2 igual
        //tambe mante el step 3 igual

        //si no m'equivoco, com al principi totes les velocitats son 0, no cal que calculi primer
        //RU ni RV, perque son 0 igualment, i ja les he inicialitzat a 0.
        
        
        PredUfill(M,N, MatU, MatV, PredU, RU, rho, dx, dy,dt,mu);
        PredVfill(M,N, MatU, MatV, PredV, RV, Temp, nuTemp, dt, beta, g, rho, dx,dy,mu);
        

        //Step 2//
        //Gauss-Seidel to solve for the pressure at each node//
        bool B = 1;

        Pressure[1][1] = 1.0;  //Fixing pressure at this point because in this version of NS
                                //with incompressible flow, pressure is only determined up to a constant
                                //and otherwise Gauss Seidel may diverge (it is in fact diverging)
        Pressure[0][1] = 1.0;
        Pressure[1][0] = 1.0;//Neumann boundary condition implies these pressures are also 0
        //Amb el continue dins el for puc treure aixo del do while entenc

        GaussSeidel(M,N,dy,dx,rho,dt,Pressure,PredU,PredV,maxiter,Omega,epsilon,StationaryError,B,CoefP
        ,CoefN,CoefE,CoefS,CoefW,k);
        //Here used to lay conjugate gradient for an ephemeral instant
        
        


        
        //VEIEM SI TENIR AIXO FORA DEL WHILE EM DONA ALGUN PROBLEMA (CG)
        //now to fill in the boundary condition dP/dn = 0.
        //Due to these rows and columns being identical to their adjacent ones
        //they are already checked in the |Pressure[i][j]-P|<epsilon
        

        //gotta put the good values into the matrix after going through the iterations
        //with the vector
        // for(int i = 1; i < M+1; i++){
        //     for(int j = 1; j < N+1; j++){
        //         Pressure[i][j] = PressureVec[idx(i-1,j-1,N)];
        //     }
        // }

        
        if(k%20==0){
            ErrorPrint(Pressure, M, N, rho, dx, dt, PredU, PredV, PressureError, q,k,StationaryError);
            StationaryError = 0;
        }
        //Step 3//
        //Energy equation, temperature at next timestep//

        TempChange(Temp,N,M,MatV,MatU,nuTemp,dt,lambda,dy,dx,rho,cp,StationaryError);

        //Step 4//
        //Velocity at next time step//

        //saving whats about to be the previous timestep into the matrix.
        vmax = 0;
        umax = 0;
        MatPreU = MatU;
        MatPreV = MatV;

        NewVelocity(MatV,MatU,PredU,PredV,Vmax,vmax,Umax,umax,k,dt,M,N,rho,Pressure,L,TSteps,StationaryError
                    ,MatPreU,MatPreV);

        //std::cout << "k: " << k << " vmax: " << vmax << " umax: " << umax << std::endl;

        
        if(StationaryError < StationaryEpsilon && k > 1000){
            std::cout << "Stationary state reached. It took: " << k*dt << " seconds to reach it." << std::endl;
            break;
        }//vale per ara ho deixo pero el pep deia que aquest error havia de ser amb
         //les velocitats no la pressio de nou
        //commenting it for now because I want the program to run for a while just in case
        //to see if anything wonky happens
        //Aixo amb velocitat no Pressio
    } 

    //Temp = nuTemp;//cal fer aixo al posar el Temp = nuTemp abans del loop aixi que encara li queda
                    //la ultima actualitzacio a Temp
                    // veient el que ve despres ho trec i poso 
                    // simplement nuTemp en lloc de Temp


    //matrix modifications to not have the virtual nodes    (?)
    std::vector<std::vector<double>> UV(2*M+2, std::vector<double>(2*N+2, 0)); //matrix for both velocities 
                                                                               //Same size as for temperature
    



    for(int i=1; i < M+2; i++){
        for(int j = 1; j < N+2; j++){
            UV[2*i-1][2*(j-1)] = MatU[i][j]; //in UV the coordinates (odd,even) are u velocity
            UV[2*(i-1)][2*j-1] = MatV[i][j]; //in UV the coordinates (even,odd) are v velocity
        }                                    //incidentally (odd,odd) would be pressure
    }


    //computing the flux thing//
    double Sum = 1; //suma telescopica hauria de cancelarse tot per molt malament q estigui
    //Sum = MatU[1][1]/4+MatU[N][1]/4; //Aqui la temperatura es 1, pero la velocitat 0 lol
    // for(int i=2;i<N;i++){
    //     Sum += MatU[i][0]/2; //Aqui tambe es 1, aixo es paret vertical esquerra
    // }
    // for(int i = 2; i < N+2 ;i++){
    //     Sum += MatU[1][i]*Temp[1][i]+;   //paret adiabatica aixi que no fem mitjanes
    // }
    for(int i = 1; i < M+1; i++){
        for(int j = 1; j < N+1; j++){
            Sum += nuTemp[i][j]*avg(MatU[i][j],MatU[i+1][j])*dx*dy;
        }
    }

    writeInnerMatrixToFile(MatU, "MatU.txt");
    writeInnerMatrixToFile(MatV, "MatV.txt");
    writeInnerMatrixToFile(Pressure, "Pressure.txt");
    write_centralcolumn(MatU,N,dx);
    write_centralfile(MatV,N,dx);
    write_full_matrix_to_file(UV);
    write_full_matrix_to_fileTemp(nuTemp);
    std::cout << "Matrix UV written to 'UV.txt'" << std::endl;
    std::cout << "Average Nusselt number: " << Sum << std::endl;
    std::cout << "Maximum U velocity is " << Umax << " at height " << I << std::endl;
    std::cout << "Maximum V velocity is " << Vmax << " at length " << J << std::endl;

    Vmax = 0;
    Umax = 0;
    for(int i = 1; i < N+1; i++){
        if(MatV[1+floor(M/2)][i] > Vmax){
            Vmax = MatV[floor(M/2)][i];
            J = (double)(i-0.5)*dx;
        }
    }
    for(int i = 1;i< M+1; i++){
        if(MatU[i][1+floor(N/2)] > Umax){
            Umax = MatU[i][floor(N/2)];
            I = L-(double)(i-0.5)*dx;
        }
    }
    std::cout << "Maximum U velocity with better method is " << Umax << " at height " << I << std::endl;
    std::cout << "Maximum V velocity with better method is " << Vmax << " at length " << J << std::endl;
    auto t5 = high_resolution_clock::now();
    std::cout << "The entire program took " << duration_cast<seconds>(t5 - t1).count() << " s\n";

    return 0;
}

