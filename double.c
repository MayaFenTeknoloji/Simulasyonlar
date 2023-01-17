#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <quadmath.h>
int FCOUNT; 

void mekanik_enerji(__float128 r[3],__float128 v[3], __float128 *E)  {
    __float128 mu= 2.9619474286664206E-04q;
    __float128 v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    __float128 r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    __float128 norm_r=sqrtq(r2);
    __float128 ek = 0.5*v2;         // kinetic energy 
    __float128 epot = -mu/norm_r;    // potential energy
          *E = ek + epot;       // mechanical energy                           
    } 

// equation of motion's right-hand side evaluation as function:
void f(__float128 r[3], __float128 v[3], __float128 *r_i, __float128 *v_i)  {
     __float128 mu= 2.9619474286664206E-04q;
     __float128 dot = r[0]*r[0]+ r[1]*r[1]+ r[2]*r[2];
     __float128 r3=dot*sqrtq(dot); 

        for(int i=0; i<3; i++)                {
         v_i[i] = -r[i]*mu/r3; 
         r_i[i] = v[i];                       }               
                                              }    
                                          
                                                           
// 4th Order Runge Kutta Method with quadruple-precision:
void RK4(__float128 r[3], __float128 v[3], __float128 dt) {
    __float128 k_r1[3], k_v1[3], k_r2[3], k_v2[3];
    __float128 k_r3[3], k_v3[3], k_r4[3], k_v4[3];
    __float128 v_temp1[3], v_temp2[3], v_temp3[3];
    __float128 r_temp1[3], r_temp2[3], r_temp3[3];
    __float128 kr[3], kv[3];

    f(r, v, k_r1, k_v1);
        for(int i=0; i<3; i++)               {
            v_temp1[i]= v[i]+k_v1[i]*0.5*dt;
            r_temp1[i]= r[i]+k_r1[i]*0.5*dt; 
                                             }

    f(r_temp1, v_temp1, k_r2, k_v2);
        for(int i=0; i<3; i++)               {
            v_temp2[i]= v[i]+k_v2[i]*0.5*dt;
            r_temp2[i]= r[i]+k_r2[i]*0.5*dt;
                                             }

    f(r_temp2, v_temp2, k_r3, k_v3);
        for(int i=0; i<3; i++)               {
            v_temp3[i]= v[i]+k_v3[i]*dt;
            r_temp3[i]= r[i]+k_r3[i]*dt;    
                                             }

    f(r_temp3, v_temp3, k_r4, k_v4);

        for(int i=0; i<3; i++)               {
           kv[i] = (k_v1[i] +2*(k_v2[i]+k_v3[i])+ k_v4[i])/6.0;
           kr[i] = (k_r1[i] +2*(k_r2[i]+k_r3[i])+ k_r4[i])/6.0;
           v[i] += kv[i]*dt;
           r[i] += kr[i]*dt;                 }
    }  
                                  

int main() {

    double t, dt, T;
     __float128 energy;     
     __float128 r[]= {6.931720357474623E+08q, -2.745889195466235E+08q, -1.436878793848155E+07q}; // position
     __float128 v[]= {4.653807892682983E+00q, 1.276054417367929E+01q, -1.571074085212407E-01q}; // velocity
     __float128 E= 92.25704912411518333633168582252779905q;                                      // mechanical energy of system

      char buf1[128], buf2[128];
    //  quadmath_snprintf (buf1, sizeof buf1, "%*.33Qf", 20, E);
    
    T=  4500;                                
    dt= 0.01;
    t=0.0 ;   

    double time_spent = 0.0;
    clock_t begin = clock(); 
    int Count=0.0;

    FILE *tt;
    tt = fopen("RK4_quadruple_hata.txt","w");
    while ( t <= T )  {          // ( t <= tmax/dt )
       RK4(r, v, dt);  
       mekanik_enerji(r, v, &energy);   


        __float128 energy_1= (E-energy)/E;
        quadmath_snprintf (buf2, sizeof buf2, "%*.35Qf", 10, energy_1);
 
        fprintf(tt,"%.3lf %s  \n", t, buf2);
        t +=dt;       
        Count++; }

    clock_t end = clock();
    fclose(tt);

    time_spent += (double)(end - begin) / CLOCKS_PER_SEC; 
 
    printf("The elapsed time is %f seconds \n", time_spent);
    return 0;  }
    //gcc -o quadruple.exe quadruple.c -lm -lquadmath 
    //  0.203000
 
