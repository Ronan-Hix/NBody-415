#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nr.h"
#include "nrutil.h"


#define nbodies 2
#define noutputs 1000
#define eps 0

void f(double t,double *y, double *dydt){
    double dx,dy,dz,d;
    for(int i=0;i<nbodies;i++){
        dydt[1+7*i] = y[4+7*i];
        dydt[2+7*i] = y[5+7*i];
        dydt[3+7*i] = y[6+7*i];
        dydt[4+7*i] = 0;
        dydt[5+7*i] = 0;
        dydt[6+7*i] = 0;
        for(int j=0;j<nbodies;j++){
            if (i != j){ //Omits i==j
                dx = y[1+7*i]-y[1+7*j];
                dy = y[2+7*i]-y[2+7*j];
                dz = y[3+7*i]-y[3+7*j];
                d = dx*dx+dy*dy+dz*dz;
                dydt[4+7*i] -= y[0+7*j]*y[0+7*i]*dx/pow(d+eps,1.5);
                dydt[5+7*i] -= y[0+7*j]*y[0+7*i]*dy/pow(d+eps,1.5);
                dydt[6+7*i] -= y[0+7*j]*y[0+7*i]*dz/pow(d+eps,1.5);
            }
        }
    }
}

void frk(float t,float *y, float *dydt){
    float dx,dy,dz,d;
    for(int i=0;i<nbodies;i++){
        dydt[1+6*i] = y[4+6*i];
        dydt[2+6*i] = y[5+6*i];
        dydt[3+6*i] = y[6+6*i];
        dydt[4+6*i] = 0;
        dydt[5+6*i] = 0;
        dydt[6+6*i] = 0;
        for(int j=0;j<nbodies;j++){
            if (i != j){ //Omits i==j
                dx = y[1+6*i]-y[1+6*j];
                dy = y[2+6*i]-y[2+6*j];
                dz = y[3+6*i]-y[3+6*j];
                d = dx*dx+dy*dy+dz*dz;
                dydt[4+6*i] -= dx/pow(d+eps,1.5);
                dydt[5+6*i] -= dy/pow(d+eps,1.5);
                dydt[6+6*i] -= dz/pow(d+eps,1.5);
            }
        }
    }
}


void energy(double *y,double *E){
    double d,vsq,dx,dy,dz;
    *E = 0;
    for(int i=0;i<nbodies;i++){
        vsq = y[4+7*i]*y[4+7*i]+y[5+7*i]*y[5+7*i]+y[6+7*i]*y[6+7*i];
        *E += 0.5*vsq;
        for(int j=0;j<nbodies;j++){
            if (i < j){
                dx = y[1+7*i]-y[1+7*j];
                dy = y[2+7*i]-y[2+7*j];
                dz = y[3+7*i]-y[3+7*j];
                d = dx*dx+dy*dy+dz*dz;
                *E -= 1/sqrt(d);
            }
        }
    }
}

void energyrk(float *y,double *E){
    float d,vsq,dx,dy,dz;
    *E = 0;
    for(int i=0;i<nbodies;i++){
        vsq = y[4+6*i]*y[4+6*i]+y[5+6*i]*y[5+6*i]+y[6+6*i]*y[6+6*i];
        *E += 0.5*vsq;
        for(int j=0;j<nbodies;j++){
            if (i < j){
                dx = y[1+6*i]-y[1+6*j];
                dy = y[2+6*i]-y[2+6*j];
                dz = y[3+6*i]-y[3+6*j];
                d = dx*dx+dy*dy+dz*dz;
                *E -= 1/sqrt(d);
            }
        }
    }
}

void leapfrog(double tn, double *y,double *dydt, double E,int N, double h, double nstepmax, double *yout,void (*f)(double, double [], double [])){
//    double vn1_2 = y[1]*(0.5*h);//Kick from tn=0
//    double xn = y[0];
    int i;
    f(tn,y,dydt);
    for(i=0;i<nbodies;i++){ //Initialize
        yout[0+7*i] = y[0+7*i];
        yout[1+7*i] = y[1+7*i];
        yout[2+7*i] = y[2+7*i];
        yout[3+7*i] = y[3+7*i];
        yout[4+7*i] = y[4+7*i];
        yout[5+7*i] = y[5+7*i];
        yout[6+7*i] = y[6+7*i];
    }
    f(tn,yout,dydt);
    for(i=0;i<nbodies;i++){
        yout[4+7*i] += (0.5*h)*dydt[4+7*i]; //Kick
        yout[5+7*i] += (0.5*h)*dydt[5+7*i];
        yout[6+7*i] += (0.5*h)*dydt[6+7*i];
    }
    for(i=0;i<nbodies;i++){
        yout[1+7*i] += h*yout[4+7*i]; //Drift
        yout[2+7*i] += h*yout[5+7*i];
        yout[3+7*i] += h*yout[6+7*i];
    }
    
    f(tn,yout,dydt); //f(xn+1)
    for(i=0;i<nbodies;i++){
        yout[4+7*i] += (0.5*h)*dydt[4+7*i]; //Resync
        yout[5+7*i] += (0.5*h)*dydt[5+7*i];
        yout[6+7*i] += (0.5*h)*dydt[6+7*i];
    }
    //printf("%f,%f,%f,%f,%f,%f,%f\n",h,tn,yout[0],yout[1],yout[2],yout[3],E);
    
}

int main(int nargs, char **argv){
    int nstepmax, outputfreq, n,m,i,j;
    double t,h,tmax,E;
    double x[7*nbodies],dxdt[7*nbodies],yout[7*nbodies];
    float trk,hrk,tmaxrk;
    float mass;
    float *xrk,*dxdtrk,*youtrk;
    char fileoutname[100];
    
    FILE *data;
    data = fopen(argv[1],"r+");
    

    if (nargs != 5){
        printf("Error: Wrong Number of arguments. Please supply, 'file path', 'method', 'tmax' and 'h' \n");
        exit(-1);
    }
    
    if (strcmp("rk4",argv[2])==0){
        xrk=vector(1,6*nbodies);
        dxdtrk=vector(1,6*nbodies);
        youtrk=vector(1,6*nbodies);
        hrk = atof(argv[4]);
        trk = 0;
        tmaxrk=atof(argv[3]);
        nstepmax = (tmaxrk/hrk);
        for (n=0; n < nbodies; n++){
            fscanf(data, "%f, %f, %f, %f, %f, %f, %f",&mass,&xrk[1+6*n],&xrk[2+6*n],&xrk[3+6*n],&xrk[4+6*n],&xrk[5+6*n],&xrk[6+6*n]);
        }
    }

    else{
        h = atof(argv[4]);
        t = 0;
        tmax=atof(argv[3]);
        nstepmax = (tmax/h);
        for (n=0; n < nbodies; n++){
            for (m=0; m < 7 && fscanf(data, "%lf,",&x[m+7*n]) != EOF; m++);
        }
    }

    outputfreq = nstepmax/noutputs;
    
    if (strcmp("rk4",argv[2])==0){
        for(n = 0;n<=nstepmax;n++){
            (*frk)(trk,xrk,dxdtrk);
            rk4(xrk,dxdtrk,6*nbodies,trk,hrk,youtrk,frk);
            //printf("[%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",xrk[1],xrk[2],xrk[3],xrk[4],xrk[5],xrk[6],xrk[7],xrk[8],xrk[9]);
            trk += hrk;
            for(i=0;i<nbodies;i++){ //Initialize
                xrk[1+6*i] = youtrk[1+6*i];
                xrk[2+6*i] = youtrk[2+6*i];
                xrk[3+6*i] = youtrk[3+6*i];
                xrk[4+6*i] = youtrk[4+6*i];
                xrk[5+6*i] = youtrk[5+6*i];
                xrk[6+6*i] = youtrk[6+6*i];
            }
            energyrk(youtrk,&E);
            if(n%outputfreq == 0){
                
                sprintf(fileoutname,"outputs/out_%04d.dat",n/outputfreq);
                FILE *fileout;
                fileout = fopen(fileoutname, "w+");
                /*
                printf("%f,%f",hrk,trk);
                for(i=0;i<nbodies;i++){
                    printf(",%d",1);
                    for(j=1;j<7;j++){
                        printf(",%f",youtrk[j+6*i]);
                    }
                    
                }
                printf(",%f\n",E);
                */
                for(i=0;i<nbodies;i++){
                    //fprintf(fileout,"%f,%f",h,t);
                    for(j=1;j<7;j++){
                        fprintf(fileout,"%f,",yout[6*i+j]);
                    }
                    fprintf(fileout,"\n");
                    //fprintf(fileout,",%f\n",E);
                }
                
                fclose(fileout);
                
            }
        }
    }
    else if (strcmp("leapfrog",argv[2])==0){
        for(n =0;n<=nstepmax;n++){
            //printf("[%f,%f,%f,%f,%f,%f,%f,%f]\n",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]);
            leapfrog(t,x,dxdt,E,4,h,nstepmax,yout,&f);
            t += h;
            for(i=0;i<nbodies;i++){ //Initialize
                x[1+7*i] = yout[1+7*i];
                x[2+7*i] = yout[2+7*i];
                x[3+7*i] = yout[3+7*i];
                x[4+7*i] = yout[4+7*i];
                x[5+7*i] = yout[5+7*i];
                x[6+7*i] = yout[6+7*i];
            }
            energy(yout,&E);
            
            if(n%outputfreq == 0){
                
                /*
                printf("%f,%f",h,t);
                for(i=0;i<7*nbodies;i++){
                    printf(",%f",yout[i]);
                }
                printf(",%f\n",E);
                */
                sprintf(fileoutname,"outputs/out_%04d.dat",n/outputfreq);
                FILE *fileout;
                fileout = fopen(fileoutname, "w+");
                
                for(i=0;i<nbodies;i++){
                    //fprintf(fileout,"%f,%f",h,t);
                    for(j=0;j<7;j++){
                        fprintf(fileout,"%lf,",yout[7*i+j]);
                    }
                    fprintf(fileout,"\n");
                    //fprintf(fileout,",%f\n",E);
                }
                
                fclose(fileout);
            }
             
        }
    }
    else{
        printf("Error: '%s' is not a valid method name. Use euler, rk4, or leapfrog\n",argv[1]);
        exit(-1);
    }
}
