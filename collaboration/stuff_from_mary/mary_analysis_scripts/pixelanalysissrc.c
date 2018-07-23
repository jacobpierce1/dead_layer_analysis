

#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TDatime.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TGraphErrors.h"
#include "selectionSort.h"
#include "swap.h"

void pixelanalysissrc()
    {
        //will create a graph using text files from runs documenting error and line center of each of 3 peaks. Shows error bar vs length of data taking for a strip and error bar over 3 different runs.
        gStyle->SetLabelSize(0.047,"y");
        gStyle->SetLabelSize(0.047,"x");
        gStyle->SetTitleFontSize(0.2);
        
        FILE *CENTER, *RIGHT, *COS;
        char center[50],right[50],canvasname[50],graphname[50];
        int i, r[900], pixel, A,j,k,z=0,y=0,can=0,q=0;
        float pu240ec[900], pu240er[900], pu240erc[900], pu240err[900], pu240sc[900], pu240sr[900], pu240serc[900], pu240serr[900], pu240cc[900], pu240cr[900];
        
        float pu238ec[900], pu238er[900], pu238erc[900], pu238err[900], pu238sc[900], pu238sr[900], pu238serc[900], pu238serr[900], pu238cc[900], pu238cr[900];
        
        float cf249ec[900], cf249er[900], cf249erc[900], cf249err[900], cf249sc[900], cf249sr[900], cf249serc[900], cf249serr[900], cf249cc[900], cf249cr[900];
        
        float a, L[900], dEdx[3], dEdxe[3], Epu240, Epu238, Ecf249, Cpu[900], Ccf[900], Ccent[900], Crt[900], Cpuerr[900], Ccferr[900], Ccenterr[900], Crterr[900], pu240x, pu240y, pu238cx, pu238cy, pu238rx, pu238ry, cf249x, cf249y, pu240xerr, pu240yerr, pu238cxerr, pu238cyerr, pu238rxerr, pu238ryerr, cf249xerr, cf249yerr,Pu240x[900], Pu240y[900], Pu238cx[900], Pu238cy[900], Pu238rx[900], Pu238ry[900], Cf249x[900], Cf249y[900], Hpu240, Hpu238, Hcf249, Herr, X[30], Xerr[30],Ediff[30], Lavg=0.0, tot=0.0, diff[900], prevdl, prevdlerr, mc[900], mcerr[900], mr[900], mrerr[900], bc[900], bcerr[900], br[900], brerr[900], Energyerr[30],PU238c[900], PU238r[900],Sl[30],Sle[30],Ict[30],Icte[30],icount[30],Sla=0.0,Sls=0.0,Icta=0.0,Icts=0.0,Sln=0.0,Ictn=0.0,Hpu238a,Herra,nmlength,nmlengtherr,nmtopx,nmtopy,nmtopye,nmtopxe,Hypoth[900],Hypothe[900],jcount[30];
        float Ccentd[30],Crtd[30],Ccented[30],Crted[30],invCc[30],invCr[30],invCce[30],invCre[30],values[150],valuese[150],mean=0.0,stdev=0.0;
        
        sprintf(center,"flatfinalsorted.txt");
        sprintf(right,"anglefinalsorted.txt");
        
        //ANGLE NOW CORRESPONDS TO RIGHT
        
        CENTER=fopen(center,"r");
        RIGHT=fopen(right,"r");
        COS=fopen("cosines.txt","w");
        fprintf(COS,"pixel (fstrip+32*(bstrip-1), Cos(Center), Cos(angled), 1/Cos(Center), 1/Cos(angled), Angle Difference, Angle Diff error, Energy difference, Energy Diff error \n");
        //j= diagnostic pixel
        k=104;
        
        
       // reading out numbers into arrays
        for(i=0; i<900; i++){
            r[i]=0;
            fscanf(CENTER,"%d", &A);
            fscanf(RIGHT,"%d", &A);
            
            //Pu240 energy
            fscanf(CENTER,"%f", &a);
            pu240ec[i]=a;
            fscanf(RIGHT,"%f", &a);
            pu240er[i]=a;
            
            //Pu240 energy error
            fscanf(CENTER,"%f", &a);
            pu240erc[i]=sqrt(a*a+1.231688*1.231688/4.0); //+1.231688 0.615844
            fscanf(RIGHT,"%f", &a);
            pu240err[i]=sqrt(a*a+1.231688*1.231688/4.0);
            
            //Pu240 sigma
            fscanf(CENTER,"%f", &a);
            pu240sc[i]=a;
            fscanf(RIGHT,"%f", &a);
            pu240sr[i]=a;
            
            //Pu240 sigma error
            fscanf(CENTER,"%f", &a);
            pu240serc[i]=a;
            fscanf(RIGHT,"%f", &a);
            pu240serr[i]=a;
            
            //Pu 240 chi^2
            fscanf(CENTER,"%f", &a);
            pu240cc[i]=a;
            fscanf(RIGHT,"%f", &a);
            pu240cr[i]=a;
            if(pu240cc[i]>3.0||pu240cr[i]>3.0) r[i]=r[i]+1;

            //Pu238 energy
            fscanf(CENTER,"%f", &a);
            pu238ec[i]=a;
            fscanf(RIGHT,"%f", &a);
            pu238er[i]=a;
            
            //Pu238 energy error
            fscanf(CENTER,"%f", &a);
            pu238erc[i]=sqrt(a*a+0.973391*0.973391/4.0); //0.973391 0.486696
            fscanf(RIGHT,"%f", &a);
            pu238err[i]=sqrt(a*a+0.973391*0.973391/4.0);
            
            //Pu238 sigma
            fscanf(CENTER,"%f", &a);
            pu238sc[i]=a;
            fscanf(RIGHT,"%f", &a);
            pu238sr[i]=a;
            
            //Pu238 sigma error
            fscanf(CENTER,"%f", &a);
            pu238serc[i]=a;
            fscanf(RIGHT,"%f", &a);
            pu238serr[i]=a;
            
            //Pu 238 chi^2
            fscanf(CENTER,"%f", &a);
            pu238cc[i]=a;
            fscanf(RIGHT,"%f", &a);
            pu238cr[i]=a;
            if(pu238cc[i]>3.0||pu238cr[i]>3.0) r[i]=r[i]+1;
          
            //Cf249 energy
            fscanf(CENTER,"%f", &a);
            cf249ec[i]=a;
            fscanf(RIGHT,"%f", &a);
            cf249er[i]=a;
            
            //Cf249 energy error
            fscanf(CENTER,"%f", &a);
            cf249erc[i]=sqrt(a*a+0.783275*0.783275/4.0); //0.391637
            fscanf(RIGHT,"%f", &a);
            cf249err[i]=sqrt(a*a+0.783275*0.783275/4.0);
            
            //Cf249 sigma
            fscanf(CENTER,"%f", &a);
            cf249sc[i]=a;
            fscanf(RIGHT,"%f", &a);
            cf249sr[i]=a;
            
            //Cf249 sigma error
            fscanf(CENTER,"%f", &a);
            cf249serc[i]=a;
            fscanf(RIGHT,"%f", &a);
            cf249serr[i]=a;
            
            //Cf249 chi^2
            fscanf(CENTER,"%f", &a);
            cf249cc[i]=a;
            fscanf(RIGHT,"%f", &a);
            cf249cr[i]=a;
            if(cf249cc[i]>3.0||cf249cr[i]>3.0) r[i]=r[i]+1;
            
            
        }
        
        
        //gStyle->SetOptStat(0);
        gStyle->SetLabelSize(0.047,"y");
        gStyle->SetLabelSize(0.047,"x");
        gStyle->SetTitleFontSize(0.06);
        /*
        TCanvas *c1 = new TCanvas("c1", "c1", 1200, 600);
        c1->Divide(4, 2);
        TCanvas *c2 = new TCanvas("c2", "c2", 1200, 600);
        c2->Divide(4, 2);
        TCanvas *c3 = new TCanvas("c3", "c3", 1200, 600);
        c3->Divide(4, 2);
        TCanvas *c4 = new TCanvas("c4", "c4", 1200, 600);
        c4->Divide(4, 2);

        */
        //diagnostics
        //printf("center: pu240=%f pu238=%f cf249=%f  \n", pu240ec[k],pu238ec[k],cf249ec[k]);
        //printf("center chi^2 : pu240=%f pu238=%f cf249=%f  \n", pu240cc[k],pu238cc[k],cf249cc[k]);
        //printf("right: pu240=%f pu238=%f cf249=%f  \n", pu240er[k],pu238er[k],cf249er[k]);
       // printf("right chi^2: pu240=%f pu238=%f cf249=%f  \n", pu240cr[k],pu238cr[k],cf249cr[k]);

        
        
        //INPUTS WITH ERRORS!!!
        //------------------------------------------------------------------------------
        
        //Energies of alphas
        Epu240=5168.17; //(15) from Matt's report, error negligible
        Epu238=5499.03; //(10)
        Ecf249=5813.3; //(10)
      
      float  Epu240e=0.15;
      float  Epu238e=0.10;
      float  Ecf249e=0.10;
       /*
        //Values from SRIM: stopping power of helium in silicon
        dEdx[0]=215.91-14.916*Epu240/1000.0; // units = kev/um
        dEdx[1]=215.91-14.916*Epu238/1000.0;
        dEdx[2]=215.91-14.916*Ecf249/1000.0;
        
        //errors (data found in SRIM spreadsheet)
        dEdxe[0]=sqrt(0.636*0.636*(Epu240/1000.0)*(Epu240/1000.0)+3.507*3.507);
        dEdxe[1]=sqrt(0.636*0.636*(Epu238/1000.0)*(Epu238/1000.0)+3.507*3.507);
        dEdxe[2]=sqrt(0.636*0.636*(Ecf249/1000.0)*(Ecf249/1000.0)+3.507*3.507);
        */
        //printf("stopping power from SRIM and errors: Pu240=%f +/- %f, Pu238=%f +/- %f, Cf249=%f +/- %f\n", dEdx[0], dEdxe[0], dEdx[1], dEdxe[1], dEdx[2], dEdxe[2]);
       
        
        //heights
        Hpu238=5.764-0.1260-0.267-3.003-0.050; //inches
        Hpu240=5.764-0.1260-0.267-3.0025-0.032;//inches
        Hcf249=5.764-0.1260-0.267-3.006-0.0522; //inches
        Hpu238a=5.764-0.1260-0.267-3.003-0.050;
        
        //convert to um
        Hpu238=(Hpu238*25.40-1.0)*1000.0;
        Hpu240=(Hpu240*25.40-1.0)*1000.0;
        Hcf249=(Hcf249*25.40-1.0)*1000.0;
        Hpu238a=(Hpu238a*25.40-1.0+0.85)*1000.0;

        //errors of heights are the same for all positions except for the dl one.
        Herr=1000.0; //um
        Herra=2000.0;
        
        //Using the previously measured dead layer to handle systematics when finding energy shifts (taken from Adrian's note on pixel dead layer measurement)
        prevdl=0.1243; //um
        prevdlerr=0.0083; //um
        
        
        //from bottom edge of the plate (looking down from the stool) the center of the first pixel is 7.475 p/m 0.001 mm
        //from the left edge of the plate (looking down from the stool) 151.696 pm 0.003 mm
        //initial values
        //see placement_caluclations (mathematica) for calculations
        //good means double check ended up a-ok!
        

        pu240x=(87.10-2.00)*1000.0; // for this one, we don't even know where the source dot is. :( approximating it as the center of the rod (good)
        pu240y=(9.31-2.00)*1000.0;//good
        pu238cx=(32.95-2.00)*1000.0;//good
        pu238cy=(31.45-2.00)*1000.0;//good
        pu238rx=(32.95-0.35-2.00)*1000.0;//good
        pu238ry=(31.45-2.00)*1000.0;//good
        cf249x=(85.35-2.00)*1000.0; // the dot has an 8.33 mm diameter. Fuck //good
        cf249y=(46.06-2.00)*1000.0;//good
        nmlength=sqrt(2.0)*Hpu238a;
        nmtopx=89.005*1000.0;
        nmtopy=23.44*1000.0;
        
        pu240xerr=3.0*1000.0;// don't know where the source dot is, but we're going to assume it is in the middle 0.25 inch diameter circle.
        pu240yerr=3.0*1000.0;
        pu238cxerr=1.0*1000.0;
        pu238cyerr=1.0*1000.0;//the dot is sufficently small that we can assume 1 mm error
        pu238rxerr=2.0*1000.0;
        pu238ryerr=1.0*1000.0;
        cf249xerr=4.0*1000.0; // taking the approximate radius of the cf dot as the error
        cf249yerr=4.0*1000.0;
        nmlengtherr=sqrt(2.0)*Herra;
        nmtopxe=1000.0;
        nmtopye=1500.0;
        
//Change the i one by one to go strip by strip
        
        for(i=3; i<4; i++){
            memset(Ediff,0.0,30);
            memset(Energyerr,0.0,30);
            memset(Xerr,0.0,30);
            memset(X,0.0,30);
           
            can=can+1;
            if ( can < 9)  c1->cd(can);
            if (can > 8 && can < 17)  c2->cd(can - 8);
            if (can > 16 && can < 25 ) c3->cd(can - 16);
            if (can > 24 && can < 33 ) c4->cd(can - 24);
            
            icount[i]=i+2.0;
            
            for(j=0; j<30; j++){
                
                
                float square(float b)
                {
                    float a;
                    a = b*b;
                    //printf("%f %f",a ,b);
                    return(a);
                }
                
                jcount[j]=j+2.0;
                
                //printf("i=%d\n",i+2);
                //placements of sources
                Pu240x[i+30*j]=pu240x-j*2000.0;
                Pu240y[i+30*j]=pu240y-i*2000.0;
                Pu238cx[i+30*j]=pu238cx-j*2000.0;
                Pu238cy[i+30*j]=pu238cy-i*2000.0;
                Pu238rx[i+30*j]=pu238rx-j*2000.0;
                Pu238ry[i+30*j]=pu238ry-i*2000.0;
                Cf249x[i+30*j]=cf249x-j*2000.0;
                Cf249y[i+30*j]=cf249y-i*2000.0;
                
                Hypoth[i+30*j]=sqrt(square(Hpu238a)+square(Pu238rx[i+30*j])+square(Pu238ry[i+30*j]));
        
                Hypothe[i+30*j]=sqrt((square(Hpu238a)/(square(Hpu238a)+square(Pu238rx[i+30*j])+square(Pu238ry[i+30*j])))*square(Herra)+(square(Pu238rx[i+30*j])/(square(Hpu238a)+square(Pu238rx[i+30*j])+square(Pu238ry[i+30*j])))*square(pu238rxerr)+(square(Pu238ry[i+30*j])/(square(Hpu238a)+square(Pu238rx[i+30*j])+square(Pu238ry[i+30*j])))*square(pu238ryerr));
                //printf("hypothenuse %f \n",Hypothe[i+30*j]);
                
                Cpu[i+30*j]=Hpu240/sqrt(Hpu240*Hpu240+(pu240y-i*2000.0)*(pu240y-i*2000.0)+(pu240x-j*2000.0)*(pu240x-j*2000.0));
                Ccf[i+30*j]=Hcf249/sqrt(Hcf249*Hcf249+(cf249y-i*2000.0)*(cf249y-i*2000.00)+(cf249x-j*2000.0)*(cf249x-j*2000.0));
                Ccent[i+30*j]=Hpu238/sqrt(Hpu238*Hpu238+(pu238cy-i*2000.0)*(pu238cy-i*2000.0)+(pu238cx-j*2000.0)*(pu238cx-j*2000.0));
                
                Crt[i+30*j]=((square(nmtopy-i*2000.0)+square(nmtopx-j*2000.0))-square(nmlength)-square(Hypoth[i+30*j]))/(-2.0*nmlength*sqrt(Hpu238a*Hpu238a+(pu238ry-i*2000.0)*(pu238ry-i*2000.0)+(pu238rx-j*2000.0)*(pu238rx-j*2000.0)));
                
                
                //printf("isthistheproblem: %f \n", Crt[i+30*j]);
                //Hpu238a/sqrt(Hpu238a*Hpu238a+(pu238ry-i*2000.0)*(pu238ry-i*2000.0)+(pu238rx+j*2000.0)*(pu238rx+j*2000.0));
               // printf("Cpu=%f, Ccf=%f, Ccent=%f, Crt=%f \n", Cpu[i+30*j], Ccf[i+30*j], Ccent[i+30*j], Crt[i+30*j]);
                
                //adding the fluctuation of the system to the line center errors
                
                //fprintf(COS,"%d %f %f %f %f %f %f %f \n", i+2+(j+1)*32,Cpu[i+30*j], Ccf[i+30*j],Ccent[i+30*j],Crt[i+30*j],(1/Ccent[i+30*j]),(1/Crt[i+30*j]),(-1/Ccent[i+30*j]+1/Crt[i+30*j]));
            
                        
                
                //errors for the cosine terms
                Ccenterr[i+30*j]=sqrt(square((sqrt(square(Hpu238)+square(Pu238cy[i+30*j])+square(Pu238cx[i+30*j]))-square(Hpu238)/sqrt(square(Hpu238)+square(Pu238cy[i+30*j])+square(Pu238cx[i+30*j])))/(square(Hpu238)+square(Pu238cy[i+30*j])+square(Pu238cx[i+30*j])))*square(Herr)+square((Hpu238*Pu238cy[i+30*j])/pow(square(Hpu238)+square(Pu238cy[i+30*j])+square(Pu238cx[i+30*j]),1.5))*square(pu238cyerr)+square((Hpu238*Pu238cx[i+30*j])/pow(square(Hpu238)+square(Pu238cy[i+30*j])+square(Pu238cx[i+30*j]),1.5))*square(pu238cxerr));
                
                Crterr[i+30*j]=sqrt((square(nmtopx-j*2000.0)/square(nmlength*Hypoth[i+30*j]))*square(nmtopxe)+(square(nmtopy-i*2000.0)/square(nmlength*Hypoth[i+30*j]))*square(nmtopye)+square((4.0*square(nmlength)*Hypoth[i+30*j]+2.0*Hypoth[i+30*j]*(square(nmtopx-j*2000.0)+square(nmtopy-i*2000.0)-square(nmlength)-square(Hypoth[i+30*j])))/(4.0*square(nmlength)*square(Hypoth[i+30*j])))*square(nmlengtherr)+square((4.0*nmlength*square(Hypoth[i+30*j])+2.0*nmlength*(square(nmtopx-j*2000.0)+square(nmtopy-i*2000.0)-square(nmlength)-square(Hypoth[i+30*j])))/(4.0*square(nmlength)*square(Hypoth[i+30*j])))*square(Hypothe[i+30*j]));
                                     
               // printf("hypoth contribution %f \n",square((4.0*nmlength*square(Hypoth[i+30*j]))/(4.0*square(nmlength)*square(Hypoth[i+30*j])))*square(Hypothe[i+30*j]));
                       
                Cpuerr[i+30*j]=sqrt(square((sqrt(square(Hpu240)+square(Pu240y[i+30*j])+square(Pu240x[i+30*j]))-square(Hpu240)/sqrt(square(Hpu240)+square(Pu240y[i+30*j])+square(Pu240x[i+30*j])))/(square(Hpu240)+square(Pu240y[i+30*j])+square(Pu240x[i+30*j])))*square(Herr)+square((Hpu240*Pu240y[i+30*j])/pow(square(Hpu240)+square(Pu240y[i+30*j])+square(Pu240x[i+30*j]),1.5))*square(pu240yerr)+square((Hpu240*Pu240x[i+30*j])/pow(square(Hpu240)+square(Pu240y[i+30*j])+square(Pu240x[i+30*j]),1.5))*square(pu240xerr));
                
                Ccferr[i+30*j]=sqrt(square((sqrt(square(Hcf249)+square(Cf249y[i+30*j])+square(Cf249x[i+30*j]))-square(Hcf249)/sqrt(square(Hcf249)+square(Cf249y[i+30*j])+square(Cf249x[i+30*j])))/(square(Hcf249)+square(Cf249y[i+30*j])+square(Cf249x[i+30*j])))*square(Herr)+square((Hcf249*Cf249y[i+30*j])/pow(square(Hcf249)+square(Cf249y[i+30*j])+square(Cf249x[i+30*j]),1.5))*square(cf249yerr)+square((Hcf249*Cf249x[i+30*j])/pow(square(Hcf249)+square(Cf249y[i+30*j])+square(Cf249x[i+30*j]),1.5))*square(cf249xerr));
                
               //printf("Errors: Cpu=%f, Ccf=%f, Ccent=%f, Crt=%f \n", Cpuerr[i+30*j], Ccferr[i+30*j], Ccenterr[i+30*j], Crterr[i+30*j]);
                
                //X axis values and errors for end graph that will determine the dead layer (slope)
                
                X[j]=(-1/Ccent[i+30*j]+1/Crt[i+30*j]);
                
                Xerr[j]=sqrt(square(1.0/square(Ccent[i+30*j]))*square(Ccenterr[i+30*j])+square(1.0/square(Crt[i+30*j]))*square(Crterr[i+30*j]));
                
                //printf("X[j]=%f +/- %f\n",X[j],Xerr[j]);

                
                //slopes and intercepts and errors
                
                mc[i+30*j]=((Ecf249/*-prevdl*dEdx[2]/Ccf[i+30*j]*/)-(Epu240/*-prevdl*dEdx[0]/Cpu[i+30*j]*/))/(cf249ec[i+30*j]-pu240ec[i+30*j]);
                
                mr[i+30*j]=((Ecf249/*-prevdl*dEdx[2]/Ccf[i+30*j]*/)-(Epu240/*-prevdl*dEdx[0]/Cpu[i+30*j]*/))/(cf249er[i+30*j]-pu240er[i+30*j]);

                //printf("Energies and prevdl: Ecf249=%f Epu240=%f, prevdl=%f\n",Ecf249,Epu240,prevdl);
                //printf("stopping powers: cf249=%f pu240=%f\n",dEdx[2],dEdx[0]);
                //printf("linecenters(c): cf249=%f p/m %f pu240=%f p/m %f\n",cf249ec[i+30*j],cf249erc[i+30*j], pu240ec[i+30*j],pu240erc[i+30*j]);
                //printf("linecenters(r): cf249=%f p/m %f pu240=%f p/m %f\n",cf249er[i+30*j],cf249err[i+30*j],pu240er[i+30*j],pu240err[i+30*j]);
                //printf("pu238 line centers: c=%f +/- %f  r=%f +/- %f \n",pu238ec[i+30*j],pu238erc[i+30*j],pu238er[i+30*j],pu238err[i+30*j]);
                
                //printf("slopes: mc=%f and mr=%f \n",mc[i+30*j],mr[i+30*j]);
                
                mcerr[i+30*j]=sqrt(/*square((dEdx[0]/Cpu[i+30*j]-dEdx[2]/Ccf[i+30*j])/(cf249ec[i+30*j]-pu240ec[i+30*j]))*square(prevdlerr)+square(prevdl/Ccf[i+30*j]/(cf249ec[i+30*j]-pu240ec[i+30*j]))*square(dEdxe[2])+square(prevdl/Cpu[i+30*j]/(cf249ec[i+30*j]-pu240ec[i+30*j]))*square(dEdxe[0])+*/square((Ecf249/*-prevdl*dEdx[2]/Ccf[i+30*j]*/-Epu240/*+prevdl*dEdx[0]/Cpu[i+30*j]*/)/square(cf249ec[i+30*j]-pu240ec[i+30*j]))*(square(cf249erc[i+30*j])+square(pu240erc[i+30*j]))/*+square(prevdl*dEdx[2]/(square(Ccf[i+30*j])*(cf249ec[i+30*j]-pu240ec[i+30*j])))*square(Ccferr[i+30*j])+square(prevdl*dEdx[0]/(square(Cpu[i+30*j])*(cf249ec[i+30*j]-pu240ec[i+30*j])))*square(Cpuerr[i+30*j]))*/+square(1.0/(cf249ec[i+30*j]-pu240ec[i+30*j]))*(square(Ecf249e)+square(Epu240e)));
                
                mrerr[i+30*j]=sqrt(/*square((dEdx[0]/Cpu[i+30*j]-dEdx[2]/Ccf[i+30*j])/(cf249er[i+30*j]-pu240er[i+30*j]))*square(prevdlerr)+square(prevdl/Ccf[i+30*j]/(cf249er[i+30*j]-pu240er[i+30*j]))*square(dEdxe[2])+square(prevdl/Cpu[i+30*j]/(cf249er[i+30*j]-pu240er[i+30*j]))*square(dEdxe[0])+*/square((Ecf249/*-prevdl*dEdx[2]/Ccf[i+30*j]*/-Epu240/*+prevdl*dEdx[0]/Cpu[i+30*j]*/)/square(cf249er[i+30*j]-pu240er[i+30*j]))*(square(cf249err[i+30*j])+square(pu240err[i+30*j]))/*+square(prevdl*dEdx[2]/(square(Ccf[i+30*j])*(cf249er[i+30*j]-pu240er[i+30*j])))*square(Ccferr[i+30*j])+square(prevdl*dEdx[0]/(square(Cpu[i+30*j])*(cf249er[i+30*j]-pu240er[i+30*j])))*square(Cpuerr[i+30*j]))*/+square(1.0/(cf249er[i+30*j]-pu240er[i+30*j]))*(square(Ecf249e)+square(Epu240e)));
                
                //printf("slope errors mcerr=%f, mrerr=%f\n", mcerr[i+30*j],mrerr[i+30*j]);
                
                if(pu240erc[i+30*j]<cf249erc[i+30*j]) {
                    bc[i+30*j]=(Epu240-prevdl*dEdx[0]/Cpu[i+30*j])-mc[i+30*j]*pu240ec[i+30*j];
                }
                else {
                  bc[i+30*j]=(Ecf249-prevdl*dEdx[2]/Ccf[i+30*j])-mc[i+30*j]*cf249ec[i+30*j];
                }
                if(pu240err[i+30*j]<cf249err[i+30*j]){
                    br[i+30*j]=(Epu240-prevdl*dEdx[0]/Cpu[i+30*j])-mr[i+30*j]*pu240er[i+30*j];
                }
                else{
                    br[i+30*j]=(Ecf249-prevdl*dEdx[2]/Ccf[i+30*j])-mr[i+30*j]*cf249er[i+30*j];
                }
                
                //printf("intercepts: bc=%f and br=%f \n",bc[i+30*j],br[i+30*j]);
                
                //bcerr[i+30*j]=sqrt(square(dEdx[0]/Cpu[i+30*j])*square(prevdlerr)+square(prevdl/Cpu[i+30*j])*square(dEdxe[0])+square((prevdl*dEdx[0])/square(Cpu[i+30*j]))*square(Cpuerr[i+30*j])+square(pu240ec[i+30*j])*square(mcerr[i+30*j])+square(mc[i+30*j])*square(pu240erc[i+30*j]));
                
                //brerr[i+30*j]=sqrt(square(dEdx[0]/Cpu[i+30*j])*square(prevdlerr)+square(prevdl/Cpu[i+30*j])*square(dEdxe[0])+square((prevdl*dEdx[0])/square(Cpu[i+30*j]))*square(Cpuerr[i+30*j])+square(pu240er[i+30*j])*square(mrerr[i+30*j])+square(mr[i+30*j])*square(pu240err[i+30*j]));
                
                //printf("intercept errors: bc=%f, br=%f \n",bcerr[i+30*j],brerr[i+30*j]);
                
                //Energyerr[j]=sqrt(square(pu238ec[i+30*j])*square(mcerr[i+30*j])+square(pu238er[i+30*j])*square(mrerr[i+30*j])+square(bcerr[i+30*j])+square(brerr[i+30*j])+square(mc[i+30*j])*square(pu238erc[i+30*j])+square(mr[i+30*j])*square(pu238err[i+30*j]));
                if(pu240erc[i+30*j]+pu240err[i+30*j]<cf249erc[i+30*j]+cf249err[i+30*j]){
                Energyerr[j]=sqrt(square(mc[i+30*j])*square(pu238erc[i+30*j])+square(mc[i+30*j])*square(pu240erc[i+30*j])+square(mcerr[i+30*j])*square(pu238ec[i+30*j]-pu240ec[i+30*j])+square(mr[i+30*j])*square(pu238err[i+30*j])+square(mr[i+30*j])*square(pu240err[i+30*j])+square(mrerr[i+30*j])*square(pu238er[i+30*j]-pu240er[i+30*j]));
                }
                else{
                Energyerr[j]=sqrt(square(mc[i+30*j])*square(pu238erc[i+30*j])+square(mc[i+30*j])*square(cf249erc[i+30*j])+square(mcerr[i+30*j])*square(pu238ec[i+30*j]-cf249ec[i+30*j])+square(mr[i+30*j])*square(pu238err[i+30*j])+square(mr[i+30*j])*square(cf249err[i+30*j])+square(mrerr[i+30*j])*square(pu238er[i+30*j]-cf249er[i+30*j]));
                }
                //if(r[i+30*j]>0){
                   // PU238c[i+30*j]=0.0;
                   // PU238r[i+30*j]=0.0;
                //}
                
                //printf("energy error %f \n", Energyerr[j]);
               // else {
                    PU238c[i+30*j]=mc[i+30*j]*pu238ec[i+30*j]+bc[i+30*j]; //both of these are in kev
                    PU238r[i+30*j]=mr[i+30*j]*pu238er[i+30*j]+br[i+30*j];
                Ediff[j]=PU238c[i+30*j]-PU238r[i+30*j];
                    z=1+z;
               // }
                //printf("Pu238: center:%f right:%f, #=%d\n",PU238c[i+30*j],PU238r[i+30*j],z);
               // printf("Energy  %f +/- %f \n",Ediff[j],Energyerr[j]);
                //float Ediff[z],Cdiff[z];
               // if(PU238c[i+30*j]!=0){
                 //   Ediff[y]=PU238c[i+30*j]-PU238r[i+30*j];
                  //  Cdiff[y]=dEdx[1]*(1.0/Ccent-1.0/Crt);
                   // y=y+1;
                
                 //}
                
                if(j==0){
                    fprintf(COS,"%d %f %f %f %f %f %f %f %f \n", i+2+(j+1)*32,Ccent[i+30*j],Crt[i+30*j],(1./Ccent[i+30*j]),(1./Crt[i+30*j]),(-1./Ccent[i+30*j]+1./Crt[i+30*j]),Xerr[j],Ediff[j],Energyerr[j]);
                }
                if(j<5 && j>=0){
                    values[q]=Ediff[j]/X[j];
                    valuese[q]=sqrt(square(1./X[j])*square(Energyerr[j])+square(Energyerr[j]/square(X[j]))*square(Xerr[j]));
                    valuese[q]=1./square(valuese[q]);
                    printf("Values we need %f +/- %f \n",values[q],valuese[q]);
                    mean=values[q]*valuese[q]+mean;
                    stdev=valuese[q]+stdev;
                    q=q+1;
                }

            }
            TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
            printf("fstrip=%d\n",i+2);
            sprintf(graphname,"Strip # %d",i+2);
            gr = new TGraphErrors(30,X,Ediff,Xerr,Energyerr);
            gr->Fit("pol1");
            TF1 *fx=gr->GetFunction("pol1");
            gr->SetTitle(graphname);
            gr->GetXaxis()->SetTitle("(#frac{1}{Cos(#theta_{2})} - #frac{1}{Cos(#theta_{1})}) ");
            gr->GetXaxis()->CenterTitle();
             gr->GetYaxis()->CenterTitle();
            gr->GetXaxis()->SetTitleOffset(1.5);
            gr->GetYaxis()->SetTitle("#Delta E (keV)");
            
            gr->Draw("AP");
            gr->SetMarkerStyle(kFullCircle);
            gr->SetMarkerSize(0.8);
            
            Sl[i]=fx->GetParameter(1);
            Sle[i]=fx->GetParError(1);
            Ict[i]=fx->GetParameter(0);
            Icte[i]=fx->GetParError(0);
            printf("slope %f +/- %f int %f +/- %f \n",Sl[i],Sle[i],Ict[i],Icte[i]);
/*
            if(i==17){
                for(k=0; k<30; k++){
                    Ccentd[k]=Ccent[i+30*k];
                    Crtd[k]=Crt[i+30*k];
                    Ccented[k]=Ccenterr[i+30*k];
                    Crted[k]=Crterr[i+30*k];
                    invCc[k]=1.0/Ccent[i+30*k];
                    invCr[k]=1.0/Crt[i+30*k];
                    invCce[k]=Ccenterr[i+30*k]/square(Ccent[i+30*k]);
                    invCre[k]=Crterr[i+30*k]/square(Crt[i+30*k]);
                    
                
                    
                }
                
            
                c15 = new TCanvas("c15","wooo",700,500);
            
            
                gr9 = new TGraphErrors(30,jcount,Ccentd,NULL,Ccented);
                gr9->SetTitle("Flat Cosine values");
                gr9->GetXaxis()->SetTitle("bstrip");
                gr9->GetXaxis()->CenterTitle();
                gr9->GetXaxis()->SetTitleOffset(1.2);
                gr9->GetYaxis()->SetTitle("Cos(#theta f)");
                gr9->GetYaxis()->CenterTitle();
                gr9->GetYaxis()->SetTitleOffset(1.4);
                gr9->Draw("AP");
                gr9->SetMarkerStyle(kFullCircle);
                gr9->SetMarkerSize(0.8);
                c16 = new TCanvas("c16","wooo",700,500);

                gr10 = new TGraphErrors(30,jcount,Crtd,NULL,Crted);
                gr10->SetTitle("Angled Cosine values");
                gr10->GetXaxis()->SetTitle("bstrip");
                gr10->GetXaxis()->CenterTitle();
                gr10->GetXaxis()->SetTitleOffset(1.2);
                gr10->GetYaxis()->SetTitle("Cos(#theta a)");
                gr10->GetYaxis()->CenterTitle();
                gr10->GetYaxis()->SetTitleOffset(1.0);
                gr10->Draw("AP");
                gr10->SetMarkerStyle(kFullCircle);
                gr10->SetMarkerSize(0.8);
            
                c17 = new TCanvas("c17","wooo",700,500);

                gr11 = new TGraphErrors(30,jcount,invCc,NULL,invCce);
                gr11->SetTitle("Flat inverse Cosine values");
                gr11->GetXaxis()->SetTitle("bstrip");
                gr11->GetXaxis()->CenterTitle();
                gr11->GetXaxis()->SetTitleOffset(1.2);
                gr11->GetYaxis()->SetTitle("1/Cos(#theta f)");
                gr11->GetYaxis()->CenterTitle();
                gr11->GetYaxis()->SetTitleOffset(1.5);
                gr11->Draw("AP");
                gr11->SetMarkerStyle(kFullCircle);
                gr11->SetMarkerSize(0.8);
                
                c18 = new TCanvas("c18","wooo",700,500);
                
                gr12 = new TGraphErrors(30,jcount,invCr,NULL,invCre);
                gr12->SetTitle("Angled inverse Cosine values");
                gr12->GetXaxis()->SetTitle("bstrip");
                gr12->GetXaxis()->CenterTitle();
                gr12->GetXaxis()->SetTitleOffset(1.2);
                gr12->GetYaxis()->SetTitle("1/Cos(#theta a)");
                gr12->GetYaxis()->CenterTitle();
                gr12->GetYaxis()->SetTitleOffset(1.5);
                gr12->Draw("AP");
                gr12->SetMarkerStyle(kFullCircle);
                gr12->SetMarkerSize(0.8);
                
                c19 = new TCanvas("c19","wooo",700,500);
                
                gr13 = new TGraphErrors(30,jcount,X,NULL,Xerr);
                gr13->SetTitle("Difference in inverse Cosines");
                gr13->GetXaxis()->SetTitle("bstrip");
                gr13->GetXaxis()->CenterTitle();
                gr13->GetXaxis()->SetTitleOffset(1.2);
                gr13->GetYaxis()->SetTitle("(#frac{1}{Cos(#theta_{2})} - #frac{1}{Cos(#theta_{1})}) ");
                gr13->GetYaxis()->CenterTitle();
                gr13->GetYaxis()->SetTitleOffset(1.2);
                gr13->Draw("AP");
                gr13->SetMarkerStyle(kFullCircle);
                gr13->SetMarkerSize(0.8);


                c20 = new TCanvas("c20","wooo",700,500);
                
                gr14 = new TGraphErrors(30,jcount,Ediff,NULL,Energyerr);
                gr14->SetTitle("Difference in energy");
                gr14->GetXaxis()->SetTitle("bstrip");
                gr14->GetXaxis()->CenterTitle();
                gr14->GetXaxis()->SetTitleOffset(1.2);
                gr14->GetYaxis()->SetTitle("Energy difference (kev)");
                gr14->GetYaxis()->CenterTitle();
                gr14->GetYaxis()->SetTitleOffset(1.0);
                gr14->Draw("AP");
                gr14->SetMarkerStyle(kFullCircle);
                gr14->SetMarkerSize(0.8);
            }
*/
        }
            
         c7 = new TCanvas("c7","hope this works",700,500);
        
        
        gr1 = new TGraphErrors(30,icount,Ict,NULL,Icte);
        gr1->SetTitle("Intercept by strip");
        gr1->GetXaxis()->SetTitle("fstrip");
        gr1->GetXaxis()->CenterTitle();
        gr1->GetXaxis()->SetTitleOffset(1.2);
        gr1->GetYaxis()->SetTitle("Intercept (kev)");
        gr1->GetYaxis()->CenterTitle();
        gr1->GetYaxis()->SetTitleOffset(1.0);
        gr1->Draw("AP");
        gr1->SetMarkerStyle(kFullCircle);
        gr1->SetMarkerSize(0.8);

         TH1F *h2 = new TH1F("h2", "Source Dead Layer Distribution", 20, -5.0, 12.0);
        for(i=0; i<150; i++){
            h2->Fill(values[i],valuese[i]);
        }
        c18 = new TCanvas("c18","histogram",700,500);
        h2->GetXaxis()->SetTitle("Dead Layer (kev/(#frac{1}{Cos(#theta_{2})} - #frac{1}{Cos(#theta_{1})}))");
        h2->GetXaxis()->CenterTitle();
        h2->GetXaxis()->SetTitleOffset(1.6);
        h2->Fit("gaus");
        h2->Draw();
        
        mean=mean/stdev;
        stdev=sqrt(1./stdev);
        printf("mean: %f+/-%f \n",mean,stdev);
       /* c8 = new TCanvas("c8","hope this ",700,500);
        
        
        gr1 = new TGraphErrors(30,icount,Sl,NULL,Sle);
        gr1->SetTitle("Dead layer by strip");
        gr1->GetXaxis()->SetTitle("fstrip");
        gr1->GetXaxis()->CenterTitle();
        gr1->GetXaxis()->SetTitleOffset(1.2);
        gr1->GetYaxis()->SetTitle("dead layer (#mu m))");
        gr1->GetYaxis()->CenterTitle();
        gr1->GetYaxis()->SetTitleOffset(1.5);
        gr1->Draw("AP");
        gr1->SetMarkerStyle(kFullCircle);
        gr1->SetMarkerSize(0.8);
*/
        c9 = new TCanvas("c9","hope this ",700,500);
        TH1F *h1 = new TH1F("h1", "Source Dead Layer offset", 26, 0.0, 6.5);
        for(i=0; i<30; i++){
            h1->Fill(Sl[i]);
        }
        h1->GetXaxis()->SetTitle("Energy offset/change in angle");
        h1->GetXaxis()->CenterTitle();
        h1->GetXaxis()->SetTitleOffset(1.2);
        h1->Draw();
        c9->Modified();
        c9->Update();
        fclose(COS);
        //printf("size of X %d size of E %d /n",sizeof(X)/sizeof(X[0]),sizeof(Ediff)/sizeof(Ediff[0]));
      
        
        for(i=0; i<30; i++){
            Sln=Sln+Sl[i]/square(Sle[i]);
            Ictn=Ictn+Ict[i]/square(Icte[i]);
            Sls=Sls+1.0/square(Sle[i]);
            Icts=Icts+1.0/square(Icte[i]);
        }
        Sla=Sln/Sls;
        Icta=Ictn/Icts;
        Sls=sqrt(1.0/Sls);
        Icts=sqrt(1.0/Icts);
        
        printf("slope = %f +/- %f intercept = %f +/- %f \n", Sla,Sls,Icta, Icts);
        
        for(i=0; i<30; i++){
            Sl[i]=(Sl[i]-Sla)/Sla;
            Sle[i]=Sle[i]/Sla;
        }
        
        c8 = new TCanvas("c8","hope this ",700,500);
        
        
        gr1 = new TGraphErrors(30,icount,Sl,NULL,Sle);
        gr1->SetTitle("Dead layer by strip");
        gr1->GetXaxis()->SetTitle("fstrip");
        gr1->GetXaxis()->CenterTitle();
        gr1->GetXaxis()->SetTitleOffset(1.2);
        gr1->GetYaxis()->SetTitle("fractional deviation");
        gr1->GetYaxis()->CenterTitle();
        gr1->GetYaxis()->SetTitleOffset(1.5);
        gr1->Draw("AP");
        gr1->SetMarkerStyle(kFullCircle);
        gr1->SetMarkerSize(0.8);

        
        //printf("Cpu=%f, Ccf=%f, Ccent=%f, Crt=%f \n", Cpu[k], Ccf[k], Ccent[k], Crt[k]);
        /*
        for(i=0; i<30; i++){
            if(r[i]>0) L[i]=0.0;
            else {
                tot=tot+1.0;
                L[i]=((Epu240-Ecf249)*(cf249ec[i]* (pu240er[i] - pu238er[i] ) + pu240ec[i]* (pu238er[i] - cf249er[i]) + pu238ec[i]* (-pu240er[i] + cf249er[i]))* Cpu[i])/(dEdx[0]* (cf249ec[i] *pu240er[i] + pu240ec[i]*pu238er[i] - cf249ec[i] *pu238er[i] - pu240ec[i]* cf249er[i] + pu238ec[i]* (-pu240er[i] + cf249er[i])) +Cpu[i]* (dEdx[2]* (pu238ec[i] *pu240er[i] - cf249ec[i]* pu240er[i] - pu240ec[i]* pu238er[i] + cf249ec[i]* pu238er[i] +pu240ec[i]* cf249er[i] - pu238ec[i]* cf249er[i])/Ccf[i] +dEdx[1]* (pu240ec[i] - cf249ec[i]) *(pu240er[i] - cf249er[i])*(-1/Ccent[i] + 1/Crt[i])));
                diff[i]=pu238ec[i]-pu238er[i];
                
                Lavg=Lavg+L[i];
            }
             //printf("L[%d]= %f, diff=%f\n", i,L[i],diff[i]);
        }
        
        Lavg=Lavg/tot;
        printf("average dead layer = %f \n", Lavg);
        */
                
        
            
        }
        
