

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

void pixelanalysisnew()
    {
        //will create a graph using text files from runs documenting error and line center of each of 3 peaks. Shows error bar vs length of data taking for a strip and error bar over 3 different runs.
        gStyle->SetLabelSize(0.047,"y");
        gStyle->SetLabelSize(0.047,"x");
        gStyle->SetTitleFontSize(0.2);
        gStyle->SetTitleXSize(0.042);
        
        FILE *CENTER, *RIGHT, *COS;
        char center[50],right[50],canvasname[50],graphname[50];
        int i, r[900], pixel, A,j,k,z=0,y=0,can=0;
        float pu240ec[900], pu240er[900], pu240erc[900], pu240err[900], pu240sc[900], pu240sr[900], pu240serc[900], pu240serr[900], pu240cc[900], pu240cr[900];
        
        float pu238ec[900], pu238er[900], pu238erc[900], pu238err[900], pu238sc[900], pu238sr[900], pu238serc[900], pu238serr[900], pu238cc[900], pu238cr[900];
        
        float cf249ec[900], cf249er[900], cf249erc[900], cf249err[900], cf249sc[900], cf249sr[900], cf249serc[900], cf249serr[900], cf249cc[900], cf249cr[900];
        
        float a, L[900], dEdx[3], dEdxe[3], Epu240, Epu238, Ecf249, Cpu[900], Ccf[900], Ccent[900], Crt[900], Cpuerr[900], Ccferr[900], Ccenterr[900], Crterr[900], pu240x, pu240y, pu238cx, pu238cy, pu238rx, pu238ry, cf249x, cf249y, pu240xerr, pu240yerr, pu238cxerr, pu238cyerr, pu238rxerr, pu238ryerr, cf249xerr, cf249yerr,Pu240x[900], Pu240y[900], Pu238cx[900], Pu238cy[900], Pu238rx[900], Pu238ry[900], Cf249x[900], Cf249y[900], Hpu240, Hpu238, Hcf249, Herr, X[30], Xerr[30],Ediff[30], Lavg=0.0, tot=0.0, diff[900], prevdl, prevdlerr, mc[900], mcerr[900], mr[900], mrerr[900], bc[900], bcerr[900], br[900], brerr[900], Energyerr[30],PU238c[900], PU238r[900],Sl[30],Sle[30],Ict[30],Icte[30],icount[30],Sla=0.0,Sls=0.0,Icta=0.0,Icts=0.0,Sln=0.0,Ictn=0.0,srcdl,srcdle,Epu240e,Epu238e,Ecf249e,Xs[30],meanc,meanr,mcrms,mrrms,bcmean,brmean,bcrms,brrms,MC[30],MR[30],BC[30],BR[30],cmpu3c,cmpu4c,cmpu3r,cmpu4r,cmcfc,cmcfr;
        
        sprintf(center,"det1centcat.txt");
        sprintf(right,"det1rtcaterr.txt");
        
        
        
        CENTER=fopen(center,"r");
        RIGHT=fopen(right,"r");
        //COS=fopen("cosines.txt","w");
        //fprintf(COS,"pixel (fstrip+32*(bstrip-1), Cos(Pu240), Cos(Cf249), Cos(Center), Cos(Right), 1/Cos(Center), 1/Cos(Right), Difference \n");
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
            pu240erc[i]=sqrt(a*a);//+1.231688*1.231688/4.0); //+1.231688 0.615844
            fscanf(RIGHT,"%f", &a);
            pu240err[i]=sqrt(a*a);//+1.231688*1.231688/4.0);
            
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
            pu238erc[i]=sqrt(a*a);//+0.973391*0.973391/4.0); //0.973391 0.486696
            fscanf(RIGHT,"%f", &a);
            pu238err[i]=sqrt(a*a);//+0.973391*0.973391/4.0);
            
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
            cf249erc[i]=sqrt(a*a);//+0.783275*0.783275/4.0); //0.391637
            fscanf(RIGHT,"%f", &a);
            cf249err[i]=sqrt(a*a);//+0.783275*0.783275/4.0);
            
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
        
        
        gStyle->SetOptStat(0);
        gStyle->SetLabelSize(0.047,"y");
        gStyle->SetLabelSize(0.040,"x");
        gStyle->SetTitleFontSize(0.06);
        TCanvas *c1 = new TCanvas("c1", "c1", 1200, 600);
        c1->Divide(4, 2);
        TCanvas *c2 = new TCanvas("c2", "c2", 1200, 600);
        c2->Divide(4, 2);
        TCanvas *c3 = new TCanvas("c3", "c3", 1200, 600);
        c3->Divide(4, 2);
        TCanvas *c4 = new TCanvas("c4", "c4", 1200, 600);
        c4->Divide(4, 2);

        
        //diagnostics
        printf("center: pu240=%f pu238=%f cf249=%f  \n", pu240ec[k],pu238ec[k],cf249ec[k]);
        //printf("center chi^2 : pu240=%f pu238=%f cf249=%f  \n", pu240cc[k],pu238cc[k],cf249cc[k]);
        //printf("right: pu240=%f pu238=%f cf249=%f  \n", pu240er[k],pu238er[k],cf249er[k]);
        //printf("right chi^2: pu240=%f pu238=%f cf249=%f  \n", pu240cr[k],pu238cr[k],cf249cr[k]);

        
        
        //INPUTS WITH ERRORS!!!
        //------------------------------------------------------------------------------
        
        //Energies of alphas
        Epu240=5168.17; //(15) from Matt's report, error negligible
        Epu238=5499.03; //(10)
        Ecf249=5813.3; //(10)
      
        Epu240e=0.15;
        Epu238e=0.10;
        Ecf249e=0.10;
        
        //Values from SRIM: stopping power of helium in silicon
        dEdx[0]=215.91-14.916*Epu240/1000.0; // units = kev/um
        dEdx[1]=215.91-14.916*Epu238/1000.0;
        dEdx[2]=215.91-14.916*Ecf249/1000.0;
        
        //errors (data found in SRIM spreadsheet)
        dEdxe[0]=sqrt(0.636*0.636*(Epu240/1000.0)*(Epu240/1000.0)+3.507*3.507);
        dEdxe[1]=sqrt(0.636*0.636*(Epu238/1000.0)*(Epu238/1000.0)+3.507*3.507);
        dEdxe[2]=sqrt(0.636*0.636*(Ecf249/1000.0)*(Ecf249/1000.0)+3.507*3.507);
        
        //printf("stopping power from SRIM and errors: Pu240=%f +/- %f, Pu238=%f +/- %f, Cf249=%f +/- %f\n", dEdx[0], dEdxe[0], dEdx[1], dEdxe[1], dEdx[2], dEdxe[2]);
       
        
        //heights
        Hpu238=5.764-0.1260-0.267-3.003-0.050; //inches
        Hpu240=5.764-0.1260-0.267-3.0025-0.032;//inches
        Hcf249=5.764-0.1260-0.267-3.006-0.0522; //inches
        
        
        //convert to um
        Hpu238=(Hpu238*25.40-1.0)*1000.0;
        Hpu240=(Hpu240*25.40-1.0)*1000.0;
        Hcf249=(Hcf249*25.40-1.0)*1000.0;
        
        //errors of heights are the same for all positions.
        Herr=1000.0; //um

        
        //Using the previously measured dead layer to handle systematics when finding energy shifts (taken from Adrian's note on pixel dead layer measurement)
        prevdl=0.1243; //um
        prevdlerr=0.0083; //um
        
        srcdl=2.70;
        srcdle=0.14;
        
        //from bottom edge of the plate (looking down from the stool) the center of the first pixel is 7.475 p/m 0.001 mm
        //from the left edge of the plate (looking down from the stool) 151.696 pm 0.003 mm
        //initial values
        //see placement_caluclations (mathematica) for calculations
        //good means double check ended up a-ok!
        

        pu240x=(87.10-2.00)*1000.0; // for this one, we don't even know where the source dot is. :( approximating it as the center of the rod (good)
        pu240y=(9.31-2.00)*1000.0;//good
        pu238cx=(32.95-2.00)*1000.0;//good
        pu238cy=(31.45-2.00)*1000.0;//good
        pu238rx=(44.59-2.00)*1000.0;//good
        pu238ry=(28.58-2.00)*1000.0;//good
        cf249x=(85.35-2.00)*1000.0; // the dot has an 8.33 mm diameter. Fuck //good
        cf249y=(46.06-2.00)*1000.0;//good
        
        pu240xerr=3.0*1000.0;// don't know where the source dot is, but we're going to assume it is in the middle 0.25 inch diameter circle.
        pu240yerr=3.0*1000.0;
        pu238cxerr=1.0*1000.0;
        pu238cyerr=1.0*1000.0;//the dot is sufficently small that we can assume 1 mm error
        pu238rxerr=1.0*1000.0;
        pu238ryerr=1.0*1000.0;
        cf249xerr=4.0*1000.0; // taking the approximate radius of the cf dot as the error
        cf249yerr=4.0*1000.0;
        

//Change the i one by one to go strip by strip
        
        for(i=0; i<30; i++){
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
                //printf("i=%d\n",i+2);
                //placements of sources
                Pu240x[i+30*j]=pu240x-j*2000.0;
                Pu240y[i+30*j]=pu240y-i*2000.0;
                Pu238cx[i+30*j]=pu238cx-j*2000.0;
                Pu238cy[i+30*j]=pu238cy-i*2000.0;
                Pu238rx[i+30*j]=pu238rx+j*2000.0;
                Pu238ry[i+30*j]=pu238ry-i*2000.0;
                Cf249x[i+30*j]=cf249x-j*2000.0;
                Cf249y[i+30*j]=cf249y-i*2000.0;
                
                Cpu[i+30*j]=Hpu240/sqrt(Hpu240*Hpu240+(pu240y-i*2000.0)*(pu240y-i*2000.0)+(pu240x-j*2000.0)*(pu240x-j*2000.0));
                Ccf[i+30*j]=Hcf249/sqrt(Hcf249*Hcf249+(cf249y-i*2000.0)*(cf249y-i*2000.00)+(cf249x-j*2000.0)*(cf249x-j*2000.0));
                Ccent[i+30*j]=Hpu238/sqrt(Hpu238*Hpu238+(pu238cy-i*2000.0)*(pu238cy-i*2000.0)+(pu238cx-j*2000.0)*(pu238cx-j*2000.0));
                Crt[i+30*j]=Hpu238/sqrt(Hpu238*Hpu238+(pu238ry-i*2000.0)*(pu238ry-i*2000.0)+(pu238rx+j*2000.0)*(pu238rx+j*2000.0));
                //printf("Cpu=%f, Ccf=%f, Ccent=%f, Crt=%f \n", Cpu[i+30*j], Ccf[i+30*j], Ccent[i+30*j], Crt[i+30*j]);
                
                //adding the fluctuation of the system to the line center errors
                
                //fprintf(COS,"%d %f %f %f %f %f %f %f \n", i+2+(j+1)*32,Cpu[i+30*j], Ccf[i+30*j],Ccent[i+30*j],Crt[i+30*j],(1/Ccent[i+30*j]),(1/Crt[i+30*j]),(-1/Ccent[i+30*j]+1/Crt[i+30*j]));
            
                        
               
                
                float square(float b)
                {
                    float a;
                    a = b*b;
                    //printf("%f %f",a ,b);
                    return(a);
                }
                
                //errors for the cosine terms
                Ccenterr[i+30*j]=sqrt(square((sqrt(square(Hpu238)+square(Pu238cy[i+30*j])+square(Pu238cx[i+30*j]))-square(Hpu238)/sqrt(square(Hpu238)+square(Pu238cy[i+30*j])+square(Pu238cx[i+30*j])))/(square(Hpu238)+square(Pu238cy[i+30*j])+square(Pu238cx[i+30*j])))*square(Herr)+square((Hpu238*Pu238cy[i+30*j])/pow(square(Hpu238)+square(Pu238cy[i+30*j])+square(Pu238cx[i+30*j]),1.5))*square(pu238cyerr)+square((Hpu238*Pu238cx[i+30*j])/pow(square(Hpu238)+square(Pu238cy[i+30*j])+square(Pu238cx[i+30*j]),1.5))*square(pu238cxerr));
                
                Crterr[i+30*j]=sqrt(square((sqrt(square(Hpu238)+square(Pu238ry[i+30*j])+square(Pu238rx[i+30*j]))-square(Hpu238)/sqrt(square(Hpu238)+square(Pu238ry[i+30*j])+square(Pu238rx[i+30*j])))/(square(Hpu238)+square(Pu238ry[i+30*j])+square(Pu238rx[i+30*j])))*square(Herr)+square((Hpu238*Pu238ry[i+30*j])/pow(square(Hpu238)+square(Pu238ry[i+30*j])+square(Pu238rx[i+30*j]),1.5))*square(pu238ryerr)+square((Hpu238*Pu238rx[i+30*j])/pow(square(Hpu238)+square(Pu238ry[i+30*j])+square(Pu238rx[i+30*j]),1.5))*square(pu238rxerr));
                
                Cpuerr[i+30*j]=sqrt(square((sqrt(square(Hpu240)+square(Pu240y[i+30*j])+square(Pu240x[i+30*j]))-square(Hpu240)/sqrt(square(Hpu240)+square(Pu240y[i+30*j])+square(Pu240x[i+30*j])))/(square(Hpu240)+square(Pu240y[i+30*j])+square(Pu240x[i+30*j])))*square(Herr)+square((Hpu240*Pu240y[i+30*j])/pow(square(Hpu240)+square(Pu240y[i+30*j])+square(Pu240x[i+30*j]),1.5))*square(pu240yerr)+square((Hpu240*Pu240x[i+30*j])/pow(square(Hpu240)+square(Pu240y[i+30*j])+square(Pu240x[i+30*j]),1.5))*square(pu240xerr));
                
                Ccferr[i+30*j]=sqrt(square((sqrt(square(Hcf249)+square(Cf249y[i+30*j])+square(Cf249x[i+30*j]))-square(Hcf249)/sqrt(square(Hcf249)+square(Cf249y[i+30*j])+square(Cf249x[i+30*j])))/(square(Hcf249)+square(Cf249y[i+30*j])+square(Cf249x[i+30*j])))*square(Herr)+square((Hcf249*Cf249y[i+30*j])/pow(square(Hcf249)+square(Cf249y[i+30*j])+square(Cf249x[i+30*j]),1.5))*square(cf249yerr)+square((Hcf249*Cf249x[i+30*j])/pow(square(Hcf249)+square(Cf249y[i+30*j])+square(Cf249x[i+30*j]),1.5))*square(cf249xerr));
                
                //printf("Errors: Cpu=%f, Ccf=%f, Ccent=%f, Crt=%f \n", Cpuerr[i+30*j], Ccferr[i+30*j], Ccenterr[i+30*j], Crterr[i+30*j]);
                
                //X axis values and errors for end graph that will determine the dead layer (slope)
                
                X[j]=dEdx[1]*(-1/Ccent[i+30*j]+1/Crt[i+30*j]);
                
                Xerr[j]=sqrt(square(1/Ccent[i+30*j]-1/Crt[i+30*j])*square(dEdxe[1])+square(dEdx[1]/square(Ccent[i+30*j]))*square(Ccenterr[i+30*j])+square(dEdx[1]/square(Crt[i+30*j]))*square(Crterr[i+30*j]));
                
                Xs[j]=sqrt(square(1/Ccent[i+30*j]-1/Crt[i+30*j])*square(srcdle)+square(srcdl/square(Ccent[i+30*j]))*square(Ccenterr[i+30*j])+square(srcdl/square(Crt[i+30*j]))*square(Crterr[i+30*j]));
                
                //printf("X[j]=%f +/- %f\n",X[j],Xerr[j]);

                
                //slopes and intercepts and errors
                
                mc[i+30*j]=((Ecf249/*-prevdl*dEdx[2]/Ccf[i+30*j]*/)-(Epu240/*-prevdl*dEdx[0]/Cpu[i+30*j]*/))/(cf249ec[i+30*j]-pu240ec[i+30*j]);
                MC[j]=mc[i+30*j];
                
                mr[i+30*j]=((Ecf249/*-prevdl*dEdx[2]/Ccf[i+30*j]*/)-(Epu240/*-prevdl*dEdx[0]/Cpu[i+30*j]*/))/(cf249er[i+30*j]-pu240er[i+30*j]);
                MR[j]=mr[i+30*j];

                //if(i==6) printf("Energy difference C: %f R: %f \n",cf249ec[i+30*j]-pu240ec[i+30*j],cf249er[i+30*j]-pu240er[i+30*j]);
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
                    //if(i==6) printf("%d C: slope part of intercept %f \n",j,mc[i+30*j]*pu240ec[i+30*j]);
                }
                else {
                  bc[i+30*j]=(Ecf249-prevdl*dEdx[2]/Ccf[i+30*j])-mc[i+30*j]*cf249ec[i+30*j];
                    //if(i==6) printf("C: slope part of intercept %f \n",mc[i+30*j]*cf249ec[i+30*j]);
                }
                if(pu240err[i+30*j]<cf249err[i+30*j]){
                    br[i+30*j]=(Epu240-prevdl*dEdx[0]/Cpu[i+30*j])-mr[i+30*j]*pu240er[i+30*j];
                    //if(i==6) printf("%d R: slope part of intercept %f \n",j,mr[i+30*j]*pu240er[i+30*j]);
                }
                else{
                    br[i+30*j]=(Ecf249-prevdl*dEdx[2]/Ccf[i+30*j])-mr[i+30*j]*cf249er[i+30*j];
                   //if(i==6) printf("R: slope part of intercept %f \n",mr[i+30*j]*cf249er[i+30*j]);
                }
                BC[j]=bc[i+30*j];
                BR[j]=br[i+30*j];
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
               // else {
                    PU238c[i+30*j]=mc[i+30*j]*pu238ec[i+30*j]+bc[i+30*j]; //both of these are in kev
                    PU238r[i+30*j]=mr[i+30*j]*pu238er[i+30*j]+br[i+30*j];
                Ediff[j]=PU238c[i+30*j]-PU238r[i+30*j]-(srcdl*X[j]/dEdx[1]);
                    z=1+z;
                Energyerr[j]=sqrt(square(Energyerr[j])+square(Xs[j]));
               // }
                //printf("Pu238: center:%f right:%f, #=%d\n",PU238c[i+30*j],PU238r[i+30*j],z);
               //printf("Energy diff %f +/- %f \n",Ediff[j],Energyerr[j]);
                //float Ediff[j];
               // if(PU238c[i+30*j]!=0){
                 //   Ediff[y]=PU238c[i+30*j]-PU238r[i+30*j];
                  //  Cdiff[y]=dEdx[1]*(1.0/Ccent-1.0/Crt);
                   // y=y+1;
                
                 //}
                
             
                

            }
            
            meanc=TMath::Mean(30,&MC[0]);
            meanr=TMath::Mean(30,&MR[0]);
            
            mcrms=TMath::RMS(30,&MC[0]);
            mrrms=TMath::RMS(30,&MR[0]);
            
            bcmean=TMath::Mean(30,&BC[0]);
            brmean=TMath::Mean(30,&BR[0]);
            
            bcrms=TMath::RMS(30,&BC[0]);
            brrms=TMath::RMS(30,&BR[0]);
            
            //printf("slope means C: %f +/- %f R: %f +/- %f\n", meanc, mcrms,meanr,mrrms);
            
            //printf("intercept means C: %f +/- %f R: %f +/- %f\n", bcmean, bcrms,brmean,brrms);

            
            printf("fstrip=%d\n",i+2);
            sprintf(graphname,"Strip # %d",i+2);
            gr = new TGraphErrors(30,X,Ediff,Xerr,Energyerr);
            TF1 *fx = new TF1("fx","[0]+[1]*x", 0.0, 140.0);
            //fx->SetParameter(0,0.001);
            //fx->SetParLimits(0,0.0009,0.0011);
            gr->Fit(fx,"","",0.0,140.0);
            //gr->Fit("pol1");
            //TF1 *fx=gr->GetFunction("pol1");
            gr->SetTitle(graphname);
            gr->GetXaxis()->SetTitle("#frac{dE}{dx}(Pu238)(#frac{1}{Cos(#theta_{2})} - #frac{1}{Cos(#theta_{1})}) (keV/#mu m)");
            gr->GetXaxis()->CenterTitle();
            gr->GetXaxis()->SetTitleOffset(1.5);
            gr->GetYaxis()->SetTitle("Pu238 (Center pos energy-Right pos energy) (keV)");
            gr->GetYaxis()->SetTitleOffset(1.5);
            
            gr->Draw("AP");
            gr->SetMarkerStyle(kFullCircle);
            gr->SetMarkerSize(0.8);
            
            Sl[i]=fx->GetParameter(1);
            Sle[i]=fx->GetParError(1);
            Ict[i]=fx->GetParameter(0);
            Icte[i]=fx->GetParError(0);
            printf("slope %f +/- %f int %f +/- %f \n",Sl[i],Sle[i],Ict[i],Icte[i]);

        }
        
        cmpu3c=TMath::Mean(900,&pu238cc[0]);
        cmpu3r=TMath::Mean(900,&pu238cr[0]);
        cmpu4c=TMath::Mean(900,&pu240cc[0]);
        cmpu4r=TMath::Mean(900,&pu240cr[0]);
        cmcfc=TMath::Mean(900,&cf249cc[0]);
        cmcfr=TMath::Mean(900,&cf249cr[0]);
        printf("Average Chi^2 for each: pu240: %f %f pu238: %f %f cf249: %f %f \n",cmpu4c,cmpu4r,cmpu3c,cmpu3r,cmcfc,cmcfr);
        
        c9 = new TCanvas("c9","hope this ",700,500);
        TH1F *h1 = new TH1F("h1", "Det 1 Dead Layer distribution", 15, 40.0, 180.0);
        for(i=0; i<30; i++){
            h1->Fill(Sl[i]*1000.0);
        }
        h1->GetXaxis()->SetTitle("thickness (nm)");
        h1->GetXaxis()->CenterTitle();
        h1->GetXaxis()->SetTitleOffset(1.2);
        h1->SetLineWidth(2.0);
        h1->Draw();
        c9->Modified();
        c9->Update();
       
        for(i=0; i<30; i++){
            Sln=Sln+Sl[i]/square(Sle[i]);
            Ictn=Ictn+Ict[i]/square(Icte[i]);
            Sls=Sls+1.0/square(Sle[i]);
            Icts=Icts+1.0/square(Icte[i]);
            //printf("%f %f \n",Sl[i],Sle[i]);
        }
        Sla=Sln/Sls;
        Icta=Ictn/Icts;
        Sls=sqrt(1.0/Sls);
        Icts=sqrt(1.0/Icts);
        
        for(i=0; i<30; i++){
            Sl[i]=Sl[i]*1000.0;
           // Ict[i]=(Ict[i]-Icta)/Icta;
            Sle[i]=Sle[i]*1000.0;
           // Icte[i]=Icte[i]/Icta;
        }
      
         c7 = new TCanvas("c7","hope this works",700,500);
        
        
        gr1 = new TGraphErrors(30,icount,Ict,NULL,Icte);
        gr1->SetTitle("Det 3 Intercept by strip");
        gr1->GetXaxis()->SetTitle("fstrip");
        gr1->GetXaxis()->CenterTitle();
        gr1->GetXaxis()->SetTitleOffset(1.2);
        gr1->GetYaxis()->SetTitle("Intercept (kev)");
        gr1->GetYaxis()->CenterTitle();
        gr1->GetYaxis()->SetTitleOffset(1.0);
        gr1->Draw("AP");
        gr1->SetMarkerStyle(kFullCircle);
        gr1->SetMarkerSize(0.8);

        c8 = new TCanvas("c8","hope this ",700,500);
        
        
        gr1 = new TGraphErrors(30,icount,Sl,NULL,Sle);
        gr1->SetTitle("Det 1 Dead layer by strip");
        gr1->GetXaxis()->SetTitle("fstrip");
        gr1->GetXaxis()->SetLimits(0.0,32.0);
        gr1->GetXaxis()->CenterTitle();
        gr1->GetXaxis()->SetTitleOffset(1.2);
        gr1->GetYaxis()->SetTitle("Dead layer (nm)");
        gr1->GetYaxis()->CenterTitle();
        gr1->GetYaxis()->SetTitleOffset(1.5);
        double L1,L2,L3,L4,L5,L6;
        L1=119.6;
        L2=L1+0.7;
        L3=L1-0.7;
        L4=0.0;
        L5=32.0;
        TLine *l1 = new TLine(L4, L1,L5, L1);
        TLine *l2 = new TLine(L4, L2, L5, L2);
        TLine *l3 = new TLine(L4, L3, L5, L3);
        
        l1->SetLineWidth(2);
        l2->SetLineWidth(2);
        l3->SetLineWidth(2);
        l2->SetLineColor(2);
        l3->SetLineColor(2);
        gr1->Draw("AP");
        
        gr1->SetMarkerStyle(kFullCircle);
        gr1->SetMarkerSize(0.8);
        l1->Draw();
        l2->Draw();
        l3->Draw();
        
        c11 = new TCanvas("c11","hope this ",900,400);
        c11->Divide(2,1);
        c11->cd(1);
        TH1F *h5 = new TH1F("h5", "#chi^{2}'s of Center Histogram fits", 18, 0.0, 6.0);
        TH1F *h6 = new TH1F("h6", "#chi^{2}'s of Right Histogram fits", 18, 0.0, 6.0);
        TH1F *h7 = new TH1F("h7", "#chi^{2}'s of Histogram fits", 18, 0.0, 6.0);
        TH1F *h8 = new TH1F("h8", "#chi^{2}'s of Histogram fits", 18, 0.0, 6.0);
        TH1F *h9 = new TH1F("h9", "#chi^{2}'s of Center Histogram fits", 18, 0.0, 6.0);
        TH1F *h10 = new TH1F("h10", "#chi^{2}'s of Histogram fits", 18, 0.0, 6.0);
        for(i=0; i<900; i++){
            h5->Fill(pu240cc[i]);
            h6->Fill(pu240cr[i]);
            h7->Fill(pu238cc[i]);
            h8->Fill(pu238cr[i]);
            h9->Fill(cf249cc[i]);
            h10->Fill(cf249cr[i]);
            
        }
        h5->SetLineColor(kRed);
        h5->SetLineWidth(2.0);
        h6->SetLineColor(kRed);
        h6->SetLineWidth(2.0);
        h7->SetLineColor(kGreen);
        h7->SetLineWidth(2.0);
        h8->SetLineColor(kGreen);
        h8->SetLineWidth(2.0);
        h9->SetLineColor(kBlue);
        h9->SetLineWidth(2.0);
        h10->SetLineColor(kBlue);
        h10->SetLineWidth(2.0);
        h5->GetXaxis()->SetTitle("Reduced #chi^{2}");
        h5->GetXaxis()->CenterTitle();
        h5->GetXaxis()->SetTitleOffset(1.0);
        
        h6->GetXaxis()->SetTitle("Reduced #chi^{2}");
        h6->GetXaxis()->CenterTitle();
        h6->GetXaxis()->SetTitleOffset(1.0);
        
        h5->Draw();
        h7->Draw("sames");
        h9->Draw("sames");
        c11->cd(2);
        h6->Draw();
        h8->Draw("sames");
        h10->Draw("sames");
        
        
        
        /*
        h1->GetXaxis()->SetTitle("thickness (#mu m)");
        h1->GetXaxis()->CenterTitle();
        h1->GetXaxis()->SetTitleOffset(1.2);
        h1->Draw();
        */
        c11->Modified();
        c11->Update();


                //fclose(COS);
        //printf("size of X %d size of E %d /n",sizeof(X)/sizeof(X[0]),sizeof(Ediff)/sizeof(Ediff[0]));
      
        
        
        printf("slope = %f +/- %f intercept = %f +/- %f \n", Sla,Sls,Icta, Icts);
        
        
        
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
        /*sprintf(canvasname,"strip %d",strip);
        
        TCanvas *c5=new TCanvas("c5",canvasname,1000,900);
        c5->Divide(3, 3);
        c5->cd(1);
        TMultiGraph *mg1 = new TMultiGraph();
        TGraphErrors *gr1 = new TGraphErrors(3,runcount,Pu240r,NULL,Pu240er);
        gr1->SetMarkerColor(kBlue);
        gr1->SetMarkerStyle(kFullCircle);
        gr1->SetMarkerSize(.9);
        mg1->Add(gr1);
        TGraphErrors *gr2 = new TGraphErrors(3,runcount,Pu238r,NULL,Pu238er);
        gr2->SetMarkerColor(kRed);
        gr2->SetMarkerStyle(kFullCircle);
        gr2->SetMarkerSize(.8);
        mg1->Add(gr2);
        TGraphErrors *gr3 = new TGraphErrors(3,runcount,Cf249r,NULL,Cf249er);
        gr3->SetMarkerColor(kGreen);
        gr3->SetMarkerStyle(kFullCircle);
        gr3->SetMarkerSize(.7);
        
        mg1->Add(gr3);
        mg1->Draw("AP");
        mg1->SetTitle("line center vs run ");
        mg1->GetXaxis()->SetTitle("1st, 2nd, and 3rd 6hr run");
        mg1->GetYaxis()->SetTitle("offset from original value (channels)");
        mg1->GetXaxis()->CenterTitle();
        mg1->GetYaxis()->CenterTitle();
        mg1->GetXaxis()->SetTitleSize(0.05);
        mg1->GetYaxis()->SetTitleSize(0.05);
        mg1->GetYaxis()->SetTitleOffset(1.0);
        
        leg = new TLegend(0.1,0.7,0.48,0.9);
        leg->AddEntry(gr1,"lowest energy peak","ap");
        leg->AddEntry(gr2,"middle peak","ap");
        leg->AddEntry(gr3,"highest energy peak","ap");
        leg->Draw();


        c5->cd(2);
        TMultiGraph *mg2 = new TMultiGraph();
        TGraphErrors *gr12 = new TGraphErrors(3,runcount,Pu240sr,NULL,Pu240ser);
        gr12->SetMarkerColor(kBlue);
        gr12->SetMarkerStyle(kFullCircle);
        gr12->SetMarkerSize(.9);
        mg2->Add(gr12);
        TGraphErrors *gr22 = new TGraphErrors(3,runcount,Pu238sr,NULL,Pu238ser);
        gr22->SetMarkerColor(kRed);
        gr22->SetMarkerStyle(kFullCircle);
        gr22->SetMarkerSize(.8);
        mg2->Add(gr22);
        TGraphErrors *gr32 = new TGraphErrors(3,runcount,Cf249sr,NULL,Cf249ser);
        gr32->SetMarkerColor(kGreen);
        gr32->SetMarkerStyle(kFullCircle);
        gr32->SetMarkerSize(.7);
        mg2->Add(gr32);
        mg2->Draw("ap");
        mg2->SetTitle("sigma vs run ");
        mg2->GetXaxis()->SetTitle("1st, 2nd, and 3rd 6hr run");
        mg2->GetYaxis()->SetTitle("offset from original sigma value (channels)");
        mg2->GetXaxis()->CenterTitle();
        mg2->GetYaxis()->CenterTitle();
        mg2->GetXaxis()->SetTitleSize(0.05);
        mg2->GetYaxis()->SetTitleSize(0.05);
        mg2->GetYaxis()->SetTitleOffset(1.0);

        mg2->Draw("ap");
        
        c5->cd(3);
        TMultiGraph *mg5 = new TMultiGraph();
        TGraphErrors *gr15 = new TGraphErrors(3,runcount,Pu240cr,NULL,NULL);
        gr15->SetMarkerColor(kBlue);
        gr15->SetMarkerStyle(kFullCircle);
        gr15->SetMarkerSize(.9);
        mg5->Add(gr15);
        TGraphErrors *gr25 = new TGraphErrors(3,runcount,Pu238cr,NULL,NULL);
        gr25->SetMarkerColor(kRed);
        gr25->SetMarkerStyle(kFullCircle);
        gr25->SetMarkerSize(.8);
        mg5->Add(gr25);
        TGraphErrors *gr35 = new TGraphErrors(3,runcount,Cf249cr,NULL,NULL);
        gr35->SetMarkerColor(kGreen);
        gr35->SetMarkerStyle(kFullCircle);
        gr35->SetMarkerSize(.7);
        mg5->Add(gr35);
        mg5->Draw("ap");
        mg5->SetTitle("chi^2 vs run ");
        mg5->GetXaxis()->SetTitle("1st, 2nd, and 3rd 6hr run");
        mg5->GetYaxis()->SetTitle("reduced chi^2 of each fit");
        mg5->GetXaxis()->SetTitleSize(0.05);
        mg5->GetYaxis()->SetTitleSize(0.05);
        mg5->GetXaxis()->CenterTitle();
        mg5->GetYaxis()->CenterTitle();
        mg5->GetYaxis()->SetTitleOffset(1.0);

        
        mg5->Draw("ap");
        
        c5->cd(4);
        TMultiGraph *mg3 = new TMultiGraph();
        TGraphErrors *gr13 = new TGraphErrors(10,hrcount,Pu240t,NULL,Pu240et);
        gr13->SetMarkerColor(kBlue);
        gr13->SetMarkerStyle(kFullCircle);
        gr13->SetMarkerSize(.9);
        mg3->Add(gr13);
        TGraphErrors *gr23 = new TGraphErrors(10,hrcount,Pu238t,NULL,Pu238et);
        gr23->SetMarkerColor(kRed);
        gr23->SetMarkerStyle(kFullCircle);
        gr23->SetMarkerSize(.8);
        mg3->Add(gr23);
        TGraphErrors *gr33 = new TGraphErrors(10,hrcount,Cf249t,NULL,Cf249et);
        gr33->SetMarkerColor(kGreen);
        gr33->SetMarkerStyle(kFullCircle);
        gr33->SetMarkerSize(.7);
        mg3->Add(gr33);
        mg3->Draw("ap");
        mg3->SetTitle("line center vs time ");
        mg3->GetXaxis()->SetTitle("data over # hours");
        mg3->GetYaxis()->SetTitle("offset from original value (channels)");
        mg3->GetXaxis()->CenterTitle();
        mg3->GetYaxis()->CenterTitle();
        mg3->GetXaxis()->SetTitleSize(0.05);
        mg3->GetYaxis()->SetTitleSize(0.05);
        mg3->GetYaxis()->SetTitleOffset(1.0);

        
        mg3->Draw("ap");
        
        c5->cd(5);
        TMultiGraph *mg4 = new TMultiGraph();
        TGraphErrors *gr14 = new TGraphErrors(10,hrcount,Pu240st,NULL,Pu240set);
        gr14->SetMarkerColor(kBlue);
        gr14->SetMarkerStyle(kFullCircle);
        gr14->SetMarkerSize(.9);
        mg4->Add(gr14);
        TGraphErrors *gr24 = new TGraphErrors(10,hrcount,Pu238st,NULL,Pu238set);
        gr24->SetMarkerColor(kRed);
        gr24->SetMarkerStyle(kFullCircle);
        gr24->SetMarkerSize(.8);
        mg4->Add(gr24);
        TGraphErrors *gr34 = new TGraphErrors(10,hrcount,Cf249st,NULL,Pu238set);
        gr34->SetMarkerColor(kGreen);
        gr34->SetMarkerStyle(kFullCircle);
        gr34->SetMarkerSize(.7);
        mg4->Add(gr34);
        
        mg4->Draw("AP");
        mg4->SetTitle("sigma vs time ");
        mg4->GetXaxis()->SetTitle("data over # hours");
        mg4->GetYaxis()->SetTitle("offset from original value (channels)");
        mg4->GetXaxis()->SetTitleSize(0.05);
        mg4->GetYaxis()->SetTitleSize(0.05);
        mg4->GetXaxis()->CenterTitle();
        mg4->GetYaxis()->CenterTitle();
        mg4->GetYaxis()->SetTitleOffset(1.0);

        mg4->Draw("ap");
        
        c5->cd(6);
        TMultiGraph *mg6 = new TMultiGraph();
        TGraphErrors *gr16 = new TGraphErrors(10,hrcount,Pu240ct,NULL,NULL);
        gr16->SetMarkerColor(kBlue);
        gr16->SetMarkerStyle(kFullCircle);
        gr16->SetMarkerSize(.9);
        mg6->Add(gr16);
        TGraphErrors *gr26 = new TGraphErrors(10,hrcount,Pu238ct,NULL,NULL);
        gr26->SetMarkerColor(kRed);
        gr26->SetMarkerStyle(kFullCircle);
        gr26->SetMarkerSize(.8);
        mg6->Add(gr26);
        TGraphErrors *gr36 = new TGraphErrors(10,hrcount,Cf249ct,NULL,NULL);
        gr36->SetMarkerColor(kGreen);
        gr36->SetMarkerStyle(kFullCircle);
        gr36->SetMarkerSize(.7);
        mg6->Add(gr36);
        mg6->Draw("AP");
        mg6->SetTitle("chi^2 vs time ");
        mg6->GetXaxis()->SetTitle("data over # hours");
        mg6->GetYaxis()->SetTitle("reduced chi^2 of each fit");
        mg6->GetXaxis()->CenterTitle();
        mg6->GetYaxis()->CenterTitle();
        mg6->GetXaxis()->SetTitleSize(0.05);
        mg6->GetYaxis()->SetTitleSize(0.05);
        mg6->GetYaxis()->SetTitleOffset(1.0);


        mg6->Draw("ap");
        
        c5->cd(7);
        TMultiGraph *mg7 = new TMultiGraph();
        TGraphErrors *gr17 = new TGraphErrors(10,hrcount,Pu240tc,NULL,Pu240etc);
        gr17->SetMarkerColor(kBlue);
        gr17->SetMarkerStyle(kFullCircle);
        gr17->SetMarkerSize(.9);
        mg7->Add(gr17);
        TGraphErrors *gr27 = new TGraphErrors(10,hrcount,Pu238tc,NULL,Pu238etc);
        gr27->SetMarkerColor(kRed);
        gr27->SetMarkerStyle(kFullCircle);
        gr27->SetMarkerSize(.8);
        mg7->Add(gr27);
        TGraphErrors *gr37 = new TGraphErrors(10,hrcount,Cf249tc,NULL,Cf249etc);
        gr37->SetMarkerColor(kGreen);
        gr37->SetMarkerStyle(kFullCircle);
        gr37->SetMarkerSize(.7);
        mg7->Add(gr37);
        mg7->Draw("ap");
        mg7->SetTitle("line center vs time ");
        mg7->GetXaxis()->SetTitle("data over of cumulative hours");
        mg7->GetYaxis()->SetTitle("offset from original value (channels)");
        mg7->GetXaxis()->SetTitleSize(0.05);
        mg7->GetYaxis()->SetTitleSize(0.05);
        mg7->GetXaxis()->CenterTitle();
        mg7->GetYaxis()->CenterTitle();
        mg7->GetYaxis()->SetTitleOffset(1.0);
        
        
        mg7->Draw("ap");
        
        c5->cd(8);
        TMultiGraph *mg8 = new TMultiGraph();
        TGraphErrors *gr18 = new TGraphErrors(10,hrcount,Pu240stc,NULL,Pu240setc);
        gr18->SetMarkerColor(kBlue);
        gr18->SetMarkerStyle(kFullCircle);
        gr18->SetMarkerSize(.9);
        mg8->Add(gr18);
        TGraphErrors *gr28 = new TGraphErrors(10,hrcount,Pu238stc,NULL,Pu238setc);
        gr28->SetMarkerColor(kRed);
        gr28->SetMarkerStyle(kFullCircle);
        gr28->SetMarkerSize(.8);
        mg8->Add(gr28);
        TGraphErrors *gr38 = new TGraphErrors(10,hrcount,Cf249stc,NULL,Pu238setc);
        gr38->SetMarkerColor(kGreen);
        gr38->SetMarkerStyle(kFullCircle);
        gr38->SetMarkerSize(.7);
        mg8->Add(gr38);
        
        mg8->Draw("AP");
        mg8->SetTitle("sigma vs time ");
        mg8->GetXaxis()->SetTitle("data over # of cumulative hours");
        mg8->GetYaxis()->SetTitle("offset from original value (channels)");
        mg8->GetXaxis()->SetTitleSize(0.05);
        mg8->GetYaxis()->SetTitleSize(0.05);
        mg8->GetXaxis()->CenterTitle();
        mg8->GetYaxis()->CenterTitle();
        mg8->GetYaxis()->SetTitleOffset(1.0);
        
        mg8->Draw("ap");
        
        c5->cd(9);
        TMultiGraph *mg9 = new TMultiGraph();
        TGraphErrors *gr19 = new TGraphErrors(10,hrcount,Pu240ctc,NULL,NULL);
        gr19->SetMarkerColor(kBlue);
        gr19->SetMarkerStyle(kFullCircle);
        gr19->SetMarkerSize(.9);
        mg9->Add(gr19);
        TGraphErrors *gr29 = new TGraphErrors(10,hrcount,Pu238ctc,NULL,NULL);
        gr29->SetMarkerColor(kRed);
        gr29->SetMarkerStyle(kFullCircle);
        gr29->SetMarkerSize(.8);
        mg9->Add(gr29);
        TGraphErrors *gr39 = new TGraphErrors(10,hrcount,Cf249ctc,NULL,NULL);
        gr39->SetMarkerColor(kGreen);
        gr39->SetMarkerStyle(kFullCircle);
        gr39->SetMarkerSize(.7);
        mg9->Add(gr39);
        mg9->Draw("AP");
        mg9->SetTitle("chi^2 vs time ");
        mg9->GetXaxis()->SetTitle("data over # of cumulative hours");
        mg9->GetYaxis()->SetTitle("reduced chi^2 of each fit");
        mg9->GetXaxis()->SetTitleSize(0.05);
        mg9->GetYaxis()->SetTitleSize(0.05);
        mg9->GetXaxis()->CenterTitle();
        mg9->GetYaxis()->CenterTitle();
        mg9->GetYaxis()->SetTitleOffset(1.0);
        
    */
       // mg9->Draw("ap");
    
         
        //TGraphErrors *g2 = new TGraphErrors(32, Icount, sig, NULL, NULL);
    }
