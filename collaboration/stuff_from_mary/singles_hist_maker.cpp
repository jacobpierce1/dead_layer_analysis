//Si detectors calibration by alpha source, and Phillips 7164.

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TObject.h"
#include "TSystem.h"


// int get_strip_num(int x);

using namespace std;


const char *tree_file = "deadlayerdet3rt.root";


int  singles_hist_maker( void )
{
    TFile *fp = new TFile( tree_file );
    TTree *tree = (TTree *) fp->Get( "Singles_Tree" );
    //    tree->Process( "singles_selector.C" );

    // tree->Draw( "Singles_Energy" );  // plot all singles energy

    // see link for details: https://root.cern.ch/how/how-read-tree
    tree->Print();
    const int start_bin = 2500;
    const int stop_bin  = 3300;
    const int num_bins = stop_bin - start_bin;

    
    // TH1I *hist = new TH1I( "hist", "hist", num_bins, start_bin, stop_bin  );  // why root?
    //  hframe->GetXaxis()->SetRange( start_bin, stop_bin );
   
    tree->Draw( "Singles_Energy>>hnew(800,2500,3300)", "fstrip==20 && bstrip==20" ); //, "same" ); //  , "same" );
    // TH1I *hnew = (TH1I *) tree->GetHistogram();
    //hnew->GetXaxis()->SetRange( start_bin, stop_bin );
    //hnew->Draw();
    
// hist->Draw();
    
    TH1I *hnew = (TH1I *) gDirectory->Get( "hnew" );  // https://root.cern.ch/root/htmldoc/guides/users-guide/Trees.html#accessing-the-histogram-in-batch-mode
    hnew->Draw();
    
// hist->GetXaxis()->SetRange( start_bin, stop_bin );

//    hist->Draw();
    
    // gPad->SetLogy();
    


    /* 
    TH1I *hist = new TH1I( "test", "Counts vs. Channel: 0 degrees at Strip 1, 1", 5000, 0, 5000 );


    // declare vars to hold the data and set each one to hold data from a branch
    char global_scalers[128];
    
    
    tree->SetBranchAddress( "Global_Scalers", &global_scalers );

    // read in the array
    int num_entries = tree->GetEntries();
    for( int i=0; i<num_entries; i++)
    {
	tree->GetEntry(i);
	
    }

    

    hist->Draw();

    */

    
    fp->Close();
    return 0;
}



	


/*
  int main(int argc, char** argv)
  {
  char filename[100];
  int first_run=603, last_run=610;
  // argc=3;

  // string mypath( "/Users/edward/savard_group/2017_08_28_dead_layer_analysis/data_from_mary" );
  
  if(argc <= 1){
		
  sprintf(filename, "root_data/deadlayerdet3rt.root");
  //	sprintf(filename,"/Users/adrianperezgalvan/Desktop/physics/2014/ANL/8B/data/calibration/sources/5.27.2014/lincalinoRF.root");
    
  //	sprintf(filename,"/Users/mtburkey/Desktop/Beta-Correlation/selectorfiles/run%d.root", 24);

  } else{
  if(argc == 2){
  cout << "You are missing the second argument!" << endl;
  } else{
  if(argc == 3){
  first_run = atoi(argv[1]);
  last_run = atoi(argv[2]);
  sprintf(filename,"/Users/mtburkey/Desktop/Beta-Correlation/Deadlayers/Working_things/bypixel/deadlayerdet3rt.root");
		
  } else{
  return 0;
  }
  }
  }
	
  cout << "The root file is : " << filename << endl;
  TFile *f = new TFile(filename);
  TTree *Li_Tree = (TTree*)f->Get("Singles_Tree");

*/
/*char file_name[20];
  sprintf(file_name,"file_name%d",1);
  string outfile;
  outfile = (string)file_name;
  cout << outfile.c_str() << endl;*/

/*
  Li_Tree->Process( "singles_selector.C" );
  //    Li_Tree->Process( (mypath + "singles_selector.C").c_str() );
	
  return 0;
  }
*/


