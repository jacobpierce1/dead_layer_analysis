//Si detectors calibration by alpha source, and Phillips 7164.

#include <ctime>
#include <sys/stat.h>
#include <cstdio>
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
#include <cstdint>



#define LEN( array ) ( sizeof(array) / sizeof((array)[0]) )
#define PRINT_INT( var ) printf( "%s: %d\n", #var, (var) ) 
#define PRINT_STR( var ) printf( "%s: %s\n", #var, (var) )
#define PRINT_DOUBLE( var ) printf( "%s: %f\n", #var, (var) )
#define PRINT_LONG( var ) printf( "%s: %ld\n", #var, (var) )


// #define DEBUG  1
// #define USE_TEST_TREE  1
#define WRITE_FILES  1
// #define WRITE_ONE_FILE  1


const int num_channels = 32;  // we have a 32x32 array of pixels on the detector.

using namespace std;

const char *test_tree = "deadlayerdet3rt.root";

int open_all_files( const char *rootfile, const char *dirname, ofstream outfiles[32][32][2] );
int write_all_files( const char *rootfile, const char *dirname );
int close_all_files( ofstream outfiles[32][32][2]);
int clear_all_files( const char *rootfile, const char *dirname );
void debug( TTree *tree, Int_t *singles_scalers, double *singles_energies );
int get_prefix( char *prefix, const char *rootfile );


// int get_strip_num(int x);



int ttree_to_bin( const char *rootfile, const char *dirname )
{
    int ret = 0;
#ifdef USE_TEST_TREE
    ret = write_all_files( test_tree, dirname );
#else
    ret = write_all_files( rootfile, dirname );
#endif
    return ret ? 0 : -1; /// convert to unix standard.
}







int write_all_files( const char *rootfile, const char *dirname )
{
    // get time 
    time_t start_time;
    time( &start_time );
    time_t current_time = start_time;
    
    ofstream outfile;

    // this shall hold the name of the file
    char tmp_file_name[128];

    // remove the .root suffix
    char prefix[64];
    get_prefix( prefix, rootfile );

    // open the tree for reading    
    TFile *fp = new TFile( rootfile );
    if( !fp )
    {
	printf( "ERROR: unable to open file: %s\n", rootfile );
	return 0;
    }
    TTree *tree = (TTree *) fp->Get( "Singles_Tree" );
    if( ! tree )
    {
	puts( "ERROR: unable to obtain the requested tree from this file." );
	return 0;
    }
    

    Int_t global_scalers[9];
    Int_t singles_scalers[5];
    double singles_energy[2];
    
    tree->SetBranchAddress( "Global_Scalers", &(global_scalers[0]) );
    tree->SetBranchAddress( "Singles_Scalers", &(singles_scalers[0]) );
    tree->SetBranchAddress( "Singles_Energy", &(singles_energy[0]) );
    
#ifdef DEBUG
    debug( tree, singles_scalers, singles_energy );
#endif

    
#ifdef WRITE_FILES
    // write a header to a dedicated file
    const char *header = "efront eback t_eject t_cap cap_state n_cap rf_phase clock clock_total run_n event_n mult time fstrip bstrip det \n";
    sprintf( tmp_file_name, "%s%s", dirname, "header.txt" );
    outfile.open( tmp_file_name );
    outfile << header << endl;

    
    // vars we need 
    int num_hundredths_finished = 0;
    int num_entries = tree->GetEntries();    

    //  static bool file_opened[32][32]; //static init to 0


    clear_all_files( rootfile, dirname );
    
    // initiate for loop through tree 
    for( int j=0; j<num_entries; j++)
    {
	// populate the branch addresses 
	tree->GetEntry(j);

	// construct file  name and open stream 
	int fstrip = singles_scalers[2];
	int bstrip = singles_scalers[3];

	// subtract 1 from fstrip and bstrip in order to correct original start-at-1 indexing
	sprintf( tmp_file_name, "%s%s_%d_%d.bin", dirname, prefix, fstrip-1, bstrip-1 );
#ifdef WRITE_ONE_FILE
	PRINT_STR( tmp_file_name );
#endif
	    
	// open flie. check if has bene opened first. 
/*
  if( file_opened[fstrip][bstrip] )
  {
  puts("second");
  outfile.open( tmp_file_name, ios::out | ios::binary | ios::app );

  } else
  {
  puts("first");
  outfile.open( tmp_file_name, ios::out | ios::binary | ios::trunc );
  file_opened[fstrip][bstrip] = 1;
  }

  PRINT_INT( outfile.is_open() );
*/

	// keep trying to open file until it is good.
	while (1)
	{
	    outfile.open( tmp_file_name, ios::out | ios::binary | ios::app );
	    if ( outfile.good() )  break;
	    printf( "WARNING: failed to open %s, trying again...\n", tmp_file_name ); 
	    outfile.close();
	}

	
	// state progress
	if( j*1.0 / num_entries >= num_hundredths_finished / 100.0 )
	{
	    time( &current_time );
	    printf( "%d/100 complete, ETA: %f mins\n", num_hundredths_finished, ( difftime(current_time, start_time) / 60.0 ) * (100 - num_hundredths_finished ) / num_hundredths_finished );
	    ++num_hundredths_finished;
	}


	// write the data
#ifdef DEBUG_FILE_OPEN
	if( ! outfile.good() )
	{
	    puts("ERROR: outfile bad before write" );
	}
	PRINT_LONG( (long) outfile.tellp() );
#endif
	outfile.write( (char *) singles_energy, 2*sizeof(double) );
#ifdef DEBUG_FILE_OPEN
	PRINT_LONG( (long) outfile.tellp() );
#endif
	if( ! outfile.good() )
	{
	    puts("ERROR: write error.");
	}
	
// outfile.write( (char *) singles_scalers, 2*sizeof(Int_t) );
	// outfile.write( (char *) singles_scalers + 4, 
	
	// close file 
	outfile.close();
#ifdef WRITE_ONE_FILE
	if( j>=10) break;
#endif
    }
#endif  // WRITE_FILES

    fp->Close(); // close the root file
    return 1;
}



// eliminate very strange segfault error when resetting files.
int clear_all_files( const char *rootfile, const char *dirname )
{
    puts( "INFO: clearing all files." );
    char tmp_file_name[128];

    // remove the .root suffix
    char prefix[64];
    get_prefix( prefix, rootfile );
    
    for( int i=0; i<num_channels; i++)
    {
	for( int j=0; j<num_channels; j++)
	{
	    sprintf( tmp_file_name, "%s%s_%d_%d.bin", dirname, prefix, i, j );
	    ofstream current_file;
	    current_file.open( tmp_file_name, ios::out | ios::binary | ios::trunc );
	    
	    if( ! current_file.good() )
	    {
		printf( "ERROR: unable to trunacte file: %s\n", tmp_file_name );
	    }
	    current_file.close();
	    if( current_file.is_open() )
	    {
		printf( "ERROR: unable to close file: %s\n", tmp_file_name );
	    }
	}
    }
    puts( "INFO: successfully cleared all files." );
    return 1;
}



void debug( TTree *tree, Int_t *singles_scalers, double *singles_energies )
{
    tree->Print();
    tree->Show(0);
    tree->GetEntry(0);

    int fstrip = singles_scalers[2];
    int bstrip = singles_scalers[3];

    double efront = singles_energies[0];
    double eback = singles_energies[1];

    puts( "\nData we read: ") ;
    PRINT_INT( (int) sizeof(double) );
    PRINT_INT( fstrip );
    PRINT_INT( bstrip );
    PRINT_DOUBLE( efront );
    PRINT_DOUBLE( eback );
    puts("\n\n");
}



// populate char * with prefix of the .root file or a default if in test mode 
int get_prefix( char *prefix, const char *rootfile )
{
#ifndef USE_TEST_TREE
    const char *last_slash = strrchr( rootfile, '/' );
    if( ! last_slash ) last_slash = rootfile;
    strcpy( prefix, last_slash + 1 );
    * strstr( prefix, "." ) = '\0';
#else
    strcpy( prefix, "test" );
#endif
    return 1;
}
