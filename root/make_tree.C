//
// program to convert a simc ascii n-tuplw to a root tree
//
// read data ntuple from stdin

#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"

using namespace std;

// create a tree from the simc ascii file

using namespace std;

int main()
{

  TString NTtitle; // NT title
  TString NTname ; // NT name
  
  Int_t Ntags;  // number of tags

  Int_t i;
  

  // Tree variables and branch names
  Float_t val[100];
  Float_t Normfac;
  
  TString tag_name[100]; // tag names
  TString var_branch[100]; // branch name string
  TString var_name[100]; 
  TString var_type("/F"); // all variables of the same type

  // read file
  // header information
  // this factor MUST be added by hand or by script
  cin >> Normfac;
  // this is the default output
  cin >> NTtitle;
  cin >> NTname;
  cin >> Ntags;

  //  cout << NTtitle << endl;
  //  cout << NTname << endl;
  //  cout << Ntags << endl;

  // read tag names
  for (i = 0; i<Ntags; i++) {
    cin >> tag_name[i];
    //    cout << tag_name[i] << endl;
  }

  // create output file
  TFile *f = new TFile("simc.root", "RECREATE", NTtitle);

  // Create the tree

  TTree *stree = new TTree(NTname,NTtitle);
  
  // the tree variable 
  // 1st brang contains the norm factor
  
  var_branch[0] = "Normfac"; // get rid of blanks
  var_name[0] = "Normfac"; // get rid of blanks
  var_name[0].Append( var_type );
  stree->Branch(var_branch[0].Data(), &Normfac, var_name[0]);
  cout << "adding branch : " << var_branch[0].Data() << endl;

  Int_t offset = 1;
  for (i=0; i<Ntags;i++){
    var_branch[i+offset] = tag_name[i].Strip(); // get rid of blanks
    var_name[i+offset] = tag_name[i].Strip(); // get rid of blanks
    var_name[i+offset].Append( var_type );

    // add the branches to the tree

    stree->Branch(var_branch[i+offset].Data(), &val[i], var_name[i+offset]);

    cout << "adding branch : " << var_branch[i+offset].Data() << endl;
  }

  
  // read the values and fill the tree

  Int_t n = 0;
  while(! cin.eof() ){         // loop until eof
    n++;

    if (n%1000 == 0) cout << "Event : " << n << endl;

    for (i=0; i<Ntags; i++) {  // read the data for each event
      cin >> val[i];          
                              
    }
    // fill the tree's branches
    stree->Fill();
  }
  cout << "Read " << n << " events" << endl;

  f->Write(); // write output file
  return(0);

}
