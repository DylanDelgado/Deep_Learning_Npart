#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <pthread.h>
#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "src/TUrQMD.h"

//#define ASYNC_WRITER

using TURQMD::TUrQMD;

const int maxv = 20000;
const int vzero = 1;
struct RecordField{
  int Npart, mul;
  float b, atime;
  int pid[maxv];
  float px[maxv], py[maxv], pz[maxv];
} RecField;

void NewTree(TTree* &t, int c);
void* async_write(void*);
bool onwrite;
TTree *wtree;
TFile *EventRec;
std::string rname;

int main(void) {
  using namespace std;
  int i;
  const int event_num = (int) 6250;
  // how to split events into single root file
  const int step = 8000;
  // load random number
  ifstream fin("rdlist");
  int rd;
  fin>>rd;
  fin.close();

  // ready to emit UrQMD
  TUrQMD u;
  u.use_external_seed = false;
  u.quiet_out = true;
  u.skip_empty_event = false;
  u.rsd << Form("%d", rd);		// Random number seed
  u.nev << Form("%d", event_num);	//Event amount
  u.IMP << "0 16";    // Min and Max impact parameter
  u.ecm << "3.2";	    //Incident energy (GeV), Ecm here
  u.pro << "197 79";  // Project
  u.tar << "197 79";  // Target
  u.tim << "50 1";   // Total time and step in fm/c
  //u.CTO(13) = 1; //enhance f16 output
  //u.CTO(18) = 1; //1: unstable particles do not decay
  //u.eos << "1"; //skyrme potential, default is 0. casacde mode
  //u.f[13] = true;     // enable text output file13
  //u.f[14] = true;
  //u.f[15] = true;
  //u.f[16] = true;
  // emit UrQMD
  u.init();

  int start, end;
  static char RootName[255], TreeName[255];
  TTree *EventTree = NULL;
  EventRec = NULL;
  const int maxv=2000;
  TDirectory *cur = gDirectory;
  pthread_t wpid = 0;
  int pstat;
  bool isFinish;

  int mul, iev;

  NewTree(EventTree, 1);
  for(int iev = 1; iev <= u.nevents; ++iev) {
    //u.evt(iev, rdlist[iev-1]);
    u.evt(iev);

    RecField.mul = u.mul;
    RecField.Npart = u.Npart;
    RecField.b = u.bimp;
    // particle Loop
    u.write_pdata(RecField.pid,RecField.px,RecField.py,RecField.pz);
    EventTree->Fill();

    if (iev % step == 0 || iev == u.nevents) { // write a root file
      cout << "######### Event " << iev << " ###" << endl;
      isFinish = iev >= u.nevents;
#ifdef ASYNC_WRITER
      // waite writing thread
      if (wpid != 0) {
        if (onwrite) puts("Waiting for IO");
        pstat = 0;
        pstat = pthread_join(wpid, NULL);
        if (pstat) {
          cout << "Failed to finish IO thread on " << rname;
          cout << " Data has been dropped" << endl;
        }
      }
#endif
      // write to .root
      wtree = EventTree;
      if (!isFinish) {
        EventTree = NULL;
        NewTree(EventTree, iev+1);
      }
      rname = Form("Event_%d.root", iev);
      EventRec = new TFile(rname.c_str(), "RECREATE");
      cur->cd();

      // start writing thread
#ifdef ASYNC_WRITER
      if (!isFinish) {
        pstat = 0;
        pstat = pthread_create(&wpid, NULL, async_write, NULL);
        if (pstat) {
          cout << "Failed to start writing thread" << rname;
          cout << " Data has been dropped" << endl;
        }
      } else {
        async_write(NULL);
      }
#else
      async_write(NULL);
#endif
    }
  }// end Event Loop

  // finish simulation, output
  u.finish();
  return 0;
}

void NewTree(TTree* &t, int c) {
  // new tree
  if (t == NULL) {
    t = new TTree(Form("urqmd%d", c), "urqmd");
    // init Branch
    t->Branch("mul", &RecField.mul, "mul/I");
    t->Branch("b", &RecField.b, "b/F");
    t->Branch("Npart", &RecField.Npart, "Npart/I");
    t->Branch("pid", RecField.pid, "pid[mul]/I");
    t->Branch("px", RecField.px, "px[mul]/F");
    t->Branch("py", RecField.py, "py[mul]/F");
    t->Branch("pz", RecField.pz, "pz[mul]/F");
  }
}

void* async_write(void*) {
  using namespace std;
  onwrite = true;
  cout << "Async writing to " << rname << endl;
  EventRec->WriteTObject(wtree, "urqmd");
  EventRec->Close();
  delete EventRec;
  delete wtree;
  wtree = NULL;
  EventRec = NULL;
  cout << "Finish writing on " << rname << endl;
  onwrite = false;
  return nullptr;
}
