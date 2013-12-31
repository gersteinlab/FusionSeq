extern "C" {
#include <bios/log.h>
#include <bios/format.h>
#include <bios/linestream.h>
}


#ifndef __CINT__

//--- ROOT includes ---
#include <root/TString.h>
#include <root/TROOT.h>
#include <root/TTree.h>
#include <root/TFile.h>
#include <root/TH1.h>
#include <root/TCanvas.h>

TROOT root("Rint","The ROOT Interactive Interface");

#endif



int main (int argc, char *argv[])
{ 
	LineStream ls;
	char *line;
	char *pos;
	Stringa buffer;

	if (argc != 2) {
		usage ("%s <file.intraOffsets>");
	}

	TH1 *his = new TH1D ("","Intra-read distribution",1000,0,1000);
	TCanvas *canv = new TCanvas("","canvas",1200,400);
	ls = ls_createFromFile (argv[1]);
	while (line = ls_nextLine (ls)) {
		his->Fill (atoi (line));
	}
	ls_destroy (ls);
	his->Draw();
	his->GetXaxis()->SetLabelSize (0.04);
	his->GetYaxis()->SetLabelSize (0.04);
	buffer = stringCreate (100);
	pos = strchr (argv[1],'.');
	if (pos == NULL) {
		die ("Expected <file.intraOffsets>: %s",argv[1]);
	}
	*pos = '\0';
	stringPrintf (buffer,"%s_intraDistribution.jpg",argv[1]);
	canv->Print (string (buffer),"jpg");
	stringDestroy (buffer);
	return 0;
}
