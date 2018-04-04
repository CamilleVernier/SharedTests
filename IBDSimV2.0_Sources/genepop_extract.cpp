/***************************************************************************
© R. Leblois 2001-2004
© F. Rousset 2005-2008

rousset@isem.univ-montp2.fr

This file is part of IBDSim. This software is a computer program
whose purpose is to perform population genetic simulations.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

 ***************************************************************************/
#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream> // cout...
#include <fstream> // for fstream
#include <sstream> // for stringstream
#include "genepop_extract.h"
#include "lattice_s.h"

using namespace std;



string EOLtype="";
const string fichierIn="fichier.in";

#ifdef WIN32
#include <windows.h>
void _gotoxy(int x,int y) {
    COORD mescoord;
    mescoord.X=x;
    mescoord.Y=y;
    SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE),mescoord);
}/**/
void effacer_ecran() {
    system("cls");
}
#else //vt100 terminals
void _gotoxy(int x,int y) {
//    printf("\033[%d;%dH",y+1,x+1);
      cout<<"\033["<<y+1<<";"<<x+1<<"H"<<flush; //equiv printf+flush ?
}
void effacer_ecran() {
//    system("clear"); //hum
//    printf("\033[2J");
      cout<<"\033[2J"<<flush; //equiv printf+flush ?
      _gotoxy(0,0);
}
#endif

using namespace std;

// utilitaire de comparaison de chaines insensible a la casse
//attention sortie non intuitive! voir des usages precedents
int cmp_nocase(const string& s, const string& s2) {
	string::const_iterator p = s.begin();
	string::const_iterator p2 = s2.begin();

	while(p != s.end() && p2 != s2.end()) {
		if (toupper(*p) != toupper(*p2)) return((toupper(*p)<toupper(*p2)) ? -1 : 1);
		++p; ++p2;
	}
	return((s2.size()==s.size()) ? 0 : (s2.size()<s.size()) ? -1 : 1);
}



// vire les blancs a droite
void rtrim(string *s) {
	while ((s->length()>0)  && (s->substr(s->length()-1,1)) == " ") {
			s->erase(s->length()-1,s->length());
	}
}

int set_eof_check_EOLtype(string EOLFileName, bool set_eof){
// CBR: previously known as CFichier_genepop::set_eof_check_EOLtype(bool set_eof /*=true*/)
//cout<<"debut set_eof_check_EOLtype()";getchar();
// no system call <=> not writing any new file <=> non change of length of file
	long pos,lngth;
	char c, buf;
EOLtype="";
fstream ORIfile(EOLFileName.c_str(),ios::in|ios::out);
    if (!ORIfile.is_open()) {
        remove(fichierIn.c_str()); // otherwise, gets stuck if fichier.in contains an incorrect file name
        cerr << "From set_eof_check_EOLtype(): Cannot open file!" << endl;
		cerr << "Check input file name " << EOLFileName.c_str() << endl;
		cin.ignore(); // vire le \n restant apr[`e]s une lecture dans cin
//	cerr<<"(Return) to show working directory and list of `*.' files:"<<endl;if (cinGetOnError) cin.get();
#ifdef WIN32
//	system("dir *. /p");
#else
//	system("ls *. | more");
#endif
        return(-1);
    }
    while ( ! ORIfile.eof() ) {
         c=ORIfile.get();
         if (c=='\r' || c=='\n') break;
    }
    if (ORIfile.eof()) { //there was no \r nor \n
        cerr<<"No line terminator in the file!";
        cerr<<"I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    } else if (c=='\n') {
        EOLtype="LinuxLF";
    } else /* we found \r first*/ if (ORIfile.get()=='\n'){
        EOLtype="WindowsCRLF"; //Windows \r\n
    } else { //\r not followed by \n
        cerr<<"(!) The file appears to contain a CR line terminator."<<endl;
        cerr<<"    This was the old MacIntosh end-of-line character."<<endl;
        cerr<<"    (i) If you are using Mac OS X, this should no longer occur (line terminator should be LF as in Unix)."<<endl;
        cerr<<"        However, some Microsoft editors for Mac OS X still use CR (sigh)."<<endl;
        cerr<<"        Genepop does not handle CR line terminators. Use a standard-compliant text editor."<<endl;
        cerr<<"    (ii) If you are using Windows/Linux, convert your input file to the correct format for your OS."<<endl;
        cerr<<"I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }

    if( ! set_eof )
    	return 0;

    /** ELSE **/

    ORIfile.clear();
    ORIfile.seekg(0,ios::end);
    lngth=ORIfile.tellg(); //longueur du fichier originel
// manitenant va lire charcater par char [`a] partir de la fin jusqu'[`a] trouver un chiffre
    pos=-1;
    ORIfile.seekg(pos,ios::end);
    ORIfile.read((char *) &buf, sizeof buf);
    while ( ! ( ((buf>47) && (buf<58)) || pos<= - lngth) ) { // stops if beginning of file reached [pos<= - lngth] or if  (numberfound && EOLtype.size()>0)
        pos--;
        ORIfile.seekg(pos,ios::end);
        ORIfile.read((char *) &buf, sizeof buf);
    }
    if (pos==-lngth) {
        cerr<<"No number, hence no genotype, in the file!";
        cerr<<"Exiting";
        if (cinGetOnError) cin.get();
        exit(-1);
    }
    // now we replace anything between the final number and the EOF by space chars
    ORIfile.seekp(lngth+pos+1);
    while (ORIfile.tellp()<lngth) ORIfile.put(' ');
    ORIfile.close();
return 1;
}

void set_UNIX_EOLtype(string EOLFileName) {
// CBR: previously part of CFichier_genepop::parseFile()
#ifdef WIN32
   // Unix-alike -> Windows
   // do nothing, because file 'as is' works anyway, while the conversion code below creates a blank line...
//    if (EOLtype=="LinuxLF") { // determined by ::set_eof_check_EOLtype !
//       cout << endl << "Input file seems to follow Unix-like end-of-line format. Here converted to Windows format." << endl;
//  	   if (cinGetOnError) getchar(); // no better condition at this point
//  	   {stringstream cmdline;
//        cmdline<<"TYPE "<<EOLFileName.c_str()<<" | FIND \"\" /V > tmp.tmp.tmp.txt"<<std::endl;
//        system(cmdline.str().c_str());
//	    remove(EOLFileName.c_str());
//        rename("tmp.tmp.tmp.txt",EOLFileName.c_str());
//       }
//    } // else nothing to do
#else
    if (EOLtype=="WindowsCRLF") { // determined by ::set_eof_check_EOLtype !
       // Windows -> Unix
        cout << endl << "Input file seems to follow Windows end-of-line format. Attempting conversion to Unix format." << endl;
        cerr << "(if this fails, find some way to convert the file before running Genepop)\n";
  	    if (cinGetOnError) cin.get();
        // plus simple mais moins portable ?
        {
  	        stringstream cmdline;
            cmdline<<"dos2unix -k -q "<<EOLFileName.c_str()<<"\n";
    		system(cmdline.str().c_str());
        }
		// don't use system() reurn call... may bear no relation to the cmdline
		// the return code may come from somthing else than the cmdline, and perror is equally useless
		 set_eof_check_EOLtype(EOLFileName, false); // ONLY recheck EOL
		 if (EOLtype=="WindowsCRLF") {
			 cerr<<"Problem calling 'dos2unix' on this OS. Trying 'fromdos' for conversion...\n";
			 // ubuntu's
			{
				 stringstream cmdline;
				 cmdline<<"fromdos "<<EOLFileName.c_str()<<"\n";
				 system(cmdline.str().c_str());
			}
			 set_eof_check_EOLtype(EOLFileName, false); // ONLY recheck EOL
			 if (EOLtype=="WindowsCRLF") {
				 int bidon=0;
				 cerr<<"Problem calling 'fromdos' on this OS";
				 cerr<<"Trying 'tr' for conversion...\n";
				 char masque[] = "fileXXXXXX";
				 bidon=mkstemp(masque);
				 if (bidon==-1) {
					 cerr<<"(!) From parseFile(): Problem calling mkstemp() in this directory. Exiting...";
					 if (cinGetOnError) cin.get();
					 exit(-1);
				 }
				 string tempname=masque;
				{
					 stringstream cmdline;
					 cmdline<<"tr -d '\\r' <"<<EOLFileName.c_str()<<" > "<<tempname.c_str()<<endl;
					 system(cmdline.str().c_str());
				}
				 set_eof_check_EOLtype(EOLFileName, false); // ONLY recheck EOL
				 if (EOLtype=="WindowsCRLF") {
					 cerr<<"(!) From parseFile(): Problem calling 'tr' on this OS. Exiting...";
					 if (cinGetOnError) cin.get();
					 exit(-1);
				 } else {
					 remove(EOLFileName.c_str());
					 rename(tempname.c_str(),EOLFileName.c_str());
				 }
			 }
			 else {
				 cout << "\nConversion from Windows EOL to Unix EOL successful! Press any key to continue...";
				 if (cinGetOnError)
					 cin.get();
			 }
		 }
		 else {
			 cout << "\nConversion from Windows EOL to Unix EOL successful! Press any key to continue...";
			 if (cinGetOnError)
				 cin.get();
		 }
    }
#endif
}



