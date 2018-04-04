#include <cstdio>
#include <iostream> // cout...
#include <fstream> // ofstream
#include <sstream> //stringstream
#include <cstdlib> //needed for vector too!!
#include <cstring>
#include <cmath>
#include <algorithm> // required for min_element with ubuntu linux g++...
#include <vector>
//#ifdef WIN32
//#include <windows.h> // for Sleep();
//#endif
//#include <unistd.h> // required for sleep, not used anymore, replaced by mysleep()
//#include <ctime>
//#include <limits>

#include "isoutils2.h"
#include "lattice_s.h"
#include "settings.h"
#include "TimeVarPars.h"


using namespace std;

string boundaryType="Absorbing";

bool warningBool = false;
vector<bool> FoundDispSpatiallyHeterogInSettings(1,false);/*T/F si dispersion homogene*/


void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}


int cmp_nocase_no__(const string& cs, const string& cs2) {

    string s(cs); // makes local copy that one could modify. Aletrnatively, (string& s, string& s2) [no const] would modify its arguments
    string s2(cs2);
    replaceAll(s,"_","");
    replaceAll(s2,"_","");

	return(cmp_nocase(s,s2));
}


int evaluateBool(bool &boolean, string buf) { // safe assignment of value `buf' to 'boolean'
	stringstream strstr(buf);
    string locstring;
    strstr>>locstring;
    if(cmp_nocase(locstring,"")==0 || cmp_nocase(locstring,"T")==0 ||
          cmp_nocase(locstring,"True")==0 || cmp_nocase(locstring,"Yes")==0 ||
          cmp_nocase(locstring,"Y")==0)
       boolean=true;
    else if(cmp_nocase(locstring,"F")==0 || cmp_nocase(locstring,"False")==0 ||
          cmp_nocase(locstring,"No")==0 || cmp_nocase(locstring,"N")==0)
       boolean=false;
    else {
        cout<<"(!) Suspicious specification for a boolean: "<< buf <<endl;
		cout<<"(!) Only \"\", \"T\", \"True\", \"Yes\", \"Y\", \"F\", \"False\", \"No\", and \"N\" are allowed"<<endl;
		cout<<"I exit..."<< endl;
		if (cinGetOnError) cin.get();
		exit(-1);
    }
return(0);
}

//	Initializing boolean locus/marker parameter vectors from the settings file "SettingsFilename"
void evaluateBoolVector(string buf, vector<bool>& varBoolVector, string displayVar) {
    if (cmp_nocase(buf,"") == 0 || cmp_nocase(buf,"T") == 0 ||
          cmp_nocase(buf,"True") == 0 || cmp_nocase(buf,"Yes") == 0 ||
          cmp_nocase(buf,"Y") == 0)
    	varBoolVector.push_back(true);
    else if (cmp_nocase(buf,"F") == 0 || cmp_nocase(buf,"False") == 0 ||
          cmp_nocase(buf,"No") == 0 || cmp_nocase(buf,"N") == 0)
    	varBoolVector.push_back(false);
    else {
		cerr << "\nUnrecognized boolean argument(s) for " << displayVar << " in your settings file" << SettingsFilename;
		cerr << "\nOnly \"\", \"T\", \"True\", \"Yes\", \"Y\", \"F\", \"False\", \"No\", and \"N\" are allowed.\n";
		cerr << "\nAborting IBDSim...Press any key to exit.\n";
		if (cinGetOnError)
			cin.get();
		exit(-1);
    }
}

//  transforms unspecified/partly specified/marker-wise specified parameter values into a locus-wise specification
template <typename T>
void expandParamVector(vector<T> &paramVector, string dispVar) {

//	IF the parameter vector is of type "int" type cast the void pointer correspondingly.
//	The following code is the same for "double", "string" and "bool" cases below.

	size_t vectorSize = paramVector.size();

	//	IF the particular parameter isn't declared in the settings file, use the respective default value found at the zeroth index.
	if (vectorSize == 1)
		paramVector.assign(n_locus + 1, paramVector[0]);
	//	ELSE IF it contains a single value, use the same for the total number of simulated loci irrespective of their marker type.
	else if (vectorSize == 2)
		paramVector.assign(n_locus + 1, paramVector[1]);
	//	ELSE IF the parameter values have been declared marker-wise, use them for all loci of the respective marker.
	else if ((vectorSize - 1) == (size_t)nMarker)
		for (size_t i = vectorSize; i >= 2; i--){
            const T partToCopy=paramVector[i-1];
            if(locusNumVector[i - 1]>0) paramVector.insert((paramVector.begin() + i), size_t (locusNumVector[i - 1] - 1), partToCopy);
        }
	//	ELSE unspecified format, report it and quit IBDSim.
	else {
		cerr << "\nWrong number of arguments for \"" << dispVar << "\". \nPlease check your settings file " << SettingsFilename
				<< endl << "Aborting IBDSim...Press any key to exit.";
		if (cinGetOnError)
			cin.get();
		exit(-1);
	}

	//	FOR DEBUGGING PURPOSES ONLY: code for outputting parameter contents locus-wise
    if(false) {
        cout << endl << dispVar << endl;
        cout << "First model is the default (=KAM with 10 states) and is not used for simulation" << endl;
        for (typename vector<T>::iterator it = paramVector.begin(); it < paramVector.end(); it++)
            cout << *it << " ";
        cout << endl;
    }
}


//Lecture des paramtres dans fichier
int read_settings_file(const char filename[]) {
    string buf,var;
//stringstream good mostly for string/numbers conversions
    size_t pos;
    int tempo;
    char bidon;
    bool hacktempo=true;

    ifstream settings(filename,ios::in);

    if(!settings.is_open()) {
        cerr << "Unable to open file "<<filename << " in read_settings_file();" << endl;
        cerr<<" Check settings. I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
	} else do {
		getline(settings,buf);
//        cout << "buf=" << buf << endl;
		if(buf.length()==0) goto nextline;
		while((buf[0]==' ')||(buf[0]=='\t')||(buf[0]=='\r')) {buf.erase(0,1);}//vire les blancs initiaux
        while((buf[buf.length()-1]==' ')||(buf[buf.length()-1]=='\t')||(buf[buf.length()-1]=='\r')) {buf.erase(buf.length()-1,buf.length());}//vire les blancs finaux
		pos=std::min(buf.find('='),std::min(buf.find('\t'),buf.length()));
//		pos=std::min(buf.find('='),buf.length());
        if ((buf[pos])=='=') while (buf[pos-1]==' '||(buf[pos-1]=='\t')) {buf.erase(pos-1,1); pos--;}// vire les blancs et les tabulations avant le =
        if ((buf[pos])=='=') while (buf[pos+1]==' '||(buf[pos+1]=='\t')) {buf.erase(pos+1,1);}// vire les blancs et les tabulations apres le =
		var=buf.substr(0,pos).c_str();
		if(var.length()==0) goto nextline;
		if(var[0]=='/'||var[0]=='%'||var[0]=='#') goto nextline;//RL pour ne pas afficher les commentaires

#ifdef GOTO
_gotoxy(0,17);
#endif
        if (var[var.size()-1]=='0') {
            cerr<<"Automatically removing final '0' from "<<var<< ".\nPress any key to resume." << endl;
            var=var.substr(0,var.size()-1);
            if (cinGetOnError) cin.get();
        }
		if(cmp_nocase_no__(var,"Pause")==0) {
/* Pause determines only two correct contexts for cin.get():
       cerr<< error message + if(cinGetOnError) cin.get() + exit
and
       cout<< some info + if(pauseGP) cin.get() + execution continues
cinGetOnError is true at declaration in genepop.cpp but may then be set to false in latin.cpp
*/
			string locstring;
			stringstream strstr(buf.substr(pos+1));
			strstr>>locstring;
			if(cmp_nocase_no__(buf.substr(pos+1),"Final")==0) pauseGP=true;
			if(cmp_nocase_no__(locstring,"OnError")==0) cinGetOnError=true;
			if( cmp_nocase_no__(locstring,"Default")==0 || cmp_nocase_no__(locstring,"Never")==0) {cinGetOnError=false; pauseGP=false;}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Data_File_Name")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>fichier_genepop;
			goto nextline;
		}
		if(cmp_nocase_no__(var,".txt_extension")==0) {
			evaluateBool(txt_extensionbool,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Genepopfile_extension")==0 || cmp_nocase_no__(var,"Genepop_File_Extension")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>Genepopfile_extension;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Run_Number")==0 || cmp_nocase_no__(var,"Dataset_Number")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>repet;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Random_Seeds")==0 || cmp_nocase_no__(var,"Random_Seed")==0) {
			string reste=buf.substr(pos+1);
			stringstream strstr(reste);
			//int i=0;
			while (!strstr.eof()) {
				strstr>>tempo;
				/*strstr.get(); to skip any one-character separator which is necessarily here
				is not sufficient is there are whitespaces after the last value. If so, eof is not reached
				but no new value is read -> the last value is duplicated */
                 while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.') strstr.get();
                 Seeds=tempo;
                 //i++;
			}
			strstr.clear(); // c'est le truc essentiel pour le reutiliser... snif
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Locus_Number") == 0) {
			buf = buf.substr(pos + 1);
			string arg;
			pos = 0;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !isdigit(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && isdigit(buf[pos]))
				   pos++;
			   arg = buf.substr(0,pos);
			   locusNumVector.push_back(atoi(arg.c_str()));
               if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Mutation_Rate") == 0) {
			buf = buf.substr(pos + 1);
			pos = 0;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !(buf[pos] == '.') && !isdigit(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && (buf[pos] == '.' || isdigit(buf[pos])))
				   pos++;
			   string strstr(buf.substr(0,pos));
			   double stringToFloat = atof(strstr.c_str());
			   mutRateVector.push_back(stringToFloat);
               if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
		if (cmp_nocase_no__(var,"Mutation_Model") == 0) {
			buf = buf.substr(pos + 1);
			string arg;
			pos = 0;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !isalnum(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && isalnum(buf[pos]))
				   pos++;
			   arg = buf.substr(0,pos);
			   mutModelIDVector.push_back(arg);
			   if (!((cmp_nocase_no__(mutModelIDVector.back(),"SMM") == 0)
				   || (cmp_nocase_no__(mutModelIDVector.back(),"IAM") == 0)
				   || (cmp_nocase_no__(mutModelIDVector.back(),"KAM") == 0)
				   || (cmp_nocase_no__(mutModelIDVector.back(),"TPM") == 0)
				   || (cmp_nocase_no__(mutModelIDVector.back(),"GSM") == 0)
				   || (cmp_nocase_no__(mutModelIDVector.back(),"JC69") == 0)
				   || (cmp_nocase_no__(mutModelIDVector.back(),"K80") == 0)
                   || (cmp_nocase_no__(mutModelIDVector.back(),"K2P") == 0)
				   || (cmp_nocase_no__(mutModelIDVector.back(),"F81") == 0)
				   || (cmp_nocase_no__(mutModelIDVector.back(),"HKY85") == 0)
				   || (cmp_nocase_no__(mutModelIDVector.back(),"TN93") == 0)
				   || (cmp_nocase_no__(mutModelIDVector.back(),"SNP") == 0)
				   || (cmp_nocase_no__(mutModelIDVector.back(),"ISM") == 0))) {
				   cerr << "\nUnrecognized Mutation_Model argument(s) in " << filename;
				   cerr << "\nCurrently only a combination of SMM, IAM, KAM, TPM, GSM, JC69, K80/K2P, "
						   "F81, HKY85, TN93, SNP or ISM models are allowed.\n";
				   cerr << "\nAborting IBDSim...Press any key to exit.\n";
				   if (cinGetOnError)
					   cin.get();
				   exit(-1);
			   }
			   if (cmp_nocase_no__(mutModelIDVector.back(),"ISM") == 0
					|| cmp_nocase_no__(mutModelIDVector.back(),"SNP") == 0
					|| cmp_nocase_no__(mutModelIDVector.back(),"JC69") == 0
					|| cmp_nocase_no__(mutModelIDVector.back(),"K80") == 0
                    || cmp_nocase_no__(mutModelIDVector.back(),"K2P") == 0
					|| cmp_nocase_no__(mutModelIDVector.back(),"F81") == 0
					|| cmp_nocase_no__(mutModelIDVector.back(),"HKY85") == 0
					|| cmp_nocase_no__(mutModelIDVector.back(),"TN93") == 0)
				   specifyLocusSeq = true;
                if(cmp_nocase_no__(mutModelIDVector.back(),"K2P") == 0) mutModelIDVector.back()="K80";

                if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Min_Allele_Number") == 0) {
			buf = buf.substr(pos + 1);
			string arg;
			pos = 0;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !isalnum(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && isalnum(buf[pos]))
				   pos++;
			   arg = buf.substr(0,pos);
			   if (cmp_nocase_no__(arg,"NA") == 0)
				   arg = "1";
			   minAlleleVector.push_back(atoi(arg.c_str()));
               if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
        if(cmp_nocase_no__(var,"Max_Mutation_Number") == 0 || cmp_nocase_no__(var,"Max_Mutations_Number") == 0) {
			buf = buf.substr(pos + 1);
			string arg;
			pos = 0;
            while (buf.length() > 0) {
                while (buf.length() > pos && !isalnum(buf[pos])) {
                    if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
                        cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
                        << filename << "! Ignored." << endl;
                        if (!warningBool)
                            warningBool = true;
                    }
                    pos++;
                }
                if (buf.length() == pos)
                    goto nextline;
                buf = buf.substr(pos);
                pos = 0;
                while (buf.length() > pos && isalnum(buf[pos]))
                    pos++;
                arg = buf.substr(0,pos);
                if (cmp_nocase_no__(arg,"NA") == 0)
                    arg = "1";
                maxMutationsVector.push_back(atoi(arg.c_str()));
                if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
                } else goto nextline;
			}
			goto nextline;
		}

		if(cmp_nocase_no__(var,"Polymorphic_Loci_only") == 0 || cmp_nocase_no__(var,"Polymorphic_Locus_only") == 0
				|| cmp_nocase_no__(var,"Polymorphic_Loci") == 0 || cmp_nocase_no__(var,"Polymorphic_Locus") == 0) {
			buf = buf.substr(pos + 1);
			string arg;
			pos = 0;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !isalpha(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && isalpha(buf[pos]))
				   pos++;
			   arg = buf.substr(0,pos);
			   evaluateBoolVector(arg, polyLociBoolVector, var);
               if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Minor_Allele_Frequency") == 0) {
			buf = buf.substr(pos + 1);
			pos = 0;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !(buf[pos] == '.') && !isalnum(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && (buf[pos] == '.' || isalnum(buf[pos])))
				   pos++;
			   string strstr(buf.substr(0,pos));
			   if (cmp_nocase_no__(strstr,"NA") == 0)
				   strstr = "0";
			   double stringToFloat = atof(strstr.c_str());
			   minorAlleleFreqVector.push_back(stringToFloat);
               if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Variable_Mutation_Rate") == 0) {
			buf = buf.substr(pos + 1);
			string arg;
			pos = 0;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !isalpha(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && isalpha(buf[pos]))
				   pos++;
			   arg = buf.substr(0,pos);
			   evaluateBoolVector(arg, varMutRateBoolVector, var);
               if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Allelic_Lower_Bound") == 0 || cmp_nocase_no__(var,"Min_Allele") == 0) {
			buf = buf.substr(pos + 1);
			string arg;
			pos = 0;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !isalnum(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && isalnum(buf[pos]))
				   pos++;
			   arg = buf.substr(0,pos);
			   if (cmp_nocase_no__(arg,"NA") == 0)
				   arg = "-1";
			   kMinVector.push_back(atoi(arg.c_str()));
               if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Allelic_Upper_Bound") == 0 || cmp_nocase_no__(var,"Max_Allele") == 0) {
			buf = buf.substr(pos + 1);
			string arg;
			pos = 0;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !isalnum(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && isalnum(buf[pos]))
				   pos++;
			   arg = buf.substr(0,pos);
			   if (cmp_nocase_no__(arg,"NA") == 0)
				   arg = "-1";
			   kMaxVector.push_back(atoi(arg.c_str()));
               if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Allelic_State_MRCA") == 0) {
			buf = buf.substr(pos + 1);
			string arg;
			pos = 0;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !isalnum(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && isalnum(buf[pos]))
				   pos++;
			   arg = buf.substr(0,pos);
			   if (cmp_nocase_no__(arg,"NA") == 0)
				   arg = "-1";
			   kIniVector.push_back(atoi(arg.c_str()));
               if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Repeated_motif_size") == 0) {
			buf = buf.substr(pos + 1);
			string arg;
			pos = 0;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !isalnum(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && isalnum(buf[pos]))
				   pos++;
			   arg = buf.substr(0,pos);
			   if (cmp_nocase_no__(arg,"NA") == 0)
				   arg = "-1";
			   motifSizeVector.push_back(atoi(arg.c_str()));
               if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"SMM_Probability_In_TPM") == 0) {
			buf = buf.substr(pos + 1);
			pos = 0;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !(buf[pos] == '.') && !isalnum(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && (buf[pos] == '.' || isalnum(buf[pos])))
				   pos++;
			   string strstr(buf.substr(0,pos));
			   if (cmp_nocase_no__(strstr,"NA") == 0)
				   strstr = "-1";
			   double stringToFloat = atof(strstr.c_str());
			   pSMMVector.push_back(stringToFloat);
               if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Geometric_Variance_In_TPM") == 0) {
			buf = buf.substr(pos + 1);
			pos = 0;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !(buf[pos] == '.') && !isalnum(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && (buf[pos] == '.' || isalnum(buf[pos])))
				   pos++;
			   string strstr(buf.substr(0,pos));
			   if (cmp_nocase_no__(strstr,"NA") == 0)
				   strstr = "-1";
			   double stringToFloat = atof(strstr.c_str());
			   geomTPMVector.push_back(stringToFloat);
               if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Geometric_Variance_In_GSM") == 0) {
			buf = buf.substr(pos + 1);
			pos = 0;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !(buf[pos] == '.') && !isalnum(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && (buf[pos] == '.' || isalnum(buf[pos])))
				   pos++;
			   string strstr(buf.substr(0,pos));
			   if (cmp_nocase_no__(strstr,"NA") == 0)
				   strstr = "-1";
			   double stringToFloat = atof(strstr.c_str());
			   geomGSMVector.push_back(stringToFloat);
               if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"MRCA_Sequence")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>defMRCASequence;
			if (defMRCASequence.size())
				seqSize = (int) defMRCASequence.size();
            goto nextline;
		}
		if(cmp_nocase_no__(var,"Sequence_Size")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>seqSize;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Transition_Transversion_ratio")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>ratio_TITV;
		    if (ratio_TITV < 0) {
		    	cerr << "\nThe \"Transition_Transversion_ratio\" should always be positive. \nPlease check "
					 << SettingsFilename << "\nAborting IBDSim...Press any key to exit.\n";
				if (cinGetOnError)
					cin.get();
				exit(-1);
		    }
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Transition1_Transition2_ratio")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>ratio_TITI;
		    if (ratio_TITI < 0) {
		    	cerr << "\nThe \"Transition1_Transition2_ratio\" should always be positive. \nPlease check "
					 << SettingsFilename << "\nAborting IBDSim...Press any key to exit.\n";
				if (cinGetOnError)
					cin.get();
				exit(-1);
		    }
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Equilibrium_Frequencies") == 0) {
			int arrayIndex = 0;
			buf = buf.substr(pos + 1);
			pos = 0;
			specifyVarBaseFreqs = true;
            while (buf.length() > 0) {
			   while (buf.length() > pos && !(buf[pos] == '.') && !isdigit(buf[pos])) {
				   if ((buf[pos] != ' ') && (buf[pos] != ',') && (buf[pos] != '\t')) {
					   cerr << "Unidentified specification \"" << buf[pos] << "\" for \"" << var << "\" in "
					   	    << filename << "! Ignored." << endl;
					   if (!warningBool)
						   warningBool = true;
				   }
				   pos++;
			   }
			   if (buf.length() == pos)
				   goto nextline;
			   buf = buf.substr(pos);
			   pos = 0;
			   while (buf.length() > pos && (buf[pos] == '.' || isdigit(buf[pos])))
				   pos++;
			   string strstr(buf.substr(0,pos));
			   double stringToFloat = atof(strstr.c_str());
			   if (arrayIndex == 4) {
					cerr << "\nError in the number of parameters for \"Equilibrium_Frequencies\". Please check "
						 << SettingsFilename << "\nAborting IBDSim...Press any key to exit.\n";
					if (cinGetOnError)
						cin.get();
					exit(-1);
			   }
			   varBaseFreqs[arrayIndex] = stringToFloat;
			   ++arrayIndex;
               if (buf.length() > pos) {
    				buf = buf.substr(pos+1);
    			   	pos = 0;
			   } else goto nextline;
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Nexus_file_format") == 0) {
			string locstring;
			stringstream strstr(buf.substr(pos+1));
			strstr >> locstring;
			if(cmp_nocase_no__(locstring,"Haplotypes_and_Individuals") == 0)
				nexusFormat = locstring;
			else if(cmp_nocase_no__(locstring,"Haplotypes_only") == 0)
				nexusFormat = locstring;
			else if(cmp_nocase(locstring,"F")==0 || cmp_nocase(locstring,"False")==0 ||
			          cmp_nocase(locstring,"No")==0 || cmp_nocase(locstring,"N")==0)
				nexusFormat="F";
		    else {
				cerr << "\nUnrecognized argument for Haplotype_format in " << SettingsFilename;
				cerr << "\nOnly \"Haplotypes_only\", \"Haplotypes_and_Individuals\" and \"False\"are allowed.\n";
				cerr << "\nAborting IBDSim...Press any key to exit.\n";
				if (cinGetOnError)
					cin.get();
				exit(-1);
		    }
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Ploidy")==0) {
			string locstring;
			stringstream strstr(buf.substr(pos+1));
			strstr>>locstring;
			if(cmp_nocase_no__(locstring,"Diploid")==0) ploidy=2;
			else if(cmp_nocase_no__(locstring,"Haploid")==0) ploidy=1;
				else {cerr<<"\nBad Argument for Ploidy in " << SettingsFilename <<  " :"<<endl;
					  cerr<<"Only Diploid or Haploid are implemented"<<endl;
					  if(cinGetOnError) cin.get();
					  exit(-1);
				}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"genepop")==0 || cmp_nocase_no__(var,"genepopfile")==0) {
			evaluateBool(genepopoui,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"genepop_no_coord")==0) {
			evaluateBool(genepopNoCoordbool,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"geneland")==0 || cmp_nocase_no__(var,"genelandfile")==0) {
			evaluateBool(genelandoui,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Migraine")==0) {
			evaluateBool(DG2002,buf.substr(pos+1)); //FR->RL obsolete... à maintenir pour certaines comparaisons.
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Migraine_Settings")==0) {
			evaluateBool(migraineSettingsbool,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Migraine_AllStates")==0) {
			evaluateBool(AllStates,buf.substr(pos+1));  // colonne des 0 aux format DG2002 (obsolete)
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Migrate")==0 || cmp_nocase_no__(var,"MigrateFile")==0) {
			evaluateBool(migrate_oui,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Migrate_Letter")==0 || cmp_nocase_no__(var,"Migrate_Lettre")==0) {
			evaluateBool(migrate_lettre_oui,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Generic_Computations")==0) {
			evaluateBool(calculoui,buf.substr(pos+1));  //FR->RL remove or rename
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Hexp")==0 || cmp_nocase_no__(var,"Hexp_Nei")==0) {
			evaluateBool(HexpNeioui,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Seq_stats")==0) {
			evaluateBool(seqStatsOui,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"DeltaH")==0) {
            cerr<<"DeltaH is obsolete..."<<endl;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Allelic_Variance")==0) { // calcul variance...
			evaluateBool(Varoui,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Iterative_Identity_Probability")==0) {
			evaluateBool(suiviQ,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Iterative_Statistics")==0) {
			evaluateBool(iterativeStats,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Allelic_Number_Folder")==0) { // obsolete
			string locstring;
			stringstream strstr(buf.substr(pos+1));
			strstr>>locstring;

			if(cmp_nocase(locstring,"true")==0) {
                cout<<"\nAllelic_Number_Folder is no more implemented :"<<endl;
				cout<<"Contact R. leblois if you really need this option"<<endl;
				getchar();
				exit(-1);
                //dossier_nbre_allele=true;
            } else if(cmp_nocase(locstring,"false")==0) dossier_nbre_allele=false;
				else {cout<<"\nBad Argument for Allelic_Number_Folder in " << SettingsFilename <<  " :"<<endl;
					  cout<<"Only true or false are possible options"<<endl;
				 	  getchar();
					  exit(-1);
				}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Old_Version_Computation")==0) { // obsolete
			string locstring;
			stringstream strstr(buf.substr(pos+1));
			strstr>>locstring;
			if(cmp_nocase(locstring,"true")==0) calcul_anc=true;
			else if(cmp_nocase(locstring,"false")==0) calcul_anc=false;
				else {cout<<"\nBad Argument for Old_Version_Computation in " << SettingsFilename <<  " :"<<endl;
					  cout<<"Only true or false are possible options"<<endl;
				 	  getchar();
					  exit(-1);
				}
			goto nextline;
		}

		if(cmp_nocase_no__(var,"DiagnosticTables")==0 || cmp_nocase_no__(var,"DiagnosticTable")==0) { // FR->RL new mechanism to pass options to R in a more economical way that all the booleans in Roptions
			buf=buf.substr(pos+1);
			string arg;
			pos=0;
            while (buf.length()>0) {
			   while (buf.length()>pos && !(buf[pos]=='_') && !isalnum(buf[pos])) pos++;
			   if (buf.length()==pos) goto nextline;
			   buf=buf.substr(pos);
			   pos=0;
			   while (buf.length()>pos && (buf[pos]=='_' || isalnum(buf[pos]))) pos++;
			   arg=buf.substr(0,pos);
			   if (cmp_nocase_no__(arg,"Effective_Dispersal")==0 || cmp_nocase_no__(arg,"Realized_Dispersal")==0) {effective_disp=true;} //
			   else if (cmp_nocase_no__(arg,"Prob_Id_Matrix")==0) {Prob_Id_Matrix=true;} //
			   else if (cmp_nocase_no__(arg,"Iterative_Statistics")==0) {iterativeStats=true;} //
			   else if (cmp_nocase_no__(arg,"Iterative_Identity_Probability")==0) {suiviQ=true;} //
			   else if (cmp_nocase_no__(arg,"Allelic_Variance")==0) {Varoui=true;} //
			   else if (cmp_nocase_no__(arg,"Hexp")==0 || cmp_nocase_no__(arg,"Hexp_Nei")==0) {HexpNeioui=true;} //
			   else if (cmp_nocase_no__(arg,"Seq_stats")==0) {seqStatsOui=true;} //
			   else if (cmp_nocase_no__(arg,"Fis")==0) {Fisoui=true;} //

			   else {
                   cerr<<"Unknown DiagnosticTables arguments'"<<arg<<"'. Check settings. I exit."<<endl;
                   if (cinGetOnError) cin.get();
                   exit(-1);
                }
               if (buf.length()>pos) {
    				buf=buf.substr(pos+1);
    			   	pos=0;
			   } else goto nextline;
			}
			goto nextline;
		}
		
		if(cmp_nocase_no__(var,"Prob_Id_Matrix")==0) {
			evaluateBool(Prob_Id_Matrix,buf.substr(pos+1));
			goto nextline;
		}

		if(cmp_nocase_no__(var,"Effective_Dispersal")==0 || cmp_nocase_no__(var,"Realized_Dispersal")==0) {
			evaluateBool(effective_disp,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Constant_Dispersal")==0) {
			evaluateBool(const_disp,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Immigration_Control")==0 || cmp_nocase_no__(var,"2D_Dispersal")==0 || cmp_nocase_no__(var,"2D_Disp")==0) {
		    //FR->FR FR->RL trouver autre nom car option a un sens en 1D
			string locstring;
			stringstream strstr(buf.substr(pos+1));
			strstr>>locstring;
			if(cmp_nocase_no__(locstring,"Simple1DProduct")==0) {}// was twoD_disp="Simple1DProduct";
			else if(cmp_nocase_no__(locstring,"1DProductWithoutm0")==0) Simple1DProductDispBool=false; // was twoD_disp="1DProductWithoutm0";
			//else if(cmp_nocase(locstring,"X+YdistProb")==0) twoD_disp="X+YDistProb";
			else {cerr<<"\nBad Argument for 2D_Dispersal in " << SettingsFilename <<  " :"<<endl;
				cerr<<"Only Simple1DProduct and 1DProductWithoutm0 are possible options. I exit."<<endl;
			    if(cinGetOnError) cin.get();
				exit(-1);
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Total_Range_Dispersal")==0) {
			evaluateBool(compareMath,buf.substr(pos+1));
			/** if true: dx_max est contraint par la fonction de dispersion
			 else dx_max et la dim max du reseau**/
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Lattice_Boundaries")==0 || cmp_nocase_no__(var,"Edge_Effect")==0 || cmp_nocase_no__(var,"Edge_Effects")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>boundaryType;
            goto nextline;
		}
		if(cmp_nocase_no__(var,"Max_Lattice_SizeX")==0) {
            cerr<<"Max_Lattice_SizeX is obsolete..."<<endl;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Max_Lattice_SizeY")==0) {
            cerr<<"Max_Lattice_SizeY is obsolete..."<<endl;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Sample_SizeX")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>dim_sample1[0];
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Sample_SizeY")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>dim_sample2[0];
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Min_Sample_CoordinateX")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>xmin_sample[0];
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Min_Sample_CoordinateY")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>ymin_sample[0];
			goto nextline;
		}
		if(cmp_nocase_no__(var,"SpecificSampleDesign_SampleSize")==0) {
            cerr<<"SpecificSampleDesign_SampleSize is obsolete..."<<endl;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Specific_Sample_Design")==0) {
            cerr<<"Specific_Sample_Design is obsolete..."<<endl;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Sample_Coordinates_X")==0) {
            string reste=buf.substr(pos+1);
			stringstream strstr(reste);
			Spec_Sample_Coord[0][x].resize(0);
			while (!strstr.eof()) {
					strstr>>tempo;
					/*strstr.get(); to skip any one-character separator which is necessarily here
					is not sufficient is there are whitespaces after the last value. If so, eof is not reached
					but no new value is read -> the last value is duplicated */
	                 while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.') strstr.get();
	                 Spec_Sample_Coord[0][x].push_back(tempo);
			}
			vector<int>::iterator it=min_element(Spec_Sample_Coord[0][x].begin(), Spec_Sample_Coord[0][x].end());
			xmin_sample[0]=(*it);
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Sample_Coordinates_Y")==0) {
				string reste=buf.substr(pos+1);
				stringstream strstr(reste);
				Spec_Sample_Coord[0][y].resize(0);
				while (!strstr.eof()) {
					strstr>>tempo;
					/*strstr.get(); to skip any one-character separator which is necessarily here
					is not sufficient is there are whitespaces after the last value. If so, eof is not reached
					but no new value is read -> the last value is duplicated */
	                 while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.') strstr.get();
	                 Spec_Sample_Coord[0][y].push_back(tempo);
				}
				vector<int>::iterator it=min_element(Spec_Sample_Coord[0][y].begin(), Spec_Sample_Coord[0][y].end());
				ymin_sample[0]=(*it);
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Group_All_Samples")==0 || cmp_nocase_no__(var,"Group_Samples")==0) {
			evaluateBool(groupAllSamplesbool,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Ind_Per_Pop_Sampled")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>dens_sample[0];
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Void_Sample_Node")==0 || cmp_nocase_no__(var,"Void_Sample_Nodes")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>vide_sample[0];
			goto nextline;
		}
        // pour échantillon prédispersion
        if(cmp_nocase_no__(var,"Predisp_Sample_SizeX")==0) {
            stringstream strstr(buf.substr(pos+1));
            strstr>>dim_sample1[1];
            goto nextline;
        }
        if(cmp_nocase_no__(var,"Predisp_Sample_SizeY")==0) {
            stringstream strstr(buf.substr(pos+1));
            strstr>>dim_sample2[1];
            goto nextline;
        }
        if(cmp_nocase_no__(var,"Predisp_Min_Sample_CoordinateX")==0) {
            stringstream strstr(buf.substr(pos+1));
            strstr>>xmin_sample[1];
            goto nextline;
        }
        if(cmp_nocase_no__(var,"Predisp_Min_Sample_CoordinateY")==0) {
            stringstream strstr(buf.substr(pos+1));
            strstr>>ymin_sample[1];
            goto nextline;
        }
        if(cmp_nocase_no__(var,"Predisp_SpecificSampleDesign_SampleSize")==0) {
            cerr<<"Predisp_SpecificSampleDesign_SampleSize is obsolete..."<<endl;
            goto nextline;
        }
        if(cmp_nocase_no__(var,"Predisp_Specific_Sample_Design")==0) {
            cerr<<"Predisp_Specific_Sample_Design is obsolete..."<<endl;
            goto nextline;
        }
        if(cmp_nocase_no__(var,"Predisp_Sample_Coordinates_X")==0) {
            string reste=buf.substr(pos+1);
            stringstream strstr(reste);
            Spec_Sample_Coord[1][x].resize(0);
            while (!strstr.eof()) {
                strstr>>tempo;
                /*strstr.get(); to skip any one-character separator which is necessarily here
                 is not sufficient is there are whitespaces after the last value. If so, eof is not reached
                 but no new value is read -> the last value is duplicated */
                while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.') strstr.get();
                Spec_Sample_Coord[1][x].push_back(tempo);
            }
            vector<int>::iterator it=min_element(Spec_Sample_Coord[1][x].begin(), Spec_Sample_Coord[1][x].end());
            xmin_sample[1]=(*it);
            goto nextline;
        }
        if(cmp_nocase_no__(var,"Predisp_Sample_Coordinates_Y")==0) {
            string reste=buf.substr(pos+1);
            stringstream strstr(reste);
            Spec_Sample_Coord[1][y].resize(0);
            while (!strstr.eof()) {
                strstr>>tempo;
                /*strstr.get(); to skip any one-character separator which is necessarily here
                 is not sufficient is there are whitespaces after the last value. If so, eof is not reached
                 but no new value is read -> the last value is duplicated */
                while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.') strstr.get();
                Spec_Sample_Coord[1][y].push_back(tempo);
            }
            vector<int>::iterator it=min_element(Spec_Sample_Coord[1][y].begin(), Spec_Sample_Coord[1][y].end());
            ymin_sample[1]=(*it);
            goto nextline;
        }
        if(cmp_nocase_no__(var,"Predisp_Ind_Per_Pop_Sampled")==0) {
            stringstream strstr(buf.substr(pos+1));
            strstr>>dens_sample[1];
            goto nextline;
        }
        if(cmp_nocase_no__(var,"Predisp_Void_Sample_Node")==0 || cmp_nocase_no__(var,"Predisp_Void_Sample_Nodes")==0) {
            stringstream strstr(buf.substr(pos+1));
            strstr>>vide_sample[1];
            goto nextline;
        }
        //fin échantillon prédispersion
        
        if(cmp_nocase_no__(var,"NewDemographicPhaseAt")==0) {
			unsigned long int locMinGen;
			CTimeVaryingParams locCTP;
			stringstream strstr(buf.substr(pos+1));
			strstr>>locMinGen;
			if (locMinGen==0) {
                cerr<<"Found New demographic phase starting at 0." <<endl;
                cerr<<"generations in IOBDSim starts at 1, Timing of the New demographic phase is changed to 1 instead of 0." <<endl;
                if(cinGetOnError) cin.get();
			}
			if (locMinGen>1) {
			    if (locMinGen<=TVpars.back().minGen) {
                    cerr<<"New demographic phase must begin at later time than previous one." <<endl;
                    cerr<<"Check use of NewDemographicPhaseAt setting. I exit." <<endl;
                    if(cinGetOnError) cin.get();
                    exit(-1);
                }
                TVpars.push_back(CTimeVaryingParams(&(TVpars.back())));
                TVpars.back().minGen=locMinGen;
                FoundDispSpatiallyHeterogInSettings.resize(TVpars.size());
			}
            //remplissage du tableau densite plus loin car besoin de DimResX0 et Y0
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Ind_Per_Pop")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().initialDens;
			TVpars.back().dens=TVpars.back().initialDens;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Lattice_SizeX")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().initialDimRes1;
			TVpars.back().dimRes1=TVpars.back().initialDimRes1;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Lattice_SizeY")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().initialDimRes2;
			TVpars.back().dimRes2=TVpars.back().initialDimRes2;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Void_Nodes")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().vide;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Zone")==0) {
            evaluateBool(TVpars.back().zone,buf.substr(pos+1));
            if(TVpars.size()>1 && !(TVpars.back().zone) && TVpars[TVpars.size()-2].zone) {TVpars.back().dens_zone=-666;TVpars.back().vide_zone=-666;}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Min_Zone_CoordinateX")==0) { //FR->RL on laisse sans default pour l'instant jusqua Specific_Density_Design (RL??)
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().xmin_zone;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Max_Zone_CoordinateX")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().xmax_zone;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Min_Zone_CoordinateY")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().ymin_zone;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Max_Zone_CoordinateY")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().ymax_zone;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Specific_Density_Design")==0 || cmp_nocase_no__(var,"SpecificDensityDesign")==0) {
			evaluateBool(Specific_Density_Designbool,buf.substr(pos+1));
            //remplissage du tableau densite plus loin car besoiun de DimResX0 et Y0
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Void_Nodes_Zone")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().vide_zone;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Ind_Per_Pop_Zone")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().dens_zone;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"RealForwardBarrier")==0 || cmp_nocase_no__(var,"ForwardBarrier")==0) {
            evaluateBool(TVpars.back().forwardBarrier,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Barrier")==0 || cmp_nocase_no__(var,"BackwardBarrier")==0) {
            evaluateBool(TVpars.back().barrier,buf.substr(pos+1));
			goto nextline;
		}
		if(cmp_nocase_no__(var,"X1_Barrier")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().x1_barrier;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"X2_Barrier")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().x2_barrier;
			goto nextline;
		}if(cmp_nocase_no__(var,"Y1_Barrier")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().y1_barrier;
			goto nextline;
		}if(cmp_nocase_no__(var,"Y2_Barrier")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().y2_barrier;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"BarrierCrossingRate")==0 || cmp_nocase_no__(var,"CrossingRate")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().barrierCrossingRate;
			goto nextline;
		}		
		if(cmp_nocase_no__(var,"Dispersal_Distribution")==0) {
			stringstream strstr(buf.substr(pos+1));
			if(strstr.str().length()==1) {
				string locstring;
				//stringstream strstr1;
				//strstr1 << strstr.rdbuf();
				strstr>>locstring;
				if(cmp_nocase_no__(locstring,"g")==0) TVpars.back().mod='g';
				else if(cmp_nocase_no__(locstring,"s")==0) TVpars.back().mod='S';
				else if(cmp_nocase_no__(locstring,"p")==0) TVpars.back().mod='P';
				else TVpars.back().mod=locstring.c_str()[0];
			} else {
				string locstring;
				strstr>>locstring;
				if(cmp_nocase_no__(locstring,"steppingStone")==0 ||
                   cmp_nocase_no__(locstring,"SS")==0) TVpars.back().mod='b';
				else if(cmp_nocase_no__(locstring,"Geometric")==0) TVpars.back().mod='g';
				else if(cmp_nocase_no__(locstring,"Sichel")==0) TVpars.back().mod='S';
				else if(cmp_nocase_no__(locstring,"Pareto")==0) TVpars.back().mod='P';
				
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Total_Emigration_Rate")==0 || cmp_nocase_no__(var,"Total_Emmigration_Rate")==0) {
            stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().Mig;
					goto nextline;
		}
		if(cmp_nocase_no__(var,"dist_max")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().dist_max;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Pareto_Shape")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().Pareto_Shape;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Geometric_Shape")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().GeoG;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Sichel_Gamma")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().Sichel_Gamma;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Sichel_Xi")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().Sichel_Xi;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Sichel_Omega")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().Sichel_Omega;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Dispersal_Distribution_zone")==0) {
			stringstream strstr(buf.substr(pos+1));
			if(strstr.str().length()==1) {
				string locstring;
				//stringstream strstr1;
				//strstr1 << strstr.rdbuf();
				strstr>>locstring;
				if(cmp_nocase_no__(locstring,"g")==0) TVpars.back().mod_zone='g';
				else if(cmp_nocase_no__(locstring,"s")==0) TVpars.back().mod_zone='S';
				else if(cmp_nocase_no__(locstring,"p")==0) TVpars.back().mod_zone='P';
				else TVpars.back().mod_zone=locstring.c_str()[0];
			} else {
				string locstring;
				strstr>>locstring;
				if(cmp_nocase_no__(locstring,"steppingStone")==0) TVpars.back().mod_zone='b';
				else if(cmp_nocase_no__(locstring,"Geometric")==0) TVpars.back().mod_zone='g';
				else if(cmp_nocase_no__(locstring,"Sichel")==0) TVpars.back().mod_zone='S';
				else if(cmp_nocase_no__(locstring,"Pareto")==0) TVpars.back().mod_zone='P';
				
			}
            FoundDispSpatiallyHeterogInSettings[TVpars.size()-1]=true;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Total_Emigration_Rate_zone")==0 || cmp_nocase_no__(var,"Total_Emmigration_Rate_zone")==0) {
            stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().Mig_zone;
            FoundDispSpatiallyHeterogInSettings[TVpars.size()-1]=true;
            goto nextline;
		}
		if(cmp_nocase_no__(var,"dist_max_zone")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().dist_max_zone;
            FoundDispSpatiallyHeterogInSettings[TVpars.size()-1]=true;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Pareto_Shape_zone")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().Pareto_Shape_zone;
            FoundDispSpatiallyHeterogInSettings[TVpars.size()-1]=true;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Geometric_Shape_zone")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().GeoG_zone;
            FoundDispSpatiallyHeterogInSettings[TVpars.size()-1]=true;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Sichel_Gamma_zone")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().Sichel_Gamma_zone;
            FoundDispSpatiallyHeterogInSettings[TVpars.size()-1]=true;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Sichel_Xi_zone")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().Sichel_Xi_zone;
            FoundDispSpatiallyHeterogInSettings[TVpars.size()-1]=true;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Sichel_Omega_zone")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().Sichel_Omega_zone;
            FoundDispSpatiallyHeterogInSettings[TVpars.size()-1]=true;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"ContinuousDemeSizeVariation")==0 || cmp_nocase_no__(var,"ContinuousDemeSizeChange")==0) {
			string locstring;
			stringstream strstr(buf.substr(pos+1));
			strstr>>locstring;
			if(cmp_nocase_no__(locstring,"Linear")==0) {
			    TVpars.back().ContDemeSizeChange="Linear";
			} else if(cmp_nocase_no__(locstring,"Exponential")==0) {
                TVpars.back().ContDemeSizeChange="Exponential";
            } else if(cmp_nocase_no__(locstring,"Logistic")==0) {
                TVpars.back().ContDemeSizeChange="Logistic";
            } else if(cmp_nocase(locstring,"False")==0 || cmp_nocase(locstring,"None")==0 || cmp_nocase(locstring,"F")==0) {
                TVpars.back().ContDemeSizeChange="None";
            } else {cerr<<"\nBad Argument " << locstring << "  for ContinuousDemeSizeVariation in " << SettingsFilename <<  " :"<<endl;
					  cerr<<"Only Linear, Exponential, Logistic or False/None/F are implemented. I exit."<<endl;
					  if (cinGetOnError) cin.get();
					  exit(-1);
				}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"ContinuousPopSizeVariation")==0 || cmp_nocase_no__(var,"ContinuousPopSizeChange")==0) {
			cerr<<"ContinuousPopSizeVariation, ContinuousPopSizeChange, ContinuousVariation and ContinuousChange are obsolete..."<<endl;
			cerr<<" Use ContinuousDemeSizeVariation or ContinuousLatticeSizeVariation, depending on what you want to simulate"<<endl;
			cerr<<" I exit."<<endl;			
			goto nextline;
		}
		if(cmp_nocase_no__(var,"ContinuousLatticeSizeVariation")==0 ||cmp_nocase_no__(var,"ContinuousLatticeSizeChange")==0) {
			string locstring;
			stringstream strstr(buf.substr(pos+1));
			strstr>>locstring;
			if(cmp_nocase_no__(locstring,"Linear")==0) {
			    TVpars.back().ContLatticeSizeChange="Linear";
			} else if(cmp_nocase_no__(locstring,"Exponential")==0) {
                TVpars.back().ContLatticeSizeChange="Exponential";
            } else if(cmp_nocase_no__(locstring,"Logistic")==0) {
                TVpars.back().ContLatticeSizeChange="Logistic";
            } else if(cmp_nocase(locstring,"False")==0 || cmp_nocase(locstring,"None")==0 || cmp_nocase(locstring,"F")==0 || cmp_nocase(locstring,"N")==0) {
                TVpars.back().ContLatticeSizeChange="None";
            } else {cerr<<"\nBad Argument " << locstring << "  for ContinuousLatticeSizeVariation in " << SettingsFilename <<  " :"<<endl;
				cerr<<"Only Linear, Exponential, Logistic or False/None/F are implemented. I exit."<<endl;
				if (cinGetOnError) cin.get();
				exit(-1);
			}
			goto nextline;
		}
		if(cmp_nocase_no__(var,"DensLogisticGrowthRate")==0 || cmp_nocase_no__(var,"LogisticGrowthRate")==0
		 || cmp_nocase_no__(var,"Dens_Logistic_Growth_Rate")==0 || cmp_nocase_no__(var,"Logistic_Growth_Rate")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().densLogisticGrowthRate;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"LatticeLogisticGrowthRate")==0 || cmp_nocase_no__(var,"Lattice_Logistic_Growth_Rate")==0) {
			stringstream strstr(buf.substr(pos+1));
			strstr>>TVpars.back().latticeLogisticGrowthRate;
			goto nextline;
		}
		if(cmp_nocase_no__(var,"Random_Translation")==0) {
			evaluateBool(random_translation,buf.substr(pos+1));
			goto nextline;
		}
        cerr<<"\n(!) Unknown keyword '"<<var<<"' in file "<<filename<<endl;
        if (cinGetOnError) {
            if (hacktempo && cinGetOnError) {
                cerr<<"\n Proceed to next keyword, for testing purposes ? (Yes = return / No-Abort = control-C) "<<endl;
                cin.get();
            }
			hacktempo=false;
        } else {
            cerr<<"I exit.\n"<<endl;
            exit(-1);
        }


nextline: ;	// the pleasure of sin :-)
    } while(!settings.eof());
settings.close();
return 0;
}

/***** verif options ***************************/

int apply_final_settings() {
static bool executed=false;
using namespace NS_translation;


    if (executed) {
       cerr<<"apply_final_settings() has been called twice... Program aborted. ";
       if (cinGetOnError) cin.get();
       exit(-1);
    }
    executed=true; // vraiement utile a mettre uniquement after som "new" pointer allocation in this function

    //	pauses so that the user takes notice of the typos in the settings file
    if (warningBool && cinGetOnError) {
    	cerr << "\nPress any key to resume... " << endl;
    	cin.get();
    }

//    ***********************************************************************
//    A FINAL CHECK OF LOCUS/MARKER SPECIFIC INPUT SETTINGS BEFORE SIMULATION
//    ***********************************************************************

    //	IF the total number of loci has not been specified but the mutational models have been specified
	if (locusNumVector.size() == 1 && mutModelIDVector.size() > 2) {
		locusNumVector.assign(mutModelIDVector.size(), 10);
		nMarker = (int) locusNumVector.size() - 1;
		n_locus = nMarker * 10;
	}
	//	ELSE IF the total number of loci not specified, then use the default value stored at the zeroth index
	else if (locusNumVector.size() == 1) {
		n_locus = locusNumVector[0];
		nMarker = 1;
		locusNumVector.push_back(locusNumVector[0]);
	}
	//	ELSE sum across all the specified loci numbers to obtain the total number of loci to be simulated.
	else {
		n_locus = 0;
		nMarker = (int) locusNumVector.size() - 1;
		for (int i = 1; i <= (int) locusNumVector.size() - 1; i++)
			n_locus += locusNumVector[i];
	}

	//	Expands unspecified/partly specified/marker-wise specified parameter values into a locus-wise specification
	if ((int) mutModelIDVector.size() != (n_locus + 1))
		expandParamVector(mutModelIDVector, "Mutation_Model");

	if ((int) mutRateVector.size() != (n_locus + 1))
		expandParamVector(mutRateVector, "Mutation_Rate");

	if ((int) minAlleleVector.size() != (n_locus + 1))
		expandParamVector(minAlleleVector, "Min_Allele_Number");

    if ((int) maxMutationsVector.size() != (n_locus + 1))
		expandParamVector(maxMutationsVector, "Max_Mutation_Number");

	if ((int) polyLociBoolVector.size() != (n_locus + 1))
		expandParamVector(polyLociBoolVector, "Polymorphic_Loci_Only");

	if ((int) minorAlleleFreqVector.size() != (n_locus + 1))
		expandParamVector(minorAlleleFreqVector, "Minor_Allele_Frequency");

	if ((int) varMutRateBoolVector.size() != (n_locus + 1))
		expandParamVector(varMutRateBoolVector, "Variable_Mutation_Rate");

	if ((int) kMinVector.size() != (n_locus + 1))
		expandParamVector(kMinVector, "Allelic_Lower_Bound");

	if ((int) kMaxVector.size() != (n_locus + 1))
		expandParamVector(kMaxVector, "Allelic_Upper_Bound");

	if ((int) kIniVector.size() != (n_locus + 1))
		expandParamVector(kIniVector, "Allelic_State_MRCA");

	if ((int) motifSizeVector.size() != (n_locus + 1))
		expandParamVector(motifSizeVector, "Repeated_motif_size");

	if ((int) pSMMVector.size() != (n_locus + 1))
		expandParamVector(pSMMVector, "SMM_Probability_In_TPM");

	if ((int) geomTPMVector.size() != (n_locus + 1))
		expandParamVector(geomTPMVector, "Geometric_Variance_In_TPM");

	if ((int) geomGSMVector.size() != (n_locus + 1))
		expandParamVector(geomGSMVector, "Geometric_Variance_In_GSM");

	for (int i = 0; i < n_locus; i++) {
		// tous les XXXVector[] comporte le modèle par default = KAM en position 0, donc i+1 correspond aux locus simulés,
        // == [0 -> n_locus] 0= default=KAM not used for simulation, 1->n_locus = markers simulated
		//	eliminating all "minor" incompatibility issues for sequence data and/or SNPs
		if (specifyLocusSeq) {
			//	checks if the specified Equilibrium_Frequencies are all positive and sum up to 1
			if (specifyVarBaseFreqs) {
				if ((varBaseFreqs[0] < 0.0) || (varBaseFreqs[1] < 0.0) || (varBaseFreqs[2] < 0.0) || (varBaseFreqs[3] < 0.0)) {
					cerr << "\nThe values of \"Equilibrium_Frequencies\" must be positive. \nPlease re-check the file "
						 << SettingsFilename << endl << "Aborting IBDSim...Press any key to exit.";
					if (cinGetOnError)
					   cin.get();
					exit(-1);
				}
				if ((varBaseFreqs[0] + varBaseFreqs[1] + varBaseFreqs[2] + varBaseFreqs[3]) != 1.0) {
					cerr << "\nThe values of \"Equilibrium_Frequencies\" must sum up to 1. \nPlease re-check the file "
						 << SettingsFilename << endl << "Aborting IBDSim...Press any key to exit.";
					if (cinGetOnError)
					   cin.get();
					exit(-1);
				}
			}
			//	checks if Sequence_Size and MRCA_Sequence correspond when they have been user-specified
			if (defMRCASequence.size() && ((int)defMRCASequence.size() != seqSize)) {
				cerr << "\nThe specified \"Sequence_Size\" and \"MRCA_Sequence\" do not correspond. \nPlease re-check the file "
					 << SettingsFilename << endl << "Aborting IBDSim...Press any key to exit.";
				if (cinGetOnError)
				   cin.get();
				exit(-1);
			}
			//	eliminating all "minor" incompatibility issues for sequence data
			if ((mutModelIDVector[i+1] == "ISM") || (mutModelIDVector[i+1] == "SNP") || (mutModelIDVector[i+1] == "JC69")
					|| (mutModelIDVector[i+1] == "K80") || (mutModelIDVector[i+1] == "K2P") || (mutModelIDVector[i+1] == "F81")
					|| (mutModelIDVector[i+1] == "HKY85") || (mutModelIDVector[i+1] == "TN93")) {
                if(mutModelIDVector[i+1] == "K2P") mutModelIDVector[i+1] = "K80";
				if (kMinVector[i+1] != -1)
					kMinVector[i+1] = -1;
				if (kMaxVector[i+1] != -1)
					kMaxVector[i+1] = -1;
				if (kIniVector[i+1] != -1)
					kIniVector[i+1] = -1;
				if (motifSizeVector[i+1] != -1)
					motifSizeVector[i+1] = -1;
				if (pSMMVector[i+1] != double(-1))
					pSMMVector[i+1] = -1;
				if (geomTPMVector[i+1] != double(-1))
					geomTPMVector[i+1] = -1;
				if (geomGSMVector[i+1] != double(-1))
					geomGSMVector[i+1] = -1;
			}
			//	eliminating all "minor" incompatibility issues for SNPs
			if (mutModelIDVector[i+1] == "SNP") {
				if (mutRateVector[i+1] != -1)
					mutRateVector[i+1] = -1;
				if (!varMutRateBoolVector[i+1])
					varMutRateBoolVector[i+1] = false;
				if (minAlleleVector[i+1] > 2) {
					cerr << "\n\"Min_Allele_Number\" for SNPs cannot be bigger than 2. \nPlease re-check the settings file "
						 << SettingsFilename << endl << "Aborting IBDSim...Press any key to exit.";
					if (cinGetOnError)
					   cin.get();
					exit(-1);
			    }
			}
			else if (minorAlleleFreqVector[i+1] != 0) {
				cerr << "\nThe \"Minor_Allele_Frequency\" option is only available for SNPs. \nPlease re-check the settings file "
					 << SettingsFilename << endl << "Aborting IBDSim...Press any key to exit.";
				if (cinGetOnError)
				   cin.get();
				exit(-1);
			}
		}
		//	ONLY for microsatellites
		//	checks if the specified allelic bounds contain the state of the MRCA
		if ((mutModelIDVector[i+1] == "KAM") && (mutModelIDVector[i+1] == "SMM")
				&& (mutModelIDVector[i+1] == "TPM") && (mutModelIDVector[i+1] == "GSM")
				&& (kIniVector[i+1]) && (kIniVector[i+1] > kMaxVector[i+1] || kIniVector[i+1] < kMinVector[i+1])) {
			cerr << "\n\"Allelic_State_MRCA\" not between \"Allelic_Upper_Bound\" and \"Allelic_Lower_Bound\" at locus number "
					<< i << ". \nPlease check the file " << SettingsFilename << endl
					<< "Aborting IBDSim...Press any key to exit.";
			if (cinGetOnError)
			   cin.get();
			exit(-1);
	    }
		//	checks for 0 <= minorAlleleFreqVector < 0.5
		if ((minorAlleleFreqVector[i+1] < 0) && (minorAlleleFreqVector[i+1] > 0.5)) {
			cerr << "\n\"Minor_Allele_Frequency\" needs to be between 0 and 0.5 for locus number "
					<< i << ". \nPlease check the file " << SettingsFilename << endl
					<< "Aborting IBDSim...Press any key to exit.";
			if (cinGetOnError)
			   cin.get();
			exit(-1);
		}
		//	overrides the minAlleleVector settings if polyLociBoolVector has been set or if
		//	minorAlleleFreqVector has been specified
		if ((minAlleleVector[i+1] < 2)
				&& (polyLociBoolVector[i+1] || (minorAlleleFreqVector[i+1] > 0)))
			minAlleleVector[i+1] = 2;
	}

    if(DG2002 && n_locus > 1) {
    	printf("\nCannot write multilocus files for DG2002 (=> Set n_locus=1)");
        cerr << "Check settings. I exit." << endl;
        if (cinGetOnError)
        	cin.get();
        exit(-1);
    }
    if(migraineSettingsbool && nMarker>1) {
            cout << endl;cout << endl;cout << endl;cout << endl;cout << endl;
            cout << "nMarker=" << nMarker << endl;
            cerr << "IBDSim is not ready for writing Migraine settingsfile with multiple markers." << endl;
            cerr << "The file will not be written but the simulations continue." << endl;
            if (cinGetOnError) cin.get();
        }


//    ***********************************************************************


    if (fichier_genepop.length()==0) fichier_genepop="gp";
    fichier_geneland_geno=fichier_genepop+"_gld";
    fichier_migrate=fichier_genepop+"_m";
	fichier_stepsim=fichier_genepop+"_dg";
	if (repet!=1) {
        fichier_geneland_geno+="_";
        fichier_genepop+="_";
    }

    if(TVpars[0].dimRes1<TVpars[0].dimRes2) {
        cerr<<"\n Found Lattice_SizeX<lattice_SizeY,";
        if (TVpars.size()>1) cerr<<" in first demographic phase.";
        cerr << endl;
        cerr << "X dimension should always be greater than Y dimension." << endl;
        cerr << "Swap X and Y lattice dimension in SettingsFile.\n I exit." << endl;
        if (cinGetOnError)
            cin.get();
        exit(-1);
        //cerr<<".\nThese two dimensions are swapped.";
        //if (cinGetOnError) cin.get();
        //int bidon=TVpars[0].dimRes1;
        //TVpars[0].dimRes1=TVpars[0].dimRes2;
        //TVpars[0].dimRes2=bidon;
    }

	for (vector<CTimeVaryingParams>::iterator it=TVpars.begin();it!=TVpars.end();it++)
	   if (it->dist_max<0) it->dist_max=max(it->dimRes1,it->dimRes2);


    if(cmp_nocase(boundaryType,"Circular")==0) {
                EdgeEffect="circ";
                SortieSupfnPtr1 = &sortie_sup_circ_abs1;
                SortieSupfnPtr2 = &sortie_sup_circ_abs2;
                SortieInffnPtr1 = &sortie_inf_circ_abs1;
                SortieInffnPtr2 = &sortie_inf_circ_abs2;
    } else if(cmp_nocase(boundaryType,"Absorbing")==0 || cmp_nocase(boundaryType,"Absorbing1")==0) {
                    EdgeEffect="abs";
                    SortieSupfnPtr1 = &sortie_sup_circ_abs1;
                    SortieSupfnPtr2 = &sortie_sup_circ_abs2;
                    SortieInffnPtr1 = &sortie_inf_circ_abs1;
                    SortieInffnPtr2 = &sortie_inf_circ_abs2;
    } else if(cmp_nocase(boundaryType,"Reflecting")==0) {
                        EdgeEffect="refl";
                        SortieSupfnPtr1 = &sortie_sup_refl1;
                        SortieSupfnPtr2 = &sortie_sup_refl2;
                        SortieInffnPtr1 = &sortie_inf_refl1;
                        SortieInffnPtr2 = &sortie_inf_refl2;
    } else {cout<<"\nBad Argument for Lattice_Boundaries in " << SettingsFilename <<  " :"<<endl;
					  cerr<<"Only Circular, Absorbing or Reflecting are possible options. I exit."<<endl;
				 	  if (cinGetOnError) getchar();
					  exit(-1);
    }

    if (Spec_Sample_Coord[0][x].size()!=Spec_Sample_Coord[0][y].size()) {
                cerr<<"SampleCoordinatesX and SampleCoordinatesY are of unequal size. I exit." <<endl;
			    if(cinGetOnError) cin.get();
				exit(-1);
    } else if (Spec_Sample_Coord[0][x].size()>0) {
        Specific_Sample_Designbool[0]=true;
        Spec_SampleSize[0]= (int) Spec_Sample_Coord[0][x].size();
    }
    if (Spec_Sample_Coord[1][x].size()!=Spec_Sample_Coord[1][y].size()) {
        cerr<<"Predisp_SampleCoordinatesX and Predisp_SampleCoordinatesY are of unequal size. I exit." <<endl;
        if(cinGetOnError) cin.get();
        exit(-1);
    } else if (Spec_Sample_Coord[1][x].size()>0) {
        Specific_Sample_Designbool[1]=true;
        predispbool=true;
        Spec_SampleSize[1]= (int) Spec_Sample_Coord[1][x].size();
    }
    if( ( (Specific_Sample_Designbool[0] && Spec_SampleSize[0]>1) && (!Specific_Sample_Designbool[1] && Spec_SampleSize[1]>1) ) ||
       ( (!Specific_Sample_Designbool[0] && Spec_SampleSize[0]>1) && (Specific_Sample_Designbool[1] && Spec_SampleSize[1]>1) ) ) {
        cerr<<"Specific sampling design must be set for both pre and postdisp samples, or not." << endl;
        cerr<<"Please modify the settings file." << endl;
        cerr<<"I exit." <<endl;
        if(cinGetOnError) cin.get();
        exit(-1);
    }
    if(Specific_Sample_Designbool[0] && dim_sample1[0]*dim_sample2[0]*dens_sample[0]>0) {
        cerr<<"A specific sample design is set for predisp" << endl;
        cerr<<"but classical rectangular sample dimension (predisp_sample_sizeX & Y) and non-null pre_dens_sample" << endl;
        cerr<<"are also specified in the settingfile=" << SettingsFilename << "." << endl;
        cerr<<"Please modify the setting file to remove this ambiguity." << endl;
        cerr<<"I exit." <<endl;
        if(cinGetOnError) cin.get();
        exit(-1);
    }
    if(Specific_Sample_Designbool[1] && dim_sample1[1]*dim_sample2[1]*dens_sample[1]>0) {
        cerr<<"A specific sample design is set for postdisp" << endl;
        cerr<<"but classical rectangular sample dimension (predisp_sample_sizeX & Y) and non-null pre_dens_sample" << endl;
        cerr<<"are also specified in the settingfile=" << SettingsFilename << "." << endl;
        cerr<<"Please modify the setting file to remove this ambiguity." << endl;
        cerr<<"I exit." <<endl;
        if(cinGetOnError) cin.get();
        exit(-1);
    } else if(!Specific_Sample_Designbool[1] && dim_sample1[1]*dim_sample2[1]*dens_sample[1]>0) {
        predispbool=true;
    }

    if( (Specific_Sample_Designbool[0] && ploidy*dens_sample[0]*Spec_SampleSize[0]<2)
       || (!Specific_Sample_Designbool[0] && ploidy*dens_sample[0]*dim_sample1[0]*dim_sample2[0]<2) ) {
        cerr<<"\n (!) Number of sampled genes (postdisp) is less than 2." << endl;
        cerr<<" (!) Please modify the setting file because IBDSim can not simulate samples with less than 2 genes." << endl;
        cerr<<"     I exit." <<endl;
        if(cinGetOnError) cin.get();
        exit(-1);
    }
    if(predispbool && ( (Specific_Sample_Designbool[1] && ploidy*dens_sample[1]*Spec_SampleSize[1]<2)
                       || (!Specific_Sample_Designbool[1] && ploidy*dens_sample[1]*dim_sample1[1]*dim_sample2[1]<2) ) ) {
        cerr<<"\n (!) Number of predisp sampled genes is less than 2." << endl;
        cerr<<" (!) Please modify the setting file because IBDSim can not simulate samples with less than 2 genes." << endl;
        cerr<<"     I exit." <<endl;
        if(cinGetOnError) cin.get();
        exit(-1);
    }

    if(cmp_nocase(TVpars.back().ContDemeSizeChange,"Linear")==0 || cmp_nocase(TVpars.back().ContDemeSizeChange,"Exponential")==0) {
                cerr<<"Requested Linear or Exponential change in pop size in the last demographic phase. " <<endl;
                cerr<<"These changes only make sense if a further demographic phase exists. I exit." <<endl;
				TVpars.back().PrintTVpars();
			    if(cinGetOnError) cin.get();
				exit(-1);
    }
	
	if(cmp_nocase(TVpars.back().ContLatticeSizeChange,"Linear")==0 || cmp_nocase(TVpars.back().ContLatticeSizeChange,"Exponential")==0) {
		cerr<<"Requested Linear or Exponential change in lattice size in the last demographic phase. " <<endl;
		cerr<<"These changes only make sense if a further demographic phase exists. I exit." <<endl;
		TVpars.back().PrintTVpars();
		if(cinGetOnError) cin.get();
		exit(-1);
    }

    int maxdimRes=0;
	for (vector<CTimeVaryingParams>::iterator it=TVpars.begin();it!=TVpars.end();it++) {
        maxdimRes=max(maxdimRes,it->dimRes1);
        maxdimRes=max(maxdimRes,it->dimRes2);
    }



    if( ! const_disp && (TVpars.size()==1 || maxdimRes==1)) {
        cerr<<"\nConstant_Dispersal is set to false but there is no change in time, or there is a single population:"<<endl;
        cerr<<"Constant dispersal is now set to true"<<endl;
        const_disp=true;
    }

    if(Kini!=0&&(Kini>Kmax||Kini<Kmin)) {
        cerr<<"\nProblem : Kini is not between Kmin and Kmax. Check these settings. I exit.";
        if (cinGetOnError) cin.get();
        exit(-1);
    }

    if(DG2002&&n_locus>1)
    {printf("\nCannot write multilocus files for DG2002 (=> Set n_locus=1)");
        cerr<<"Check settings. I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }
	if(TVpars[0].dimRes2==1) vide_sampleY[0]=1; else vide_sampleY[0]=vide_sample[0];
	if(TVpars[0].dimRes1==1) vide_sampleX[0]=1; else vide_sampleX[0]=vide_sample[0];
    if(predispbool){
        if(TVpars[0].dimRes2==1) vide_sampleY[1]=1; else vide_sampleY[1]=vide_sample[1];
        if(TVpars[0].dimRes1==1) vide_sampleX[1]=1; else vide_sampleX[1]=vide_sample[1];
    }
    if( ! Specific_Sample_Designbool[0]) {
        if (xmin_sample[0]==0) // no explicit value in settings
            xmin_sample[0]=TVpars[0].vide;
        if (ymin_sample[0]==0) // no explicit value in settings
            ymin_sample[0]=TVpars[0].vide;
        if(TVpars[0].dimRes1<(dim_sample1[0]-1)*vide_sampleX[0]+1) {
            cerr<<"Habitat dimension LatticeSizeX< sample dimension ="<<(dim_sample1[0]-1)*vide_sampleX[0]+1;
            cerr<<"\n    as implied by SampleSizeX="<<dim_sample1[0]<<" and void_sample_node="<<vide_sampleX[0];
            cerr<<"Check settings. I exit."<<endl;
            if (cinGetOnError) cin.get();
            exit(-1);
        }
        if(TVpars[0].dimRes2<(dim_sample2[0]-1)*vide_sampleY[0]+1) {
            cerr<<"Habitat dimension LatticeSizeY< sample dimension ="<<(dim_sample2[0]-1)*vide_sampleY[0]+1;
            cerr<<"\n    as implied by SampleSizeY="<<dim_sample2[0]<<" and void_sample_node="<<vide_sampleY[0];
            cerr<<"Check settings. I exit."<<endl;
            if (cinGetOnError) cin.get();
            exit(-1);
        }
    }
    if(predispbool && ! Specific_Sample_Designbool[1]) {
        if (xmin_sample[1]==0) // no explicit value in settings
            xmin_sample[1]=TVpars[0].vide;
        if (ymin_sample[1]==0) // no explicit value in settings
            ymin_sample[1]=TVpars[0].vide;
        if(TVpars[0].dimRes1<(dim_sample1[1]-1)*vide_sampleX[1]+1) {
            cerr<<"Habitat dimension LatticeSizeX< Predisp_sample dimension ="<<(dim_sample1[1]-1)*vide_sampleX[1]+1;
            cerr<<"\n    as implied by Predisp_Sample_SizeX="<<dim_sample1[1]<<" and Predisp_void_sample_node="<<vide_sampleX[1];
            cerr<<"Check settings. I exit."<<endl;
            if (cinGetOnError) cin.get();
            exit(-1);
        }
        if(TVpars[0].dimRes2<(dim_sample2[1]-1)*vide_sampleY[1]+1) {
            cerr<<"Habitat dimension LatticeSizeY< sample dimension ="<<(dim_sample2[1]-1)*vide_sampleY[1]+1;
            cerr<<"\n    as implied by SampleSizeY="<<dim_sample2[1]<<" and void_sample_node="<<vide_sampleY[1];
            cerr<<"Check settings. I exit."<<endl;
            if (cinGetOnError) cin.get();
            exit(-1);
        }
    }

    /** check that samples fit the population design when there is NO specific sampling  design:
    (i) test first coordinate
    (ii) test total size
    (iii) test steps
    **/

    if( ! Specific_Sample_Designbool[0] && !((xmin_sample[0] % TVpars[0].vide)==0)) {
        cerr<<"\nX Sample Coordinates are misspecified:"<<endl;
        cerr<<"The first X coordinate does not match the population design (Void_Nodes)"<<endl;
        cerr<<"Check settings. I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }

    if( ! Specific_Sample_Designbool[0] && !((ymin_sample[0] % TVpars[0].vide)==0)) {
        cerr<<"\nY Sample Coordinates are misspecified:"<<endl;
        cerr<<"The first Y coordinate does not match the population design (Void_Nodes)"<<endl;
        cerr<<"Check settings. I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }

    if(predispbool && ! Specific_Sample_Designbool[1] && !((xmin_sample[1] % TVpars[0].vide)==0)) {
        cerr<<"\nX Predisp_Sample Coordinates are misspecified:"<<endl;
        cerr<<"The first X coordinate does not match the population design (Void_Nodes)"<<endl;
        cerr<<"Check settings. I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }
    
    if(predispbool && ! Specific_Sample_Designbool[1] && !((ymin_sample[1] % TVpars[0].vide)==0)) {
        cerr<<"\nY Predisp_Sample Coordinates are misspecified:"<<endl;
        cerr<<"The first Y coordinate does not match the population design (Void_Nodes)"<<endl;
        cerr<<"Check settings. I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }


    if(!Specific_Sample_Designbool[0]) {
        if(TVpars[0].dimRes1<(xmin_sample[0]+(dim_sample1[0]-1)*vide_sampleX[0])){
            cerr<<"Habitat dimension LatticeSizeX< maximum sample coordinate ="<<xmin_sample[0]+(dim_sample1[0]-1)*vide_sampleX[0];
            cerr<<"\n    as implied by MinSampleCoordinateX="<<xmin_sample[0]<<" ,SampleSizeX="<<dim_sample1[0]<<" and void_sample_node="<<vide_sampleX[0];
            cerr<<"Check settings. I exit."<<endl;
            if (cinGetOnError) cin.get();
            exit(-1);
        }
        if(TVpars[0].dimRes2<(ymin_sample[0]+(dim_sample2[0]-1)*vide_sampleY[0])){
            cerr<<"Habitat dimension LatticeSizeY< maximum sample coordinate ="<<ymin_sample[0]+(dim_sample2[0]-1)*vide_sampleY[0];
            cerr<<"\n    as implied by MinSampleCoordinateY="<<ymin_sample[0]<<" ,SampleSizeY="<<dim_sample2[0]<<" and void_sample_node="<<vide_sampleY[0];
            cerr<<"Check settings. I exit."<<endl;
            if (cinGetOnError) cin.get();
            exit(-1);
        }
        if(xmin_sample[0]<=0 || ymin_sample[0]<=0){
            cerr<<"\nProblem with sample position on the lattice: X/Ymin_sample < 0";
            cerr<<"Check settings. I exit."<<endl;
            if (cinGetOnError) cin.get();
            exit(-1);
        }
    }
    if( ! Specific_Sample_Designbool[0] && !((vide_sample[0] % TVpars[0].vide)==0)) {
        cerr<<"\nY Sample Coordinates are misspecified:"<<endl;
        cerr<<"The spacing between samples (void_samples) is not a multiple of the spacing between occupied nodes (void_nodes)"<<endl;
        cerr<<"Check settings. I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }

    if(predispbool && !Specific_Sample_Designbool[1]) {
        if(TVpars[0].dimRes1<(xmin_sample[1]+(dim_sample1[1]-1)*vide_sampleX[1])){
            cerr<<"Habitat dimension LatticeSizeX < maximum Predisp_sample coordinate ="<<xmin_sample[1]+(dim_sample1[1]-1)*vide_sampleX[1];
            cerr<<"\n    as implied by Predisp_MinSampleCoordinateX="<<xmin_sample[1]<<" ,Predisp_SampleSizeX="
                <<dim_sample1[1]<<" and Predisp_void_sample_node="<<vide_sampleX[1];
            cerr<<"Check settings. I exit."<<endl;
            if (cinGetOnError) cin.get();
            exit(-1);
        }
        if(TVpars[0].dimRes2<(ymin_sample[1]+(dim_sample2[1]-1)*vide_sampleY[1])){
            cerr<<"Habitat dimension LatticeSizeY < maximum Predisp_sample coordinate ="<<ymin_sample[1]+(dim_sample2[1]-1)*vide_sampleY[1];
            cerr<<"\n    as implied by Predisp_MinSampleCoordinateY="<<ymin_sample[1]<<" ,Predisp_SampleSizeY="<<dim_sample2[1]
                <<" and Predisp_void_sample_node="<<vide_sampleY[1];
            cerr<<"Check settings. I exit."<<endl;
            if (cinGetOnError) cin.get();
            exit(-1);
        }
        if(xmin_sample[1]<=0 || ymin_sample[1]<=0){
            cerr<<"\nProblem with negative Predisp_sample position on the lattice: Predisp_xmin_sample="
            << xmin_sample[1] << "< 0; and Predisp_ymin_sample=" << ymin_sample[1] << "." << endl;
            cerr<<"Check settings. I exit."<<endl;
            if (cinGetOnError) cin.get();
            exit(-1);
        }
    }
    if(predispbool && ! Specific_Sample_Designbool[1] && !((vide_sample[1] % TVpars[0].vide)==0)) {
        cerr<<"\nY Predisp_Sample Coordinates are misspecified:"<<endl;
        cerr<<"The spacing between Predisp_samples (Predisp_void_samples) is not a multiple of the spacing between occupied nodes (void_nodes)"<<endl;
        cerr<<"Check settings. I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }

    //RL -> RL verfier les quatres tests suivants...
    if(Specific_Sample_Designbool[0] && !((Spec_Sample_Coord[0][x][0] % vide_sampleX[0])==0)) {
        cerr<<"\nX Sample Coordinates are misspecified:"<<endl;
        cerr<<"The first X coordinate (=" << Spec_Sample_Coord[0][x][0] << ") does not match void_sample_node=" << vide_sampleX[0] << ";" <<endl;
        cerr<<"Check settings. I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }
    if(Specific_Sample_Designbool[0] && !((Spec_Sample_Coord[0][y][0] % vide_sampleY[0])==0)) {
        cerr<<"\nY Sample Coordinates are misspecified:"<<endl;
        cerr<<"The first Y coordinate does not match void_sample_node"<<endl;
        cerr<<"Check settings. I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }
    if(predispbool && Specific_Sample_Designbool[1] && !((Spec_Sample_Coord[1][x][0] % vide_sampleX[1])==0)) {
        cerr<<"\nX Predisp_Sample Coordinates are misspecified:"<<endl;
        cerr<<"The first X coordinate (=" << Spec_Sample_Coord[1][x][0] << ") does not match Predisp_void_sample_node=" << vide_sampleX[1] << ";" <<endl;
        cerr<<"Check settings. I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }
    if(predispbool && Specific_Sample_Designbool[1] && !((Spec_Sample_Coord[1][y][0] % vide_sampleY[1])==0)) {
        cerr<<"\nY Predisp_Sample Coordinates are misspecified:"<<endl;
        cerr<<"The first Y coordinate does not match Predisp_void_sample_node"<<endl;
        cerr<<"Check settings. I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }

    /** END check that samples (post and predisp) fit the population design when there is NO specific sampling  design **/



    {vector<CTimeVaryingParams>::iterator it=TVpars.begin();
    int locit=1;
	while (it!=TVpars.end()) {
        
		if ((it->dimRes1 % it->vide) !=0 || (it->dimRes2 % it->vide) != 0) {
            cerr<<"\n In demographic Phase "<<locit<<" :";
            cerr<<"\n Any of lattice_size_X or lattice_size_Y is not a multiple or void_nodes. I exit";
            if (cinGetOnError) cin.get();
            exit(-1);
        }
		if(it->barrier || it->forwardBarrier) {
			if ((it->x1_barrier < 1) || (it->x2_barrier < 1) || (it->y1_barrier < 1) || (it->y2_barrier < 1)) {
				cerr<<"\n In demographic Phase "<<locit<<" :";
				cerr<<"\n At least one of the barrier coordinates is < 1. Check your settings. I exit."<< endl;
				if (cinGetOnError) cin.get();
				exit(-1);
			}
			if ((it->x1_barrier > it->dimRes1) || (it->x2_barrier > it->dimRes1)) {
				cerr<<"\n In demographic Phase "<<locit<<" :";
				cerr<<"\n At least one of the barrier X coordinates is larger than lattice size X. Check your settings. I exit."<< endl;
				if (cinGetOnError) cin.get();
				exit(-1);
			}
			if ((it->y1_barrier > it->dimRes2) || (it->y2_barrier > it->dimRes2)) {
				cerr<<"\n In demographic Phase "<<locit<<" :";
				cerr<<"\n At least one of the barrier Y coordinates is larger than lattice size Y. Check your settings. I exit."<< endl;
				if (cinGetOnError) cin.get();
				exit(-1);
			}
			if (it->barrierCrossingRate > 1.0) {
				cerr<<"\n In demographic Phase "<<locit<<" :";
				cerr<<"\n (!) barrierCrossingRate is > 1.0 (!) . Check your settings. I exit."<< endl;
				if (cinGetOnError) cin.get();
				exit(-1);
			}
			if ((it->y1_barrier != it->y2_barrier) && (it->x1_barrier != it->x2_barrier)) {
				cerr<<"\n In demographic Phase "<<locit<<" :";
				cerr<<"\n The barrier is nor a horizontal or vertical line. The current version of IBDSim can not consider such barriers. I exit";
				if (cinGetOnError) cin.get();
				exit(-1);
			}
			if ((it->y1_barrier == it->y2_barrier) && (it->x1_barrier > it->x2_barrier)) {
				cerr<<"\n In demographic Phase "<<locit<<" :";
				cerr<<"\n exchanging x1_barrier and x2_barrier values, because x1 > x2 in settings file " << SettingsFilename << endl;
				int temp=it->x2_barrier;
				it->x2_barrier=it->x1_barrier;
				it->x1_barrier=temp;
				if (cinGetOnError) {
					cerr << "Press any key to resume." << endl;
					cin.get();
				}
			}
			if ((it->x1_barrier == it->x2_barrier) && (it->y1_barrier > it->y2_barrier)) {
				cerr<<"\n In demographic Phase "<<locit<<" :";
				cerr<<"\n exchanging y1_barrier and y2_barrier values, because y1 > y2 in settings file " << SettingsFilename << endl;
				int temp=it->y2_barrier;
				it->y2_barrier=it->y1_barrier;
				it->y1_barrier=temp;
				if (cinGetOnError) {
					cerr << "Press any key to resume." << endl;
					cin.get();
				}
			}
			if (((it->x1_barrier == it->x2_barrier) && (it->y1_barrier == it->y2_barrier)) || (it->barrierCrossingRate > 0.99999)) {
				cerr<<"\n In demographic Phase "<<locit<<" :";
				cerr<<"\n Barrier settings will not have any effect because x1=x2 ";
                cerr << "and y1=y2  or barrierCrossingRate > 0.99999 in settings file " << SettingsFilename<< endl;
				it->barrier=false;
				if (cinGetOnError) {
					cerr << "Press any key to resume." << endl;
					cin.get();
				}
			}
			if ( (it->barrierCrossingRate < 0.0000000001) && ( ((it->x1_barrier == 1) && (it->x2_barrier == it->dimRes1)) || ((it->y1_barrier == 1) && (it->y2_barrier == it->dimRes2)) ) ) {
                int FoundIncompleteBarrier=false;
                vector<CTimeVaryingParams>::iterator itloc=(it+1);
                while (itloc!=TVpars.end()) {
                    if( (itloc->barrier || itloc->forwardBarrier) && ! (itloc->barrierCrossingRate < 0.0000000001) && ( ((itloc->x1_barrier == 1) && (itloc->x2_barrier == itloc->dimRes1)) || ((itloc->y1_barrier == 1) && (itloc->y2_barrier == itloc->dimRes2)) )) FoundIncompleteBarrier=true;
                        else if(!itloc->barrier && !itloc->forwardBarrier) FoundIncompleteBarrier=true;
                    itloc++;
                }
                if(!FoundIncompleteBarrier) {
                    cerr<<"\n In demographic Phase "<<locit<<" :";
                    cerr<<"\n (!) barrierCrossingRate is " << it->barrierCrossingRate << " (< 0.0000000001) with a complete vertical or horizontal barrier (!) ." << endl;
                    cerr << " with such configuration, either IBDSim will never stop a run because there is no possible MRCA if sample is taken across the barrier" << endl;
                    cerr << " or the barrier will only cut the lattice into a reacheable part and a non-reacheable one." << endl;
                    cerr << " We suspect it is not the desired configuration." << endl;
                    cerr << " Check your settings. I exit."<< endl;
                    if (cinGetOnError) cin.get();
                    exit(-1);
                }
			}
			if (it->x1_barrier == it->x2_barrier) it->barrierX=true;
			if (it->y1_barrier == it->y2_barrier) it->barrierY=true;
		}
		
        it++;locit++;
    }
    }//end bloc local var


    if((Specific_Sample_Designbool[0] || Specific_Sample_Designbool[1]) && (Prob_Id_Matrix || calculoui || HexpNeioui || Varoui || iterativeStats || seqStatsOui)) {
        cerr<<"\n\nSpecific sample design is not compatible with"<<endl;
        cerr<<"Hexp=true,Allelic_Variance=true, Pob_id_Matrix=true or Iterative_Statistics=true "<<endl;
        cerr<<"All computations of genetic statistics will be set to false...."<<endl;
        if (cinGetOnError) {
            cout<<"Press any key to continue computations...."<<endl;
            cin.get();
        }
        calculoui=false;
        HexpNeioui=false;
        Varoui=false;
        iterativeStats=false;
        Prob_Id_Matrix=false;
        seqStatsOui=false;
	//exit(-1);
    } else/* if(predispbool && (vide_sample[0]!=vide_sample[1] || dens_sample[0]!=dens_sample[1]) && (Prob_Id_Matrix || calculoui || HexpNeioui || Varoui || iterativeStats || seqStatsOui)) {
        cerr<<"\n\nConsidering different predisp and postdisp values for Void_Sample_Node or Ind_Per_Pop_Sampled is not compatible with"<<endl;
        cerr<<"Hexp=true,Allelic_Variance=true, Pob_id_Matrix=true or Iterative_Statistics=true "<<endl;
        cerr<<"All computations of genetic statistics will be set to false...."<<endl;
        if (cinGetOnError) {
            cout<<"Press any key to continue computations...."<<endl;
            cin.get();
        }
        calculoui=false;
        HexpNeioui=false;
        Varoui=false;
        iterativeStats=false;
        Prob_Id_Matrix=false;
        seqStatsOui=false;
        //exit(-1);
    } else*/ {
        if((HexpNeioui ||Varoui ||iterativeStats || suiviQ || Fisoui) && !calculoui) {
            calculoui=true;
        }
		//exit(-1);
	}
	if(!Simple1DProductDispBool) {//RL072014 was cmp_nocase(twoD_disp,"1DProductWithoutm0")==0) {
        vector<CTimeVaryingParams>::iterator it=TVpars.begin();
        int locit=1;
        while (it!=TVpars.end()) {
            if (it->mod != 'g' || it->mod != 'S') {
                cerr<<"\n twoD_disp==\"1DProductWithoutm0\", with incompatible dispersal model (not geometric or Sichel) in demographic Phase "<<locit<<" :";
                cerr<<"\n 1DProductWithoutm0 works only with geometric and Sichel dispersal. twoD_disp is set to Simple1DProduct";
                //twoD_disp="Simple1DProduct";
                Simple1DProductDispBool=true;
                if (cinGetOnError) cin.get();
            }
            it++;locit++;
        }
    }

    for (unsigned int it=0;it<TVpars.size();it++) {
        dim_reseau1=max(dim_reseau1,TVpars[it].dimRes1);
        dim_reseau2=max(dim_reseau2,TVpars[it].dimRes2);
    }
    grande_dim=max(dim_reseau1,dim_reseau2);

    
    for (unsigned int it=0;it<TVpars.size();it++) {
        if (( ! TVpars[it].zone) && (TVpars[it].vide==1) && (! TVpars[it].forwardBarrier)) { // pas de zone, ni de noeuds vide, ni de barrier forward
            if (it>0 || !Specific_Density_Designbool ) {
                TVpars[it].fixe=true;
            } else { //it==0 and Specific_Density_Designbool
                TVpars[it].fixe=false;
                GlobalDensSpatiallyHeterog=true;
                TVpars[it].DensSpatiallyHeterog=true;
            }
        } else {// ajout 1 RL 07 2014, il y a une zone ou des noeuds vide ou une barrier forward
            TVpars[it].fixe=false;
            if((TVpars[it].zone) && (TVpars[it].xmin_zone==-666 || TVpars[it].xmax_zone==-666
                                          || TVpars[it].ymin_zone==-666 || TVpars[it].ymax_zone==-666) ) {
                cerr << "A special zone is specified in the settings but has missing coordinates." << endl;
                cerr << "Check your settings file " << SettingsFilename << endl;
                cerr << " I exit" << endl;
                if(cinGetOnError) cin.get();
                exit(-1);
            }
            if(TVpars[it].zone && TVpars[it].vide_zone==-666 && TVpars[it].dens_zone==-666
                                && TVpars[it].dist_max_zone==-666 && TVpars[it].mod_zone=='n' && TVpars[it].Mig_zone==-666
                                && TVpars[it].GeoG_zone==-666 && TVpars[it].Pareto_Shape_zone==-666 && TVpars[it].Sichel_Gamma_zone==-666
                                &&   TVpars[it].Sichel_Xi_zone==-666 && TVpars[it].Sichel_Omega_zone==-666 ) {
                cerr << "A special zone is specified in the settings but no special parameters have been specified for the zone." << endl;
                cerr << "Check your settings file " << SettingsFilename << endl;
                cerr << " I exit" << endl;
                if(cinGetOnError) cin.get();
                exit(-1);
            }
            //heterogeneités spatiales de densités
            if( (TVpars[it].zone && ( (TVpars[it].dens_zone!=-666) || (TVpars[it].vide_zone!=-666) )
                 && ((TVpars[it].dens != TVpars[it].dens_zone) || (TVpars[it].vide != TVpars[it].vide_zone)) )
               || (TVpars[it].vide != 1) ) {
                GlobalDensSpatiallyHeterog=true;
                TVpars[it].DensSpatiallyHeterog=true;
                if(TVpars[it].dens_zone==-666) TVpars[it].dens_zone=TVpars[it].dens;
                if(TVpars[it].vide_zone==-666) TVpars[it].vide_zone=TVpars[it].vide;
            }
            //heterogeneités spatiales de dispersion
            if( TVpars[it].zone && FoundDispSpatiallyHeterogInSettings[it] ) {
                 GlobalDispSpatiallyHeterog=true;
                 TVpars[it].DispSpatiallyHeterog=true;
                if(TVpars[it].vide_zone==-666) TVpars[it].vide_zone=TVpars[it].vide;
                if(TVpars[it].dens_zone==-666) TVpars[it].dens_zone=TVpars[it].dens;
                if(TVpars[it].dist_max_zone==-666) TVpars[it].dist_max_zone=TVpars[it].dist_max;
                if(TVpars[it].mod_zone=='n') TVpars[it].mod_zone=TVpars[it].mod;
                if(TVpars[it].Mig_zone==-666) TVpars[it].Mig_zone=TVpars[it].Mig;
                if(TVpars[it].GeoG_zone==-666) TVpars[it].GeoG_zone=TVpars[it].GeoG;
                if(TVpars[it].Pareto_Shape_zone==-666) TVpars[it].Pareto_Shape_zone=TVpars[it].Pareto_Shape;
                if(TVpars[it].Sichel_Gamma_zone==-666) TVpars[it].Sichel_Gamma_zone=TVpars[it].Sichel_Gamma;
                if(TVpars[it].Sichel_Xi_zone==-666) TVpars[it].Sichel_Xi_zone=TVpars[it].Sichel_Xi;
                if(TVpars[it].Sichel_Omega_zone==-666) TVpars[it].Sichel_Omega_zone=TVpars[it].Sichel_Omega;
            }
        }
        if( (it+1<TVpars.size()) && ( (TVpars[it].mod != TVpars[it+1].mod) || (TVpars[it].dist_max != TVpars[it+1].dist_max) || (TVpars[it].Mig != TVpars[it+1].Mig)
            || (TVpars[it].GeoG != TVpars[it+1].GeoG) || (TVpars[it].Pareto_Shape != TVpars[it+1].Pareto_Shape)
            || (TVpars[it].Sichel_Gamma != TVpars[it+1].Sichel_Gamma) || (TVpars[it].Sichel_Xi != TVpars[it+1].Sichel_Xi)
            || (TVpars[it].Sichel_Omega != TVpars[it+1].Sichel_Omega) ) ) {
            //cout << "\nTo Be Removed... const_disp set to FALSE..." << endl;
            const_disp=false;
        }
        //TVpars[it].PrintTVpars();
    }

//    {int nonfixes=0;
//        for (unsigned int it=0;it<TVpars.size();it++) if ( ! TVpars[it].fixe) nonfixes++;
//        if (nonfixes>1) {
//            cerr << "IBDSim is not designed to take into account temporal changes with heterogeneous density configurations." << endl;
//            cerr << "I exit." << endl;
//            if (cinGetOnError) cin.get();
//            exit(-1);
//        }
//    }

    if(GlobalDensSpatiallyHeterog){ // a adapter pour temporal changes with spatial heterogeneities : J'ai abandoné sym... prends juste plus de mmoire dans certains cas (e.g. vide=2, vide=3)
        if((!Specific_Density_Designbool) && (!TVpars[0].zone) && (!TVpars[0].fixe)) {
            if(TVpars[0].vide==1) {
                cerr << "Problem with vide=" << TVpars[0].vide << " and fixe=" << TVpars[0].fixe << " and zone=" << TVpars[0].zone << "." << endl;
                cerr << "Contact The authors. I exit." << endl;
                if (cinGetOnError) cin.get();
                exit(-1);
            } else {
                if(GlobalDispSpatiallyHeterog) {
                    cerr << "IBDSim is not yet designed to take into account void nodes with spatially heterogeneous dispersal." << endl;
                    cerr << " I exit." << endl;
                    if (cinGetOnError) cin.get();
                    exit(-1);
                }
            }
        }
        for(unsigned int it=0;it<TVpars.size();it++) /*if(TVpars[it].sym==0)*/ TVpars[it].sym=max(TVpars[it].dimRes1,TVpars[it].dimRes2); /// essai en abandonnant sym car trop compliqué avec dispSpatialHetero....22052014 marche.
    }
    
    if(GlobalDispSpatiallyHeterog) for(unsigned int it=1;it<TVpars.size();it++) TVpars[it].sym=max(TVpars[it].dimRes1,TVpars[it].dimRes2);


       /** ancien code
    if(DensSpatiallyHeterog){
        sym=0;
        if((!Specific_Density_Designbool) && (!zone0) && (!fixe0)) sym=vide0;
        if((!zone1)&&(!fixe1)&&(vide1>sym)&&(Gn1<2147483646)) sym=vide1;
        if((!zone2)&&(!fixe2)&&(vide2>sym)&&(Gn2<2147483646)) sym=vide2;
        if((!zone3)&&(!fixe3)&&(vide3>sym)&&(Gn3<2147483646)) sym=vide3;
        if(sym==0) sym=max(dim_reseau1,dim_reseau2);
    }
        **/
    if(GlobalDispSpatiallyHeterog && !Simple1DProductDispBool) {//RL072014 was cmp_nocase(twoD_disp,"1DProductWithoutm0")==0) {
        cout << "option 1DproductWithoutm0 is not implemented with spatially heterogeneous dispersal" << endl;
        cout << "You must change your settings file " << SettingsFilename << endl;
        cout << "I exit" << endl;
        if(cinGetOnError) cin.get();
        exit(-1);
    }

//    {bool anyContDemeSizeChange=false;
//        for (unsigned int it=0;it<TVpars.size()-1;it++)
//            if( cmp_nocase(TVpars[it].ContDemeSizeChange,"None") !=0 ) anyContDemeSizeChange=true;
//        if((GlobalDensSpatiallyHeterog || Specific_Density_Designbool) && anyContDemeSizeChange) {
//            cerr << "IBDSim is not yet designed to take into account continuous deme size changes with spatially heterogeneous density." << endl;
//            cerr << "Check ContinuousDemeSizeVariation, Void_nodes and Specific_Density_Design specifications in settings";
//            cerr << "I exit." << endl;
//            if (cinGetOnError) cin.get();
//            exit(-1);
//        }
//    }

	{bool anyContLatticeSizeChange=false;
        for (unsigned int it=0;it<TVpars.size()-1;it++)
            if( cmp_nocase(TVpars[it].ContLatticeSizeChange,"None") !=0 ) anyContLatticeSizeChange=true;
        if((GlobalDensSpatiallyHeterog || Specific_Density_Designbool) && anyContLatticeSizeChange) {
            cerr << "IBDSim has not yet been checked for simulating continuous lattice size changes with spatially heterogeneous density." << endl;
            cerr << "Check ContinuousLatticeSizeVariation, zone and Specific_Density_Design specifications in settings";
            cerr << "I exit." << endl;
            if (cinGetOnError) cin.get();
            exit(-1);
        }
    }
	
    if(cinGetOnError || pauseGP) { // no error occurred, but the test being true means that the session is interactive
        cout << "\n\n All settings checked.\n\n" << endl;
        mysleep(2000);
        
        // but not cin.get, no exit since no error
        
    }

    effacer_ecran();

//cout<<"fin apply_final_settings";
return 0;
}
/*****************/

