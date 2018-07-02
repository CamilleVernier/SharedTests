/***************************************************************************
© R. Leblois 2005- for code collection
leblois@supagro.inra.fr

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

// 28012009 : many changes and bug corrected (allocation and free memory for cum2, unaccurrate computations of dy_min & max in SetForwardDispersalDistributions())
// bugs were problematic only for cases with void_node>1
// delete procedures in isoutils2.cpp were also change to delete[] after valgrind warnings, same modifications in main()
// bug fixed in allocation of the deffective[][] mattrix
// diverse new intiatilisations (problems found with valgrind)

// V1.0 -> V1.1 : 23032009
// incorporation of a specific density matrix giving the number of individuals in each lattice node
// information is read in a file named DensityMatrix.txt
// this is yet implemented only for the first generation before a potential change but not after

//V1.1->v1.2 : 24032009
//incorporation of continuous time change of deme size : linear, exponential and logistic
//16042009 : corrected bug : problem in free(tabdis) in new_disp() not conditionned on DensSpatiallyHeterog

//31032010 : bug for void_sample_nodes>2 with a linear lattice. vide_sample replaced by vide_sampleX and vide_sampleY

//03052010 : bug with vide_sampleX and Y that did not take the vide_sample value with Specific_Sample_design=TRUE, L.2072


// regarder tout les commentaires FR et RL et les appliquer


//#define STL_DEBUG
//#define DEBUG
#ifdef DEBUG
//int maxprev=0; source de bug possible (mais trËs improbable sur nombre descendants
#endif


#include <cstdio>
#include <iostream> // cout...
#include <fstream> // ofstream
#include <sstream> //stringstream
#include <cstdlib> //needed for vector too!!
#include <cstring>
#include <cmath>
#include <vector>
#include <set>
#include <ctime>
#include <limits>
#include <iomanip>
//#include <unistd.h> // required for sleep, not used anymore, replaced by mysleep()
#include <algorithm> // for sort

// FR->RL tout a nettoyer dans ce ficheir :-)

using namespace std;


#define GOTO //comment this line to have nicer outputs under Xcode for example.


//uncomment the following line if compiling with wxDev-C++ or MinGW
//#include <windows.h> //for mkdir(), system() but also for _gotoxy
#ifdef WIN32

       #ifdef DEVCPP //user-supplied -DDEVCPP
              #include <Windows.h>
              //#define sleep Sleep //replaced by mysleep()
//              #include <values.h>
/*              #ifndef DBL_MIN //RL j'ai ajoutÈ ca car ca ne marchait pas avec DEVC++ et sous unix... prb dans BesselK.c DBL_MAX undeaclared
              #define DBL_MIN 4.94065645841246544e-324
              #endif
              #ifndef DBL_MAX
              #define DBL_MAX 1.79769313486231470e+308
              #endif */
       #else
            #include <sys/stat.h> //mkdir
       #endif //if not DEVCPP
#else
      #include <sys/stat.h> //mkdir
#endif

#include "MersenneTwister.h"
#include "genepop_extract.h"
#include "isoutils2.h"
#include "charfunc.h"
#include "Bessel.h"
#include "BesselK.h"
#include "lattice_s.h"
#include "settings.h"
#include "TimeVarPars.h"


//#define NO_MODULES
#ifdef NO_MODULES
#include "genepop_extract.cpp"
#include "charfunc.cpp"
#include "BesselK.cpp"
#include "TimeVarPars.cpp"
#include "isoutils2.cpp"
#include "settings.cpp"
#endif //NO_MODULES


#define petit 10E-10
#define t1 1
#define haut 1
#define bas 2
//#define x 1 todo
//#define y 2 tofo
#define sq(t) ({    \
auto v = (t);  \
v*v;                \
})

const int x=1,y=2;
bool cinGetOnError=false; //to pause on cerr messages in batchDebug mode; overridden by explicit call to Batch mode
bool pauseGP=false;
char spname[]="simul_pars.txt"; // initialisation de longueur suffisante
ofstream simulpars(spname,ios::out);


int (*SortieSupfnPtr1)(int coord, int dim);
int (*SortieSupfnPtr2)(int coord, int dim, int trans);
int (*SortieInffnPtr1)(int coord, int dim);
int (*SortieInffnPtr2)(int coord, int dim, int trans);

//ofstream ffishet;
bool GlobalSimulparsWrittenBool=false;
string SettingsFilename="IbdSettings.txt";//fichiers d'entree avec valeurs des parametres
string SpecificDensityFilename="DensityMatrix.txt";//fichiers avec valeurs de densites a chaque point du reseau


/*CONSTANTES MARQUEURS*/
//	Defines vector datatypes for the specification multi-locus markers (see below)
//	and initializes their zeroth index (i.e. for multi-locus markers) with their respective default values
//all those vector data types are [0 -> n_locus] 0= default=KAM not used for simulation, 1->n_locus = markers simulated
// here it fills the default KAM model, and other defaults values for other mutation models
vector<int> locusNumVector(1,10);
vector<string> mutModelIDVector(1,"KAM");
vector<double> mutRateVector(1,0.0005);
vector<int> minAlleleVector(1,1);
vector<int> maxMutationsVector(1,2147483647);
vector<bool> polyLociBoolVector(1,false);
vector<double> minorAlleleFreqVector(1,0);
vector<bool> varMutRateBoolVector(1,false);
vector<int> kMinVector(1,1);
vector<int> kMaxVector(1,10);
vector<int> kIniVector(1,0);
vector<int> motifSizeVector(1,1);
vector<double> pSMMVector(1,0.8);
vector<double> geomTPMVector(1,10);
vector<double> geomGSMVector(1,0.36);

//	Variables indexing most of the vectors defined above
double mu, _mu;//0.0000005;
string  model_mut;/*'1'=SMM, '2'=IAM, '3'=KAM, '4'=TPM, '5'=GSM*/
int min_allele;/*1 ou 2 si on veut que loci polymorphes, ou plus...*/
int max_mutations;/* si on veut limiter les mutations, e.g. pour simuler des SNPs...*/
bool variable_mubool;
int Kmin, Kmax;/* bornes pour tailles alleliques*/
int	Kini;//si 0 automatique...sinon cette valeur
int MotifSize;
double pSMM;/*proportion de SMM sous TPM*/
double var_geom_TPM;/* variance de la loi geometrique du TPM*/
double var_geom_GSM;/* variance de la loi geometrique du GSM*/

int n_locus; // total number of loci simulated for each replicate sample /*+s pr fichiers Genepop*/
int nMarker; //	variable containing the total number of markers specified by the user
// vector containing the total number of UNIQUE genetic markers specified by the user at the zeroth index
vector<int> nUniqMarker(2,1);//[0 -> nUniqMarker[0]], 0 -> nb of unique markers, then 1 ... nUniqMarkers = value of the "first mut model" for each uniq marker
vector<string> dispMut; // vector containing the names of the UNIQUE genetic markers specified by the user [1 -> NUniqMarkers[0]], for display purpose only
//	nexusFormat takes the values: "F" (default), "Haplotypes_only" or ""Haplotypes_and_Individuals"
string nexusFormat = "F";

/*PARAMETERS FOR THE JC69, K80/K2P, F81, HKY85, TN93 substitution models*/
enum {A, G, C, T};
//	the nucleotide substitution rate matrix
double subsRateMatrix[4][4];
//	the MRCA sequence is random by default unless specified by the user
string defMRCASequence = "";
//	default size of the simulated sequence
int seqSize = 20;
//	the transition-transversion ratio
double ratio_TITV = 0.5;
//	the transition(A<->G) vs. transition(C<->T) ratio (only for TN93)
double ratio_TITI = 1.0;
//	the default fixed and variable equilibrium base frequencies in the respective order 'AGCT'
double fixedBaseFreqs[4] = {0.25, 0.25, 0.25, 0.25}, varBaseFreqs[4] = {0.2, 0.3, 0.2, 0.3};

bool  genepopoui=true;/*si 1 ecrit les fichier genepop*/
bool  genepopNoCoordbool=false;/*si T n'écrit pas les coordonnées dans les fichier genepop*/
bool  migFilebool=false;
bool  graFilebool=false;
bool  commonSSwInArNumAndDenombool=true;
bool  DG2002=false;//1... si 1 ecrit les fichiers DG2002 !! que pour 1 locus!!!!!!!!
bool  genelandoui=false;/*si true ecrit les fichier geneland*/
bool  migrate_oui=false,migrate_lettre_oui=false;//si true ecrit les fichiers migrate, si true ecrit des lettres pour les alleles (limitÈ a 36 alleles)
bool  dossier_nbre_allele=false;//1 si dossiers spÈciale par nbre d'allÈles pour DG2002
bool  AllStates=false;//pour affichier tous les Ètats allÈliques, meme sans Èchantillons pour DG2002KAM = equivalent de AllStates dans migraine
int   repet=10;/*600... nbre de repetitions*/
int   ploidy=1;// 1 pour haploides 2 pour diploides


int		Seeds=14071789;//seeds
vector< vector < vector<int> > > Spec_Sample_Coord;//stocke les coord de l'echantillon quand user defined
bool	calcul_anc=false;//lance proba_d'id_anc et pas new
bool	HexpNeioui=false;//si 1 calcul l'heterozygotie attendue de Nei (somme des produit des frequences)
bool	Varoui=false;//si 1 calcul la variance de taille allelique
bool	Fisoui=false;//si 1 calcul le Fis et l'Heterozygotie observée
bool	noHeader=false;//si 1 n'écrit pas les header dans les fichiers stats et probid
bool    arRegression=false; // si 1 fait la regression de ar=f(dist)
bool    erRegression=false; // si 1 fait la regression de er=f(dist)
bool    moranIRegression=false; // si 1 fait la regression de er=f(dist)
bool  	effective_disp=false;// si 1 ecrit fichier de distribution des distance de dispersion "efficaces"
bool	suiviQ=false;/*si 1 suit les Q a chaque repet pour graphe*/
bool	iterativeStats=false; // si 1 ecrit les fichiers Fis, Het, !!!!! calcul_oui necessaire !!!!!
bool  	calculoui=false; /*si 1 lance calculs proba d'id et temps de coa, rajoute environ 10% de temps de calcul*/
bool 	const_disp=true; // permet Èviter calcul forw ‡ chaque rep
//string  twoD_disp="Simple1DProduct";
bool    Simple1DProductDispBool=true;
bool	compareMath=true; /** si true, dispersion jusqu'a dx_max et non limit'e par taille du reseau (ajout RL necessaire pour tester = pour comparer avec mathematica).
                                important for display of moments of the dispersal distribution **/
bool  	Prob_Id_Matrix=false; // estimation des probas d'identitÈ ‡ chaque position ÈchantillonnÈe
bool	seqStatsOui = false; // si 1 calcul nbre de sites segrégés, pairwise differences


/*CONSTANTES DEMOGRAPHIQUES */

string   EdgeEffect;//types of edge effects : 0=circular,1=absorbing (disperse until a ancestor is found inside the lattice),2=reflecting,3=FRAbsorbing (disperse if first ancestor found inside the lattice, or do not disperse)
int   dim_reseau1=0;/* dimensions maximale (=toutes generations confondues) en x*/
int   dim_reseau2=0;/* 1... dimensions maximale en y*/

vector< int >   dim_sample1(2,1);/* axial nb noeuds echantillonnÈs*/
vector< int >   dim_sample2(2,1);/*1... """""""""""""""*/
vector< int >   Spec_SampleSize(2,-1);//sample size for specific design
vector< int >   xmin_sample(2,1),ymin_sample(2,1);//1... coord minimales de l'Èchantillon
                          //ATTENTION si vide sample!=1 ajuster x et y min
vector< int >	dens_sample(2,1);/*10... was 30 or 20*//* nbre d'ind echantillonne par noeud*/
vector< int >	vide_sample(2,1);/*was 11 or 2*/
vector< int >	vide_sampleX(2,1);/*was 11 or 2*/
						/*1 si pas de noeuds vides
						 2,3 si il existe des noeuds vide ( 1sur 2, 1 sur 3)
						 dans l'Èchantillon*/
vector< int >	vide_sampleY(2,1);

//pour échantillon predisp
//int   dim_sample1[1]=0;/* axial nb noeuds echantillonnÈs*/
//int   dim_sample2[1]=0;/*1... """""""""""""""*/
//int   Spec_SampleSize[1]=-1;//sample size for specific design
//int   xmin_sample[1]=0,ymin_sample[1]=0;//1... coord minimales de l'Èchantillon
////ATTENTION si vide sample!=1 ajuster x et y min
//int	dens_sample[1]=0;/*10... was 30 or 20*//* nbre d'ind echantillonne par noeud*/
//int	vide_sample[1]=1;/*was 11 or 2*/
//int	vide_sampleX[1];/*was 11 or 2*/
///*1 si pas de noeuds vides
// 2,3 si il existe des noeuds vide ( 1sur 2, 1 sur 3)
// dans l'Èchantillon*/
//int	vide_sampleY[1];

/*CARACTERISTIQUES DE DISPERSION + TAILLES ET NBRE DE POPS*/

// *9=>sig≤=4 sur 48 pas* *0=>sig≤=4 sur 15 pas* *1=>step stone m=0.666*  *2=>sig≤=1*
// *3=>sig≤=100* *4=>sig≤=1 1noeud sur 2*  *5=>sig≤=1 1noeud sur 3*/
// *6=> sig≤=20* *7=>sig≤=10*  *8=>sig≤=4 1noeud sur 3* *a=>step stone m=0.333*
// *b=> step stone m=0.01* * g=geometric mig=0.01* S=>Sichel mixture, inverse gamma mixture

bool 	random_translation=false;//1 si random, 0 si au milieu
vector<CTimeVaryingParams> TVpars(1);
CTimeVaryingParams* currTVparsPtr=NULL;

/*variables globales*/
FILE      *fmoyenne,*fsuivi,*fDG2002,*fdisp,*fdisp_moy/*,*ffishet*/;/* pour fichier genepop,moyennes,...*/

ofstream fgenepop, fmigrate, fgeneland1, fgeneland2;
bool txt_extensionbool=true,predispbool=false,Specific_Density_Designbool=false,
        groupAllSamplesbool=false,individualizeAllSamplesbool=false,migraineSettingsbool=false;
vector< bool > Specific_Sample_Designbool(2,false);

MTRand alea;
int n_genes_total=0;// nbre de genes total des différents echantillons,
vector< int > n_genes(2,0);/* nbre de genes des l'échantillon post et predisp*/

string fichier_geneland_geno,fichier_geneland_coord,fichier_migrate,fichier_genepop,fichier_stepsim;
string Genepopfile_extension="";
int locus,rep;

namespace NS_translation {
    int grande_dim;
    int translationx=0,translationy=0;//facteur de translation au moment de l'augmentation de surface
               //pour ne pas tjs mettre l'echantillon dans un coin de la pop
}

namespace NS_coal_tree {
    int nbre_noeud_restant;
#ifdef STL_DEBUG
    vector<vector<int> > coord_individus;
    vector<vector<int> >coord_noeud,coord_ori;
    vector<int>no_noeud_restant,aleat_noeud;
#else
    int 	  **coord_individus/*[1601][3]*/;//tableau des coord dans l'echantillon
    int 		**coord_noeud,**coord_ori;
    int 		*no_noeud_restant,*aleat_noeud;
#endif
}

namespace NS_diagnostic_tables {
    long long int        **effective/*[60][3]*/,**effective_moy/*[60][3]*/,cum_moyx,cum_moyy;//tab des distances de dispersion "efficaces" : compteur,moy sur repets
//long long MAx=9223372036854775807
    double     **deffective/*[60][3]*/;//tab des distances de dispersion "efficaces" : frÈquences
    long long int   **effectiveImmigration,**cumulEffectiveImmigration,**effectiveImmigration_moy,**cumulEffectiveImmigration_moy;
    //SS for each sample (i.e. post & pre disp), for each sample pair (pairs of ind or pop), and for each locus
    //cf definitions in Rousset(2000)JEB IBD between individuals
    vector<vector<vector<double> > > SSb,SSw,SSb_SumForEachIndOverOthers,Qb_pair,Qw_pair,Qw_ind,QiOverAllOtherInd,Qb_distClass;
    //SS moyennes for each sample (i.e. pre & post disp) and for each sample pair (pairs of ind or pop),
    vector<vector<double> > SSw_SumOverAllPairs,Qb_meanAllPairs,Qw_meanAllInd;
    //SS & Q moyennes, moy=Mean over loci, and moy_glob=mean over repetations (= over data sets)
    vector<vector<double> > Qb_pair_moy,Qw_pair_moy,Qw_ind_moy,QiOverAllOtherInd_moy,Qb_distClass_moy;
    vector<double> Qb_meanAllPairs_moy,Qw_meanAllInd_moy,Qb_meanAllPairs_moy_glob,Qw_meanAllInd_moy_glob;
    //Geographic distance for each sample (i.e. pre & post disp) and for each sample pair (pairs of ind or pop),
    vector<vector<double> > indGeoDist,popGeoDist;
    vector<double> maxDistSample;//max distance over all sample pairs
    //Fstatistics and results of the regressions
    vector<vector<double> > Fstat_ar_fromSS,Fstat_ar_fromQ,Fstat_er_fromQ,Fstat_moranI_fromQ_forEachIndPairs,Fstat_moranI_fromQ_forEachDistClass,regression_ar,regression_er,regression_moranI;
    double minDistReg=1E-9;// minimal distance to take into account for the regression of ar and er against distance/log

}

int nbre_noeud_existant;
int no_noeud_coalesce;
int n_mut,compte_mut_loc;//,compte_mut_glob;//RL 022017 not used
//double mut_moy_glob,mut_moy_glob_var;//RL 022017 not used
int allelesAtCurrentLocus,*allele_nmbr;
//	vector containing the allele ID's (node.etat_allele) per run of the locus loop
vector<int> alleleID(1,0);
//	pointer to the allele frequency spectrum for a given locus
double *alleleFreq;
int *allele_max,*allele_min;
float mean_allele_nmbr;
char model_mig,model_mig_zone;
bool GlobalDensSpatiallyHeterog=false;/*T/F si densite constante en tout point du reseau*/ //DensSpatiallyHeterog mis dans TVpars pour chaque phase demo
bool GlobalDispSpatiallyHeterog=false;/*T/F si dispersion constante en tout point du reseau*/ //DsipSpatiallyHeterog mis dans TVpars pour chaque phase demo
int dx_min,dx_max,dy_min,dy_max,dispxx,dispyy  // dy_min,dy_max FER
			  ,deja /*0/1 si tab migra deja rempli*/
			  ,deja1/*0/1 si cumul deja fait*/,deja3;
long double     *migra,*migra_zone,*cum_fixex,*cum_fixey;
int 	   **densite;
double	   **deffectiveImmigration;
double     moy_dispx,sig_dispx,kurto_dispx,skew_dispx;//stat sur ditrib disp efficace en x
double     moy_dispy,sig_dispy,kurto_dispy,skew_dispy;//stat sur ditrib disp efficace en y
double     moy_demidispx,sig_demidispx,kurto_demidispx,skew_demidispx;//stat sur demi ditrib disp efficace en x
double     moy_demidispy,sig_demidispy,kurto_demidispy,skew_demidispy;//stat sur demi ditrib disp efficace en y
double     denom;
long double     norm;
long double* Sicheltable;
long double* Sicheltable_zone;
unsigned long int  mrca_max,*mrca/*[501]*/;
double mrca_moy,mrca_moy_glob;
unsigned long int currentGeneration;
//int 	  **coord_noeud/*[1601][3]*/;//tableau des coord des noeuds pdt la construction de l'arbre
//int	  	  **coord_ori/*[1601][3]*/;//a virer?
int       coord_origine[4];
//int   	  *aleat_noeud/*[1601]*/;
//int   	  *no_noeud_restant/*[1601]*/;

int       pas;
//long int  temps_coa;
int  temps_coa;
double    temps_coa_moy,coa_moy_glob;
vector< vector< vector < vector < vector<double> > > > > QMatrix_Moy_Glob,QMatrix_Moy;//moyenne globale sur les repets, moyenne sur tous les locus pour une repet
vector< vector< vector< vector< vector< vector<double> > > > > > QMatrix;//par locus pour une repet

vector<vector<double> >backwardPhilo; //matrix of backward non-immigration rate // FR 0610 essential for controlling backward rates

//long int  compteur3;
int  compteur3;
vector<double> Qind_moy,Qind_moy_glob,Qr_mean_moy,Qr_mean_moy_glob,hetero_moy,fis_moy_glob,hetero_moy_glob,
                fis_moy,fis_moy2,fis_moy_denom,fis_moy_numer,fis_moy_denom_glob,fis_moy_numer_glob,
                HexpNei_moy, HexpNei_moy_glob,HexpNei_moy_glob_var, Var_moy,Var_moy_glob,Var_moy_glob_var,
                M_moy,M_moy_glob,M_moy_glob_var;
vector< vector<double> > Q1_moy,Q1_moy_glob,Qind2,hetero,fis,HexpNei,Var,M;
vector< vector< vector<double> > > Q1,freq;
// maximum_allele_number is needed for allocating and initializing freq
int maximum_allele_number;

//	specifies if the the locus is of sequence type or SNP
bool specifyLocusSeq = false, specifyVarBaseFreqs = false;

vector<int> numSegSites;	//	total number of segregating sites since the MRCA per locus [0 -> n_locus-1], not per sample because it is computed on the whole tree
vector< vector<double> > numPairMismatch;	//	average number of pairwise mismatches for all sequences in the sample [SampleNb][0 -> n_locus-1]
vector< vector<double> > empTITVRatio;	//	average empirical (i.e. pairwise) transitions upon transversions ratio [SampleNb][0 -> n_locus-1]

double	  Hexp, Vexp;//attendu theoriques pour l'heterozygotie et la variance de taille allelique !! que pour certain mod mut
int 	  comptLoc=0;/*compteur pour proportion de locus monomorphes*/
double	  mono;/*proportion de locus monomorphes, ou ayant moins de MinAllelNb allelles*/
double	  tooManyMut;/*proportion de locus  ayant trop de mutation*/

struct    node{vector<int> descendant;
//					long int generation;
					unsigned long int generation;
					int etat_allele;
					char etat_allele_lettre;//pour migrate
					int ancetre;
					string sequence;	// DNA sequence information (variable used only for sequence loci and SNP)
//					int nbre_descendant;
				  };
struct    node **noeud;//[1601][501]
//vector <vector<node> > noeud;

struct    disdis{
    //double cum2[100][100]; /*
#ifdef DEBUG
	vector<vector<double> > cum2;
#else
	double** cum2;
#endif
};
struct    disdis **tabdis;




/*****************************************************************************/
/**********************programme principal ***********************************/
int main(int argc, char *argv[]) {

	int noeudi,noeudj,noeudk,i_ini,sampleNb;

    clock_t start, end;
    double temps_ecoule = 0.0;

    using namespace NS_coal_tree;
    using namespace NS_translation;
    using namespace NS_diagnostic_tables;

    unsigned int phaseDemo=0;
    //pour échantillon predisp, on met tout a 0...
    dim_sample1[1]=0,dim_sample2[1]=0,Spec_SampleSize[1]=-1,xmin_sample[1]=0,ymin_sample[1]=0,dens_sample[1]=0,vide_sample[1]=1;
    //et aussi pour le postdisp
    dim_sample1[0]=0,dim_sample2[0]=0,Spec_SampleSize[0]=-1,xmin_sample[0]=0,ymin_sample[0]=0,dens_sample[0]=0,vide_sample[0]=1;

#ifdef GOTO
    effacer_ecran();
#endif

/*cout<<numeric_limits<int>::max()<<endl;
cout<<numeric_limits<long int>::max()<<endl;
cout<<numeric_limits<long long int>::max()<<endl;*/

    Spec_Sample_Coord.resize(2);
    for(int a=0;a<2;a++) {
        Spec_Sample_Coord[a].resize(3); // FR->RL "3" pcq x=1 et y=2 et non x=0,y=1 Faire gaffe si on remplace x et y un jour...
        Spec_Sample_Coord[a][x].resize(0); // defaults implying no specific sampling design
        Spec_Sample_Coord[a][y].resize(0);
    }
    if (argc>1) {
        for(int i=1;i<argc;i++){
            string buf(argv[i]);
            string::size_type pos=std::min(buf.find('='),std::min(buf.find('\t'),buf.length()));
            string var=buf.substr(0,pos).c_str();
            if(cmp_nocase_no__(var,"IBDSettings")==0 or cmp_nocase_no__(var,"Settings")==0
               or cmp_nocase_no__(var,"SettingsFile")==0 or cmp_nocase_no__(var,"SettingFile")==0) SettingsFilename=buf.substr(pos+1);
            else if(cmp_nocase_no__(var,"densityFile")==0 or cmp_nocase_no__(var,"densityMatrix")==0 or cmp_nocase_no__(var,"densityFileName")==0 or cmp_nocase_no__(var,"densityMatrixFileName")==0) SpecificDensityFilename=buf.substr(pos+1);
            else if(cmp_nocase_no__(var,"RandomSeeds")==0 or cmp_nocase_no__(var,"Random_Seeds")==0 or cmp_nocase_no__(var,"seeds")==0) Seeds=atoi(buf.substr(pos+1).c_str());
            else {
                cerr<<"\n(!) Unknown keyword '"<<var<<"' in commandline arguments."<<endl;
                cerr<<" Check settings. I exit."<<endl;
                exit(-1);
            }
        }
    } else {SettingsFilename="IbdSettings.txt";SpecificDensityFilename="DensityMatrix.txt";}


	/*
	 //	CBR: same as the code above without the "SettingsFile" command line prefix
	 if (argc>1) {
	 string buf(argv[1]);
	 SettingsFilename=buf;
	 }
	 else
	 SettingsFilename="IbdSettings.txt";
	 */

	//	Checks for the End-of-Line characters of the input file
	int rc = set_eof_check_EOLtype(SettingsFilename, false);
	//rc indicates file not found
	if (rc == -1) {
		cerr<<"... now exiting";
		if (cinGetOnError)
			cin.get();
		exit(-1);
	}

	//	Attempts to convert the input file EOL format to that of the UNIX environment
	set_UNIX_EOLtype(SettingsFilename);

	read_settings_file(SettingsFilename.c_str());
    apply_final_settings();
#ifdef GOTO
	effacer_ecran();
#endif

//	code for determining the number and the names of UNIQUE genetic markers specified by the user (used in affichage_ecran() and simulpars)
	dispMut.assign(2,mutModelIDVector[1]);
	for (int i = 2; i <= (int) mutModelIDVector.size() - 1; i++) {
		int j = 1;
		while ((j != i-1) && (mutModelIDVector[j] != mutModelIDVector[i]))
			j++;
		if (mutModelIDVector[j] != mutModelIDVector[i]) {
			nUniqMarker[0]++; // increase the nb of UniqMarkers
			nUniqMarker.push_back(i); // store the Nb corresponding to its first occurence in mutModelIDVector
			dispMut.push_back(mutModelIDVector[i]); // store its name for display purpose
		}
	}
	affichage_ecran();
    fflush(stdin);
    fflush(stdout);

    start=clock();
    alea.seed(Seeds);
    if(Specific_Sample_Designbool[0]) n_genes[0]=ploidy*dens_sample[0]*Spec_SampleSize[0];
    else n_genes[0]=ploidy*dens_sample[0]*dim_sample1[0]*dim_sample2[0];
    if(predispbool){
        if(Specific_Sample_Designbool[1]) n_genes[1]=ploidy*dens_sample[1]*Spec_SampleSize[1];
        else n_genes[1]=ploidy*dens_sample[1]*dim_sample1[1]*dim_sample2[1];
    }
    
    n_genes_total=n_genes[0]+n_genes[1];
    grande_dim=max(dim_reseau1,dim_reseau2);

/*-----Allocation vecteurs/pointeurs utilises pour moyennes tout le long des repetitions---------*/
    if( !(Specific_Sample_Designbool[0] || Specific_Sample_Designbool[1])) {
        Q1_moy_glob.resize(2);
        for(int i=0;i<2;i++) Q1_moy_glob[i].resize(/*max(TVpars[0].dimRes1, */max( max(dim_sample1[0],dim_sample2[0]), max(dim_sample1[1],dim_sample2[1])/*)*/),0.0);
        Qr_mean_moy_glob.resize(2,0.0);
        Qb_meanAllPairs_moy_glob.resize(2,0.0);
        Qw_meanAllInd_moy_glob.resize(2,0.0);
        Qind_moy_glob.resize(2,0.0);
    }
    HexpNei_moy_glob.resize(2,0.0);
    HexpNei_moy_glob_var.resize(2,0.0);
    hetero_moy_glob.resize(2,0.0);
    M_moy_glob.resize(2,0.0);
    M_moy_glob_var.resize(2,0.0);
    Var_moy_glob.resize(2,0.0);
    Var_moy_glob_var.resize(2,0.0);
    fis_moy_glob.resize(2,0.0);
    fis_moy_denom_glob.resize(2,0.0);
    fis_moy_numer_glob.resize(2,0.0);
    densite=imatrix(dim_reseau1+1,dim_reseau2+1);
    if(effective_disp) {
        effective_moy=llimatrix(2*grande_dim+2,4);//tab des occurences de dispersion pour distribution de disp efficace
        deffective=dmatrix(2*grande_dim+2,4);//tab des frÈquences pour distribution de disp efficace
        effectiveImmigration_moy=llimatrix(dim_reseau1+1,dim_reseau2+1);//tab des occurences d'immigration efficaces
        cumulEffectiveImmigration_moy=llimatrix(dim_reseau1+1,dim_reseau2+1);//tab des occurences d'immigration efficaces
        deffectiveImmigration=dmatrix(dim_reseau1+1,dim_reseau2+1);
    }
/*---------------------------------------------------------------------------*/
/*----dimensionnement des vecteurs utilisÈs tout le long des repetitions-----*/
    if(Prob_Id_Matrix) {
        int SampleNb;
        if(predispbool) SampleNb=2; else SampleNb=1;
        QMatrix_Moy_Glob.resize(SampleNb);
        for(int a=0;a<SampleNb;a++) if(!Specific_Sample_Designbool[a]) {
            QMatrix_Moy_Glob[a].resize(dim_sample1[a]);
            for(int i=0;i<dim_sample1[a];i++) {
                QMatrix_Moy_Glob[a][i].resize(dim_sample2[a]);
                for(int j=0;j<dim_sample2[a];j++) {
                    QMatrix_Moy_Glob[a][i][j].resize(dim_sample1[a]);
                    for(int i2=0;i2<dim_sample1[a];i2++)
                        QMatrix_Moy_Glob[a][i][j][i2].resize(dim_sample2[a],0.0);
                }
            }
        }
    }
/*---------------------------------------------------------------------------*/

    ini_moy_glob();//initialisation pour moyennes sur repets

/*---initialisation de structure **noeud---------*/
//noeud=(struct node **) malloc((unsigned) (2*n_genes+2)*sizeof(struct node*));
    noeud=new node*[2*n_genes_total+2];
//noeud.resize(2*n_genes+2);
    if (!noeud) nrerror("allocation failure 1 in structure node()");
    for (int i=0; i<=(2*n_genes_total+1); i++){//noeud[i].resize(n_locus);
//  noeud[i]=(struct node *) malloc((unsigned) (n_locus)*sizeof(struct node));
        noeud[i]=new node[n_locus];
        if (!noeud[i]) nrerror("allocation failure 2 in structure node()");
    }

    //allocation de tabdis maintenant dans AllocNotFixDisp(), desaoull'e dans set_new_disp() et a la fin de chaque locus, RL 052015, a enlever
/*---fin initialisation de **noeud-------*/

    
    // les tableaux cum seront redimensionnÈs par AllocNotFixDisp appel'e dans forw

    /*pour vérifs matrice de mutation GSM
    int delta,iter=10000000,size=1;
    double sum=0.0;

    vector< int > locusNumVector2(1,10);
    vector< int > deltaCountTable(kMaxVector[1]+1,0);
    vector<double> deltaFreqTable(kMaxVector[1]+1,0.0);

    var_geom_GSM = geomGSMVector[1];
    cout << "var_geom_GSM=" << var_geom_GSM << endl;
    cout << "kMaxVector[1]=" << kMaxVector[1] << endl;

    for(int i=0;i<iter;i++){
        ici:
        delta=loi_geom(var_geom_GSM);
        //cout << "delta=" << delta << endl;
        double aleat=alea();
        if(aleat<0.5)
            delta=size+delta;
        else delta=size-delta;
        //printf("etat apres mut:%d",noeud[noeud[i][locus].descendant[k]][locus].etat_allele);
        if(delta>kMaxVector[1] || delta<=0) goto ici;
            else deltaCountTable[delta]++;
    }
    for(int i=1;i<=kMaxVector[1];i++) cout << deltaCountTable[i] << " ";
    cout << endl;
    for(int i=1;i<=kMaxVector[1];i++) {
        deltaFreqTable[i]=(double) deltaCountTable[i]/iter;
        if(i<=kMaxVector[1])  sum+=deltaFreqTable[i];
    }
    for(int i=1;i<=kMaxVector[1];i++) cout << deltaFreqTable[i] << " ";
    cout << "sum=" << sum << endl;
    //for(int i=1;i<=kMaxVector[1];i++) deltaFreqTable[i]=deltaFreqTable[i]/sum;
    //for(int i=1;i<=kMaxVector[1];i++) cout << deltaFreqTable[i] << " ";
    //cout << endl;

    fflush(stdout);//marche pas bien
    fflush(stdin);

    cin.get();
    */
#ifdef GOTO
    _gotoxy(0,16);
#endif
    cout<<"Dataset:          Locus:";//"\n nbre_noeud_restant"<<ends;
    fflush(stdout);//marche pas bien
    fflush(stdin);

/** @@@@@@@@@@@@@@@@@@@@@ BOUCLE PRINCIPALE @@@@@@@@@@@@@@@@@ **/
    for(rep=1;rep<=repet;rep++) {/*debut boucle sur repetition*/
#ifdef GOTO
        _gotoxy(8,16);
#endif
        cout << rep << "   ";
        fflush(stdout);//marche pas bien
        fflush(stdin);

        if(predispbool) sampleNb=2;  else sampleNb=1; // RL 062018 a deplacer juste apres readsettingsfile car commun a toute les repets

        //Initializing vectors containing the info on segregating sites for all loci with zeros for each repetition/run
        numSegSites.assign(n_locus, 0);// pas par échantillon car calculé sur tout l'arbre au moment de l'ajout des mutations
        numPairMismatch.resize(sampleNb);
        empTITVRatio.resize(sampleNb);
        for(int a=0;a<sampleNb;a++) {
            numPairMismatch[a].assign(n_locus, 0);
            empTITVRatio[a].assign(n_locus, 0);
        }
/*------pointeurs utilises pdt une boucle de repetition---*/

        if(effective_disp) {
            effective=llimatrix(2*grande_dim+2,4);//tab pour distribution de dispersion efficace
            effectiveImmigration=llimatrix(dim_reseau1+1,dim_reseau2+1);//tab des taux d'immigration efficaces par demes
            cumulEffectiveImmigration=llimatrix(dim_reseau1+1,dim_reseau2+1);//tab des taux d'immigration efficaces par demes
        }
        allele_nmbr=ivector(n_locus);//tableau du nbre d'allele pour chaque locus
        allele_max=ivector(n_locus);//pour range de taille
        allele_min=ivector(n_locus);//pour range de taille
        HexpNei.resize(2);
        for(int i=0;i<2;i++) HexpNei[i].resize(n_locus,0.0);
        hetero.resize(2);
        for(int i=0;i<2;i++) hetero[i].resize(n_locus,0.0);
        M.resize(2);
        for(int i=0;i<2;i++) M[i].resize(n_locus,0.0);
        Var.resize(2);
        for(int i=0;i<2;i++) Var[i].resize(n_locus,0.0);
        fis.resize(2);
        for(int i=0;i<2;i++) fis[i].resize(n_locus,0.0);
        Qind2.resize(2);
        for(int i=0;i<2;i++) Qind2[i].resize(n_locus,0.0);
        if(!Specific_Sample_Designbool[0] || !Specific_Sample_Designbool[1])
            Q1.resize(2);
            for(int i=0;i<2;i++) {
                Q1[i].resize(n_locus);
                for(int j=0;j<n_locus;j++)
                    Q1[i][j].resize( /*max(TVpars[0].dimRes1,*/ max( max(dim_sample1[0],dim_sample2[0]), max(dim_sample1[1],dim_sample2[1]) ) /*)*/);
            }
        if( !(Specific_Sample_Designbool[0] || Specific_Sample_Designbool[1]) ) {
            Q1_moy.resize(2);
            for(int i=0;i<2;i++) Q1_moy[i].resize(/*max(TVpars[0].dimRes1,*/max(max(dim_sample1[0],dim_sample2[0]), max(dim_sample1[1],dim_sample2[1]) /*)*/ ),0.0);
            Qr_mean_moy.resize(2,0.0);
            Qind_moy.resize(2,0.0);
        }
        HexpNei_moy.resize(2,0.0);
        hetero_moy.resize(2,0.0);
        M_moy.resize(2,0.0);
        Var_moy.resize(2,0.0);
        fis_moy2.resize(2,0.0);
        fis_moy.resize(2,0.0);
        fis_moy_numer.resize(2,0.0);
        fis_moy_denom.resize(2,0.0);


        mrca=luivector(n_locus);/*age mrca de tous les individus pour chaque locus*/
#ifdef STL_DEBUG
        coord_individus.resize(2*n_genes_total+2);
        for (vector<int>::iterator ii=coord_individus.begin();ii<coord_individus.end();ii++) (*ii).resize(4);
#else
        coord_individus=imatrix(2*n_genes_total+2,4);/*coordonnees des individus de l'echantillon de depart*/
#endif
/*------------------------------------------*/

/*-----dimensionnement des vecteurs utilisÈ pdt une repetition-------*/
        if(Prob_Id_Matrix && !(Specific_Sample_Designbool[0] || Specific_Sample_Designbool[1])) {
            QMatrix.resize(0);
            int SampleNb;
            if(predispbool) SampleNb=2; else SampleNb=1;
            QMatrix.resize(SampleNb);
            for(int a=0;a<SampleNb;a++) {
                QMatrix[a].resize(n_locus);
                for(int l=0;l<n_locus;l++) {
                    QMatrix[a][l].resize(dim_sample1[a]);
                    for(int i=0;i<dim_sample1[a];i++) {
                        QMatrix[a][l][i].resize(dim_sample2[a]);
                        for(int j=0;j<dim_sample2[a];j++) {
                            QMatrix[a][l][i][j].resize(dim_sample1[a]);
                            for(int i2=0;i2<dim_sample1[a];i2++)
                                QMatrix[a][l][i][j][i2].resize(dim_sample2[a],0.0);
                        }
                    }
                }
            }
            QMatrix_Moy.resize(0);//pour reinitialiser a chaque fois
            QMatrix_Moy.resize(SampleNb);
            for(int a=0;a<SampleNb;a++) {
                QMatrix_Moy[a].resize(dim_sample1[a]);
                for(int i=0;i<dim_sample1[a];i++) {
                    QMatrix_Moy[a][i].resize(dim_sample2[a]);
                    for(int j=0;j<dim_sample2[a];j++) {
                        QMatrix_Moy[a][i][j].resize(dim_sample1[a]);
                        for(int i2=0;i2<dim_sample1[a];i2++)
                            QMatrix_Moy[a][i][j][i2].resize(dim_sample2[a],0.0);
                    }
                }
            }
        }
/*-------------------------------------------------------------------*/



        ini_struct_noeud();
/*for(i=0;i<n_locus;i++)
{for(j=1;j<=(2*n_genes_total);j++)
  {printf("1noeud: %d;\n generation: %ld;\n etat allele: %d;ancetre: %d;
			  nbre de descendants: %d;\n"
			  ,j, noeud[j][i].generation,noeud[j][i].etat_allele,
			  noeud[j][i].ancetre,noeud[j][i].nbre_descendant);
  getchar();
  }
}
clrscr();*/

        ini_moy();/*initialisations pour calculs moyennes sur locus*/

        if(rep==1 && ! GlobalSimulparsWrittenBool) {
            const std::string datestring=__DATE__;
            const std::string timestring=__TIME__;
            simulpars<<"IBDsim built on "+datestring+" at "+timestring+"."<<endl<< endl; //FR->FR here maybe a version number issue
            simulpars<<"Random Seeds : "<<Seeds<<endl;
            simulpars<<"nbr of simulated samples : " << repet << " dataset(s) x " << n_locus << " locus = " << repet*n_locus << endl;
            simulpars<<"generic data file name: "<<fichier_genepop<<endl;
            if (EdgeEffect=="circ") simulpars<<"Closed (circular or toroidal) lattice"<<endl;
            if (EdgeEffect=="abs") simulpars<<"Lattice with absorbing boundaries"<<endl;
            if (EdgeEffect=="refl") simulpars<<"Lattice with reflecting boundaries"<<endl;
// (unreachable dispersal probability mass reported on the whole dispersal distribution including P(dx=0)=(1-mig))
            simulpars<<"nbr of loci per sample file (n_locus): "<<n_locus<<endl;
            if(!Specific_Sample_Designbool[0]) simulpars<<"Postdisp Sample: nbr of sampled nodes (dim_sample1 X dim_sample2): "<<dim_sample1[0]<<" x "<<dim_sample2[0]<<endl;
               else simulpars<<"Postdisp Sample: nbr of sampled nodes: "<<Spec_SampleSize[0] << ";" << endl;
            if(predispbool) {
                if(!Specific_Sample_Designbool[1])
                    simulpars<<"Predisp Sample: nbr of sampled nodes (dim_sample1[1] X dim_sample2[1]): "
                        <<dim_sample1[1]<<" x "<<dim_sample2[1]<<endl;
                    else simulpars<<"Predisp Sample: nbr of sampled nodes: "<<Spec_SampleSize[1]<< ";" << endl;
            }
            if(vide_sample[0]!=1) simulpars<<"Postdisp sample: Unsampled nodes between sampled nodes (void_sample_node-1): "<<vide_sample[0]-1<<endl;
            if(vide_sample[1]!=1) simulpars<<"Predisp sample: Unsampled nodes between sampled nodes (predisp_void_sample_node-1): "
                                                <<vide_sample[1]-1<<endl;
            simulpars<<"Postdisp sample: nbr of sampled individuals per sampled node (IndPerPopSampled): "<<dens_sample[0]<<endl;
            simulpars<<"Postdisp sample: x/y coordinate of \"first\" (bottom left) sampled node (xmin_sample/ymin_sample): "<<xmin_sample[0] << "/" << ymin_sample[0]<<endl<<endl;
            if(predispbool) {
                simulpars<<"Predisp Sample: nbr of sampled individuals per sampled node (Predisp_Ind_Per_Pop_Sampled): "<<dens_sample[1]<<endl;
                simulpars<<"Predisp Sample: x/y coordinate of \"first\" (bottom left) sampled node (xmin_sample[1]/ymin_sample[1]): "
                            <<xmin_sample[1] << "/" << ymin_sample[1]<<endl<<endl;
            }
            simulpars<<"=============="<<endl;
            simulpars << "Mutation model(s) (if multiple markers/loci have been specified, only the first locus of each marker is mentioned below):"
            		  << endl << endl;
            for (int i = 1; i <= nUniqMarker[0]; i++) {
            	if (dispMut[i] == "KAM" && dispMut[i] == "TPM" && dispMut[i] == "GSM" && dispMut[i] == "SMM")
            		simulpars << dispMut[i] << " (" << kMaxVector[nUniqMarker[i]] - kMinVector[nUniqMarker[i]]
							  << " alleles)";
            	else
            		simulpars << dispMut[i];

                if(mutRateVector[nUniqMarker[i]] <= 0)
                	simulpars << " Mutation rate not specified for this marker" << endl;
                else {
                    if(varMutRateBoolVector[nUniqMarker[i]])
                    	simulpars << " Mutation rate (mu) is variable for this marker : distribution is gamma with mean "
    							  <<  mutRateVector[nUniqMarker[i]] << endl;
                    else
                    	simulpars << " Mutation rate (mu) is constant for this marker : "
                    			  << mutRateVector[nUniqMarker[i]] << endl;
                }
            }
            simulpars<<"=============="<<endl<<endl;
            simulpars<<"Demographic phase: 0" <<endl;
            simulpars<<"At G=0 (sampling time ==  present) : "<<endl;
            simulpars<<"lattice dimensions are (LatticeSizeX * LatticeSizeY): "<<TVpars[0].dimRes1<<" x "<<TVpars[0].dimRes2<<endl;
            if(TVpars[0].vide!=1)  simulpars<<"with "<< (TVpars[0].vide-1) << " nodes over " <<TVpars[0].vide
                << " being empty, so that real gene density is : "<< ((double) ploidy*TVpars[0].initialDens/TVpars[0].vide)<<endl;
            if(ploidy==2) simulpars<<"Diploid individuals : nbr of genes per non empty lattice node (2 * IndsPerPop): "<<2*TVpars[0].initialDens<<endl;
            else simulpars<<"Haploid individuals : nbr of genes per lattice node (1* IndsPerPop): "<<TVpars[0].initialDens<<endl;
            if(TVpars[0].zone && TVpars[0].DensSpatiallyHeterog) {
                simulpars<< "xmin_zone=" << TVpars[0].xmin_zone << " and xmax_zone=" << TVpars[0].xmax_zone << endl;
                simulpars<< "ymin_zone=" << TVpars[0].ymin_zone << " and ymax_zone=" << TVpars[0].ymax_zone << endl;
                simulpars<< "voidNode in Zone =" << TVpars[0].vide_zone << endl;
                simulpars<< "Density in Zone (nbr of genes) =" << ploidy*TVpars[0].dens_zone << endl;
            }
            if(TVpars[0].backwardBarrier || TVpars[0].forwardBarrier) {
                simulpars << endl;
                if(TVpars[0].forwardBarrier) simulpars<<"Presence of a real forward barrier to gene flow with coordinates :" <<endl;
                    else simulpars<<"Presence of an approximate (i.e. backward) barrier to gene flow with coordinates :" <<endl;
                simulpars << "x1_barrier=" << TVpars[0].x1_barrier << " and x2_barrier=" << TVpars[0].x2_barrier << endl;
                simulpars << "y1_barrier=" << TVpars[0].y1_barrier << " and y2_barrier=" << TVpars[0].y2_barrier << endl;
                if (TVpars[0].barrierCrossingRateUp == TVpars[0].barrierCrossingRateDown) simulpars << "and a BarrierCrossingRate =" << TVpars[0].barrierCrossingRateUp << endl;
                else if(TVpars[0].forwardBarrier) {
                    simulpars << "with forward BarrierCrossingRateUp =" << TVpars[0].barrierCrossingRateUp << endl;
                    simulpars << "and forward BarrierCrossingRateDown =" << TVpars[0].barrierCrossingRateDown << endl;
                } else if(TVpars[0].backwardBarrier) {
                    simulpars << "with backward BarrierCrossingRateUp =" << TVpars[0].barrierCrossingRateUp << endl;
                    simulpars << "and backward BarrierCrossingRateDown =" << TVpars[0].barrierCrossingRateDown << endl;
                }
                    simulpars << "Be carefull with forward and backward barrier settings" << endl;
                    simulpars << "and the interpretation of asymetrical BarrierCrossingRateUp and BarrierCrossingRateDown" << endl;
             }
            if(TVpars[0].ContDemeSizeChange!="None") {
                simulpars<< "ContDemeSizeChange=" << TVpars[0].ContDemeSizeChange << endl;
            }
            if(TVpars[0].ContLatticeSizeChange!="None") {
                simulpars<< "ContLatticeSizeChange=" << TVpars[0].ContLatticeSizeChange << endl;
            }

            flush(cout);
            mysleep(2000);//so that the user can check all messages

            simulpars<<endl<<endl;
        }

//cout<<"ola\n"<<flush;

        for(locus=0;locus<n_locus;)/*boucle sur locus*/ {


//			Initializing the multi-locus marker parameters with their respective vector datatypes
            _mu = mutRateVector[locus+1];
            model_mut = mutModelIDVector[locus+1];
            min_allele = minAlleleVector[locus+1];
            max_mutations = maxMutationsVector[locus+1];
            variable_mubool = varMutRateBoolVector[locus+1];
            Kmin = kMinVector[locus+1];
            Kmax = kMaxVector[locus+1];
            Kini = kIniVector[locus+1];
            MotifSize = motifSizeVector[locus+1];
            pSMM = pSMMVector[locus+1];
            var_geom_TPM = geomTPMVector[locus+1];
            var_geom_GSM = geomGSMVector[locus+1];

            /*----pointeurs utilisÈ pdt une boucle sur un locus---*/
            //cout<<"huma"<<flush;
            //coord_noeud.resize(2*n_genes_total+2);
            //for (vector<vector<int> >::iterator ii=coord_noeud.begin();ii<coord_noeud.end();ii++) ii->resize(4);
            coord_noeud=imatrix(2*n_genes_total+2,4);/*position geographique des noeuds*/
            //int **coord_noeud=new int*[2*n_genes_total+2];
            //for (int ii=0;ii<(2*n_genes_total+2);ii++) coord_noeud[ii]=new int[4];
            //cout<<"humb"<<flush;
            //coord_ori.resize(2*n_genes_total+2);
            //for (vector<vector<int> >::iterator ii=coord_ori.begin();ii<coord_ori.end();ii++) ii->resize(4);
            coord_ori=imatrix(2*n_genes_total+2,4);/* position des lignees ancestrales*/
            //cout<<"humc"<<flush;
            //aleat_noeud.resize(2*n_genes_total+2);
            aleat_noeud=ivector(2*n_genes_total+2);/*identifie "les deux genes d'un individu"*/
            //cout<<"humd"<<flush;
            //no_noeud_restant.resize(n_genes_total+2);
            no_noeud_restant=ivector(n_genes_total+2);/*numero noeuds pas encore coalesces*/
            //cout<<"hume"<<flush;
            /*-------------------------------------------------*/

            //cout<<"\nlocus"<<locus+1<<endl<<flush;
            //getchar();
            //gotoxy(17,17);
            /*printf("locus: %d ;",locus+1);*/
            if(n_locus<100 || (n_locus<1000 && fmod(1.*(locus+1),10.0)==0.0) || (fmod(1.*(locus+1),100.0)==0.0)) {
#ifdef GOTO
                _gotoxy(25,16);
#endif
                cout<<locus+1<<"   ";
                fflush(stdout);
                fflush(stdin);
            }
            currentGeneration=1;/*initialisation du nbre de generations*/
            phaseDemo=0;
            //ini_tableaux_spe_FR_migrate();
            if( Specific_Sample_Designbool[0] || Specific_Sample_Designbool[1]) ini_tableaux_specifique();
                else ini_tableaux();/*initialise tableaux utilises pendant boucle sur locus*/
            nbre_noeud_existant=n_genes_total;/*nbre de noeud total existant a Gn=0*/
            nbre_noeud_restant=n_genes_total;/*nbre de noeuds n'ayant pas coalesce a GN=0*/
            compte_mut_loc=0;

            currTVparsPtr=&(TVpars[0]);
            currTVparsPtr->phase=0;
//            cout << "\ncurrTVparsPtr before set_new_disp(0)" << endl;
//            currTVparsPtr->PrintTVpars();
			set_new_disp(0);
//            cout << "\ncurrTVparsPtr after set_new_disp(0)" << endl;
//            currTVparsPtr->PrintTVpars();
            {//int anterieur=n_genes_total; // pour afficher la rÈduction progressive de l'arbre ancestral (rep==1)
                while (nbre_noeud_restant>1){ //boucle ppale pour chaque locus
                    //if (rep<3 && nbre_noeud_restant<anterieur-99) {
                    //  _gotoxy(20,20);
                    //  cout<<nbre_noeud_restant<<"   "<<flush;

                    // cout<<" "<<nbre_noeud_restant<<" "<<n_mut<<" "<<compte_mut_glob<<flush;
                    // anterieur=nbre_noeud_restant;
                    //} //if
                    //cout<<"\nola1";

                /*creation arbre de coalescence*/
                    //if((fmod(currentGeneration,10))==0.0) {
                    //if(currentGeneration==Gn1) {
                    //	printf("\ngeneration: %ld; ",currentGeneration);
                    //	cout<<"dens="<<dens<<flush;
                    //	cout<<"noeuds restants="<<nbre_noeud_restant<<flush;
                    //	getchar();
                    //}
                    if((phaseDemo+1<TVpars.size()) && (currentGeneration==TVpars[phaseDemo+1].minGen)) {
                        phaseDemo++;
//                        cout << "\nBefore Changing demo phase" << endl;
//                        currTVparsPtr->PrintTVpars();
                        currTVparsPtr=&(TVpars[phaseDemo]);// attention certains parametres ne sont changé que dans reset_demo_params() appelé dans set_new_disp()
                        currTVparsPtr->phase=phaseDemo;
                        if(  ( ( currTVparsPtr->dimRes1*currTVparsPtr->dimRes2 ) == 1 ) &&
                           ( ( TVpars[phaseDemo-1].dimRes1*TVpars[phaseDemo-1].dimRes2 ) > 1) )
                            currTVparsPtr->realise_disp(n_genes);  //fait la dispersion des noeuds restants pour tous les mettre dans la pop finale, car quand il n'y a qu'une pop, la migration ne se fait plus dans main()...

//                        cout << "\nBefore set_new_disp()" << endl;
//                        currTVparsPtr->PrintTVpars();
                        set_new_disp(phaseDemo);//appelle SetForwardDispersalDistributions() = choix et calucls de la disp; + appelle reset_demo_params()
//                        cout << "\nAfter set_new_disp() for a nex demo phase" << endl;
//                        currTVparsPtr->PrintTVpars();
//                        cout << endl;
                    } else {
                        translationx=translationy=0;//que quand changement de taille
                    }
					if( ( currTVparsPtr->dimRes1*currTVparsPtr->dimRes2 ) > 1) currTVparsPtr->realise_disp(n_genes);  //fait la dispersion des noeuds restants

					if( (cmp_nocase(currTVparsPtr->ContDemeSizeChange,"None") != 0) || (cmp_nocase(currTVparsPtr->ContLatticeSizeChange,"None") != 0) )
                        currTVparsPtr->currentDensAndLatticeSize(currentGeneration); //Compute current density;
//					cout << "\nMain : Gn=" << currentGeneration << "; phaseDemo=" << phaseDemo << endl;
//					cout << "\tContDemeSizeChange=" << currTVparsPtr->ContDemeSizeChange << "; ContLatticeSizeChange=" << currTVparsPtr->ContLatticeSizeChange << endl;
//					cout << "\tinitialDens=" <<currTVparsPtr->initialDens << ";dens=" <<currTVparsPtr->dens << endl;
//					cout << "\tinitialDimRes1=" <<currTVparsPtr->initialDimRes1 << "; dimRes1=" <<currTVparsPtr->dimRes1 << endl;
//					cout << "\tinitialDimRes2=" <<currTVparsPtr->initialDimRes2 << "; dimRes2=" <<currTVparsPtr->dimRes2 << endl;
                    nbre_aleat_noeud();/*identifie les deux genes de chaque individus pour coalescence*/
                    i_ini=2;
                    boucle1:
                    //processus de coalescence
                    for(noeudi=i_ini;noeudi<=nbre_noeud_restant;noeudi++) {/*comparaison noeuds 2 a 2*/
                       if (nbre_noeud_restant>1) {
                          for(noeudj=1;noeudj<noeudi;noeudj++) {
                             if((coord_noeud[no_noeud_restant[noeudi]][x]==coord_noeud[no_noeud_restant[noeudj]][x])//verifie tout d'abord s'ils sont danns le meme deme
									&&(coord_noeud[no_noeud_restant[noeudi]][y]==coord_noeud[no_noeud_restant[noeudj]][y])
									&&(aleat_noeud[no_noeud_restant[noeudi]]==aleat_noeud[no_noeud_restant[noeudj]])//puis si ils ont le meme parent
                               ) {
                                  if (noeud[no_noeud_restant[noeudi]][locus].generation!=currentGeneration)/*coa non multiples*/ {
                                     no_noeud_coalesce=nbre_noeud_existant+1;
                                     nbre_noeud_existant++;
                                     aleat_noeud[no_noeud_coalesce]=aleat_noeud[no_noeud_restant[noeudj]];
                                     coord_noeud[no_noeud_coalesce][x]=coord_noeud[no_noeud_restant[noeudj]][x];
                                     coord_noeud[no_noeud_coalesce][y]=coord_noeud[no_noeud_restant[noeudj]][y];
                                     noeud[no_noeud_coalesce][locus].generation=currentGeneration;
                                     noeud[no_noeud_restant[noeudi]][locus].ancetre=no_noeud_coalesce;
                                     noeud[no_noeud_restant[noeudj]][locus].ancetre=no_noeud_coalesce;
                                     //noeud[no_noeud_coalesce][locus].nbre_descendant++;
    #ifdef DEBUG
                                    // if (noeud[no_noeud_coalesce][locus].nbre_descendant++>maxprev)
                                    //  	{maxprev=noeud[no_noeud_coalesce][locus].nbre_descendant++; cout<<"!!"<<maxprev<<"!!"<<flush;}
    #endif
                                     noeud[no_noeud_coalesce][locus].descendant.push_back(no_noeud_restant[noeudi]);
                                     noeud[no_noeud_coalesce][locus].descendant.push_back(no_noeud_restant[noeudj]);
                                     //noeud[no_noeud_coalesce][locus].descendant[noeud[no_noeud_coalesce][locus]
                                    //	.nbre_descendant]=no_noeud_restant[noeudi];
                                     //noeud[no_noeud_coalesce][locus].nbre_descendant++;
                                     //noeud[no_noeud_coalesce][locus].descendant[noeud[no_noeud_coalesce][locus]
                                    //	.nbre_descendant]=no_noeud_restant[noeudj];
                                     for(noeudk=noeudj;noeudk<nbre_noeud_restant;noeudk++)
                                        no_noeud_restant[noeudk]=no_noeud_restant[noeudk+1];/*enleve le noeud j des noeud restants*/
                                     no_noeud_restant[nbre_noeud_restant]=no_noeud_coalesce;
                                     for(noeudk=noeudi-1;noeudk<nbre_noeud_restant;noeudk++)
                                        no_noeud_restant[noeudk]=no_noeud_restant[noeudk+1];/*enleve le noeud i des noeud restants*/
                                     no_noeud_restant[nbre_noeud_restant]=0;
                                     nbre_noeud_restant--;
                                     i_ini=noeudi-1;
                                     goto boucle1;
                                  }/*fin coalescence non multiple*/
                                  else {/* pour coalescences multiples*/
                                     no_noeud_coalesce=no_noeud_restant[noeudi];
                                     //noeud[no_noeud_coalesce][locus].nbre_descendant++;
                                     noeud[no_noeud_coalesce][locus].descendant.push_back(no_noeud_restant[noeudj]);
                                     //noeud[no_noeud_coalesce][locus].descendant[noeud[no_noeud_coalesce][locus].nbre_descendant]=no_noeud_restant[noeudj];
                                     noeud[no_noeud_restant[noeudj]][locus].ancetre=no_noeud_coalesce;
                                     for(noeudk=noeudj;noeudk<=nbre_noeud_restant;noeudk++)
                                        no_noeud_restant[noeudk]=no_noeud_restant[noeudk+1];
                                     nbre_noeud_restant--;
                                     i_ini=noeudi;
                                     goto boucle1;
                                  }/*fin coalescence multiple*/
                             }/*fin if coalescence possible*/
                          }/*fin boucle sur noeudj*/
                       }/*fin if sur noeud restant*/
                    }/*fin boucle sur noeudi=fin comparaison noeuds 2 a 2*/
                    currentGeneration++;/*incremente le numero de la generation*/
                }/*fin while: creation arbre de coalescence*/
            } //bloc int anterieur
//cout<<"fini\n"<<flush;


//verif
//printf("\nn_genes_total= %d  ;",n_genes_total);
//for(atej=1;atej<=(2*n_genes_total);atej++){
//   printf("\nancetre de %d= %d",atej,noeud[atej][locus].ancetre);
//   getchar();
//   }

            if(variable_mubool) mu=gamma_mut(_mu/2); else mu=_mu;//taux de mutation variable ou fixe
            comptLoc+=1;
            if (model_mut.compare("SMM")==0) {
               mutate_under_SMM();
            }
            else       if (model_mut.compare("IAM")==0) mutate_under_IAM();
            else       if (model_mut.compare("KAM")==0) mutate_under_KAM();
            else       if (model_mut.compare("TPM")==0) mutate_under_TPM();
            else       if (model_mut.compare("GSM")==0) mutate_under_GSM();
            else       if (model_mut.compare("JC69")==0) sequence_mutation_model();
            else       if (model_mut.compare("K80")==0) sequence_mutation_model();
            else       if (model_mut.compare("F81")==0) sequence_mutation_model();
            else       if (model_mut.compare("HKY85")==0) sequence_mutation_model();
            else       if (model_mut.compare("TN93")==0) sequence_mutation_model();
            else       if (model_mut.compare("SNP")==0) mutate_as_SNP(nbre_noeud_existant);
            else       if (model_mut.compare("ISM")==0) mutate_under_ISM();

            if(compte_mut_loc>max_mutations){// alors on resimule un nouveau locus
                tooManyMut += 1;
                ini_struct_noeud2();
            } else if (allelesAtCurrentLocus >= min_allele) {
            	if (minorAlleleFreqVector[locus+1] > 0) {
            		if (checkMinorAlleleFreq()) {
                        mrca[locus] = currentGeneration;
                        mrca_moy += 1.0*currentGeneration/n_locus;
                        if(calculoui && ploidy == 2)
                        	calcul_coa_individuel();//calcul temps de coa intra-individus
                        if(mrca[locus] > mrca_max)
                        	mrca_max = mrca[locus];
                        locus++;
            		}
                    else {
                        mono += 1;
                        ini_struct_noeud2();
                    }
            	}
            	else {
                    mrca[locus] = currentGeneration;
                    mrca_moy += 1.0*currentGeneration/n_locus;
                    if(calculoui && ploidy == 2)
                    	calcul_coa_individuel();//calcul temps de coa intra-individus
                    if(mrca[locus] > mrca_max)
                    	mrca_max = mrca[locus];
                    locus++;
            	}
            }
            else {
                mono += 1;
                ini_struct_noeud2();
            }
            //cout << "mono=" << mono << "; tooManyMut=" << tooManyMut << "; comptLoc=" << comptLoc << endl;
            //cout << "mono/comptLoc=" << mono/comptLoc << "; tooManyMut/comptLoc=" << tooManyMut/comptLoc << endl;


/*------liberation des pointeurs utilisÈs pdt boucle sur locus-----*/
//cout<<"hum1"<<flush;
            free_imatrix(coord_noeud,2*n_genes_total+2);
//for (int ii=0;ii<(2*n_genes_total+2);ii++) delete(coord_noeud[ii]);delete(coord_noeud);

//cout<<"hum2"<<flush;
            free_imatrix(coord_ori,2*n_genes_total+2);
//cout<<"hum3"<<flush;
            free_ivector(aleat_noeud);
//cout<<"hum4"<<flush;
            free_ivector(no_noeud_restant);
//cout<<"hum5"<<flush;
            free_ldvector(migra);
            free_ldvector(cum_fixex);
            free_ldvector(cum_fixey);
            if(GlobalDispSpatiallyHeterog) free_ldvector(migra_zone);


            if( (GlobalDensSpatiallyHeterog || GlobalDispSpatiallyHeterog || TVpars[phaseDemo].forwardBarrier) && !TVpars[phaseDemo].fixe) {
                //cum2 alloue dans disdis qui n'est appele qu'une fois a la premiere repet, au premier locuset a la premiere generation ou lors de changements demo si necessaire
                for(int i=TVpars[phaseDemo].dimRes1+1; i>=0; i--) {
                    for(int j=TVpars[phaseDemo].dimRes2+1; j>=0; j--) {
#ifdef DEBUG
#else
                        //cout << "free memory for tabids[" << i << "][" << j << "]." << endl;

                        for(int k=2*dx_max+1; k>=0; k--) {
                            //cout << "and cum2[" << k << "]" << endl;
                            delete[] (tabdis[i][j].cum2[k]);
                        }
                        //cout << "free cum2" << endl;
                        //getchar();
                        delete[] (tabdis[i][j].cum2);
#endif
                    }
                    //tabdis est maintenant allou'e pour chaque locus si il y a des heterogeneites dans AllocNotFixDisp(),
                    //il est desaoull'e dans set_new_disp() que lors d'une nouvelle phase, il faut donc le desallouer ici aussi
                    delete[] tabdis[i];
                }
                delete[] tabdis;

            }

/*-----------------------------------------------*/

        }/**fin boucle locus**/
//cout<<"fin boucle locus\n"<<flush;

//      maximum_allele_number identifies the maximum number of alleles found across all loci and for each run
        maximum_allele_number = 0;
    	for (int i = 0; i < n_locus; i++) {
    		maximum_allele_number = max(maximum_allele_number, allele_max[i]);
    	}
//      freq is allocated here as maximum_allele_number can be determined only at the end of the locus loop
//    	and freq is potentially needed when calcul_VariousStats() is called
        freq.resize(2);
        for(int i=0;i<2;i++) {
            freq[i].resize(n_locus);
            for(int j=0;j<n_locus;j++)
                freq[i][j].resize(maximum_allele_number+1,0.0);
        }

//file_tester(2);

        if(genepopoui) ecriture_fichier_genepop();
        if(genelandoui) {
            ecriture_fichier_genotypes();
            ecriture_fichier_coordonnees();
        }
        if(migrate_oui) ecriture_fichier_migrate();
        if(DG2002) {
            if(model_mut.compare("IAM")==0 || model_mut.compare("KAM")==0) ecriture_fichiers_DG2002KAM(); else ecriture_fichiers_DG2002SMM();
        }

//calcul la distribution de dispersion "efficace"
        if(effective_disp) {calcul_stat_disp();ecriture_disp();}

//fait divers calculs
        if(rep==repet) {
            mono=mono/comptLoc;
            tooManyMut=tooManyMut/comptLoc;
        }
        // Computation of various Proba identity, and regressions genet against distance/log
        if(arRegression || erRegression || moranIRegression) {
            computeSSForFstats();
            computeGeoDist();
            calcul_proba_identite_par_paires_individus();
            computeFstats_ar_er_moranI(commonSSwInArNumAndDenombool);// true if for bool commonSSwInNumAndDenom
            if(arRegression) {
                ecriture_fichier_MIG("ar");// true is for individual data, false for pops, RL 062018 remplacer par la stat = ar, er, Fst, FST/(1-Fst),C
                regression_ar.resize(sampleNb);
                for(int a=0;a<sampleNb;a++) regression_ar[a].resize(2,0.0);
                regression_ar=indRegression_ar_er_moranI(minDistReg, Fstat_ar_fromSS,"ar");
            }
            if(erRegression) {
                ecriture_fichier_MIG("er");// true is for individual data, false for pops, RL 062018 remplacer par la stat = ar, er, Fst, FST/(1-Fst),C
                regression_er.resize(sampleNb);
                for(int a=0;a<sampleNb;a++) regression_er[a].resize(2,0.0);
                regression_er=indRegression_ar_er_moranI(minDistReg, Fstat_er_fromQ,"er");
            }
            if(moranIRegression) {
                ecriture_fichier_MIG("moranI");// true is for individual data, false for pops, RL 062018 remplacer par la stat = ar, er, Fst, FST/(1-Fst),C
                regression_moranI.resize(sampleNb);
                for(int a=0;a<sampleNb;a++) regression_moranI[a].resize(2,0.0);
                regression_moranI=indRegression_ar_er_moranI(minDistReg, Fstat_moranI_fromQ_forEachIndPairs,"moranI");
            }
        }
        
        if(calculoui) {
            if( !(Specific_Sample_Designbool[0] || Specific_Sample_Designbool[1]) ) calcul_proba_identite_pour_distances_axiales();/*file_tester(21);*/
                else if(rep==repet) cout << "\n\n\n\n!! Calcul_Proba_Id() can not be run with a specific sample design, and will be ignored." << endl;
        //file_tester(3);
         if(!Specific_Sample_Designbool[0]  && !Specific_Sample_Designbool[1]) {
                //file_tester(4);
                calcul_VariousStats();
                if(iterativeStats) ecriture_fichier_iterativeStats();
                //ecriture_fichier_iterativeStats_cpp(); //!!!!RL : essai FR inutilisé
            } else if(rep==repet) cout << "\n!! Iterative_Statistics, Hexp and Fis can not be run with a specific sample design, and will be ignored." << endl;
            mrca_moy_glob+=mrca_moy/repet;
            ecriture_fichier_moyennes();
        }//en if(calculoui)
        
//calcul de la matrice des probabilitÈ d'identitÈ pour toutes les paires de coordonnÈes ÈchantillonnÈes
//indÈpendant des calculs de proba d'identitÈ par pas.
        if (Prob_Id_Matrix) {
            if( !(Specific_Sample_Designbool[0] || Specific_Sample_Designbool[1]) ) calcul_proba_identite_matrice();
             else {
                if(rep==repet) cout << "\n!! Prob_Id_Matrix can not be run with a specific sample design, and will be ignored." << endl;
             }
        //cout << ":apres calcul_proba_identitÈ_matrix()" << flush;
        }


/*-----liberation des pointeurs utilises pdt une boucle de repetition---*/
//free_livectorl(mrca);
        free_luivector(mrca);
#ifdef STL_DEBUG
#else
        free_imatrix(coord_individus,2*n_genes_total+2);
#endif
/*
        if( !(Specific_Sample_Designbool[0] || Specific_Sample_Designbool[1]) ) {
            free_dtab3(Q1,2,n_locus);
            free_dmatrix(Q1_moy,2);
            free_dmatrix(Qind2,2);
            free_dvector(Qr_mean_moy);
            free_dvector(Qind_moy);
        }
        free_dvector(HexpNei_moy);
        free_dvector(hetero_moy);
        free_dvector(Var_moy);
        free_dvector(M_moy);
        free_dvector(fis_moy2);
        free_dvector(fis_moy);
        free_dvector(fis_moy_numer);
        free_dvector(fis_moy_denom);
        free_dmatrix(HexpNei,2);
        free_dmatrix(hetero,2);
        free_dmatrix(Var,2);
        free_dmatrix(M,2);
        free_dmatrix(fis,2);
        free_dtab3(freq,2,n_locus);*/
        free_ivector(allele_max);
        free_ivector(allele_min);
        free_ivector(allele_nmbr);
        if(effective_disp) {
            free_llimatrix(effective,2*grande_dim+2);
            free_llimatrix(effectiveImmigration,dim_reseau1+1);
            free_llimatrix(cumulEffectiveImmigration,dim_reseau1+1);
        }
    /*------------------------------------------------------*/
    }/** @@@@@ fin boucle principale repetition @@@@@ **/

    if(effective_disp) {calcul_stat_disp_moy();ecriture_disp_moy();}
    if(Prob_Id_Matrix) if(!Specific_Sample_Designbool[0]  && !Specific_Sample_Designbool[1]) ecriture_matrice_proba_id();

/*-----liberation pointeurs utilises sur toutes les repetitions multilocus----------------------*/
/*    if( !(Specific_Sample_Designbool[0] || Specific_Sample_Designbool[1]) ) {
        free_dmatrix(Q1_moy_glob,2);
        free_dvector(Qr_mean_moy_glob);
        free_dvector(Qind_moy_glob);
    }
    free_dvector(HexpNei_moy_glob);
    free_dvector(HexpNei_moy_glob_var);
    free_dvector(hetero_moy_glob);
    free_dvector(M_moy_glob);
    free_dvector(M_moy_glob_var);
    free_dvector(Var_moy_glob);
    free_dvector(Var_moy_glob_var);
    free_dvector(fis_moy_glob);
    free_dvector(fis_moy_denom_glob);
    free_dvector(fis_moy_numer_glob);
*/
    free_imatrix(densite,dim_reseau1+1);
    if(effective_disp) {
        free_llimatrix(effective_moy,2*grande_dim+2);
        free_dmatrix(deffective,2*grande_dim+2);//crÈÈ dans calcul stat disp moy
        free_llimatrix(effectiveImmigration_moy,dim_reseau1+1);
        free_llimatrix(cumulEffectiveImmigration_moy,dim_reseau1+1);
        free_dmatrix(deffectiveImmigration, dim_reseau1+1);
    }
//for(i=(2*n_genes_total+1); i>=0; i--) free((noeud[i]));
//free((noeud));
    for(int i=(2*n_genes_total+1); i>=0; i--) delete[] noeud[i];
    delete[] noeud;
/*-----------------------------------------------------------------*/


    simulpars.close();
    
    if(min_allele>1) printf("\n\nProportion of discarded loci with less than %d alleles = %f",min_allele,mono);
    if(max_mutations<2147483647) printf("\n\nProportion of discarded loci with tooManyMut = %f",tooManyMut);

    end=clock();
    temps_ecoule= (double)(end - start) /CLOCKS_PER_SEC;
	printf("\n\n\n\nNormal ending of IBDSim (Total computation time is %f seconds)\n\n", temps_ecoule);

	if (pauseGP || cinGetOnError) {
		printf("\n...Press any key to stop the program and close the window...\n");
        cin.get();
    }
fflush(stdout);
fflush(stdin);
return 0;
}/***  fin programme principal main() ***/

/************************fin programme principal**********************/
/*********************************************************************/


//circle/torus = no edge effects : turn around the lattica in each direction
//also for absorbing edges : have to come from inside the lattice (loop while inside dispersal code)
int sortie_sup_circ_abs1(int coord, int dim) {
if(coord>dim) return((coord - (int) floor( (double) coord/ (double) dim - 1E-5)*dim));
else return(coord);
}
int sortie_inf_circ_abs1(int coord, int dim) {
if(coord<=0) return( (coord + (int) ceil(-((double) coord/ (double) dim) + 1E-5)*dim) );
else return(coord);
}
int sortie_sup_circ_abs2(int coord, int dim, int trans) {
if(coord>dim) return( (coord - (int) floor( (double) coord/ (double) dim - 1E-5)*dim + (trans) ) );
else return(coord + (trans));
}
int sortie_inf_circ_abs2(int coord, int dim, int trans) {
if(coord<=0) return( (coord + (int) ceil(- ((double) coord/ (double) dim) + 1E-5)*dim + (trans)) );
else return(coord + (trans));
}

//reflective edges = miror effect on edges
int sortie_sup_refl1(int coord, int dim) {
if(coord>dim) return( (2*dim-(coord)) );
else return(coord);
}
int sortie_inf_refl1(int coord, int dim) {
if(coord<=0) return( (2-(coord)) );
else return(coord);
}
int sortie_sup_refl2(int coord, int dim, int trans) {
if(coord>dim) return( (2*dim-(coord) + (trans) ) );
else return(coord + (trans));
}
int sortie_inf_refl2(int coord, int dim, int trans) {
if(coord<=0) return( (2-(coord) + (trans)) );
else return(coord + (trans));
}

/*****************/


/****************/
/*int file_tester(const int i) {
char Fishet[] = "atuer";
ffishet.open(Fishet,ios::out|ios::app);
ffishet<<i<<endl;
cout<<i<<endl<<flush;
ffishet.close();
return 0;
}*/
/****************/


/*****************/
//cree et/ou alloue la memoire pour cum2[][] de tabdis[][].cum2[][]
// avec dx_max et dy_max potentiellement diffÈrents
// a Fusionner avec
int AllocNotFixDisp(const int dxmax,const int dymax) {
if(currTVparsPtr->DispSpatiallyHeterog || currTVparsPtr->DensSpatiallyHeterog  || currTVparsPtr->forwardBarrier) {// si la densite ou la dispersion ne sont pas les meme partout dans cette phase demo


    tabdis=new disdis*[currTVparsPtr->dimRes1 + 2];//RL: remplacer par +1?
    if (!tabdis) nrerror("allocation failure 1 in structure disdis()");
    for (int i=0; i<=currTVparsPtr->dimRes1+1; i++) {
        tabdis[i]=new disdis[currTVparsPtr->dimRes2+2];
        if (!tabdis[i]) {
            printf("boucle %d",i);
            nrerror("allocation failure 2 in structure disdis()");
            getchar();
        }
    }//for i
    

    for (int i=0; i<=currTVparsPtr->dimRes1+1; i++) {
         for (int j=0; j<=currTVparsPtr->dimRes2+1; j++) {
#ifdef DEBUG
            tabdis[i][j].cum2.resize(2*dxmax+2);
            for (int k=0; k<2*dxmax+2; k++)	tabdis[i][j].cum2[k].resize(2*dymax+2);
#else
            //cout << "Allocate memory for tabdis[" << i << "][" << j << ".cum2" << ends;
            (tabdis[i][j].cum2)=new double*[2*dxmax+2];
            if (!tabdis[i][j].cum2) {
                printf("boucle %d",j);
                nrerror("allocation failure 3 in structure disdis()");
                getchar();
            }
            for (int k=0; k<2*dxmax+2; k++)	{
                //cout << "k=" << k << "(2*dymax+2=" << 2*dymax+2 << "); " << ends;
                tabdis[i][j].cum2[k]=new double[2*dymax+2];
                if (!tabdis[i][j].cum2[k]) {
                    printf("boucle %d",j);
                    nrerror("allocation failure 4 in structure disdis()");
                    getchar();
                }
            }//for k
            //cout << endl;
            //cout << endl;
#endif
        }//for  j
    }// for i
}//if
return 0;
}
/*****************/


/*****************/
//divers affichage Ècran au dÈbut du run a ameillorer...
void affichage_ecran()
{
#ifdef GOTO
	_gotoxy(0,0);
#endif

//printf("\n");
printf("================================================================\n");
printf("This is IBDSim v2.0 (Leblois et al. 2009 MolEcol Ressources) for the\n");
printf("simulation of genetic data under isolation by distance models\n");
printf("based in the algorithms described in Leblois et al. 2003 MBE,\n");
printf("Leblois et al. 2004 Genetics, Leblois et al. 2006 MolEcol\n");
printf("================================================================\n");
printf("Summary settings: generic file name is %s;\n",fichier_genepop.c_str());
if (nUniqMarker[0] < 2)
	cout << "Mutational model is ";
else
	cout << "The " << nUniqMarker[0] << " mutational models are ";
for (int i = 1; i <= nUniqMarker[0]; i++)
	cout << dispMut[i] << " ";
printf("\nSimulation of %d data set(s) each containing %d locus;",repet,n_locus);
if(!Specific_Sample_Designbool[0]) printf("\nSample: (%d x %d) * %d = %d ",dim_sample1[0],dim_sample2[0],dens_sample[0],dim_sample1[0]*dim_sample2[0]*dens_sample[0]);
    else printf("\nSample: (%d) * %d = %d ",Spec_SampleSize[0],dens_sample[0],Spec_SampleSize[0]*dens_sample[0]);
if(predispbool) {
    if(!Specific_Sample_Designbool[1])
        printf("\nPredisp Sample: (%d x %d) * %d = %d ",dim_sample1[1],dim_sample2[1],dens_sample[1],dim_sample1[1]*dim_sample2[1]*dens_sample[1]);
      else printf("\nPredisp Sample: (%d) * %d = %d ",Spec_SampleSize[1],dens_sample[1],Spec_SampleSize[1]*dens_sample[1]);
}
if(ploidy==2) printf("diploid individuals\n"); else printf("haploid individuals\n");
if (EdgeEffect=="circ") printf("evolving at G=0 on a %d x %d torus",TVpars[0].dimRes1,TVpars[0].dimRes2);
if (EdgeEffect=="abs") printf("evolving at G=0 on a %d x %d lattice with absorbing boundaries",TVpars[0].dimRes1,TVpars[0].dimRes2);
if (EdgeEffect=="refl") printf("evolving at G=0 on a %d x %d lattice with reflecting boundaries",TVpars[0].dimRes1,TVpars[0].dimRes2);
if(Specific_Density_Designbool) printf("\nwhere the # of individuals at each node is \ngiven in the '%s' file.",SpecificDensityFilename.c_str());
if(!Specific_Density_Designbool && TVpars[0].vide==1) printf("\nwhere each node has %d individuals",TVpars[0].initialDens);
if(!Specific_Density_Designbool && TVpars[0].vide==2) printf("\nwhere one every second node has %d individuals",TVpars[0].initialDens);
if(!Specific_Density_Designbool && TVpars[0].vide==3) printf("\nwhere one every third node has %d individuals",TVpars[0].initialDens);
printf("\nDispersal settings are summarized in the 'simul_pars.txt' file ");
printf("\n================================================================\n");
}/*fin affichage ecran*/



/*---------initialisations pour deux dimensions---------------------------*/


/*****************/
//initialisation des position geographiques des individus ÈchantillonnÈs
void ini_tableaux()
{int p,k,o,d;// j et d virÈs FER 28/08/2005 rajoutÈ RL car importants :
// d= si l'on Èchantillonne plusieurs individus sur un meme deme on veut qu'ils aient les meme coordonnÈes
//j= pour prendre en compte la deuxieme dimension de l'Èchantillon = la position en x est incrÈmentÈ uniquement
// quand on a placÈ tout les individus voulus en diffÈrents y

using namespace NS_coal_tree;

for(int i=1;i<=(2*n_genes_total);i++) coord_ori[i][x]=coord_ori[i][y]=0;
for(int i=1;i<=(n_genes_total);i++)  no_noeud_restant[i]=i;

p=1; k=xmin_sample[0];
for(int i=1;i<=(2*n_genes_total);)
 {if(i<=n_genes[0])
	{o=ymin_sample[0];
	for(int j=1;j<=dim_sample2[0];j++)
	 {for(d=1;d<=dens_sample[0];d++)
		{coord_noeud[p][x]=coord_individus[p][x]=k;
		 if(ploidy==2) coord_noeud[p+1][x]=coord_individus[p+1][x]=k;
		 coord_noeud[p][y]=coord_individus[p][y]=o;
		 if(ploidy==2) coord_noeud[p+1][y]=coord_individus[p+1][y]=o;
		 p+=ploidy;
		}
	 o+=vide_sampleY[0];
	 }
	i=p;
	}
 else if(i<=n_genes_total)
	{if(i==(n_genes[0]+1)) k=xmin_sample[1];o=ymin_sample[1];
        for(int j=1;j<=dim_sample2[1];j++) {
            for(d=1;d<=dens_sample[1];d++) {
                //cout << "i=" << i << "; p=" << p << "; k=" << k << endl;
                coord_noeud[p][x]=coord_individus[p][x]=k;
                if(ploidy==2) coord_noeud[p+1][x]=coord_individus[p+1][x]=k;
                coord_noeud[p][y]=coord_individus[p][y]=o;
                if(ploidy==2) coord_noeud[p+1][y]=coord_individus[p+1][y]=o;
                p+=ploidy;
            }
            o+=vide_sampleY[1];
        }
        i=p;
    } else {
        coord_noeud[i][x]=coord_individus[i][x]=0;
        if(ploidy==2) coord_noeud[i+1][x]=coord_individus[i+1][x]=0;
        coord_noeud[i][y]=coord_individus[i][y]=0;
        if(ploidy==2) coord_noeud[i+1][y]=coord_individus[i+1][y]=0;
        i=i+ploidy;
    }
  if(i<=n_genes[0])
      k+=vide_sampleX[0];
   else if(i<=n_genes_total)
            k+=vide_sampleX[1];
 }
//verif
//for(int i=1;i<=(n_genes_total);i++) {
//    printf("\n coord_individu[%d]=[%d,%d]",i,coord_individus[i][x],coord_individus[i][y]);
//    getchar();
//}
}

/*****************/
//procedure speciale pour positionnement specifiques des individus sur le reseau
void ini_tableaux_specifique()
{int i,p,d,k; //j, d vires FER 28/08/2005 RL rajoutÈ  :meme chose que plus haut
using namespace NS_coal_tree;

for(i=1;i<=2*n_genes_total;i++) coord_ori[i][x]=coord_ori[i][y]=0;
for(i=1;i<=n_genes_total;i++)  no_noeud_restant[i]=i;

//cout << "\n\n\n\nBeginning of ini_tableaux(), n_genes[0]=" << n_genes[0] << "; Specific Sample Coordinates are :" << endl;
//for(int i=0;i<Spec_SampleSize;i++)	cout << Spec_Sample_Coord[0][x][i] << " " << Spec_Sample_Coord[0][y][i] << endl;
//cout << "n_genes_predisp=" << n_genes[1] << "; Predisp_Specific_Sample_Coordinates are :" << endl;
//for(int i=0;i<Spec_SampleSize[1];i++)	cout << Spec_Sample_Coord[1][x][i] << " " << Spec_Sample_Coord[1][y][i] << endl;


for(i=0,p=1;p<=2*n_genes_total;) {
    if(p<=n_genes[0]) {
        //cout << "p=" << p << "i=" << i << endl;
            for(d=1;d<=dens_sample[0];d++) {
                coord_noeud[p][x]=coord_individus[p][x]=Spec_Sample_Coord[0][x][i];
                if(ploidy==2) coord_noeud[p+1][x]=coord_individus[p+1][x]=Spec_Sample_Coord[0][x][i];
                coord_noeud[p][y]=coord_individus[p][y]=Spec_Sample_Coord[0][y][i];
                if(ploidy==2) coord_noeud[p+1][y]=coord_individus[p+1][y]=Spec_Sample_Coord[0][y][i];
                //cout << "coord_individu[" << p << "]=[" << coord_individus[p][x] << "," << coord_individus[p][y] << "]" << endl;
                p+=ploidy;
    	   }
    	i++;
    }
    else if(p<=n_genes_total) {
        k=i-n_genes[0];
        //cout << "p=" << p << "; i=" << i << "; k=" << k << endl;
        for(d=1;d<=dens_sample[1];d++) {
            coord_noeud[p][x]=coord_individus[p][x]=Spec_Sample_Coord[1][x][k];
            if(ploidy==2) coord_noeud[p+1][x]=coord_individus[p+1][x]=Spec_Sample_Coord[1][x][k];
            coord_noeud[p][y]=coord_individus[p][y]=Spec_Sample_Coord[1][y][k];
            if(ploidy==2) coord_noeud[p+1][y]=coord_individus[p+1][y]=Spec_Sample_Coord[1][y][k];
            //cout << "coord_individu[" << p << "]=[" << coord_individus[p][x] << "," << coord_individus[p][y] << "]" << endl;
            p+=ploidy;
    	   }
        i++;
    } else {
        coord_noeud[p][x]=coord_individus[p][x]=0;
	    if(ploidy==2) coord_noeud[p+1][x]=coord_individus[p+1][x]=0;
	    coord_noeud[p][y]=coord_individus[p][y]=0;
	    if(ploidy==2) coord_noeud[p+1][y]=coord_individus[p+1][y]=0;
	    p+=ploidy;
	}
}
//cout << "\nEnd of ini_tableaux_specifiques(), we have :" << endl;
//for(i=1;i<=2*n_genes_total;i++) {
//	cout << "coord_individu[" << i << "]=[" << coord_individus[i][x] << "," << coord_individus[i][y] << "]" << endl;
//}
//getchar();
}

/*****************/
//initialisation des stats utilisÈes sur plusieurs repetitions
void ini_moy_glob() {
    using namespace NS_diagnostic_tables;
    using namespace NS_translation;
    cum_moyx=cum_moyy=0;
    int sampleNb=1;

if(predispbool) sampleNb=2;

if(effective_disp) {
	for(int i=0;i<=2*grande_dim;i++) effective_moy[i][x]=effective_moy[i][y]=0;
	for(int i=0;i<dim_reseau1+1;i++) for(int j=0;j<dim_reseau2+1;j++) effectiveImmigration_moy[i][j]=cumulEffectiveImmigration_moy[i][j]=0;
}
mono=0;
tooManyMut=0;
//compte_mut_glob=0;//RL 022017 not used
//mut_moy_glob=0.0;//RL 022017 not used
//mut_moy_glob_var=0.0;//RL 022017 not used
mrca_moy_glob=0.0;
coa_moy_glob=0.0;

for(int a=0;a<sampleNb;a++) {
    HexpNei_moy_glob[a]=0.0;
    HexpNei_moy_glob_var[a]=0.0;
    hetero_moy_glob[a]=0.0;
    M_moy_glob[a]=0.0;
    M_moy_glob_var[a]=0.0;
    Var_moy_glob[a]=0.0;
    Var_moy_glob_var[a]=0.0;
    fis_moy_glob[a]=0.0;
    fis_moy_denom_glob[a]=0.0;
    fis_moy_numer_glob[a]=0.0;
    if( !Specific_Sample_Designbool[0]  && !Specific_Sample_Designbool[1]) {
        for(int j=0;j</*max(TVpars[0].dimRes1,*/max(max(dim_sample1[0],dim_sample2[0]), max(dim_sample1[1],dim_sample2[1]) ) /*)*/;j++) Q1_moy_glob[a][j]=0.0;
        Qr_mean_moy_glob[a]=0.0;
        Qb_meanAllPairs_moy_glob[a]=0.0;
        Qw_meanAllInd_moy_glob[a]=0.0;
        Qind_moy_glob[a]=0.0;
    }
}
}

/*****************/
//initialisation des stats utilisÈes sur une repetition multilocus
void ini_moy() {
    using namespace NS_diagnostic_tables;
    using namespace NS_translation;
    int sampleNb=1;

if(predispbool) sampleNb=2;
mean_allele_nmbr=0.0;
if(effective_disp) {
	for(int i=0;i<=2*grande_dim;i++) effective[i][x]=effective[i][y]=0;
	for(int i=0;i<dim_reseau1+1;i++) for(int j=0;j<dim_reseau2+1;j++) effectiveImmigration[i][j]=cumulEffectiveImmigration[i][j]=0;
}
for(int i=0;i<n_locus;i++){
	mrca[i]=0;
	allele_nmbr[i]=0;
	allele_max[i]=-1000;
	allele_min[i]=1000;
}

mrca_moy=0.0;
temps_coa=compteur3=0;
temps_coa_moy=0.0;

for(int a=0;a<sampleNb;a++) {
    for(int i=0;i<n_locus;i++){
        HexpNei[a][i]=0.0;
        hetero[a][i]=0.0;
        Var[a][i]=0.0;
        M[a][i]=0.0;
        fis[a][i]=0.0;
    }
    HexpNei_moy[a]=0.0;
    hetero_moy[a]=0.0;
    M_moy[a]=0.0;
    Var_moy[a]=0.0;
    fis_moy2[a]=0.0;
    fis_moy[a]=0.0;
    fis_moy_denom[a]=0.0;
    fis_moy_numer[a]=0.0;
    if( !(Specific_Sample_Designbool[0] || Specific_Sample_Designbool[1]) ) {
        for(int i=0;i<max( max(dim_sample1[0],dim_sample2[0]), max(dim_sample1[1],dim_sample2[1]) );i++) {
            for(int j=0;j<n_locus;j++) Q1[a][j][i]=0.0;
            Q1_moy[a][i]=0.0;
        }
        Qind_moy[a]=0.0;
        Qr_mean_moy[a]=0.0;
    }

}
}

/*****************/
//initialisation de la structure noeud
void ini_struct_noeud()
{ for(int i=1;i<=2*n_genes_total;i++)
  {for(int j=0;j<n_locus;j++)
	{noeud[i][j].etat_allele=111;
	noeud[i][j].generation=0;
	noeud[i][j].descendant.resize(0);
    noeud[i][j].ancetre=0;
    noeud[i][j].sequence = "";
	//for(int k=0;k<10;k++) noeud[i][j].descendant[k]=0;
	}
  }
}

/*****************/
//re-initialisation de la structure noeud pour le locus en cours de simulation
//effectuÈe si le locus venant d'etre simulÈ est monophorphe ou a trop de mutation et ne nous interesse pas
void ini_struct_noeud2()
{int j;
 for(int i=1;i<=2*n_genes_total;i++)
  {j=locus;
  noeud[i][j].etat_allele=111;

	noeud[i][j].generation=0;
	noeud[i][j].descendant.resize(0);
    noeud[i][j].ancetre=0;
    noeud[i][j].sequence = "";

//  noeud[i][j].generation=noeud[i][j].
//  nbre_descendant=noeud[i][j].ancetre=0;
//  for(int k=0;k<10;k++) noeud[i][j].descendant[k]=0;
  }
}


/*-------------------------------------------------------------------------*/
/*--------------------dispersion-------------------------------------------*/

/*****************/
//initialisation des diffÈrents parametres dÈmographiques en fonction des changement temporels planifiÈs
//prÈcise aussi si l'on doit ou non changer certains parametres de migration (remet deja et deja1 a 0)
//crÈÈ une translation en x et y si les genes suivis doivent etre placÈs sur un nouveau reseau plus grand


int reset_demo_params(int newPhase) { //FR-> this should completely replace migra_tps...
using namespace NS_translation;
    norm=0.0;
    deja=deja1=0;
    deja3=1;

    model_mig=currTVparsPtr->mod;
    if(currTVparsPtr->DispSpatiallyHeterog) {
        model_mig_zone=currTVparsPtr->mod_zone;
        //cout << "DispSptiallyHeterog=T; model_mig_zone set to " << model_mig_zone << endl;
    }

   	TVpars[newPhase].dimRes1=TVpars[newPhase].initialDimRes1;TVpars[newPhase].dimRes2=TVpars[newPhase].initialDimRes2;

    TVpars[newPhase].SetForwardDispersalDistributions();// recalcule distr disp,

    if( (cmp_nocase(currTVparsPtr->ContDemeSizeChange,"None") != 0 ) || (cmp_nocase(currTVparsPtr->ContLatticeSizeChange,"None") != 0 ) ) {
            currTVparsPtr->setGrowthPars(&(TVpars[newPhase+1]));
    }
    if(currTVparsPtr->DensSpatiallyHeterog || currTVparsPtr->DispSpatiallyHeterog || currTVparsPtr->forwardBarrier) currTVparsPtr->computeCumuls();

    if (newPhase==0) {
		 translationx=translationy=0;
    } else {
        if(currTVparsPtr->dimRes1>TVpars[newPhase-1].dimRes1) {//translation sur x
                    translationx=(currTVparsPtr->dimRes1-TVpars[newPhase-1].dimRes1)/2;
                    if(random_translation==1) translationx=int(floor(alea()*translationx));
                    translationx-=translationx % currTVparsPtr->vide;
        }
        if(currTVparsPtr->dimRes2>TVpars[newPhase-1].dimRes2) {//translation sur y
                    translationy=(currTVparsPtr->dimRes2-TVpars[newPhase-1].dimRes2)/2;
                    if(random_translation==1) translationy=int(floor(alea()*translationy));
                    translationy-=translationy % currTVparsPtr->vide;
        }
    }
return(0);
}





/*****************/
//dÈfini les probabilitÈ de migration en fonction du modele de dispersion choisi
//dÈfini dx_max et dy_max = distance de dispersion maximale consideree
//accessoirement calcul sigma et kurtosis pour vÈrification
//appelÈe dans new_disp()
void CTimeVaryingParams::SetForwardDispersalDistributions() {
    vector<long double>SichelDispPars(4);
    float Dhs2=numeric_limits<float>::quiet_NaN();
    float Nhm=numeric_limits<float>::quiet_NaN();
//cout<<"debut forw"<<flush;
    //double mig,kpar;
    int pas=1;
/*gotoxy(1,20);
printf("model mig:%c",model_mig);*/
//if (locus ==2) cout<<"deja="<<deja<<flush;

    if(model_mig=='a') {model_mig='g';GeoG=0;Mig=1./3.;}
    if(model_mig=='b') {model_mig='g';GeoG=0;}
    if(model_mig=='1') {model_mig='g';GeoG=0;Mig=2./3.;}


    if(deja==0) {// remis a 0 en de multiples endroits eg pour Gn=1 dans reset_demo_params
        norm=0.0;
        if (model_mig=='P' || (model_mig=='g' && GeoG>0) || model_mig=='S') {
            if(compareMath) {
                dx_max=dist_max;
            } else {
                dx_max=dimRes1-1;
            }
        } else {
            if (model_mig=='g' && GeoG==0) {
                dx_max=1;
            } else if (model_mig=='0') {/*pareto modifiee sigma≤=4 pour analyses a petite distance cf FR 2000 (c) (d)*/
                dx_max=15;// parameters Mig & Pareto_Shape specified below because of specific mig[1] and mig[2] values
            } else if (model_mig=='2') {/*kurtosis forte cf waser dipodomys sigma2=1*/
                //dx_max=49;model_mig='P';Mig=0.599985;Pareto_Shape=3.79809435;
                dx_max=49;model_mig='P';Mig=0.6;Pareto_Shape=3.7980943485447316; //checked with mathematica
            } else if (model_mig=='3') {/*sigma2=100*/
                dx_max=48;model_mig='P';Mig=0.6;Pareto_Shape=1.246085;
            } else if (model_mig=='4') {/*sigma=1 pour un noeud sur deux*/
                dx_max=48;model_mig='P';pas=2;Mig=1.-0.824095;Pareto_Shape=4.1078739681;
            } else if (model_mig=='5') {/*sigma=1 pour un noeud sur trois*/
                dx_max=48;model_mig='P';pas=3;Mig=1.-0.913378;Pareto_Shape=4.43153111547;
            } else if (model_mig=='6') {/*sigma2=20*/
                dx_max=48;model_mig='P';Mig=0.719326;Pareto_Shape=2.0334337244;
            } else if (model_mig=='7') {/*sigma2=10*/
                dx_max=49;model_mig='P';Mig=0.702504;Pareto_Shape=2.313010658;
            } else if (model_mig=='8') {/*sigma≤=4 pour un noeud sur trois*/
                dx_max=48;model_mig='P';pas=3;Mig=1.-0.678842;Pareto_Shape=4.1598694692;
            } else if (model_mig=='9') {/*"sigma2=4 mais sur 48 pas"*/
                dx_max=48;model_mig='P';Mig=0.700013;Pareto_Shape=2.74376568017753;
            }
            if( ! compareMath) dx_max=dimRes1-1;
        }
        // setting dymax
        if (dimRes2==1) dy_max=0; else if(compareMath) {
            dy_max=dx_max;
        } else {
            dy_max=min(dimRes2-1,dx_max);
        }
		dx_min= - dx_max;
		dy_min= - dy_max; //obligatoire
        migra=ldvector(dx_max+1);
        for(int i=0;i<=dx_max;i++) migra[i]=0.0;
        AllocNotFixDisp(dx_max,dy_max);// allocation tabdis[][] (new 052015 RL) et tabdis[][].cum2[...] useful if spatial heterogeneities ("fixe=0")
        //FR->RL : ...ELSE... :
        cum_fixex=ldvector(2*dx_max+2); // backward distribs //FR -> RL: I would rather say it's the numerator of the expresion for the backward prob.
        cum_fixey=ldvector(2*dy_max+2);
        if(!Simple1DProductDispBool) {//RL072014 was cmp_nocase(twoD_disp,"1DProductWithoutm0")==0
            backwardPhilo.resize(dimRes1);
            for (std::vector<vector<double> >::iterator it=backwardPhilo.begin();it!=backwardPhilo.end();it++) {
                it->resize(dimRes2/vide); // FR 0610
            }
        }
		/** @@@@@@@@@@@@@@@@@@ SWITCH :**/
        switch(model_mig){
            case'g': {/*geometric*/
		//*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		// NOUVELLE VERSION
//cout<<dimRes1<<" "<<dx_max<<flush;getchar();
		  //printf("\n distrib geometric m=%f; alpha=%f",1-migra[0],alpha);
                if(!Simple1DProductDispBool) {//RL072014 wascmp_nocase(twoD_disp,"1DProductWithoutm0")==0) {
                    for(int i=0;i<=dx_max;i++) norm+=(double)pow(GeoG,i); //+1 if kpar=1
                    for(int i=0;i<=dx_max;i++) {
                        migra[i]=(double) pow(GeoG,i)/(2.0*norm-1.0); //1/(2*(dx_max+1)-1)= 1/(2*dx_max+1)
                  // AT THIS STEP migra["-dx_max"] to [dx_max] ([0] included) should sum to 1; migra[0] is NOT 1-mig.
                  /** before correction and after correction
                    m11 / m01 = g (deplacements de (1,1) vs (1,0)
                    car migra[1]=g migra[0]

                    dans le cas suivant ce n'est pas vrai, meme si on fait tourner oneDProductWithoutm0() apres
                  **/
                    }
                } else {
                    migra[0]=(double) (1.0-Mig);
                    if (GeoG==1) { //disp en ile ATTENTION PAS CONTINUITE ENTRE LES 2 CAS
                        //migra[0]+=(double) Mig/(dx_max);//RL 122012 : pour avoir un vrai model en ile, RL 092014 enleve car vraiment bizarre...
                        //for(int i=1;i<=dx_max;i++) migra[i]=(double) Mig/(dx_max-1); //aucune perte sur bord !!
						for(int i=1;i<=dx_max;i++) migra[i]=(double) Mig/(2*dx_max); //aucune perte sur bord !! modifs RL 15062011, voir avec FR   RL->FR
                    } else {
                        for(int i=1;i<=dx_max;i++) norm+=(double) pow(GeoG,i-1.0);
                        for(int i=1;i<=dx_max;i++)
                            migra[i]=(double) pow(GeoG,i-1.0)*Mig/(2.0*norm); //so migra[<>0] sum to 1-migra[0]=mig. "Cross migration pattern" for low mig
                    }
                }
        /****** in the original version of IBDsim there is no further change on migra[] ******/
        /****** However, if we want to control the maximal dispersal rate, we will change the migra[]: ******/
                if(!Simple1DProductDispBool) {//RL072014 wascmp_nocase(twoD_disp,"1DProductWithoutm0")==0) {
		        // this writes in simul_pars...
                    oneDProductWithoutm0(); /*** also meaningful in 1D !!!!!!!!!!!!!!!!! ***/
                } else { //computation of maximal dispersal rate in original version of IBDsim
// rien ?
                }
                deja=1; //FR->RL "deja=1" avant chanque break: on peut mettre ce code juste apres la sortie du switch + renommer cette variable...
            }
            break;
            case'0': {/*pareto sigma≤=4 pour analyses a petite distance cf FR (c) (d)*/
		//AllocNotFixDisp(dx_max);
                Mig=0.3;
                Pareto_Shape=2.51829;
                migra[0]=(double)(1-Mig);
                migra[1]=(double) 0.06;
                migra[2]=(double) 0.03;
                for(int i=3;i<=dx_max;i++) norm+=((double) (pow(1.0/i,Pareto_Shape)));
                for(int i=3;i<=dx_max;i++)
                    migra[i]=((double) ((Mig-2*migra[1]-2*migra[2])*pow(1.0/i,Pareto_Shape)/(2*norm)));
                deja=1;
            }
            break;
            case'P': {/*Pareto standard*/
                migra[0]=(double)(1.0-Mig);
                for(int i=pas;i<=dx_max;i+=pas) norm+=pow(1.0/(i/pas),Pareto_Shape);
                for(int i=pas;i<=dx_max;i+=pas) migra[i]=((double) (Mig*pow(1.0/(i/pas),Pareto_Shape)/(2*norm)));
                deja=1;
            }
            break;
            case'S': {/*Sichel and inverse gamma mixtures*/
		//cout<<"apres AllocNotFixDisp"<<flush;
                vector<long double>costable;
                unsigned int lPSON;
                long double temp;
		//SichelDispPars[0]: Sichel gamma
                SichelDispPars[0]=Sichel_Gamma;
                if (SichelDispPars[0]>0) {cout<<"gamma arg. must be negative";getchar();exit(-1);}
		//SichelDispPars[1]: Sichel xi (or reciprocal gamma's kappa\sim xi omega)
                SichelDispPars[1]=Sichel_Xi;
		//SichelDispPars[2]: Sichel omega (or something negative for reciprocal gamma)
                SichelDispPars[2]=Sichel_Omega;
		/*SichelDispPars[3]: island (total) dispersal rate: with proba SichelDispPars[3]
		the propagule disperses somewhere (staying excluded)
		with proba 1-SichelDispPars[3] the propagule follows the Sichel dist.*/
                SichelDispPars[3]=0.; //non charfunc fraction
		//what follows somehow assumes dim_reseau1=dim_reseau2 ?
                if (EdgeEffect=="circ") {lPSON=dimRes1;} /*ca donne des transformÈes sur un habitat fermÈ
		=> Sicheltable[->lPSON]->Sicheltable[1]... pas terrible quand on veut calculer var et kurtosis;
		ou bien il faut que lPSON>2 dxmax...*/
                else {lPSON=100*dx_max;} // pour Èviter effet repliement
                long double DeuxPiSurNx=2.*PI/lPSON;
                costable.resize(lPSON);
                for(long int ii=0;ii<lPSON;ii++) costable[ii]=cos(DeuxPiSurNx*ii);
                Sicheltable=ldvector(lPSON/2+1);
                for(long int ii=0;ii<lPSON/2+1;ii++) {  //la valeur en lPSON/2+1 est utilisÈe
                    Sicheltable[ii]=Sichelcharfunc((long double)(2.*PI*ii/lPSON),SichelDispPars);
                }
                for(int _dx=0;_dx<=dx_max;_dx++) {
                    temp=1.0; // c'est charfn[0]
                    for(long int ii=1;ii<lPSON/2;ii++) {
				//2 parce que les 2 termes sym de la transformÈe. mais c'est la valeur en |_dx|
                        temp+=2*Sicheltable[ii]*costable[(_dx*ii)%lPSON];
				//if ((_dx==3500) && (ii%500==0)) {cout<<ii<<" "<<temp<<endl;getchar();}
                    }
                    if(lPSON%2==0) temp+=Sicheltable[lPSON/2]*costable[(_dx%2)*lPSON/2];
                    temp/=lPSON;
                    migra[_dx]=temp;  // migration selon Sichel
                }
                    // SichelDispPars[3] is a proba of axial uniform dispersal, ce n'est pas Mig !
                for(long int _dx=0;_dx<=dx_max;_dx++) {
                    temp=migra[_dx]*(1.-SichelDispPars[3]);
                    if (_dx>0) temp+=SichelDispPars[3]/lPSON;
                    norm+=(_dx==0?1.:2.)*temp;
                    migra[_dx]=temp;
                } // for _dx ...dx_max
                for(int _dx=0;_dx<=dx_max;_dx++) migra[_dx]/=norm;
                // at this point (if there is no uniform dispersal) the axial dispersal terms (0 included) are still in proportion to those of a true Sichel distribution
                if(!Simple1DProductDispBool) {//RL072014 wascmp_nocase(twoD_disp,"1DProductWithoutm0")==0) {
                    // only then a Mig value is used
                    oneDProductWithoutm0(); /*** also meaningful in 1D !!!!!!!!!!!!!!!!! ***/
                }
                free_ldvector(Sicheltable);
                deja=1;
            } // end case S
            break;
            cerr<<"(!) From SetForwardDispersalDistributions(): Unknown migration model (nicknamed '"<<model_mig<<"'. I exit."<<endl;
            if (cinGetOnError) cin.get();
            exit(-1);
        } // END OF SWITCH

        /** ECRITURE SIMULPARS **/
        if( ! simulparsWrittenBool) {
            if(phase>0){
                simulpars<<"=============="<<endl<<endl;
                simulpars<<"Demographic phase: "<< phase <<endl;
                simulpars<<"At G=" << currentGeneration << ": "<<endl;
                simulpars<<"lattice dimensions are (LatticeSizeX * LatticeSizeY): "<<dimRes1<<" x "<<dimRes2<<endl;
                if(vide!=1)  simulpars<<"with "<< (vide-1) << " nodes over " <<vide
                    << " being empty, so that real gene density is : "<< ((double) 2.0*initialDens/vide)<<endl;
                if(ploidy==2) simulpars<<"Diploid individuals : nbr of genes per non empty lattice node (2 * IndsPerPop): "<<2*initialDens<<endl;
                else simulpars<<"Haploid individuals : nbr of genes per lattice node (1* IndsPerPop): "<<initialDens<<endl;
                if(zone && DensSpatiallyHeterog) {
                    simulpars<< "xmin_zone=" << xmin_zone << " and xmax_zone=" << xmax_zone << endl;
                    simulpars<< "ymin_zone=" << ymin_zone << " and ymax_zone=" << ymax_zone << endl;
                    simulpars<< "voidNode in Zone =" << vide_zone << endl;
                    simulpars<< "Density in Zone (nbr of genes) =" << ploidy*dens_zone << endl;
                }

                if(backwardBarrier || forwardBarrier) {
                    if(forwardBarrier) simulpars<<"Presence of a real forward barrier to gene flow with coordinates :" <<endl;
                    else simulpars<<"Presence of an approximate (i.e. backward) barrier to gene flow with coordinates :" <<endl;
                    simulpars << "x1_barrier=" << x1_barrier << " and x2_barrier=" << x2_barrier << endl;
                    simulpars << "y1_barrier=" << y1_barrier << " and y2_barrier=" << y2_barrier << endl;
                    if (barrierCrossingRateUp == barrierCrossingRateDown) simulpars << "and a BarrierCrossingRate =" << barrierCrossingRateUp << endl;
                    else if(forwardBarrier) {
                        simulpars << "with forward BarrierCrossingRateUp =" << barrierCrossingRateUp << endl;
                        simulpars << "and forward BarrierCrossingRateDown =" << barrierCrossingRateDown << endl;
                    } else if(backwardBarrier) {
                        simulpars << "with backward BarrierCrossingRateUp =" << barrierCrossingRateUp << endl;
                        simulpars << "and backward BarrierCrossingRateDown =" << barrierCrossingRateDown << endl;
                    }
                    simulpars << "Be carefull with forward and backward barrier settings" << endl;
                    simulpars << "and the interpretation of asymetrical BarrierCrossingRateUp and BarrierCrossingRateDown" << endl;
                }
                if(ContDemeSizeChange!="None") {
                    simulpars<< "ContDemeSizeChange=" << ContDemeSizeChange << endl;
                }
                if(ContLatticeSizeChange!="None") {
                    simulpars<< "ContLatticeSizeChange=" << ContLatticeSizeChange << endl;
                }

                simulpars<<endl<<endl;
            }
            if(!Simple1DProductDispBool) {//RL072014 wascmp_nocase(twoD_disp,"1DProductWithoutm0")==0) {
                simulpars<<"Dispersal proba (if 2D: NOT axial !) = TotalEmmigrationRate: "<<Mig<<endl<<endl;
                simulpars<<"Maximal immigrant gene number among demes, 2N_h m "<<2*ploidy*initialDens*Mig<<endl;
            } else {
                simulpars<<"Dispersal model: "<< mod <<endl;
                simulpars<<"Axial dispersal proba (TotalEmmigrationRate): "<< Mig <<endl;
            }

            switch(model_mig){
                case'g': {/*geometric*/
                    if (GeoG==0) {
                        simulpars<<"Stepping stone model "<<endl ;
                    } else if (GeoG==1) {
                        simulpars<<"Island model of dispersal "<<endl ;
                    } else {
                        simulpars<<"Geometric dispersal model "<<endl ;
                        simulpars<<"GeometricShape: "<<GeoG<<endl;
                    }
                }
                    break;
                case 'S': {
                    if (SichelDispPars[2]<0) {
                        simulpars<<"Inverse Gamma mixture (Chesson & Lee 2005) "<<endl ;
                        if (SichelDispPars[3]>0) simulpars<<"also including a fraction "<<SichelDispPars[3]<<" of random dispersal [excluding focal patch!]: "<<endl;
                        simulpars<<"gamma="<<SichelDispPars[0]<<endl ;
                        simulpars<<"kappa (i.e. lim (omega xi))= "<<SichelDispPars[1]<<endl ;
                    } else {
                        simulpars<<"Sichel dispersal distribution model (Chesson & Lee 2005)"<<endl;
                        simulpars<<"gamma="<<SichelDispPars[0]<<endl ;
                        simulpars<<"xi="<<SichelDispPars[1]<<endl ;
                        simulpars<<"omega="<<SichelDispPars[2]<<endl ;
                    }
                }
                    break;
                case 'P': {
                    simulpars<<"Pareto model "<<endl ;
                    simulpars<<"ParetoShape: "<<Pareto_Shape<<endl;
                }
                    break;
            }

            long double sig2;
            long double boundedC=0,unboundedC=0;
            if (model_mig=='S') {
		//cout<<"non charfunc fraction: "<<SichelDispPars[3]<<endl;
                if (SichelDispPars[2]<0) {
                    sig2=-(1.-SichelDispPars[3])*SichelDispPars[1]/(2.*(1.+SichelDispPars[0]));
                } else {
                    sig2=-(1.-SichelDispPars[3])*SichelDispPars[1]*bessel_k(1.+SichelDispPars[0],SichelDispPars[2],1)
                                                    /(2.*bessel_k(SichelDispPars[0],SichelDispPars[2],1));
                }
                Dhs2=ploidy*TVpars[0].initialDens*sig2;
            } else if (model_mig=='g') {
                if(!Simple1DProductDispBool) {//RL072014 wascmp_nocase(twoD_disp,"1DProductWithoutm0")==0) {
                    long double condsig2;
                    if (dimRes2/vide>1) { //2D
			          /*cf code R
			            condaxialS2fromg<-function(gv) {\n\
                            return((1+gv)/((2-gv)*(1-gv)^2))\n\
                        }\n\                      */
                        condsig2=(1.+GeoG)/((2-GeoG)*(1.-GeoG)*(1.-GeoG));
                    } else {
                        condsig2=(1.+GeoG)/((1.-GeoG)*(1.-GeoG));
                    }
                    Nhm=ploidy*TVpars[0].initialDens*Mig;
                    Dhs2=Nhm*condsig2;
                    sig2=condsig2*Mig;
                } else {
                    sig2=Mig*(1.+GeoG)/((1.-GeoG)*(1.-GeoG));
                    Dhs2=ploidy*TVpars[0].initialDens*sig2;
                }
            }

            int maxdd=max(dimRes1,dimRes2)-1;
            long double meanInddmax=0;
            long double sigInddmax=0;
            long double kurtInddmax=0;
            if(!Simple1DProductDispBool) {//RL072014 wascmp_nocase(twoD_disp,"1DProductWithoutm0")==0) {
                simulpars<<endl<<"Marginal terms and checksum "<<endl;
                if(dimRes2/vide>1) {
                    simulpars<<"   The following computation assumes that the habitat is square-shaped"<<endl;
                    simulpars<<"   and reconstructs the marginal backward distribution into the most central demes(s)"<<endl;
                    simulpars<<"   The following computation assumes that the habitat is square-shaped."<<endl;
                    simulpars<<"   The bounded cumulative distribution should then reach 1."<<endl;
                }
                long double sumNonzero;
                long double proba;
                boundedC=unboundedC= 1-Mig;
            simulpars<<"\nBounds (dx_max,dy_max)=("<<dx_max<<","<<dy_max<<")"<<endl<<endl ;
                simulpars<<"                          cumulative distributions"<<endl;
                simulpars<<"distances     probability    bounded  unbounded"<<endl;
                simulpars<<"x=y=0:         "<<setw(10)<<1-Mig<<" "<<setw(10)<<boundedC<<" "<<setw(10)<<unboundedC<<endl;
                if(dimRes2/vide>1) { // this test matches the one in the sampling algorithm
			         /**Now the maximal \sum\sum_{x,y\neq0,0}=Mig
			         =>\sum\sum_{x,y}-migra[0]^2=(\sum_x)^2-migra[0]^2=Mig
			         =>\sum_{y<>0}= sqrt[Mig+migra[0]^2]-migra[0]
			         **/
                    sumNonzero=sqrt(Mig+migra[0]*migra[0])-migra[0];
                    boundedC+=migra[0]*sumNonzero;
                    unboundedC=boundedC;
                    simulpars<<"x=0,y<>0:      "<<setw(10)<< migra[0]*sumNonzero<<" "<<setw(10)<<boundedC<<" "<<setw(10)<<unboundedC<<endl;
                } //else sumNonzero=1; //RL 062018, unused
                for(int i=1;i<=maxdd;i++) { // up to maxdd=dimRes1-1
                        /** by construction,
                            (1D)    \sum_x\neq(0) migra[x] =Mig

                            (2D)    \sum_x,y\neq(0,0) migra[x] migra[y]=Mig
                                      \sum_x,y migra[x] migra[y]=Mig+migra[0]^2
                                      \sum_y migra[y]=sqrt(Mig+migra[0]^2)

                            hence the marginal term in x!=0,
                            \sum_y migra[x] migra[y] = migra[x] \sum_y migra[y] = migra[x] sqrt(Mig+migra[0]^2)
                        **/
                    double euh;
                    if(dimRes2/vide>1) {
                        euh=sqrt(Mig+migra[0]*migra[0]);
                    } else
                        euh=1;
                    proba=migra[i]*euh;
                    unboundedC+=2*proba;
                    if (i<=maxdd/2) /*integer div*/ {
                        boundedC=unboundedC;
                        meanInddmax+=2*i*proba;
                        sigInddmax+=2*i*i*proba;
                        kurtInddmax+=2*i*i*i*i*proba;
                    } else if (i==(maxdd+1)/2) /*integer div*/ {
                        int fac=2-(maxdd%2);
                        boundedC+=fac*proba;
                        meanInddmax+=fac*i*proba;
                        sigInddmax+=fac*i*i*proba;
                        kurtInddmax+=fac*i*i*i*i*proba;
                    }
                        /**ELSE  boundedC not further increased*/
                    simulpars<<"x="<<setw(5)<<i<<":      "<<setw(10)<<proba<<" "<<setw(10)<<boundedC<<" "<<setw(10)<<unboundedC<<endl;
                }
            } else { //pas 1D...withoutm0
				simulpars<<"Bounds (dx_max,dy_max)=("<<dx_max<<","<<dy_max<<")"<<endl;
                simulpars << "Marginal dispersal distribution and checksum:" << endl;
                double boundedC= - migra[0];
                int DistMax;
                if(compareMath) DistMax=dx_max; else DistMax=maxdd;
                for(int i=0;i<=DistMax;i++) { // up to maxdd=dimRes1-1 ou dx_max
                    if(i<=dx_max) {//ajout RL car migra[] était parcouru en dehors de ces bornes (dx_max)
						boundedC+=2*migra[i];
						if(i!=0) simulpars<<i<<" "<<2*migra[i]<<" "<<boundedC<<endl; else simulpars<<i<<" "<<migra[i]<<" "<<boundedC<<endl;
						meanInddmax+=2*i*migra[i];
						sigInddmax+=2*i*i*migra[i];
						kurtInddmax+=2*i*i*i*i*migra[i];
					} else {
						boundedC+=0;
						simulpars<<i<<" "<<0<<" "<<boundedC<<endl;
						meanInddmax+=0;
						sigInddmax+=0;
						kurtInddmax+=0;
					}
                }
            }

            kurtInddmax=kurtInddmax/(sigInddmax*sigInddmax)-3.;

      		simulpars<<"\naxial mean dispersal (within bounds): "<<meanInddmax<<endl;
      		simulpars<<"axial sigma2 (within bounds): "<<sigInddmax<<endl;
      		simulpars<<"axial kurtosis (within bounds): "<<kurtInddmax<<endl;
      		simulpars<<"\nIf kurtosis appears wrong there may be many reasons\n";
  	        simulpars<<"   It may be because the dispersal distribution is not\n";
            simulpars<<"   computed with enough precision."<<endl;
	 		{   float Ds2InDxmax;
                Ds2InDxmax=ploidy*TVpars[0].initialDens*sigInddmax;
                simulpars<<"\nD sigma2 (within bounds): "<<Ds2InDxmax<<endl;
                simulpars<<"Neighborhood  (within bounds): "<<(2*Ds2InDxmax)<<"(1D) or "<<(2*PI*Ds2InDxmax)<<" (2D)."<<endl;
                simulpars<<"expected slope (within bounds): "<<1./(2*Ds2InDxmax)<<"(1D) or "<<1./(2*PI*Ds2InDxmax)<<" (2D)."<<endl;
            }
            simulpars << endl;
            if (model_mig_zone=='S' || model_mig_zone=='g') {
                simulpars<<"axial sigma2 (no bounds): "<<sig2<<endl;
                simulpars<<"D_h sigma2 (no bounds): "<<Dhs2<<endl;
                simulpars<<"Neighborhood  (no bounds): "<<(2*Dhs2)<<"(1D) or "<<(2*PI*Dhs2)<<" (2D)."<<endl;
                simulpars<<"expected slope  (no bounds): "<<1./(2*Dhs2)<<"(1D) or "<<1./(2*PI*Dhs2)<<" (2D)."<<endl;
            }
        /** AFFICHAGE ECRAN RL viré car non necessaire et moche
  //NOTA  CHECKSUM differe de 1 pour modËle en Óle...
            if(rep==1 && locus==0) { // pcq on peut changer de distrib dans le temps=> repasser dans cette fonction plusieurs fois par rep et par locus
                cout << endl;
                cout<<"model_mig: "<<model_mig<<"; checksum: "<<boundedC<<";\n sigma2 within bounds: "<<sigInddmax<<";\n kurtosis within bounds: "<<kurtInddmax<<flush;
  	//for(i=0;i<=dx_max;i++) printf("migra(%d): %g \n",i,migra[i]);
  	//getchar();
            }**/

/*            ********* A crude migraine.txt file (imperfect, but better than nothing)*/
			if(migraineSettingsbool) { //only for IBD and OnePop and OnePopVarSize
                //only for sibngle marker simulations
                if(!(nMarker>1)) {
                    // A integrer pour ecrire correctement quand multimarkers
    //                _mu = mutRateVector[locus+1];
    //                model_mut = mutModelIDVector[locus+1];
    //                min_allele = minAlleleVector[locus+1];
    //                max_mutations = maxMutationsVector[locus+1];
    //                variable_mubool = varMutRateBoolVector[locus+1];
    //                Kmin = kMinVector[locus+1];
    //                Kmax = kMaxVector[locus+1];
    //                Kini = kIniVector[locus+1];
    //                MotifSize = motifSizeVector[locus+1];
    //                pSMM = pSMMVector[locus+1];
    //                var_geom_TPM = geomTPMVector[locus+1];
    //                var_geom_GSM = geomGSMVector[locus+1];

                    //for IBD first
                    if(dimRes1*dimRes2>4 && Mig>0.0000001) {
                        float loctwoNmu=2*ploidy*TVpars[0].initialDens*mutRateVector[1];
                        ofstream migrainetxt("migraine.txt",ios::out);
                        if (dimRes2>1) {
                              migrainetxt<<"DemographicModel=PlanarIBD"<<endl;
                        } else
                              migrainetxt<<"DemographicModel=LinearIBD"<<endl;
                        migrainetxt<<"habitatPars= 0.5 0.5 "<<dimRes1<<" "<<dimRes2<<" 0"<<endl;
                        if ( ! isnan(Dhs2) ) {
                            if (dimRes2/vide>1) { //2D
                              migrainetxt<<"testPoint=Nb="<<2*PI*Dhs2<<endl;
                            } else
                              migrainetxt<<"testPoint=Nb="<<2*Dhs2<<endl;
                        }
                        if ( ! isnan(Nhm) ) {
                            migrainetxt<<"testPoint=2Nm="<<2*Nhm<<endl;
                        }
                        migrainetxt<<"testPoint=2Nmu="<<loctwoNmu<<endl;
                        migrainetxt<<"testPoint=g="<<GeoG<<endl;
                        if(model_mut.compare("KAM")==0) migrainetxt<<"GivenK="<<kMaxVector[1] - kMinVector[1] +1<<endl;
                        migrainetxt<<"AxialDistanceBins="<<dimRes1<<endl;
                        migrainetxt<<"Jobmin=1"<<endl;
                        migrainetxt<<"jobmax="<<repet<<endl;
                        migrainetxt<<"outputType=Perlocus"<<endl;
                        migrainetxt<<"genepoprootName=are_"<<endl;
                        migrainetxt<<"PointNumber=512"<<endl;
                        migrainetxt<<"NRunsPerPoint=30"<<endl;
                        if((dimRes1*dimRes2)<26) {
                            migrainetxt<<"statistic=IS"<<endl;
                        } else {
                            migrainetxt<<"statistic=PAC"<<endl;
                        }
                        if((dimRes1*dimRes2)<100) {
                            migrainetxt<<"WriteSequence=Over,Append"<<endl;
                        } else {
                            migrainetxt<<"WriteSequence=Over"<<endl;
                        }
                        migrainetxt<<"//Writesequence=WriteRKrig"<<endl;
                        migrainetxt<<"//writesequence=readpoints,Append"<<endl;
                        migrainetxt<<"Lowerbound= "<<loctwoNmu/3.<<", "<<Nhm*2./3.<<",  0"<<endl;
                        migrainetxt<<"UpperBound= "<<loctwoNmu*3.<<",  "<<Nhm*6<<",  0.999"<<endl;
                        migrainetxt<<"NextBounds=PointsFromR"<<endl;
                        migrainetxt.close();
                    } else if(dimRes1*dimRes2==1 && TVpars.size()>1) { //one pop Var Size
                        ofstream migrainetxt("migraine.txt",ios::out);
                        migrainetxt<<"Jobmin=1"<<endl;
                        migrainetxt<<"jobmax="<<repet<<endl;
                        migrainetxt<<"outputType=Perlocus"<<endl;
                        migrainetxt<<"genepopRootName="<<fichier_genepop<<endl;
                        migrainetxt<<"//Writesequence=WriteRKrig"<<endl;
                        migrainetxt<<"//writesequence=readpoints,Append"<<endl;
                        if(dispMut[1]!="GSM" && dispMut[1]!="TPM") {
                            migrainetxt<<"WriteSequence=Over,Append"<<endl;
                            migrainetxt<<"PointNumber=600"<<endl;
                        } else {
                            migrainetxt<<"WriteSequence=Over,Append,Append"<<endl;
                            migrainetxt<<"PointNumber=2400,2400,1200"<<endl;
                        }
                        migrainetxt<<"NRunsPerPoint=2000"<<endl;
                        migrainetxt<<"statistic=IS"<<endl;
                        migrainetxt<<"DemographicModel=OnePopVarSize"<<endl;
                        if(TVpars[0].ContDemeSizeChange=="Exponential")
                            migrainetxt<<"VarSizeFunction=Exponential"<<endl;
                        else if(TVpars[0].ContDemeSizeChange=="Logistic")
                            migrainetxt<<"VarSizeFunction=Logistic"<<endl;
                        else if(TVpars[0].ContDemeSizeChange=="Linear")
                            migrainetxt<<"VarSizeFunction=Linear"<<endl;
                        else migrainetxt<<"VarSizeFunction=Discrete"<<endl;
                        if(dispMut[1]=="KAM") // code pas terrible...
                            migrainetxt<<"MutationModel=KAM"<<endl;
                        else if(dispMut[1]=="SMM")
                            migrainetxt<<"MutationModel=SMM"<<endl;
                        else if(dispMut[1]=="GSM" || dispMut[1]=="TPM" )
                            migrainetxt<<"MutationModel=GSM"<<endl;
                        if(dispMut[1]=="GSM" || dispMut[1]=="TPM" || dispMut[1]=="SMM")
                            migrainetxt<<"StepSizes="<< MotifSize <<endl;
                        migrainetxt<<"//GivenK="<<kMaxVector[1] - kMinVector[1] + 1<<endl;

                        migrainetxt<<"TestPoint=2Nmu="<<(float) 2.*ploidy*TVpars[0].initialDens*mutRateVector[1]<<endl;
                        migrainetxt<<"TestPoint=Dchange="<<(float) 1.*TVpars[1].minGen/(2.*ploidy*TVpars[0].initialDens)<<endl;
                        migrainetxt<<"TestPoint=twoNancmu="<<(float) 2.*ploidy*TVpars[1].initialDens*mutRateVector[1]<<endl;
                        migrainetxt<<"TestPoint=Nratio="<<(float) 1.*TVpars[0].initialDens/TVpars[1].initialDens<<endl;

                        migrainetxt<<"NextBounds=PointsFromR"<<endl;
                        migrainetxt<<"RArguments=AutomatedCleaning"<<endl;
                        migrainetxt<<"GridSteps=40"<<endl;
                        if(dispMut[1]=="GSM" || dispMut[1]=="TPM" )
                            migrainetxt<<"1DCI=pGSM"<<endl;
                        migrainetxt<<"1DCI=2Nmu,Dchange,2Nancmu,Nratio"<<endl;
                        migrainetxt<<"RArguments=AutomatedCleaning"<<endl;

                        if(dispMut[1]=="GSM" || dispMut[1]=="TPM") {
                            migrainetxt<<"Lowerbound=0.1," << (float) 2*ploidy*TVpars[0].initialDens*mutRateVector[1]/10.;
                            migrainetxt<<",0.0,"<<(float) 1.*TVpars[1].minGen/(2.*ploidy*TVpars[0].initialDens)/50.<<","<< (float) 2*ploidy*TVpars[0].initialDens*mutRateVector[1]/10.<<endl;
                            migrainetxt<<"Upperbound=0.8," << (float) 2.*ploidy*TVpars[0].initialDens*mutRateVector[1]*10.;
                            migrainetxt<<",0.0,"<< (float) TVpars[1].minGen/(2.*ploidy*TVpars[0].initialDens)*10.<<","<< (float) 2.*ploidy*TVpars[0].initialDens*mutRateVector[1]*5.<<endl;
                            migrainetxt<<"SamplingScale=,logscale,,logscale,logscale"<<endl;
                        } else {
                            migrainetxt<<"Lowerbound=" << (float) 2*ploidy*TVpars[0].initialDens*mutRateVector[1]/10.;
                            migrainetxt<<",0.0,"<<(float) 1.*TVpars[1].minGen/(2.*ploidy*TVpars[0].initialDens)/50.<<","<< (float) 2*ploidy*TVpars[0].initialDens*mutRateVector[1]/10.<<endl;
                            migrainetxt<<"Upperbound=" << (float) 2.*ploidy*TVpars[0].initialDens*mutRateVector[1]*10.;
                            migrainetxt<<",0.0,"<< max((float) TVpars[1].minGen/(2.*ploidy*TVpars[0].initialDens)*10.,3.0)<<","<< (float) 2.*ploidy*TVpars[0].initialDens*mutRateVector[1]*5.<<endl;
                            migrainetxt<<"SamplingScale=logscale,,logscale,logscale"<<endl;
                        }
                        migrainetxt<<"Extrascale=Nratio=logscale"<<endl;

                        migrainetxt.close();
                    }
                }
            }
        } // endif !simulparsWrittenBool
        if(DispSpatiallyHeterog) SetForwardDispersalDistributions_zone();
         else simulparsWrittenBool=true;
    }//fin if deja
} // end def forw

/*****************/
//RL 072014
// nouvelle fonction pour les caractéristiques de dispersion sur la "zone" (dispersion heterogène dans l'espace)
//defini les probabilite de migration en fonction du modele de dispersion choisi
//defini dx_max et dy_max = distances de dispersion maximale consideree
//accessoirement calcul sigma et kurtosis pour verification
//appelee dans new_disp()
void CTimeVaryingParams::SetForwardDispersalDistributions_zone() {
    vector<long double>SichelDispPars(4);
    float Dhs2=numeric_limits<float>::quiet_NaN();
    int pas=1,dx_max_zone=dimRes1-1,dx_min_zone,dy_max_zone,dy_min_zone;
    //printf("model mig zone:%c",model_mig_zone);
    //if (locus ==2) cout<<"deja="<<deja<<flush;

if(!Simple1DProductDispBool) {//RL072014 wascmp_nocase(twoD_disp,"1DProductWithoutm0")==0) {
    cout << "option 1DproductWithoutm0 is not implemented with spatially heterogeneous dispersal" << endl;
    cout << "You must change your settings file " << SettingsFilename << endl;
    cout << "I exit" << endl;
    if(cinGetOnError) cin.get();
    exit(-1);
}

if(model_mig_zone=='a') {model_mig_zone='g';GeoG_zone=0;Mig_zone=1./3.;}
if(model_mig_zone=='b') {model_mig_zone='g';GeoG_zone=0;}
if(model_mig_zone=='1') {model_mig_zone='g';GeoG_zone=0;Mig_zone=2./3.;}


long double norm=0.0;
if (model_mig_zone=='P' || (model_mig_zone=='g' && GeoG_zone>0) || model_mig_zone=='S') {
    if(compareMath) {
        dx_max_zone=dist_max_zone;
    } else {
        dx_max_zone=dimRes1-1;
    }
} else {
    if (model_mig_zone=='g' && GeoG_zone==0) {
        dx_max_zone=1;
    } else if (model_mig_zone=='0') {/*pareto modifiee sigma≤=4 pour analyses a petite distance cf FR 2000 (c) (d)*/
        dx_max_zone=15; // parameters Mig & Pareto_Shape specified below because of specific mig[1] and mig[2] values
    } else if (model_mig_zone=='2') {/*kurtosis forte cf waser dipodomys sigma2=1*/
        dx_max_zone=49;model_mig_zone='P';Mig_zone=0.599985;Pareto_Shape_zone=3.79809435;
    } else if (model_mig_zone=='3') {/*sigma2=100*/
        dx_max_zone=48;model_mig_zone='P';Mig_zone=0.6;Pareto_Shape_zone=1.246085;
    } else if (model_mig_zone=='4') {/*sigma=1 pour un noeud sur deux*/
        cout << "Dispersal distribution 4, with void=2 is not implemented for simulations with spatially heterogeneous dispersal" << endl;
        cout << "I exit" << endl;
        if(cinGetOnError) cin.get();
        exit(-1);
    } else if (model_mig_zone=='5') {/*sigma=1 pour un noeud sur trois*/
        cout << "Dispersal distribution 5, with void=3 is not implemented for simulations with spatially heterogeneous dispersal" << endl;
        cout << "I exit" << endl;
        if(cinGetOnError) cin.get();
        exit(-1);
    } else if (model_mig_zone=='6') {/*sigma2=20*/
        dx_max_zone=48;model_mig_zone='P';Mig_zone=0.719326;Pareto_Shape_zone=2.0334337244;
    } else if (model_mig_zone=='7') {/*sigma2=10*/
        dx_max_zone=49;model_mig_zone='P';Mig_zone=0.702504;Pareto_Shape_zone=2.313010658;
    } else if (model_mig_zone=='8') {/*sigma≤=4 pour un noeud sur trois*/
        cout << "Dispersal distribution 8, with void=3 is not implemented for simulations with spatially heterogeneous dispersal" << endl;
        cout << "I exit" << endl;
        if(cinGetOnError) cin.get();
        exit(-1);
    } else if (model_mig_zone=='9') {/*"sigma2=4 mais sur 48 pas"*/
        dx_max_zone=48;model_mig_zone='P';Mig_zone=0.700013;Pareto_Shape_zone=2.74376568;
    }
    if( ! compareMath) dx_max_zone=dimRes1-1;
}
// setting dymax
if (dimRes2==1) dy_max_zone=0; else if(compareMath) {
    dy_max_zone=dx_max_zone;
} else {
    dy_max_zone=min(dimRes2-1,dx_max_zone);
}
dx_min_zone= - dx_max_zone;
dy_min_zone= - dy_max_zone;// RL : pas vraiment implémenté
if(dx_min_zone!=dx_min || dy_min_zone != dy_min) {
    cerr << "dx_min_zone!=dx_zone || dy_min_zone != dy_min is not implemented" << endl;
    cerr << "I must exit" << endl;
    if(cinGetOnError) cin.get();
    exit(1);
}
dist_max_zone=dx_max_zone;
if(dx_max_zone>dx_max) {
    cout << "dist_max_zone > dist_max (" << dx_max_zone << " vs " << dx_max << ")  but IBDSim can only consider dist_max_zone <= dist_max." << endl;
    cout << "You must change your settings file " << SettingsFilename << "to clearly specify dist_max_zone or use a build in function with smaller dx_max." << endl;
    cout << "I exit." << endl;
    if(cinGetOnError) cin.get();
    exit(-1);

}
migra_zone=ldvector(dx_max+1);
for(int i=0;i<=dx_max;i++) migra_zone[i]=0.0;

/************************************************************/
if( ! simulparsWrittenBool) {
    simulpars<<"\n*** Dispersal settings for the \"zone\" only ***"<<endl;
    simulpars<<"axial Dispersal proba (TotalEmigrationRate): "<<Mig<<endl;

    switch(model_mig_zone){
        case'g': {/*geometric*/
            if (GeoG_zone==0) {
                simulpars<<"Stepping stone model "<<endl ;
            } else if (GeoG_zone==1) {
                simulpars<<"Island model of dispersal "<<endl ;
            } else {
                simulpars<<"Geometric dispersal model "<<endl ;
                simulpars<<"GeometricShape: "<<GeoG_zone<<endl;
            }
        }
            break;
        case 'S': {
            if (SichelDispPars[2]<0) {
                simulpars<<"Inverse Gamma mixture (Chesson & Lee 2005) "<<endl ;
                if (SichelDispPars[3]>0) simulpars<<"also including a fraction "<<SichelDispPars[3]<<" of random dispersal [excluding focal patch!]: "<<endl;
                simulpars<<"gamma="<<SichelDispPars[0]<<endl ;
                simulpars<<"kappa (i.e. lim (omega xi))= "<<SichelDispPars[1]<<endl ;
            } else {
                simulpars<<"Sichel dispersal distribution model (Chesson & Lee 2005)"<<endl;
                simulpars<<"gamma="<<SichelDispPars[0]<<endl ;
                simulpars<<"xi="<<SichelDispPars[1]<<endl ;
                simulpars<<"omega="<<SichelDispPars[2]<<endl ;
            }
        }
            break;
        case 'P': {
            simulpars<<"Pareto model "<<endl ;
            simulpars<<"ParetoShape: "<<Pareto_Shape_zone<<endl;
        }
            break;
    }
}
/** @@@@@@@@@@@@@@@@@@ SWITCH :**/
switch(model_mig_zone){
    case'g': {/*geometric*/
            migra_zone[0]=(double) (1.0-Mig_zone);
            if (GeoG_zone==1) { //disp en ile ATTENTION PAS CONTINUITE ENTRE LES 2 CAS
                //for(int i=1;i<=dx_max;i++) migra[i]=(double) Mig/(dx_max-1); //aucune perte sur bord !!
                migra_zone[0]+=(double) Mig_zone/(dx_max_zone);//RL 122012 : pour avoir un vrai model en ile
                for(int i=1;i<=dx_max;i++) if(i<=dx_max_zone) migra_zone[i]=(double) Mig_zone/(2*dx_max_zone); else migra_zone[i]=0.0; //aucune perte sur bord !! modifs RL 15062011, voir avec FR   RL->FR
            } else {
                for(int i=1;i<=dx_max;i++) if(i<=dx_max_zone) norm+=(double) pow(GeoG_zone,i-1.0);
                for(int i=1;i<=dx_max;i++)
                    if(i<=dx_max_zone) migra_zone[i]=(double) pow(GeoG_zone,i-1.0)*Mig_zone/(2.0*norm); else migra_zone[i]=0.0;//so migra[<>0] sum to 1-migra[0]=mig. "Cross migration pattern" for low mig
            }
    }
        break;
    case'0': {/*pareto sigma≤=4 pour analyses a petite distance cf FR (c) (d)*/
        //migra=ldvector(dx_max+1);
        //AllocNotFixDisp(dx_max);
        Mig_zone=0.3;
        Pareto_Shape_zone=2.51829;
        migra_zone[0]=(double)(1-Mig_zone);
        migra_zone[1]=(double) 0.06;
        migra_zone[2]=(double) 0.03;
        for(int i=3;i<=dx_max;i++) if(i<=dx_max_zone) norm+=((double) (pow(1.0/i,Pareto_Shape_zone)));
        for(int i=3;i<=dx_max;i++)
            if(i<=dx_max_zone) migra_zone[i]=((double) ((Mig_zone-2*migra_zone[1]-2*migra_zone[2])*pow(1.0/i,Pareto_Shape_zone)/(2*norm))); else migra_zone[i]=0.0;
    }
        break;
    case'P': {/*Pareto standard*/
        migra_zone[0]=(double)(1.0-Mig_zone);
        for(int i=pas;i<=dx_max;i+=pas) if(i<=dx_max_zone) norm+=pow(1.0/(i/pas),Pareto_Shape_zone);
        for(int i=pas;i<=dx_max;i+=pas) if(i<=dx_max_zone) migra_zone[i]=((double) (Mig_zone*pow(1.0/(i/pas),Pareto_Shape_zone)/(2*norm))); else migra_zone[i]=0.0;
    }
        break;
    case'S': {/*Sichel and inverse gamma mixtures*/
        //cout<<"apres AllocNotFixDisp"<<flush;
        vector<long double>costable;
        unsigned int lPSON;
        long double temp;
        //SichelDispPars[0]: Sichel gamma
        SichelDispPars[0]=Sichel_Gamma_zone;
        if (SichelDispPars[0]>0) {cout<<"gamma arg. for Sichel distribution must be negative";getchar();exit(-1);}
        //SichelDispPars[1]: Sichel xi (or reciprocal gamma's kappa\sim xi omega)
        SichelDispPars[1]=Sichel_Xi_zone;
        //SichelDispPars[2]: Sichel omega (or something negative for reciprocal gamma)
        SichelDispPars[2]=Sichel_Omega_zone;
        /*SichelDispPars[3]: island (total) dispersal rate: with proba SichelDispPars[3]
         the propagule disperses somewhere (staying excluded)
         with proba 1-SichelDispPars[3] the propagule follows the Sichel dist.*/
        SichelDispPars[3]=0.; //non charfunc fraction
        //what follows somehow assumes dim_reseau1=dim_reseau2 ?
        if (EdgeEffect=="circ") {lPSON=dimRes1;} /*ca donne des transformÈes sur un habitat fermÈ
                                                  => Sicheltable[->lPSON]->Sicheltable[1]... pas terrible quand on veut calculer var et kurtosis;
                                                  ou bien il faut que lPSON>2 dxmax...*/
        else {lPSON=100*dx_max;} // pour Èviter effet repliement
        long double DeuxPiSurNx=2.*PI/lPSON;
        costable.resize(lPSON);
        for(long int ii=0;ii<lPSON;ii++) costable[ii]=cos(DeuxPiSurNx*ii);
        Sicheltable_zone=ldvector(lPSON/2+1);
        for(long int ii=0;ii<lPSON/2+1;ii++) {  //la valeur en lPSON/2+1 est utilisÈe
            Sicheltable_zone[ii]=Sichelcharfunc((long double)(2.*PI*ii/lPSON),SichelDispPars);
        }
        for(int _dx=0;_dx<=dx_max;_dx++) {
            temp=1.0; // c'est charfn[0]
            for(long int ii=1;ii<lPSON/2;ii++) {
                //2 parce que les 2 termes sym de la transformÈe. mais c'est la valeur en |_dx|
                temp+=2*Sicheltable_zone[ii]*costable[(_dx*ii)%lPSON];
                //if ((_dx==3500) && (ii%500==0)) {cout<<ii<<" "<<temp<<endl;getchar();}
            }
            if(lPSON%2==0) temp+=Sicheltable_zone[lPSON/2]*costable[(_dx%2)*lPSON/2];
            temp/=lPSON;
            if(_dx<=dx_max_zone) migra_zone[_dx]=temp; else migra_zone[_dx]=0.0; // migration selon Sichel
        }
        // SichelDispPars[3] is a proba of axial uniform dispersal, ce n'est pas Mig !
        for(long int _dx=0;_dx<=dx_max;_dx++) {
            if(_dx<=dx_max_zone) {
                temp=migra_zone[_dx]*(1.-SichelDispPars[3]);
                if (_dx>0) temp+=SichelDispPars[3]/lPSON;
                norm+=(_dx==0?1.:2.)*temp;
                migra_zone[_dx]=temp;
            }
        } // for _dx ...dx_max
        for(int _dx=0;_dx<=dx_max;_dx++) if(_dx<=dx_max_zone) migra_zone[_dx]/=norm; else migra_zone[_dx]=0.0;
        // at this point (if there is no uniform dispersal) the axial dispersal terms (0 included) are still in proportion to those of a true Sichel distribution
    	free_ldvector(Sicheltable);
    } // end case S
        break;
        cerr<<"(!) From SetForwardDispersalDistributions_zone(): Unknown migration model in the \"zone\" (nicknamed '"<<model_mig_zone<<"'. I exit."<<endl;
        if (cinGetOnError) cin.get();
        exit(-1);
} // END OF SWITCH

/** SUITE ECRITURE SIMULPARS **/

if (! simulparsWrittenBool) {
    long double sig2;
    long double boundedC=0;
    if (model_mig_zone=='S') {
        //cout<<"non charfunc fraction: "<<SichelDispPars[3]<<endl;
        if (SichelDispPars[2]<0) {
            sig2=-(1.-SichelDispPars[3])*SichelDispPars[1]/(2.*(1.+SichelDispPars[0]));
        } else {
            sig2=-(1.-SichelDispPars[3])*SichelDispPars[1]*bessel_k(1.+SichelDispPars[0],SichelDispPars[2],1)
            /(2.*bessel_k(SichelDispPars[0],SichelDispPars[2],1));
        }
        Dhs2=ploidy*TVpars[0].dens_zone*sig2;
    } else if (model_mig_zone=='g') {
            sig2=Mig_zone*(1.+GeoG_zone)/((1.-GeoG_zone)*(1.-GeoG_zone));
            Dhs2=ploidy*TVpars[0].dens_zone*sig2;
    }

    int maxdd=max(dimRes1,dimRes2)-1;
    long double sigInddmax=0;
    long double kurtInddmax=0;
    //pas 1D...withoutm0
    simulpars<<"Bounds (dx_max,dy_max)=("<<dx_max_zone<<","<<dy_max_zone<<")"<<endl;
    simulpars << "Marginal dispersal distribution and checksum:" << endl;
    boundedC= - migra_zone[0];
    int DistMax;
    if(compareMath) DistMax=dx_max; else DistMax=maxdd;
    for(int i=0;i<=DistMax;i++) { // up to maxdd=dimRes1-1 ou dx_max
        if(i<=dx_max) {//ajout RL car migra[] était parcouru en dehors de ces bornes (dx_max)
            boundedC+=2*migra_zone[i];
            if(i!=0) simulpars<<i<<" "<<2*migra_zone[i]<<" "<<boundedC<<endl; else simulpars<<i<<" "<<migra_zone[i]<<" "<<boundedC<<endl;
            sigInddmax+=2*i*i*migra_zone[i];
            kurtInddmax+=2*i*i*i*i*migra_zone[i];
        } else {
            boundedC+=0;
            simulpars<<i<<" "<<0<<" "<<boundedC<<endl;
            sigInddmax+=0;
            kurtInddmax+=0;
        }
    }

    kurtInddmax=kurtInddmax/(sigInddmax*sigInddmax)-3.;

    simulpars<<"\naxial sigma2 (within bounds): "<<sigInddmax<<endl;
    simulpars<<"axial kurtosis (within bounds): "<<kurtInddmax<<endl;
    simulpars<<"\nIf kurtosis appears wrong there may be many reasons\n";
    simulpars<<"   It may be because the dispersal distribution is not\n";
    simulpars<<"   computed with enough precision."<<endl;
    {   float Ds2InDxmax;
        Ds2InDxmax=ploidy*TVpars[0].dens_zone*sigInddmax;
        simulpars<<"\nD sigma2 (within bounds): "<<Ds2InDxmax<<endl;
        simulpars<<"Neighborhood  (within bounds): "<<(2*Ds2InDxmax)<<"(1D) or "<<(2*PI*Ds2InDxmax)<<" (2D)."<<endl;
        simulpars<<"expected slope (within bounds): "<<1./(2*Ds2InDxmax)<<"(1D) or "<<1./(2*PI*Ds2InDxmax)<<" (2D)."<<endl;
    }
    simulpars << endl;
    if (model_mig_zone=='S' || model_mig_zone=='g') {
        simulpars<<"axial sigma2 (no bounds): "<<sig2<<endl;
        simulpars<<"D_h sigma2 (no bounds): "<<Dhs2<<endl;
        simulpars<<"Neighborhood  (no bounds): "<<(2*Dhs2)<<"(1D) or "<<(2*PI*Dhs2)<<" (2D)."<<endl;
        simulpars<<"expected slope  (no bounds): "<<1./(2*Dhs2)<<"(1D) or "<<1./(2*PI*Dhs2)<<" (2D)."<<endl;
    }
    simulpars<<"\n**** End of dispersal settings for the \"zone\"****\n";
    simulparsWrittenBool=true;
} // endif !simulparsWrittenBool
} // end def forw_zone


/*****************/
//procedure principale de migration des noeuds restants


void set_new_disp(int locPhase) {
//	if ((rep==1 && locus==0 && currentGeneration==1) // calcule tjrs la distrib ‡ la premiËre rÈpÈt... sauf s constant disp
//        ||(currTVparsPtr->DispSpatiallyHeterog || currTVparsPtr->DensSpatiallyHeterog || !const_disp)) {//et aux autres repetition et moments de changement si on lui dit qu'elle n'est pas constante
    
    //si on est dans un changement temporel, libere les vecteurs avant de les reallouer, notament dans SetForwardDispersalDistributions(); et AllocNotFixDisp();
    if(locPhase > 0)
        if(TVpars[locPhase-1].DispSpatiallyHeterog || TVpars[locPhase-1].DensSpatiallyHeterog || TVpars[locPhase-1].forwardBarrier || !const_disp) {
                free_ldvector(migra);
                free_ldvector(cum_fixex);
                free_ldvector(cum_fixey);
                if(TVpars[locPhase-1].DispSpatiallyHeterog || TVpars[locPhase-1].DensSpatiallyHeterog || TVpars[locPhase-1].forwardBarrier) {// si on sort d'une phase avec heterogeneites spatiales
                    if(TVpars[locPhase-1].DispSpatiallyHeterog) free_ldvector(migra_zone);
#ifdef DEBUG
#else
                    for(int i=TVpars[locPhase-1].dimRes1+1; i>=0; i--) {
                        for(int j=TVpars[locPhase-1].dimRes2+1; j>=0; j--) {
                            //cout << "set_new_disp() free memory for tabdis[" << i << "][" << j << "].cum2" << ends;
                            for(int k=2*dx_max+1; k>=0; k--) {
                                delete[] tabdis[i][j].cum2[k];
                                //cout << "k=" << k << "; " << ends;
                            }
                            delete[] tabdis[i][j].cum2;
                        }
                        //cout << endl;
                        //cout << endl;
                    }
#endif
                    for(int i=TVpars[locPhase-1].dimRes1+1; i>=0; i--) delete[] tabdis[i];
                    delete[] tabdis;
                    

                }//end if(currTVparsPtr->DispSpatiallyHeterog || currTVparsPtr->DensSpatiallyHeterog)
    } //if( TVpars[locPhase-1].DispSpatiallyHeterog ||...
    
    //currTVparsPtr->SetForwardDispersalDistributions();// recalcule distr disp, deplacé dans reset_demo_params()
//	}//end if ((rep==1 && locus==0) || (!const_disp))


    /***** ancienne fonction
	model_mig=migra_tps(currentGeneration);//change les parametres demographiques si planifiÈ
    ******/
    reset_demo_params(locPhase); //compute specific parameters that are not direct variables of TVpars
	//if(!isinf(Gn1))


} // end set_new_disp


/*****************/
//calcul le denominateur des probabilitÈ cumulÈes de migration quand reseau non homogene (= fixe=0)
//appelÈe par int CTimeVaryingParams::computeCumuls() (avant : appelé par disp_new()) uniquement quand necessaire (= quand il y a un changement dÈmo)
/** RL->RL Total_Range_Dispersal=true, dispersion geomtrique, ABSORBING boundaries
    la distribution de distance de dispersion ne doit pas dépendre de dx_max, mais les tirages aleatoires en dependent
    EmpDisp semble OK : tirages differents dans la meme distribution, mais un code plus efficace
    donnerait exactement les memes tirages selon que Total_Range_Dispersal=true ou false. A revoir...
**/
int cumul1(int coordOrigineX, int coordOrigineY);
void cumul2(int coordOrigineX, int coordOrigineY);

int cumul1(int coordOrigineX, int coordOrigineY) {


int coordx,coordy;
double barrierCorrectionFactor;


 denom=0.0;
 for(int ddy=dy_min;ddy<=dy_max;ddy++) for(int ddx=dx_min;ddx<=dx_max;ddx++) {

     barrierCorrectionFactor=1.0;
     if(EdgeEffect!="abs") {
      
         if(currTVparsPtr->forwardBarrierX && ( (coordOrigineY >= currTVparsPtr->y1_barrier) && (coordOrigineY <= currTVparsPtr->y2_barrier) ) ) {
             if( (coordOrigineX < currTVparsPtr->x1_barrier) && ( (coordOrigineX+ddx) >= currTVparsPtr->x1_barrier) ) barrierCorrectionFactor *= currTVparsPtr->barrierCrossingRateUp;
             else if( (coordOrigineX >= currTVparsPtr->x1_barrier) && ( (coordOrigineX+ddx) < currTVparsPtr->x1_barrier) ) barrierCorrectionFactor *=  currTVparsPtr->barrierCrossingRateDown;
         }
         coordx=int(SortieSupfnPtr1(coordOrigineX+ddx,currTVparsPtr->dimRes1));
         coordx=int(SortieInffnPtr1(coordx,currTVparsPtr->dimRes1));
         if(EdgeEffect=="refl") {//pour reflection multiples
             if(currTVparsPtr->forwardBarrierX && ( (coordOrigineY >= currTVparsPtr->y1_barrier) && (coordOrigineY <= currTVparsPtr->y2_barrier) ) ) {
                 if( ( (coordOrigineX+ddx) < currTVparsPtr->x1_barrier) && ( coordx >= currTVparsPtr->x1_barrier) ) barrierCorrectionFactor *= currTVparsPtr->barrierCrossingRateUp;
                 else if( ( (coordOrigineX+ddx) >= currTVparsPtr->x1_barrier) && ( coordx < currTVparsPtr->x1_barrier) ) barrierCorrectionFactor *=  currTVparsPtr->barrierCrossingRateDown;
             }
             while((coordx<=0)||(coordx>currTVparsPtr->dimRes1)) {
                 coordx=int(SortieInffnPtr1(coordx,currTVparsPtr->dimRes1));
                 coordx=int(SortieSupfnPtr1(coordx,currTVparsPtr->dimRes1));
             }
         }
         
         if(currTVparsPtr->forwardBarrierY && ( (coordOrigineX >= currTVparsPtr->x1_barrier) && (coordOrigineX <= currTVparsPtr->x2_barrier) ) ) {
             if( (coordOrigineY < currTVparsPtr->y1_barrier) && ( (coordOrigineY+ddy) >= currTVparsPtr->y1_barrier) ) barrierCorrectionFactor *= currTVparsPtr->barrierCrossingRateUp;
             else if( (coordOrigineY >= currTVparsPtr->y1_barrier) && ( (coordOrigineY+ddy) < currTVparsPtr->y1_barrier) ) barrierCorrectionFactor *=  currTVparsPtr->barrierCrossingRateDown;
         }

         coordy=int(SortieSupfnPtr1(coordOrigineY+ddy,currTVparsPtr->dimRes2));
         coordy=int(SortieInffnPtr1(coordy,currTVparsPtr->dimRes2));
         if(EdgeEffect=="refl") {//pour reflection multiples
             if(currTVparsPtr->forwardBarrierY && ( (coordOrigineX >= currTVparsPtr->x1_barrier) && (coordOrigineX <= currTVparsPtr->x2_barrier) ) ) {
                 if( ( (coordOrigineY+ddy) < currTVparsPtr->y1_barrier) && ( coordy >= currTVparsPtr->y1_barrier) ) barrierCorrectionFactor *= currTVparsPtr->barrierCrossingRateUp;
                 else if( ( (coordOrigineY+ddy) >= currTVparsPtr->y1_barrier) && ( coordy < currTVparsPtr->y1_barrier) ) barrierCorrectionFactor *=  currTVparsPtr->barrierCrossingRateDown;
             }
             while((coordy<=0)||(coordy>currTVparsPtr->dimRes2)) {
                 coordy=int(SortieInffnPtr1(coordy,currTVparsPtr->dimRes2));
                 coordy=int(SortieSupfnPtr1(coordy,currTVparsPtr->dimRes2));
             }
         }
//      cout << endl;
//      cout << "cumul1->(coordOrigineX=" << coordOrigineX << "+ddx=" << ddx << ")=coordx=" << coordx << endl;
//      cout << "cumul1->(coordOrigineY=" << coordOrigineY << "+ddy=" << ddy << ")=coordy=" << coordy << endl;
//      cout << "cumul1->densite[" << coordx << "][" << coordy << "]=" << densite[coordx][coordy] << endl;

         if(currTVparsPtr->DispSpatiallyHeterog && ( (currTVparsPtr->xmin_zone<=coordx) && (coordx<=currTVparsPtr->xmax_zone)
                            && (currTVparsPtr->ymin_zone<=coordy) && (coordy<=currTVparsPtr->ymax_zone) ) ){
             denom+= ((double) 1.*barrierCorrectionFactor*densite[coordx][coordy]*migra_zone[abs(ddx)]*migra_zone[abs(ddy)]);
         } else denom+= ((double) 1.*barrierCorrectionFactor*densite[coordx][coordy]*migra[abs(ddx)]*migra[abs(ddy)]);
//      cout << "cumul1->migra[abs(ddx)=" << abs(ddx) << "]=" << migra[abs(ddx)] << endl;
//      cout << "cumul1->migra[abs(ddy)=" << abs(ddy) << "]=" << migra[abs(ddy)] << endl;
//      cout << "cumul1->denom+" << ((double) 1.*densite[coordx][coordy]*migra[abs(ddx)]*migra[abs(ddy)]) << "->" << denom << endl;
     } else {
        if(((coordOrigineX+ddx)>0)&&((coordOrigineX+ddx)<=currTVparsPtr->dimRes1)&&((coordOrigineY+ddy)>0)&&((coordOrigineY+ddy)<=currTVparsPtr->dimRes2)) { // if the new coordinates are inside the habitat

            if(currTVparsPtr->forwardBarrierX && ( (coordOrigineY >= currTVparsPtr->y1_barrier) && (coordOrigineY <= currTVparsPtr->y2_barrier) ) ) {
                if( (coordOrigineX < currTVparsPtr->x1_barrier) && ( (coordOrigineX+ddx) >= currTVparsPtr->x1_barrier) ) barrierCorrectionFactor *= currTVparsPtr->barrierCrossingRateUp;
                else if( (coordOrigineX >= currTVparsPtr->x1_barrier) && ( (coordOrigineX+ddx) < currTVparsPtr->x1_barrier) ) barrierCorrectionFactor *=  currTVparsPtr->barrierCrossingRateDown;
            }
            if(currTVparsPtr->forwardBarrierY && ( (coordOrigineX >= currTVparsPtr->x1_barrier) && (coordOrigineX <= currTVparsPtr->x2_barrier) ) ) {
                if( (coordOrigineY < currTVparsPtr->y1_barrier) && ( (coordOrigineY+ddy) >= currTVparsPtr->y1_barrier) ) barrierCorrectionFactor *= currTVparsPtr->barrierCrossingRateUp;
                else if( (coordOrigineY >= currTVparsPtr->y1_barrier) && ( (coordOrigineY+ddy) < currTVparsPtr->y1_barrier) ) barrierCorrectionFactor *=  currTVparsPtr->barrierCrossingRateDown;
            }
    //	     cout << endl;
    //	     cout << endl;
    //         cout << "cumul1->(coordOrigineX=" << coordOrigineX << "+ddx=" << ddx << ")=" << coordOrigineX+ddx << endl;
    //	     cout << "cumul1->(coordOrigineY=" << coordOrigineY << "+ddy=" << ddy << ")=" << coordOrigineY+ddy << endl;
    //	     cout << "cumul1->densite[][]=" << densite[coordOrigineX+ddx][coordOrigineY+ddy] << endl;
            if(currTVparsPtr->DispSpatiallyHeterog && ( (currTVparsPtr->xmin_zone<=(coordOrigineX+ddx)) && ((coordOrigineX+ddx)<=currTVparsPtr->xmax_zone)
                                        && (currTVparsPtr->ymin_zone<=(coordOrigineY+ddy)) && ((coordOrigineY+ddy)<=currTVparsPtr->ymax_zone) ) ){
                denom+= ((double) 1.*barrierCorrectionFactor*densite[coordOrigineX+ddx][coordOrigineY+ddy]*migra_zone[abs(ddx)]*migra_zone[abs(ddy)]);
    //                cout << "cumul1->migra_zone[abs(ddx)=" << abs(ddx) << "]=" << migra_zone[abs(ddx)] << endl;
    //                cout << "cumul1->migra_zone[abs(ddy)=" << abs(ddy) << "]=" << migra_zone[abs(ddy)] << endl;
    //                cout << "cumul1->denom+" << ((double) 1.*densite[coordOrigineX+ddx][coordOrigineY+ddy]*migra_zone[abs(ddx)]*migra_zone[abs(ddy)]) << "->" << denom << endl;
            } else {
                denom+= ((double) 1.*barrierCorrectionFactor*densite[coordOrigineX+ddx][coordOrigineY+ddy]*migra[abs(ddx)]*migra[abs(ddy)]);
    //                cout << "cumul1->migra[abs(ddx)=" << abs(ddx) << "]=" << migra[abs(ddx)] << endl;
    //                cout << "cumul1->migra[abs(ddy)=" << abs(ddy) << "]=" << migra[abs(ddy)] << endl;
    //                cout << "cumul1->denom+" << ((double) 1.*densite[coordOrigineX+ddx][coordOrigineY+ddy]*migra[abs(ddx)]*migra[abs(ddy)]) << "->" << denom << endl;
            }
        } // if inside the habitat
      } // else !ABS
    }//boucle ddx, ddy
if(isnan(denom)) {
   cout << "Probleme with the computation of tabdis denominator, denom is NaN." << endl;
   cout << "This should never happen, please contact the authors, and report this bug." << endl;
    cout << "I exit. Press any key to resume." << endl;
   if(cinGetOnError) cin.get();
}
return 0;
}/*fin cumul1*/

/*****************/
//calcul le numerateur et finalise les probabilitÈ cumulÈes de migration pour chaque noeud du reseau ou pour une "symetrie" quand reseau non homogene (= fixe=0)
//appelÈe par new_disp() uniquement quand c'est necessaire (= quand il y a un changement dÈmo)
void cumul2(int coordOrigineX, int coordOrigineY) {
    
double numer,frac;
int coordx,coordy;
double barrierCorrectionFactor;


 numer=0.0;
 for(int ddy=dy_min;ddy<=dy_max;ddy++) for(int ddx=dx_min;ddx<=dx_max;ddx++) {

     barrierCorrectionFactor=1.0;

     if(EdgeEffect!="abs") {
         if(currTVparsPtr->forwardBarrierX && ( (coordOrigineY >= currTVparsPtr->y1_barrier) && (coordOrigineY <= currTVparsPtr->y2_barrier) ) ) {
             if( (coordOrigineX < currTVparsPtr->x1_barrier) && ( (coordOrigineX+ddx) >= currTVparsPtr->x1_barrier) ) barrierCorrectionFactor *= currTVparsPtr->barrierCrossingRateUp;
             else if( (coordOrigineX >= currTVparsPtr->x1_barrier) && ( (coordOrigineX+ddx) < currTVparsPtr->x1_barrier) ) barrierCorrectionFactor *=  currTVparsPtr->barrierCrossingRateDown;
         }

         coordx=int(SortieSupfnPtr1(coordOrigineX+ddx,currTVparsPtr->dimRes1));
         coordx=int(SortieInffnPtr1(coordx,currTVparsPtr->dimRes1));
         if(EdgeEffect=="refl") {//pour reflection multiples
             if(currTVparsPtr->forwardBarrierX && ( (coordOrigineY >= currTVparsPtr->y1_barrier) && (coordOrigineY <= currTVparsPtr->y2_barrier) ) ) {
                 if( ( (coordOrigineX+ddx) < currTVparsPtr->x1_barrier) && ( coordx >= currTVparsPtr->x1_barrier) ) barrierCorrectionFactor *= currTVparsPtr->barrierCrossingRateUp;
                 else if( ( (coordOrigineX+ddx) >= currTVparsPtr->x1_barrier) && ( coordx < currTVparsPtr->x1_barrier) ) barrierCorrectionFactor *=  currTVparsPtr->barrierCrossingRateDown;
             }
             while((coordx<=0)||(coordx>currTVparsPtr->dimRes1)) {
                 coordx=int(SortieInffnPtr1(coordx,currTVparsPtr->dimRes1));
                 coordx=int(SortieSupfnPtr1(coordx,currTVparsPtr->dimRes1));
             }
         }
         
         if(currTVparsPtr->forwardBarrierY && ( (coordOrigineX >= currTVparsPtr->x1_barrier) && (coordOrigineX <= currTVparsPtr->x2_barrier) ) ) {
             if( (coordOrigineY < currTVparsPtr->y1_barrier) && ( (coordOrigineY+ddy) >= currTVparsPtr->y1_barrier) ) barrierCorrectionFactor *= currTVparsPtr->barrierCrossingRateUp;
             else if( (coordOrigineY >= currTVparsPtr->y1_barrier) && ( (coordOrigineY+ddy) < currTVparsPtr->y1_barrier) ) barrierCorrectionFactor *=  currTVparsPtr->barrierCrossingRateDown;
         }
         
         coordy=int(SortieSupfnPtr1(coordOrigineY+ddy,currTVparsPtr->dimRes2));
         coordy=int(SortieInffnPtr1(coordy,currTVparsPtr->dimRes2));
         if(EdgeEffect=="refl") {//pour reflection multiples
             if(currTVparsPtr->forwardBarrierY && ( (coordOrigineX >= currTVparsPtr->x1_barrier) && (coordOrigineX <= currTVparsPtr->x2_barrier) ) ) {
                 if( ( (coordOrigineY+ddy) < currTVparsPtr->y1_barrier) && ( coordy >= currTVparsPtr->y1_barrier) ) barrierCorrectionFactor *= currTVparsPtr->barrierCrossingRateUp;
                 else if( ( (coordOrigineY+ddy) >= currTVparsPtr->y1_barrier) && ( coordy < currTVparsPtr->y1_barrier) ) barrierCorrectionFactor *=  currTVparsPtr->barrierCrossingRateDown;
             }
             while((coordy<=0)||(coordy>currTVparsPtr->dimRes2)) {
                 coordy=int(SortieInffnPtr1(coordy,currTVparsPtr->dimRes2));
                 coordy=int(SortieSupfnPtr1(coordy,currTVparsPtr->dimRes2));
             }
         }
         if(currTVparsPtr->DispSpatiallyHeterog && ( (currTVparsPtr->xmin_zone<=coordx) && (coordx<=currTVparsPtr->xmax_zone)
                                     && (currTVparsPtr->ymin_zone<=coordy) && (coordy<=currTVparsPtr->ymax_zone) ) ){
             numer+= ((double) 1.*barrierCorrectionFactor*densite[coordx][coordy]*migra_zone[abs(ddx)]*migra_zone[abs(ddy)]);
         } else numer+= ((double) 1.*barrierCorrectionFactor*densite[coordx][coordy]*migra[abs(ddx)]*migra[abs(ddy)]);
//      cout << "cumul2->migra[abs(ddx)=" << abs(ddx) << "]=" << migra[abs(ddx)] << endl;
//      cout << "cumul2->migra[abs(ddy)=" << abs(ddy) << "]=" << migra[abs(ddy)] << endl;
//      cout << "cumul2->densite[" << coordx << "][" << coordy << "]=" << densite[coordx][coordy] << endl;
//      cout << "denom=" << denom << endl;
//      cout << "numer+" << ((double) 1.*densite[coordx][coordy]*migra[abs(ddx)]*migra[abs(ddy)]) << "->" << numer << endl;
         if(denom!=0.0) frac= ((double) numer/denom);
         else if(numer!=0.0) {printf("probleme : denom=0,numer=%f \n	the program must terminate",numer);getchar();exit(1);}
			 else frac=0.0;
//       cout << "frac=" << frac << endl;

         tabdis[coordOrigineX][coordOrigineY].cum2[ddx+dx_max][ddy+dy_max]=frac;
//      cout << "tabdis[coordOrigineX][coordOrigineY].cum2[ddx+dx_max][ddy+dy_max]=" << tabdis[coordOrigineX][coordOrigineY].cum2[ddx+dx_max][ddy+dy_max] << endl;
	 } else {
         //cout << endl;
         //cout << "cumul2->(coordOrigineX=" << coordOrigineX << "+ddx=" << ddx << ")=" << coordOrigineX+ddx << endl;
         //cout << "cumul2->(coordOrigineY=" << coordOrigineY << "+ddy=" << ddy << ")=" << coordOrigineY+ddy << endl;
         if( ((coordOrigineX+ddx)>0) &&((coordOrigineX+ddx)<=currTVparsPtr->dimRes1)&&((coordOrigineY+ddy)>0) && ((coordOrigineY+ddy)<=currTVparsPtr->dimRes2) ){
             
             if(currTVparsPtr->forwardBarrierX && ( (coordOrigineY >= currTVparsPtr->y1_barrier) && (coordOrigineY <= currTVparsPtr->y2_barrier) ) ) {
                 if( (coordOrigineX < currTVparsPtr->x1_barrier) && ( (coordOrigineX+ddx) >= currTVparsPtr->x1_barrier) ) barrierCorrectionFactor *= currTVparsPtr->barrierCrossingRateUp;
                 else if( (coordOrigineX >= currTVparsPtr->x1_barrier) && ( (coordOrigineX+ddx) < currTVparsPtr->x1_barrier) ) barrierCorrectionFactor *=  currTVparsPtr->barrierCrossingRateDown;
             }
             if(currTVparsPtr->forwardBarrierY && ( (coordOrigineX >= currTVparsPtr->x1_barrier) && (coordOrigineX <= currTVparsPtr->x2_barrier) ) ) {
                 if( (coordOrigineY < currTVparsPtr->y1_barrier) && ( (coordOrigineY+ddy) >= currTVparsPtr->y1_barrier) ) barrierCorrectionFactor *= currTVparsPtr->barrierCrossingRateUp;
                 else if( (coordOrigineY >= currTVparsPtr->y1_barrier) && ( (coordOrigineY+ddy) < currTVparsPtr->y1_barrier) ) barrierCorrectionFactor *=  currTVparsPtr->barrierCrossingRateDown;
             }

             if(currTVparsPtr->DispSpatiallyHeterog && ( (currTVparsPtr->xmin_zone<=(coordOrigineX+ddx)) && ((coordOrigineX+ddx)<=currTVparsPtr->xmax_zone)
                                         && (currTVparsPtr->ymin_zone<=(coordOrigineY+ddy)) && ((coordOrigineY+ddy)<=currTVparsPtr->ymax_zone) ) ){
                 numer+= ((double) 1.*barrierCorrectionFactor*densite[coordOrigineX+ddx][coordOrigineY+ddy]*migra_zone[abs(ddx)]*migra_zone[abs(ddy)]);
                 //cout << "cumul2->migra_zone[abs(ddx)=" << abs(ddx) << "]=" << migra_zone[abs(ddx)] << endl;
                 //cout << "cumul2->migra_zone[abs(ddy)=" << abs(ddy) << "]=" << migra_zone[abs(ddy)] << endl;
                 //cout << "cumul2->densite[][]=" << densite[coordOrigineX+ddx][coordOrigineY+ddy] << endl;
                 //cout << "denom=" << denom << endl;
                 //cout << "numer+" << ((double) 1.*densite[coordOrigineX+ddx][coordOrigineY+ddy]*migra_zone[abs(ddx)]*migra_zone[abs(ddy)]) << "->" << numer << endl;
             } else {
                 numer+= ((double) 1.*barrierCorrectionFactor*densite[coordOrigineX+ddx][coordOrigineY+ddy]*migra[abs(ddx)]*migra[abs(ddy)]);
                 //cout << "cumul2->migra[abs(ddx)=" << abs(ddx) << "]=" << migra[abs(ddx)] << endl;
                 //cout << "cumul2->migra[abs(ddy)=" << abs(ddy) << "]=" << migra[abs(ddy)] << endl;
                 //cout << "cumul2->densite[][]=" << densite[coordOrigineX+ddx][coordOrigineY+ddy] << endl;
                 //cout << "denom=" << denom << endl;
                 //cout << "numer+" << ((double) 1.*densite[coordOrigineX+ddx][coordOrigineY+ddy]*migra[abs(ddx)]*migra[abs(ddy)]) << "->" << numer << endl;
             }
             if(denom!=0.0) frac= ((double) numer/denom);
               else if(numer!=0.0) {printf("probleme : denom=0,numer=%f \n	the program must terminate",numer);getchar();exit(1);}
                     else frac=0.0;
                 //cout << "frac=" << frac << endl;
             tabdis[coordOrigineX][coordOrigineY].cum2[ddx+dx_max][ddy+dy_max]=frac;
         } else {
          if(ddx!=-dx_max) tabdis[coordOrigineX][coordOrigineY].cum2[ddx+dx_max][ddy+dy_max]=tabdis[coordOrigineX][coordOrigineY].cum2[ddx-1+dx_max][ddy+dy_max]; //RL 072014 was =0.0;
          else tabdis[coordOrigineX][coordOrigineY].cum2[ddx+dx_max][ddy+dy_max]=0.0;
         }
         //cout << "Filling tabdis[" << coordOrigineX << "][" << coordOrigineY << "].cum2[" << ddx+dx_max << "][" << ddy+dy_max<< "]=" << tabdis[coordOrigineX][coordOrigineY].cum2[ddx+dx_max][ddy+dy_max] << endl;
	  }
 }
if(isnan(numer)) {
    cout << "Probleme with the computation of tabdis numerator, numer is NaN." << endl;
    cout << "This should never happen, please contact the authors, and report this bug." << endl;
    cout << "I exit. Press any key to resume." << endl;
    if(cinGetOnError) cin.get();
}
}/*fin cumul2*/

/********************/
//for Specific_Density_Designbool : get density for each point of the lattice from file DensityMatrix.txt
int getDensityFromFile() {
    int pos;

ifstream inFile(SpecificDensityFilename.c_str(),ios::in);
if(!inFile.is_open()) {
	cerr << "Unable to open file "<<SpecificDensityFilename << " in read_settings_file();" << endl;
	cerr << "Abnormal termination...Press any key to close the window..." << endl;
	if(cinGetOnError) cin.get();
	exit(-1);
	}
else {
    string locstring;
	int i=1;
	do {
        if(i>TVpars[0].dimRes2) {
			cout << endl << endl << endl;
			cout << "Problem with column dimension in DensityMatrix.txt" << endl;
			cout << "There are " << i<< " rows in " << SpecificDensityFilename << endl;
			cout << "But Y dimension of the lattice is " << TVpars[0].dimRes2 << endl;
			cout << "Programm aborted" << endl;
			if(cinGetOnError) cin.get();
            exit(-1);
        }
        getline(inFile,locstring);
		if(locstring.length()==0) break;
		int j=1;
		do {
			if(j>TVpars[0].dimRes1) {
				  cout << endl << endl << endl;
                  cout << "Problem with row dimension in DensityMatrix.txt" << endl;
    			  cout << "There are " << j << " column in " << SpecificDensityFilename << endl;
    			  cout << "But X dimension of the lattice is " << TVpars[0].dimRes1 << endl;
    			  cout << "Programm aborted" << endl;
                if(cinGetOnError) cin.get();
                exit(-1);
            }
            while((locstring[0]==' ')||(locstring[0]=='\t')) {locstring.erase(0,1);}//vire les blancs
			pos=std::min(locstring.find(' '),std::min(locstring.find('\t'),locstring.length()));
			densite[j][i]=atoi(locstring.substr(0,pos).c_str());
			//cout << "densite[" << j << "][" << i << "]" << densite[j][i] << endl;
			//getchar();
			locstring.erase(0,pos);
			j++;
			}while((locstring.length()>0)&&(locstring[0]!=';')&&(locstring[0]!='\0'));

		if((j-1)!=TVpars[0].dimRes1) {
			cout << endl << endl << endl;
			cout << "Problem with row dimension in DensityMatrix.txt" << endl;
			cout << "There are " << j-1 << " column in " << SpecificDensityFilename << endl;
			cout << "But X dimension of the lattice is " << TVpars[0].dimRes1 << endl;
			cout << "Programm aborted" << endl;
			if(cinGetOnError) cin.get();
            exit(-1);
			}
		i++;
		}while(!inFile.eof());

		inFile.close();

	    if((i-1)!=TVpars[0].dimRes2) {
		    cout << endl << endl << endl;
			cout << "Problem with column dimension in DensityMatrix.txt" << endl;
			cout << "There are " << i-1 << " rows in " << SpecificDensityFilename << endl;
			cout << "But Y dimension of the lattice is " << TVpars[0].dimRes2 << endl;
			cout << "Programm aborted" << endl;
			if(cinGetOnError) cin.get();
            exit(-1);
		}

}//else = file open

if(Specific_Sample_Designbool[0]) {
    for(int i=0;i<Spec_SampleSize[0];i++) if(densite[Spec_Sample_Coord[0][x][i]][Spec_Sample_Coord[0][y][i]]==0) {
        cout << endl << endl << endl;
		cout<<"Deme (x=" << Spec_Sample_Coord[0][x][i] << ")(y=" << Spec_Sample_Coord[0][y][i] << ") is present in Specific_Sample_Designbool"<<endl;
	    cout<<"But density at this point is " << densite[Spec_Sample_Coord[0][x][i]][Spec_Sample_Coord[0][y][i]] << ";"<<endl;
	    cout<<"Be carrefull that both specific designs are compatible.\nProgram aborted"<<endl;
	    if(cinGetOnError) cin.get();
        exit(-1);
    }
}
if(Specific_Sample_Designbool[1]) {
    for(int i=0;i<Spec_SampleSize[1];i++) if(densite[Spec_Sample_Coord[1][x][i]][Spec_Sample_Coord[1][y][i]]==0) {
        cout << endl << endl << endl;
        cout<<"Deme (x=" << Spec_Sample_Coord[1][x][i] << ")(y=" << Spec_Sample_Coord[1][y][i] << ") is present in Predisp_Specific_Sample_Coord"<<endl;
        cout<<"But density at this point is " << densite[Spec_Sample_Coord[1][x][i]][Spec_Sample_Coord[1][y][i]] << ";"<<endl;
        cout<<"Be carrefull that both specific designs are compatible.\nProgram aborted"<<endl;
        if(cinGetOnError) cin.get();
        exit(-1);
    }
}
return(0);
}
/********************/

/*****************/
//calcul dispxx et dispyy, les distance de migration quand le reseau est non homogene (= fixe=0) sinon c'est dispx()
void dispxy()
{double alea1;
 int ddx,ddy,dy=0,dx=0;

Disp1:

alea1=alea();
ddy=0;
//cout << endl << endl <<  "Disp from coord_origine[x]="  << coord_origine[x] << "; coord_origine[y]="  << coord_origine[y] << endl;
//for(int j=-dy_max;j<=dy_max;j++) for(int i=-dx_max;i<=dx_max;i++) {
//    //if(isnan(tabdis[coord_origine[x]][coord_origine[y]].cum2[i+dx_max][j+dy_max]))
////    if(tabdis[coord_origine[x]][coord_origine[y]].cum2[i+dx_max][j+dy_max]>0.0
////       && tabdis[coord_origine[x]][coord_origine[y]].cum2[i+dx_max][j+dy_max]<1.0)
//        printf("\t tabdis[%d][%d].cum2[%d][%d]= %f\n",coord_origine[x],coord_origine[y]
//              ,i,j,tabdis[coord_origine[x]][coord_origine[y]].cum2[i+dx_max][j+dy_max]);
////        cin.get();
//}
fflush(stdout);

// si on est en une dimension, alors dy=0 et on bouge que sur x
if(dy_max==0) {dy=0;goto fin1;}


if((alea1<tabdis[coord_origine[x]][coord_origine[y]]
	  .cum2[2*dx_max][ddy+dy_max])
   && (alea1>tabdis[coord_origine[x]][coord_origine[y]]
       .cum2[2*dx_max][ddy-1+dy_max]) )
		 {dy=0;goto fin1;}
if(alea1>tabdis[coord_origine[x]][coord_origine[y]]
	  .cum2[2*dx_max][ddy+dy_max]) goto fin1plus;
if(alea1<tabdis[coord_origine[x]][coord_origine[y]]
	  .cum2[2*dx_max][ddy-1+dy_max]) goto fin1moins;
fin1moins:
if(alea1<tabdis[coord_origine[x]][coord_origine[y]]
	  .cum2[2*dx_max][0]) {dy=dy_min;goto fin1;}
for(ddy=-1;ddy>dy_min;ddy--)
{
//printf("tabdis[%d][%d].cum2[%d][%d]= %f\n",coord_origine[x],coord_origine[y]
//  ,dx_max,ddy,tabdis[coord_origine[x]][coord_origine[y]].cum2[2*dx_max][ddy+dy_max]);
 //getchar();
 if(alea1>tabdis[coord_origine[x]][coord_origine[y]]
	  .cum2[2*dx_max][ddy-1+dy_max]) {dy=ddy;goto fin1;}
}
fin1plus:
for(ddy=1;ddy<=dy_max;ddy++)
{
//    printf("tabdis[%d][%d].cum2[%d][%d]= %f\n",coord_origine[x],coord_origine[y]
//			  ,dx_max,ddy,tabdis[coord_origine[x]][coord_origine[y]].cum2[2*dx_max][ddy+dy_max]);
 //getchar();
 if(alea1<tabdis[coord_origine[x]][coord_origine[y]]
	  .cum2[2*dx_max][ddy+dy_max]) {dy=ddy;goto fin1;}
}

fin1:
ddx=0;
if((alea1<tabdis[coord_origine[x]]
		  [coord_origine[y]].cum2[ddx+dx_max][dy+dy_max])
  &&(alea1>tabdis[coord_origine[x]]
		  [coord_origine[y]].cum2[ddx-1+dx_max][dy+dy_max]))
		 {dx=0;goto fin2;}
if(alea1>tabdis[coord_origine[x]][coord_origine[y]]
			.cum2[ddx+dx_max][dy+dy_max]) goto fin2plus;
if(alea1<tabdis[coord_origine[x]][coord_origine[y]]
			.cum2[ddx-1+dx_max][dy+dy_max]) goto fin2moins;
fin2moins:
if(alea1<tabdis[coord_origine[x]][coord_origine[y]]
			 .cum2[0][dy+dy_max]) {dx=dx_min;goto fin2;}

for(ddx=-1;ddx>dx_min;ddx--)
{
//    printf("tabdis[%d][%d].cum2[%d][%d]= %f\n",coord_origine[x],coord_origine[y]
//			  ,ddx-1,ddy,tabdis[coord_origine[x]][coord_origine[y]].cum2[2*dx_max][ddy+dy_max]);
 //getchar();
 if(alea1>tabdis[coord_origine[x]][coord_origine[y]]
			 .cum2[ddx-1+dx_max][dy+dy_max]) {dx=ddx;goto fin2;}
}
fin2plus:
for(ddx=1;ddx<=dx_max;ddx++)
{
//    printf("tabdis[%d][%d].cum2[%d][%d]= %f\n",coord_origine[x],coord_origine[y]
//			  ,dx_max,ddy,tabdis[coord_origine[x]][coord_origine[y]].cum2[2*dx_max][ddy+dy_max]);
 //getchar();
 if(alea1<tabdis[coord_origine[x]][coord_origine[y]]
			 .cum2[ddx+dx_max][dy+dy_max]) {dx=ddx;goto fin2;}
}
fin2:


// RL 072014, a virer a terme
if(EdgeEffect=="abs" && ( ((coord_origine[x]+dx)>currTVparsPtr->dimRes1) || ((coord_origine[x]+dx)<=0)
                         || ((coord_origine[y]+dy)>currTVparsPtr->dimRes2) || ((coord_origine[y]+dy)<=0)) ) {
    cerr << "Dispersal occurred outside the lattice with absorbing boundaries." << endl;
    cerr << "CurrentGeneration=" << currentGeneration << endl;
    cerr << "This should never happen...There is probably a bug. I exit." << endl;
    if(cinGetOnError) cin.get();
    exit(-1);
    goto Disp1;
}

coord_origine[x]+=dx;
coord_origine[y]+=dy;

    if(false && (dx>0 || dy>0)) {// && dx==-5) {
        cout << endl << endl <<  "Dispxy() from coord_origine[x]="  << coord_origine[x]-dx << "; coord_origine[y]="  << coord_origine[y]-dy << endl;
        printf("nbre noeuds restants %d aleat %f",NS_coal_tree::nbre_noeud_restant,alea1);
        printf(";dx=  %d",dx);
        printf(";dy=  %d",dy);
        cout << "New coord[x]="  << coord_origine[x] << "; coord[y]="  << coord_origine[y] << endl;
        //    getchar();
    }


} /*fin dispxy*/


/*****************/
//procedure principale de migration quand le reseau est homogene
//appellee par disp_new() si fixe!=0
void disp_fixe() {
using namespace NS_coal_tree;
using namespace NS_translation;
    using namespace NS_diagnostic_tables;

	int i,dispx=0,dispy=0;
	int originex,originey,nbre_noeud_restant_loc;
	double aleat1;
    double migra0=1.0-currTVparsPtr->Mig;  //CONSTANT local variable compared to random deviate many times

//	static int cptUp=0,cptDown=0,cpt1Up=0,cpt1Down=0;

	if (currTVparsPtr==NULL) {
	    cerr<<"(!!!) From disp_fixe: currTVparsPtr not valued ";
	    if (cinGetOnError) cin.get();
	    exit(-1);
    }

//lance le calcul des probabilitÈ cumulÈs de migration si ils ne sont pas deja calculÈ
if(deja1==0) {
 for(i=dx_min;i<=dx_max;i++) cumul_fixex(i);
 if(currTVparsPtr->dimRes2/currTVparsPtr->vide>1) //RL072014 was && (cmp_nocase(twoD_disp,"Simple1DProduct")==0 || cmp_nocase(twoD_disp,"1DProductWithoutm0")==0))
      for(i=dy_min;i<=dy_max;i++) cumul_fixey(i);
 deja1=1;
}

//coord_origine=ivector(2);
//migration de chaque noeud restant
if(currentGeneration==1) nbre_noeud_restant_loc=n_genes[0]; else nbre_noeud_restant_loc=nbre_noeud_restant;
for (i=1;i<=nbre_noeud_restant_loc;i++) {
    //cout << endl << "noeud " << i << ":" <<  coord_noeud[no_noeud_restant[i]][x] << " , " << coord_noeud[no_noeud_restant[i]][y] << endl;

	 if(EdgeEffect=="abs") {//si bord dispersifs/absorbants avec masse de probabilité "inateignable" reportée sur toute la distribution de dispersion
						// = on migre (ou non) jusqu'a ce que l'on se retrouve dans le reseau

		//determination de l'origine, translation si necessaire
		originex=coord_origine[x]=int(SortieSupfnPtr2(coord_noeud[no_noeud_restant[i]][x],currTVparsPtr->dimRes1,translationx));// RL -> RL modifiable : pour gagner du temps,
		originey=coord_origine[y]=int(SortieSupfnPtr2(coord_noeud[no_noeud_restant[i]][y],currTVparsPtr->dimRes2,translationy));//  ne le faire qu'en cas de changement de reseau ??vérifier??

		//dispersion
		if(Simple1DProductDispBool) {//RL072014 wascmp_nocase(twoD_disp,"Simple1DProduct")==0) { // default value
			//migration jsuq'a ce que l'on soit dans la pop pour bords dispersifs/absorbant
			do {coord_origine[x]=originex; coord_origine[x]+=dispx_fixe();} while((coord_origine[x]>currTVparsPtr->dimRes1)||(coord_origine[x]<=0));
			if(currTVparsPtr->backwardBarrierX && ( (originey >= currTVparsPtr->y1_barrier) && (originey <= currTVparsPtr->y2_barrier) ) ) {
				double alea2=alea();
//                if((originex < currTVparsPtr->x1_barrier) && (coord_origine[x] >= currTVparsPtr->x1_barrier)) cptUp++;
//                if((originex >= currTVparsPtr->x1_barrier) && (coord_origine[x] < currTVparsPtr->x1_barrier)) cptDown++;
				if( ((originex < currTVparsPtr->x1_barrier) && (coord_origine[x] >= currTVparsPtr->x1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateUp) ) {
//                    if((originex < currTVparsPtr->x1_barrier) && (coord_origine[x] >= currTVparsPtr->x1_barrier)) cpt1Up++;
						coord_origine[x]=originex; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
                } else if( ((originex >= currTVparsPtr->x1_barrier) && (coord_origine[x] < currTVparsPtr->x1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateDown) ) {
//                    if((originex >= currTVparsPtr->x1_barrier) && (coord_origine[x] < currTVparsPtr->x1_barrier)) cpt1Down++;
                    coord_origine[x]=originex; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
                }
			}

				  /*if(dimRes2/currTVparsPtr->vide!=1)*/
			do {coord_origine[y]=originey; coord_origine[y]+=dispy_fixe();} while((coord_origine[y]>currTVparsPtr->dimRes2)||(coord_origine[y]<=0));
			if(currTVparsPtr->backwardBarrierY &&  ( (originex >= currTVparsPtr->x1_barrier) && (originex <= currTVparsPtr->x2_barrier) ) ) {
				double alea2=alea();
//                if((originey < currTVparsPtr->y1_barrier) && (coord_origine[y] >= currTVparsPtr->y1_barrier)) cptUp++;
//                if((originey >= currTVparsPtr->y1_barrier) && (coord_origine[y] < currTVparsPtr->y1_barrier)) cptDown++;
                if( ((originey < currTVparsPtr->y1_barrier) && (coord_origine[y] >= currTVparsPtr->y1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateUp) ) {
//                    if((originey < currTVparsPtr->y1_barrier) && (coord_origine[y] >= currTVparsPtr->y1_barrier)) cpt1Up++;
                    coord_origine[y]=originey; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
                } else if( ((originey >= currTVparsPtr->y1_barrier) && (coord_origine[y] < currTVparsPtr->y1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateDown) ) {
//                    if((originey >= currTVparsPtr->y1_barrier) && (coord_origine[y] < currTVparsPtr->y1_barrier)) cpt1Down++;
                    coord_origine[y]=originey; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
                }
			}

		} else if(!Simple1DProductDispBool) {//RL072014 was cmp_nocase(twoD_disp,"1DProductWithoutm0")==0) {
			  aleat1=alea();
			  if(aleat1>backwardPhilo[originex-1][originey-1]) { //FR0610 pas sur que originex-1 et originey-1 soient les bonnes coord ?
				  //cout << "mig" << endl;
				  //migration jsuq'a ce que l'on soit dans la pop pour bords dispersifs/absorbant
				  do {
					  if(currTVparsPtr->dimRes2/currTVparsPtr->vide>1/*if true 2D*/) do {dispx=dispx_fixe();dispy=dispy_fixe();} while(dispx==0 && dispy==0);
					  else do {dispx=dispx_fixe();} while(dispx==0);
					  //cout << "dispx=" << dispx << "dispy=" << dispy << endl;
					  //cin.get();
					  coord_origine[x]=originex;
					  coord_origine[x]+=dispx;
				      coord_origine[y]=originey;
					  coord_origine[y]+=dispy;
				  } while((coord_origine[x]>currTVparsPtr->dimRes1)||(coord_origine[x]<=0)||(coord_origine[y]>currTVparsPtr->dimRes2)||(coord_origine[y]<=0));
				  
                  if(currTVparsPtr->backwardBarrierX && ( (originey >= currTVparsPtr->y1_barrier) && (originey <= currTVparsPtr->y2_barrier) ) ) {
					  double alea2=alea();
					  if( ((originex < currTVparsPtr->x1_barrier) && (coord_origine[x] >= currTVparsPtr->x1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateUp) )
                          coord_origine[x]=originex; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
                      else if( ((originex >= currTVparsPtr->x1_barrier) && (coord_origine[x] < currTVparsPtr->x1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateDown) )
                          coord_origine[x]=originex; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
				  }
				  
                  if(currTVparsPtr->backwardBarrierY &&  ( (originex >= currTVparsPtr->x1_barrier) && (originex <= currTVparsPtr->x2_barrier) ) ){
					  double alea2=alea();
					  if( ((originey < currTVparsPtr->y1_barrier) && (coord_origine[y] >= currTVparsPtr->y1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateUp) )
                          coord_origine[y]=originey; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
                      else if( ((originey >= currTVparsPtr->y1_barrier) && (coord_origine[y] < currTVparsPtr->y1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateDown) )
                          coord_origine[y]=originey; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
				  }
			  } //else cout << "no mig;" << endl;
		} /*else if(cmp_nocase(twoD_disp,"X+YDistProb")==0) {
			do {disp=abs(dispx_fixe());
			    if(disp>0) {
			      if(dimRes2/currTVparsPtr->vide>1) dispx= (int) floor(alea()*(disp+1));
				  else dispx=disp;
					if(dispx>0) {
			            if(alea()>0.5) dispx=-dispx;
				        coord_origine[x]=originex;
				        coord_origine[x]+=dispx;
			        }
				    dispy=disp-abs(dispx);
			        if(dispy>0) {
					     if(alea()>0.5) dispy=-dispy;
						 coord_origine[y]=originey;
				         coord_origine[y]+=dispy;
			        }
				}
			} while((coord_origine[y]+dispy>dimRes2)||(coord_origine[y]+dispy<=0)||(coord_origine[x]+dispx>currTVparsPtr->dimRes1)||(coord_origine[x]+dispx<=0));
		}*/
		coord_ori[no_noeud_restant[i]][x]=coord_origine[x];
		coord_ori[no_noeud_restant[i]][y]=coord_origine[y];
	 } else if(EdgeEffect=="absFR") {
	     cerr<<"obsolete option. Look for versions pre Avril 2011 or pre 2.0 for orginal code"<<endl;
	} else {//si torre ou reflectifs//pour remettre les individus a la bonne place sur le reseau quand reduction de surface
     originex=int(SortieSupfnPtr2(coord_noeud[no_noeud_restant[i]][x],currTVparsPtr->dimRes1,translationx));
     coord_origine[x]=int(SortieSupfnPtr2(coord_noeud[no_noeud_restant[i]][x],currTVparsPtr->dimRes1,translationx));
     originey=int(SortieSupfnPtr2(coord_noeud[no_noeud_restant[i]][y],currTVparsPtr->dimRes2,translationy));
     coord_origine[y]=int(SortieSupfnPtr2(coord_noeud[no_noeud_restant[i]][y],currTVparsPtr->dimRes2,translationy));
     /*if(currentGeneration==1 || currentGeneration==Gn1 || currentGeneration==Gn2)
	   	{printf("\n at genaration %d coord_x=%d ;  coord_y=%d ",Gn,originex,originey);
	   	 getchar();
	   	}*/
     if(originex>currTVparsPtr->dimRes1 || originey>currTVparsPtr->dimRes2)
	   	{printf("\n problem : at generation %ld originex=%d ;  originey=%d ",currentGeneration,originex,originey);
	   	 printf("\n coord_noeud_x=%d ;  coord_noeud_y=%d ",coord_noeud[no_noeud_restant[i]][x],coord_noeud[no_noeud_restant[i]][y]);
	   	 fflush(stdout);
	   	 getchar();
		}
	if(Simple1DProductDispBool) {//RL072014 wascmp_nocase(twoD_disp,"Simple1DProduct")==0) {
		coord_origine[x]+=dispx_fixe();
		/*if(dimRes2/currTVparsPtr->vide!=1)*/ coord_origine[y]+=dispy_fixe();
	} else if(!Simple1DProductDispBool) {//RL072014 was cmp_nocase(twoD_disp,"1DProductWithoutm0")==0) { /** this is where this option makes a difference on dispersal**/
			aleat1=alea();
			if(aleat1>migra0) {
				if(currTVparsPtr->dimRes2/currTVparsPtr->vide>1) do {dispx=dispx_fixe();dispy=dispy_fixe();} while(dispx==0 && dispy==0);
				else do {dispx=dispx_fixe();} while(dispx==0);
				coord_origine[x]+=dispx;
				coord_origine[y]+=dispy;
			} //else cout << "no mig;" << endl;
	} /*else if(cmp_nocase(twoD_disp,"X+YDistProb")==0) {
			disp=abs(dispx_fixe());
			if(disp>0) {
				if(dimRes2/currTVparsPtr->vide>1) dispx=(int) floor(alea()*(disp+1));
				else dispx=disp;
				if(dispx>0) if(alea()>0.5) dispx=-dispx;

				dispy=disp-abs(dispx);
				if(dispy>0) if(alea()>0.5) dispy=-dispy;
				//cout << "dispx dispy = " << dispx << dispy << endl;
			}
			coord_origine[x]+=dispx;
			coord_origine[y]+=dispy;
	}*/
    
	  //Effets de bord reflectifs et circulaires pour x
     //if((coord_origine[x]<=0)||(coord_origine[x]>dim_reseau1)) printf("\n coord X avant=%d",coord_origine[x]);
     coord_ori[no_noeud_restant[i]][x]=int(SortieSupfnPtr1(coord_origine[x],currTVparsPtr->dimRes1));
     coord_ori[no_noeud_restant[i]][x]=int(SortieInffnPtr1(coord_ori[no_noeud_restant[i]][x],currTVparsPtr->dimRes1));
     if(EdgeEffect=="refl")//pour reflection multiples
       {while((coord_ori[no_noeud_restant[i]][x]<=0)||(coord_ori[no_noeud_restant[i]][x]>currTVparsPtr->dimRes1))
		 {coord_ori[no_noeud_restant[i]][x]=int(SortieInffnPtr1(coord_ori[no_noeud_restant[i]][x],currTVparsPtr->dimRes1));
          coord_ori[no_noeud_restant[i]][x]=int(SortieSupfnPtr1(coord_ori[no_noeud_restant[i]][x],currTVparsPtr->dimRes1));
         }
       }
     //
     //if((coord_origine[x]<=0)||(coord_origine[x]>dim_reseau1)) printf("coord X apres effets de bord =%d",coord_ori[no_noeud_restant[i]][x]);
     //printf("transfo donne=%d",(int) ((coord_ori[no_noeud_restant[i]][x])-(floor(((double) coord_ori[no_noeud_restant[i]][x]/(double) dim_reseau1)-(1E-5))*dim_reseau1)));
     //getchar();

        //if((coord_origine[y]<=0)||(coord_origine[y]>dim_reseau2)) {printf("coord Y apres effets bords=%d",coord_ori[no_noeud_restant[i]][y]);getchar();}
        if(currTVparsPtr->backwardBarrierX && ( (originey >= currTVparsPtr->y1_barrier) && (originey <= currTVparsPtr->y2_barrier) ) ) {
            double alea2=alea();
            //        if((originex < currTVparsPtr->x1_barrier) && (coord_ori[no_noeud_restant[i]][x] >= currTVparsPtr->x1_barrier)) cptUp++;
            //        if((originex >= currTVparsPtr->x1_barrier) && (coord_ori[no_noeud_restant[i]][x] < currTVparsPtr->x1_barrier)) cptDown++;
            if( ((originex < currTVparsPtr->x1_barrier) && (coord_ori[no_noeud_restant[i]][x] >= currTVparsPtr->x1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateUp) ) {
                //            if((originex < currTVparsPtr->x1_barrier) && (coord_ori[no_noeud_restant[i]][x] >= currTVparsPtr->x1_barrier)) cpt1Up++;
                coord_ori[no_noeud_restant[i]][x]=originex; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
            } else if( ((originex >= currTVparsPtr->x1_barrier) && (coord_ori[no_noeud_restant[i]][x] < currTVparsPtr->x1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateDown) ) {
                //            if((originex >= currTVparsPtr->x1_barrier) && (coord_ori[no_noeud_restant[i]][x] < currTVparsPtr->x1_barrier)) cpt1Down++;
                coord_ori[no_noeud_restant[i]][x]=originex; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
            }
        }
        
        if(currTVparsPtr->backwardBarrierY &&  ( (originex >= currTVparsPtr->x1_barrier) && (originex <= currTVparsPtr->x2_barrier) ) ){
            double alea2=alea();
            //        if((originey < currTVparsPtr->y1_barrier) && (coord_ori[no_noeud_restant[i]][y] >= currTVparsPtr->y1_barrier)) cptUp++;
            //        if((originey >= currTVparsPtr->y1_barrier) && (coord_ori[no_noeud_restant[i]][y] < currTVparsPtr->y1_barrier)) cptDown++;
            if( ((originey < currTVparsPtr->y1_barrier) && (coord_ori[no_noeud_restant[i]][y] >= currTVparsPtr->y1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateUp) ) {
                //            if((originey < currTVparsPtr->y1_barrier) && (coord_ori[no_noeud_restant[i]][y] >= currTVparsPtr->y1_barrier)) cpt1Up++;
                coord_ori[no_noeud_restant[i]][y]=originey; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
            } else if( ((originey >= currTVparsPtr->y1_barrier) && (coord_ori[no_noeud_restant[i]][y] < currTVparsPtr->y1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateDown) ) {
                //            if((originey >= currTVparsPtr->y1_barrier) && (coord_ori[no_noeud_restant[i]][y] < currTVparsPtr->y1_barrier)) cpt1Down++;
                coord_ori[no_noeud_restant[i]][y]=originey; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
            }
        }
        
        //printf("coord X apres barriere =%d",coord_ori[no_noeud_restant[i]][x]);
        //printf("coord Y apres barriere=%d",coord_ori[no_noeud_restant[i]][y]);
        //getchar();
     //Effets de bord reflectifs et circulaires pour y
     //if((coord_origine[y]<=0)||(coord_origine[y]>dim_reseau2)) printf("\n coord Y avant=%d",coord_origine[y]);
     coord_ori[no_noeud_restant[i]][y]=int(SortieSupfnPtr1(coord_origine[y],currTVparsPtr->dimRes2));
     coord_ori[no_noeud_restant[i]][y]=int(SortieInffnPtr1(coord_ori[no_noeud_restant[i]][y],currTVparsPtr->dimRes2));
     if(EdgeEffect=="refl")//pour reflection multiples
       {while((coord_ori[no_noeud_restant[i]][y]<=0)||(coord_ori[no_noeud_restant[i]][y]>currTVparsPtr->dimRes2))
		 {coord_ori[no_noeud_restant[i]][y]=int(SortieInffnPtr1(coord_ori[no_noeud_restant[i]][y],currTVparsPtr->dimRes2));
          coord_ori[no_noeud_restant[i]][y]=int(SortieSupfnPtr1(coord_ori[no_noeud_restant[i]][y],currTVparsPtr->dimRes2));
         }
        //if((coord_origine[y]<=0)||(coord_origine[y]>dim_reseau2)) {printf("coord Y apres effets bords=%d",coord_ori[no_noeud_restant[i]][y]);getchar();}
        if(currTVparsPtr->backwardBarrierX && ( (originey >= currTVparsPtr->y1_barrier) && (originey <= currTVparsPtr->y2_barrier) ) ) {
            double alea2=alea();
            //        if((originex < currTVparsPtr->x1_barrier) && (coord_ori[no_noeud_restant[i]][x] >= currTVparsPtr->x1_barrier)) cptUp++;
            //        if((originex >= currTVparsPtr->x1_barrier) && (coord_ori[no_noeud_restant[i]][x] < currTVparsPtr->x1_barrier)) cptDown++;
            if( ((originex < currTVparsPtr->x1_barrier) && (coord_ori[no_noeud_restant[i]][x] >= currTVparsPtr->x1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateUp) ) {
                //            if((originex < currTVparsPtr->x1_barrier) && (coord_ori[no_noeud_restant[i]][x] >= currTVparsPtr->x1_barrier)) cpt1Up++;
                coord_ori[no_noeud_restant[i]][x]=originex; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
            } else if( ((originex >= currTVparsPtr->x1_barrier) && (coord_ori[no_noeud_restant[i]][x] < currTVparsPtr->x1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateDown) ) {
                //            if((originex >= currTVparsPtr->x1_barrier) && (coord_ori[no_noeud_restant[i]][x] < currTVparsPtr->x1_barrier)) cpt1Down++;
                coord_ori[no_noeud_restant[i]][x]=originex; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
            }
        }
        
        if(currTVparsPtr->backwardBarrierY &&  ( (originex >= currTVparsPtr->x1_barrier) && (originex <= currTVparsPtr->x2_barrier) ) ){
            double alea2=alea();
            //        if((originey < currTVparsPtr->y1_barrier) && (coord_ori[no_noeud_restant[i]][y] >= currTVparsPtr->y1_barrier)) cptUp++;
            //        if((originey >= currTVparsPtr->y1_barrier) && (coord_ori[no_noeud_restant[i]][y] < currTVparsPtr->y1_barrier)) cptDown++;
            if( ((originey < currTVparsPtr->y1_barrier) && (coord_ori[no_noeud_restant[i]][y] >= currTVparsPtr->y1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateUp) ) {
                //            if((originey < currTVparsPtr->y1_barrier) && (coord_ori[no_noeud_restant[i]][y] >= currTVparsPtr->y1_barrier)) cpt1Up++;
                coord_ori[no_noeud_restant[i]][y]=originey; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
            } else if( ((originey >= currTVparsPtr->y1_barrier) && (coord_ori[no_noeud_restant[i]][y] < currTVparsPtr->y1_barrier)) && (alea2 > currTVparsPtr->barrierCrossingRateDown) ) {
                //            if((originey >= currTVparsPtr->y1_barrier) && (coord_ori[no_noeud_restant[i]][y] < currTVparsPtr->y1_barrier)) cpt1Down++;
                coord_ori[no_noeud_restant[i]][y]=originey; // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL ameillorer en reflectif.
            }
        }

        //printf("coord X apres barriere et reflexion =%d",coord_ori[no_noeud_restant[i]][x]);
        //printf("coord Y apres barriere=%d",coord_ori[no_noeud_restant[i]][y]);
        //getchar();
       }
    }
  if(effective_disp)//calcul de la dispersion "efficace"
     {effective[grande_dim+ (coord_ori[no_noeud_restant[i]][x]-originex)][x]++;
      effective[grande_dim+ (coord_ori[no_noeud_restant[i]][y]-originey)][y]++;
	  if((coord_ori[no_noeud_restant[i]][x]-originex)!=0 || (coord_ori[no_noeud_restant[i]][y]-originey)!=0) effectiveImmigration[originex][originey]++;
	  cumulEffectiveImmigration[originex][originey]++;
      /*if((coord_ori[no_noeud_restant[i]][x]-originex>dx_max)||(coord_ori[no_noeud_restant[i]][y]-originey>dy_max)) {
      	printf("\n ****CurrentGn=%d,no_noeud_restant=%d",currentGeneration,no_noeud_restant[i]);
      	printf("\noriginex=%d new coordx=%d",originex,coord_ori[no_noeud_restant[i]][x]);
      	printf("\noriginey=%d new coordy=%d",originey,coord_ori[no_noeud_restant[i]][y]);
      	getchar();
      	}*/
      //printf("\n \n effective[%d]X=%lld ",grande_dim+ (coord_ori[no_noeud_restant[i]][x]-originex),effective[grande_dim+ (coord_ori[no_noeud_restant[i]][x]-originex)][x]);
      //printf("\n effective[%d]Y=%lld",grande_dim+ (coord_ori[no_noeud_restant[i]][y]-originey),effective[grande_dim+ (coord_ori[no_noeud_restant[i]][y]-originey)][y]);
	  //if((coord_ori[no_noeud_restant[i]][x]-originex)!=0 || (coord_ori[no_noeud_restant[i]][y]-originey)!=0) cout << "EffectiveImmigration[originex=" << originex << "][originey" << originey << "]++=" << effectiveImmigration[originex][originey] << endl;
      //cout << "cumulEffectiveImmigration[originex=" << originex << "][originey" << originey << "]++=" << cumulEffectiveImmigration[originex][originey] << endl;
	  //getchar();

     }
//cout << "apres disp : " <<  coord_ori[no_noeud_restant[i]][x] << " , " << coord_ori[no_noeud_restant[i]][y] << endl;
//getchar();

 }//fin boucle sur nbre noeud restant

//cout << "crossing rate Up =" << 1. - (1.*cpt1Up)/(1.*cptUp) /*<< ", cpt1Up=" << cpt1Up << ", cptUp=" << cptUp*/;
//cout << " & Down =" << 1. - (1.*cpt1Down)/(1.*cptDown) /* << ", cpt1Down=" << cpt1Down << ", cptDown=" << cptDown*/ << endl;

for(i=1;i<=nbre_noeud_restant_loc;i++)
 {coord_noeud[no_noeud_restant[i]][x]=coord_ori[no_noeud_restant[i]][x];
  coord_noeud[no_noeud_restant[i]][y]=coord_ori[no_noeud_restant[i]][y];
 }

//free_ivector(coord_origine);
} //fin de disp_fixe()

/*****************/
//calcul les probailitÈs cumulÈes de migration pour reseau homogene (= fixe!=0)
void cumul_fixex(int k) {

 int ddx;
 double numer2,denom2,frac2;

//pour x
denom2=0.0;
numer2=0.0;
for (ddx=dx_min;ddx<=dx_max;ddx++)
{denom2+=migra[abs(ddx)];
 if(ddx==k) numer2=denom2;
//cout<<numer2<<" "<<denom2<<endl;
}
frac2= ((double) numer2/denom2);
//cout<<k<<" "<<dx_max<<" "<<frac2<<endl;
cum_fixex[k+dx_max]=frac2;
//cout<<" alo";

}/*fin cumul_fixe*/

/*****************/
//calcul les probailitÈs cumulÈes de migration pour reseau homogene (= fixe!=0)
void cumul_fixey(int k)
{int ddy;
 double numer2,denom2,frac2;

//pour y
denom2=0.0;
numer2=0.0;
for (ddy=dy_min;ddy<=dy_max;ddy++)
{denom2+=migra[abs(ddy)];
 if(ddy==k) numer2=denom2;
}
frac2= ((double) numer2/denom2);
cum_fixey[k+dy_max]=frac2;

}/*fin cumul_fixe*/



/*****************/
//procedure calcul la distance de migration en x quand reseau homogene (= fixe!=0)
int dispx_fixe()
{float alea1;
int ddx=0;
alea1=alea();
if((alea1<cum_fixex[dx_max])&&(alea1>cum_fixex[dx_max-1])) {ddx=0;goto finx;}
if(alea1>cum_fixex[dx_max])
{for(int k=1;k<=dx_max;k++)
  if(alea1<cum_fixex[k+dx_max])
	 {ddx=k;goto finx;
	 }
}
if( alea1 < cum_fixex[dx_max-1])
{if(alea1<cum_fixex[0]) {ddx=dx_min;goto finx;}
for(int k=-1;k>dx_min;k--)
  if(alea1>cum_fixex[k-1+dx_max]) {ddx=k;goto finx;}
}
finx:
	//cout << "dx=" << ddx << endl;
return ddx;
} /*fin dispx_fixe*/

/*****************/
//procedure calcul la distance de migration en y quand reseau homogene (= fixe!=0)
int dispy_fixe()
{float alea1;
int ddy=0;

alea1=alea();
if(currTVparsPtr->dimRes2/currTVparsPtr->vide==1 || ((alea1<cum_fixey[dy_max])&&(alea1>cum_fixey[dy_max-1]))) {ddy=0;goto finy;}
if(alea1>cum_fixey[dy_max])
{for(int k=1;k<=dy_max;k++)
  if(alea1<cum_fixey[k+dy_max])
	 {ddy=k;goto finy;
	 }
}
if( alea1 < cum_fixey[dy_max-1])
{if(alea1<cum_fixey[0]) {ddy=dy_min;goto finy;}
for(int k=-1;k>dy_min;k--)
  if(alea1>cum_fixey[k-1+dy_max]) {ddy=k;goto finy;}
}
finy:
	//cout << "ddy=" << ddy << endl;
return ddy;
} /*fin dispy_fixe*/

/*---------------------procedure nombre aleatoire des noeuds------------------*/

void nbre_aleat_noeud()
{ int n;
  float aleat;
using namespace NS_coal_tree;

	if (currTVparsPtr==NULL) {
	    cerr<<"(!!!) From nbre_aleat_noeud: currTVparsPtr not valued ";
	    if (cinGetOnError) cin.get();
	    exit(-1);
    }

for(int j=1;j<=nbre_noeud_restant;j++)
 {if(currTVparsPtr->DensSpatiallyHeterog) n=densite[coord_noeud[no_noeud_restant[j]][x]]
		  [coord_noeud[no_noeud_restant[j]][y]];
    else n=currTVparsPtr->dens;

//printf("cordx:%d, coordy:%d; "
//          ,coord_noeud[no_noeud_restant[j]][x]
//			 ,coord_noeud[no_noeud_restant[j]][y]);
// printf("densi=n:%d\n",n);

     aleat=1.*ploidy*n*alea();
 //pb : remplacer : doit etre bcp plus rapide quand grands demes
 //aleat_noeud[no_noeud_restant[j]]=( (int) aleat) +1;
 aleat_noeud[no_noeud_restant[j]]=( (int) aleat) + 1;
 //for(jj=1;jj<=2*n;jj++)
 // {if (aleat<=jj) {aleat_noeud[no_noeud_restant[j]]=jj; break;/* goto fin;*/}
 // }
 /*fin:*/
//printf("aleat: %f %d nbre alea noeud[noeud_restant[%d]]: %d\n"
//			  ,aleat,( (int) aleat) + 1,j,aleat_noeud[no_noeud_restant[j]]);
//cin.get();
 }
}
/*--------------------------------------------------------------------------*/
/*------------------procedure gamma_mut pour tx de mut variable-------------*/
/*------------------et loi binormale allant avec----------------------------*/

double gamma_mut(double rate)/*de parameres alpha=2,lambda=0.00025)*/
{
 float aleat1,aleat2;
 double x1,x2,mut,mutl;
do
 {mut=0;
 for(long int i=1;i<=2/*alpha*/;i++)
  {aleat1=alea();
  aleat2=alea();
  x1=sqrt(-2.0*log(aleat1))*cos(2*PI*aleat2);/*loi normale(0.1)*/
  x2=sqrt(-2.0*log(aleat2))*cos(2*PI*aleat1);/*idem*/
  mut+=x1*x1+x2*x2;
  }
 mutl=mut*rate/*lambda*//2;
 }while(mutl<=0);
//printf("mutl= %6f",mutl);
//getchar();
return(mutl);
}

/*-------------counts the number of alleles in the total simulated genetic sample after having added mutations----------------*/
/*-------------also stores the allele ID's----------------*/

void allele_count()
{
	int j;
	allelesAtCurrentLocus = 1;

	//	IF "simulated data = sequence type" then count alleles separately
	if ((model_mut.compare("ISM") == 0) || (model_mut.compare("SNP") == 0)
			|| (model_mut.compare("JC69") == 0) || (model_mut.compare("K80") == 0)
			|| (model_mut.compare("F81") == 0)
			|| (model_mut.compare("HKY85") == 0) || (model_mut.compare("TN93") == 0)) {
		//	loop for counting the total number of haplotypes in the sample
		noeud[1][locus].etat_allele = allelesAtCurrentLocus;
		alleleID.push_back(allelesAtCurrentLocus);
		for (int i = 2; i <= n_genes_total; i++) {
			j = 1;
			while ((j != i-1) && (noeud[j][locus].sequence != noeud[i][locus].sequence))
				j++;
			if (noeud[j][locus].sequence != noeud[i][locus].sequence) {
				allelesAtCurrentLocus++;
				noeud[i][locus].etat_allele = allelesAtCurrentLocus;
				alleleID.push_back(allelesAtCurrentLocus);
			}
			else
				noeud[i][locus].etat_allele = noeud[j][locus].etat_allele;
		}
	}
	else {
		alleleID.push_back(noeud[1][locus].etat_allele);
		for (int i = 2; i <= n_genes_total; i++) {
			j = 1;
			while ((j != i-1) && (noeud[j][locus].etat_allele != noeud[i][locus].etat_allele))
				j++;
			if (noeud[j][locus].etat_allele != noeud[i][locus].etat_allele) {
				allelesAtCurrentLocus++;
				alleleID.push_back(noeud[i][locus].etat_allele);
			}
		}
	}

	allele_max[locus] = 0;
	allele_min[locus] = 2147483647;
	for (int i = 1; i <= n_genes_total; i++){
		if (noeud[i][locus].etat_allele > allele_max[locus])
			allele_max[locus] = noeud[i][locus].etat_allele;
		if (noeud[i][locus].etat_allele < allele_min[locus])
			allele_min[locus] = noeud[i][locus].etat_allele;
	}
	allele_nmbr[locus] = allelesAtCurrentLocus;
	mean_allele_nmbr += ((double) allelesAtCurrentLocus) / (double) n_locus;
}



/*---------------------------------------------------------------------------*/

/*-------------procedure Kini pour donner l'etat allelic du MRCA-------------*/
int K_ini()
{double aleat;
aleat=alea();
//printf("\n aleat=%f,fraction=%f",aleat,(((double)i)/((double)(Kmax-Kmin+1))));
if(Kini==0){
for(int i=Kmin;i<=Kmax;i++) if(aleat <= (((double)i)/((double)(Kmax-Kmin+1)))) {return i;}
 }
else return Kini;
return(-1);
}
/*---------------------------------------------------------------------------*/

/*------------------generates KAM mutations for the simulated coalescent tree------------------*/

void mutate_under_KAM()
{int aleat;
long int temps;
char lettre;

noeud[nbre_noeud_existant][locus].etat_allele=K_ini();/*etat du noeud ancestral*/
for(int i=nbre_noeud_existant;i>0;i--)
	for(unsigned int k=0;k<noeud[i][locus].descendant.size();k++)
	{if(noeud[noeud[i][locus].descendant[k]][locus].etat_allele==111)
		{noeud[noeud[i][locus].descendant[k]][locus].etat_allele
						 =noeud[i][locus].etat_allele;
		 temps=noeud[i][locus].generation
		 -noeud[noeud[i][locus].descendant[k]][locus].generation;/*nbre generation
																	entre ancetre et descendant*/
		 n_mut=mutation(temps);
		 /*printf("%d", n_mut);
		 getchar();*/
		 if (n_mut!=0)
			{for(int j=1;j<=n_mut;j++)
				{do {aleat=(int) (Kmax*alea() +1);}
						while(aleat
							==noeud[noeud[i][locus].descendant[k]][locus].etat_allele);
				 noeud[noeud[i][locus].descendant[k]][locus].etat_allele=aleat;
				}
			}
		 }
	}
allele_count();
if(migrate_lettre_oui)
 for(int i=1;i<=nbre_noeud_existant;i++)
  { if(noeud[i][locus].etat_allele==1) lettre='a';
	else if(noeud[i][locus].etat_allele==2) lettre='b';
	else if(noeud[i][locus].etat_allele==3) lettre='c';
	else if(noeud[i][locus].etat_allele==4) lettre='d';
	else if(noeud[i][locus].etat_allele==5) lettre='e';
	else if(noeud[i][locus].etat_allele==6) lettre='f';
	else if(noeud[i][locus].etat_allele==7) lettre='g';
	else if(noeud[i][locus].etat_allele==8) lettre='h';
	else if(noeud[i][locus].etat_allele==9) lettre='i';
	else if(noeud[i][locus].etat_allele==10) lettre='j';
	else if(noeud[i][locus].etat_allele==11) lettre='k';
	else if(noeud[i][locus].etat_allele==12) lettre='l';
	else if(noeud[i][locus].etat_allele==13) lettre='m';
	else if(noeud[i][locus].etat_allele==14) lettre='n';
	else if(noeud[i][locus].etat_allele==15) lettre='o';
	else if(noeud[i][locus].etat_allele==16) lettre='p';
	else if(noeud[i][locus].etat_allele==17) lettre='q';
	else if(noeud[i][locus].etat_allele==18) lettre='r';
	else if(noeud[i][locus].etat_allele==19) lettre='s';
	else if(noeud[i][locus].etat_allele==20) lettre='t';
	else if(noeud[i][locus].etat_allele==21) lettre='u';
	else if(noeud[i][locus].etat_allele==22) lettre='v';
	else if(noeud[i][locus].etat_allele==23) lettre='w';
	else if(noeud[i][locus].etat_allele==24) lettre='x';
	else if(noeud[i][locus].etat_allele==25) lettre='y';
	else if(noeud[i][locus].etat_allele==26) lettre='z';
	else if(noeud[i][locus].etat_allele==27) lettre='0';
	else if(noeud[i][locus].etat_allele==28) lettre='1';
	else if(noeud[i][locus].etat_allele==29) lettre='2';
	else if(noeud[i][locus].etat_allele==30) lettre='3';
	else if(noeud[i][locus].etat_allele==31) lettre='4';
	else if(noeud[i][locus].etat_allele==32) lettre='5';
	else if(noeud[i][locus].etat_allele==33) lettre='6';
	else if(noeud[i][locus].etat_allele==34) lettre='7';
	else if(noeud[i][locus].etat_allele==35) lettre='8';
	else if(noeud[i][locus].etat_allele==36) lettre='9';
	else {lettre=-1;printf("probleme : too many alleles (>36) to translate to letter code");getchar();}
				 //printf("\n lettre=%c",lettre);
    noeud[i][locus].etat_allele_lettre=lettre;

  }


}
/*---------------------------------------------------------------------------*/

/*------------------generates SMM mutations for the simulated coalescent tree------------------*/


void mutate_under_SMM()
{
int etat_precedent;
long int temps;
double aleat;

for(int i=1;i<=nbre_noeud_existant;i++) noeud[i][locus].etat_allele=111;
noeud[nbre_noeud_existant][locus].etat_allele=K_ini();
				/*etat du noeud ancestral*/
for(int i=nbre_noeud_existant;i>0;i--)
	  { for(unsigned int k=0;k<noeud[i][locus].descendant.size();k++)
		  {if(noeud[noeud[i][locus].descendant[k]][locus].etat_allele==111)
		{noeud[noeud[i][locus].descendant[k]][locus].etat_allele
			=noeud[i][locus].etat_allele;
					  temps=noeud[i][locus].generation
							  -noeud[noeud[i][locus].descendant[k]][locus].generation;
					  n_mut=mutation(temps);
		 if (n_mut!=0)
			{for(int j=1;j<=n_mut;j++)
				{etat_precedent=noeud[noeud[i][locus].descendant[k]][locus].etat_allele;
				precedent1:
				 aleat=alea();
				 if (aleat<0.5) noeud[noeud[i][locus].descendant[k]][locus].etat_allele=(etat_precedent-MotifSize);
				  else noeud[noeud[i][locus].descendant[k]][locus].etat_allele=(etat_precedent+MotifSize);

				 if (noeud[noeud[i][locus].descendant[k]][locus].etat_allele<Kmin) {
                     noeud[noeud[i][locus].descendant[k]][locus].etat_allele=etat_precedent;
                     goto precedent1;/*refait les mutations depassant les bornes, apriori comme dans migraine*/
//                      aleat=alea();// sinon réparti la mutation entre rien et mutation dans l'autre sens, comme dans mathematica
//                      if(aleat<0.5) noeud[noeud[i][locus].descendant[k]][locus].etat_allele=etat_precedent;
//                       else noeud[noeud[i][locus].descendant[k]][locus].etat_allele=etat_precedent+1;
					} else if (noeud[noeud[i][locus].descendant[k]][locus].etat_allele>Kmax) {
                      noeud[noeud[i][locus].descendant[k]][locus].etat_allele=etat_precedent;
                      goto precedent1;/*refait les mutations depassant les bornes, apriori comme dans migraine*/
//                         aleat=alea();// sinon réparti la mutation entre rien et mutation dans l'autre sens, comme dans mathematica
//                         if(aleat<0.5) noeud[noeud[i][locus].descendant[k]][locus].etat_allele=etat_precedent;
//                         else noeud[noeud[i][locus].descendant[k]][locus].etat_allele=etat_precedent-1;
					 }

				}
			}
		}
	}
	  }
allele_count();
}
/*----------------------------------------------------------------------------*/

/*------------------generates IAM mutations for the simulated coalescent tree------------------*/

void mutate_under_IAM()
{
int nbre_etat;
char lettre;
long int temps;

for(int i=1;i<=nbre_noeud_existant;i++) noeud[i][locus].etat_allele=111;
nbre_etat=1;
noeud[nbre_noeud_existant][locus].etat_allele
	=nbre_etat;/*etat du noeud ancestral*/
for(int i=nbre_noeud_existant;i>0;i--)
 {for(unsigned int k=0;k<noeud[i][locus].descendant.size();k++)
	{if(noeud[noeud[i][locus].descendant[k]][locus].etat_allele==111)
		{temps=(noeud[i][locus].generation
		  -noeud[noeud[i][locus].descendant[k]][locus].generation);
		 n_mut=mutation(temps);
		 if (n_mut!=0)
			{for(int j=1;j<=n_mut;j++)
				{nbre_etat++;
				 noeud[noeud[i][locus].descendant[k]][locus].etat_allele=nbre_etat;
				}
			 }
			else noeud[noeud[i][locus].descendant[k]][locus].etat_allele
				=noeud[i][locus].etat_allele;
		 }
	}
 }
allele_count();
if(migrate_lettre_oui)
 for(int i=1;i<=nbre_noeud_existant;i++)
  { if(noeud[i][locus].etat_allele==1) lettre='a';
	else if(noeud[i][locus].etat_allele==2) lettre='b';
	else if(noeud[i][locus].etat_allele==3) lettre='c';
	else if(noeud[i][locus].etat_allele==4) lettre='d';
	else if(noeud[i][locus].etat_allele==5) lettre='e';
	else if(noeud[i][locus].etat_allele==6) lettre='f';
	else if(noeud[i][locus].etat_allele==7) lettre='g';
	else if(noeud[i][locus].etat_allele==8) lettre='h';
	else if(noeud[i][locus].etat_allele==9) lettre='i';
	else if(noeud[i][locus].etat_allele==10) lettre='j';
	else if(noeud[i][locus].etat_allele==11) lettre='k';
	else if(noeud[i][locus].etat_allele==12) lettre='l';
	else if(noeud[i][locus].etat_allele==13) lettre='m';
	else if(noeud[i][locus].etat_allele==14) lettre='n';
	else if(noeud[i][locus].etat_allele==15) lettre='o';
	else if(noeud[i][locus].etat_allele==16) lettre='p';
	else if(noeud[i][locus].etat_allele==17) lettre='q';
	else if(noeud[i][locus].etat_allele==18) lettre='r';
	else if(noeud[i][locus].etat_allele==19) lettre='s';
	else if(noeud[i][locus].etat_allele==20) lettre='t';
	else if(noeud[i][locus].etat_allele==21) lettre='u';
	else if(noeud[i][locus].etat_allele==22) lettre='v';
	else if(noeud[i][locus].etat_allele==23) lettre='w';
	else if(noeud[i][locus].etat_allele==24) lettre='x';
	else if(noeud[i][locus].etat_allele==25) lettre='y';
	else if(noeud[i][locus].etat_allele==26) lettre='z';
	else if(noeud[i][locus].etat_allele==27) lettre='0';
	else if(noeud[i][locus].etat_allele==28) lettre='1';
	else if(noeud[i][locus].etat_allele==29) lettre='2';
	else if(noeud[i][locus].etat_allele==30) lettre='3';
	else if(noeud[i][locus].etat_allele==31) lettre='4';
	else if(noeud[i][locus].etat_allele==32) lettre='5';
	else if(noeud[i][locus].etat_allele==33) lettre='6';
	else if(noeud[i][locus].etat_allele==34) lettre='7';
	else if(noeud[i][locus].etat_allele==35) lettre='8';
	else if(noeud[i][locus].etat_allele==36) lettre='9';
	else {lettre=-1;
          printf("probleme : too many alleles (>36) to translate to letter code");
		  printf("\n allele number %d can not be translate",noeud[i][locus].etat_allele);
		  getchar();
		 }

	noeud[i][locus].etat_allele_lettre=lettre;

  }

}

/*----------------------------------------------------------------------------*/

/*------------------generates TPM mutations for the simulated coalescent tree------------------*/

 void mutate_under_TPM()
{
int delta_mut,etat_precedent;
double aleat;
long int temps;

for(int i=1;i<=nbre_noeud_existant;i++) noeud[i][locus].etat_allele=111;
noeud[nbre_noeud_existant][locus].etat_allele
	=K_ini();/*etat du noeud ancestral*/
for(int i=nbre_noeud_existant;i>0;i--)
 {for(unsigned int k=0;k<noeud[i][locus].descendant.size();k++)
	{if(noeud[noeud[i][locus].descendant[k]][locus].etat_allele==111)
		{noeud[noeud[i][locus].descendant[k]][locus].etat_allele
			=noeud[i][locus].etat_allele;
		 temps=(noeud[i][locus].generation
			-noeud[noeud[i][locus].descendant[k]][locus].generation);
		 n_mut=mutation(temps);
		 if (n_mut!=0)
			{for(int j=1;j<=n_mut;j++)
				{etat_precedent
					  =noeud[noeud[i][locus].descendant[k]][locus].etat_allele;
				 precedent1:
				 aleat=alea();
				 if (aleat<pSMM) delta_mut=1;/*strict SMM*/
				  else delta_mut=loi_geom(var_geom_TPM);/*loi geometrique*/
				  /*gotoxy(1,17);
				  printf("delta_mut= %d\n",delta_mut);
				  getchar();*/
				 aleat=alea();
				 if(aleat<0.5)
					 noeud[noeud[i][locus].descendant[k]][locus].etat_allele
						 =(etat_precedent-delta_mut*MotifSize);
				  else noeud[noeud[i][locus].descendant[k]][locus].etat_allele
					=(etat_precedent+delta_mut*MotifSize);
				 if (noeud[noeud[i][locus].descendant[k]][locus].etat_allele<Kmin)
					{noeud[noeud[i][locus].descendant[k]][locus].etat_allele
					 =etat_precedent;
					 goto precedent1; // comme dans migraine (à vérifier), voir SMM pour comparaison avec mathematica
					}
				 if (noeud[noeud[i][locus].descendant[k]][locus].etat_allele>Kmax)
					 {noeud[noeud[i][locus].descendant[k]][locus].etat_allele
						=etat_precedent;
                         goto precedent1; // comme dans migraine (à vérifier), voir SMM pour comparaison avec mathematica
					 }
				 /*elimine mutations depassant les bornes*/
				}
			}
		}
	}
 }
allele_count();
}

/*----------------------------------------------------------------------------*/

/*------------------generates GSM mutations for the simulated coalescent tree------------------*/

 void mutate_under_GSM()
{
int delta_mut,etat_precedent;
double aleat;
long int temps;

for(int i=1;i<=nbre_noeud_existant;i++) noeud[i][locus].etat_allele=111;
noeud[nbre_noeud_existant][locus].etat_allele
	=K_ini();/*etat du noeud ancestral*/
for(int i=nbre_noeud_existant;i>0;i--)
 {for(unsigned int k=0;k<noeud[i][locus].descendant.size();k++)
	{if(noeud[noeud[i][locus].descendant[k]][locus].etat_allele==111)
		{noeud[noeud[i][locus].descendant[k]][locus].etat_allele
			=noeud[i][locus].etat_allele;
		 temps=(noeud[i][locus].generation
			-noeud[noeud[i][locus].descendant[k]][locus].generation);
		 n_mut=mutation(temps);
		 if (n_mut!=0)
			{for(int j=1;j<=n_mut;j++)
				{etat_precedent
					  =noeud[noeud[i][locus].descendant[k]][locus].etat_allele;
				 //printf("\n etat initial:%d",etat_precedent);
				 precedent1:
				 delta_mut=loi_geom(var_geom_GSM);/*loi geometrique*/
				  /*gotoxy(1,17);
				  printf("delta_mut= %d\n",delta_mut);
				  getchar();*/
				 aleat=alea();
				 if(aleat<0.5)
					 noeud[noeud[i][locus].descendant[k]][locus].etat_allele
						 =(etat_precedent-delta_mut*MotifSize);
				  else noeud[noeud[i][locus].descendant[k]][locus].etat_allele
					=(etat_precedent+delta_mut*MotifSize);
				 //printf("etat apres mut:%d",noeud[noeud[i][locus].descendant[k]][locus].etat_allele);
				 if (noeud[noeud[i][locus].descendant[k]][locus].etat_allele<Kmin)
					{noeud[noeud[i][locus].descendant[k]][locus].etat_allele
					 =etat_precedent;
					 goto precedent1;// comme dans migraine (à vérifier), voir SMM pour comparaison avec mathematica
					}
				 if (noeud[noeud[i][locus].descendant[k]][locus].etat_allele>Kmax)
					 {noeud[noeud[i][locus].descendant[k]][locus].etat_allele
						=etat_precedent;
					  goto precedent1;// comme dans migraine (à vérifier), voir SMM pour comparaison avec mathematica
					 }
				 //printf("etat apres elimination:%d",noeud[noeud[i][locus].descendant[k]][locus].etat_allele);
				 //getchar();
				 /*elimine mutations depassant les bornes*/
				}
			}
		}
	}
 }
allele_count();
}
/*--------------------------------------------------------------------------*/

 /*------------------nucleotide substitution models for the simulated coalescent tree------------------*/
void sequence_mutation_model()
{
	long int temps;
	numSegSites[locus] = 0;

	//	constituting the nucleotide substitution rate matrix
	setSubsRateMatrix();

	//	IF the user hasn't specified a sequence for the MRCA
	if (!defMRCASequence.size()) {
		//	Initializing the ancestral node (i.e. MRCA) based on the equilibrium frequencies of the mutational model
		if ((mutModelIDVector[locus+1] == "K80") || (mutModelIDVector[locus+1] == "JC69"))
			for (int i = 0; i < seqSize; i++)
				noeud[nbre_noeud_existant][locus].sequence.push_back(random_base(fixedBaseFreqs));
		else
			for (int i = 0; i < seqSize; i++)
				noeud[nbre_noeud_existant][locus].sequence.push_back(random_base(varBaseFreqs));
	}
	//	ELSE use the user-specified MRCA
	else
		noeud[nbre_noeud_existant][locus].sequence = defMRCASequence;

	//	looping over every branch in the coalescent starting with the MRCA and its immediate descendants
	for(int i = nbre_noeud_existant; i > 0; i--) {
		for(unsigned int k = 0; k < noeud[i][locus].descendant.size(); k++) {
			//	no. of mutations based on the time in generations between the immediate ancestor and the selected descendant
			temps = (noeud[i][locus].generation - noeud[noeud[i][locus].descendant[k]][locus].generation);
			n_mut = mutation(temps);
			//	Initializing the sequence of the selected descendant with that of its immediate ancestor
			noeud[noeud[i][locus].descendant[k]][locus].sequence = noeud[i][locus].sequence;

			//	Substituting a nucleotide of the selected descendant for each mutation on the current branch of the tree
			for(int j=1;j<=n_mut;j++) {
				//	Uniform sampling between [0,seqSize - 1]
				int basePos = floor(alea() * seqSize);
				if (basePos == seqSize)
					basePos = seqSize - 1;

				if (noeud[noeud[i][locus].descendant[k]][locus].sequence[basePos] == 'A')
					noeud[noeud[i][locus].descendant[k]][locus].sequence[basePos] = random_base(subsRateMatrix[A]);
				else if (noeud[noeud[i][locus].descendant[k]][locus].sequence[basePos] == 'G')
					noeud[noeud[i][locus].descendant[k]][locus].sequence[basePos] = random_base(subsRateMatrix[G]);
				else if (noeud[noeud[i][locus].descendant[k]][locus].sequence[basePos] == 'C')
					noeud[noeud[i][locus].descendant[k]][locus].sequence[basePos] = random_base(subsRateMatrix[C]);
				else if (noeud[noeud[i][locus].descendant[k]][locus].sequence[basePos] == 'T')
					noeud[noeud[i][locus].descendant[k]][locus].sequence[basePos] = random_base(subsRateMatrix[T]);
			}
		}
	}

	//	counts the number of segregating sites w.r.t. the MRCA
	for (int j = 0; j < seqSize; j++)
		for (int i = 1; i <= n_genes_total; i++)
			if (noeud[nbre_noeud_existant][locus].sequence[j] != noeud[i][locus].sequence[j]) {
				++numSegSites[locus];
				break;
			}

	//	counts the number of alleles
	allele_count();

	//	writes a NEXUS type file containing only haplotypes and/or another containing all of the individuals in the sample
	write_haplotype_file();
}

/*--------------randomly draws a nucleotide from the discrete base frequency distribution--------------*/
int random_base(double *baseFreq)
{
	double randNum = alea();

	if (randNum < baseFreq[A])
		return 'A';
	else if (randNum < (baseFreq[A] + baseFreq[G]))
		return 'G';
	else if (randNum < (baseFreq[A] + baseFreq[G] + baseFreq[C]))
		return 'C';
	else
		return 'T';
}

/*------------------sets the nucleotide substitution matrix where rows sum to unity------------------*/
void setSubsRateMatrix()
{
	double *baseFreq;
/*
	The base ratios are used as follows: If the ratio R = X/Y has been specified by the user, then
	Y is assigned a multiplicative factor of 1 and X is assigned R
*/
	double tmp_ratio_TITV = ratio_TITV / 0.5, tmp_ratio_TITI = ratio_TITI;

	//	initializing the 4x4 matrix with zeros
	for (int i = A; i <= T; i++)
		for (int j = A; j <= T; j++)
			subsRateMatrix[i][j] = 0.0;

	//	setting the fixed(i.e. 1/4) and variable base frequencies
	if ((mutModelIDVector[locus+1] == "JC69") || (mutModelIDVector[locus+1] == "K80"))
		baseFreq = fixedBaseFreqs;
	else
		baseFreq = varBaseFreqs;

	//	setting up the transition probabilities along with the respective
	//	transition-transversion bias and/or the transition-transition bias
	subsRateMatrix[A][G] = tmp_ratio_TITI * tmp_ratio_TITV * baseFreq[G];
	subsRateMatrix[A][C] = baseFreq[C];
	subsRateMatrix[A][T] = baseFreq[T];

	subsRateMatrix[G][A] = tmp_ratio_TITI * tmp_ratio_TITV * baseFreq[A];
	subsRateMatrix[G][C] = baseFreq[C];
	subsRateMatrix[G][T] = baseFreq[T];

	subsRateMatrix[C][A] = baseFreq[A];
	subsRateMatrix[C][G] = baseFreq[G];
	subsRateMatrix[C][T] = tmp_ratio_TITV * baseFreq[T];

	subsRateMatrix[T][A] = baseFreq[A];
	subsRateMatrix[T][G] = baseFreq[G];
	subsRateMatrix[T][C] = tmp_ratio_TITV * baseFreq[C];

	//	normalizing by row sums
	for (int i = A; i <= T; i++) {
		double rowSum = 0.0;

		for (int j = A; j <= T; j++)
			rowSum += subsRateMatrix[i][j];

		for (int j = A; j <= T; j++)
			subsRateMatrix[i][j] /= rowSum;
	}
}

 /*----------------------------------------------------------------------------*/

/*------------------generates SNP's based on the coalescent ------------------*/
/*CBR: N.B. Here, the ancestral allele is denoted 1 and the derived allele is denoted 2
as a zero in a GENEPOP file may be read as missing data*/
void mutate_as_SNP(int mutNode)
{
	//	check if it is the first run of mutate_as_SNP
	if (mutNode == nbre_noeud_existant) {
		long int totTime = 0, mutTime;
		//	looping over every branch in the coalescent starting with the MRCA and its immediate descendants
		for (int i = nbre_noeud_existant; i > 0; i--) {
			for (unsigned int k = 0; k < noeud[i][locus].descendant.size(); k++) {
				//	summing the time in generations between the immediate ancestor and the selected descendant
				totTime += (noeud[i][locus].generation - noeud[noeud[i][locus].descendant[k]][locus].generation);
			}
			noeud[i][locus].sequence = "1";	//	Initializing every node on the tree by its ancestral state
		}
		//	generating a random mutation between 0 and totTime
		mutTime = ceil(alea.rand(totTime));
		totTime = 0;

		/*//	looping over every branch in the coalescent starting with the MRCA and its immediate descendants
		for (int i = nbre_noeud_existant; i > 0; i--) {
			for (unsigned int k = 0; k < noeud[i][locus].descendant.size(); k++) {
				//	summing the time in generations between the immediate ancestor and the selected descendant
				totTime += (noeud[i][locus].generation - noeud[noeud[i][locus].descendant[k]][locus].generation);

				//	checks if the unique SNP mutation occurred on the current branch
				if (totTime >= mutTime) {
					mutNode = noeud[i][locus].descendant[k];
					//	recursively add SNP polymorphism in the sub tree with the root at mutNode
					mutate_as_SNP(mutNode);
					break;
				}
			}
		}
		*/
		int i = nbre_noeud_existant;
		unsigned int k=0;
		while(totTime < mutTime && i>0) {
			k=0;
			while( totTime < mutTime && k < noeud[i][locus].descendant.size()) {
				//	summing the time in generations between the immediate ancestor and the selected descendant
				totTime += (noeud[i][locus].generation - noeud[noeud[i][locus].descendant[k]][locus].generation);
				k++;
			}
			i--;
		}
		//	the unique SNP mutation occurred on the current branch
		mutNode = noeud[i+1][locus].descendant[k-1];
		//	recursively add SNP polymorphism in the sub tree with the root at mutNode
		mutate_as_SNP(mutNode);

		//	counts the SNP polymorphism
		allele_count();
	}
	//	once a mutation has been detected, add SNP to the specified subtree recursively
	else {
		noeud[mutNode][locus].sequence = "2";
		for (unsigned int k = 0; k < noeud[mutNode][locus].descendant.size(); k++) {
			//	recursively add SNP polymorphism in the sub tree with the root at mutNode
			mutate_as_SNP(noeud[mutNode][locus].descendant[k]);
		}
	}
}

  /*------------------generates ISM mutations (only DNA sequence samples) for the simulated coalescent tree------------------*/
  //	CBR		NB.  Please be familiar with string::resize before trying to understand the code below

void mutate_under_ISM()
{
	long int temps;
	noeud[nbre_noeud_existant][locus].sequence = "";	//	Initializing the ancestral node (i.e. MRCA)

	//	additionally prudent here in order to ensure a zero number of segregated sites before generating any mutations
	numSegSites[locus] = 0;

	//	looping over every branch in the coalescent starting with the MRCA and its immediate descendants
	for (int i = nbre_noeud_existant; i > 0; i--) {
		for (unsigned int k = 0; k < noeud[i][locus].descendant.size(); k++) {

			//	no. of mutations based on the time in generations between the immediate ancestor and the selected descendant
			temps = (noeud[i][locus].generation - noeud[noeud[i][locus].descendant[k]][locus].generation);
			n_mut = mutation(temps);

			//	updates the descendant locus with its immediate ancestors sequence followed by the
			//	number of segregated sites not having occurred in the ancestors sequence
			noeud[noeud[i][locus].descendant[k]][locus].sequence = noeud[i][locus].sequence;
			noeud[noeud[i][locus].descendant[k]][locus].sequence.resize(numSegSites[locus],'0');

			//	updating the total number of segregating sites where n_mut >= 0
			numSegSites[locus] += n_mut;

			//	adds n_mut (>= 0) number of ADDITIONAL segregated sites to the selected descendant locus (cf. string::resize)
			noeud[noeud[i][locus].descendant[k]][locus].sequence.resize(numSegSites[locus],'1');

			//	adds n_mut (>= 0) number of ADDITIONAL UNsegregated sites to the immediate ancestor locus (cf. string::resize)
			noeud[i][locus].sequence.resize(numSegSites[locus],'0');
		}
	}

	//	adds unsegregated sites to all the genes whose sequence.size() < numSegSites in the sample
	for (int i = 1; i <= n_genes_total; i++)
		noeud[i][locus].sequence.resize(numSegSites[locus],'0');

	//	IF no mutations have occurred in the tree then all genes are set for the MRCA state i.e. 0
	if (numSegSites[locus] == 0)
		for (int i = 1; i <= n_genes_total; i++)
			noeud[i][locus].sequence.resize(1,'0');

	//	converts the number code into nucleotide code
	for (int i = 1; i <= n_genes_total; i++)
		noeud[i][locus].sequence = convert_binary_format(noeud[i][locus].sequence);

	//	counts the number of alleles
	allele_count();

	//	writes a NEXUS type file containing only haplotypes and/or another containing all of the individuals in the sample
	write_haplotype_file();
}

/* writes a NEXUS type file containing only haplotypes and/or another containing all of the individuals in the sample */
 //	CBR:		see also http://en.wikipedia.org/wiki/Nexus_file for the basic NEXUS file format details
void write_haplotype_file()
{
	if (nexusFormat != "F") {
		ofstream hapFile;
		{
			string fichier;
			stringstream stst;
			if (n_locus > 1) {
				if (repet != 1)
					stst << fichier_genepop << "Haplotypes" << "_rep" << rep << "_loc" << locus+1 << ".nex";
				else
					stst << fichier_genepop << "_Haplotypes" << "_loc" << locus+1 << ".nex";
			}
			else {
				if (repet != 1)
					stst << fichier_genepop << "Haplotypes" << "_rep" << rep << ".nex";
				else
					stst << fichier_genepop << "_Haplotypes.nex";
			}
			fichier += stst.str();
			hapFile.open(fichier.c_str(), std::ios::out);
		}
		hapFile << "#NEXUS\n" << endl
				<< "begin data;" << endl
				<< "dimensions ntax=" << allele_nmbr[locus] + 1 << " nchar=" << noeud[1][locus].sequence.size() << ";" << endl
				<< "format datatype=dna symbols=\"ACTG\";" << endl
				<< "matrix" << endl;

		if (model_mut.compare("ISM") == 0) {
			string tmpStr;
			if (numSegSites[locus] == 0)
				tmpStr.assign(1, '0');
			else
				tmpStr.assign(numSegSites[locus], '0');
			hapFile << "Anc\t" << convert_binary_format(tmpStr) << endl;

		}
		else
			hapFile << "Anc\t" << noeud[nbre_noeud_existant][locus].sequence << endl;

		for (int i = 1; i <= allele_nmbr[locus]; i++) {
			for (int j = 1; j <= n_genes_total; j++) {
				if (noeud[j][locus].etat_allele == i) {
					{
						string locstring;
						stringstream stst_num;
						stst_num << "00" << i;
						stst_num >> locstring;
						locstring = locstring.substr(locstring.size() - 3, locstring.size());
						hapFile << locstring << '\t';
					}
					hapFile << noeud[j][locus].sequence << endl;
					break;
				}
			}
		}
		hapFile << ';' << endl
				<< "end;" << endl;
		hapFile.close();
	}

	if (nexusFormat == "Haplotypes_and_Individuals") {
		ofstream allIndFile;
		{
			string fichier;
			stringstream stst;
			if (n_locus > 1) {
				if (repet != 1)
					stst << fichier_genepop << "AllIndividuals" << "_rep" << rep << "_loc" << locus+1 << ".nex";
				else
					stst << fichier_genepop << "_AllIndividuals" << "_loc" << locus+1 << ".nex";
			}
			else {
				if (repet != 1)
					stst << fichier_genepop << "AllIndividuals" << "_rep" << rep << ".nex";
				else
					stst << fichier_genepop << "_AllIndividuals.nex";
			}
			fichier += stst.str();
			allIndFile.open(fichier.c_str(), std::ios::out);
		}
		allIndFile << "#NEXUS\n" << endl
				   << "begin data;" << endl
				   << "dimensions ntax=" << n_genes_total + 1 << " nchar=" << noeud[1][locus].sequence.size() << ";" << endl
				   << "format datatype=dna symbols=\"ACTG\";" << endl
				   << "matrix" << endl;

		if (model_mut.compare("ISM") == 0) {
			string tmpStr;
			if (numSegSites[locus] == 0)
				tmpStr.assign(1, '0');
			else
				tmpStr.assign(numSegSites[locus], '0');
			allIndFile << "Anc\t" << convert_binary_format(tmpStr) << endl;

		}
		else
			allIndFile << "Anc\t" << noeud[nbre_noeud_existant][locus].sequence << endl;

		for (int i = 1; i <= n_genes_total; i++) {
			{
				string locstring;
				stringstream stst_num;
				stst_num << "00" << i;
				stst_num >> locstring;
				locstring = locstring.substr(locstring.size() - 3, locstring.size());
				allIndFile << locstring << '\t';
			}
			allIndFile << noeud[i][locus].sequence << endl;
		}
		allIndFile << ';' << endl
				   << "end;" << endl;
		allIndFile.close();
	}
}

/*-------------------coverts a binary sequence into an equivalent ATGC format---------------------------*/
string convert_binary_format(string binSequence)
{
	string atgcSequence;
	string::iterator it;
    {
		stringstream stst;
    	for (it = binSequence.begin(); it < binSequence.end(); it++) {
    		if (*it == '0')
    			stst << 'A';
    		else
    			stst << 'T';
    	}
    	atgcSequence += stst.str();
    }
    return atgcSequence;
}

/*-------------------coverts a quartery sequence into an equivalent ATGC format---------------------------*/
string convert_quartery_format(string quartSequence)
{
	string atgcSequence;
	string::iterator it;
    {
		stringstream stst;
    	for (it = quartSequence.begin(); it < quartSequence.end(); it++) {
    		if (*it == '0')
    			stst << 'A';
    		else if (*it == '1')
    			stst << 'G';
    		else if (*it == '2')
    			stst << 'C';
    		else if (*it == '3')
    			stst << 'T';
    		else {
    			cout << "Unrecognized character in sequence data! Aborting IBDSim."
    				 << "\nPress any key to continue..." << endl;
    			cin.get();
    		}
    	}
    	atgcSequence += stst.str();
    }
    return atgcSequence;
}

/*-------------------returns the number of differences between two strings---------------------------*/
int seqPairDiff(string str1, string str2)
{
	string::iterator iter1, iter2;
	int numPairDiff = 0;

	for (iter1 = str1.begin(), iter2 = str2.begin(); (iter1  < str1.end() || iter2  < str2.end()); iter1++, iter2++)
		if (*iter1 != *iter2)
			++numPairDiff;

	return numPairDiff;
}

/*----------returns the number of transitions(TI) and transversions(TV) between two strings----------*/
double seqPairNucSubDiff(string str1, string str2)
{
	string::iterator iter1, iter2;
	int numTI = 0;
	int numTV = 0;

	for (iter1 = str1.begin(), iter2 = str2.begin(); (iter1  < str1.end() || iter2  < str2.end()); iter1++, iter2++) {
		switch (*iter1) {
		case 'A':
			if ((*iter2 == 'C') || (*iter2 == 'T'))
				++numTV;
			else if (*iter2 == 'G')
				++numTI;
			break;
		case 'G':
			if ((*iter2 == 'C') || (*iter2 == 'T'))
				++numTV;
			else if (*iter2 == 'A')
				++numTI;
			break;
		case 'C':
			if ((*iter2 == 'A') || (*iter2 == 'G'))
				++numTV;
			else if (*iter2 == 'T')
				++numTI;
			break;
		case 'T':
			if ((*iter2 == 'A') || (*iter2 == '1'))
				++numTV;
			else if (*iter2 == 'C')
				++numTI;
			break;
		}
	}

	if (!numTV)
		return 0;
	else
		return (double) numTI / numTV;
}

/*-----------------ecriture entetes fichiers genepop-------------------------*/

void entete_fichier_genepop(int a/*for predisp samples*/)
{
	string fichier = fichier_genepop;
    if (repet!=1) {
		stringstream stst;
        stst << rep;
        fichier += stst.str();
    }
    if(a>0) fichier += "_predisp";
	if (Genepopfile_extension.size()>0) fichier += Genepopfile_extension;
	else if (txt_extensionbool) fichier += ".txt";


    fgenepop.open(fichier.c_str(), std::ios::out);

	if (!fgenepop.is_open()) {
		cerr << "Could not open the file " << fichier << endl;
		cerr << "\nAborting IBDSim...Press any key to exit.\n";
		if(cinGetOnError)
			cin.get();
		exit(-1);
	}
	if(a>0) fgenepop << "This predisp sample file has been generated by the IBDSim program." << endl;
        else fgenepop << "This sample file has been generated by the IBDSim program." << endl;
	for (int i = 0; i < n_locus; i++)
		fgenepop << "locus" << i+1 << "_" << mutModelIDVector[i+1] << endl;


/*fgenepop ferme dans procedure ecriture_fichier*/
}
/*---------------------------------------------------------------------------*/
/*-------------------ecriture fichiers fgenepop-------------------------------*/
void ecriture_fichier_genepop(void)
{
	using namespace NS_coal_tree;
    int sampleNb=1;
    if(predispbool) sampleNb=2;

    for(int a=0;a<sampleNb;a++){

        entete_fichier_genepop(a);/*fgenepop ouvert dans procedure entete_fichier*/

        string locstring;

        int startIndex=1,endIndex=n_genes[0];
        if(a==1) {startIndex=n_genes[0]+1;endIndex=n_genes_total;}

        for (int i = startIndex; i <= endIndex; i = i + ploidy) {
            if ( individualizeAllSamplesbool ) fgenepop << "pop \n";
             else if ((!(groupAllSamplesbool && i!=1) && ((i-1) % (ploidy*dens_sample[a])) == 0) || i == 1  ) fgenepop << "pop \n";
            if(!genepopNoCoordbool) fgenepop << coord_individus[i][x] << " " << coord_individus[i][y] << " , ";
             else if(ploidy==2) fgenepop << (int) (i+1)/ploidy << ", "; else fgenepop << (int) (i+1)/ploidy << ", ";
            for (int j = 0; j < n_locus; j++) {
                for (int pl = 0; pl < ploidy; pl++) {
                    if (noeud[i+pl][j].etat_allele > 999) {
    #ifdef GOTO
                        _gotoxy(0,17);
    #endif
                        cerr << "An allelic state has been found to be larger than 999," << endl;
                        cerr << "IBDSim cannot write Genepop files with more than 3 digits" << endl;
                        cerr << "Choose a lower mutation rate or a mutation model with bounds," << endl;
                        cerr << "\nand start the program again...." << endl;
                        if (cinGetOnError)
                            cin.get();
                        exit(-1);
                    }
                    else {
                        stringstream stst;
                        stst << "00" << noeud[i+pl][j].etat_allele;
                        stst >> locstring;
                        locstring = locstring.substr(locstring.size() - 3, locstring.size());
                        fgenepop << locstring;
                    }
                } //ploidy
                fgenepop << " ";
            } // end loop over loci
            fgenepop << endl;
        }
        fgenepop.close();
    }
}

/*-------------------ecriture fichiers fgeneland1-------------------------------*/
void ecriture_fichier_genotypes(void)
{
    using namespace NS_coal_tree;
    string locstring;
    string nomfich;
    int sampleNb=1;
    if(predispbool) sampleNb=2;

    for(int a=0;a<sampleNb;a++){

        nomfich=fichier_geneland_geno;
        if (repet!=1) {stringstream stst;
            stst<<rep;
            nomfich+=stst.str();
        }
        if(a>0) nomfich+="_predisp";
        if(txt_extensionbool) nomfich+=".txt";

        fgeneland1.open(nomfich.c_str(),std::ios::out);

        if (!fgeneland1.is_open()) {cerr<<"Could not open  file "<<nomfich<<". I exit."; if(cinGetOnError) cin.get();exit(-1);};

        int startIndex=1,endIndex=n_genes[0];
        if(a==1) {startIndex=n_genes[0]+1;endIndex=n_genes_total;}

        for (int i = startIndex; i <= endIndex; i = i + ploidy) {
            for (int j = 0; j < n_locus; j++) {
                for (int pl = 0; pl < ploidy; pl++) {
                    if (noeud[i+pl][j].etat_allele > 999) {
    #ifdef GOTO
                        _gotoxy(0,17);
    #endif
                        cerr << "An allelic state has been found to be larger than 999," << endl;
                        cerr << "IBDSim cannot write Geneland files with more than 3 digits" << endl;
                        cerr << "Choose a lower mutation rate or a mutation model with bounds," << endl;
                        cerr << "\nand start the program again...." << endl;
                        if (cinGetOnError) cin.get();
                        exit(-1);
                    } else {
                        stringstream stst;
                        stst << "00" << noeud[i+pl][j].etat_allele;
                        stst >> locstring;
                        locstring = locstring.substr(locstring.size() - 3, locstring.size());
                        fgeneland1 << locstring << " ";
                    }
                } //end loop over ploidy
            } // end loop over loci
            fgeneland1 << endl;
        } // end loop over ngenes
        fgeneland1.close();
    }
}

/*-------------------ecriture fichiers fgeneland1-------------------------------*/
void ecriture_fichier_coordonnees(void)
{
    using namespace NS_coal_tree;
    string locstring;
    string nomfich;
    int sampleNb=1;
    if(predispbool) sampleNb=2;

    for(int a=0;a<sampleNb;a++){

        nomfich=fichier_geneland_geno;
        if (repet!=1) {stringstream stst;
            nomfich+="coord_";
            stst<<rep;
            nomfich+=stst.str();
        } else nomfich+=+"_coord";
        if(a>0) nomfich+=+"_predisp";
        if(txt_extensionbool) nomfich+=".txt";

        fgeneland2.open(nomfich.c_str(),std::ios::out);

        int startIndex=1,endIndex=n_genes[0];
        if(a==1) {startIndex=n_genes[0]+1;endIndex=n_genes_total;}

        for (int i = startIndex; i <= endIndex; i = i + ploidy)
            fgeneland2 << coord_individus[i][x] << " " << coord_individus[i][y] << endl;

        fgeneland2.close();
    }
}

/*-------------------ecriture fichiers fmigrate-------------------------------*/
void ecriture_fichier_migrate(void)
{
using namespace NS_coal_tree;
    string locstring;
    string nomfich=fichier_migrate;
    int sampleNb=1;
    if(predispbool) sampleNb=2;

    for(int a=0;a<sampleNb;a++){
        if (repet!=1) {stringstream stst;
            stst<<rep;
            nomfich+=stst.str();
        }
        if(a>0) nomfich+="_predisp";
        if(txt_extensionbool) nomfich+=".txt";
        fmigrate.open(nomfich.c_str(),std::ios::out);
        if (!fmigrate.is_open()) {cerr<<"Could not open  file "<<nomfich<<". I exit."; if(cinGetOnError) cin.get();exit(-1);};
        if(!Specific_Sample_Designbool[a]) fmigrate<<dim_sample1[a]/vide_sampleX[a]*dim_sample2[a]/vide_sampleY[a];//nbre de pop echantillonnÈes
           else fmigrate<<Spec_SampleSize[a];//nbre de pop echantillonnÈes
        fmigrate<<" "<<n_locus<<". "<<nomfich.c_str()<<" under "<<model_mut<<"; ";
        if(!Specific_Sample_Designbool[a]) fmigrate<<dim_sample1[a]<<" x "<<dim_sample2[a];
            else fmigrate<<Spec_SampleSize[a];
        if(ploidy==2) fmigrate<<" diploid individuals ";
        else fmigrate<<" haploid idividuals ";
        fmigrate<<" evolving at G=0 on a "<<TVpars[0].dimRes1<<" x "<<TVpars[0].dimRes2<<" lattice";
        if (EdgeEffect=="abs") fmigrate<<" with absorbing boundaries";
        if (EdgeEffect=="refl") fmigrate<<" with reflecting boundaries";

        fmigrate<<"; mutation proba="<<_mu;
        fmigrate<<"; replicate # "<<rep;
        fmigrate<<"; Mrca Moy= "<<mrca_moy<<"; MRCA MAX="<<mrca_max<<"\n";

        int startIndex=1,endIndex=n_genes[0];
        if(a==1) {startIndex=n_genes[0]+1;endIndex=n_genes_total;}

        for (int i = startIndex; i <= endIndex; i = i + ploidy) {
            if (((i-1) % (ploidy*dens_sample[a]))==0 || i==1   )
                    fmigrate<<dens_sample[a]<<" pop"<<coord_individus[i][x]<<coord_individus[i][y]<<endl;
            fmigrate<<(i+1)/2<<"     ,";
            for(int j=0;j<n_locus;j++) {
               for (int pl=0;pl<ploidy;pl++) {
                if(pl==1) fmigrate<<".";
                if(!migrate_lettre_oui) fmigrate<<" "<<noeud[i+pl][j].etat_allele;
                else fmigrate<<" "<<noeud[i+pl][j].etat_allele_lettre;
               }
             }
             fmigrate<<endl;
        }
        fmigrate.close();
    }
}


/*--------------ecriture fichiers DG2002------------------------------------*/
void ecriture_fichiers_DG2002KAM(void)
{int pop,PopNbr,***nb_alleles_etat,*petit2,*grand,nb_etats;
using namespace NS_coal_tree;
//double ***nb_alleles_etat;
string fichier1;
int sampleNb=1,local_sample_size=0;
if(predispbool) sampleNb=2;


    for(int a=0;a<sampleNb;a++){

        nb_etats=Kmax-Kmin+1;

       if(!Specific_Sample_Designbool[a]) local_sample_size=dim_sample1[a]*dim_sample2[a];
           else local_sample_size=Spec_SampleSize[a];
        nb_alleles_etat=(int ***) itab3(nb_etats+1,n_locus,local_sample_size+1);
        petit2=(int *) ivector(n_locus);
        grand=(int *) ivector(n_locus);

        fichier1=fichier_stepsim;
        if (repet!=1) {stringstream stst;
            stst<<rep;
            fichier1+=stst.str();
        }
        if(a>0) fichier1+="_predisp";
        if(txt_extensionbool) fichier1+=".txt";

    /*if(dossier_nbre_allele) {
        //a=malloc(100*sizeof(char *));
    //#ifdef WIN32
    //	sprintf(chaine, ".\\%d alleles", allele_nmbr[1]);
    //	mkdir(chaine);
    //    sprintf(chaine, ".\\%d alleles\\", allele_nmbr[1]);
    //#else //speculative but should allow compilation at least
    #ifdef WIN32
    //not really tested but avoids Dev-C++ error on mkdir()
        sprintf(chaine, "mkdir ./%d alleles", allele_nmbr[1]);
        system(chaine);
    #else
        sprintf(chaine, "./%d alleles", allele_nmbr[1]);
        mkdir(chaine,0777);
    #endif
        sprintf(chaine, "./%d alleles/", allele_nmbr[1]);
    //#endif
        strcat(chaine,fichier1);
        strcpy(fichier1,chaine);

    }*/


    if((fDG2002=fopen(fichier1.c_str(),"w"))==0)
        {printf("ecriture_fichiers_DG2002KAM() cannot open file: %s\n",fichier1.c_str());
         getchar();
         //exit(1);
        }

    for(int j=0;j<n_locus;j++) {
        petit2[j]=nb_etats;
        grand[j]=0;
    }
    pop=0;

    if(!Specific_Sample_Designbool[a]) {
        for(int j=0;j<n_locus;j++)
            for(int k=1;k<=dim_sample1[a];k++)
                for(int k2=1;k2<=dim_sample2[a];k2++){
                    pop++;
                    for(int l=1;l<=nb_etats;l++) nb_alleles_etat[l][j][pop]=0;
                }
    } else {
        for(int j=0;j<n_locus;j++)
            for(int k=1;k<=Spec_SampleSize[a];k++) {
                for(int l=1;l<=nb_etats;l++) nb_alleles_etat[l][j][k]=0;
            }
    }

    int startIndex=1,endIndex=n_genes[0];
    if(a==1) {startIndex=n_genes[0]+1;endIndex=n_genes_total;}

     for(int j=0;j<n_locus;j++)
         for (int i = startIndex; i <= endIndex; i++) {
            if(!Specific_Sample_Designbool[a]) pop=(coord_individus[i][y]-ymin_sample[a])*dim_sample1[a]
                                                    +(coord_individus[i][x]-xmin_sample[a]+vide_sampleX[a])/vide_sampleY[a];
            else pop=(int) (floor(i/ploidy/dens_sample[a])+1);
            //printf("\n locus=%d; gene=%d;pop=%d;etat=%d",j+1,i,pop,noeud[i][j].etat_allele);
            //getchar();
            for(int l=1;l<=nb_etats;l++) {
                if(noeud[i][j].etat_allele==(Kmin+l-1)) {
                    nb_alleles_etat[l][j][pop]+=1;
                    //printf("nb_allele_etat[%d][%d][%d]=%d",l,j,pop,nb_alleles_etat[l][j][pop]);
                    if(l<petit2[j]) petit2[j]=l;
                    if(l>grand[j]) grand[j]=l;
                }
            }
        }
    //for(int j=0;j<n_locus;j++) printf("\n Locus %d griffiths2() OK alleles between %d and %d\n",j+1,petit2[j],grand[j]);
    //getchar();
    if(!Specific_Sample_Designbool[a]) PopNbr=dim_sample1[a]*dim_sample2[a];
    else PopNbr=Spec_SampleSize[a];
    for(int j=0;j<n_locus;j++) {
     for(int i=1;i<=PopNbr;i++) {
      for(int l=1;l<=nb_etats;l++){
        if(nb_alleles_etat[l][j][i]>0 || AllStates) {
                fprintf(fDG2002,"%d",nb_alleles_etat[l][j][i]);
                if(l!=nb_etats) fprintf(fDG2002," ");
        }
      }
      fprintf(fDG2002,";\n");
     }
    }

    free_itab3(nb_alleles_etat,nb_etats+1,n_locus);
    free_ivector(petit2);
    free_ivector(grand);

    if(fclose(fDG2002)) printf("file close error.\n");
}
}
/*----------------------------------------------------------------------------*/

/*--------------ecriture fichiers DG2002, Stepsim et migraine pour SMM------------------------------------*/
void ecriture_fichiers_DG2002SMM(void)
{
using namespace NS_coal_tree;
    //FR->RL voir si on peut virer cette fn ou la fusionner avec la precedente

    int pop,PopNbr,***nb_alleles_etat,*petit2,*grand,nb_etats;
//double ***nb_alleles_etat;
string fichier1;
int sampleNb=1,local_sample_size=0;
if(predispbool) sampleNb=2;

nb_etats=Kmax-Kmin+1;


    for(int a=0;a<sampleNb;a++){

        if(!Specific_Sample_Designbool[a]) local_sample_size=dim_sample1[a]*dim_sample2[a];
        else local_sample_size=Spec_SampleSize[a];
        nb_alleles_etat=(int ***) itab3(nb_etats+1,n_locus,local_sample_size+1);
        petit2=(int *) ivector(n_locus);
        grand=(int *) ivector(n_locus);

        fichier1=fichier_stepsim;
        if (repet!=1) {stringstream stst;
            stst<<rep;
            fichier1+=stst.str();
        }
        if(a>0) fichier1+="_predisp";
        if(txt_extensionbool) fichier1+=".txt";

        /*if(dossier_nbre_allele) {
            //a=malloc(100*sizeof(char *));
        //#ifdef WIN32
        //	sprintf(chaine, ".\\%d alleles", allele_nmbr[1]);
        //	mkdir(chaine);
        //    sprintf(chaine, ".\\%d alleles\\", allele_nmbr[1]);
        //#else //speculative but should allow compilation at least
        #ifdef WIN32
        //not really tested but avoids Dev-C++ error on mkdir()
        //marche pas avec wxdevcpp... donc virer dans les nouvelles versions
            sprintf(chaine, "mkdir ./%d alleles", allele_nmbr[1]);
            system(chaine);
        #else
            sprintf(chaine, "./%d alleles", allele_nmbr[1]);
            mkdir(chaine,0777);
        #endif
            sprintf(chaine, "./%d alleles/", allele_nmbr[1]);
        //#endif
            strcat(chaine,fichier1);
            strcpy(fichier1,chaine);

        }*/


        if((fDG2002=fopen(fichier1.c_str(),"w"))==0)
            {printf("ecriture_fichiers_DG2002SMM() cannot open file: %s\n",fichier1.c_str());
             getchar();
             //exit(1);
            }

        for(int j=0;j<n_locus;j++) {
            petit2[j]=nb_etats;
            grand[j]=0;
        }

        pop=0;

        if(!Specific_Sample_Designbool[a]) {
            for(int j=0;j<n_locus;j++)
                for(int k=1;k<=dim_sample1[a];k++)
                    for(int k2=1;k2<=dim_sample2[a];k2++){
                        pop++;
                        for(int l=1;l<=nb_etats;l++) nb_alleles_etat[l][j][pop]=0;
                    }
        } else {
            for(int j=0;j<n_locus;j++)
                for(int k=1;k<=Spec_SampleSize[a];k++) {
                    for(int l=1;l<=nb_etats;l++) nb_alleles_etat[l][j][k]=0;
                }
        }
        int startIndex=1,endIndex=n_genes[0];
        if(a==1) {startIndex=n_genes[0]+1;endIndex=n_genes_total;}

        for(int j=0;j<n_locus;j++)
            for (int i = startIndex; i <= endIndex;i++){
            if(!Specific_Sample_Designbool[a]) pop=(coord_individus[i][y]-ymin_sample[a])*dim_sample1[a]+(coord_individus[i][x]-xmin_sample[a]+vide_sampleX[a])/vide_sampleY[a];
                else pop=(int) (floor(i/ploidy/dens_sample[a])+1);
                //printf("\n locus=%d; gene=%d;pop=%d;etat=%d",j+1,i,pop,noeud[i][j].etat_allele);
                //getchar();
                for(int l=1;l<=nb_etats;l++)
                  {	if(noeud[i][j].etat_allele==(Kmin+l-1)) {
                        nb_alleles_etat[l][j][pop]+=1;
                        //printf("nb_allele_etat[%d][%d][%d]=%d",l,j,pop,nb_alleles_etat[l][j][pop]);
                        if(l<petit2[j]) petit2[j]=l;
                        if(l>grand[j]) grand[j]=l;
                    }
                  }
                }
        //for(int j=0;j<n_locus;j++) printf("\n Locus %d griffiths2() OK alleles between %d and %d\n",j+1,petit2[j],grand[j]);
        //getchar();
        if(!Specific_Sample_Designbool[a]) PopNbr=dim_sample1[a]*dim_sample2[a];
        else PopNbr=Spec_SampleSize[a];
        for(int j=0;j<n_locus;j++) {
            if(j!=1) fprintf(fDG2002,"\n");
         //nbre d'ind ÈchantillonnÈs par pop
         for(int i=1;i<=PopNbr;i++) fprintf(fDG2002,"%d",ploidy*dens_sample[a]);
         fprintf(fDG2002,"\n");
         for(int i=1;i<=PopNbr;i++) {
                fprintf(fDG2002,"pop\n");
          for(int l=petit2[j];l<=grand[j];l++){
            if(nb_alleles_etat[l][j][i]>0 || AllStates) {
                    fprintf(fDG2002,"%d %d\n",l/*-petit2[j]*/,nb_alleles_etat[l][j][i]);
            }
          }
         }
        }


        free_itab3(nb_alleles_etat,nb_etats+1,n_locus);
        free_ivector(petit2);
        free_ivector(grand);

        if(fclose(fDG2002)) printf("file close error.\n");
    }
}
/*----------------------------------------------------------------------------*/

/*---------ecriture des moyennes des calculs ds fichiers moy_ et suiviQ_--------------*/

void ecriture_fichier_moyennes()
{string moy,suiv;
 //long int  i;

int sampleNb=1;
if(predispbool) sampleNb=2;

for(int a=0;a<sampleNb;a++){


    if(a==1) {
     moy="Various_Statistics_predisp.txt";
     suiv="Iterative_IdProb_predisp.txt";
    } else {
        moy="Various_Statistics_postdisp.txt";
        suiv="Iterative_IdProb_postdisp.txt";
    }


    if(rep==1) {
        if((fmoyenne=fopen(moy.c_str(),"w"))==NULL) {
            cerr << "ecriture fichiers cannot open file: " << moy << endl;
            cerr << "I exit" << endl;
            if(cinGetOnError) cin.get();
            exit(1);
        }

        if(suiviQ && !Specific_Sample_Designbool[a]) {
            if((fsuivi=fopen(suiv.c_str(),"w"))==NULL) {
                cerr << "ecriture fichiers cannot open file: " << suiv << endl;
                cerr << "I exit" << endl;
                if(cinGetOnError) cin.get();
                exit(1);
            }
            if(!noHeader) {
                fprintf(fsuivi,"#   ");
                for(int i=0;i<(/*max(TVpars[0].dimRes1-1,*/max(dim_sample1[a],dim_sample2[a])/*)*//*/vide_sample[a]*/);i++)
                    fprintf(fsuivi,"  Qr(%2d)     meanQr(%2d)   ",i,i);
                fprintf(fsuivi,"\n");
            }
        }
        fprintf(fmoyenne,"run: %s with ",fichier_genepop.c_str());
        fprintf(fmoyenne,"%d simulated data files of %d loci;\n",repet,n_locus);
    } else {/*fin if(rep==1)*/
        if((fmoyenne=fopen(moy.c_str(),"a"))==NULL) {
            cerr << "ecriture fichiers cannot open file: " << moy << endl;
            cerr << "I exit" << endl;
            if(cinGetOnError) cin.get();
            exit(1);
        }
        if(suiviQ && !Specific_Sample_Designbool[a])
            if((fsuivi=fopen(suiv.c_str(),"a"))==NULL) {
                cerr << "ecriture fichiers cannot open file: " << suiv << endl;
                cerr << "I exit" << endl;
                if(cinGetOnError) cin.get();
                exit(1);
            }

    }

    if(suiviQ  && !Specific_Sample_Designbool[a]) {
        fprintf(fsuivi,"%d ",rep);
        for(int i=0;i<(/*max(TVpars[0].dimRes1,*/max(dim_sample1[a],dim_sample2[a])/*)*/);i++) {
            if(rep!=repet)fprintf(fsuivi,"%.10f %.10f ",Q1_moy[a][i],Q1_moy_glob[a][i]*repet/rep);
             else fprintf(fsuivi,"%.10f %.10f ",Q1_moy[a][i],Q1_moy_glob[a][i]);
        }
        fprintf(fsuivi,"\n");
        fflush(fsuivi);
    }
    
    
    if(rep==repet) {
        fprintf(fmoyenne,"MRCA MOY (total sample) = %8.6e\n",mrca_moy_glob);
        fprintf(fmoyenne,"MRCA MAX (total sample) = %lu\n",mrca_max);
        fprintf(fmoyenne,"Monomorphs (total sample) = %8.6e\n",mono);
        fprintf(fmoyenne,"tooManyMut (total sample) = %8.6e\n",tooManyMut);
        if(!Specific_Sample_Designbool[a]) {
            if(ploidy==2 || TVpars[0].initialDens>1) fprintf(fmoyenne,"Mean within-deme identity prob. (within-individual for cont. pop.) Q0= %12.10e;\n",Qind_moy_glob[a]);
            fprintf(fmoyenne,"Mean between-deme identity prob. (within-pop for cont. pop.) Q1= %12.10e \n",NS_diagnostic_tables::Qb_meanAllPairs_moy_glob[a]); //RL 062018 was Qr_mean_moy_glob, anciennenement mal nommé Qpop_moy_glob
            for(int i=0;i<(/*max(TVpars[0].dimRes1,*/max(dim_sample1[a],dim_sample2[a])/*)*/);i++)
                fprintf(fmoyenne,"Mean identity prob. for genes at %2d steps (Qr%2d) = %12.10e;\n",i,i,Q1_moy_glob[a][i]);
            if(ploidy==2 || TVpars[0].initialDens>1) {
                fprintf(fmoyenne,"Mean intra-deme coalescence time (total sample) = %.8f\n",coa_moy_glob);
                fprintf(fmoyenne,"Mean observed heterozygosity Ho = %12.10e\n",hetero_moy_glob[a]);
            }
            if(HexpNeioui) fprintf(fmoyenne,"Mean genetic diversity He (Nei Expected Heterozygosity) and SE = %12.10e (+- %12.10e))\n",HexpNei_moy_glob[a],1.96*HexpNei_moy_glob_var[a]);
            if(Varoui) {
                fprintf(fmoyenne,"Mean Variance of allelic size and SE = %12.10e (+- %12.10e))\n",Var_moy_glob[a],1.96*Var_moy_glob_var[a]);
                fprintf(fmoyenne,"Mean M statistic (Garza & Williamson 2001) and SE = %12.10e (+- %12.10e))\n",M_moy_glob[a],1.96*M_moy_glob_var[a]);
            }
            fprintf(fmoyenne,"Mean Fst (Fis for cont. pop.) = %12.10e \n",fis_moy_glob[a]);
            fprintf(fmoyenne,"Mean Fst (Fis for cont. pop.) calculated as mean_num/mean_denom = %12.10e\n",fis_moy_numer_glob[a]/fis_moy_denom_glob[a]);
            if (nMarker == 1 && TVpars[0].dimRes1*TVpars[0].dimRes2==1 && (model_mut.compare("IAM")==0 || model_mut.compare("KAM")==0 || model_mut.compare("SMM")==0)) {
                fprintf(fmoyenne,"Theoretical expectations  for a wright fisher populaton of size: %d\n",TVpars[0].initialDens);
                fprintf(fmoyenne,"Genetic diversity Hexp = %12.10e\n",Hexp);
                fprintf(fmoyenne,"Variance in allelic size Vexp = %12.10e\n",Vexp);
            }
        }
    }/*fin if(rep==repet)*/
    fclose(fmoyenne);
    if(suiviQ && !Specific_Sample_Designbool[a]) fclose(fsuivi);
} //fin for a in sampleNb
}

/*---------------------------------------------------------------------------*/

/*---------ecriture des moyennes Fis et Het --------------*/

/*void ecriture_fichier_iterativeStats_cpp() //not used RL, obsolete
{char Fishet[] = "Iterative_Statistics.txt"; // initialisation de longueur suffisante
 int j;

ffishet.open(Fishet,ios::out);

if(rep==1)
{if(!ffishet.is_open())
	{cout<<"Cannot open file: "<<Fishet<<endl<<flush;
//	getchar(); //not good for batch debugging...
	 exit(1);
	}

//ffishet.setf(ios_base::left);
ffishet.width(3);
 if(ploidy==2) for(int j=0;j<n_locus;j++)
  ffishet<<"Het_loc"<<j<<"  ";
 if(HexpNeioui) for(int j=0;j<n_locus;j++)
  ffishet<<"Hexp_loc"<<j<<" ";
 if(DeltaHoui) for(int j=0;j<n_locus;j++)
  ffishet<<"DeltaH_loc"<<j<<" ";
 if(Varoui) for(int j=0;j<n_locus;j++)
  ffishet<<"  Var_loc"<<j<<"    ";
 if(Varoui) for(int j=0;j<n_locus;j++)
  ffishet<<"  MGW_loc"<<j<<"    ";
 for(int j=0;j<n_locus;j++)
  ffishet<<"Fis_loc"<<j<<"  ";
 for(int j=0;j<n_locus;j++)
  ffishet<<"nb_alleles"<<j<<" ";
 if(ploidy==2) ffishet<<"Het_moy     ";
 if(HexpNeioui) ffishet<<"Hexp_moy    ";
 if(DeltaHoui) ffishet<<"DeltaH_moy  ";
 if(Varoui) ffishet<<" Var_moy        ";
 if(Varoui) ffishet<<" MGW_moy        ";
 ffishet<<"fis_moy    ";
 ffishet<<"nb_allele_moy \n";
 //cout<<"zappe un morceau"<<flush;

 ffishet.close();
}//fin if(rep==1)

else

	{
	j=0;
	ffishet.open(Fishet,ios::app);
	while(!(ffishet.is_open()) && j<40){
cout<<"Cannot open file: "<<Fishet<<" ("<<j<<"th attempt)."<<endl<<flush;
		j++;
		mysleep(50); //milliseconds
		//getchar();
	 	//exit(1);
	ffishet.open(Fishet,ios::app);
	}
	if(j>=40){
//		printf("\n ecriture fichiers cannot open file: %sb after 20 tries",Fishet);
//	 	printf("\npress any key to exit");
//		getchar();
 		exit(1);
	}
 cout<<"zappe un autre morceau"<<flush;

 	for(int j=0;j<n_locus;j++)
		fprintf(ffishet,"%8.6e   ",hetero[j]);
	if(HexpNeioui) for(int j=0;j<n_locus;j++)
		fprintf(ffishet,"%8.6e   ",HexpNei[j]);
	if(DeltaHoui) for(int j=0;j<n_locus;j++)
		fprintf(ffishet,"%8.6e     ",DeltaH[j]);
	if(Varoui) for(int j=0;j<n_locus;j++)
		fprintf(ffishet,"% 8.6e   ",Var[j]);
	if(Varoui) for(int j=0;j<n_locus;j++)
		fprintf(ffishet,"% 8.6e   ",M[j]);
	for(int j=0;j<n_locus;j++)
		fprintf(ffishet,"%8.6e   ",fis[j]);
	for(int j=0;j<n_locus;j++)
		fprintf(ffishet,"%10d    ",allele_nmbr[j]);

	fprintf(ffishet,"%8.6e   ",hetero_moy);
	if(HexpNeioui) fprintf(ffishet,"%8.6e   ",HexpNei_moy);
	if(DeltaHoui) fprintf(ffishet,"%8.6e   ",DeltaH_moy);
	if(Varoui) fprintf(ffishet,"%8.6e   ",Var_moy);
	if(Varoui) fprintf(ffishet,"%8.6e   ",M_moy);
	fprintf(ffishet,"%8.6e   ",fis_moy);
	fprintf(ffishet,"%8.6e   \n",mean_allele_nmbr);


	ffishet.close();
}//fin if rep!=1

}
*/

void ecriture_fichier_iterativeStats() {
string Fishet,Fishet2;
FILE* ffishet;
FILE* ffishet2;
int sampleNb=1;

if(predispbool) sampleNb=2;

for(int a=0;a<sampleNb;a++){

    if(a==1) Fishet="Iterative_Statistics_predisp_PerLocus.txt"; else Fishet="Iterative_Statistics_postdisp_PerLocus.txt";
    if(a==1) Fishet2="Iterative_Statistics_predisp_MeanLoc.txt"; else Fishet2="Iterative_Statistics_postdisp_MeanLoc.txt";
    

    if(rep==1) { //entete
        if((ffishet=fopen(Fishet.c_str(),"w"))==NULL) {
            cerr << "ecriture fichiers cannot open file:" << Fishet << endl;
            if(cinGetOnError) cin.get();
            exit(1);
        }
        if((ffishet2=fopen(Fishet2.c_str(),"w"))==NULL) {
            cerr << "ecriture fichiers cannot open file:" << Fishet2 << endl;
            if(cinGetOnError) cin.get();
            exit(1);
        }
        if(!noHeader) {
            //valeurs par locus
            if(Fisoui && (ploidy==2 || TVpars[0].initialDens>1) ) for(int j=0;j<n_locus;j++)
            fprintf(ffishet,"Het_loc%3d  ",j+1);
            if(HexpNeioui) for(int j=0;j<n_locus;j++)
            fprintf(ffishet,"Hexp_loc%3d  ",j+1);
            if(Varoui) for(int j=0;j<n_locus;j++)
            fprintf(ffishet,"Var_loc%3d  ",j+1);
            if(Varoui) for(int j=0;j<n_locus;j++)
            fprintf(ffishet,"MGW_loc%3d  ",j+1);
            if(Fisoui && (ploidy==2 || TVpars[0].initialDens>1)) for(int j=0;j<n_locus;j++)
            fprintf(ffishet,"Fis_loc%3d  ",j+1);
            if(seqStatsOui) {
             for(int j=0;j<n_locus;j++) {
                 if((mutModelIDVector[j+1]=="ISM") || (mutModelIDVector[j+1]=="JC69")
                         || (mutModelIDVector[j+1]=="K80") || (mutModelIDVector[j+1]=="F81")
                         || (mutModelIDVector[j+1]=="HKY85") || (mutModelIDVector[j+1]=="TN93")) {
                     fprintf(ffishet,"SegSitesTotalSample_loc%3d  ",j+1);
                     fprintf(ffishet,"NumPairDiff_loc%3d  ",j+1);
                     if(mutModelIDVector[j+1]!="ISM")
                         fprintf(ffishet,"EmpTITV_ratio_loc%3d  ",j+1);
                 }
             }
            }
            for(int j=0;j<n_locus;j++)
                fprintf(ffishet,"nb_allelesTotalSample_loc%3d  ",j+1);
            
            //moyennes sur locus
            if(Fisoui && (ploidy==2 || TVpars[0].initialDens>1)) {fprintf(ffishet,"Het_moy     ");fprintf(ffishet2,"Het_moy     ");}
            if(HexpNeioui) {fprintf(ffishet,"Hexp_moy  ");fprintf(ffishet2,"Hexp_moy  ");}
            if(Varoui) {fprintf(ffishet,"Var_moy     ");fprintf(ffishet2,"Var_moy     ");}
            if(Varoui) {fprintf(ffishet,"MGW_moy     ");fprintf(ffishet2,"MGW_moy     ");}
            if(Fisoui && (ploidy==2 || TVpars[0].initialDens>1)) {fprintf(ffishet,"fis_moy  ");fprintf(ffishet2,"fis_moy  ");}
            fprintf(ffishet,"nb_allele_moyTotalSample ");fprintf(ffishet2,"nb_allele_moyTotalSample ");
            for(int i=0;i<(/*max(TVpars[0].dimRes1-1,*/max(dim_sample1[a],dim_sample2[a])/*)*//*/vide_sample[a]*/);i++) {
                fprintf(ffishet,"Qr(%2d)  ",i);fprintf(ffishet2,"Qr(%2d)  ",i);
            }
            if(arRegression) {fprintf(ffishet,"ar_slope ar_intercept ");fprintf(ffishet2,"ar_slope ar_intercept ");}
            if(erRegression) {fprintf(ffishet,"er_slope er_intercept ");fprintf(ffishet2,"er_slope er_intercept ");}
            if(moranIRegression) {fprintf(ffishet,"moranI_slope moranI_intercept ");fprintf(ffishet2,"moranI_slope moranI_intercept ");}
            fprintf(ffishet,"\n");
            fprintf(ffishet2,"\n");
        }
    }/*fin if(rep==1)*/
    else {
        int jj=0;
        while((ffishet=fopen(Fishet.c_str(),"a"))==NULL && jj<40){
            printf("\necriture fichiers cannot open file: %s; try=%d\n",Fishet.c_str(),jj);
            jj++;
            mysleep(50);//milliseconds
            //getchar();
            //exit(1);
        }
        if(jj>=40){
            printf("\n ecriture fichiers cannot open file: %sb after 40 tries",Fishet.c_str());
            printf("\npress any key to exit");
            if(cinGetOnError) cin.get();
            exit(1);
        }
        while((ffishet2=fopen(Fishet2.c_str(),"a"))==NULL && jj<40){
            printf("\necriture fichiers cannot open file: %s; try=%d\n",Fishet2.c_str(),jj);
            jj++;
            mysleep(50);//milliseconds
            //getchar();
            //exit(1);
        }
        if(jj>=40){
            printf("\n ecriture fichiers cannot open file: %sb after 40 tries",Fishet2.c_str());
            printf("\npress any key to exit");
            if(cinGetOnError) cin.get();
            exit(1);
        }
    }
    if(Fisoui && (ploidy==2 || TVpars[0].initialDens>1)) for(int j=0;j<n_locus;j++)
        fprintf(ffishet,"%8.6e   ",hetero[a][j]);
    if(HexpNeioui) for(int j=0;j<n_locus;j++)
        fprintf(ffishet,"%8.6e   ",HexpNei[a][j]);
    if(Varoui) for(int j=0;j<n_locus;j++)
        fprintf(ffishet,"%8.6e   ",Var[a][j]);
    if(Varoui) for(int j=0;j<n_locus;j++)
        fprintf(ffishet,"%8.6e   ",M[a][j]);
    if(Fisoui && (ploidy==2 || TVpars[0].initialDens>1)) for(int j=0;j<n_locus;j++)
        fprintf(ffishet,"%8.6e   ",fis[a][j]);
    if(seqStatsOui) {
        for(int j=0;j<n_locus;j++) {
             if((mutModelIDVector[j+1]=="ISM") || (mutModelIDVector[j+1]=="JC69")
                     || (mutModelIDVector[j+1]=="K80") || (mutModelIDVector[j+1]=="F81")
                     || (mutModelIDVector[j+1]=="HKY85") || (mutModelIDVector[j+1]=="TN93")) {
                 fprintf(ffishet,"%10d    ",numSegSites[j]);
                 fprintf(ffishet,"%8.6e    ",numPairMismatch[a][j]);
                 if(mutModelIDVector[j+1]!="ISM")
                     fprintf(ffishet,"%8.6e    ",empTITVRatio[a][j]);
             }
        }
    }
    for(int j=0;j<n_locus;j++)
        fprintf(ffishet,"%8.6e    ", (double) allele_nmbr[j]);

    if(Fisoui && (ploidy==2 || TVpars[0].initialDens>1)) {fprintf(ffishet,"%8.6e   ",hetero_moy[a]);fprintf(ffishet2,"%8.6e   ",hetero_moy[a]);}
    if(HexpNeioui) {fprintf(ffishet,"%8.6e   ",HexpNei_moy[a]);fprintf(ffishet2,"%8.6e   ",HexpNei_moy[a]);}
    if(Varoui) {fprintf(ffishet,"%8.6e   ",Var_moy[a]);fprintf(ffishet2,"%8.6e   ",Var_moy[a]);}
    if(Varoui) {fprintf(ffishet,"%8.6e   ",M_moy[a]);fprintf(ffishet2,"%8.6e   ",M_moy[a]);}
    if(Fisoui && (ploidy==2 || TVpars[0].initialDens>1)) {fprintf(ffishet,"%8.6e  ",fis_moy[a]);fprintf(ffishet2,"%8.6e  ",fis_moy[a]);}
    fprintf(ffishet,"%8.6e   ",mean_allele_nmbr);fprintf(ffishet2,"%8.6e   ",mean_allele_nmbr);
    for(int i=0;i<(/*max(TVpars[0].dimRes1,*/max(dim_sample1[a],dim_sample2[a])/*)*/);i++) {
        fprintf(ffishet,"%8.6e   ",Q1_moy[a][i]);fprintf(ffishet2,"%8.6e   ",Q1_moy[a][i]);
    }
    if(arRegression) {
        fprintf(ffishet,"%8.6e   %8.6e   ",NS_diagnostic_tables::regression_ar[a][0],NS_diagnostic_tables::regression_ar[a][1]);
        fprintf(ffishet2,"%8.6e   %8.6e   ",NS_diagnostic_tables::regression_ar[a][0],NS_diagnostic_tables::regression_ar[a][1]);
    }
    if(erRegression) {
        fprintf(ffishet,"%8.6e   %8.6e   ",NS_diagnostic_tables::regression_er[a][0],NS_diagnostic_tables::regression_er[a][1]);
        fprintf(ffishet2,"%8.6e   %8.6e   ",NS_diagnostic_tables::regression_er[a][0],NS_diagnostic_tables::regression_er[a][1]);
    }
    if(moranIRegression) {
        fprintf(ffishet,"%8.6e   %8.6e   ",NS_diagnostic_tables::regression_moranI[a][0],NS_diagnostic_tables::regression_moranI[a][1]);
        fprintf(ffishet2,"%8.6e   %8.6e   ",NS_diagnostic_tables::regression_moranI[a][0],NS_diagnostic_tables::regression_moranI[a][1]);
    }
    fprintf(ffishet,"\n");fclose(ffishet);
    fprintf(ffishet2,"\n");fclose(ffishet2);
}//end a<sampleNb
}
/*****************************************************************************/

/*-------------------ecriture fichiers Mig-------------------------------*/
void ecriture_fichier_MIG(string Fstat_str) {
string migFileName;
ofstream fmig;
int sampleNb=1,pairNb;
    
if(predispbool) sampleNb=2;

for(int a=0;a<sampleNb;a++){
    
    size_t startGene=1,endGene=n_genes[0],indPairNb=n_genes[0]/ploidy*(n_genes[0]/ploidy-1)/2,popPairNb;
    if(a==1) {startGene=n_genes[0]+1;endGene=n_genes_total;indPairNb=n_genes[1]/ploidy*(n_genes[1]/ploidy-1)/2;}
    
    //computation of geographic distances between all sampled population pairs
    if(!Specific_Sample_Designbool[a]) popPairNb=dim_sample1[a]*dim_sample2[a]*(dim_sample1[a]*dim_sample2[a]-1)/2;
    else popPairNb=Spec_SampleSize[a]*(Spec_SampleSize[a]-1)/2;

    if(cmp_nocase(Fstat_str, "ar") || cmp_nocase(Fstat_str, "er")) pairNb=indPairNb; else pairNb=popPairNb;

    stringstream stst,stst2;
    if(a>0) stst2 << "_predisp_" << Fstat_str; else stst2 << "_" << Fstat_str;;
    

    if (repet!=1) {
        stst << fichier_genepop << rep << stst2.str() << ".MIG";
    } else stst << fichier_genepop << stst2.str() << ".MIG";
    migFileName=stst.str();
    
    fmig.open(migFileName.c_str(), std::ios::out);
    
    if (!fmig.is_open()) {
        cerr << "Could not open the file " << migFileName << endl;
        cerr << "\nAborting IBDSim...Press any key to exit.\n";
        if(cinGetOnError)
            cin.get();
        exit(-1);
    }
    
    fmig << "This MIG file was generated by IBDSim" << endl;
    fmig << (endGene-startGene+1)/ploidy << " populations/individuals, thus " <<  pairNb << "pairs;" << endl;
    fmig << "Genetic statistic (" << Fstat_str << "):" << endl;
    
    size_t pair=0;
    //loop over each individuals
    for(size_t i=startGene+ploidy;i<=(endGene-ploidy+1);i+=ploidy) {
        //loop over "all<" other individuals consider all individual pair
        for(size_t ii=startGene;ii<i;ii+=ploidy) {
            if(cmp_nocase(Fstat_str, "ar")==0)
                fmig << fixed << setprecision(15) << NS_diagnostic_tables::Fstat_ar_fromSS[a][pair] << " ";
            if(cmp_nocase(Fstat_str, "er")==0)
                fmig << fixed << setprecision(15) << NS_diagnostic_tables::Fstat_er_fromQ[a][pair] << " ";
            if(cmp_nocase(Fstat_str, "moranI")==0)
                fmig << fixed << setprecision(15) << NS_diagnostic_tables::Fstat_moranI_fromQ_forEachIndPairs[a][pair] << " ";
            pair++;
        }
        fmig << endl;
    }

    fmig << "distances:" << endl;
    
    pair=0;
    //loop over each individuals
    for(size_t i=startGene+ploidy;i<=(endGene-ploidy+1);i+=ploidy) {
        //loop over "all<" other individuals consider all individual pair
        for(size_t ii=startGene;ii<i;ii+=ploidy) {
            fmig << setprecision(15) << NS_diagnostic_tables::indGeoDist[a][pair] << " ";
            pair++;
        }
        /*if(i<(endGene-ploidy+1))*/ fmig << endl;
    }

    fmig.close();
}
}



/*************calcul du temps de coalescence ***********
 ****************pour deux genes d'un individu**************/

void calcul_coa_individuel(void)
{int j,k,o,coa,mauvais;
// int noeud_gene1[1601],noeud_gene2[1601];

vector<int>noeud_gene1(2*n_genes_total+1); //FR 08/2005 le +1 est important
vector<int>noeud_gene2(2*n_genes_total+1);

for(int i=1;i<n_genes_total;i+=2)
 {coa=0;
 mauvais=0;
 retour:
  for(int oo=1;oo<=2*n_genes_total;oo++)
	{noeud_gene1[oo]=0;
	 noeud_gene2[oo]=0;
	}
  j=2;
  k=2;
  //stock la lignÈe du gene i
  noeud_gene1[1]=i;
  noeud_gene1[2]=noeud[i][locus].ancetre;
  while(noeud[noeud_gene1[j]][locus].ancetre>0){
	 j++;
	 noeud_gene1[j]=noeud[noeud_gene1[j-1]][locus].ancetre;
	 }
  //stock la lignÈe du gene i+1
  noeud_gene2[1]=i+1;
  noeud_gene2[2]=noeud[i+1][locus].ancetre;
  while(noeud[noeud_gene2[k]][locus].ancetre>0){
	 k++;
	 noeud_gene2[k]=noeud[noeud_gene2[k-1]][locus].ancetre;
	 }
  //affichage des lignÈes
  if(mauvais==1){
  printf("\nn_genes_total= %d  ;2*n_genes_total= %d   ;\n",n_genes_total,2*n_genes_total);
  for(int oo=1;oo<=2*n_genes_total;oo++){
   printf("\nancetre de %d= %d",oo,noeud[oo][locus].ancetre);
   getchar();
   }
  printf("\nlignÈe gene %d  ;lignÈe gene %d   ;\n",i,i+1);
  o=1;
  while((noeud_gene1[o]!=0)&&(noeud_gene2[o]!=0))
    {printf("    %d       ;   %d   ;\n",noeud_gene1[o],noeud_gene2[o]);
    getchar();
    o++;
    }
  }
  //comparaison des  deux lignÈes pour trouver le premier ancetre commun
  for(int oo=1;oo<=2*n_genes_total;oo++)
	{if(noeud_gene1[oo]==0) break;
	 for(int m=1;m<=2*n_genes_total;m++)
	  {if((coa==1)||(noeud_gene2[m]==0)) break;
		if(noeud_gene1[oo]==noeud_gene2[m])
			{temps_coa+=noeud[noeud_gene1[oo]][locus].generation;
			 compteur3+=1;
			 coa=1;
			 //printf("coa apres %d generations pour gene %d, %d;\n"
			 //			  ,temps_coa,i,i+1);
			 //printf("nbre de coalescence individuelles %d\n",compteur3);
			 //getchar();
			}
	  }
	}
	if(coa==0)
	{printf("\n\n\n\n !!probleme de calcul de temps de coa intra individu au locus: %d\n",locus+1);
	  printf("\nn_genes_total= %d  ;2*n_genes_total= %d   ;\n",n_genes_total,2*n_genes_total);
	  printf("\nlignée gene %d  ;lignée gene %d   ;\n",i,i+1);
	  o=1;
	  while((noeud_gene1[o]!=0)&&(noeud_gene2[o]!=0))
	    {printf("    %d       ;   %d   ;\n",noeud_gene1[o],noeud_gene2[o]);
	    getchar();
	    o++;
	    }
		getchar();
		mauvais=1;
		goto retour;
	  /*exit(1);*/
	}
 }
if(locus==(n_locus-1))
{temps_coa_moy=(double) temps_coa/compteur3;
coa_moy_glob+=temps_coa_moy/repet;
}
//free_ivector(noeud_gene1);
//free_ivector(noeud_gene2);
}

/*****************************************************************************/

/******************ecriture des distances de dispersion "efficaces"***********/
void ecriture_disp(void) {
    using namespace NS_diagnostic_tables;
    using namespace NS_translation;

    string Immig,disp;
    /*char disp[50]="EmpDisp_",ch1[100];
char Immig[50]="EmpImmigRate_";*/
	int i,j;

    disp="EmpDisp_"+fichier_genepop;
    Immig="EmpImmigRate_"+fichier_genepop;
    if (repet!=1) {stringstream stst;
        stst<<rep;
        disp+=stst.str();
        Immig+=stst.str();
    }
    if(txt_extensionbool) {disp+=".txt";Immig+=".txt";}

ofstream fImmig(Immig.c_str(),ios::out);
if(!fImmig.is_open())
{printf("ecriture fichiers cannot open file: %s\n",Immig.c_str());
	getchar();
	exit(1);
}
for(j=dim_reseau2;j>=1;j--) {
	for(i=1;i<dim_reseau1+1;i++) if(i!=dim_reseau1) {
		fImmig << (float) effectiveImmigration[i][j]/cumulEffectiveImmigration[i][j] << " ";
		//cout << effectiveImmigration[i][j] << " / " << cumulEffectiveImmigration[i][j] << "-> " << (float) effectiveImmigration[i][j]/cumulEffectiveImmigration[i][j] << endl;
		//getchar();
		effectiveImmigration_moy[i][j]+=effectiveImmigration[i][j];
		cumulEffectiveImmigration_moy[i][j]+=cumulEffectiveImmigration[i][j];
	} else {
		fImmig << (float) effectiveImmigration[i][j]/cumulEffectiveImmigration[i][j] << endl;
		//cout << effectiveImmigration[i][j] << " / " << cumulEffectiveImmigration[i][j] << "-> " << (float) effectiveImmigration[i][j]/cumulEffectiveImmigration[i][j] << endl;
		//getchar();
		effectiveImmigration_moy[i][j]+=effectiveImmigration[i][j];
		cumulEffectiveImmigration_moy[i][j]+=cumulEffectiveImmigration[i][j];
	}
}

fImmig.close();


if((fdisp=fopen(disp.c_str(),"w"))==NULL)
	{printf("ecriture fichiers cannot open file: %s\n",disp.c_str());
	getchar();
	 exit(1);
	}
for(i=0;i<=2*grande_dim;i++)
  {fprintf(fdisp," %4d  %8ld %8ld %12.10e %12.10e \n",i-grande_dim,(long int) effective[i][x],(long int) effective[i][y],deffective[i][x],deffective[i][y]);
  }

fprintf(fdisp,"Various statistics on the whole empirical distribution (non-central moments)\n");
fprintf(fdisp,"mean = %12.10e %12.10e \n",moy_dispx,moy_dispy);
fprintf(fdisp,"non-central sigma≤  = %12.10e %12.10e \n",sig_dispx,sig_dispy);
fprintf(fdisp,"non-central kurtosis= %12.10e %12.10e \n",kurto_dispx,kurto_dispy);
fprintf(fdisp,"non-central skewness= %12.10e %12.10e \n",skew_dispx,skew_dispy);
fprintf(fdisp,"Various statistics on the semi empirical distribution (one side, axial dispersal, central moments)\n");
fprintf(fdisp,"mean axial dispersal = %12.10e %12.10e \n",moy_demidispx,moy_demidispy);
fprintf(fdisp,"central/axial sigma≤  = %12.10e %12.10e \n",sig_demidispx,sig_demidispy);
fprintf(fdisp,"central/axial kurtosis= %12.10e %12.10e \n",kurto_demidispx,kurto_demidispy);
fprintf(fdisp,"central/axial skewness= %12.10e %12.10e",skew_demidispx,skew_demidispy);

fclose(fdisp);
}

void ecriture_disp_moy()
{
using namespace NS_diagnostic_tables;
    using namespace NS_translation;
string disp_moy="MeanEmpDisp_"+fichier_genepop+".txt";
string Immig_moy="MeanEmpImmigRate_"+fichier_genepop+".txt";

ofstream fImmig_moy(Immig_moy,ios::out);
if(!fImmig_moy.is_open())
{printf("ecriture fichiers cannot open file: %s\n",Immig_moy.c_str());
	getchar();
	exit(1);
}
for(int j=dim_reseau2;j>=1;j--) {
	for(int i=1;i<dim_reseau1+1;i++) if(i!=dim_reseau1) fImmig_moy << (float) effectiveImmigration_moy[i][j]/cumulEffectiveImmigration_moy[i][j] << " ";
	else fImmig_moy << (float) effectiveImmigration_moy[i][j]/cumulEffectiveImmigration_moy[i][j] << endl;
}
fImmig_moy.close();



if((fdisp_moy=fopen(disp_moy.c_str(),"w"))==NULL)
	{printf("ecriture fichiers cannot open file: %s\n",disp_moy.c_str());
	getchar();
	 exit(1);
	}

for(int i=0;i<=2*grande_dim;i++)
  {fprintf(fdisp_moy," %4d %8ld %8ld %12.10e %12.10e \n",i-grande_dim,(long int) effective_moy[i][x],(long int)effective_moy[i][y],((double) effective_moy[i][x])/((double) cum_moyx), ((double) effective_moy[i][y]) / ((double) cum_moyy) );
  }
fprintf(fdisp_moy,"Various statistics on the whole empirical distribution (non-central moments)\n");
fprintf(fdisp_moy,"mean = %12.10e %12.10e \n",moy_dispx,moy_dispy);
fprintf(fdisp_moy,"non-central sigma2 = %12.10e %12.10e \n",sig_dispx,sig_dispy);
fprintf(fdisp_moy,"non-central kurtosis= %12.10e %12.10e \n",kurto_dispx,kurto_dispy);
fprintf(fdisp_moy,"non-central skewness= %12.10e %12.10e \n",skew_dispx,skew_dispy);
fprintf(fdisp_moy,"Various statistics on the semi empirical distribution (one side, axial dispersal, central moments)\n");
fprintf(fdisp_moy,"mean axial dispersal = %12.10e %12.10e \n",moy_demidispx,moy_demidispy);
fprintf(fdisp_moy,"central/axial sigma2 = %12.10e %12.10e \n",sig_demidispx,sig_demidispy);
fprintf(fdisp_moy,"central/axial kurtosis= %12.10e %12.10e \n",kurto_demidispx,kurto_demidispy);
fprintf(fdisp_moy,"central/axial skewness= %12.10e %12.10e",skew_demidispx,skew_demidispy);
fclose(fdisp_moy);
}


/*****************************************************************************/
/*******************************************************************************/
void calcul_proba_identite_pour_distances_axiales(void)
{/*totalement changee le 10072002,
  mais toujours pas en diagonale
  modifs le 13082002 pour avpoir Qindividu
  */
using namespace NS_coal_tree;

int 	compteur_ind,Qind=0,_pas,nSample=1;
//int compteur[2000],identity[2000];
//int *compteur,*identity;

fflush(stdout);
fflush(stdin);

if(predispbool) nSample=2;
for(int a=0;a<nSample;a++){

    int startIndex=1,endIndex=n_genes[0];
    if(a==1) {startIndex=n_genes[0]+1;endIndex=n_genes_total;}

    for(int iLoc=0;iLoc<n_locus;iLoc++) {//pour chaque locus
        if(ploidy==2 || TVpars[0].initialDens>1) {Qind=0;compteur_ind=0;}
        vector<int>compteur(max(TVpars[0].dimRes1,max(dim_sample1[a],dim_sample2[a]))+1,0);
        vector<int>identity(max(TVpars[0].dimRes1,max(dim_sample1[a],dim_sample2[a]))+1,0);
        for(int d1=startIndex+1;d1<=endIndex;d1++)//on compare tous les genes de l'echantillon deux a deux
            for(int d2=startIndex;d2<d1;d2++) {//calcul sur dim1
                if(coord_individus[d1][y]==coord_individus[d2][y]) {
                     _pas=(coord_individus[d1][x]-coord_individus[d2][x])/vide_sampleX[a];
                     compteur[_pas]++;
                     if(noeud[d1][iLoc].etat_allele==noeud[d2][iLoc].etat_allele) identity[_pas]++;
                     //printf("\n dim1 1compteur[%d]=%d,identity[%d]=%d",_pas,compteur[_pas],_pas,identity[_pas]);
                     //getchar();
                     if( (ploidy==2 || TVpars[0].initialDens>1) &&(coord_individus[d1][x]==coord_individus[d2][x])&&((d1 % 2)==0))//car un ind=(impair,pair)
                        {compteur_ind++;
                         if(noeud[d1][iLoc].etat_allele==noeud[d2][iLoc].etat_allele) Qind++;
                         //printf("\n dim1 compteur_ind=%d,Qind=%d",compteur_ind,Qind);
                         //getchar();
                        }
                     if(TVpars[0].dimRes1!=vide_sampleX[a]
                            && int(SortieInffnPtr1(-_pas*vide_sampleX[a],TVpars[0].dimRes1)/vide_sampleX[a])!=_pas && EdgeEffect=="circ")//distance en faisant le tour du torre
                       {_pas=int(SortieInffnPtr1(-_pas*vide_sampleX[a],TVpars[0].dimRes1))/vide_sampleX[a];
                        compteur[_pas]++;
                        if(noeud[d1][iLoc].etat_allele==noeud[d2][iLoc].etat_allele) identity[_pas]++;
                        //printf("\n dim1 2compteur[%d]=%d,identity[%d]=%d",_pas,compteur[_pas],_pas,identity[_pas]);
                        //getchar();
                       }
                  }
                  //else//calcul sur dim2
                  if(coord_individus[d1][x]==coord_individus[d2][x])
                    {_pas=(coord_individus[d1][y]-coord_individus[d2][y])/vide_sampleY[a];
                     if(_pas!=0)
                        {compteur[_pas]++;
                         if(noeud[d1][iLoc].etat_allele==noeud[d2][iLoc].etat_allele) identity[_pas]++;
                         //printf("\n dim2 1compteur[%d]=%d,identity[%d]=%d",_pas,compteur[_pas],_pas,identity[_pas]);
                         //getchar();
                        }
                     if(TVpars[0].dimRes2!=vide_sampleY[a] && int(SortieInffnPtr1(-_pas*vide_sampleY[a],TVpars[0].dimRes2))/vide_sampleY[a]!=_pas
                                       && EdgeEffect=="circ")//distance en faisant le tour du torre
                        {_pas=int(SortieInffnPtr1(-_pas*vide_sampleY[a],TVpars[0].dimRes2))/vide_sampleY[a];
                         compteur[_pas]++;
                         if(noeud[d1][iLoc].etat_allele==noeud[d2][iLoc].etat_allele) identity[_pas]++;
                         //printf("\n dim2 2compteur[%d]=%d,identity[%d]=%d",_pas,compteur[_pas],_pas,identity[_pas]);
                         //getchar();
                        }
                    }
                }//fin boucle sur d2
         for(_pas=0;_pas<max(dim_sample1[a],dim_sample2[a]);_pas++)
            {if(identity[_pas]>0 )
                 {//printf("\nfinal compteur[%d]=%d,identity[%d]=%d",_pas,compteur[_pas],_pas,identity[_pas]);
                  //getchar();
                  Q1[a][iLoc][_pas]=(double) identity[_pas]/compteur[_pas];
                  //printf("\nfinal Q1[%d][loc=%d][_pas=%d]=%f",a,iLoc,_pas,Q1[a][iLoc][_pas]);
                  //getchar();
                  Q1_moy[a][_pas]+=Q1[a][iLoc][_pas]/n_locus;
                 }
            }
        if(ploidy==2 || TVpars[0].initialDens>1) if(Qind>0) {
            Qind2[a][iLoc]= (double) Qind/compteur_ind;
            Qind_moy[a]+=Qind2[a][iLoc]/n_locus;
            //cout << "Qind2=" << Qind2[a][iLoc] << "=Qind/compteur_ind=" << Qind << "/" << compteur_ind << endl;
        }

        }//fin boucle sur iLoc

    for(_pas=0;_pas<max(dim_sample1[a],dim_sample2[a]);_pas++) if(Q1_moy[a][_pas]!=0.0) Q1_moy_glob[a][_pas]+=Q1_moy[a][_pas]/repet;
    if(ploidy==2 || TVpars[0].initialDens>1) if(Qind_moy[a]!=0.0) Qind_moy_glob[a]+=Qind_moy[a]/repet;
}
}

/*****************************************************************************/
/*******************************************************************************/
void calcul_proba_identite_par_paires_individus(void)
{
using namespace NS_coal_tree;
using namespace NS_diagnostic_tables;


size_t nSample=1,indNb,genePairPerIndPair=4;
vector<size_t> comptPairDistClass;


if(ploidy==1) genePairPerIndPair=1;
double inc=1.0/(1.0*genePairPerIndPair);
if(predispbool) nSample=2;
Qb_pair.resize(nSample);Qw_pair.resize(nSample);Qw_ind.resize(nSample);QiOverAllOtherInd.resize(nSample);
Qb_meanAllPairs.resize(nSample);Qw_meanAllInd.resize(nSample);Qb_distClass.resize(nSample);

Qb_pair_moy.resize(nSample);Qw_pair_moy.resize(nSample);Qw_ind_moy.resize(nSample);QiOverAllOtherInd_moy.resize(nSample);
Qb_meanAllPairs_moy.resize(nSample);Qw_meanAllInd_moy.resize(nSample);Qb_distClass_moy.resize(nSample);

for(int a=0;a<nSample;a++){
    
    size_t startGene=1,endGene=n_genes[0],pairNb=n_genes[0]/ploidy*(n_genes[0]/ploidy-1)/2;
    if(a==1) {startGene=n_genes[0]+1;endGene=n_genes_total;pairNb=n_genes[1]/ploidy*(n_genes[1]/ploidy-1)/2;}
    indNb=(endGene-startGene+1)/ploidy;


    Qb_pair[a].resize(n_locus);Qw_pair[a].resize(n_locus);Qw_ind[a].resize(n_locus);QiOverAllOtherInd[a].resize(n_locus);
    Qb_meanAllPairs[a].resize(n_locus,0.0);Qw_meanAllInd[a].resize(n_locus,0.0);Qb_distClass[a].resize(n_locus);

    Qb_pair_moy[a].resize(pairNb,0.0);Qw_pair_moy[a].resize(pairNb,0.0);Qw_ind_moy[a].resize(pairNb,0.0);QiOverAllOtherInd_moy[a].resize(indNb);
    Qb_distClass_moy[a].resize((int) floor(maxDistSample[a]) + 1,0.0);
    
    for(size_t iLoc=0;iLoc<n_locus;iLoc++) {//pour chaque locus
        Qb_pair[a][iLoc].resize(pairNb,0.0);Qw_pair[a][iLoc].resize(pairNb,0.0);Qw_ind[a][iLoc].resize(indNb,0.0);QiOverAllOtherInd[a][iLoc].resize(indNb);
        Qb_distClass[a][iLoc].resize((int) floor(maxDistSample[a]) + 1,0.0);comptPairDistClass.resize((int) floor(maxDistSample[a]) + 1,0);

        
        size_t pair=0;
        //loop over each individuals
        for(size_t i=startGene+ploidy;i<=(endGene-ploidy+1);i+=ploidy) {
            //loop over "all<" other individuals consider all individual pair
            for(size_t ii=startGene;ii<i;ii+=ploidy) {
                //loop over each gene of each individual
                for(size_t j=0;j<ploidy;j++) {
                    for(size_t k=0;k<ploidy;k++) {
                        if(noeud[i+j][iLoc].etat_allele==noeud[ii+k][iLoc].etat_allele) Qb_pair[a][iLoc][pair]+=inc;
                        //cout << i+j << "vs" << ii+k << " -> allele=" << noeud[i+j][iLoc].etat_allele << "&" << noeud[ii+k][iLoc].etat_allele << endl;
                        //cout << "Qb_pair[a][iLoc][pair]=" << Qb_pair[a][iLoc][pair] << endl;
                    }
                }
                
                if(noeud[i][iLoc].etat_allele==noeud[i+1][iLoc].etat_allele) Qw_pair[a][iLoc][pair]+=1.0/2.0;
                if(noeud[ii][iLoc].etat_allele==noeud[ii+1][iLoc].etat_allele) Qw_pair[a][iLoc][pair]+=1.0/2.0;
                Qw_pair_moy[a][pair]+=Qw_pair[a][iLoc][pair]/n_locus;
                
                Qb_meanAllPairs[a][iLoc]+=Qb_pair[a][iLoc][pair]/pairNb;
                Qb_pair_moy[a][pair]+=Qb_pair[a][iLoc][pair]/n_locus;
                Qb_meanAllPairs_moy[a]+=Qb_pair[a][iLoc][pair]/pairNb/n_locus;
                
                QiOverAllOtherInd[a][iLoc][(i-1)/ploidy]+=Qb_pair[a][iLoc][pair]/(indNb-1);
                QiOverAllOtherInd[a][iLoc][(ii-1)/ploidy]+=Qb_pair[a][iLoc][pair]/(indNb-1);
                QiOverAllOtherInd_moy[a][(i-1)/ploidy]+=Qb_pair[a][iLoc][pair]/(indNb-1)/n_locus;
                QiOverAllOtherInd_moy[a][(ii-1)/ploidy]+=Qb_pair[a][iLoc][pair]/(indNb-1)/n_locus;
                
                if(ii==startGene) { // pour boucle sur individu = une seule fois par ind (i-1)/ploidy
                    if(ploidy==2 && noeud[i][iLoc].etat_allele==noeud[i+1][iLoc].etat_allele) {
                        Qw_ind[a][iLoc][(i-1)/ploidy]=1.0;
                        Qw_ind_moy[a][(i-1)/ploidy]+=1.0/n_locus;
                        Qw_meanAllInd[a][iLoc]+=1.0/(indNb);
                        Qw_meanAllInd_moy[a]+=1.0/(indNb)/n_locus;
                    }
                    //et une fois pour l'individu (startGene-1)/ploidy qui n'est pas dans la boucle i
                    if(i==(startGene+ploidy) ) {
                        if(ploidy==2 && noeud[ii][iLoc].etat_allele==noeud[ii+1][iLoc].etat_allele ) {
                            Qw_ind[a][iLoc][(ii-1)/ploidy]=1.0;
                            Qw_ind_moy[a][(ii-1)/ploidy]+=Qw_ind[a][iLoc][(ii-1)/ploidy]/n_locus;
                            Qw_meanAllInd[a][iLoc]+=1.0/(indNb);
                            Qw_meanAllInd_moy[a]+=1.0/(indNb)/n_locus;
                        }
                    }
                }
                size_t dist=floor(sqrt( sq( (double) (coord_individus[i][x] - coord_individus[ii][x]) )
                                       + sq( (double) (coord_individus[i][y] - coord_individus[ii][y]) ) ) );
                Qb_distClass[a][iLoc][dist]+=Qb_pair[a][iLoc][pair];
                Qb_distClass_moy[a][dist]+=Qb_pair[a][iLoc][pair]/n_locus;

                comptPairDistClass[dist]++;
                pair++;

            }
        }
        for(size_t d=0;d<=(int) floor(maxDistSample[a]);d++) if(comptPairDistClass[d] != 0) {
            Qb_distClass[a][iLoc][d]/=1.0*comptPairDistClass[d];
            Qb_distClass_moy[a][d]/=1.0*comptPairDistClass[d];
        } else {
            Qb_distClass[a][iLoc][d]=0.0;
            Qb_distClass_moy[a][d]=0.0;

        }
    }//fin boucle sur iLoc
    
    
    //affichage debug
//    cout << endl << endl;
//    for(size_t iLoc=0;iLoc<n_locus;iLoc++) {
//        size_t pair=0;
//        for(size_t i=startGene+ploidy;i<=(endGene-ploidy+1);i+=ploidy) {
//            for(size_t ii=startGene;ii<i;ii+=ploidy) {
//                cout << "p=" << pair << "->";
//                cout << "Qb=" << Qb_pair[a][iLoc][pair] << "; Qw=" << Qw_pair[a][iLoc][pair] << endl;
//                pair++;
//            }
//        }
//            
//        cout << endl;
//        for(size_t i=startGene;i<=(endGene-ploidy+1);i+=ploidy) {
//            cout << "ind=" << (i-1)/ploidy << "->";
//            cout << "Qw=" << Qw_ind[a][iLoc][(i-1)/ploidy] << "; ";
//            cout << "QiOverAllOtherInd=" << QiOverAllOtherInd[a][iLoc][(i-1)/ploidy] << endl;
//        }
//        cout << endl;
//        cout << "Qw_meanAllInd=" << Qw_meanAllInd[a][iLoc] << endl;
//        cout << "Qb_meanAllPairs=" << Qb_meanAllPairs[a][iLoc] << endl;
//
//
//        for(size_t d=0;d<=(int) floor(maxDistSample[a]);d++) {
//            cout << "d=" << d << "->" << Qb_distClass[a][iLoc][d] << endl;;
//        }
//    }//fin boucle sur iLoc
//
//    cout << endl << endl;
//    size_t pair=0;
//    for(size_t i=startGene+ploidy;i<=(endGene-ploidy+1);i+=ploidy) {
//        for(size_t ii=startGene;ii<i;ii+=ploidy) {
//            cout << "p=" << pair << "->";
//            cout << "Qb_moy=" << Qb_pair_moy[a][pair] << endl;
//            pair++;
//        }
//    }
//    
//    cout << endl;
//    for(size_t i=startGene;i<=(endGene-ploidy+1);i+=ploidy) {
//        cout << "ind=" << (i-1)/ploidy << "->";
//        cout << "Qw_moy=" << Qw_ind_moy[a][(i-1)/ploidy] << "; ";
//        cout << "QiOverAllOtherInd_moy=" << QiOverAllOtherInd_moy[a][(i-1)/ploidy] << endl;
//    }
//    cout << endl;
//    cout << "Qw_meanAllInd_moy=" << Qw_meanAllInd_moy[a] << endl;
//    cout << "Qb_meanAllPairs_moy=" << Qb_meanAllPairs_moy[a] << endl;
//    
//    
//    for(size_t d=0;d<=(int) floor(maxDistSample[a]);d++) {
//        cout << "d=" << d << "->Qb_distClass_moy=" << Qb_distClass_moy[a][d] << endl;;
//    }

}//fin boucle sur sample a
}

/*****************************************************************************/


/*******************************************************************************/
void calcul_proba_identite_matrice(void)
{/*created 12/2005, adapted from calcul_proba_identitÈ_par_pas()*/
//calcul les probilitÈ d'identitÈ de paires de genes intra individus compris
//facilement modifiable si necessaire
using namespace NS_coal_tree;

int iLoc,i,j,i2,j2,d1,d2;
int sampleNb=1;

if(predispbool) sampleNb=2;

for(int a=0;a<sampleNb;a++){

    vector< vector< vector< vector<int> > > > compteur;
    vector< vector< vector< vector<int> > > > identity;

    int startIndex=1,endIndex=n_genes[0];
    if(a==1) {startIndex=n_genes[0]+1;endIndex=n_genes_total;}

    //dimensionement des matrices locales
    compteur.resize(dim_sample1[a]);
    identity.resize(dim_sample1[a]);
    for(i=0;i<dim_sample1[a];i++) {
        compteur[i].resize(dim_sample2[a]);
        identity[i].resize(dim_sample2[a]);
        for(j=0;j<dim_sample2[a];j++) {
            compteur[i][j].resize(dim_sample1[a]);
            identity[i][j].resize(dim_sample1[a]);
            for(i2=0;i2<dim_sample1[a];i2++) {
                compteur[i][j][i2].resize(dim_sample2[a],0);
                identity[i][j][i2].resize(dim_sample2[a],0);
            }
        }
    }


    for(iLoc=0;iLoc<n_locus;iLoc++) {//pour chaque locus
        //initialisation des matrices locales
        for(i=0;i<dim_sample1[a];i++)
            for(j=0;j<dim_sample2[a];j++)
                for(i2=0;i2<dim_sample1[a];i2++)
                    for(j2=0;j2<dim_sample2[a];j2++) {
                        compteur[i][j][i2][j2]=0;
                        identity[i][j][i2][j2]=0;
                    }

        //cout << "\n\n\n\n**For locus=" << iLoc+1 << endl;
        for(d1=startIndex+1;d1<=endIndex;d1++) //on compare tous les genes de l'echantillon deux a deux
            for(d2=startIndex;d2<d1;d2++) {
                 //cout << "\n\nComparison between genes " << d1 << "(" << coord_individus[d1][x] << "," << coord_individus[d1][y];
                 //cout <<  ") and " << d2 << "(" << coord_individus[d2][x] << "," << coord_individus[d2][y] << ")"<< endl ;
                 //cout << "new coord genes " << d1 << "(" << (coord_individus[d1][x]-xmin_sample)/vide_sampleX << "," << (coord_individus[d1][y]-ymin_sample)/vide_sampleY;
                 //cout <<  ") and " << d2 << "(" << (coord_individus[d2][x]-xmin_sample)/vide_sampleX << "," << (coord_individus[d2][y]-ymin_sample)/vide_sampleY << ")"<< endl ;
                 compteur[(coord_individus[d1][x]-xmin_sample[a])/vide_sampleX[a]][(coord_individus[d1][y]-ymin_sample[a])/vide_sampleY[a]][(coord_individus[d2][x]-xmin_sample[a])/vide_sampleX[a]][(coord_individus[d2][y]-ymin_sample[a])/vide_sampleY[a]]++;
                 //cout << "compteur=" << compteur[(coord_individus[d1][x]-xmin_sample)/vide_sample][(coord_individus[d1][y]-ymin_sample)/vide_sample][(coord_individus[d2][x]-xmin_sample)/vide_sample][(coord_individus[d2][y]-ymin_sample)/vide_sample];
                 if(noeud[d1][iLoc].etat_allele==noeud[d2][iLoc].etat_allele)
                     identity[(coord_individus[d1][x]-xmin_sample[a])/vide_sampleX[a]][(coord_individus[d1][y]-ymin_sample[a])/vide_sampleY[a]][(coord_individus[d2][x]-xmin_sample[a])/vide_sampleX[a]][(coord_individus[d2][y]-ymin_sample[a])/vide_sampleY[a]]++;
            }//fin boucle sur d2

        for(i=0;i<dim_sample1[a];i++)
            for(j=0;j<dim_sample2[a];j++)
                for(i2=0;i2<dim_sample1[a];i2++)
                    for(j2=0;j2<dim_sample2[a];j2++) {

                        if(compteur[i][j][i2][j2]>0) QMatrix[a][iLoc][i][j][i2][j2]= (double) identity[i][j][i2][j2]/compteur[i][j][i2][j2];
                        else QMatrix[a][iLoc][i][j][i2][j2]=-1.0;

                        QMatrix_Moy[a][i][j][i2][j2]+=QMatrix[a][iLoc][i][j][i2][j2]/n_locus;
                        //if(QMatrix_Moy[a][i][j][i2][j2]>=0.0) cout << "\n<< a << ", " << i << ", " << j << ", " << i2 << ", " << j2 << "->" << QMatrix_Moy[a][i][j][i2][j2]*n_locus;
                    }

    }//fin boucle sur iLoc

    //computation of the global stat = pour toutes mes repets
    //cout << endl;
    //cout << endl;
    for(i=0;i<dim_sample1[a];i++)
        for(j=0;j<dim_sample2[a];j++)
            for(i2=0;i2<dim_sample1[a];i2++)
                for(j2=0;j2<dim_sample2[a];j2++) {
                    //if(QMatrix_Moy[a][i][j][i2][j2]>=0.0) cout << "\n" << a << ", " << i << "," << j << "," << i2 << "," << j2 << "->" <<QMatrix_Moy[a][i][j][i2][j2];
                    QMatrix_Moy_Glob[a][i][j][i2][j2]+=QMatrix_Moy[a][i][j][i2][j2]/repet;
                }
}// for a<sampleNb
}

/*****************************************************************************/

/******************ecriture Matrice Proba d'identitÈ**************************/
void ecriture_matrice_proba_id(void) {
string FileName="Matrix_ProbId";
int i,j,i2,j2,sampleNb=1;//,count;
ofstream File;


if(predispbool) sampleNb=2;

for(int a=0;a<sampleNb;a++){

    if(a==1) FileName+="_predisp";
    FileName+=".txt";
    File.open(FileName.c_str(),ios::out);
    if(!File.is_open()) {
        cout << "Cannot open file :" << FileName << endl << flush;
        getchar();
        exit(1);
    }
    //ecriture "colone" x1,y1,x2,y2,Q((x1,y1),(x2,y2))  (x1 = coord en x du gene 1... y2 = coord en y du gene2)
    //count=0;
    //cout << endl;
    //cout << endl;
    for(i=0;i<dim_sample1[a];i++)
        for(j=0;j<dim_sample2[a];j++)
            for(i2=0;i2<dim_sample1[a];i2++)
                for(j2=0;j2<dim_sample2[a];j2++) {
                    if(QMatrix_Moy_Glob[a][i][j][i2][j2]>=0.0) { //FER
                        //count++;
                        //cout << endl;
                        //cout << count << ";" <<  i*vide_sample[a]+xmin_sample[a] << "," << j*vide_sample[a]+ymin_sample[a] << "," << i2*vide_sample[a]+xmin_sample[a];
                        //cout << "," << j2*vide_sample[a]+ymin_sample[a] << "->" <<QMatrix_Moy_Glob[a][i][j][i2][j2] << endl;
                        File <<  i*vide_sampleX[a]+xmin_sample[a] << " " << j*vide_sampleY[a]+ymin_sample[a] << " " << i2*vide_sampleX[a]+xmin_sample[a];
                        File << " " << j2*vide_sampleY[a]+ymin_sample[a] << " " <<QMatrix_Moy_Glob[a][i][j][i2][j2] << endl;
                    }
                }
    /*
     //ecriture sous forme "matricielle" x1 = coord en x du gene 1; y2 = coord en y du gene2
     //             x1-x2               x1-x3     ...    x2-x1.    ..
     //y1-y2 Q((x1,y1),(x2,y2))
     //y1-y3
     //...
     cout << "\n\nMatrice des proba d'identitÈ"<< endl;
     cout << setw(12) << "j\\i     ";
     File << setw(12) << "j\\i     ";
     for(int i=0;i<dim_sample1[a];i++)
     for(int i2=0;i2<dim_sample1[a];i2++) {
     cout << setw(4) << i*vide_sampleX[a]+xmin_sample[a] << "-" << setw(4) << i2*vide_sampleX[a]+xmin_sample[a] << " ";
     File << setw(4) << i*vide_sampleX[a]+xmin_sample[a] << "-" << setw(4) << i2*vide_sampleX[a]+xmin_sample[a] << " ";
     }

     cout << endl;
     File << endl;

     for(int j=0;j<dim_sample2[a];j++)
     for(int j2=0;j2<dim_sample2[a];j2++) {
     cout << setw(4) << j*vide_sampleY[a]+ymin_sample[a] << "-" << setw(4) << j2*vide_sampleY[a]+ymin_sample[a] << " ";
     File << setw(4) << j*vide_sampleY[a]+ymin_sample[a] << "-" << setw(4) << j2*vide_sampleY[a]+ymin_sample[a] << " ";
     for(int i=0;i<dim_sample1[a];i++)
     for(int i2=0;i2<dim_sample1[a];i2++) {
     cout << setw(9) << QMatrix_Moy_Glob[a][i][j][i2][j2] << " " ;
     File << setw(9) << QMatrix_Moy_Glob[a][i][j][i2][j2] << " " ;
     }
     cout << endl;
     File << endl;
     }
     */

//    cout << endl;
    File.close();
}//e dn a<sampleNb
}
/*****************************************************************************/

/*****************************************************************************/
void calcul_VariousStats() {
//modifiÈe le 31072002
//prends comme base les proba d'identitÈ.
//11022003 pb? comparaison avec mathematica
//pb car Fis plus grand pour petit torre alors que ca devrait etre l'inverse (cf mathematica)
//22022003 pb car on prends que sur l"echantillon pour calculer Qpop! faudrait toute la pop
    
// RL 06 2018 : Fis était complement faux car basé sur Qr_mean_moy...ex Qpop... qui n'était calculé que sur les pas et pas sur toutes les paires d'individus.

using namespace NS_diagnostic_tables;

double storeTITV;
vector<double> Qr_mean(2,0.0);// RL 062018 a virer car ne sert a rien
int sampleNb=1;

if(predispbool) sampleNb=2;

for(int a=0;a<sampleNb;a++){

int startIndex=1,endIndex=n_genes[0];
if(a==1) {startIndex=n_genes[0]+1;endIndex=n_genes_total;}

//	IF sequence loci hace been specified
if (seqStatsOui) {
	//	loop over total loci
	for (int loc = 0; loc < n_locus; loc++)
		//	IF current locus is of type sequence
		if((mutModelIDVector[loc+1]=="ISM") || (mutModelIDVector[loc+1]=="JC69")
				|| (mutModelIDVector[loc+1]=="K80") || (mutModelIDVector[loc+1]=="F81")
				|| (mutModelIDVector[loc+1]=="HKY85")	|| (mutModelIDVector[loc+1]=="TN93")) {
			//	double loop over all the genes in the sample for pairwise comparisons
			for (int i = startIndex; i <= endIndex - 1; i++)
				for (int j = i+1; j <= endIndex; j++) {
					//	IF the two genes selected from the sample are different from each each other
					if (noeud[j][loc].sequence != noeud[i][loc].sequence) {
						numPairMismatch[a][loc] += seqPairDiff(noeud[i][loc].sequence, noeud[j][loc].sequence);
						//	IF the current locus isn't an ISM locus then calculate transition/transversion ratio
						if (mutModelIDVector[loc+1]!="ISM") {
							storeTITV = seqPairNucSubDiff(noeud[i][loc].sequence, noeud[j][loc].sequence);
//	Only for DEBUG purposes
/*
							if (!storeTITV)
								cerr << "\nDivision by zero! Found no transversions on run " << rep << " at locus " << loc+1
									 << " between sequence " << i << " and " << j << endl
									 << "Please check the Nexus file containing the individuals if it has been generated."
									 << endl;
*/
							empTITVRatio[a][loc] += storeTITV;
						}
					}
				}
			numPairMismatch[a][loc] *= (double) 2 / (n_genes[a] * (n_genes[a] - 1));
			empTITVRatio[a][loc] *= (double) 2 / (n_genes[a] * (n_genes[a] - 1));
		}
}

//initialisations dans ini_moyenne

if(HexpNeioui && a==0) freq_sample();

for(int loc=0;loc<n_locus;loc++)
	{Qr_mean[a]=0.0;// RL 062018 a virer car ne sert a rien
        if(ploidy==2 || TVpars[0].initialDens>1) {
            if(  dens_sample[a]==1 && ( fabs(Qind2[a][loc] - Qw_meanAllInd[a][loc]) > 0.000001 ) ) {
                cout << "Qind2[a][loc] != Qw_meanAllInd[a][loc], I exit" << endl;
                if(cinGetOnError) cin.get();
                exit(-1);
            }
            hetero[a][loc]=1.0-Qind2[a][loc];
            //cout << "Het[" << a << "][" << loc << "]=" << hetero[a][loc] << " Qind2=" << Qind2[a][loc] << endl;
            //cin.get();

        }
        if(HexpNeioui){
            for(int i=allele_min[loc];i<=allele_max[loc];i++)
                HexpNei[a][loc]+=freq[a][loc][i]*(1.0-freq[a][loc][i]);
            HexpNei[a][loc]*=n_genes[a]/(n_genes[a]-1);
        }

        //pas clair ce calcul de Fis, je pense qu'a l'origine c'était pour calculer un Fis en pop continue
//	CBR->RL _pas <= 0:max(dim_sample1,dim_sample2) ou bien <= 1:max(dim_sample1,dim_sample2)?
//  for(int _pas=0;_pas<=max(dim_sample1,dim_sample2);_pas++)

    for(int _pas=1;_pas<max(dim_sample1[a],dim_sample2[a]);_pas++)// proba d'identité inter "pop" ou inter individus en pop cont, a vérifier
        //RL 062018 : faux c'est la moyenne des Qr, et pas la moyenne sur toutes les paires d'individus
    	Qr_mean[a]+=Q1[a][loc][_pas];// RL 062018 a virer car ne sert a rien

	if(max(dim_sample1[a],dim_sample2[a])>1) Qr_mean[a]=Qr_mean[a]/( (double) max(dim_sample1[a],dim_sample2[a])-1.0);
    if(Fisoui && (ploidy==2 || TVpars[0].initialDens>1) && max(dim_sample1[a],dim_sample2[a])>1) {
        if((1.0 - Qb_meanAllPairs[a][loc])>petit) fis[a][loc]=(Qind2[a][loc]-Qb_meanAllPairs[a][loc])/(1.0-Qb_meanAllPairs[a][loc]); else fis[a][loc]=0.0;
    }
	if(ploidy==2 || TVpars[0].initialDens>1) hetero_moy[a]+=hetero[a][loc];
	if(HexpNeioui) {
		HexpNei_moy[a]+=HexpNei[a][loc];
		HexpNei_moy_glob_var[a]+=pow(HexpNei[a][loc],2);
	}
	if(Fisoui && (ploidy==2 || TVpars[0].initialDens>1) && max(dim_sample1[a],dim_sample2[a])>1) {
		fis_moy_numer[a]+=Qind2[a][loc]-Qb_meanAllPairs[a][loc];
		fis_moy_denom[a]+=1.0-Qb_meanAllPairs[a][loc];
	}
	Qr_mean_moy[a]+=Qr_mean[a];// RL 062018 a virer car ne sert a rien
	}
//moyennes sur locus
if(ploidy==2 || TVpars[0].initialDens>1) hetero_moy[a]=hetero_moy[a]/n_locus;
if(HexpNeioui) HexpNei_moy[a]/=n_locus;
if(Fisoui && (ploidy==2 || TVpars[0].initialDens>1) && max(dim_sample1[a],dim_sample2[a])>1) {
    if(fis_moy_denom[a]>petit) fis_moy[a]=fis_moy_numer[a]/fis_moy_denom[a]; else fis_moy[a]=0.0;
} else fis_moy[a]=0.0;
Qr_mean_moy[a]=Qr_mean_moy[a]/n_locus; // RL 062018 a virer car ne sert a rien
if(  dens_sample[a]==1 && ( fabs(Qind_moy[a] - Qw_meanAllInd_moy[a]) > 0.000001 ) ) {
    cout << "Qind_moy[a] != Qw_meanAllInd_moy[a], I exit" << endl;
    if(cinGetOnError) cin.get();
    exit(-1);
}


if(Varoui && a==0) var_allele();


//moyenne sur repet
if(ploidy==2 || TVpars[0].initialDens>1) hetero_moy_glob[a]+=hetero_moy[a];
if(HexpNeioui) HexpNei_moy_glob[a]+=HexpNei_moy[a];
if(ploidy==2 || TVpars[0].initialDens>1) fis_moy_numer_glob[a]+=fis_moy_numer[a];
if(ploidy==2 || TVpars[0].initialDens>1) fis_moy_denom_glob[a]+=fis_moy_denom[a];
Qr_mean_moy_glob[a]+=Qr_mean_moy[a];// RL 062018 a virer car ne sert a rien
Qb_meanAllPairs_moy_glob[a]+=Qb_meanAllPairs_moy[a];
Qw_meanAllInd_moy_glob[a]+=Qw_meanAllInd_moy[a];
if(Fisoui && (ploidy==2 || TVpars[0].initialDens>1)) fis_moy_glob[a]+=fis_moy[a];//le bon
if( dens_sample[a]==1 && ( fabs(Qind_moy_glob[a] - Qw_meanAllInd_moy_glob[a]) > 0.000001) ) {
    cout << "Qind_moy_glob[a] != Qw_meanAllInd_moy_glob[a], I exit" << endl;
    if(cinGetOnError) cin.get();
    exit(-1);
}

    
    
if(rep==repet) {
  if(ploidy==2 || TVpars[0].initialDens>1) hetero_moy_glob[a]/=repet;
  if(HexpNeioui) {
  	HexpNei_moy_glob[a]/=repet;
  	if(n_locus*repet>1) HexpNei_moy_glob_var[a]=sqrt(HexpNei_moy_glob_var[a]/(n_locus*repet-1)-pow(HexpNei_moy_glob[a],2) )/sqrt(double(n_locus*repet));
    else HexpNei_moy_glob_var[a]=0.0;
  }
  if(Varoui) {
  	Var_moy_glob[a]/=repet;
  	Var_moy_glob_var[a]=sqrt(Var_moy_glob_var[a]/(n_locus*repet-1)-pow(Var_moy_glob[a],2))/sqrt(double(n_locus*repet));
  	M_moy_glob[a]/=repet;
  	M_moy_glob_var[a]=sqrt(M_moy_glob_var[a]/(n_locus*repet-1)-pow(M_moy_glob[a],2))/sqrt(double(n_locus*repet));
//      cout << "\nFinal M_moy_glob[" << a << "]=" << M_moy_glob[a] << endl;
//      cout << "Final M_moy_glob_var[" << a << "]=" << M_moy_glob_var[a] << endl;
//      fflush(stdout);
  }
  Qr_mean_moy_glob[a]/=repet;
  Qb_meanAllPairs_moy_glob[a]/=repet;
  Qw_meanAllInd_moy_glob[a]/=repet;
  if(Fisoui && (ploidy==2 || TVpars[0].initialDens>1)) {
	  fis_moy_glob[a]/=repet;//le bon
	  fis_moy_numer_glob[a]/=repet;
	  fis_moy_denom_glob[a]/=repet;
  }

    printf("\n");
  if(a==0) printf("\nVarious statistics on the genetic sample (postdisp)");
      else printf("\nVarious statistics on the predisp genetic sample");
     if(ploidy==2 || TVpars[0].initialDens>1) printf("\nMean within-deme identity prob. (within-individual for cont. pop.) =%12.10e\nMean observed heterozygosity within demes (within-individual for cont. pop.) =%12.10e",Qind_moy_glob[a],hetero_moy_glob[a]);
    printf("\nMean between-deme identity prob. (within-pop for a cont. pop.)=%12.10e \n ",Qb_meanAllPairs_moy_glob[a]);
  if(HexpNeioui) printf("\nHexpNei=%12.10e (+- %12.10e)",HexpNei_moy_glob[a],1.96*HexpNei_moy_glob_var[a]);
  if(HexpNeioui && a==0) printf("\nAllele number=%12.10e",mean_allele_nmbr);
  if(Varoui)     printf("\nVar=%12.10e (+- %12.10e)",Var_moy_glob[a],1.96*Var_moy_glob_var[a]);
  if(Varoui)     printf("\nM=%12.10e (+- %12.10e)",M_moy_glob[a],1.96*M_moy_glob_var[a]);
  if(Fisoui && (ploidy==2 || TVpars[0].initialDens>1)) {
	  printf("\nMean Fst (Fis for cont. pop.)=%12.10e ",fis_moy_glob[a]);
	  //printf("\nFis moy2=Fis multirepet=%12.10e ",fis_moy_numer_glob/fis_moy_denom_glob);
  }
 }
}

    //calcul des Hexp theoriques
    if (nMarker == 1) {	// calculating Hexp only in the case of a single user specified marker
        if (model_mut.compare("SMM")==0) {
            Hexp=1.0-1.0/sqrt(1.0+2*4*_mu*TVpars[0].initialDens);
            Vexp=4*_mu*TVpars[0].initialDens/2;
        } else if (model_mut.compare("IAM")==0) {
            Hexp=4*_mu*TVpars[0].initialDens/(1.0+4*_mu*TVpars[0].initialDens);
            Vexp=numeric_limits<double>::infinity();
        } else if (model_mut.compare("KAM")==0) {
            Hexp=1.0*(1.0 - (1+4*_mu*TVpars[0].initialDens/(Kmax-Kmin-1)) / (1.0+4*_mu*TVpars[0].initialDens*(Kmax-Kmin)/(Kmax-Kmin-1)));
            Vexp=0.25*4*_mu*TVpars[0].initialDens*(Kmax-Kmin)*(Kmax-Kmin-1)*(Kmax-Kmin+1)/6.0;
        } else {
            Hexp=numeric_limits<double>::infinity();
            Vexp=numeric_limits<double>::infinity();
        }
        if (rep==repet && (TVpars[0].dimRes1*TVpars[0].dimRes2==1 && (model_mut.compare("IAM")==0 || model_mut.compare("KAM")==0 || model_mut.compare("SMM")==0) ) ) {
            if(min_allele>1) printf("\nProportion of discarded loci with less than %d alleles = %f",min_allele,mono);
            if(max_mutations<2147483647)printf("\nProportion of loci with tooManyMut in total sample = %f",tooManyMut);
            printf("\nSome theoretical expectations for a wright fisher populaton of size: %d",TVpars[0].initialDens);
            printf("\nHexp=%12.10e",Hexp);
            printf("\nVexp=%12.10e (Infinity expected under IAM)",Vexp);
        }

    }

}
/*******************************************************************/

/***********calcul de la moyenne, sigma≤, kurtosis et skewness de disp efficcace**/

void calcul_stat_disp()
{int cumx,cumy;
    using namespace NS_diagnostic_tables;
    using namespace NS_translation;

cumx=cumy=0;
moy_dispx=moy_dispy=0.0;
sig_dispx=sig_dispy=0.0;
kurto_dispx=kurto_dispy=0.0;
skew_dispx=skew_dispy=0.0;
moy_demidispx=moy_demidispy=0.0;
sig_demidispx=sig_demidispy=0.0;
kurto_demidispx=kurto_demidispy=0.0;
skew_demidispx=skew_demidispy=0.0;

for(int i=0;i<=2*grande_dim;i++)
 {cumx+=effective[i][x];
  cumy+=effective[i][y];
  effective_moy[i][x]+=effective[i][x];
  effective_moy[i][y]+=effective[i][y];
 }
//printf("\n cumx: %d; cumy: %d",cumx,cumy);
//printf("\n cum_moyx: %d; cum_moyy: %d",cum_moyx,cum_moyy);


for(int i=0;i<=2*grande_dim;i++)
{deffective[i][x]=double(effective[i][x])/ double(cumx);
 deffective[i][y]=double(effective[i][y])/ double(cumy);
 moy_dispx+=deffective[i][x]*(i-grande_dim);
 moy_dispy+=deffective[i][y]*(i-grande_dim);
 sig_dispx+=deffective[i][x]*pow(double(i-grande_dim),double(2));
 sig_dispy+=deffective[i][y]*pow(double(i-grande_dim),double(2));
 kurto_dispx+=deffective[i][x]*pow(double(i-grande_dim),double(4));
 kurto_dispy+=deffective[i][y]*pow(double(i-grande_dim),double(4));
 skew_dispx+=deffective[i][x]*pow(double(i-grande_dim),double(3));
 skew_dispy+=deffective[i][y]*pow(double(i-grande_dim),double(3));
 moy_demidispx+=deffective[i][x]*fabs( double(i-grande_dim));
 moy_demidispy+=deffective[i][y]*fabs( double(i-grande_dim));
}

for(int i=0;i<=2*grande_dim;i++)
{sig_demidispx+=deffective[i][x]*pow(fabs( double(i-grande_dim))-moy_demidispx,double(2));
 sig_demidispy+=deffective[i][y]*pow(fabs( double(i-grande_dim))-moy_demidispy,double(2));
 kurto_demidispx+=deffective[i][x]*pow(fabs( double(i-grande_dim))-moy_demidispx,double(4));
 kurto_demidispy+=deffective[i][y]*pow(fabs( double(i-grande_dim))-moy_demidispy,double(4));
 skew_demidispx+=deffective[i][x]*pow(fabs( double(i-grande_dim))-moy_demidispx,double(3));
 skew_demidispy+=deffective[i][y]*pow(fabs( double(i-grande_dim))-moy_demidispy,double(3));
}

if(sig_dispx!=0.0)
  {kurto_dispx=kurto_dispx/pow(sig_dispx,double(2))-3.0;
   skew_dispx=skew_dispx/pow(sig_dispx,double(3/2));
   kurto_demidispx=kurto_demidispx/pow(sig_demidispx,double(2))-3.0;
   skew_demidispx=skew_demidispx/pow(sig_demidispx,double(3/2));
  }
 else
   {kurto_dispx=0.0;
    skew_dispx=0.0;
    kurto_demidispx=0.0;
    skew_demidispx=0.0;
   }
if(sig_dispy!=0.0)
  {kurto_dispy=kurto_dispy/pow(sig_dispy,double(2))-3.0;
   skew_dispy=skew_dispy/pow(sig_dispy,double(3/2));
   kurto_demidispy=kurto_demidispy/pow(sig_demidispy,double(2))-3.0;
   skew_demidispy=skew_demidispy/pow(sig_demidispy,double(3/2));
  }
 else
   {kurto_dispy=0.0;
    skew_dispy=0.0;
    kurto_demidispy=0.0;
    skew_demidispy=0.0;
   }

}

void calcul_stat_disp_moy()
{
    using namespace NS_diagnostic_tables;
    using namespace NS_translation;

moy_dispx=moy_dispy=0.0;
sig_dispx=sig_dispy=0.0;
kurto_dispx=kurto_dispy=0.0;
skew_dispx=skew_dispy=0.0;
moy_demidispx=moy_demidispy=0.0;
sig_demidispx=sig_demidispy=0.0;
kurto_demidispx=kurto_demidispy=0.0;
skew_demidispx=skew_demidispy=0.0;

for(int i=0;i<=2*grande_dim;i++)
 {cum_moyx+=effective_moy[i][x];
  cum_moyy+=effective_moy[i][y];
  deffective[i][x]=0.0;
  deffective[i][y]=0.0;

 }
for(int i=0;i<=2*grande_dim;i++)
{deffective[i][x]=((double) effective_moy[i][x])/ double(cum_moyx);
 deffective[i][y]=((double) effective_moy[i][y])/ double(cum_moyy);
 moy_dispx+=deffective[i][x]*(i-grande_dim);
 moy_dispy+=deffective[i][y]*(i-grande_dim);
 sig_dispx+=deffective[i][x]*pow(double(i-grande_dim),double(2));
 sig_dispy+=deffective[i][y]*pow(double(i-grande_dim),double(2));
 kurto_dispx+=deffective[i][x]*pow(double(i-grande_dim),double(4));
 kurto_dispy+=deffective[i][y]*pow(double(i-grande_dim),double(4));
 skew_dispx+=deffective[i][x]*pow(double(i-grande_dim),double(3));
 skew_dispy+=deffective[i][y]*pow(double(i-grande_dim),double(3));
 moy_demidispx+=deffective[i][x]*fabs( double(i-grande_dim));
 moy_demidispy+=deffective[i][y]*fabs( double(i-grande_dim));
}

for(int i=0;i<=2*grande_dim;i++)
{sig_demidispx+=deffective[i][x]*pow(fabs( (double) (i-grande_dim))-moy_demidispx,double(2));
 sig_demidispy+=deffective[i][y]*pow(fabs( (double) (i-grande_dim))-moy_demidispy,double(2));
 kurto_demidispx+=deffective[i][x]*pow(fabs( (double) (i-grande_dim))-moy_demidispx,double(4));
 kurto_demidispy+=deffective[i][y]*pow(fabs( (double) (i-grande_dim))-moy_demidispy,double(4));
 skew_demidispx+=deffective[i][x]*pow(fabs( (double) (i-grande_dim))-moy_demidispx,double(3));
 skew_demidispy+=deffective[i][y]*pow(fabs( (double) (i-grande_dim))-moy_demidispy,double(3));
}

if(sig_dispx!=0.0)
  {kurto_dispx=kurto_dispx/pow(sig_dispx,double(2))-3.0;
   skew_dispx=skew_dispx/pow(sig_dispx,double(3/2));
   kurto_demidispx=kurto_demidispx/pow(sig_demidispx,double(2))-3.0;
   skew_demidispx=skew_demidispx/pow(sig_demidispx,double(3/2));
  }
 else
   {kurto_dispx=0.0;
    skew_dispx=0.0;
    kurto_demidispx=0.0;
    skew_demidispx=0.0;
   }
if(sig_dispy!=0.0)
  {kurto_dispy=kurto_dispy/pow(sig_dispy,double(2))-3.0;
   skew_dispy=skew_dispy/pow(sig_dispy,double(3/2));
   kurto_demidispy=kurto_demidispy/pow(sig_demidispy,double(2))-3.0;
   skew_demidispy=skew_demidispy/pow(sig_demidispy,double(3/2));
  }
 else
   {kurto_dispy=0.0;
    skew_dispy=0.0;
    kurto_demidispy=0.0;
    skew_demidispy=0.0;
   }
}
/******************************************************************************/

/******* Calculates the allele frequency spectrum of the current locus ******/
bool checkMinorAlleleFreq()
{
	//	allocating, initializing and normalzing the allele frequency spectrum
	double *alleleFreq = dvector(allelesAtCurrentLocus + 1);
	for (int i = 1; i <= allelesAtCurrentLocus; i++)
			*(alleleFreq + i) = 0.0;

	for (int i = 1; i <= allelesAtCurrentLocus; i++)
		for (int k = 1; k <= n_genes_total; k++)
			if (alleleID[i] == noeud[k][locus].etat_allele)
				++(*(alleleFreq + i));

	for (int i = 1; i <= allelesAtCurrentLocus; i++)
			*(alleleFreq + i) /= n_genes_total;

	//	CBR: N.B. sort() ignores the last index of the specified range delimiters!
	sort(alleleFreq + 1, alleleFreq + allelesAtCurrentLocus + 1);

	//	check if the second largest allele frequency meets the specified limit
	if (*(alleleFreq + allelesAtCurrentLocus - 1) >= minorAlleleFreqVector[locus+1]) {
	    free_dvector(alleleFreq);
		return true;
	}
	else {
	    free_dvector(alleleFreq);
		return false;
	}
}

/*********calcul des frequences alleliques dans l'echantillon***/
void freq_sample() {
int sampleNb=1;

if(predispbool) sampleNb=2;

for(int a=0;a<sampleNb;a++){
//	freq is initialized here as maximum_allele_number is determined only at the end of the locus loop
	for (int i=0;i<n_locus;i++) {
		for (int j = 1; j <= maximum_allele_number; j++) {
			freq[a][i][j] = 0.0;
		}
	}

int startIndex=1,endIndex=n_genes[0];
if(a==1) {startIndex=n_genes[0]+1;endIndex=n_genes_total;}

	for(int j=0;j<n_locus;j++){
		for(int i=startIndex;i<=endIndex;i++) {
			for(int k=allele_min[j];k<=allele_max[j];k++) {
				if(noeud[i][j].etat_allele==k){
					freq[a][j][k]+=1.0;
				}
			}
		}
		for(int k=allele_min[j];k<=allele_max[j];k++)
            freq[a][j][k]/=n_genes[a];
	}
} //end for a<sampleNb
//return(0);
}

/****************calcul des SSb et SSw pour ensuite calculer ar/er/etc, cf Rousset 2000JEB ************/
void computeSSForFstats() {
    using namespace NS_diagnostic_tables;
size_t sampleNb=1,indVar;
vector<vector<vector<vector<double> > > > indicVar;
    
//computation of indicator variables for all loci, all individuals, all genes, and all allelic states
indicVar.resize(n_locus);//locus indexing starting at 0
//loop over each locus
for(int l=0;l<n_locus;l++){
    indicVar[l].resize(n_genes_total-ploidy+2);//To keep node/gene indexing starting at 1, and use 0 to store the mean
    //loop over each individuals
    for(int i=1;i<=(n_genes_total-ploidy+1);i+=ploidy) {
        indicVar[l][i].resize(ploidy+1);//gene indexing starting at 1, and use 0 to store the mean
        indicVar[l][i][0].resize(ploidy+1);//gene indexing starting at 1, and use 0 to store the mean
        for(int j=1;j<=ploidy;j++) {
            indicVar[l][i][0].resize(maximum_allele_number+1);
        }
        //loop over each gene of each individual
        for(int j=1;j<=ploidy;j++) {
            indicVar[l][i][j].resize(maximum_allele_number+1);
            //Loop over each possible allelic type
            for(int u=1;u<=maximum_allele_number;u++) {//allele indexing starting at 1
                if(noeud[i+j-1][l].etat_allele==u) indVar=1; else indVar=0;
                indicVar[l][i][j][u]=indVar;
                indicVar[l][i][0][u]+=1.0*indVar/ploidy;
                //indicVar[l][0][0][u]+=indVar/ploidy/(endGene-startGene+1); //RL 062018 not used to be removed
            }
        }
    }
//    for(int i=1;i<=(n_genes_total-ploidy+1);i+=ploidy) {
//        for(int j=1;j<=ploidy;j++) {
//            for(int u=1;u<=maximum_allele_number;u++) {//allele indexing starting at 1
//                cout << /*"X_" << i << "_" << j  << "_" << u << "=" <<*/ indicVar[l][i][j][u] << " ";
//            }
//        }
//        cout << endl;
//    }
//    cout << endl;
//    
//    for(int i=1;i<=(n_genes_total-ploidy+1);i+=ploidy) {
//        for(int u=1;u<=maximum_allele_number;u++) {//allele indexing starting at 1
//            cout << /*"X_" << i << "."  << "_" << u << "=" <<*/ indicVar[l][i][0][u] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;
}

    
if(predispbool) sampleNb=2;
    SSw.resize(sampleNb);SSb.resize(sampleNb);SSw_SumOverAllPairs.resize(sampleNb);
    SSb_SumForEachIndOverOthers.resize(sampleNb);

for(int a=0;a<sampleNb;a++){//lopp over pre and postdisp samples
    
    size_t startGene=1,endGene=n_genes[0],pairNb=n_genes[0]/ploidy*(n_genes[0]/ploidy-1)/2;
    if(a==1) {startGene=n_genes[0]+1;endGene=n_genes_total;pairNb=n_genes[1]/ploidy*(n_genes[1]/ploidy-1)/2;}

//    for(size_t l=0;l<n_locus;l++){
//        size_t pair=0;
//        for(size_t i=startGene+ploidy;i<=(endGene-ploidy+1);i+=ploidy) {
//            for(size_t ii=startGene;ii<i;ii+=ploidy) {
//                for(size_t u=1;u<=maximum_allele_number;u++) {//allele indexing starting at 1
//                    cout << /*"pair(" << i << "," << ii << ") X.."  << "_" << u << "=" << */(indicVar[l][i][0][u] + indicVar[l][ii][0][u])/2.0 << " ";
//                }
//                pair++;
//                cout << "  ";
//            }
//            cout << endl;
//        }
//    }
    //computation of SSw,SSb for all samples, all individual pairs, and all loci
    SSw[a].resize(pairNb);SSb[a].resize(pairNb);SSw_SumOverAllPairs[a].resize(n_locus,0.0);
    SSb_SumForEachIndOverOthers[a].resize((endGene-startGene+1)/ploidy);

    size_t pair=0;
    //loop over each individuals
    SSb_SumForEachIndOverOthers[a][(startGene-1)/ploidy].resize(n_locus,0.0);
    for(size_t i=startGene+ploidy;i<=(endGene-ploidy+1);i+=ploidy) {
        SSb_SumForEachIndOverOthers[a][(i-1)/ploidy].resize(n_locus,0.0);
        //loop over "all<" other individuals consider all individual pair
        for(size_t ii=startGene;ii<i;ii+=ploidy) {
            SSw[a][pair].resize(n_locus,0.0);SSb[a][pair].resize(n_locus,0.0);
            //loop over each locus
           for(size_t l=0;l<n_locus;l++){
               //loop over each gene of each individual
               for(size_t j=1;j<=ploidy;j++) {
                   //Loop over each possible allelic type
                   for(size_t u=1;u<=maximum_allele_number;u++) {
                       //cout << "pair i=" << i << ", ii=" << ii << " j=" << j  << " u=" << u << "; Xi_j_u=" << indicVar[l][i][j][u] << endl;
                       SSw[a][pair][l]+=sq( (double) indicVar[l][i][j][u] - indicVar[l][i][0][u] );
                       SSw[a][pair][l]+=sq( (double) indicVar[l][ii][j][u] - indicVar[l][ii][0][u] );
                       SSb[a][pair][l]+=sq( (double) indicVar[l][i][0][u] - (indicVar[l][i][0][u] + indicVar[l][ii][0][u])/2.0 );
                       SSb[a][pair][l]+=sq( (double) indicVar[l][ii][0][u] - (indicVar[l][i][0][u] + indicVar[l][ii][0][u])/2.0 );
                   }
               }
               SSw_SumOverAllPairs[a][l]+=SSw[a][pair][l];
               SSb_SumForEachIndOverOthers[a][(ii-1)/ploidy][l]+=SSb[a][pair][l];
               SSb_SumForEachIndOverOthers[a][(i-1)/ploidy][l]+=SSb[a][pair][l];
//               cout << "Loc=" << l << "; pair" << pair << " i=" << i << ", ii=" << ii << " SSb=" << SSb[a][pair][l] << " SSw=" << SSw[a][pair][l] << endl;
//               if(pair == (pairNb-1)) cout <<"Loc=" << l <<  "; Last pair, SSw_SumOverAllPairs=" << SSw_SumOverAllPairs[a][l] << endl;
           }
           pair++;
        }
    }

    

}

    
}

/****************calcul des distances entre paires d'individus, cf Rousset 2000JEB ************/
void computeGeoDist() {
    using namespace NS_diagnostic_tables;
    using namespace NS_coal_tree;
    size_t sampleNb=1;
    
    if(predispbool) sampleNb=2;
    indGeoDist.resize(sampleNb);popGeoDist.resize(sampleNb);maxDistSample.resize(sampleNb);
    
    for(int a=0;a<sampleNb;a++){//lopp over pre and postdisp samples
        
        size_t startGene=1,endGene=n_genes[0],indPairNb=n_genes[0]/ploidy*(n_genes[0]/ploidy-1)/2,popPairNb;
        if(a==1) {startGene=n_genes[0]+1;endGene=n_genes_total;indPairNb=n_genes[1]/ploidy*(n_genes[1]/ploidy-1)/2;}
        
        //computation of geographic distances between all sampled individual pairs
        indGeoDist[a].resize(indPairNb,0.0);
        
        size_t pair=0;
        //loop over each individuals
        for(size_t i=startGene+ploidy;i<=(endGene-ploidy+1);i+=ploidy) {
            //loop over "all<" other individuals consider all individual pair
            for(size_t ii=startGene;ii<i;ii+=ploidy) {
                indGeoDist[a][pair]=sqrt( sq( (double) (coord_individus[i][x] - coord_individus[ii][x]) ) + sq( (double) (coord_individus[i][y] - coord_individus[ii][y]) ) );
                pair++;
            }
        }
        maxDistSample[a]=*max_element(indGeoDist[a].begin(), indGeoDist[a].end());
        //RL 062018 a virer
        if(pair != indPairNb) cerr << "!! pair=" << pair << " is different from indPairNb=" << indPairNb << endl;

        //computation of geographic distances between all sampled population pairs, chercher pop partout, deja fait quelque part
        if(!Specific_Sample_Designbool[a]) popPairNb=dim_sample1[a]*dim_sample2[a]*(dim_sample1[a]*dim_sample2[a]-1)/2;
        else popPairNb=Spec_SampleSize[a]*(Spec_SampleSize[a]-1)/2;
        popGeoDist[a].resize(popPairNb,0.0);
        
        //set<vector<int> > popPairSet;
        pair=0;
        //loop over each individuals
        for(size_t i=startGene+ploidy;i<=(endGene-ploidy+1);i+=ploidy) {
            //loop over "all<" other individuals consider all individual pair
            for(size_t ii=startGene;ii<i;ii+=ploidy) {
                //if( (coord_individus[i][x] != coord_individus[j][x]) || (coord_individus[i][y] != coord_individus[j][y]) ) {
                    //vector<int> pairVec(4);
                    //pairVec.push_back(coord_individus[i][x]);pairVec.push_back(coord_individus[j][x]);pairVec.push_back(coord_individus[i][y]);pairVec.push_back(coord_individus[j][y]);
                    //if(popPairSet.empty() || ( popPairSet.find(pairVec) != popPairSet.end() ) ) {
                        indGeoDist[a][pair]=sqrt( sq( (double) (coord_individus[i][x] - coord_individus[ii][x]) ) + sq( (double) (coord_individus[i][y] - coord_individus[ii][y]) ) );
                        //copy( pairVec.begin(), pairVec.end(), inserter( popPairSet, popPairSet.end() ) );
                        pair++;
                    //}
                //}
            }
        }
        //RL 062018 a virer
        //if(pair != popPairNb) cerr << "!! pair=" << pair << " is different from popPairNb=" << popPairNb << endl;
        
    }
}

/****************calcul de ar et er, cf Rousset 2000JEB ************/
void computeFstats_ar_er_moranI(bool commonSSwInNumAndDenom=true) {
    using namespace NS_diagnostic_tables;
    size_t sampleNb=1;
    double sumLocDenom=0.0;
    vector <double> sumLocNum_ar,sumLocNum_er;
    
    if(predispbool) sampleNb=2;
    Fstat_ar_fromSS.resize(sampleNb);Fstat_ar_fromQ.resize(sampleNb);
    Fstat_er_fromQ.resize(sampleNb);Fstat_moranI_fromQ_forEachIndPairs.resize(sampleNb);
    Fstat_moranI_fromQ_forEachDistClass.resize(sampleNb);
    
    for(size_t a=0;a<sampleNb;a++){//lopp over pre and postdisp samples
        
        size_t startGene=1,endGene=n_genes[0],pairNb=n_genes[0]/ploidy*(n_genes[0]/ploidy-1)/2;
        if(a==1) {startGene=n_genes[0]+1;endGene=n_genes_total;pairNb=n_genes[1]/ploidy*(n_genes[1]/ploidy-1)/2;}
        //indNb=(endGene-startGene+1)/ploidy;
        sumLocNum_ar.resize(pairNb,0.0);
        sumLocNum_er.resize(pairNb,0.0);

        //loop over each locus (multilocus estimator is the ratio of the sum of locus specific num & denom)
        for(size_t l=0;l<n_locus;l++){
            //loop over each pair of individuals
            size_t pair=0;
            //loop over each individuals
            for(size_t i=startGene+ploidy;i<=(endGene-ploidy+1);i+=ploidy) {
                //loop over "all<" other individuals consider all individual pair
                for(size_t ii=startGene;ii<i;ii+=ploidy) {
                    if (! commonSSwInNumAndDenom ) {
                        sumLocNum_ar[pair] += (2*SSb[a][pair][l] - SSw[a][pair][l])*pairNb;

                    } else {
                        sumLocNum_ar[pair] += SSb[a][pair][l]*pairNb;
                    }
                    pair++;
                }
            }
            sumLocDenom+=SSw_SumOverAllPairs[a][l];
        }
        
        size_t pair=0;
       
        Fstat_ar_fromSS[a].resize(pairNb);
        Fstat_ar_fromQ[a].resize(pairNb);
        Fstat_er_fromQ[a].resize(pairNb);
        Fstat_moranI_fromQ_forEachIndPairs[a].resize(pairNb);
        Fstat_moranI_fromQ_forEachDistClass[a].resize(maxDistSample[a]);
        //loop over each pair
        pair=0;
        //loop over each individuals
        for(size_t i=startGene+ploidy;i<=(endGene-ploidy+1);i+=ploidy) {
            //loop over "all<" other individuals consider all individual pair
            for(size_t ii=startGene;ii<i;ii+=ploidy) {
                if (! commonSSwInNumAndDenom ) Fstat_ar_fromSS[a][pair] = sumLocNum_ar[pair]/(2*sumLocDenom);
                 else Fstat_ar_fromSS[a][pair] = sumLocNum_ar[pair]/sumLocDenom - 1.0/2.0;
                
                Fstat_ar_fromQ[a][pair]=(Qw_pair_moy[a][pair] - Qb_pair_moy[a][pair])/(1 - Qw_meanAllInd_moy[a]);
//                cout << endl << "ar from SS=" << Fstat_ar_fromSS[a][pair] << "=SumLocNumPair/sumLocDenom - 0.5=" << sumLocNum_ar[pair] << "/" << sumLocDenom << "-0.5";
//                cout << " (2*sumLocNum-sumLocDenom=" << 2*sumLocNum_ar[pair] - sumLocDenom << ")" << endl;
//                cout << "ar from Q=" << Fstat_ar_fromQ[a][pair] <<"=(Qw_pair_moy - Qb_pair_moy) / (1 - Qw_meanAllInd_moy[a])="
//                        << (Qw_pair_moy[a][pair] - Qb_pair_moy[a][pair]) << "/" << (1 - Qw_meanAllInd_moy[a]) << endl;
                if ((!commonSSwInNumAndDenom && ( fabs( Fstat_ar_fromQ[a][pair] - Fstat_ar_fromSS[a][pair] ) > 0.000001 ) ) ) {
                    cout << endl << "ar from SS=" << Fstat_ar_fromSS[a][pair] << "=SumLocNumPair/(2*sumLocDenom) =" << sumLocNum_ar[pair] << "/" << 2*sumLocDenom << endl;
                    cout << "ar from Q=" << Fstat_ar_fromQ[a][pair] <<"=(Qw_pair_moy - Qb_pair_moy) / (1 - Qw_meanAllInd_moy[a])="
                            << (Qw_pair_moy[a][pair] - Qb_pair_moy[a][pair]) << "/" << (1 - Qw_meanAllInd_moy[a]) << endl;
                
                }
                
                Fstat_er_fromQ[a][pair] = ( ( QiOverAllOtherInd_moy[a][ ( i-1 ) / ploidy ] + QiOverAllOtherInd_moy[a][ ( ii -1 ) / ploidy ] ) - Qb_pair_moy[a][pair] /* - Qw_meanAllInd_moy[a]*/ ) / ( 1 - Qw_meanAllInd_moy[a] );
                //cout << endl << endl << "er from Q=( Qb_pair - ( QiOverAllOtherInd + QiOverAllOtherInd ) ) / ( 1 - Qw_meanAllInd) =" << Fstat_er_fromQ[a][pair] << endl;
                //cout << "Qb_pair_moy=" << Qb_pair_moy[a][pair] << "; QiOverAllOtherInd_moy=" << QiOverAllOtherInd_moy[a][ ( i-1 ) / ploidy ] << "; QjOverAllOtherInd_moy=" << QiOverAllOtherInd_moy[a][ ( ii -1 ) / ploidy ] << "; (1 - Qw_meanAllInd) = " <<  ( 1 - Qw_meanAllInd_moy[a]) ;

                Fstat_moranI_fromQ_forEachIndPairs[a][pair] = ( Qb_pair_moy[a][pair] - Qb_meanAllPairs_moy[a] ) / ( (1 + Qw_meanAllInd_moy[a] ) / 2 - Qb_meanAllPairs_moy[a] );
                //cout << "moranI_indPairs from Q=" << Fstat_moranI_fromQ_forEachDistClass[a][pair] << endl;

                pair++;
            }
        }
        for(size_t d=0;d<=(int) floor(maxDistSample[a]);d++) {
            Fstat_moranI_fromQ_forEachDistClass[a][pair] = ( Qb_distClass_moy[a][d] - Qb_meanAllPairs_moy[a] ) / ( (1 + Qw_meanAllInd_moy[a] ) / 2 - Qb_meanAllPairs_moy[a] );
            //cout << "moranI_DistClass from Q=" << Fstat_moranI_fromQ_forEachDistClass[a][d] << endl;
        }
        
    }
}

/****************regression sur donnés individuelles de ar, er, et la distance/log cf Rousset 2000JEB ************/
vector<vector<double> > indRegression_ar_er_moranI(double minDist, const vector<vector<double> >& genetDist,string Fstat_str) {
using namespace NS_diagnostic_tables;
size_t sampleNb=1;
vector <double> locGeoDist,locGenetDist;
vector<vector<double> > regRes;

if(predispbool) sampleNb=2;
regRes.resize(sampleNb);

if(minDist<=0 && TVpars[0].dimRes2>1 ) {
    cerr << endl << "log(distance) for 2D IBD regression can not be computed because MinDistReg<=0 " <<  endl;
    cerr << "MinDistReg is set to 0.000001. Change settings if you're not satisified with this setting." <<  endl;
    minDist=0.000001;
    minDistReg=0.000001;
//    if (cinGetOnError) cin.get();
//    exit(-1);

}
    
for(int a=0;a<sampleNb;a++){//lopp over pre and postdisp samples
//    for(auto& e : indGeoDist[a]) cout << e << " ";
//    cout << endl;
//    for(auto& e : genetDist[a]) cout << e << " ";
//    cout << endl;

    for (int i=0;i<indGeoDist[a].size();i++){
        if(indGeoDist[a][i]>minDist) {
            if(TVpars[0].dimRes2>1) locGeoDist.push_back( log( indGeoDist[a][i] ) );
            else locGeoDist.push_back( indGeoDist[a][i] );
                locGenetDist.push_back( genetDist[a][i] );
        }
    }
    
    regRes[a].resize(2);
    
    //ecriture fichier .GRA with coordinnates of points used for the regression
    if(graFilebool) {
        string graFileName;
        ofstream fgra;
        stringstream stst,stst2;
        if(a>0) stst2 << "_predisp_" << Fstat_str; else stst2 << "_" << Fstat_str;;
        
        if (repet!=1) {
            stst << fichier_genepop << rep << stst2.str() << ".GRA";
        } else stst << fichier_genepop << stst2.str() << ".GRA";
        graFileName=stst.str();
        
        fgra.open(graFileName.c_str(), std::ios::out);
        
        if (!fgra.is_open()) {
            cerr << "Could not open the file " << graFileName << "for writing."  << endl;
            cerr << "\nAborting IBDSim...Press any key to exit.\n";
            if(cinGetOnError) cin.get();
            exit(-1);
        }
        
        for(int i=0;i<indGeoDist[a].size();i++) fgra << indGeoDist[a][i] << " " << genetDist[a][i] << endl;
        
        fgra.close();
    }
    regRes[a]=getLinearFit(locGeoDist,locGenetDist);
}
return regRes;
}


//*********calcul de la variance de taille allelique dans l'echantillon***/
void var_allele() {
 double *Mean,*Mean2;

int sampleNb=1;

if(predispbool) sampleNb=2;

for(int a=0;a<sampleNb;a++){

    Mean=dvector(n_locus);
    Mean2=dvector(n_locus);
    int startIndex=1,endIndex=n_genes[0];
    if(a==1) {startIndex=n_genes[0]+1;endIndex=n_genes_total;}

    //calcul la moyenne et la moyenne des carrÈs sans frequences
    for(int j=0;j<n_locus;j++){
        Mean[j]=0.0;
        Mean2[j]=0.0;

        for(int i=startIndex;i<=endIndex;i++){
            Mean[j]+=noeud[i][j].etat_allele;
            Mean2[j]+=pow(double(noeud[i][j].etat_allele),double(2));
        }
        Mean[j]/=n_genes[a];
        Mean2[j]/=n_genes[a];
        //printf("Mean[%d]=%12.10e",j,Mean[j]);
        //printf("Mean2[%d]=%12.10e",j,Mean2[j]);
        //getchar();
    }
    for(int j=0;j<n_locus;j++){
        //for(i=1;i<=n_genes;i++) Var[j]+=pow((noeud[i][j].etat_allele- Mean[j]),2);
        //Var[j]/=n_genes;
        Var[a][j]=(Mean2[j]-pow(Mean[j],double(2)))*n_genes[a]/(n_genes[a]-1);
        Var_moy[a]+=Var[a][j];
        Var_moy_glob_var[a]+=pow(Var[a][j],double(2));
    }
    Var_moy[a]/=n_locus;
    Var_moy_glob[a]+=Var_moy[a];

    free_dvector(Mean);
    free_dvector(Mean2);

    for(int j=0;j<n_locus;j++){
        M[a][j]= (double) allele_nmbr[j]/(allele_max[j]-allele_min[j]+1.0);
        //printf("M[%d]=%e, nre_allele=%d, max=%d, min=%d",j,M[j],allele_nmbr[j],allele_max[j],allele_min[j]);
        //getchar();
        M_moy[a]+=M[a][j];
        M_moy_glob_var[a]+=pow(M[a][j],double(2));
    }
    M_moy[a]/=n_locus;
    M_moy_glob[a]+=M_moy[a];
//    cout << "M_Moy[" << a << "]=" << M_moy[a] << endl;
//    cout << "M_moy_glob[" << a << "]=" << M_moy_glob[a] << endl;

}//end for a<sampleNb
}


/*-------------------------------------------------------------------------*/
/*------------------loi binomiale pour mutations----------------------------*/
double binom(double pp,int n)
{int j;
 int nold=(-1);
 double am,em,g,angle,p,bnl,sq,t,z,pold=(-1),pc,plog,pclog,en,oldg;

p=(pp<=0.5 ? pp : 1.0-pp);
pc=1.0-p;
plog=log(p);
pclog=log(pc);
am=n*p;
if(n<1000)
 {bnl=0.0;
  for(j=1;j<=n;j++) if(alea()<p) ++bnl;
 }
 else
  if(am<1.0)
	{g=exp(-am);
	 t=1.0;
	 for(j=0;j<=n;j++) {t*=alea(); if(t<g) break;}
	 bnl=(j<=n ? j : n);
	}
	else
	 {if(n!=nold)
		{en=n;
		oldg=gammln(en+1.0);
		//nold=n;
		}
	  if(p!=pold)
		{pc=1.0-p;
		 plog=log(p);
		 pclog=log(pc);
		 //pold=p;
		}
	  sq=sqrt(2.0*am*pc);
	  do
		{do
		 {angle=PI*alea();
		  z=tan(angle);
		  em=sq*z+am;
		 } while(em<0.0||em>=(en+1.0));
		 em=floor(em);
		 t=1.2*sq*(1.0+z*z)*exp(oldg-gammln(em+1.0)-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while(alea()>t);
		bnl=em;
//printf("\n binomiale=%f",bnl);
//cout<<pp<<" "<<n;
//getchar();
	 }
  if(p!=pp) bnl=n-bnl;
  //printf("\n binomiale=%f",bnl);
  //getchar();
  return bnl;
}
/*---------------------------------------------------------------------------*/

/*--------------------loi gamma pour loi de poisson--------------------------*/
double gammln(double xx)
{double xxx,tmp,ser;
 static double cof[6]={76.18009173,-86.50532033,24.01409822,
	-1.231739516,0.120858003e-2,-0.53363382e-5};
 int j;

 xxx=xx-1.0;
 tmp=xxx+5.5;
 tmp-=(xxx+0.5)*log(tmp);
 ser=1.0;
 for(j=0;j<5;j++)
	{xxx+=1.0;
	 ser+=cof[j]/xxx;
	}
 return -tmp+log(2.50662827465*ser);
}
/*---------------------------------------------------------------------------*/

/*-----------loi geometrique pour TPM et GSM---------------------------------*/
int loi_geom(const double var_geom)
{int g;
 double r,a,c,aleat;

do {
	r=((1.0+sqrt(1.0+4.0*var_geom))/2);
	c=(1.0/(r-1.0));
	a=(1.0/(1.0+c));
	do aleat=alea(); while (aleat<=0.0);
	g=int(floor((log(aleat/c))/log(a)));
	}while(g<=0);
return g;
}
/*---------------------------------------------------------------------------*/
/*----------essayer une autre procedure pour loi geometrique----------------*/


/*------------------procedure mutation---------------------------------------*/
int mutation(long int temps)
{ //double m;

if (temps==0) n_mut=0;
else
 {//m=temps*mu;
 //n_mut=poisson(m);
 //n_mut=arnaud(m);
 n_mut=int(binom(mu,temps));
 }


//if(temps!=0) {mut_moy_glob+= double(n_mut)/temps;mut_moy_glob_var+=pow(double(n_mut)/temps,double(2));}
//compte_mut_glob+=n_mut;//RL 022017 not used
compte_mut_loc+=n_mut;
return n_mut;
}
/*---------------------------------------------------------------------------*/


/*------------------loi de poisson pour procedure mutation-------------------*/
double poisson(double poiss)
{static double sq,alpoiss,g,oldm=(-1.0);
 double em,t,y1;
 float aleat;

 if(poiss<12.0)
 {if(poiss!=oldm)
	{oldm=poiss;
	 g=exp(-poiss);
	}
  em=-1;
  t=1.0;
  do{aleat=alea();
	  em+=1.0;
	  t*=aleat;
	  }while(t>g);

 }
  else
	{if(poiss!=oldm)
		{oldm=poiss;
		 sq=sqrt(2.0*poiss);
		 alpoiss=log(poiss);
		 g=poiss*alpoiss-gammln(poiss+1.0);
		}
	 do{
			do {aleat=alea();
				 y1=tan(PI*aleat);
				 em=sq*y1+poiss;
				}while(em<0.0);
			em=floor(em);
			t=0.9*(1.0+y1*y1)*exp(em*alpoiss-gammln(em+1.0)-g);
		}while(aleat>t);
	 }
	return em;
}
