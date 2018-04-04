/***************************************************************************
© R. leblois 2005- for code collection
leblois@mnhn.fr

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
#ifndef SANS_PTR_H
#define SANS_PTR_H

#include <string>
#include <vector>

#include "genepop_extract.h"

#define PI 3.14159265358979323846264338328

// comments for each function are in the main .cpp files before each function implementation
//order is the same than in the main .cpp file
int analyse(int mls,int tranche);
int analyse(int mls);
int sortie_sup_circ_abs1(int coord, int dim);
int sortie_inf_circ_abs1(int coord, int dim);
int sortie_sup_circ_abs2(int coord, int dim, int trans);
int sortie_inf_circ_abs2(int coord, int dim, int trans);
int sortie_sup_refl1(int coord, int dim);
int sortie_inf_refl1(int coord, int dim);
int sortie_sup_refl2(int coord, int dim, int trans);
int sortie_inf_refl2(int coord, int dim, int trans);
int AllocNotFixDisp(const int dxmax,const int dymax);
void affichage_ecran(void);
void ini_moy_glob(void);
void ini_tableaux(void);/*initialisation des tableaux par locus*/
void ini_tableaux_specifique(void);//initialisation pour echantillon special
void ini_struct_noeud(void);
void ini_struct_noeud2(void);/*reinitialisation si locus
									  precedent monomorphe ou a trop de mutations*/
void ini_moy(void);
void set_new_disp(int locPhase);
int getDensityFromFile(void);
int  calcul_densite(int coordx, int coordy);
//int cumul1();
//void cumul2(void);
void dispxy(void);
int  decim_tabdis(int coordx);
void disp_fixe(void);
void cumul_fixex(int k);
void cumul_fixey(int k);
int  dispx_fixe(void);
int dispy_fixe(void);
int ComputeDemeSize(void);
void nbre_aleat_noeud(void);
double gamma_mut(double rate);/*pour taux de mut variable*/
void allele_count (void);
int K_ini(void);
void mutate_under_KAM(void);
void mutate_under_SMM(void);
void mutate_under_IAM(void);
void mutate_under_TPM(void);
void mutate_under_GSM(void);
void sequence_mutation_model(void);
int random_base(double *);
void setSubsRateMatrix();
void mutate_as_SNP(int);
void mutate_under_ISM(void);
void write_haplotype_file(void);
std::string convert_binary_format(std::string);
std::string convert_quartery_format(std::string);
int seqPairDiff(std::string, std::string);
double seqPairNucSubDiff(std::string, std::string);
void entete_fichier(void);
void ecriture_fichier(void);
void ecriture_fichier_genotypes(void);
void ecriture_fichier_coordonnees(void);
void ecriture_fichier_migrate(void);
void ecriture_fichiers_DG2002KAM(void);
void ecriture_fichiers_DG2002SMM(void);
void ecriture_fichier_moyennes(void);
void ecriture_disp(void);//ecriture des distances de dispersion "efficaces"
void ecriture_disp_moy(void);//ecriture des distances de dispersion "efficaces" moyennes sur repets
void ecriture_fichier_iterativeStats(void);
void ecriture_fichier_iterativeStats_cpp();
void calcul_proba_identite(void);
void ecriture_matrice_proba_id(void);
void calcul_proba_identite_matrice(void);
void calcul_coa_individuel(void);
void calcul_stat_disp(void);
void calcul_stat_disp_moy(void);
void calcul_VariousStats(void);
bool checkMinorAlleleFreq();
void freq_sample(void);
void var_allele(void);



double binom(double pp,int n);
double gammln(double xx);/*pour procedure loi de poisson*/
int loi_geom(const double var_geom);
int  mutation(long int temps);
//void nrerror(char error_text[]);
double poisson(double poiss);


extern bool cinGetOnError,pauseGP; //to pause on cerr messages in batchDebug mode; overridden by explicit call to Batch mode

extern std::string fichier_geneland_geno,fichier_geneland_coord,fichier_migrate,fichier_genepop,fichier_stepsim,SettingsFilename,SpecificDensityFilename,EdgeEffect,Genepopfile_extension;//,twoD_disp

extern bool variable_mubool,txt_extensionbool,genepopoui,genepopNoCoordbool,genelandoui,DG2002,AllStates,migrate_oui,migrate_lettre_oui,
            calculoui,HexpNeioui,Varoui,Fisoui,suiviQ,iterativeStats,dossier_nbre_allele,Prob_Id_Matrix,calcul_anc,effective_disp,seqStatsOui,
            const_disp,compareMath,predispbool,Specific_Density_Designbool,groupAllSamplesbool,
            random_translation,GlobalDensSpatiallyHeterog,GlobalDispSpatiallyHeterog,migraineSettingsbool,Simple1DProductDispBool;

extern int repet,Seeds,n_locus,min_allele,max_mutations,Kmin,Kmax,Kini,MotifSize,ploidy,/*sym, dans TVpars*/dim_reseau1,dim_reseau2;

extern std::vector< int > dim_sample1,dim_sample2, xmin_sample,ymin_sample,Spec_SampleSize,vide_sample,
                            vide_sampleY,vide_sampleX,dens_sample;
extern std::vector< bool > Specific_Sample_Designbool;

extern const int x,y;

extern double _mu,pSMM,var_geom_TPM,var_geom_GSM;

extern std::string model_mut;

extern std::vector< std::vector < std::vector<int> > > Spec_Sample_Coord;//stocke les coord des échantillons quand user defined

extern long double     *migra;
extern int 	   **densite;
extern int       coord_origine[4];

extern int dx_max,dy_max,dx_min,dy_min;

extern unsigned long int currentGeneration;

extern std::vector<std::vector<double> >backwardPhilo; //matrix of backward non-immigration rate // FR 0610 essential for controlling backward rates

//************************************************************************************************************
//	Vector datatypes for multi-locus genetic markers [0 -> n_locus] 0= default=KAM not used for simulation, 1->n_locus = markers simulated
//************************************************************************************************************
extern std::vector<int> minAlleleVector, maxMutationsVector, locusNumVector, kMinVector, kMaxVector, motifSizeVector, kIniVector;

extern std::vector<double> mutRateVector, minorAlleleFreqVector, pSMMVector, geomTPMVector, geomGSMVector;

extern std::vector<std::string> mutModelIDVector;

extern std::vector<bool> polyLociBoolVector, varMutRateBoolVector;
//************************************************************************************************************
// other markers related variables, but not vector datatype

//specifies if the the locus is of sequence type or SNP
extern bool specifyLocusSeq, specifyVarBaseFreqs;


extern int nMarker, seqSize;

extern std::string nexusFormat, defMRCASequence;

extern double varBaseFreqs[4], ratio_TITV, ratio_TITI;
//************************************************************************************************************

namespace NS_translation {
    extern int translationx,translationy,grande_dim;
}

namespace NS_coal_tree {
    extern int nbre_noeud_restant;
#ifdef STL_DEBUG
//    extern std::vector<std::vector<int> > coord_individus;
    extern std::vector<std::vector<int> >coord_noeud;
    extern std::vector<std::vector<int> >coord_ori;
    extern std::vector<int>no_noeud_restant;
//    extern std::vector<int>aleat_noeud;
#else
//    extern int 	  **coord_individus/*[1601][3]*/;//tableau des coord dans l'echantillon
    extern int 		**coord_noeud;
    extern int **coord_ori;
    extern int 		*no_noeud_restant;
//    extern int 		*aleat_noeud;
#endif
}

namespace NS_diagnostic_tables {
    extern long long int        **effective;
    extern long long int **cumulEffectiveImmigration, **effectiveImmigration;
}


extern int (*SortieSupfnPtr1)(int coord, int dim);
extern int (*SortieSupfnPtr2)(int coord, int dim, int trans);
extern int (*SortieInffnPtr1)(int coord, int dim);
extern int (*SortieInffnPtr2)(int coord, int dim, int trans);

extern std::ofstream simulpars;

#endif
