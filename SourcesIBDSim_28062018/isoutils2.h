/***************************************************************************
� F. Rousset 2005- for code collection
francois.rousset@univ-montp2.fr

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
#ifndef H_ISOUTILS
#define H_ISOUTILS

void mysleep(unsigned int mseconds);
void nrerror(std::string error_text);
void MyCopyFile(char acopier[80], char copie[80]);
//void readln(FILE *toto);
int batcommande(char contenu[80]);
int cmdcommande(char contenu[80]);
double **dmatrix(int nrh, int nch);
void free_dmatrix(double **m, int nrh);
int **imatrix(int nrh, int nch);
void free_imatrix(int **m, int nrh);
long long int **llimatrix(int nrh, int nch); // used; lli>i
void free_llimatrix(long long int **m, int nrh);
int *ivector(int nh);
void free_ivector(int *v);
long int *livector(int nh);
void free_livector(long int *v);
long unsigned int *luivector(int nh);
void free_luivector(long unsigned int *v);
double *dvector(int nh);
void free_dvector(double *v);
long double *ldvector(int nh);
void free_ldvector(long double *v);
int ***itab3(int nrh, int nch,int nh);
void free_itab3(int ***m,int nrh,int nch);
double ***dtab3(int nrh, int nch,int nh);
void free_dtab3(double ***m,int nrh,int nch);
int **istructure(int nr,int nc);
void free_istructure(int **struc, int nr);
std::vector< double > getLinearFit(const std::vector<double>& X,const std::vector<double>& Y);
double slopeLinReg(const std::vector<double> & X,const std::vector< double >& Y);


#endif
