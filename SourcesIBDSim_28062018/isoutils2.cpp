/***************************************************************************
This contains small bits of code from Numerical Recipes in C by Press et al.

© F. Rousset 2005- for code collection
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

// si il dit qu'il trouve pas la fonction c'est peut Ítre que (int) doit Ítre remplacÈ par (long int)

/*readln, int to ansi et gestion de pointers*/
#include <iostream> // cout...
#include <string>
#include <cstdlib> //exit()...
#include <cstdio>
#include <vector>
#include <cmath>
#include <numeric>
#include "isoutils2.h"
#include "lattice_s.h"

using namespace std;

#include <time.h>

void mysleep(unsigned int mseconds)
{
    clock_t goal = mseconds + clock();
    while (goal > clock());
}


/**********************************************************/
void nrerror(string error_text) {
    cout<<"You have probably run out of memory...\n";
    cout<<error_text<<endl;
    cout<<"...the program must terminate...\n";
    cin.get();
    exit(-1);
}
/***********************************************************/
/***********************************************************/

/**********************************************************/
/**********************************************************
void readln(FILE *toto) {
  char c; while ((c = getc(toto)) != '\n') c=c;  // bumpkin's method to avoid boring Warnings...
	  }
*********************************************************/


/************************************************************/
/************************************************************************************/
int batcommande(char contenu[80])
{FILE *command_bat;
int success;

if((command_bat=fopen("commande.bat","w"))==0)
  {printf("batcommande cannot open file commande.bat");exit(1);}

fprintf(command_bat,"%s",contenu);

if(fclose(command_bat)) printf("commande.bat file close error");

if(system("commande.bat")==1) {success=1; /*printf("%d",success);
getchar();getchar();*/
}
 else {success=0;/*printf("%d",success);
getchar();getchar();*/
}

remove("commande.bat");

return success;


}

/**********************************************************************/
/************************************************************************************/
int cmdcommande(char contenu[80])
{FILE *command_cmd;
int success;

if((command_cmd=fopen("commande.cmd","w"))==0)
  {printf("cmdcommande cannot open file commande.bat");exit(1);}

fprintf(command_cmd,"%s",contenu);

if(fclose(command_cmd)) printf("commande.cmd file close error");

if(system("commande")==1) {success=1; /*printf("%d",success);
getchar();getchar();*/
}
 else {success=0;
       //printf("\n problem with cmdcommande");
       //getchar();
       //getchar();
}

//remove("commande.cmd");

return success;


}

/**********************************************************************/

/*************************************************************/
double **dmatrix(int nrh, int nch)  {
    int i;
    double **m;

//    m=(double **) malloc((unsigned) (nrh)*sizeof(double*));
	m=new double*[nrh];
    if (!m) nrerror("allocation failure 1 in dmatrix()");

    for (i=0; i<=nrh-1; i++) {
		m[i]=new double[nch];
//        m[i]=(double *) malloc((unsigned) (nch)*sizeof(double));
        if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
    }
    return m;
}  /*  end dmatrix  */

void free_dmatrix(double **m, int nrh)
{
int i;
for(i=nrh-1;i>=0;i--) delete[] m[i];
delete[] m;
}
/*******************************************************************/

/*******************************************************************/
/*double **dmatrixl(long int nrh, long int nch); // was used for Q1
double **dmatrixl(long int nrh, long int nch)  {
    long int i;
    double **m;

    m=(double **) malloc((unsigned) (nrh)*sizeof(double*));
	 if (!m) nrerror("allocation failure 1 in dmatrixl()");

    for (i=0; i<=nrh-1; i++) {
		  m[i]=(double *) malloc((unsigned) (nch)*sizeof(double));
        if (!m[i]) nrerror("allocation failure 2 in dmatrixl()");
    }
    return m;
}  //  end dmatrixl

void free_dmatrixl(double **m, long int nrh);
void free_dmatrixl(double **m, long int nrh)
{
long int i;
for(i=nrh-1;i>=0;i--) free((m[i]));
free((m));
}*/
/**********************************************************************/

/**********************************************************************/
/*long double **ldmatrix(int nrh, int nch); //was not used
long double **ldmatrix(int nrh, int nch)  {
    int i;
    long double **m;
    void nrerror(char error_text[]);

    m=(long double **) malloc((unsigned) (nrh)*sizeof(long double*));
    if (!m) nrerror("allocation failure 1 in dmatrix()");

    for (i=0; i<=nrh-1; i++) {
        m[i]=(long double *) malloc((unsigned) (nch)*sizeof(long double));
        if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
    }
    return m;
}  //  end ldmatrix

void free_ldmatrix(long double **m, int nrh);
void free_ldmatrix(long double **m, int nrh)
{
int i;
for(i=nrh-1;i>=0;i--) free((m[i]));
free((m));
}*/
/***********************************************************************/

/****************************************************************/
/*float **fmatrix(int nrh, int nch); //was not used
float **fmatrix(int nrh, int nch)  {
    int i;
    float **m;
    void nrerror(char error_text[]);

    m=(float **) malloc((unsigned) (nrh)*sizeof(float*));
    if (!m) nrerror("allocation failure 1 in fmatrix()");

    for (i=0; i<=nrh-1; i++) {
        m[i]=(float *) malloc((unsigned) (nch)*sizeof(float));
        if (!m[i]) nrerror("allocation failure 2 in fmatrix()");
    }
    return m;
}  //  end fmatrix

void free_fmatrix(float **m, int nrh);
void free_fmatrix(float **m, int nrh){
int i;
for(i=nrh-1;i>=0;i--) free((m[i]));
free((m));
}*/
/*******************************************************************/

/*******************************************************************/
/*float **fmatrixl(long int nrh, long int nch); //was not used
float **fmatrixl(long int nrh, long int nch)  {
    long int i;
    float **m;
    void nrerror(char error_text[]);

    m=(float **) malloc((unsigned) (nrh)*sizeof(float*));
    if (!m) nrerror("allocation failure 1 in fmatrixl()");

    for (i=0; i<=nrh-1; i++) {
        m[i]=(float *) malloc((unsigned) (nch)*sizeof(float));
        if (!m[i]) nrerror("allocation failure 2 in fmatrixl()");
    }
    return m;
}  //  end fmatrix


void free_fmatrixl(float **m, long int nrh);
void free_fmatrixl(float **m, long int nrh)
{
long int i;
for(i=nrh-1;i>=0;i--) free((m[i]));
free((m));
}*/
/********************************************************************/

/****************************************************************/
int **imatrix(int nrh, int nch)  {
    int i;
    int **m;

//  m=(int **) malloc((unsigned) (nrh)*sizeof(int*));
	m=new int*[nrh];
    if (!m) nrerror("allocation failure 1 in imatrix()");

    for (i=0; i<=nrh-1; i++) {
//        m[i]=(int *) malloc((unsigned) (nch)*sizeof(int));
		  m[i]=new int[nch];
        if (!m[i]) nrerror("allocation failure 2 in imatrix()");
    }
    return m;
}  /*  end imatrix  */

void free_imatrix(int **m, int nrh)
{
int i;
for(i=nrh-1;i>=0;i--) delete[] m[i];
delete[] m;
}
/*****************************************************************/

/****************************************************************/
/*int **imatrixl(long int nrh, long int nch); //was not used
int **imatrixl(long int nrh, long int nch)  {
    long int i;
    int **m;
    void nrerror(char error_text[]);

    m=(int **) malloc((unsigned) (nrh)*sizeof(int *));
	 if (!m) nrerror("allocation failure 1 in imatrixl()");

    for (i=0; i<=nrh-1; i++) {
        m[i]=(int *) malloc((unsigned) (nch)*sizeof(int));
        if (!m[i]) nrerror("allocation failure 2 in imatrixl()");
    }
    return m;
}

void free_imatrixl(int **m, long int nrh);
void free_imatrixl(int **m, long int nrh)
{
long int i;
for(i=nrh-1;i>=0;i--) free((m[i]));
free((m));
}*/
/********************************************************************/
/****************************************************************/
/*long int **limatrix(int nrh, int nch); // was not used... li = i...
long int **limatrix(int nrh, int nch)  {
    int i;
    long int **m;
    void nrerror(char error_text[]);

    m=(long int **) malloc((unsigned) (nrh)*sizeof(long int*));
    if (!m) nrerror("allocation failure 1 in imatrix()");

    for (i=0; i<=nrh-1; i++) {
        m[i]=(long int *) malloc((unsigned) (nch)*sizeof(long int));
        if (!m[i]) nrerror("allocation failure 2 in imatrix()");
    }
    return m;
}  //  end limatrix

void free_limatrix(long int **m, int nrh);
void free_limatrix(long int **m, int nrh)
{
int i;
for(i=nrh-1;i>=0;i--) free((m[i]));
free((m));
}*/
/*****************************************************************/
/****************************************************************/
long long int **llimatrix(int nrh, int nch)  {
    int i;
    long long int **m;

//    m=(long long int **) malloc((unsigned) (nrh)*sizeof(long long int*));
	m=new long long int*[nrh];
    if (!m) nrerror("allocation failure 1 in imatrix()");

    for (i=0; i<=nrh-1; i++) {
//        m[i]=(long long int *) malloc((unsigned) (nch)*sizeof(long long int));
			m[i]=new long long int[nch];
        if (!m[i]) nrerror("allocation failure 2 in imatrix()");
    }
    return m;
}  //  end llimatrix

void free_llimatrix(long long int **m, int nrh)
{
int i;
for(i=nrh-1;i>=0;i--) delete[] m[i];
delete[] m;
}
/*****************************************************************/

/***************************************************************/
/*char *cvector(int nh); // not used
char *cvector(int nh) {
    char *v;
    void nrerror(char error_text[]);

    v = (char *)malloc((unsigned) (nh+1)*sizeof(char));
    if (!v) nrerror("allocation failure in ivector()");
    return v;
}

void free_cvector(char *v);
void free_cvector(char *v)
{
free((v));
}*/
/************************************************************/

/************************************************************/
int *ivector(int nh) {
    int *v;
	v=new int[nh];
//    v = (int *)malloc((unsigned) (nh)*sizeof(int));
    if (!v) nrerror("allocation failure in ivector()");
    return v;
}


void free_ivector(int *v)
{delete[] v;}
/*********************************************************/
/************************************************************/
long int *livector(int nh) {
    long int *v;
	v=new long int[nh];
//    v = (long int *)malloc((unsigned) (nh)*sizeof(long int));
    if (!v) nrerror("allocation failure in ivector()");
    return v;
}


void free_livector(long int *v)
{delete[] v;}
/*********************************************************/
/************************************************************/
long unsigned int *luivector(int nh) {
    long unsigned int *v;
	v=new long unsigned int[nh];
    //    v = (long int *)malloc((unsigned) (nh)*sizeof(long int));
    if (!v) nrerror("allocation failure in ivector()");
    return v;
}


void free_luivector(long unsigned int *v)
{delete[] v;}
/*********************************************************/


/********************************************************/
double *dvector(int nh) {
    double *v;
	v=new double[nh];
//  v = (double *)malloc((unsigned) (nh)*sizeof(double));
    if (!v) nrerror("allocation failure in dvector()");
    return v;
}

void free_dvector(double *v)
{delete[] v;}
/*******************************************************/

/********************************************************/
long double *ldvector(int nh) {
    long double *v;
	v=new long double[nh];
//  v = (long double *)malloc((unsigned) (nh)*sizeof(long double));
    if (!v) nrerror("allocation failure in dvector()");
    return v;
}

void free_ldvector(long double *v)
{delete[] v;}
/*******************************************************/

/*******************************************************/
/*oat *fvector(int nh);
float *fvector(int nh) {
    float *v;
    void nrerror(char error_text[]);

    v = (float *)malloc((unsigned) (nh)*sizeof(float));
	 if (!v) nrerror("allocation failure in fvector()");
	 return v;
}

void free_fvector(float *v);
void free_fvector(float *v)
{
free((v));
}*/
/********************************************************/

/********************************************************/
/*oat *fvectorl(long int nh);
float *fvectorl(long int nh) {
	 float *v;
	 void nrerror(char error_text[]);

	 v = (float *)malloc((unsigned) (nh)*sizeof(float));
	 if (!v) nrerror("allocation failure in fvectorl()");
	 return v;
}

void free_fvectorl(float *v);
void free_fvectorl(float *v)
{
free((v));
}*/
/**********************************************************/

/**************************************************/
/*long int *livectorl(long int nh);
long int *livectorl(long int nh) {
	 long int *v;
	 void nrerror(char error_text[]);
	 v = (long int *)malloc((unsigned) (nh)*sizeof(long int));
	 if (!v) nrerror("allocation failure in livectorl()");
	 return v;
}

void free_livectorl(long int *v);
void free_livectorl(long int *v)
{
free((v));
}*/
/**************************************************/

/**************************************************/
/*int *ivectorl(long int nh);
int *ivectorl(long int nh) {
	 int *v;
	 void nrerror(char error_text[]);

	 v = (int *)malloc((unsigned) (nh)*sizeof(int));
	 if (!v) nrerror("allocation failure in ivectorl()");
	 return v;
}

void free_ivectorl(int *v);
void free_ivectorl(int *v)
{
free((v));
}*/
/**************************************************/

/**************************************************/
/*double *dvectorl(long int nh);
double *dvectorl(long int nh) {
	 double *v;
	 void nrerror(char error_text[]);

	 v = (double *)malloc((unsigned) (nh)*sizeof(double));
	 if (!v) nrerror("allocation failure in dvectorl()");
	 return v;
}

void free_dvectorl(double *v);
void free_dvectorl(double *v)
{
free((v));
}*/
/**************************************************/
/*******************************************************************/
int ***itab3(int nrh, int nch,int nh)  {
	 int i,j;
	 int ***m;
//	m=new int**[nrh];
	 m=(int ***) malloc((unsigned) (nrh)*sizeof(int**));
	 if (!m) nrerror("allocation failure 1 in dtab3()");

	 for (i=0; i<=nrh-1; i++) {
	 	m[i]=new int*[nch];
//		  m[i]=(int **) malloc((unsigned) (nch)*sizeof(int*));
		  if (!m[i]) nrerror("allocation failure 2 in dtab3()");
		  for(j=0;j<=nch-1;j++)
			{//m[i][j]=(int *) malloc((unsigned) (nh)*sizeof(int));
		 	m[i][j]=new int[nh];
			 if(!m[i][j]) nrerror("allocation failure 3 in dtab3()");
			}
	 }

    return m;
}  //  end itab3

void free_itab3(int ***m,int nrh,int nch)
{
int i,j;
for(i=nrh-1;i>=0;i--)
{for(j=nch-1;j>=0;j--) delete[] (m[i][j]);
 delete[] m[i];
}
free(m);
}
/************************************************************************************/


/*******************************************************************/
double ***dtab3(int nrh, int nch,int nh)  {
	 int i,j;
	 double ***m;

	 m=(double ***) malloc((unsigned) (nrh)*sizeof(double**));
	 if (!m) nrerror("allocation failure 1 in dtab3()");

	 for (i=0; i<=nrh-1; i++) {
		  m[i]=(double **) malloc((unsigned) (nch)*sizeof(double*));
		  if (!m[i]) nrerror("allocation failure 2 in dtab3()");
		  for(j=0;j<=nch-1;j++)
			{m[i][j]=(double *) malloc((unsigned) (nh)*sizeof(double));
			 if(!m[i][j]) nrerror("allocation failure 3 in dtab3()");
			}
	 }

    return m;
} //  end dtab3

void free_dtab3(double ***m,int nrh,int nch) {
    int i,j;
    for(i=nrh-1;i>=0;i--) {
        for(j=nch-1;j>=0;j--) free((m[i][j]));
        free((m[i]));
    }
    free((m));
}
/************************************************************************************/

/*********************************************************************/
int **istructure(int nr,int nc)
	{int **struc,i;

	 struc=(int **)malloc((unsigned) (nr)*sizeof(int*));
	 if (!struc) nrerror("allocation failure in istructure()");

	 for(i=0;i<nr;i++){
		struc[i]= (int*) malloc((unsigned) (nc)*sizeof(int));
		if (!struc[i]) nrerror("allocation failure in istructure()");
		}
	 return struc;
	 }


void free_istructure(int **struc, int nr)
{
int i;
for(i=nr-1;i>=0;i--) free((struc[i]));
free((struc));
}
/**************************************************************************/

vector< double > getLinearFit(const vector<double>& X,const vector<double>& Y)
{
    double xSum = 0, ySum = 0, xxSum = 0, xySum = 0, slope, intercept;
    
    if(X.size() != Y.size()) {
        cerr << "in getLinearFit() X and Y do not have the same size..." << endl;
        if (cinGetOnError) cin.get();
        exit(-1);

    }

    for (size_t i = 0; i < X.size(); i++) {
        xSum += X[i];
        ySum += Y[i];
        xxSum += X[i] * X[i];
        xySum += X[i] * Y[i];
    }
    
    slope = (X.size() * xySum - xSum * ySum) / (X.size() * xxSum - xSum * xSum);
    intercept = (ySum - slope * xSum) / X.size();
    
    if(fabs(slope - slopeLinReg(X,Y) )> 0.000001) {
        cerr << "in getLinearFit(), n=" << X.size() << "  Two slopes are different:" << slope << " " << slopeLinReg(X,Y) <<  endl;
        if (cinGetOnError) cin.get();
        exit(-1);

    }
    

    
    vector<double> res;
    res.push_back(slope);
    res.push_back(intercept);
    return res;
}

double slopeLinReg(const vector<double> & X,const vector< double >& Y) {

    size_t n=X.size();
    double avgX = accumulate(X.begin(), X.end(), 0.0) / n;
    double avgY = accumulate(Y.begin(), Y.end(), 0.0) / n;
    
    double numerator = 0.0;
    double denominator = 0.0;
    
    for(int i=0; i<n; ++i){
        numerator += (X[i] - avgX) * (Y[i] - avgY);
        denominator += (X[i] - avgX) * (X[i] - avgX);
    }
    
    if(denominator == 0){
        cerr << "in slope(), denominator is 0" << endl;
        if (cinGetOnError) cin.get();
        exit(-1);
    }
    
    return numerator / denominator;
}
