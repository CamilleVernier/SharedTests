#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <fstream> // string ... !.
#include <iostream> // string ... !.
#include <vector>
#include <cmath>


#include "TimeVarPars.h"
#include "MersenneTwister.h" // CB->RL can be removed
#include "genepop_extract.h"
#include "lattice_s.h"

using namespace std;

int CTimeVaryingParams::setGrowthPars(CTimeVaryingParams* nextTVpars) {
		dens=initialDens;
		dimRes1=initialDimRes1;
		dimRes2=initialDimRes2;
		Dt=(nextTVpars->minGen-minGen+1);//+1 a cause de G0=1.
        Ddens=nextTVpars->initialDens-initialDens;
		Rdens=((float) nextTVpars->initialDens)/((float) initialDens);
		DdimRes1=nextTVpars->initialDimRes1-initialDimRes1;
		DdimRes2=nextTVpars->initialDimRes2-initialDimRes2;
		RdimRes1=nextTVpars->initialDimRes1/initialDimRes1;
		RdimRes2=nextTVpars->initialDimRes2/initialDimRes2;
        if(cmp_nocase(ContDemeSizeChange,"Linear") == 0) {//continuous linear change in deme size through time
            growthRate=Ddens/Dt;
        } else if(cmp_nocase(ContDemeSizeChange,"Exponential") == 0) {//
			growthRate=log((float) Rdens)/Dt;
        } else if(cmp_nocase(ContDemeSizeChange,"Logistic") == 0) {//continuous logistic change in deme size through time
            KLogistic=nextTVpars->initialDens;
        }
		if(cmp_nocase(ContLatticeSizeChange,"Linear") == 0) {//continuous linear change in deme size through time
			growthRate1=DdimRes1/Dt;
			growthRate2=DdimRes2/Dt;
		} else if(cmp_nocase(ContLatticeSizeChange,"Exponential") == 0) {//
			growthRate1=log((float) RdimRes1)/Dt;
			growthRate2=log((float) RdimRes2)/Dt;
		} else if(cmp_nocase(ContLatticeSizeChange,"Logistic") == 0) {//continuous logistic change in deme size through time
			KLogistic1=nextTVpars->initialDimRes1;
			KLogistic2=nextTVpars->initialDimRes2;
		}
return(0);
}

int CTimeVaryingParams::currentDensAndLatticeSize(long unsigned int& locCurrentGeneration) {
	
	//cout << "currentDensAndLatticeSize : Gn=" << locCurrentGeneration << endl;
	//cout << "ContDemeSizeChange=" << ContDemeSizeChange << "; ContLatticeSizeChange=" << ContLatticeSizeChange << endl;
    
	if(cmp_nocase(ContDemeSizeChange,"Linear") == 0) {
        dens=(int) floor(initialDens + growthRate*(locCurrentGeneration-minGen));
		//cout << "growthrate=" << growthRate << endl;
    } else if(cmp_nocase(ContDemeSizeChange,"Exponential")==0) {
        dens=(int) floor(initialDens*exp(growthRate*(locCurrentGeneration-minGen)));
		//cout << "growthrate=" << growthRate << endl;
    } else if(cmp_nocase(ContDemeSizeChange,"Logistic")==0) {
        dens=(int) floor(KLogistic/(1+(KLogistic-initialDens)/initialDens*exp(-densLogisticGrowthRate*(locCurrentGeneration-minGen))));
		//cout << "densLogisticGrowthRate=" << densLogisticGrowthRate << "; KLogistic="  << KLogistic << endl;
    }

	if(cmp_nocase(ContLatticeSizeChange,"Linear") == 0) {
        dimRes1=(int) floor(initialDimRes1 + growthRate1*(locCurrentGeneration-minGen));
        dimRes2=(int) floor(initialDimRes2 + growthRate2*(locCurrentGeneration-minGen));
		//cout << "growthrate1=" << growthRate1 << "; growthrate2=" << growthRate2 << endl;
    } else if(cmp_nocase(ContLatticeSizeChange,"Exponential")==0) {
        dimRes1=(int) floor(initialDimRes1*exp(growthRate1*(locCurrentGeneration-minGen)));
        dimRes2=(int) floor(initialDimRes2*exp(growthRate2*(locCurrentGeneration-minGen)));
		//cout << "growthrate1=" << growthRate1 << "; growthrate2=" << growthRate2 << endl;
    } else if(cmp_nocase(ContLatticeSizeChange,"Logistic")==0) {
        dimRes1=(int) floor(KLogistic1/(1+(KLogistic1-initialDimRes1)/initialDimRes1*exp(-latticeLogisticGrowthRate*(locCurrentGeneration-minGen))));
        dimRes2=(int) floor(KLogistic2/(1+(KLogistic2-initialDimRes2)/initialDimRes2*exp(-latticeLogisticGrowthRate*(locCurrentGeneration-minGen))));
		//cout << "LatticeLogisticGrowthRate=" << LatticeLogisticGrowthRate << endl;
		//cout << "KLogistic1=" << KLogistic1 << "; KLogistic2=" << KLogistic2 << endl;
    }	
	//cout << "Gn=" << currentGeneration << " (minGen=" << minGen << endl;
	//cout << "LatticeSizeX=" << dimRes1 << "; LatticeSizeY" << dimRes2 << endl;
	//cout << "Gn=" << currentGeneration << " -> N=" << dens << endl;

	return(0);
}

void CTimeVaryingParams::oneDProductWithoutm0() {// arg should be the maximum realized immig rate among demes
   long double norm;
           if (Mig>1) {
                    cerr<<"(!) From oneDProductWithoutm0(): maxMig argument ="<<Mig<<"; value not feasible"<<endl;
                    if (cinGetOnError) cin.get();
                    exit(-1);
            }
            long double maxcumul=0,tmp;
            /*at this point in geometric model the migra[] terms are normalized to sum to 1 and decrease as powers of g
            so the norm recomputed here should be 1. maxcumul 2D below could approach 1-migra[0]^2 (except for edge effects at dist <dxmax) */
            norm=migra[0]/2;
            for(int i=1;i<=dx_max;i++) {
                if (isnan(norm)) {
                    cerr<<"(!) From oneDProductWithoutm0(): migra["<<i-1<<"] is NaN"<<endl;
                    if (cinGetOnError) cin.get();
                    exit(-1);
                }
                norm+=migra[i];
            }
            norm*=2;
            for(int i=0;i<=dx_max;i++) migra[i]/=norm; // the sum of the axial terms from -dxmax to dxmax, 0 included

            /****** We want to control the maximal dispersal rate among focal demes, we will change the migra[]: ******/
        { // bloc a,b
		      int a=0;int b; double cumul;
            if(dimRes2/vide>1) { // this test matches the one in the sampling algorithm
                for(int xpos=0;xpos<dimRes1;xpos++) {  // double loop over focal patches
                    for(int ypos=0;ypos<dimRes2;ypos++) {
                        b=0;cumul=0;
                        for(int xxpos=0;xxpos<dimRes1;xxpos++){
                            for(int yypos=0;yypos<dimRes2;yypos++){
                                if(a!=b) { // different pops => not the migra[0]*migra[0] term
								  // distances between a and b for absorbing boundaries:
								  cumul+=migra[abs(xpos-xxpos)]*migra[abs(ypos-yypos)];
                                }
                                b++;
                            }// yypos
                        }// xxpos
                        if (cumul>maxcumul) {
                            maxcumul=cumul; // NOT the maximal immigration proba into one position if the migra[] don't depend on Mig
						  //    cout<<maxcumul;getchar();
                        }
                        a++;
                    } // for ypos ...
                } // for xpos ...
            } else { // 1D...
			    for(int xpos=0;xpos<dimRes1;xpos++) {
					  b=0;cumul=0;
					  for(int xxpos=0;xxpos<dimRes1;xxpos++){
							  if(a!=b) { // different pops => not the migra[0]*migra[0] term
								  // distances between a and b for absorbing boundaries:
								  cumul+=migra[abs(xpos-xxpos)];
							  }
							  b++;
					  }// xxpos
					  if (cumul>maxcumul) {
						  maxcumul=cumul; // NOT the maximal immigration proba into one position  if the migra[] don't depend on Mig
					  }
					  a++;
			    } // for xpos ...
            } //2D else 1D
        } // block a ,b
           /*** now we compute tmp such that migra0/(migra0+maxcumul*tmp^d)=1-maxMig:=migra0
           ***/
           tmp=Mig/maxcumul;
		   if(dimRes2/vide>1) { // this test matches the one in the sampling algorithm
		       tmp=sqrt(tmp);
		   }
           /*** now we recompute immig rates and write them as a check***/
			  for(int i=0;i<=dx_max;i++) { //[0] matters for [0][!0] terms in 2D....
				  migra[i]*=tmp; /* now the sum of the 2D terms over feasible dispersal positions must be mig for the maximal deme
				  and the sum of the migra[] terms must be > sqrt(mig) (and somewhat opaque). */
			  }
			  if (! simulparsWrittenBool) {
                int a=0;int b; double cumul;
                double migra0=1.0-Mig;  //CONSTANT should be the final non-immgration rate, NOT used in the computation of norm

			      // recomputation of the the cumul: not economical, but safe check
			      a=0;
				  simulpars<<endl<<"immigration rates in the different demes:"<<endl;
			      if(dimRes2/vide>1) { // this test matches the one in the sampling algorithm
			        for(int xpos=0;xpos<dimRes1;xpos++) {
				      for(int ypos=0;ypos<dimRes2;ypos++) {
					      b=0;cumul=0;
//cout << endl << xpos << " " << ypos;getchar();
					      for(int xxpos=0;xxpos<dimRes1;xxpos++){
						      for(int yypos=0;yypos<dimRes2;yypos++){
							      if(a!=b) { // different pops => not the migra[0]*migra[0] term
								       // distances between a and b for absorbing boundaries:
								       cumul+=migra[abs(xpos-xxpos)]*migra[abs(ypos-yypos)];
							      }
							      b++;
                              }// yypos
					      }// xxpos
                          backwardPhilo[xpos][ypos]=migra0/(cumul+migra0); //FR 0610
/** At this point, for each focal position we have
cumul= \sum g^k *maxMig/max(\sum g^k)
backward immigration rates will be of the form
[\sum g^k *maxMig/max(\sum g^k)]/[\sum g^k *maxMig/max(\sum g^k)]+1-Mig)
and we could think that the relative immigration proba in each demes are the relative
\sum g^k values, but this is clearly not so since the denominators differ for each term.

There is an independent check of this code in (isolde_stat_ML)/technical_docs/checkIBDsim.nb.
**/
						  simulpars<<1-backwardPhilo[xpos][ypos]<<" ";
					      a++;
				      } // for ypos ...
					  simulpars<<endl;
			        } // for xpos ...
			     } else { //1D...
                    int ypos=0; //FR 0610
			        for(int xpos=0;xpos<dimRes1;xpos++) {
					      b=0;cumul=0;
					      for(int xxpos=0;xxpos<dimRes1;xxpos++){
							      if(a!=b) { // different pops => not the migra[0]*migra[0] term
								       // distances between a and b for absorbing boundaries:
								       cumul+=migra[abs(xpos-xxpos)];
							      }
							      b++;
					      }// xxpos
                          backwardPhilo[xpos][ypos]=migra0/(cumul+migra0); //FR 0610
						  simulpars<<1-backwardPhilo[xpos][ypos]<<" ";
					      a++;
					  simulpars<<endl;
			        } // for xpos ...
                 }
              }// if (! simulpars...
} //end def oneDProductWithoutm0

extern int cumul1(int coordOrigineX, int coordOrigineY);
extern void cumul2(int coordOrigineX, int coordOrigineY);


int CTimeVaryingParams::computeCumuls() {//
    if(Specific_Density_Designbool) getDensityFromFile();
    else for(int cc=1;cc<=dimRes1;cc++) {
                for(int dd=1;dd<=dimRes2;dd++) {
                    densite[cc][dd]=calcul_densite(cc,dd);
		        }
            }
//    for(int cc=1;cc<=dimRes1;cc++)
//        for(int dd=1;dd<=dimRes2;dd++) {
//            cout << "densite[" << cc << "][" << dd << "]=" << densite[cc][dd] << endl;
//        }
    
	for(int cc=1;cc<=dimRes1;cc++) {
	    for(int dd=1;dd<=dimRes2;dd++){
            cumul1(cc,dd);//calcul le demonimateur des probabilite de migration
				// (pour chaque points du reseau si pas de symetrie ou pour chaque point de la "symetrie")
            cumul2(cc,dd);//calcul le nominatuer et donc finalise les probabilite de migration
		}
    }
return(0);
}




void CTimeVaryingParams::realise_disp(vector<int> n_genes_loc) {
    int originex,originey,nbre_noeud_restant_loc;
    using namespace NS_coal_tree;
    using namespace NS_translation;
    using namespace NS_diagnostic_tables;
    if(!fixe) {//si il existe une source d'heterogeneite sur le reseau = noeuds vides ou zone; sinon dispfixe()
        
        //migration de chaque noeud restant
        if(currentGeneration==1) nbre_noeud_restant_loc=n_genes_loc[0]; else nbre_noeud_restant_loc=nbre_noeud_restant;
        for(int i=1;i<=nbre_noeud_restant_loc;i++){//fait la translation si necessaire
//            cout << "Avant disp: noeud " << i << ":" <<  coord_noeud[no_noeud_restant[i]][x] << " , " << coord_noeud[no_noeud_restant[i]][y] << endl;
            // RL : a mettre en mode DEBUG uniquement
            if( currentGeneration!=minGen && ( (coord_noeud[no_noeud_restant[i]][x]<=0) || (coord_noeud[no_noeud_restant[i]][x]>dimRes1)
               || (coord_noeud[no_noeud_restant[i]][y]<=0) || (coord_noeud[no_noeud_restant[i]][y]>dimRes2) ) ) {
                cout << "Pb0 : noeud i=" << i << "coord_noeud[no_noeud_restant[i]][x]=" << coord_noeud[no_noeud_restant[i]][x];
                cout << "coord_noeud[no_noeud_restant[i]][y]=" << coord_noeud[no_noeud_restant[i]][y] << endl;
            }
            originex=coord_origine[x]=SortieSupfnPtr2(coord_noeud[no_noeud_restant[i]][x],dimRes1,translationx);
            originey=coord_origine[y]=SortieSupfnPtr2(coord_noeud[no_noeud_restant[i]][y],dimRes2,translationy);
            //dispersion en x et y en fonction de la place du neoud sur le reseau (=originex et originey), calcul dispxx et dispyy
            dispxy(); // 	increments coord_origine[x]+=dispxx; coord_origine[y]+=dispyy;
            
            // RL : a mettre en mode DEBUG uniquement
            if(EdgeEffect=="abs" && ( (coord_origine[x]<=0) || (coord_origine[x]>dimRes1)
                                     || (coord_origine[y]<=0) || (coord_origine[y]>dimRes2)) ) {
                cout << "Pb1 : coord_origine[x]=" << coord_origine[x];
                cout << "coord_origine[y]=" << coord_origine[y] << endl;
            }
            
            //ajuste les nouvelle coordonnees en fonction des effets de bords
            coord_ori[no_noeud_restant[i]][x]=SortieSupfnPtr1(coord_origine[x],dimRes1);
            coord_ori[no_noeud_restant[i]][x]=SortieInffnPtr1(coord_ori[no_noeud_restant[i]][x],dimRes1);
            if(EdgeEffect=="refl") {//pour reflection multiples
                while((coord_ori[no_noeud_restant[i]][x]<=0)||(coord_ori[no_noeud_restant[i]][x]>dimRes1)){
                    coord_ori[no_noeud_restant[i]][x]=SortieInffnPtr1(coord_ori[no_noeud_restant[i]][x],dimRes1);
                    coord_ori[no_noeud_restant[i]][x]=SortieSupfnPtr1(coord_ori[no_noeud_restant[i]][x],dimRes1);
                }
            }
            coord_ori[no_noeud_restant[i]][y]=SortieSupfnPtr1(coord_origine[y],dimRes2);
            coord_ori[no_noeud_restant[i]][y]=SortieInffnPtr1(coord_ori[no_noeud_restant[i]][y],dimRes2);
            if(EdgeEffect=="refl") {//pour reflection multiples
                while((coord_ori[no_noeud_restant[i]][y]<=0)||(coord_ori[no_noeud_restant[i]][y]>dimRes2)) {
                    coord_ori[no_noeud_restant[i]][y]=SortieInffnPtr1(coord_ori[no_noeud_restant[i]][y],dimRes2);
                    coord_ori[no_noeud_restant[i]][y]=SortieSupfnPtr1(coord_ori[no_noeud_restant[i]][y],dimRes2);
                }
            }
            // RL : a mettre en mode DEBUG uniquement
            if( (coord_ori[no_noeud_restant[i]][x]<=0) || (coord_ori[no_noeud_restant[i]][x]>dimRes1)
               || (coord_ori[no_noeud_restant[i]][y]<=0) || (coord_ori[no_noeud_restant[i]][y]>dimRes2) ) {
                cout << "Pb11 : noeud i=" << i << "ccoord_ori[no_noeud_restant[i]][x]=" << coord_ori[no_noeud_restant[i]][x];
                cout << "coord_ori[no_noeud_restant[i]][y]=" << coord_ori[no_noeud_restant[i]][y] << endl;
            }

            if(backwardBarrierX && ( (originey > y1_barrier) && (originey < y2_barrier) ) ) {
                double alea2=rand();//RL -> RL a modifier, mais alea() est galere à utiliser ici....
                if(((originex < x1_barrier) && (coord_ori[no_noeud_restant[i]][x] > x1_barrier)) && (alea2 > barrierCrossingRateUp) ) coord_ori[no_noeud_restant[i]][x]=originex;
                else if(((originex > x1_barrier) && (coord_ori[no_noeud_restant[i]][x] < x1_barrier)) && (alea2 > barrierCrossingRateDown) ) coord_ori[no_noeud_restant[i]][x]=originex;
            }
            if(backwardBarrierY &&  ( (originex > x1_barrier) && (originex < x2_barrier) ) ) {
                double alea2=rand();//RL -> RL a modifier, mais alea() est galere à utiliser ici....
                // si barrier alors on revient ou pas au point de départ en fonction de la "permeabilite" de la barriere RL -> RL a ameillorer en reflectif.
                if(((originey < y1_barrier) && (coord_ori[no_noeud_restant[i]][y] > y1_barrier)) && (alea2 > barrierCrossingRateUp) ) coord_ori[no_noeud_restant[i]][y]=originey;
                else if(((originey > y1_barrier) && (coord_ori[no_noeud_restant[i]][y] < y1_barrier)) && (alea2 > barrierCrossingRateDown) ) coord_ori[no_noeud_restant[i]][y]=originey;
            }
            
            
            if(effective_disp) {//calcul de la dispersion "efficace"
                effective[grande_dim+(coord_ori[no_noeud_restant[i]][x]-originex)][x]++;
                effective[grande_dim+ (coord_ori[no_noeud_restant[i]][y]-originey)][y]++;
                // RL : a mettre en mode DEBUG uniquement
                if(EdgeEffect=="abs" && ( (fabs(coord_ori[no_noeud_restant[i]][x]-originex)>dx_max)||(fabs(coord_ori[no_noeud_restant[i]][y]-originey)>dy_max) ) ) {
                    cout << "\
                    nPb2 : (coord_ori[no_noeud_restant[i]][x]-originex)=" << (coord_ori[no_noeud_restant[i]][x]-originex) << endl;
                    cout << "Pb2 : (coord_ori[no_noeud_restant[i]][y]-originey)=" << (coord_ori[no_noeud_restant[i]][y]-originey) << endl;
                    //getchar();
                }

                if((coord_ori[no_noeud_restant[i]][x]-originex)!=0 || (coord_ori[no_noeud_restant[i]][y]-originey)!=0)
                    effectiveImmigration[originex][originey]++;
                cumulEffectiveImmigration[originex][originey]++;
            }
//            cout << "apres disp : " <<  coord_ori[no_noeud_restant[i]][x] << " , " << coord_ori[no_noeud_restant[i]][y] << endl;
//            cout << endl;
//            cin.get();
        }//fin de la migration de chaque noeud
        //attribut les nouvelles coordonnees
        for(int i=1;i<=nbre_noeud_restant_loc;i++) {
            coord_noeud[no_noeud_restant[i]][x]=coord_ori[no_noeud_restant[i]][x];
            coord_noeud[no_noeud_restant[i]][y]=coord_ori[no_noeud_restant[i]][y];
            // RL : a mettre en mode DEBUG uniquement
            if( (coord_noeud[no_noeud_restant[i]][x]<=0) || (coord_noeud[no_noeud_restant[i]][x]>dimRes1)
               || (coord_noeud[no_noeud_restant[i]][y]<=0) || (coord_noeud[no_noeud_restant[i]][y]>dimRes2) ) {
                cout << "Pb22 : noeud i=" << i << "coord_noeud[no_noeud_restant[i]][x]=" << coord_noeud[no_noeud_restant[i]][x];
                cout << "coord_noeud[no_noeud_restant[i]][y]=" << coord_noeud[no_noeud_restant[i]][y] << endl;
            }

        }

    } else {
        disp_fixe();
    }//si reseau totalement homogene = ni vide ni zone = beaucoup plus simple et beaucoup plus rapide!!
} // end realise_disp();

/*****************/
//calcul la densite en un point du reseau selon le model demo
int CTimeVaryingParams::calcul_densite(int coordx,int coordy) {
    int densite_locale;//variable locale de la densite
    int vide_locale;//variable locale de vide

    if(zone==1) {
        if((xmin_zone<=coordx)&&(coordx<=xmax_zone)&&(ymin_zone<=coordy)&&(coordy<=ymax_zone)) {
            if(vide_zone==1) return(dens_zone);
            else {densite_locale=dens_zone;vide_locale=vide_zone;}
        } else {densite_locale=dens;vide_locale=vide;}
    } else {densite_locale=dens;vide_locale=vide;}

    if(vide_locale!=1) {
        if((coordx % vide_locale)==0 &&(coordy % vide)==0) return(densite_locale);
        else return(0);
    } else return(densite_locale);
}/*fin calcul densite*/

/*****************/
//affichage des differentes phases démo
int CTimeVaryingParams::PrintTVpars() {

	cout << "Current Density (nbr of genes)=" << ploidy*dens << endl;
	cout << "Initial Density for this phase=" << ploidy*initialDens << endl;
	cout << "dimRes1=" << dimRes1 << " and dimRes2=" << dimRes2 << endl;
	cout << "voidNode=" << vide << endl;
    cout << "zone=" << zone << endl;
	cout << "DensSpatiallyHeterog =" << DensSpatiallyHeterog << endl;
    cout << "DispSpatiallyHeterog=" << DispSpatiallyHeterog << endl;
	cout << "Specific_Density_Designbool=" << Specific_Density_Designbool << endl;

	if(zone && DensSpatiallyHeterog) {
			cout << "xmin_zone=" << xmin_zone << " and xmax_zone=" << xmax_zone << endl;
			cout << "ymin_zone=" << ymin_zone << " and ymax_zone=" << ymax_zone << endl;
			cout << "voidNode in Zone =" << vide_zone << endl;
			cout << "Density in Zone (nbr of genes) =" << ploidy*dens_zone << endl;
	}
	if(backwardBarrier || forwardBarrier) {
        if(forwardBarrier) simulpars<<"Presence of a real forward barrier to gene flow with coordinates :" <<endl;
        else simulpars<<"Presence of an approximate (i.e. backward) barrier to gene flow with coordinates :" <<endl;
		cout << "x1_barrier=" << x1_barrier << " and x2_barrier=" << x2_barrier << endl;
		cout << "y1_barrier=" << y1_barrier << " and y2_barrier=" << y2_barrier << endl;
		if(barrierCrossingRateUp==barrierCrossingRateDown) cout << "barrierCrossingRate=" << barrierCrossingRateUp << endl;
        else {
         cout << "barrierCrossingRateUp=" << barrierCrossingRateUp << endl;
         cout << "barrierCrossingRateDown=" << barrierCrossingRateDown << endl;
        }
	}
	cout << "dispersal distribution =" << mod << endl;
	if(mod=='g' || mod=='P' || mod=='S') {
        cout << "Emmigration rate=" << Mig << endl;
        cout << "dist_max=" << dist_max << endl;
    }
	if(mod=='g') cout << "GeoG=" << GeoG << endl;
	if(mod=='P') cout << "Pareto_Shape=" << Pareto_Shape << endl;
	if(mod=='S') {
		cout << "Sichel_Gamma=" << Sichel_Gamma << endl;
		cout << "Sichel_Xi=" << Sichel_Xi << endl;
		cout << "Sichel_Omega=" << Sichel_Omega << endl;
	}
    cout << "dx_max=" << dx_max << endl;
	if(zone && DispSpatiallyHeterog) {
        cout << "xmin_zone=" << xmin_zone << " and xmax_zone=" << xmax_zone << endl;
        cout << "ymin_zone=" << ymin_zone << " and ymax_zone=" << ymax_zone << endl;
        cout << "dispersal distribution in the \"zone\" =" << mod_zone << endl;
        if(mod_zone=='g' || mod_zone=='P' || mod_zone=='S') {
         cout << "Emmigration rate in the \"zone\"=" << Mig_zone << endl;
         cout << "dist_max_zone (! not dx_max...)=" << dist_max_zone << endl;
        }
        if(mod_zone=='g') cout << "GeoG_zone=" << GeoG_zone << endl;
        if(mod_zone=='P') cout << "Pareto_Shape_zone=" << Pareto_Shape_zone << endl;
        if(mod_zone=='S') {
            cout << "Sichel_Gamma_zone=" << Sichel_Gamma_zone << endl;
            cout << "Sichel_Xi_zone=" << Sichel_Xi_zone << endl;
            cout << "Sichel_Omega_zone=" << Sichel_Omega_zone << endl;
        }
        cout << "dx_max=" << dx_max << endl;
	}
	if(ContDemeSizeChange!="None") {
		cout << "ContDemeSizeChange=" << ContDemeSizeChange << endl;
		cout << "minGen=" << minGen << endl;
		cout << "Dt=" << Dt << endl;
		cout << "Ddens=" << Ddens << endl;
		cout << "Rdens=" << Rdens << endl;
		cout << "growthRate=" << growthRate << endl;
		cout << "densLogisticGrowthRate=" << densLogisticGrowthRate << endl;
	}
	if(ContLatticeSizeChange!="None") {
		cout << "ContLatticeSizeChange=" << ContLatticeSizeChange << endl;
		cout << "minGen=" << minGen << endl;
		cout << "Dt=" << Dt << endl;
		cout << "DdimRes1=" << DdimRes1 << endl;
		cout << "DdimRes2=" << DdimRes2 << endl;
		cout << "RdimRes1=" << RdimRes1 << endl;
		cout << "RdimRes2=" << RdimRes2 << endl;
		cout << "growthRate1=" << growthRate1 << endl;
		cout << "growthRate2=" << growthRate2 << endl;
		cout << "latticeLogisticGrowthRate=" << latticeLogisticGrowthRate << endl;
	}
return 0;
}/*fin PrintTVpars*/


