#ifndef TIMEVARPARS_H
#define TIMEVARPARS_H

#include <cstring>
#include <vector>


class CTimeVaryingParams {
    public:
    int phase; // phase demo
    bool simulparsWrittenBool;
    int 	dens;//nbre d'individus diploides par noeud (non vide!)
	int 	initialDens;//nbre d'individus diploides par noeud (non vide!) ayu debut de la phase demo
	int		initialDimRes1;//Lattice size at the begiining of the demographic phase
	int		initialDimRes2;//Lattice size at the begiining of the demographic phase
    int 	dimRes1; // doit etre <=  dim_reseau1
    int 	dimRes2; //1...doit etre <=  dim_reseau2
    bool 	fixe;//=0 si il existe des noeud vides sur le reseau, ou des heterogeneites spatialles de densité ou dispersion
    int 	vide;//1 si pas de noeuds vides 2,3 si il existe des noeuds vide ( 1sur 2, 1 sur 3)
    bool 	zone;//1 si il existe une zone de densite differente du reste
    bool    DensSpatiallyHeterog;/*0/1 si densité constante en tout point du reseau*/
	int		xmin_zone,xmax_zone,ymin_zone,ymax_zone;// coord de la zone
    int 	vide_zone;//si il y a des noeuds vides dans la zone
    int 	dens_zone;//nbre d'individus par noeud dans la zone
	bool	backwardBarrier,forwardBarrier,backwardBarrierX,backwardBarrierY,forwardBarrierX,forwardBarrierY;
	int		x1_barrier,x2_barrier,y1_barrier,y2_barrier;
	double	barrierCrossingRateUp,barrierCrossingRateDown;
    char 	mod; /*S ou g... Modele de dispersion *chercher cette option pour en modif les params*/
    double	Mig,GeoG,Pareto_Shape;
	int 	dist_max;//distance maximale de dispersion
    double  Sichel_Gamma,Sichel_Xi,Sichel_Omega;
    bool    DispSpatiallyHeterog;/*0/1 si dispersion constante en tout point du reseau*/
    char 	mod_zone; /*S ou g... Modele de dispersion *chercher cette option pour en modif les params*/
    double	Mig_zone,GeoG_zone,Pareto_Shape_zone;
	int 	dist_max_zone;//distance maximale de dispersion
    double Sichel_Gamma_zone,Sichel_Xi_zone,Sichel_Omega_zone;
    std::string ContDemeSizeChange,ContLatticeSizeChange;
    double densLogisticGrowthRate,latticeLogisticGrowthRate;
    unsigned long int minGen; //starting time of the demographic phase, starting from recent=0
    CTimeVaryingParams()  {
        simulparsWrittenBool=false;
        dens=1;
		initialDens=1;
        initialDimRes1=10;
        initialDimRes2=1;
        dimRes1=10;
        dimRes2=1;
        fixe=true;
        vide=1;
        zone=false;
        DensSpatiallyHeterog=false;
        vide_zone=-666;
        dens_zone=-666;
        xmin_zone=-666;
        xmax_zone=-666;
        ymin_zone=-666;
        ymax_zone=-666;
		backwardBarrier=false;
		forwardBarrier=false;
		backwardBarrierX=false;
		backwardBarrierY=false;
        forwardBarrierX=false;
        forwardBarrierY=false;
        x1_barrier=1;
        x2_barrier=1;
        y1_barrier=1;
        y2_barrier=1;
        barrierCrossingRateUp=1.0;
        barrierCrossingRateDown=1.0;
        dist_max=-1;
        mod='b';
        Mig=0.5;
        GeoG=0.5;
        Pareto_Shape=5;
        Sichel_Gamma=-2.15;
        Sichel_Xi=100;
        Sichel_Omega=-1;
        DispSpatiallyHeterog=false;
        dist_max_zone=-666;
        mod_zone='n';
        Mig_zone=-666;
        GeoG_zone=-666;
        Pareto_Shape_zone=-666;
        Sichel_Gamma_zone=-666;
        Sichel_Xi_zone=-666;
        Sichel_Omega_zone=-666;
        ContDemeSizeChange="None";
        ContLatticeSizeChange="None";		
        densLogisticGrowthRate=0.0;
        latticeLogisticGrowthRate=0.0;
        minGen=1;
    };
    CTimeVaryingParams(CTimeVaryingParams* modele)  {
		initialDens=modele->initialDens;
		initialDimRes1=modele->initialDimRes1;
		initialDimRes2=modele->initialDimRes2;
		dens=modele->dens;
        dimRes1=modele->dimRes1;
        dimRes2=modele->dimRes2;
        fixe=modele->fixe;
        vide=modele->vide;
        zone=modele->zone;
        DensSpatiallyHeterog=modele->DensSpatiallyHeterog;
        dens_zone=modele->dens_zone;
        vide_zone=modele->vide_zone;
		xmin_zone=modele->xmin_zone;
        xmax_zone=modele->xmax_zone;
        ymin_zone=modele->ymin_zone;
        ymax_zone=modele->ymax_zone;
        backwardBarrier=modele->backwardBarrier;
        forwardBarrier=modele->forwardBarrier;
        backwardBarrierX=modele->backwardBarrierX;
        backwardBarrierY=modele->backwardBarrierY;
        forwardBarrierX=modele->forwardBarrierX;
        forwardBarrierY=modele->forwardBarrierY;
        x1_barrier=modele->x1_barrier;
        x2_barrier=modele->x2_barrier;
        y1_barrier=modele->y1_barrier;
        y2_barrier=modele->y2_barrier;
        barrierCrossingRateUp=modele->barrierCrossingRateUp;
        barrierCrossingRateDown=modele->barrierCrossingRateDown;
        dist_max=modele->dist_max;
        mod=modele->mod;
        Mig=modele->Mig;
        GeoG=modele->GeoG;
        Pareto_Shape=modele->Pareto_Shape;
        Sichel_Gamma=modele->Sichel_Gamma;
        Sichel_Xi=modele->Sichel_Xi;
        Sichel_Omega=modele->Sichel_Omega;
        DispSpatiallyHeterog=modele->DispSpatiallyHeterog;
        dist_max_zone=modele->dist_max_zone;
        mod_zone=modele->mod_zone;
        Mig_zone=modele->Mig_zone;
        GeoG_zone=modele->GeoG_zone;
        Pareto_Shape_zone=modele->Pareto_Shape_zone;
        Sichel_Gamma_zone=modele->Sichel_Gamma_zone;
        Sichel_Xi_zone=modele->Sichel_Xi_zone;
        Sichel_Omega_zone=modele->Sichel_Omega_zone;
        ContDemeSizeChange="None";//il ne faut pas passer cet argument sinon on refait un changement de taille de pop continu, a revoir RL->RL
		ContLatticeSizeChange="None";//il ne faut pas passer cet argument sinon on refait un changement de taille de pop continu, a revoir RL->RL
		densLogisticGrowthRate=0.0;
		latticeLogisticGrowthRate=0.0;
        minGen=1;
    };
    ~CTimeVaryingParams() {};
    int setGrowthPars(CTimeVaryingParams* nextTVpars);
    int computeCumuls();
    int currentDensAndLatticeSize(long unsigned int& locCurrentGeneration);
    void SetForwardDispersalDistributions();
    void SetForwardDispersalDistributions_zone();
    void oneDProductWithoutm0();
    void realise_disp(std::vector<int> n_genes_loc);
    int calcul_densite(int coordx,int coordy);
	int PrintTVpars();

    private:

    long int Dt;//difference de temps entre phases demo
    int Ddens,DdimRes1,DdimRes2;//differences de taille de demes ou de tailles de lattice entre phases demo
	float Rdens,RdimRes1,RdimRes2;//ratio de taille de pops entre phases demo
    float growthRate,growthRate1,growthRate2;
    int KLogistic,KLogistic1,KLogistic2;
};

extern std::vector<CTimeVaryingParams> TVpars;


#endif
