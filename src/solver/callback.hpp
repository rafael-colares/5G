#ifndef __callback__hpp
#define __callback__hpp

/****************************************************************************************/
/*										LIBRARIES										*/
/****************************************************************************************/

/*** C++ Libraries ***/
#include <thread>
#include <mutex>

/*** CPLEX Libraries ***/
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

/*** Own Libraries ***/
#include "../instance/data.hpp"
#include "../tools/others.hpp"

/****************************************************************************************/
/*										TYPEDEFS										*/
/****************************************************************************************/

/*** Modelling matrices of variables and coefficients with CPLEX types ***/	
typedef std::vector<IloNumVar>          IloNumVarVector;
typedef std::vector<IloNumVarVector>    IloNumVarMatrix;
typedef std::vector<IloNumVarMatrix>    IloNumVar3DMatrix;
typedef std::vector<IloNumVar3DMatrix>  IloNumVar4DMatrix;
typedef std::vector<IloNumVar4DMatrix>  IloNumVar5DMatrix;

typedef std::vector<IloNum>             IloNumVector;
typedef std::vector<IloNumVector>       IloNumMatrix;
typedef std::vector<IloNumMatrix>       IloNum3DMatrix;
typedef std::vector<IloNum3DMatrix>     IloNum4DMatrix;
typedef std::vector<IloNum4DMatrix>     IloNum5DMatrix;

typedef IloCplex::Callback::Context     Context;

/****************************************************************************************/
/*										DEFINES			    							*/
/****************************************************************************************/
#define EPS 1e-4        // Large tolerance, used for float precision
#define EPSILON 1e-6    // Small tolerance, used for float precision

/************************************************************************************
 * This class implements the generic callback interface. It has two main 
 * functions: addUserCuts and addLazyConstraints.
 ************************************************************************************/
class Callback: public IloCplex::Callback::Function {

    
private:
    /*** General variables ***/
    const IloEnv&   env;    /**< IBM environment **/
    const Data&     data;   /**< Data read in data.hpp **/


    /*** LP data ***/
	const IloNumVar3DMatrix&    x;                  /**< VNF assignement variables **/
	const IloNumVarMatrix&      y;                  /**< VNF assignement variables **/
    const IloNumVarMatrix&      secAvail;			/**< Real variable between 0 and 1 representing the availability of a section**/
    const IloNumVarMatrix&      secUnavail;			/**< Real variable between 0 and 1 representing the unavailability of a section**/
    	
    IloRangeArray cutPool;                          /**< Cutpool to be checked on each node. **/

    
    /*** Solution data ***/
    IloNumMatrix        ySol;               /**< Stores the y variables from a given solution **/
    IloNum3DMatrix      xSol;               /**< Stores the x variables from a given solution **/
    double              objSol;             /**< Stores the objective function value from a given solution **/
    std::vector<double> remainingCapacity;  /**< Stores the remaining capacity of each node in the graph **/

    /*** Manage execution and control ***/
    std::mutex  thread_flag;                /**< A mutex for synchronizing multi-thread operations. **/
    int         nb_cuts_avail_heuristic;    /**< Number of availability cuts added through heuristic procedure. **/
    int         nbLazyConstraints;          /**< Number of lazy constraints added. **/
    int         nbCuts;                     /**< Total number of user cuts added. **/
    IloNum      timeAll;                    /**< Total time spent on callback. **/


public:

	/****************************************************************************************/
	/*										Constructors									*/
	/****************************************************************************************/
    /** Constructor. Initializes callback variables. **/
	Callback(const IloEnv& env_, const Data& data_, 
                const IloNumVar3DMatrix& x_, const IloNumVarMatrix& y_,
                const IloNumVarMatrix& secAvail_, const IloNumVarMatrix& secUnavail_);


    /****************************************************************************************/
	/*									Auxliary Structs     								*/
	/****************************************************************************************/
    /** Stores the section id and its availability. Used for the separation of integer solutions. **/
    struct MapAvailability { 
        int section;            /**< The section id.*/
        double availability;    /**< The section availability. */
    }; 


	/****************************************************************************************/
	/*									Main operations  									*/
	/****************************************************************************************/
    /** CPLEX will call this method during the solution process at the places that we asked for. @param context Defines on which places the method is called and we use it do define what to do in each case.**/
    void    invoke                  (const Context& context);

    /** Solves the separation problems for a given fractional solution. @note Should only be called within relaxation context.**/
	void    addUserCuts             (const Context& context); 
    
    /** Solves the separation problems for a given integer solution. @note Should only be called within candidate context.**/
    void    addLazyConstraints      (const Context& context);

    /** Launches the matheuristic procedure based on a given fractional solution. @note Should only be called within relaxation context.**/
    void    runHeuristic            (const Context& context);
    
    /** Returns the current integer solution. @note Should only be called within candidate context. **/ 
    void    getIntegerSolution      (const Context &context);
    
    /** Returns the current fractional solution. @note Should only be called within relaxation context. **/ 
    void    getFractionalSolution   (const Context &context);
    
    /** Checks whether the current solution satisfies all cuts in the pool and add the unsatisfied one. **/
    bool    checkCutPool            (const Context &context);

	/****************************************************************************************/
	/*							Heuristic Related Methods  				    			    */
	/****************************************************************************************/
    /** Launches the phase I of the matheuristic procedure.**/
    void    runHeuristic_Phase_I    (const Context& context);

    /** Launches the phase II of the matheuristic procedure. Returns true if a feasible solution was found. **/
    bool    runHeuristic_Phase_II   (const Context& context);

    /** Checks if the heuristic should be launched. **/
    bool    heuristicRule           (const Context &context);

    /** Posts an heuristic solution into the optimization procedure. **/
    void    insertHeuristicSolution (const Context &context);

    /** Returns the availability of SFC k obtained from the solution stored in xSol. @note least will store the index of the least available section of the SFC **/
    double  getSolutionAvail_k      (int k, int& least);

    /** Chooses on which node VNF f should be installed for demand k **/
    int     getNodeToInstall        (int f, int k);

	/****************************************************************************************/
	/*							Cut Pool Definition Methods  							    */
	/****************************************************************************************/
    /** Sets up the cut pool that is checked on relaxation context. On this pool, only cuts appearing in a polynomial number are added. **/
    void setCutPool                      ();

    /** Add (some) availability cover constraints to the cut pool. Only sets {j \in V : a(j) <= a(v) } for each v \in V are considered. **/
    void addAvailabilityCoverConstraints ();
    
    /** Add vnf lower bound constraints to the cut pool. @note These are chain cover constraints where the set of considered sections is the whole set of sections, i.e., Q = I. **/
    void addVnfLowerBoundConstraints ();

    /** Add section failure constraints to the cut pool. **/
    void addSectionFailureConstraints ();

	/****************************************************************************************/
	/*							Availability Separation Methods  							*/
	/****************************************************************************************/
    /** Greedly solves the separation problem associated with the availability constraints. **/
    void heuristicSeparationOfAvailibilityConstraints(const Context &context, const IloNum3DMatrix& xSol);

    /** Initializes the availability heuristic. **/
    void initiateHeuristic(const int k, std::vector< std::vector<int> >& coeff, std::vector< std::vector<int> >& sectionNodes, std::vector< double >& sectionAvailability, const IloNum3DMatrix& xSol);
	
    /** Computes the availability increment resulted from the instalation of a new vnf. @param CHAIN_AVAIL The chain required availability. @param deltaAvail The matrix to be computed. @param sectionAvail THe current section availabilities. @param coeff The matrix of coefficients storing the possible vnfs to be placed. **/
    void computeDeltaAvailability(const double CHAIN_AVAIL, std::vector< std::vector<double> >& deltaAvail, const std::vector< double >& sectionAvail, const std::vector< std::vector<int> >& coeff);
    
    /** Tries to add new vnf placements to the current solution without changing its availability violation. @param xSol The current solution for a given demand. @param availabilityRequired The SFC required availability. @param sectionAvailability The current section availabilities. @param nbSections The number of sections that can be modified. **/
    void lift(IloNumMatrix& xSol, const double& availabilityRequired, std::vector<MapAvailability>& sectionAvailability, const int& nbSections);

	/****************************************************************************************/
	/*							    Cover Separation Methods    							*/
	/****************************************************************************************/
    /** Solves the separation problem associated with the chain cover constraints. **/
    void chainCoverSeparation(const Context &context, const IloNum3DMatrix& xSol);
    
    /** Solves the separation problem associated with the generalized cover constraints. **/
    void generalizedCoverSeparation(const Context &context, const IloNum3DMatrix& xSol);

    /****************************************************************************************/
	/*							    Integer solution query methods 							*/
	/****************************************************************************************/
    /** Returns the availability of the i-th section of a SFC demand obtained from an integer solution. @param k The demand id. @param i The section id. @param xSol The current integer solution. **/
    double getAvailabilityOfSection (const int& k, const int& i, const IloNum3DMatrix& xSol) const;
    
    /** Returns the availabilities of the sections of a SFC demand obtained from an integer solution. @param k The demand id. @param xSol The current integer solution. **/
    std::vector<MapAvailability> getAvailabilitiesOfSections (const int& k, const IloNum3DMatrix& xSol) const;
    

	/****************************************************************************************/
	/*								     Other Query Methods	    	    	    		*/
	/****************************************************************************************/
    /** Returns the number of user cuts added so far. **/ 
    const int    getNbUserCuts()           const{ return nbCuts; }

    /** Returns the number of lazy constraints added so far. **/ 
    const int    getNbLazyConstraints()    const{ return nbLazyConstraints; }

    /** Returns the total time spent on callback so far. **/ 
    const IloNum getTime()                 const{ return timeAll; }

    /** Checks if all placement variables of a given SFC demand are integers. @param k The demand id. @param xSol The current solution. **/
    const bool   isIntegerAssignment (const int& k, const IloNum3DMatrix& xSol) const;
    

	/****************************************************************************************/
	/*								Thread Protected Methods			    				*/
	/****************************************************************************************/
    /** Increase by one the number of lazy constraints added. **/
    void incrementLazyConstraints();
    /** Increase by one the number of availability cuts added through the heuristic procedure. **/
    void incrementAvailabilityCutsHeuristic();
    /** Increase by one the number of user cuts added. **/
    void incrementUsercuts();
    /** Increases the total callback time. @param time The time to be added. **/
    void incrementTime(const IloNum time);

	/****************************************************************************************/
	/*										Destructors			    						*/
	/****************************************************************************************/
    /** Destructor **/
    ~Callback() {}

};

/** Checks if the availability of a is lower than the one of b. **/
bool compareAvailability(Callback::MapAvailability a, Callback::MapAvailability b);

#endif