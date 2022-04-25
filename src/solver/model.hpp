#ifndef __model__hpp
#define __model__hpp


/****************************************************************************************/
/*										LIBRARIES										*/
/****************************************************************************************/

/*** Own Libraries ***/
#include "callback.hpp"

#include <limits>
/****************************************************************************************/
/*										TYPEDEFS										*/
/****************************************************************************************/

/*** CPLEX ***/	
typedef std::vector<IloNumVar>          IloNumVarVector;
typedef std::vector<IloNumVarVector>    IloNumVarMatrix;
typedef std::vector<IloNumVarMatrix>    IloNumVar3DMatrix;
typedef std::vector<IloNumVar3DMatrix>  IloNumVar4DMatrix;
typedef std::vector<IloNumVar4DMatrix>  IloNumVar5DMatrix;

typedef std::vector<IloNum>            IloNumVector;
typedef std::vector<IloNumVector>      IloNumMatrix;
typedef std::vector<IloNumMatrix>      IloNum3DMatrix;
typedef std::vector<IloNum3DMatrix>    IloNum4DMatrix;
typedef std::vector<IloNum4DMatrix>    IloNum5DMatrix;


/********************************************************************************************
 * This class models the MIP formulation and solves it using CPLEX. 											
********************************************************************************************/

class Model
{
	private:
		/*** General features ***/
		const IloEnv&   env;    /**< IBM environment **/
	 	IloModel        model;  /**< IBM Model **/
		IloCplex        cplex;  /**< IBM Cplex **/
		const Data&     data;   /**< Data read from parameters files **/

		/*** Formulation specific ***/
		
		// Variables required for modelling the Resilient VNF placement problem
		IloNumVarMatrix 	y;              /**< VNF placement variables. y[v][f] **/
		IloNumVar3DMatrix 	x;            	/**< VNF assignement variables. x[k][i][v] **/
		
		// Variables required for including routing constraints
		IloNumVar4DMatrix 	z;            	/**< VNF pair assignement variables. z[k][i][s][t] **/
		IloNumVar5DMatrix 	r;            	/**< Routing variables. r[k][i][a][s][t] **/
		IloNumVarMatrix 	delay;          /**< Section delay variables. delay[k][i] **/
		IloNumVar3DMatrix 	arc_usage;      /**< Section arc usage variables. arc_usage[k][i][a] **/

		// Approximation related Variables
		IloNumVarMatrix secAvail;			/**< Real variable between 0 and 1 representing the availability of a section**/
		IloNumVarMatrix secUnavail;			/**< Real variable between 0 and 1 representing the unavailability of a section**/
		
		/*** Formulation general ***/
		IloObjective    	obj;            /**< Objective function **/
		IloRangeArray   	constraints;    /**< Set of constraints **/

		/*** Manage execution and control ***/
		Callback* 			callback; 		/**< User generic callback **/
		IloNum time;						/**< Time spent during the optimization **/

	public:
	/****************************************************************************************/
	/*										Constructors									*/
	/****************************************************************************************/
		/** Constructor. Builds the model (variables, objective function, constraints and further parameters). **/
		Model(const IloEnv& env, const Data& data);
		Model(const IloEnv& env, const Data&&) = delete;
		Model() = delete;

	/****************************************************************************************/
	/*									    Formulation  									*/
	/****************************************************************************************/
        /** Set up the Cplex parameters. **/
        void setCplexParameters();

        /** Set up the variables. **/
        void setVariables();
		/** Set up placement variables. **/
		void setPlacementVariables();
		/** Set up assignment variables. **/
		void setAssignmentVariables();
		/** Set up VNF pair assignement variables. **/
		void setPairAssignmentVariables();
		/** Set up routing variables. **/
		void setRoutingVariables();
		/** Set up delay variables. **/
		void setDelayVariables();
		/** Set up arc usage variables. **/
		void setArcUsageVariables();
		/** Set up availability variables. **/
		void setAvailabilityVariables();

        /** Set up the objective function. **/
        void setObjective();


        /** Set up the constraints. **/
        void setConstraints();
        /** Add up the VNF assignment constraints: At least one VNF must be assigned to each section of each demand. **/
        void setVnfAssignmentConstraints();
        /** Add up the VNF placement constraints: a VNF can only be assigned to a demand if it is already placed. **/
        void setVnfPlacementConstraints();
        /** Add up the original aggregated VNF placement constraints. **/
        void setOriginalVnfPlacementConstraints();
        /** Add up the node capacity constraints: the bandwidth treated in a node must respect its capacity. **/
        void setNodeCapacityConstraints();
        /** Add up the strong node capacity constraints. **/
        void setStrongNodeCapacityConstraints();
        /** Add up the delay constraints: The longest SFC path must not exceed its required latency. **/
        void setDelayConstraints();
        /** Add up the bandwidth constraints: The bandwidth allocated on each arc must be at most its capacity. **/
        void setBandwidthConstraints();
        /** Add up the linking constraints: Constraints linking variables. **/
        void setLinkingConstraints();
        /** Add up the routing constraints: There must be a path between any two consecutive VNFs. **/
        void setRoutingConstraints();

		void setSectionAvailabilityApproxConstraints();
		void setSFCAvailabilityApproxConstraints();

	/****************************************************************************************/
	/*										   Methods  									*/
	/****************************************************************************************/
		/** Solves the MIP. **/
		void run();

		/** Approximation related methods **/
		void buildApproximationFunctionAvail	(int demand, IloNumArray &breakpoints, IloNumArray &slopes);
		void buildApproximationFunctionUnavail	(int demand, IloNumArray &breakpoints, IloNumArray &slopes);
		void buildAvailTouchs					(int demand, IloNumArray &touchs);
		void buildAvailBreakpoints				(int demand, IloNumArray &breakpoints, IloNumArray &touchs);
		void buildUnavailTouchs					(int demand, IloNumArray &touchs);
		void buildUnavailBreakpoints			(int demand, IloNumArray &breakpoints, IloNumArray &touchs);

	/****************************************************************************************/
	/*									Solution Query  									*/
	/****************************************************************************************/
		/** Displays the obtained results **/
		void printResult				();
		void printDemand				(const int demand);
		void printSectionPlacement		(const int demand, const int section);
		void printSectionAvailability	(const int demand, const int section);
		void printRouting				(const int demand, const int section);
		void printPath					(const int demand, const int section, Graph::Node &s, Graph::Node &t);
		
		/** Auxialiary getters **/
		double 	getServiceAvail		(const int demand);
		int 	getNbAvailViolation	();
		double 	getMaxAvailViolation();

		/** Outputs the obtained results **/
		void output();

	/****************************************************************************************/
	/*										Destructors 									*/
	/****************************************************************************************/
		/** Destructor. Free dynamic allocated memory. **/
		~Model();
};


#endif // MODELS_HPP