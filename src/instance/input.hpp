#ifndef __input__hpp
#define __input__hpp

#include <iostream>
#include <fstream>
#include <string>
#include <vector>


/*****************************************************************************************
 * This class stores all the information recovered from the parameter file, that is,
 * input/output file paths, execution and control parameters.						
*****************************************************************************************/
class Input{

public: 
	/** States wheter disaggregated VNF placement constraints are activated.**/
	enum Disaggregated_VNF_Placement_Constraints {
		DISAGGREGATED_VNF_PLACEMENT_OFF = 0,  		
		DISAGGREGATED_VNF_PLACEMENT_ON  = 1 	
	};
	/** States wheter strong node capacity constraints are activated.**/
	enum Strong_Node_Capacity_Constraints {
		STRONG_NODE_CAPACITY_OFF = 0,  		
		STRONG_NODE_CAPACITY_ON  = 1 	        
	};
	/** States wheter node cover cuts are activated.**/
	enum Availability_Usercuts {
		AVAILABILITY_USERCUTS_OFF = 0,  		
		AVAILABILITY_USERCUTS_ON  = 1 	        
	};
	/** States wheter node cover cuts are activated.**/
	enum Node_Cover_Cuts {
		NODE_COVER_OFF = 0,  		
		NODE_COVER_ON  = 1 	        
	};
	/** States wheter chain cover cuts are activated.**/
	enum Chain_Cover_Cuts {
		CHAIN_COVER_OFF = 0,  		
		CHAIN_COVER_ON  = 1 	        
	};
    
	/** States wheter vnf lower bound cuts are activated.**/
	enum VNF_Lower_Bound_Cuts {
		VNF_LOWER_BOUND_CUTS_OFF = 0,  		
		VNF_LOWER_BOUND_CUTS_ON  = 1 	        
	};

	/** States wheter section failure cuts are activated.**/
	enum Section_Failure_Cuts {
		SECTION_FAILURE_CUTS_OFF = 0,  		
		SECTION_FAILURE_CUTS_ON  = 1 	        
	};
	/** States wheter routing is activated.**/
	enum Routing {
		ROUTING_OFF = 0,  		
		ROUTING_ON  = 1 	        
	};
	enum Approximation_Type {
		APPROXIMATION_TYPE_RESTRICTION  = -1,
		APPROXIMATION_TYPE_NONE			= 0,  		
		APPROXIMATION_TYPE_RELAXATION   = 1 	        
	};
	/** States wheter lazy constraints are activated.**/
	enum Lazy_Constraints {
		LAZY_OFF = 0,  		
		LAZY_ON  = 1 	        
	};
	/** States wheter heuristics are activated.**/
	enum Heuristic {
		HEURISTIC_OFF = 0,  		
		HEURISTIC_ON  = 1 	        
	};

private:
    /***** Input file paths *****/
    const std::string   parameters_file;
    std::string         node_file;
    std::string         link_file;
    std::string         demand_file;
    std::string         vnf_file;

    /***** Formulation parameters*****/
    Disaggregated_VNF_Placement_Constraints disaggregated_vnf_placement;    /**< Refers to the activation of disaggregated VNF placement constraints. **/
    Strong_Node_Capacity_Constraints        strong_node_capacity;           /**< Refers to the activation of strong node capacity constraints. **/
    Availability_Usercuts                   availability_cuts;              /**< Refers to the activation of availability cuts. **/
    Node_Cover_Cuts                         node_cover;                     /**< Refers to the activation of node cover cuts. **/
    Chain_Cover_Cuts                        chain_cover;                    /**< Refers to the activation of chain cover cuts. **/
    VNF_Lower_Bound_Cuts                    vnf_lower_bound;                /**< Refers to the activation of vnf lower bound cuts. **/
    Section_Failure_Cuts					section_failure_cuts; 			/**< Refers to the activation of section failure cuts. **/
    Routing									routing_activation; 			/**< Refers to the activation of routing decisions. **/
	Approximation_Type                      approx_type;                    /**< Refers to the type of approximation used for modeling availability constraints. **/
	Lazy_Constraints 						lazy; 							/**< Refers to the activation of lazy constraints. **/
	Heuristic								heuristic_activation;			/**< Refers to the activation of heuristics. **/

	/***** Optimization parameters*****/
    bool                linear_relaxation;
    int                 time_limit;
    int                 nb_breakpoints;


    /***** Output file paths *****/
    std::string         output_file;
    
public:
	/****************************************************************************************/
	/*				    				Constructors    									*/
	/****************************************************************************************/
	/** Constructor initializes the object with the information contained in the parameter file. @param file The address of the parameter file (usually the address of file 'parameters.txt'). **/
    Input(const std::string file);
    /** Constructor always need a parameter file. **/
    Input() = delete;



	/****************************************************************************************/
	/*				    					Getters	    									*/
	/****************************************************************************************/
    /** Returns the parameters file. */
    const std::string& getParameterFile()  const { return this->parameters_file; }

    /** Returns the node file. */
    const std::string& getNodeFile()       const { return this->node_file; }

    /** Returns the link file. */
    const std::string& getLinkFile()       const { return this->link_file; }

    /** Returns the demand file. */
    const std::string& getDemandFile()     const { return this->demand_file; }

    /** Returns the VNF file. */
    const std::string& getVnfFile()        const { return this->vnf_file; }

	/** Returns whether disaggregated vnf placement constraints are activated. **/
    const Disaggregated_VNF_Placement_Constraints & getDisaggregatedVnfPlacement()  const { return disaggregated_vnf_placement; }
	/** Returns whether strong node capacity constraints are activated. **/
    const Strong_Node_Capacity_Constraints &        getStrongNodeCapacity()         const { return strong_node_capacity; }
	/** Returns whether node cover cuts are activated. **/
    const Node_Cover_Cuts &                         getNodeCover()                  const { return node_cover; }
	/** Returns whether chain cover cuts are activated. **/
    const Chain_Cover_Cuts &                        getChainCover()                 const { return chain_cover; }
	/** Returns whether availability cuts are activated. **/
    const Availability_Usercuts &                   getAvailabilityUsercuts()       const { return availability_cuts; }
	/** Returns whether vnf lower bound cuts are activated. **/
    const VNF_Lower_Bound_Cuts &                    getVnfLowerBoundCuts()          const { return vnf_lower_bound; }
	/** Returns whether section failure cuts are activated. **/
    const Section_Failure_Cuts &                    getSectionFailureCuts()         const { return section_failure_cuts; }
	/** Returns whether routing is activated. **/
    const Routing &                    				getRoutingActivation()         	const { return routing_activation; }
	/** Returns the type of availability approximation to be used. **/ 
    const Approximation_Type &                      getApproximationType()         	const { return approx_type; }
	/** Returns whether lazy constraints are activated **/ 
    const Lazy_Constraints &                      	getLazy()         				const { return lazy; }
	/** Returns whether heuristics are used. **/ 
    const Heuristic &                      			getHeuristic()         			const { return heuristic_activation; }
	

    /** Returns true if linear relaxation is to be applied. */
    const bool&        isRelaxation()      const { return this->linear_relaxation; }

    /** Returns time limit in seconds to be applied. */
    const int&         getTimeLimit()      const { return this->time_limit; }

    /** Returns the number of breakpoints to be used in the log approximation. */
    const int&         getNbBreakpoints() const { return this->nb_breakpoints; }

    /** Returns the output file. */
    const std::string& getOutputFile()     const { return this->output_file; }

	/****************************************************************************************/
	/*				    					Methods	    									*/
	/****************************************************************************************/
    /** Returns the pattern value in the parameters file. */
    std::string getParameterValue(const std::string pattern);

	/****************************************************************************************/
	/*				    					Display	    									*/
	/****************************************************************************************/
    /** Print the parameters stored in the parameter file. */
    void print();

};

#endif