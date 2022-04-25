#include "callback.hpp"

/****************************************************************************************/
/*										CONSTRUCTOR										*/
/****************************************************************************************/
/** Callback constructor. This is called only once, before the optimization procedure is launched. **/
Callback::Callback(const IloEnv& env_, const Data& data_, 
                    const IloNumVar3DMatrix& x_, const IloNumVarMatrix& y_,
                    const IloNumVarMatrix& secAvail_, const IloNumVarMatrix& secUnavail_) :
                    env(env_), data(data_),	
                    x(x_), y(y_), 
                    secAvail(secAvail_), secUnavail(secUnavail_),
                    cutPool(env)
{	
	/*** Control ***/
    thread_flag.lock();

	nb_cuts_avail_heuristic = 0;
    nbLazyConstraints       = 0;
    nbCuts                  = 0;
	timeAll                 = 0;
    setCutPool();

	// Solution related initialization
    const int NB_NODES = lemon::countNodes(data.getGraph());
    ySol.resize(lemon::countNodes(data.getGraph()));
    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
        int v = data.getNodeId(n);
        ySol[v].resize(data.getNbVnfs());
    }
    xSol.resize(data.getNbDemands());
    for (int k = 0; k < data.getNbDemands(); k++){
        xSol[k].resize(data.getDemand(k).getNbVNFs());
        for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
            xSol[k][i].resize(NB_NODES);
        }
    }

    // heuristic related initializations
    objSol = 0.0;
    remainingCapacity.resize(NB_NODES);

    // input a big hard-coded number to serve as seed
    const int SEED = 20102019;
    srand(SEED);

    thread_flag.unlock();
}

/****************************************************************************************/
/*										Operations  									*/
/****************************************************************************************/
/** This is the main callback function and determines what to do when the callback is invoked during the optimization. **/
void Callback::invoke(const Context& context)
{
    IloNum time_start = context.getDoubleInfo(IloCplex::Callback::Context::Info::Time);
    switch (context.getId()){
        /* When fractional solution is available */
        case Context::Id::Relaxation:
        {
            int current_number_of_cuts = getNbUserCuts();
            // look up for user cuts and add them
            addUserCuts(context);
            // if no additional cut was found, launch heuristic
            if (current_number_of_cuts == getNbUserCuts())  runHeuristic(context);
            break;
        }
        /* When integer solution is available */
        case Context::Id::Candidate:
        {
            // if the candidate solution is considered feasible, check if all lazy constraints are satisfied
			if (context.isCandidatePoint()) {
	    		addLazyConstraints(context);
			}
            break;
        }
        /* Not an option for callback */
		default:
			throw IloCplex::Exception(-1, "ERROR: Unexpected context id !");
    }
    IloNum time_spent = context.getDoubleInfo(IloCplex::Callback::Context::Info::Time) - time_start;
    incrementTime(time_spent);
}

/** Launches the matheuristic procedure based on a given fractional solution. @note Should only be called within relaxation context.**/
void Callback::runHeuristic(const Context &context)
{
    if ((data.getInput().getHeuristic() == Input::HEURISTIC_OFF) || (heuristicRule(context) == false)) return;
    
    try{
        runHeuristic_Phase_I(context);
        bool isFeasible = runHeuristic_Phase_II(context);
        if ((isFeasible) && (objSol < context.getIncumbentObjective())){
            insertHeuristicSolution(context);
        }
    }
    catch (...) {
        throw;
    }
}

/** Builds (a possibly unfeasible) integer solution **/
void Callback::runHeuristic_Phase_I(const Context &context){
    //build placement y
    objSol = 0.0;
    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
        int v = data.getNodeId(n);
        for (int f = 0; f < data.getNbVnfs(); f++){
            const double RND = ((double) rand() / (RAND_MAX));
            if (RND <= context.getRelaxationPoint(y[v][f])){
                ySol[v][f] = 1;
                objSol += data.getPlacementCost(data.getNode(v), data.getVnf(f));
            }
            else{
                ySol[v][f] = 0;
            }
        }
    }
    //build assignment x
    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
        int v = data.getNodeId(n);
        remainingCapacity[v] = data.getNode(v).getCapacity();
        for (int k = 0; k < data.getNbDemands(); k++){
            for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
                int f = data.getDemand(k).getVNF_i(i);
                double req_capacity = data.getDemand(k).getBandwidth() * data.getVnf(f).getConsumption();
                if (ySol[v][f] == 1 && req_capacity <= remainingCapacity[v]){
                    // there is a chance of assigning the vnf
                    const double RND = ((double) rand() / (RAND_MAX));
                    if (RND <= context.getRelaxationPoint(x[k][i][v])){
                        xSol[k][i][v] = 1;
                        remainingCapacity[v] -= req_capacity;
                    }
                    else{
                        xSol[k][i][v] = 0;
                    }
                }
                else{
                    xSol[k][i][v] = 0;
                }
            }
        }
    }
}

/** Launches the phase II of the matheuristic procedure. Returns true if a feasible solution was found. **/
bool Callback::runHeuristic_Phase_II(const Context &context){
    for (int k = 0; k < data.getNbDemands(); k++){
        int i = 0;
        const double REQ_AVAIL = data.getDemand(k).getAvailability();
        while (getSolutionAvail_k(k, i) < REQ_AVAIL){
            int f = data.getDemand(k).getVNF_i(i);
            //choose node to install the ith vnf of sfc k
            int v = getNodeToInstall(f, k);
            if (v == -1) return false;
            
            // set x[k,i,v] to 1 and y[v, f(i,k)] also if needed
            xSol[k][i][v] = 1;
            remainingCapacity[v] -= data.getDemand(k).getBandwidth() * data.getVnf(f).getConsumption();
            if (ySol[v][f] == 0){
                ySol[v][f] = 1;
                objSol += data.getPlacementCost(data.getNode(v), data.getVnf(f));
            }
        }
    }
    return true;
}

/** Chooses on which node VNF f should be installed for demand k **/
int Callback::getNodeToInstall(int f, int k){
    const double REQ_CAPACITY   = data.getDemand(k).getBandwidth() * data.getVnf(f).getConsumption();
    IloNum maxRemainingCapacity = 0.0;
    IloNum minValue             = IloInfinity;
    int selectedNode            = -1;

    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
        int v = data.getNodeId(n);
        if (REQ_CAPACITY <= remainingCapacity[v]){
            const IloNum ADDITIONAL_COST = (data.getPlacementCost(data.getNode(v), data.getVnf(f))) * (1.0 - ySol[v][f]);
            if (ADDITIONAL_COST <= minValue + EPSILON){
                if (ADDITIONAL_COST <= minValue - EPSILON){
                    selectedNode         = v;
                    minValue             = ADDITIONAL_COST;
                    maxRemainingCapacity = remainingCapacity[v];
                }
                else{
                    //check the remaining capacity
                    if (remainingCapacity[v] >= maxRemainingCapacity + EPSILON){
                        selectedNode         = v;
                        minValue             = ADDITIONAL_COST;
                        maxRemainingCapacity = remainingCapacity[v];
                    }
                }
            }
        }
    }
    return selectedNode;
}

/** Returns the availability of SFC k obtained from the solution stored in xSol. @note least will store the index of the least available section of the SFC **/
double Callback::getSolutionAvail_k(int k, int &leastAvailableSection){
    double availability     = 1.0;
    double minSectionAvail  = 1.0;
    leastAvailableSection   = -1;
    for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
        double section_fail = 1.0;
        // compute section availability
        for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
            int v = data.getNodeId(n);
            if (xSol[k][i][v] == 1){
                section_fail *= (1.0 - data.getNode(v).getAvailability());
            }
        }
	    double section_avail = 1.0 - section_fail;
        // check if it is the least available section
        if (section_avail < minSectionAvail){
            minSectionAvail = section_avail;
            leastAvailableSection = i;
        }
        // compute global availability
        availability *= section_avail;
    }
    return availability;
}

/** Posts an heuristic solution into the optimization procedure. **/
void Callback::insertHeuristicSolution(const Context &context){
    IloNumArray vals(context.getEnv());
    try{  
      IloNumVarArray vars(context.getEnv());
    
        // build heuristic solution representation
        for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
            int v = data.getNodeId(n);
            for (int f = 0; f < data.getNbVnfs(); f++){
                vars.add(y[v][f]);
                vals.add(ySol[v][f]);
           }
        }
        /*... fill Arrays ... */
        for (int k = 0; k < data.getNbDemands(); k++){
            for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
                for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                    int v = data.getNodeId(n);
                    vars.add(x[k][i][v]);
                    vals.add(xSol[k][i][v]);
                }
            }
        }
        
        //post heuristic solution
        double objVal = 0.0;
        for(NodeIt n(data.getGraph()); n != lemon::INVALID; ++n) {
            int v = data.getNodeId(n);
            for (int i = 0; i < data.getNbVnfs(); i++){
                int f = data.getVnf(i).getId();
                double cost = data.getPlacementCost(data.getNode(v), data.getVnf(f));
                objVal += ( cost*ySol[v][f] ); 
            }
        }
        context.postHeuristicSolution(vars, vals, objVal, IloCplex::Callback::Context::SolutionStrategy::NoCheck);
        vals.end();
    }
    catch (...) {
        vals.end();
        throw;
    }
}

/** Checks if the heuristic should be launched. **/
bool Callback::heuristicRule(const Context &context){
    if (data.getInput().getApproximationType() != Input::APPROXIMATION_TYPE_NONE) return false;
    
    const double OBJ    = context.getRelaxationObjective();
    const double UB     = context.getIncumbentObjective();
    const double LIMIT  = (UB - OBJ) / UB;
    const double RND    = ((double) rand() / (RAND_MAX));

    return (RND <= LIMIT);
}

/** Sets up the cut pool that is checked on relaxation context. On this pool only cuts appearing in a polynomial number are added. **/
void Callback::setCutPool()
{
    std::cout << "\t > Setting up pool of cuts... " << std::endl;
    if (data.getInput().getNodeCover() == Input::NODE_COVER_ON){
        addAvailabilityCoverConstraints();
    }
    
    if (data.getInput().getVnfLowerBoundCuts() == Input::VNF_LOWER_BOUND_CUTS_ON){
        addVnfLowerBoundConstraints();
    }

    if (data.getInput().getSectionFailureCuts() == Input::SECTION_FAILURE_CUTS_ON){
        addSectionFailureConstraints();
    }
}


/** Add section failure constraints to the cut pool. **/
void Callback::addSectionFailureConstraints ()
{
    std::cout << "Adding section failure cuts to the pool..." << std::endl;
    /* For each SFC */
    for (int k = 0; k < data.getNbDemands(); k++){
        /* For each VNF section */
        for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
            const double REQ_AVAIL = data.getDemand(k).getAvailability();
            const double RHS       = -std::log(1.0 - REQ_AVAIL);
            IloExpr exp(env);
            /* For each node */
            for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                const int v             = data.getNodeId(n);
                const double NODE_AVAIL = data.getNode(v).getAvailability();
                double coeff            = -std::log(1.0 - NODE_AVAIL);
                exp += coeff*x[k][i][v];
            }
            std::string name = "Section_Fail(" + std::to_string(k) + "," + std::to_string(i) + ")";
            cutPool.add(IloRange(env, RHS, exp, IloInfinity, name.c_str()));
            exp.clear();
            exp.end();
        }
    }
}

/** Add vnf lower bound constraints to the cut pool. **/
void Callback::addVnfLowerBoundConstraints ()
{
    std::cout << "Adding chain cover cuts to the pool..." << std::endl;
    /* For each SFC */
    for (int k = 0; k < data.getNbDemands(); k++){
        //std::cout << "Demand " << k << ": ";
        const double REQUIRED_AVAIL = data.getDemand(k).getAvailability();
        const int    NB_SECTIONS    = data.getDemand(k).getNbVNFs();
        const int    RHS            = data.getVnfLB(REQUIRED_AVAIL, NB_SECTIONS);
        if (RHS > NB_SECTIONS*data.getVnfLB(REQUIRED_AVAIL, 1)){
            IloExpr exp(env);
            /* For each VNF section */
            for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
                /* For each node */
                for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                    int v = data.getNodeId(n);
                    exp += x[k][i][v];
                }
            }
            std::string name = "VNF_LowerBound(" + std::to_string(k) + ")";
            cutPool.add(IloRange(env, RHS, exp, IloInfinity, name.c_str()));
            exp.clear();
            exp.end();
        }
    }
}

/** Add availability cover constraints to the cut pool. **/
void Callback::addAvailabilityCoverConstraints(){
    std::cout << "Adding node cover cuts to the pool..." << std::endl;
    /* For each SFC */
    for (int k = 0; k < data.getNbDemands(); k++){
        /* Compute the coefficient c[v] for each node v. */
        std::vector<int> c;
        c.resize(data.getAvailNodeRank().size());
        for (unsigned int i = 0; i < data.getAvailNodeRank().size(); i++){
            int v = data.getAvailNodeRank()[i];
            c[v] = data.getMinNbNodes(data.getDemand(k).getAvailability(), data.getNode(v).getAvailability());
        }
        /* For each VNF section */
        for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
            /* For each node */
            for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                int v = data.getNodeId(n);
                int pos = data.getNodeRankPosition(v);
                /* Define only the non-dominated constraints. */
                if (pos > 0) {
                    if (c[v] > c[data.getAvailNodeRank()[pos-1]]) {
                        /* Define availability cover constraint for S = {j \in V : a(j) <= a(v) } */
                        IloExpr exp(env);
                        for (NodeIt it(data.getGraph()); it != lemon::INVALID; ++it){
                            int node_id = data.getNodeId(it);
                            int coeff = 0;
                            if (data.getNodeRankPosition(node_id) < data.getNodeRankPosition(v)){
                                coeff = std::max(c[v] - c[node_id] + 1, 1);
                            }
                            else{
                                coeff = 1;
                            }
                            exp += coeff*x[k][i][node_id];
                        }
                        std::string name = "NodeCover(" + std::to_string(k) + "," + std::to_string(i) + "," + std::to_string(v) + ")";
                        cutPool.add(IloRange(env, c[v], exp, IloInfinity, name.c_str()));
                        exp.clear();
                        exp.end();
                    }
                }
            }
        }
    }
}

/** Solves the separation problems for a given fractional solution. @note Should only be called within relaxation context.**/
void Callback::addUserCuts(const Context &context)
{
    try {    
        getFractionalSolution(context);
        /** If no cut in the cutpool is violated, `
         *  then, look for the violated cuts in the 
         *  exponential-sized families of valid inequalities **/
        if (checkCutPool(context) == false){
            if (data.getInput().getChainCover() == Input::CHAIN_COVER_ON){
                generalizedCoverSeparation(context, xSol);
                chainCoverSeparation(context, xSol);
            }
            if (data.getInput().getAvailabilityUsercuts() == Input::AVAILABILITY_USERCUTS_ON){
                heuristicSeparationOfAvailibilityConstraints(context, xSol);
            }
        }
    }
    catch (...) {
        throw;
    }
}


/* Checks whether the current solution satisfies all cuts in the pool and add the unsatisfied one. */
bool Callback::checkCutPool(const Context &context){
    bool found_violated_cut = false;
    for (IloInt i = 0; i < cutPool.getSize(); ++i) {
        const IloRange& cut = cutPool[i];
        const IloNum    LHS = context.getRelaxationValue(cut.getExpr());

        if ( LHS < cut.getLB() - EPS || LHS > cut.getUB() + EPS ) {
            std::cout << "Adding " << cut.getName() << std::endl;
            context.addUserCut(cut, IloCplex::UseCutForce, IloFalse);
            incrementUsercuts();
            found_violated_cut = true;
            /* Uncomment next line to add only one violated cut at a time. */
            // return true;
        }
    }
    return found_violated_cut;
}

// user cut related heuristic
void Callback::initiateHeuristic(const int k, std::vector< std::vector<int> >& coeff, std::vector< std::vector<int> >& sectionNodes, std::vector< double >& sectionAvailability, const IloNum3DMatrix& xSol)
{
    coeff.resize(data.getDemand(k).getNbVNFs());
    sectionNodes.resize(data.getDemand(k).getNbVNFs());
    sectionAvailability.resize(data.getDemand(k).getNbVNFs());

    for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
        coeff[i].resize(xSol[k][i].size());
        for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
            int v = data.getNodeId(n);
            coeff[i][v] = 1;
        }
    }

    /* Initialization of placement */
    for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
        /* Place every integer variable */
        for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
            int v = data.getNodeId(n);
            if (xSol[k][i][v] >= 1 - EPS){
                sectionNodes[i].push_back(v);
                coeff[i][v] = 0;
            }
        }
        /* If still empty, select some initial node based on the best x/availability ratio. */
        if (sectionNodes[i].empty()){
            int selectedNode = -1;
            double bestValue = -1.0;
            for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                int v = data.getNodeId(n);
                if ( (xSol[k][i][v] / data.getNode(v).getAvailability()) > bestValue){
                    bestValue = (xSol[k][i][v] / data.getNode(v).getAvailability());
                    selectedNode = v;
                }
            }
            sectionNodes[i].push_back(selectedNode);
            coeff[i][selectedNode] = 0;
        }
        /* Set initial section availability. */
        sectionAvailability[i] = 1.0 - data.getFailureProb(sectionNodes[i]);
    }
}


/* Solves the separation problem associated with the chain cover constraints. */
void Callback::chainCoverSeparation(const Context &context, const IloNum3DMatrix& xSol)
{
    /* Check VNF placement availability for each demand */
    for (int k = 0; k < data.getNbDemands(); k++){
        std::vector<double> sum_over_nodes(data.getDemand(k).getNbVNFs(), 0.0);
        for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
            sum_over_nodes[i] = 0.0;
            for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                int v = data.getNodeId(n);
                sum_over_nodes[i] += xSol[k][i][v];
            }
        }
        std::vector<int> sorted_sections = getSortedIndexes_Asc(sum_over_nodes);
        for (int nb_sections = 1; nb_sections <= data.getDemand(k).getNbVNFs(); nb_sections++){
            double lhs = 0.0;
            for (int i = 0; i < nb_sections; i++){
                lhs += sum_over_nodes[sorted_sections[i]];
            }
            double rhs = data.getVnfLB(data.getDemand(k).getAvailability(), nb_sections);
            /* If violated, build and add cut. */
            if (lhs < rhs - EPS){
                IloExpr expr(env);
                for (int i = 0; i < nb_sections; i++){
                    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                        int v = data.getNodeId(n);
                        int section = sorted_sections[i];
                        expr += x[k][section][v];
                    }
                }
                std::string name = "ChainCoverCut" + std::to_string(k) + "," + std::to_string(nb_sections) + ")";

                IloRange cut(env, rhs, expr, IloInfinity, name.c_str());
                std::cout << "Adding " << cut.getName() << std::endl;
                context.addUserCut(cut, IloCplex::UseCutForce, IloFalse);
                incrementUsercuts();
                expr.end();
                break;
            }
        }
    }
}

/* Solves the separation problem associated with the generalized cover constraints. */
void Callback::generalizedCoverSeparation(const Context &context, const IloNum3DMatrix& xSol)
{
    /* Check VNF placement availability for each demand */
    for (int k = 0; k < data.getNbDemands(); k++){
        for (NodeIt node(data.getGraph()); node != lemon::INVALID; ++node){
            int limit_node = data.getNodeId(node);
            /* Define set U */
            std::vector<bool> set_U(data.getNbNodes(), false);
            for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                int v = data.getNodeId(n);
                if (data.getNodeRankPosition(v) >= data.getNodeRankPosition(limit_node)){
                    set_U[v] = true;
                }
            }
            for (int nb_sections = 1; nb_sections <= data.getDemand(k).getNbVNFs(); nb_sections++){
                double rhs = data.getVnfLB(set_U, nb_sections, data.getDemand(k).getAvailability());
                if (rhs >= 0){
                    /* Define sum of nodes of U */
                    std::vector<double> sum_over_nodes(data.getDemand(k).getNbVNFs(), 0.0);
                    for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
                        for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                            int v = data.getNodeId(n);
                            if (set_U[v] == true){
                                sum_over_nodes[i] += xSol[k][i][v];
                            }
                            else{
                                sum_over_nodes[i] += (rhs*xSol[k][i][v]);
                            }
                        }
                    }
                    std::vector<int> sorted_sections = getSortedIndexes_Asc(sum_over_nodes);
                    /* Compute left-hand side value */
                    double lhs = 0.0;
                    for (int i = 0; i < nb_sections; i++){
                        lhs += sum_over_nodes[sorted_sections[i]];
                    }
                    /* If violated, build and add cut. */
                    if (lhs < rhs - EPS){
                        IloExpr expr(env);
                        for (int i = 0; i < nb_sections; i++){
                            int section = sorted_sections[i];
                            for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                                int v = data.getNodeId(n);
                                if (set_U[v] == true){
                                    expr += x[k][section][v];
                                }
                                else{
                                    expr += rhs*x[k][section][v];
                                }
                            }
                        }
                        std::string name = "GenCoverCut" + std::to_string(k) + "," + std::to_string(nb_sections) + "," + std::to_string(limit_node) + ")";
                        IloRange cut(env, rhs, expr, IloInfinity, name.c_str());
                        std::cout << "Adding " << cut.getName() << std::endl;
                        context.addUserCut(cut, IloCplex::UseCutForce, IloFalse);
                        incrementUsercuts();
                        expr.end();
                        return;
                    }
                }
            }
        }
    }
}

/* Greedly solves the separation problem associated with the availability constraints. */
void Callback::heuristicSeparationOfAvailibilityConstraints(const Context &context, const IloNum3DMatrix& xSol)
{
    /* Check VNF placement availability for each demand */
    for (int k = 0; k < data.getNbDemands(); k++){
        
        /* Declare auxiliary structures. */
        std::vector< std::vector<int> > coeff;          // the variable coefficient in the constraint
        std::vector< std::vector<int> > sectionNodes;   // the set of nodes placed in each section
        std::vector< double > sectionAvailability;      // the availability assoaciated with the placement
        
        initiateHeuristic(k, coeff, sectionNodes, sectionAvailability, xSol);

        double chainAvailability = data.getChainAvailability(sectionAvailability);
        const double REQUIRED_AVAIL = data.getDemand(k).getAvailability(); 
        
        if (chainAvailability < REQUIRED_AVAIL){
            std::vector< std::vector<double> > deltaAvailability;
            deltaAvailability.resize(data.getDemand(k).getNbVNFs());
            for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
                deltaAvailability[i].resize(xSol[k][i].size());
            }

            bool STOP = false;
            while (!STOP){
                computeDeltaAvailability(REQUIRED_AVAIL, deltaAvailability, sectionAvailability, coeff);  
                int nextSection = -1;
                int nextNode = -1;
                double bestRatio = -1.0;

                /* Search for next vnf to include on placement without satifying the chain availability. */
                for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
                    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                        int v = data.getNodeId(n);
                        if ((chainAvailability + deltaAvailability[i][v] < REQUIRED_AVAIL) && ((xSol[k][i][v]/deltaAvailability[i][v]) > bestRatio)){
                            bestRatio = (xSol[k][i][v]/deltaAvailability[i][v]);
                            nextSection = i;
                            nextNode = v;
                        }
                    }
                }
                /* If a vnf is found, include it. */
                if ((nextSection != -1) && (nextNode != -1)){
                    chainAvailability += deltaAvailability[nextSection][nextNode];
                    sectionAvailability[nextSection] = (1.0 - ((1.0 - sectionAvailability[nextSection])*(1.0 - data.getNode(nextNode).getAvailability())));
                    coeff[nextSection][nextNode] = 0;
                    sectionNodes[nextSection].push_back(nextNode);
                }
                /* If not, stop */
                else{
                    STOP = true;
                }
            }
            
            double lhs = 0.0;
            for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
                for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                    int v = data.getNodeId(n);
                    lhs += (coeff[i][v]*xSol[k][i][v]);
                }
            }

            if (lhs < 1 - EPS){
                IloExpr expr(env);
                for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
                    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                        int v = data.getNodeId(n);
                        if (coeff[i][v] == 1){
                            expr += x[k][i][v];
                        }
                    }
                }
                std::string name = "heurAvailabilityCut";

                IloRange cut(env, 1, expr, IloInfinity, name.c_str());
                std::cout << "Adding " << cut.getName() << std::endl;
                context.addUserCut(cut, IloCplex::UseCutForce, IloFalse);
                incrementAvailabilityCutsHeuristic();
                incrementUsercuts();
                expr.end();
            }
        }
    }
}

/** Computes the availability increment resulted from the instalation of a new vnf. @param CHAIN_AVAIL The chain required availability. @param deltaAvail The matrix to be computed. @param sectionAvail THe current section availabilities. @param coeff The matrix of coefficients storing the possible vnfs to be placed. **/
void Callback::computeDeltaAvailability(const double CHAIN_AVAIL, std::vector< std::vector<double> >& deltaAvail, const std::vector< double >& sectionAvail, const std::vector< std::vector<int> >& coeff){
    for (unsigned int i = 0; i < sectionAvail.size(); i++){
        for (unsigned int v = 0; v < coeff[i].size(); v++){
            /* If node is already placed, forbid inclusion */
            if (coeff[i][v] == 0){
                deltaAvail[i][v] = 10.0;
            }
            else{
                double newSectionAvail = (1.0 - ((1.0 - sectionAvail[i])*(1.0 - data.getNode(v).getAvailability())));
                double newChainAvail = (CHAIN_AVAIL / sectionAvail[i])*newSectionAvail;
                deltaAvail[i][v] = newChainAvail - CHAIN_AVAIL;
            }
        }
    }
}

/** Solves the separation problems for a given integer solution. @note Should only be called within candidate context.**/
void Callback::addLazyConstraints(const Context &context)
{
    try {
        /* Get current integer solution */
        getIntegerSolution(context); 

        /* Check VNF placement availability for each demand */
        for (int k = 0; k < data.getNbDemands(); k++){
            
            /* Compute sections availability and sort them by increasing order */
            std::vector<MapAvailability> sectionAvailability = getAvailabilitiesOfSections(k, xSol);
            std::sort(sectionAvailability.begin(), sectionAvailability.end(), compareAvailability);

            /* Find smallest subset of sections violating the SFC availability. */
            const double REQUIRED_AVAIL = data.getDemand(k).getAvailability(); 
            double chainAvailability = 1.0;
            int index = 0;
            int nbSelectedSections = 0;
            while ((chainAvailability >= REQUIRED_AVAIL) && (index < data.getDemand(k).getNbVNFs())){
                chainAvailability *= sectionAvailability[index].availability;
                nbSelectedSections++;
                index++;
            }
            /* If such subset is found, add lazy constraint. */
            if (chainAvailability < REQUIRED_AVAIL){
                /* Try to lift the separating inequality */
                lift(xSol[k], REQUIRED_AVAIL, sectionAvailability, nbSelectedSections);

                /* Build inequality. */
                IloExpr exp(env);
                for (int s = 0; s < nbSelectedSections; ++s){
                    int i = sectionAvailability[s].section;
                    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                        int v = data.getNodeId(n);
                        if (xSol[k][i][v] < 1 - EPS){
                            exp += x[k][i][v];
                        }
                    }
                }
                IloRange cut(env, 1.0, exp, IloInfinity);
                context.rejectCandidate(cut);
                incrementLazyConstraints();
                exp.end();
            }
        }
    }
    catch (...) {
        throw;
    }
}

/** Tries to add new vnf placements to the current solution without changing its availability violation. @param xSol The current solution for a given demand. @param availabilityRequired The SFC required availability. @param sectionAvailability The current section availabilities. @param nbSections The number of sections that can be modified. **/
void Callback::lift(IloNumMatrix& xSol, const double& availabilityRequired, std::vector<Callback::MapAvailability>& sectionAvailability, const int& nbSections)
{
    for (int s = 0; s < nbSections; ++s){
        int i = sectionAvailability[s].section;
        for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
            int v = data.getNodeId(n);
            /* If the i-th vnf is not placed on node v */
            if (xSol[i][v] < 1 - EPS){
                /* Compute the availability obtained if a i-th vnf was placed on node v*/
                double futureAvailability = 1.0;
                double futureAvailabilityOfSection = sectionAvailability[s].availability;
                for (int j = 0; j < nbSections; ++j){
                    if (s == j){
                        double newFailureRate = (1.0 - sectionAvailability[j].availability)*(1.0 - data.getNode(v).getAvailability());
                        futureAvailabilityOfSection = (1.0 - newFailureRate);
                        futureAvailability *= futureAvailabilityOfSection;
                    }
                    else{
                        futureAvailability *= sectionAvailability[j].availability;
                    }
                }
                /* If the availability would still be violated */
                if (futureAvailability < availabilityRequired){
                    /* Place vnf */
                    xSol[i][v] = 1;
                    sectionAvailability[s].availability = futureAvailabilityOfSection;
                }
            }
        }
    }
}

/** Returns the availability of the i-th section of a SFC demand obtained from an integer solution. @param k The demand id. @param i The section id. @param xSol The current integer solution. **/
double Callback::getAvailabilityOfSection(const int& k, const int& i, const IloNum3DMatrix& xSol) const
{
    double failure_prob = 1.0;
    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
        int v = data.getNodeId(n);
        if (xSol[k][i][v] >= 1 - EPS){
            failure_prob *= (1.0 - data.getNode(v).getAvailability());
        }
    }
    double availability = 1.0 - failure_prob;
    return availability;
}

/* Returns the availabilities of the sections of a SFC demand obtained from an integer solution. */
std::vector<Callback::MapAvailability> Callback::getAvailabilitiesOfSections (const int& k, const IloNum3DMatrix& xSol) const
{   
    std::vector<MapAvailability> sectionAvailability;
    for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
        MapAvailability entry;
        entry.section = i;
        entry.availability = getAvailabilityOfSection(k,i, xSol);
        sectionAvailability.push_back(entry);
    }
    return sectionAvailability;
}

/** Returns the current integer solution. @note Should only be called within candidate context. **/ 
void Callback::getIntegerSolution(const Context &context)
{
    /* Fill solution matrix */
    if (context.getId() == Context::Id::Candidate){
        if (context.isCandidatePoint()) {
            for (int k = 0; k < data.getNbDemands(); k++){
                for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
                    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                        int v = data.getNodeId(n);
                        xSol[k][i][v] = context.getCandidatePoint(x[k][i][v]);
                    }
                }
            }
        }
        else{
            throw IloCplex::Exception(-1, "ERROR: Unbounded solution within callback !");
        }
    }
    else{
        throw IloCplex::Exception(-1, "ERROR: Trying to get integer solution while not in candidate context !");
    }
}

/** Returns the current fractional solution. @note Should only be called within relaxation context. **/
void Callback::getFractionalSolution(const Context &context)    
{
    /* Fill solution matrix */
    if (context.getId() == Context::Id::Relaxation){
        for (int k = 0; k < data.getNbDemands(); k++){
            for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
                for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                    int v = data.getNodeId(n);
                    xSol[k][i][v] = context.getRelaxationPoint(x[k][i][v]);
                }
            }
        }
    }
    else{
       throw IloCplex::Exception(-1, "ERROR: Trying to get fractional solution while not in relaxation context !");
    }
    
}

/** Checks if all placement variables of a given SFC demand are integers. @param k The demand id. @param xSol The current solution. **/
const bool Callback::isIntegerAssignment(const int& k, const IloNum3DMatrix& xSol) const{
    for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
        for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
            int v = data.getNodeId(n);
            if ((xSol[k][i][v] >= EPS)  && (xSol[k][i][v] <= 1 - EPS)){
                return false;
            }
        }
    }
    return true;
}


void Callback::incrementLazyConstraints()
{
    thread_flag.lock();
    ++nbLazyConstraints;
    thread_flag.unlock();
}

void Callback::incrementAvailabilityCutsHeuristic()
{
    thread_flag.lock();
    ++nb_cuts_avail_heuristic;
    thread_flag.unlock();
}

void Callback::incrementUsercuts()
{
    thread_flag.lock();
    ++nbCuts;
    thread_flag.unlock();
}

void Callback::incrementTime(const IloNum time)
{
    thread_flag.lock();
    timeAll += time;
    thread_flag.unlock();
}

bool compareAvailability(Callback::MapAvailability a, Callback::MapAvailability b)
{
    return (a.availability < b.availability);
}