#include "model.hpp"

/** Constructor: Builds and exports the mathematical model to mip.lp file. Also sets up CPLEX parameters. **/
Model::Model(const IloEnv& env_, const Data& data_) : 
                env(env_), model(env), cplex(model), data(data_), 
                obj(env), constraints(env)
{
	std::cout << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "-                  Building optimization model.                 -" << std::endl;
    std::cout << "=================================================================" << std::endl;
    setVariables();
    setObjective();  
    setConstraints();  
    setCplexParameters();
    cplex.exportModel("mip.lp");
    std::cout << std::endl << "Model was correctly built ! " << std::endl;                 
}


/** Set up the Cplex parameters. **/
void Model::setCplexParameters(){
    std::cout << std::endl << "Setting up CPLEX optimization parameters... " << std::endl;
    // build callback
    callback = new Callback(env, data, x, y, secAvail, secUnavail);

    // define contexts on which the callback will be used
    CPXLONG contextmask = 0;
    if (data.getInput().getLazy() == Input::LAZY_ON){
	    contextmask |= IloCplex::Callback::Context::Id::Candidate;
    }
    contextmask |= IloCplex::Callback::Context::Id::Relaxation;

    // activate the callback usage
	cplex.use(callback, contextmask);

    /** Time limit definition **/
    cplex.setParam(IloCplex::Param::TimeLimit, data.getInput().getTimeLimit());
	// Treads limited to one
    cplex.setParam(IloCplex::Param::Threads, 1); 
    
    // cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1); // Uncomment to desactivate CPLEX automatic heuristics
}
/*************************************************************************/
/*                          VARIABLE DEFINITIONS                         */
/*************************************************************************/

/** Set up variables **/
void Model::setVariables(){

    std::cout << "Setting up variables... " << std::endl;
    // Define y variables
    setPlacementVariables();

    // Define x variables
    setAssignmentVariables();

    // If routing is activated, define rounting related variables
    if (data.getInput().getRoutingActivation() == Input::ROUTING_ON){
        setPairAssignmentVariables();
        setRoutingVariables();
        setDelayVariables();
        setArcUsageVariables();
    }

    // If availability approximation is used, define auxiliary availability variables
    if (data.getInput().getApproximationType() != Input::APPROXIMATION_TYPE_NONE){
        setAvailabilityVariables();
    }

    std::cout << std::endl;
}

/** Set up VNF placement variables: For any node v and VNF f, 
    y[v][f] =  1 if VNF f is installed on node v; 0 otherwise **/
void Model::setPlacementVariables(){
    std::cout << "\t > Setting up VNF placement variables... " << std::endl;
    y.resize(lemon::countNodes(data.getGraph()));
    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
        int v = data.getNodeId(n);
        y[v].resize(data.getNbVnfs());
        for (int f = 0; f < data.getNbVnfs(); f++){
            int vnf = data.getVnf(f).getId();
            std::string name = "y(" + std::to_string(v) + "," + std::to_string(vnf) + ")";
            if (data.getInput().isRelaxation()){
                y[v][f] = IloNumVar(env, 0.0, 1.0, ILOFLOAT, name.c_str());
            }
            else{
                y[v][f] = IloNumVar(env, 0.0, 1.0, ILOINT, name.c_str());
            }
            model.add(y[v][f]);
        }
    }
}

/** Set up VNF assignment variables: For any demand k, section i, and node v, 
    x[k][i][v] = 1 if the i-th VNF of SFC k can be processed on node v; 0 otherwise. **/
void Model::setAssignmentVariables(){
    std::cout << "\t > Setting up VNF assignment variables... " << std::endl;
    x.resize(data.getNbDemands());
    for (int k = 0; k < data.getNbDemands(); k++){
        x[k].resize(data.getDemand(k).getNbVNFs());
        for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
            x[k][i].resize(lemon::countNodes(data.getGraph()));
            for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                int v = data.getNodeId(n);
                std::string name = "x(" + std::to_string(v) + "," + std::to_string(i) + "," + std::to_string(data.getDemand(k).getId()) + ")";
                if (data.getInput().isRelaxation()){
                    x[k][i][v] = IloNumVar(env, 0.0, 1.0, ILOFLOAT, name.c_str());
                }
                else{
                    x[k][i][v] = IloNumVar(env, 0.0, 1.0, ILOINT, name.c_str());
                }
                model.add(x[k][i][v]);
            }
        }
    }
}

/** Set up VNF pair assignment variables: For any demand k, section i, nodes s and t,
    z[k][i][s][t] = 1 if for demand k, its i-th VNF is installed on node s and its (i+1)-th VNF is installed on node t **/
void Model::setPairAssignmentVariables(){
    std::cout << "\t > Setting up VNF pair assignment variables... " << std::endl;
    z.resize(data.getNbDemands());
    for (int k = 0; k < data.getNbDemands(); k++){
        const int NB_SECTIONS = data.getDemand(k).getNbVNFs()+1;
        z[k].resize(NB_SECTIONS);
        for (int i = 0; i < NB_SECTIONS; i++){
            const int NB_NODES = lemon::countNodes(data.getGraph());
            z[k][i].resize(NB_NODES);
            for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                int s = data.getNodeId(n);
                z[k][i][s].resize(NB_NODES);
                for (NodeIt n2(data.getGraph()); n2 != lemon::INVALID; ++n2){
                    int t = data.getNodeId(n2);
                    std::string name = "z(" + std::to_string(data.getDemand(k).getId()) + "," + std::to_string(i) + "," + std::to_string(s) + "," + std::to_string(t) + ")";
                    double lb = 0.0;
                    double ub = 1.0;
                    // specific cases where z is fixed
                    if (i == 0 && s != data.getDemand(k).getSource()){
                        ub = lb;
                    }
                    if (i == (NB_SECTIONS-1) && t != data.getDemand(k).getTarget()){
                        ub = lb;
                    }
                    
                    if (data.getInput().isRelaxation()){
                        z[k][i][s][t] = IloNumVar(env, lb, ub, ILOFLOAT, name.c_str());
                    }
                    else{
                        z[k][i][s][t] = IloNumVar(env, lb, ub, ILOINT, name.c_str());
                    }
                    model.add(z[k][i][s][t]);
                }
            }
        }
    }
}

/** Set up SFC routing variables: For any demand k, section i, arc a, and nodes s and t,
    r[k][i][a][s][t] = 1 if the arc a is used for routing the i-th section of demand k from s to t **/
void Model::setRoutingVariables(){
    std::cout << "\t > Setting up SFC routing variables... " << std::endl;
    r.resize(data.getNbDemands());
    for (int k = 0; k < data.getNbDemands(); k++){
        const int NB_SECTIONS = data.getDemand(k).getNbVNFs()+1;
        r[k].resize(NB_SECTIONS);
        for (int i = 0; i < NB_SECTIONS; i++){
            const int NB_ARCS = lemon::countArcs(data.getGraph());
            r[k][i].resize(NB_ARCS);
            for (ArcIt arc_it(data.getGraph()); arc_it != lemon::INVALID; ++arc_it){
                int a = data.getArcId(arc_it);
                const int NB_NODES = lemon::countNodes(data.getGraph());
                r[k][i][a].resize(NB_NODES);
                for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                    int s = data.getNodeId(n);
                    r[k][i][a][s].resize(NB_NODES);
                    for (NodeIt n2(data.getGraph()); n2 != lemon::INVALID; ++n2){
                        int t = data.getNodeId(n2);
                        std::string name = "r(" + std::to_string(data.getDemand(k).getId()) + "," + std::to_string(i) + "," + std::to_string(a) + "," + std::to_string(s) + "," + std::to_string(t) + ")";
                        double lb = 0.0;
                        double ub = 1.0;
                        // specific cases where r is fixed to zero
                        if (i == 0 && s != data.getDemand(k).getSource()){
                            ub = lb;
                        }
                        if (i == (NB_SECTIONS-1) && t != data.getDemand(k).getTarget()){
                            ub = lb;
                        }
                        if (s == t){
                            ub = lb;
                        }
                        // define variable
                        if (data.getInput().isRelaxation()){
                            r[k][i][a][s][t] = IloNumVar(env, lb, ub, ILOFLOAT, name.c_str());
                        }
                        else{
                            r[k][i][a][s][t] = IloNumVar(env, lb, ub, ILOINT, name.c_str());
                        }
                        model.add(r[k][i][a][s][t]);
                    }
                }
            }
        }
    }
}

/** Set up SFC delay variables: For any demand k, and section i,
    delay[k][i] corresponds to the maximal delay that can be obtained in this section **/
void Model::setDelayVariables(){
    std::cout << "\t > Setting up SFC delay variables... " << std::endl;
    delay.resize(data.getNbDemands());
    for (int k = 0; k < data.getNbDemands(); k++){
        const int NB_SECTIONS = data.getDemand(k).getNbVNFs()+1;
        delay[k].resize(NB_SECTIONS);
        for (int i = 0; i < NB_SECTIONS; i++){
            std::string name = "l(" + std::to_string(data.getDemand(k).getId()) + "," + std::to_string(i) + ")";
            delay[k][i] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT, name.c_str());
            model.add(delay[k][i]);
        }
    }
}

/** Set up SFC arc usage variables: For any demand k, section i, and arc a,
    arc_usage[k][i][a] = 1 if arc a is used for routing the i-th section of demand k **/
void Model::setArcUsageVariables(){
    std::cout << "\t > Setting up SFC arc usage variables... " << std::endl;
    arc_usage.resize(data.getNbDemands());
    for (int k = 0; k < data.getNbDemands(); k++){
        const int NB_SECTIONS = data.getDemand(k).getNbVNFs()+1;
        arc_usage[k].resize(NB_SECTIONS);
        for (int i = 0; i < NB_SECTIONS; i++){
            const int NB_ARCS = lemon::countArcs(data.getGraph());
            arc_usage[k][i].resize(NB_ARCS);
            for (ArcIt arc_it(data.getGraph()); arc_it != lemon::INVALID; ++arc_it){
                int a = data.getArcId(arc_it);
                int tail = data.getNodeId(data.getGraph().source(arc_it));
                int head = data.getNodeId(data.getGraph().target(arc_it));
                std::string name = "pi(" + std::to_string(data.getDemand(k).getId()) + "," + std::to_string(i) + "," + std::to_string(a) + ")";
                double lb = 0.0;
                double ub = 1.0;
                // specific cases where arc can be fixed to zero
                if (i == 0 && head == data.getDemand(k).getSource()){
                    ub = lb;
                }
                if (i == (NB_SECTIONS-1) && tail == data.getDemand(k).getTarget()){
                    ub = lb;
                }

                if (data.getInput().isRelaxation()){
                    arc_usage[k][i][a] = IloNumVar(env, lb, ub, ILOFLOAT, name.c_str());
                }
                else{
                    arc_usage[k][i][a] = IloNumVar(env, lb, ub, ILOINT, name.c_str());
                }
                model.add(arc_usage[k][i][a]);
            }
        }
    }
}

/** Set up availability variables: For any demand k, and section i,
    secAvail[k][i] refers to the availability of the i-th section of demand k.**/
void Model::setAvailabilityVariables(){
    std::cout << "\t > Setting up section availability variables... " << std::endl;
    const int NB_DEMANDS = data.getNbDemands();
    secAvail.resize(NB_DEMANDS);
    for (int k = 0; k < NB_DEMANDS; k++){
        const int NB_VNFS = data.getDemand(k).getNbVNFs();
        secAvail[k].resize(NB_VNFS);
        for (int i = 0; i < NB_VNFS; i++){
            std::string name = "secAvail(" + std::to_string(k) + "," + std::to_string(i) + ")";
            const double LB = data.getDemand(k).getAvailability();
            const double UB = data.getParallelAvailability(data.getAvailNodeRank());
            secAvail[k][i] = IloNumVar(env, LB, UB, ILOFLOAT, name.c_str());
            model.add(secAvail[k][i]);
        }
    }

    // secUnavail is an auxiliary variable such that secUnavail = 1 - secAvail
    secUnavail.resize(NB_DEMANDS);
    for (int k = 0; k < NB_DEMANDS; k++){
        const int NB_VNFS = data.getDemand(k).getNbVNFs();
        secUnavail[k].resize(NB_VNFS);
        for (int i = 0; i < NB_VNFS; i++){
            std::string name = "secUnavail(" + std::to_string(k) + "," + std::to_string(i) + ")";
            const double LB = 1.0 - data.getParallelAvailability(data.getAvailNodeRank());
            const double UB = 1.0 - data.getDemand(k).getAvailability();
            secUnavail[k][i] = IloNumVar(env, LB, UB, ILOFLOAT, name.c_str());
            model.add(secUnavail[k][i]);
        }
    }
}

/*************************************************************************/
/*                          OBJECTIVE FUNCTION                           */
/*************************************************************************/

/* Set up objective function. */
void Model::setObjective(){

    std::cout << "Setting up objective function... " << std::endl;

	IloExpr exp(env);
    
    std::cout << "\t > Minimizing VNF placement cost. " << std::endl;
    /* Objective: minimize VNF placement cost */
	for(NodeIt n(data.getGraph()); n != lemon::INVALID; ++n) {
        int v = data.getNodeId(n);
        for (int i = 0; i < data.getNbVnfs(); i++){
            int f = data.getVnf(i).getId();
            double cost = data.getPlacementCost(data.getNode(v), data.getVnf(f));
            exp += ( cost*y[v][f] ); 
        }
    }
    // maximize avail
    // for (int k = 0; k < data.getNbDemands(); k++){
    //     for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
    //         exp -= 10*secAvail[k][i];
    //     }
    // }
	obj.setExpr(exp);
	obj.setSense(IloObjective::Minimize);
    model.add(obj);
	exp.clear();
    exp.end();
}

/*************************************************************************/
/*                          CONSTRAINTS DEFINITION                       */
/*************************************************************************/

/* Set up constraints. */
void Model::setConstraints(){
    std::cout << std::endl << "Setting up constraints... " << std::endl;
    // Placement related constraints
    setVnfAssignmentConstraints();
    setNodeCapacityConstraints();
    switch (data.getInput().getDisaggregatedVnfPlacement()){
        case Input::DISAGGREGATED_VNF_PLACEMENT_OFF:
            setOriginalVnfPlacementConstraints();
        break;
        case Input::DISAGGREGATED_VNF_PLACEMENT_ON:
            setVnfPlacementConstraints();
        break;
        default:
            std::cout << "ERROR: DISAGGREGATED_VNF_PLACEMENT not recognized." << std::endl;
            exit(EXIT_FAILURE);
        break;
    }

    // Routing related constraints
    if (data.getInput().getRoutingActivation() == Input::ROUTING_ON){
        setDelayConstraints();
        setLinkingConstraints();
        setRoutingConstraints();
        //setBandwidthConstraints();
    }

    // Availability approx related constraints
    if (data.getInput().getApproximationType() != Input::APPROXIMATION_TYPE_NONE){
        setSFCAvailabilityApproxConstraints();
        setSectionAvailabilityApproxConstraints();
    }


    if (data.getInput().getStrongNodeCapacity() == Input::STRONG_NODE_CAPACITY_ON){
        setStrongNodeCapacityConstraints();
    }
    
    model.add(constraints);
}

/* Add up the original aggregated VNF placement constraints. */
void Model::setOriginalVnfPlacementConstraints()
{
    std::cout << "\t > Setting up aggregated VNF Placement constraints... " << std::endl;
    for (int f = 0; f < data.getNbVnfs(); f++){
        for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
            int v = data.getNodeId(n);
            IloExpr exp(env);
            // build constraint expression
            for (int k = 0; k < data.getNbDemands(); k++){
                for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
                    int f_ik = data.getDemand(k).getVNF_i(i);
                    if (f_ik == f){
                        exp += x[k][i][v];
                    }
                }
            }
            int bigM = 0;
            for (int k = 0; k < data.getNbDemands(); k++){
                bigM += data.getDemand(k).getNbVNFs();
            }
            exp -= (bigM * y[v][f]);
            // add constraint
            std::string name = "Original_VNF_Placement(" + std::to_string(f) + "," + std::to_string(v) + ")";
            constraints.add(IloRange(env, -IloInfinity, exp, 0, name.c_str()));
            exp.clear();
            exp.end();
        }
    }
}

/* Add up the VNF placement constraints: a VNF can only be assigned to a demand if it is already placed. */
void Model::setVnfPlacementConstraints()
{
    std::cout << "\t > Setting up VNF Placement constraints... " << std::endl;
    for (int k = 0; k < data.getNbDemands(); k++){
        for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
            int f = data.getDemand(k).getVNF_i(i);
            for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                int v = data.getNodeId(n);
                IloExpr exp(env);
                exp += x[k][i][v];
                exp -= y[v][f];
                std::string name = "VNF_Placement(" + std::to_string(k) + "," + std::to_string(i) + "," + std::to_string(v) + ")";
                constraints.add(IloRange(env, -IloInfinity, exp, 0, name.c_str()));
                exp.clear();
                exp.end();
            }
        }
    }
}

/* Add up the VNF assignment constraints: At least lb VNFs must be assigned to each section of each demand. */
void Model::setVnfAssignmentConstraints(){
    std::cout << "\t > Setting up VNF Assignment constraints... " << std::endl;
    for (int k = 0; k < data.getNbDemands(); k++){
        for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
            IloExpr exp(env);
            for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                int v = data.getNodeId(n);
                exp += x[k][i][v];
            }
            std::string name = "VNF_Assignment(" + std::to_string(k) + "," + std::to_string(i) + ")";
            int rhs = data.getMinNbNodes(data.getDemand(k).getAvailability());
            //std::cout << rhs << std::endl;
            if (rhs >= 1){
                constraints.add(IloRange(env, rhs, exp, IloInfinity, name.c_str()));
            }
            else{
                std::cerr << "ERROR: Error within method getMinNbNodes. \n" << std::endl;
                exit(0);
            }
            exp.clear();
            exp.end();
        }
    }
}

/* Add up the node capacity constraints: the bandwidth treated in a node must respect its capacity. */
void Model::setNodeCapacityConstraints(){
    std::cout << "\t > Setting up Node Capacity constraints... " << std::endl;
    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
        int v = data.getNodeId(n);
        IloExpr exp(env);
        double capacity = data.getNode(v).getCapacity();
        for (int k = 0; k < data.getNbDemands(); k++){
            for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
                int vnf = data.getDemand(k).getVNF_i(i);
                double coeff = data.getDemand(k).getBandwidth() * data.getVnf(vnf).getConsumption();
                exp += (coeff * x[k][i][v]);
            }
        }
        std::string name = "Node_Capacity(" + std::to_string(v) + ")";
        constraints.add(IloRange(env, 0, exp, capacity, name.c_str()));
        exp.clear();
        exp.end();
    }
}

/* Add up the strong node capacity constraints. */
void Model::setStrongNodeCapacityConstraints(){
    std::cout << "\t Setting up Strong Node Capacity constraints... " << std::endl;
    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
        int v = data.getNodeId(n);
        double capacity = data.getNode(v).getCapacity();
        for (int f = 0; f < data.getNbVnfs(); f++){
            IloExpr exp(env);
            for (int k = 0; k < data.getNbDemands(); k++){
                for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
                    int vnf = data.getDemand(k).getVNF_i(i);
                    if (vnf == f){
                        double coeff = data.getDemand(k).getBandwidth() * data.getVnf(vnf).getConsumption();
                        exp += (coeff * x[k][i][v]);
                    }
                }
            }
            exp -= ( capacity * y[v][f]);
            std::string name = "Strong_Node_Capacity(" + std::to_string(v) + "," + std::to_string(f) + ")";
            constraints.add(IloRange(env, -IloInfinity, exp, 0, name.c_str()));
            exp.clear();
            exp.end();
        }
    }
}

/* Add up the delay constraints. */
void Model::setDelayConstraints(){
    std::cout << "\t > Setting up Delay constraints... " << std::endl;
    // Section delay
    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
        int s = data.getNodeId(n);
        for (NodeIt n2(data.getGraph()); n2 != lemon::INVALID; ++n2){
            int t = data.getNodeId(n2);
            if (s != t){
                for (int k = 0; k < data.getNbDemands(); k++){
                    const int NB_SECTIONS = data.getDemand(k).getNbVNFs()+1;
                    for (int i = 0; i < NB_SECTIONS; i++){
                        IloExpr exp(env);
                        for (ArcIt arc_it(data.getGraph()); arc_it != lemon::INVALID; ++arc_it){
                            int a = data.getArcId(arc_it);
                            double arc_delay = data.getLink(a).getDelay();
                            exp += (arc_delay * r[k][i][a][s][t]);
                        }
                        exp -= delay[k][i];
                        std::string name = "Section_Delay(" + std::to_string(s) + "," + std::to_string(t) + "," + std::to_string(k) + "," + std::to_string(i) +")";
                        constraints.add(IloRange(env, -IloInfinity, exp, 0, name.c_str()));
                        exp.clear();
                        exp.end();
                    }
                }
            }
        }
    }

    // Total delay
    for (int k = 0; k < data.getNbDemands(); k++){
        IloExpr exp(env);
        double rhs = data.getDemand(k).getMaxLatency();
        const int NB_SECTIONS = data.getDemand(k).getNbVNFs()+1;
        for (int i = 0; i < NB_SECTIONS; i++){
            exp += delay[k][i];
        }
        std::string name = "Delay(" + std::to_string(k) +")";
        constraints.add(IloRange(env, -IloInfinity, exp, rhs, name.c_str()));
        exp.clear();
        exp.end();
    }
}

/* Add up the linking constraints. */
void Model::setLinkingConstraints(){
    std::cout << "\t > Setting up Linking constraints... " << std::endl;
    for (int k = 0; k < data.getNbDemands(); k++){
        const int NB_SECTIONS = data.getDemand(k).getNbVNFs()+1;
        for (int i = 0; i < NB_SECTIONS; i++){
            for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                int s = data.getNodeId(n);
                for (NodeIt n2(data.getGraph()); n2 != lemon::INVALID; ++n2){
                    int t = data.getNodeId(n2);
                    // tail link
                    if (i != 0){
                        IloExpr exp(env);
                        int vnf = i-1;
                        exp += z[k][i][s][t];
                        exp -= x[k][vnf][s];
                        std::string name = "Tail_link(" + std::to_string(s) + "," + std::to_string(t) + "," + std::to_string(k) + "," + std::to_string(i) +")";
                        constraints.add(IloRange(env, -IloInfinity, exp, 0, name.c_str()));
                        exp.clear();
                        exp.end();
                    }
                    // head link
                    if (i != NB_SECTIONS-1){
                        IloExpr exp(env);
                        exp += z[k][i][s][t];
                        exp -= x[k][i][t];
                        std::string name = "Head_link(" + std::to_string(s) + "," + std::to_string(t) + "," + std::to_string(k) + "," + std::to_string(i) +")";
                        constraints.add(IloRange(env, -IloInfinity, exp, 0, name.c_str()));
                        exp.clear();
                        exp.end();
                    }
                    // imposition link
                    if (i == 0){
                        if (s == data.getDemand(k).getSource()){
                            IloExpr exp(env);
                            exp += x[k][i][t];
                            exp -= z[k][i][s][t];
                            std::string name = "Imp_link(" + std::to_string(s) + "," + std::to_string(t) + "," + std::to_string(k) + "," + std::to_string(i) +")";
                            constraints.add(IloRange(env, 0, exp, 0, name.c_str()));
                            exp.clear();
                            exp.end();
                        }
                    }
                    else{
                        if (i == NB_SECTIONS-1){
                            if (t == data.getDemand(k).getTarget()){
                                IloExpr exp(env);
                                exp += x[k][i-1][s];
                                exp -= z[k][i][s][t];
                                std::string name = "Imp_link(" + std::to_string(s) + "," + std::to_string(t) + "," + std::to_string(k) + "," + std::to_string(i) +")";
                                constraints.add(IloRange(env, 0, exp, 0, name.c_str()));
                                exp.clear();
                                exp.end();
                            }
                        }
                        else{
                            IloExpr exp(env);
                            exp += x[k][i-1][s];
                            exp += x[k][i][t];
                            exp -= z[k][i][s][t];
                            std::string name = "Imp_link(" + std::to_string(s) + "," + std::to_string(t) + "," + std::to_string(k) + "," + std::to_string(i) +")";
                            constraints.add(IloRange(env, -IloInfinity, exp, 1.0, name.c_str()));
                            exp.clear();
                            exp.end();
                        }
                    }
                }
            }
        }
    }

    // Routing linking
    for (int k = 0; k < data.getNbDemands(); k++){
        const int NB_SECTIONS = data.getDemand(k).getNbVNFs()+1;
        for (int i = 0; i < NB_SECTIONS; i++){
            for (ArcIt arc_it(data.getGraph()); arc_it != lemon::INVALID; ++arc_it){
                int a = data.getArcId(arc_it);
                for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                    int s = data.getNodeId(n);
                    for (NodeIt n2(data.getGraph()); n2 != lemon::INVALID; ++n2){
                        int t = data.getNodeId(n2);
                        if (s != t){
                            IloExpr exp(env);
                            exp += r[k][i][a][s][t];
                            exp -= z[k][i][s][t];
                            std::string name = "Route_link(" + std::to_string(k) + "," + std::to_string(i) + "," +  std::to_string(a) + "," + std::to_string(s) + "," +  std::to_string(t) + ")";
                            constraints.add(IloRange(env, -IloInfinity, exp, 0, name.c_str()));
                            exp.clear();
                            exp.end();
                        }
                    }
                }
            }
        }
    }
}

/* Add up the routing constraints. */
void Model::setRoutingConstraints(){
    std::cout << "\t > Setting up Routing constraints... " << std::endl;
    // tail route
    for (int k = 0; k < data.getNbDemands(); k++){
        const int NB_SECTIONS = data.getDemand(k).getNbVNFs()+1;
        for (int i = 0; i < NB_SECTIONS; i++){
            for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                int s = data.getNodeId(n);
                for (NodeIt n2(data.getGraph()); n2 != lemon::INVALID; ++n2){
                    int t = data.getNodeId(n2);
                    if (s != t){
                        for (NodeIt n3(data.getGraph()); n3 != lemon::INVALID; ++n3){
                            int v = data.getNodeId(n3);
                            IloExpr exp(env);
                            for (OutArcIt arc_it(data.getGraph(), n3); arc_it != lemon::INVALID; ++arc_it){
                                int a = data.getArcId(arc_it);
                                exp += r[k][i][a][s][t];
                            }
                            for (InArcIt arc_it(data.getGraph(), n3); arc_it != lemon::INVALID; ++arc_it){
                                int a = data.getArcId(arc_it);
                                exp -= r[k][i][a][s][t];
                            }
                            if (v == s){
                                exp -= z[k][i][s][t];
                            }
                            if (v == t){
                                exp += z[k][i][s][t];
                            }
                            std::string name = "Routing_tail(" + std::to_string(s) + "," + std::to_string(t) + "," + std::to_string(k) + "," + std::to_string(i) + "," + std::to_string(v) +")";
                            constraints.add(IloRange(env, 0, exp, 0, name.c_str()));
                            exp.clear();
                            exp.end();
                        }
                    }
                }
            }
        }
    }
}

/* Add up the bandwidth constraints. */
void Model::setBandwidthConstraints(){
    std::cout << "\t Setting up Bandwidth constraints... " << std::endl;
    for (ArcIt arc_it(data.getGraph()); arc_it != lemon::INVALID; ++arc_it){
        int a = data.getArcId(arc_it);
        IloExpr exp(env);
        double rhs = data.getLink(a).getBandwidth();
        for (int k = 0; k < data.getNbDemands(); k++){
            const int NB_SECTIONS = data.getDemand(k).getNbVNFs()+1;
            for (int i = 0; i < NB_SECTIONS; i++){
                double demand_band = data.getDemand(k).getBandwidth();
                exp += (demand_band * arc_usage[k][i][a]);
            }
        }
        std::string name = "Band(" + std::to_string(a) +")";
        constraints.add(IloRange(env, -IloInfinity, exp, rhs, name.c_str()));
        exp.clear();
        exp.end();
    }

    // Linking band
    for (int k = 0; k < data.getNbDemands(); k++){
        const int NB_SECTIONS = data.getDemand(k).getNbVNFs()+1;
        for (int i = 0; i < NB_SECTIONS; i++){
            for (ArcIt arc_it(data.getGraph()); arc_it != lemon::INVALID; ++arc_it){
                int a = data.getArcId(arc_it);
                for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                    int s = data.getNodeId(n);
                    for (NodeIt n2(data.getGraph()); n2 != lemon::INVALID; ++n2){
                        int t = data.getNodeId(n2);
                        if (s != t){
                            IloExpr exp(env);
                            exp += r[k][i][a][s][t];
                            exp -= arc_usage[k][i][a];
                            std::string name = "Link_band(" + std::to_string(s) + "," + std::to_string(t) + "," + std::to_string(k) + "," + std::to_string(i) + "," + std::to_string(a) + ")";
                            constraints.add(IloRange(env, -IloInfinity, exp, 0, name.c_str()));
                            exp.clear();
                            exp.end();
                        }
                    }
                }
            }
        }
    }
}

void Model::setSFCAvailabilityApproxConstraints(){
    std::cout << "\t > Setting up section availability piecewise linear approximation constraints... " << std::endl;

    std::cout << "\t > Setting up approximated SFC availability constraints... " << std::endl;
    // The availability of a chain should be at least its required SLA
    for (int k = 0; k < data.getNbDemands(); k++){
        IloExpr exp(env);
        IloNumArray breakpoints(env);
        IloNumArray slopes(env);
        // define the piecewise linear approximation parameters
        buildApproximationFunctionAvail(k, breakpoints, slopes);
        for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
            exp += IloPiecewiseLinear(secAvail[k][i], breakpoints, slopes, 1, 0);
        }
        std::string name = "ReqAvail(" + std::to_string(k) + ")";
        double rhs = std::log(data.getDemand(k).getAvailability());
        constraints.add(IloRange(env, rhs, exp, IloInfinity, name.c_str()));
        exp.clear();
        exp.end();
    }

    // Define link between avail and unavail
    std::cout << "\t > Setting up availability linking constraints... " << std::endl;
    for (int k = 0; k < data.getNbDemands(); k++){
        for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
            IloExpr exp(env);
            exp += secUnavail[k][i];
            exp += secAvail[k][i];
            std::string name = "availLink(" + std::to_string(k) + "," + std::to_string(i) + ")";
            constraints.add(IloRange(env, 1, exp, 1, name.c_str()));
            exp.clear();
            exp.end();
        }
    }
}

void Model::setSectionAvailabilityApproxConstraints(){

    std::cout << "\t > Setting up approximated section unavailability constraints... " << std::endl;
    // The unavailability of a section should be at least the product of the unavailabilities of its nodes
    for (int k = 0; k < data.getNbDemands(); k++){ 
        IloNumArray breakpoints(env);
        IloNumArray slopes(env);
        // define the piecewise linear approximation parameters
        buildApproximationFunctionUnavail(k, breakpoints, slopes);
        for (int i = 0; i < data.getDemand(k).getNbVNFs(); i++){
            IloExpr exp(env);
            exp += IloPiecewiseLinear(secUnavail[k][i], breakpoints, slopes, 1, 0);
            for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
                int v = data.getNodeId(n);
                double coeff = std::log(1.0 - data.getNode(v).getAvailability());
                exp -= coeff * x[k][i][v];
            }
            std::string name = "SectionAvail(" + std::to_string(k) + "," + std::to_string(i) + ")";
        
            constraints.add(IloRange(env, 0, exp, 0, name.c_str()));
            exp.clear();
            exp.end();
        }
    }
}

void Model::buildApproximationFunctionAvail(int demand, IloNumArray &breakpoints, IloNumArray &slopes){
    //initialize parameters
    breakpoints.clear();
    slopes.clear();

    // get touch points and breakpoints
    IloNumArray touchs(env);
    buildAvailTouchs(demand, touchs);
    buildAvailBreakpoints(demand, breakpoints, touchs);

    // get slopes
    switch (data.getInput().getApproximationType()){
        case Input::APPROXIMATION_TYPE_RESTRICTION:
            slopes.add( 1.0/touchs[0] );
            for (int i = 0; i < touchs.getSize()-1; i++){
                double delta_X = touchs[i+1] - touchs[i];
                double delta_Y = std::log(touchs[i+1]) - std::log(touchs[i]);
                slopes.add( delta_Y/delta_X );
            }
            slopes.add( 0 );
            
            break;
        case Input::APPROXIMATION_TYPE_RELAXATION:
            for (unsigned int i = 0; i < touchs.getSize(); i++){
                double tangent = 1.0/touchs[i];
                slopes.add(tangent);
            }
            break;
        default:
            throw IloCplex::Exception(-1, "ERROR: Unexpected availability approx !");
    }
}


/* Set up the vector u for approximating log(avail). */
void Model::buildAvailTouchs(int demand, IloNumArray &touchs){
    const int NB_BREAKS  = data.getInput().getNbBreakpoints();
    int       NB_TOUCHS  = NB_BREAKS;

    // If approximation from above
    if (data.getInput().getApproximationType() ==  Input::APPROXIMATION_TYPE_RELAXATION){
        NB_TOUCHS = NB_BREAKS + 1;
    }

    // get upper and lower bounds for avail
    const double MIN_AVAIL = data.getDemand(demand).getAvailability();
    const double LB = MIN_AVAIL;
    const double UB = 1.0;
    // const double UB = data.getParallelAvailability(data.getNMostAvailableNodes(5));
    for (int t = 1; t <= NB_TOUCHS; t++){
        double exponent  = ((double) (NB_TOUCHS - t)) / (NB_TOUCHS - 1);
        double touch_val = UB * std::pow((LB / UB), exponent);
        touchs.add(touch_val);
    }
    // std::cout << "Avail touchs: " << touchs << std::endl;
}

/* Set up the breakpoints for approximating log(avail). */
void Model::buildAvailBreakpoints(int demand, IloNumArray &breakpoints, IloNumArray &touchs){
    
    if (data.getInput().getApproximationType() == Input::APPROXIMATION_TYPE_RELAXATION){
        for (int k = 1; k < touchs.getSize(); k++){
            double log_u_k1 = std::log(touchs[k]);
            double log_u_k = std::log(touchs[k-1]);
            double inv_u_k1 = 1.0/touchs[k];
            double inv_u_k = 1.0/touchs[k-1];
            
            breakpoints.add( (log_u_k1 - log_u_k) / (inv_u_k - inv_u_k1) );
        }
    }
    else{
        for (int k = 0; k < touchs.getSize(); k++){
            breakpoints.add(touchs[k]);
        }
    }
}

void Model::buildApproximationFunctionUnavail(int demand, IloNumArray &breakpoints, IloNumArray &slopes){
    //initialize parameters
    breakpoints.clear();
    slopes.clear();
    
    // get touch points and breakpoints
    IloNumArray touchs(env);
    buildUnavailTouchs(demand, touchs);
    buildUnavailBreakpoints(demand, breakpoints, touchs);

    switch (data.getInput().getApproximationType()){
        case Input::APPROXIMATION_TYPE_RESTRICTION:
            slopes.add(1.0/breakpoints[0]);
            for (int i = 0; i < breakpoints.getSize()-1; i++){
                double delta_X = breakpoints[i+1] - breakpoints[i];
                double delta_Y = std::log(breakpoints[i+1]) - std::log(breakpoints[i]);
                slopes.add( delta_Y/delta_X );
            }
            slopes.add(0);
            break;
        case Input::APPROXIMATION_TYPE_RELAXATION:
            for (int i = 0; i < touchs.getSize(); i++){
                double tangent = 1.0/touchs[i];
                slopes.add(tangent);
            }
            break;
        default:
            throw IloCplex::Exception(-1, "ERROR: Unexpected availability approx !");
    }
    // std::cout << "touchs: " << touchs << std::endl;
    // std::cout << "breakpoints: " << breakpoints << std::endl;
    // std::cout << "slopes: " << slopes << std::endl;
}
/* Set up the vector u for approximating log(unavail). */
void Model::buildUnavailTouchs(int demand, IloNumArray &touchs){

    const int NB_BREAKS  = data.getInput().getNbBreakpoints();
    int       NB_TOUCHS  = NB_BREAKS;

    // If approximation from above
    if (data.getInput().getApproximationType() ==  Input::APPROXIMATION_TYPE_RELAXATION){
        NB_TOUCHS = NB_BREAKS + 1;
    }

    const double minAvail = data.getDemand(demand).getAvailability();
    // const double minAvail = data.getParallelAvailability(data.getNMostAvailableNodes(1));
    const double maxAvail = data.getParallelAvailability(data.getAvailNodeRank());
    // const double maxAvail = data.getParallelAvailability(data.getNMostAvailableNodes(5));
    const double EPSILON_PRECISION = 1e-8;
    double UB = std::min(1.0 - minAvail, 1.0 - EPSILON_PRECISION);
    double LB = std::max(1.0 - maxAvail, EPSILON_PRECISION);
    // std::cout <<"UB : " << UB << std::endl;
    // std::cout <<"LB : " << LB << std::endl;

    /** billionnet spacement of touchs **/
    for (int k = 1; k <= NB_TOUCHS; k++){
        double expo = ((double) (NB_TOUCHS - k)) / (NB_TOUCHS - 1);
        // std::cout <<"exp : " << expo << std::endl;
        double u = (UB) * std::pow((LB/UB), expo);
        // std::cout <<"touch : " << u << std::endl;
        touchs.add(u);
    }

    /** linearly spaced touchs **/
    // double interval = (UB - LB)/(NB_TOUCHS-1);
    // for (int k = 0; k < NB_TOUCHS; k++){
    //     double u = LB + interval*k;
    //     touchs.add(u);
    // }


    touchs.add(1.0);
    // std::cout <<"Unavail touchs : " << touchs << std::endl;
}

/* Set up the breakpoints for approximating log(unavail). */
void Model::buildUnavailBreakpoints(int demand, IloNumArray &breakpoints, IloNumArray &touchs){
    // std::cout << "NBtouchs " << touchs.getSize()-1 << std::endl;
    if (data.getInput().getApproximationType() == Input::APPROXIMATION_TYPE_RELAXATION){
        for (int k = 1; k <= touchs.getSize()-1; k++){
            double log_u_k1 = std::log(touchs[k]);
            double log_u_k = std::log(touchs[k-1]);
            double inv_u_k1 = 1.0/touchs[k];
            double inv_u_k = 1.0/touchs[k-1];
            breakpoints.add( (log_u_k1 - log_u_k) / (inv_u_k - inv_u_k1) );
            
            // std::cout << "break[" << k-1 << "] = " << breakpoints[k-1] << std::endl;
        }
    }
    else{
        for (int k = 0; k < touchs.getSize(); k++){
            breakpoints.add(touchs[k]);
        }
    }
    // std::cout << "end breaks" << std::endl;
}


void Model::run()
{
	std::cout << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "-                Running optimization procedure.                -" << std::endl;
    std::cout << "=================================================================" << std::endl;
    
    time = cplex.getCplexTime();
	cplex.solve();
	time = cplex.getCplexTime() - time;
}

int Model::getNbAvailViolation(){
    int nbViolations = 0;
    for (int k = 0; k < data.getNbDemands(); k++){
        if (getServiceAvail(k) + 1e-12 < data.getDemand(k).getAvailability()){
            nbViolations++;
        }
    }
    return nbViolations;
}

double Model::getMaxAvailViolation(){
    double maxViolation = 0.0;
    for (int k = 0; k < data.getNbDemands(); k++){
        double violation = data.getDemand(k).getAvailability() - getServiceAvail(k);
        if (violation > maxViolation + 1e-12){
            maxViolation = violation;
        }
    }
    return maxViolation;
}

void Model::printResult(){
    
	std::cout << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "-                 Printing best solution found.                 -" << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "Printing VNF placement..." << std::endl;
    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
        int v = data.getNodeId(n);
        std::string vnfs;
        for (int f = 0; f < data.getNbVnfs(); f++){
            if (cplex.getValue(y[v][f]) > 1.0 - EPS){
                vnfs += data.getVnf(f).getName();
                vnfs += ", ";
            }
        }
        if (!vnfs.empty()){
            vnfs.pop_back();
            vnfs.pop_back();
            vnfs += ".";
            std::cout << "\t" << data.getNode(v).getName() << ": " << vnfs << std::endl;
        }
    }

    std::cout << std::endl << "Printing Service Function Chain deployment..." << std::endl;
    for (int k = 0; k < data.getNbDemands(); k++){
        printDemand(k);
        std::cout << std::endl;
    }

    std::cout << std::endl << "Printing optimization informations..."                   << std::endl;
    std::cout << "\t Objective value:           " << cplex.getValue(obj)                << std::endl;
    std::cout << "\t Nodes evaluated:           " << cplex.getNnodes()                  << std::endl;
    std::cout << "\t User cuts added:           " << callback->getNbUserCuts()          << std::endl;
    std::cout << "\t Lazy constraints added:    " << callback->getNbLazyConstraints()   << std::endl;
    std::cout << "\t Time on cuts:              " << callback->getTime()                << std::endl;
    std::cout << "\t Total time:                " << time << std::endl                  << std::endl;

}
void Model::printDemand(const int demand){
    std::cout << "Demand " << demand << ": " << std::endl;
    std::cout << "\t Placement: " << getServiceAvail(demand) << " > " << data.getDemand(demand).getAvailability();
    if (getServiceAvail(demand) < data.getDemand(demand).getAvailability()) std::cout << "  NOT OK";
    std::cout << std::endl;
    for (int i = 0; i < data.getDemand(demand).getNbVNFs(); i++){
        printSectionPlacement(demand, i);
        if (data.getInput().getApproximationType() != Input::APPROXIMATION_TYPE_NONE){
        printSectionAvailability(demand, i);
        }
    }
    // print routing
    if (data.getInput().getRoutingActivation() == Input::ROUTING_ON){
        std::cout << "\t Routing: " << std::endl;
        for (int i = 0; i < data.getDemand(demand).getNbVNFs()+1; i++){
            printRouting(demand, i);
        }
    }
}

void Model::printSectionAvailability(const int demand, const int section){
    std::cout << std::setprecision(17) << "\t\t Computed Avail: " << 1.0 - cplex.getValue(secUnavail[demand][section]) << "; Real Avail: ";
    double unavail = 1.0;
    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
        int v = data.getNodeId(n);
        if (cplex.getValue(x[demand][section][v]) > 1.0 - EPS ){
            unavail *= (1.0 - data.getNode(v).getAvailability());
        }
    }
    std::cout << std::setprecision(17) <<  1.0 - unavail << std::endl;
}

double Model::getServiceAvail(const int demand){
    double result = 1.0;
    for (int i = 0; i < data.getDemand(demand).getNbVNFs(); ++i){
        double unavail = 1.0;
        for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
            int v = data.getNodeId(n);
            if (cplex.getValue(x[demand][i][v]) > 1.0 - EPS ){
                unavail *= (1.0 - data.getNode(v).getAvailability());
            }
        }
        double avail = 1.0 - unavail;
        result *= avail;
    }
    return result;
}
void Model::printSectionPlacement(const int demand, const int section){
    std::cout << "\t Section " << section << ": " ;
    std::string placement = "";
    for (NodeIt n(data.getGraph()); n != lemon::INVALID; ++n){
        int v = data.getNodeId(n);
        if (cplex.getValue(x[demand][section][v]) > 1.0 - EPS ){
            placement += data.getNode(v).getName();
            placement += ", ";
        }
    }
    std::cout << placement << std::endl;
}
void Model::printRouting(const int demand, const int section){
    std::cout << "\t Section " << section << ": " ;
    for (NodeIt source(data.getGraph()); source != lemon::INVALID; ++source){
        int s = data.getNodeId(source);
        for (NodeIt target(data.getGraph()); target != lemon::INVALID; ++target){
            int t = data.getNodeId(target);
            if (cplex.getValue(z[demand][section][s][t]) > 1.0 - EPS ){
                printPath(demand, section, source, target);
            }
        }
    }
}

void Model::printPath(const int demand, const int section, Graph::Node &source, Graph::Node &target){
    Graph::Node node = source;
    int s = data.getNodeId(source);
    int t = data.getNodeId(target);
    if (section == 0){
        std::cout << "\t\t [o] "; 
    }
    else{
        std::cout << "\t\t (" << data.getVnf(data.getDemand(demand).getVNF_i(section-1)).getName() << ")";
    }
    while (node != target){
        int currentNode = data.getNodeId(node);
        std::cout << currentNode << " -- ";
        int nextNode = currentNode;
        for (OutArcIt arc(data.getGraph(), node); arc != lemon::INVALID; ++arc){
            int a = data.getArcId(arc);
            if (cplex.getValue(r[demand][section][a][s][t]) > 1.0 - EPS ){
                node = data.getGraph().target(arc);
                nextNode = data.getNodeId(node);
            }
        }
        if (currentNode == nextNode){
            std::cerr << std::endl << "ERROR: Could not find next node on path." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    if (section == data.getDemand(demand).getNbVNFs()){
        std::cout << t << "[d];" << std::endl;
    }
    else{
        std::cout << t << "(" << data.getVnf(data.getDemand(demand).getVNF_i(section)).getName() << ");" << std::endl;
    }
}
void Model::output(){
    std::cout << "Writting results to file..." << std::endl;
    std::string output_file = data.getInput().getOutputFile();
    if (output_file.empty()){
        std::cout << "Warning: There is no output file." << std::endl;
        return;
    }

	std::ofstream fileReport(output_file, std::ios_base::app); // File report
    // If file_output can't be opened 
    if(!fileReport)
    {
        std::cerr << "ERROR: Unable to access output file." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string demand_name = data.getInput().getDemandFile();
    std::string node_name   = data.getInput().getNodeFile();
    std::string link_name   = data.getInput().getLinkFile();
    std::string vnf_name    = data.getInput().getVnfFile();
    std::string relax_type  = "";
    if (data.getInput().getApproximationType() == Input::APPROXIMATION_TYPE_NONE){
        relax_type = "NO_APPROX";
    }
    if (data.getInput().getApproximationType() == Input::APPROXIMATION_TYPE_RESTRICTION){
        relax_type = "RESTRICTION_" + std::to_string(data.getInput().getNbBreakpoints());
    }
    if (data.getInput().getApproximationType() == Input::APPROXIMATION_TYPE_RELAXATION){
        relax_type = "RELAX_" + std::to_string(data.getInput().getNbBreakpoints());
    }
    std::string instance_name = node_name + "_" + demand_name;

    fileReport << link_name << ";"
    		   << node_name << ";"
    		   << demand_name << ";"
    		   << vnf_name << ";"
    		   << relax_type << ";"
    		   << time << ";"
    		   << cplex.getObjValue() << ";"
    		   << cplex.getBestObjValue() << ";"
    		   << cplex.getMIPRelativeGap()*100 << ";"
    		   << cplex.getNnodes() << ";"
    		   << cplex.getNnodesLeft() << ";" 
               << callback->getNbLazyConstraints() << ";" 
    		   << callback->getNbUserCuts() << ";" 
    		   << callback->getTime() << ";" 
               << getNbAvailViolation() << ";"
               << getMaxAvailViolation() << ";"
               << std::endl;
    		   
    // Finalization ***
    fileReport.close();

}

/****************************************************************************************/
/*										Destructors 									*/
/****************************************************************************************/
Model::~Model(){
    delete callback;
}