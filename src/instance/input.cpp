#include "input.hpp"

/****************************************************************************************/
/*										Constructor										*/
/****************************************************************************************/

/** Constructor. **/
Input::Input(const std::string filename) : parameters_file(filename){
    std::cout << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "-                   Reading input parameters.                   -" << std::endl;
    std::cout << "=================================================================" << std::endl;
    
    node_file   = getParameterValue("nodeFile=");
    link_file   = getParameterValue("linkFile=");
    demand_file = getParameterValue("demandFile=");
    vnf_file    = getParameterValue("vnfFile=");

    disaggregated_vnf_placement = (Disaggregated_VNF_Placement_Constraints)std::stoi(getParameterValue("disaggregated_VNF_Placement="));
    strong_node_capacity        = (Strong_Node_Capacity_Constraints)std::stoi(getParameterValue("strong_node_capacity="));
    node_cover                  = (Node_Cover_Cuts)std::stoi(getParameterValue("node_cover="));
    chain_cover                 = (Chain_Cover_Cuts)std::stoi(getParameterValue("chain_cover="));
    vnf_lower_bound             = (VNF_Lower_Bound_Cuts)std::stoi(getParameterValue("vnf_lower_bound="));
    section_failure_cuts        = (Section_Failure_Cuts)std::stoi(getParameterValue("section_failure="));
    routing_activation          = (Routing)std::stoi(getParameterValue("routing="));
    approx_type                 = (Approximation_Type)std::stoi(getParameterValue("availability_approx="));
    lazy                        = (Lazy_Constraints)std::stoi(getParameterValue("lazy="));
    heuristic_activation        = (Heuristic)std::stoi(getParameterValue("heuristic="));

    linear_relaxation           = std::stoi(getParameterValue("linearRelaxation="));
    time_limit                  = std::stoi(getParameterValue("timeLimit="));
    nb_breakpoints              = std::stoi(getParameterValue("nb_breakpoints="));

    output_file                 = getParameterValue("outputFile=");

    print();
}



/* Returns the pattern value in the parameters file. */
std::string Input::getParameterValue(std::string pattern){
    std::string line;
    std::string value = "";
    std::ifstream param_file (parameters_file.c_str());
    if (param_file.is_open()) {
        while ( std::getline (param_file, line) ) {
            std::size_t pos = line.find(pattern);
            if (pos != std::string::npos){
                value = line.substr(pos + pattern.size());
                if (value.empty()){
                    std::cout << "WARNING: Field '" << pattern << "' is empty." << std::endl; 
                }
                return value;
            }
        }
        param_file.close();
    }
    else {
        std::cerr << "ERROR: Unable to open parameters file '" << parameters_file << "'." << std::endl; 
        exit(EXIT_FAILURE);
    }
    std::cout << "WARNING: Did not found field '" << pattern << "' inside parameters file." << std::endl; 
    return value;
}

/** Print the info stored in the parameter file. */
void Input::print(){
    std::cout << "\t Node File:                     " << node_file    << std::endl;
    std::cout << "\t Link File:                     " << link_file    << std::endl;
    std::cout << "\t Service Chain Function File:   " << demand_file  << std::endl;
    std::cout << "\t Virtual Network Function File: " << vnf_file     << std::endl;
    std::cout << "\t Output File:                   " << output_file  << std::endl;
    std::cout << "\t Linear Relaxation:             ";
    if (linear_relaxation)  std::cout << "TRUE" << std::endl;
    else                    std::cout << "FALSE" << std::endl;
    
    std::cout << "\t Time Limit:                    " << time_limit   << " seconds"   << std::endl;
    std::cout << std::endl;

    std::cout << "\t Lazy Constraints:        " << lazy                         << std::endl;
    std::cout << "\t Heuristics:              " << heuristic_activation         << std::endl;
    std::cout << "\t Disaggregated placement: " << disaggregated_vnf_placement  << std::endl;
    std::cout << "\t Strong capacity:         " << strong_node_capacity         << std::endl;
    std::cout << "\t Node cover:              " << node_cover                   << std::endl;
    std::cout << "\t Chain cover:             " << chain_cover                  << std::endl;
}