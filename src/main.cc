
#include "cmdline.hh"
#include "Spark.hh"
#include <algorithm>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <sys/stat.h>

void seqtoRNA(std::string &sequence) {
    for (char &c : sequence) {
        if (c == 'T') c = 'U';
    }
}

void validate_structure(std::string &seq, std::string &structure) {
    cand_pos_t n = structure.length();
    std::vector<cand_pos_t> pairs;
    for (cand_pos_t j = 0; j < n; ++j) {
        if (structure[j] == '(') pairs.push_back(j);
        if (structure[j] == ')') {
            if (pairs.empty()) {
                std::cout << "Incorrect input: More left parentheses than right" << std::endl;
                exit(0);
            } else {
                cand_pos_t i = pairs.back();
                pairs.pop_back();
                if (seq[i] == 'A' && seq[j] == 'U') {
                } else if (seq[i] == 'C' && seq[j] == 'G') {
                } else if ((seq[i] == 'G' && seq[j] == 'C') || (seq[i] == 'G' && seq[j] == 'U')) {
                } else if ((seq[i] == 'U' && seq[j] == 'G') || (seq[i] == 'U' && seq[j] == 'A')) {
                } else if ((seq[i] == 'A' && seq[j] == 'T') || (seq[i] == 'T' && seq[j] == 'A')) {
                } else {
                    std::cout << "Incorrect input: " << seq[i] << " does not pair with " << seq[j] << std::endl;
                    exit(0);
                }
            }
        }
    }
    if (!pairs.empty()) {
        std::cout << "Incorrect input: More left parentheses than right" << std::endl;
        exit(0);
    }
}

bool exists(const std::string path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

void get_input(std::string file, std::string &sequence, std::string &structure) {
    if (!exists(file)) {
        std::cout << "Input file does not exist" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::ifstream in(file.c_str());
    std::string str;
    cand_pos_t i = 0;
    while (getline(in, str)) {
        if (str[0] == '>') continue;
        if (i == 0) sequence = str;
        if (i == 1) structure = str;
        ++i;
    }
    in.close();
}

/**
 * @brief Simple driver for @see Spark.
 *
 * Reads sequence from command line or stdin and calls folding and
 * trace-back methods of Spark.
 */
int main(int argc, char **argv) {

    args_info args_info;

    // get options (call gengetopt command line parser)
    if (cmdline_parser(argc, argv, &args_info) != 0) {
        exit(1);
    }

    std::string seq;
    if (args_info.inputs_num > 0) {
        seq = args_info.inputs[0];
    } else {
        if (!args_info.input_file_given) std::getline(std::cin, seq);
    }

    std::string restricted = args_info.input_structure_given ? args_info.input_structure_arg: "";

    std::string fileI = args_info.input_file_given ? args_info.input_file_arg : "";

    if (fileI != "") {

        if (exists(fileI)) {
            get_input(fileI, seq, restricted);
        }
        if (seq == "") {
            std::cout << "sequence is missing from file" << std::endl;
        }
    }
    cand_pos_t n = seq.length();
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    if (args_info.noConv_flag) seqtoRNA(seq);
    if (restricted == "") restricted = std::string(n, '.');

    if (restricted.length() != (cand_pos_tu)n) {
        std::cout << "input sequence and structure are not the same size" << std::endl;
        std::cout << seq << std::endl;
        std::cout << restricted << std::endl;
        exit(0);
    }

    if(args_info.paramFile_given){
        std::string file = args_info.paramFile_arg;
        if (exists(file)) vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
        else{
            std::cerr << "Not a valid parameter file!" << std::endl;
            exit(EXIT_FAILURE);
        }
    } else {
        if (seq.find('T') != std::string::npos) {
            vrna_params_load_DNA_Mathews2004();
        } else{
            std::string file = std::string(PARAMS_DIR) + "/rna_DirksPierce09.par";
            if (exists(file)) vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
            else{
                std::cerr << "Not a valid parameter file!" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    bool verbose = args_info.verbose_flag;

    bool mark_candidates = args_info.mark_candidates_given;

    noGU = args_info.noGU_given;
    validate_structure(seq, restricted);

    sparse_tree tree(restricted, n);

    Spark spark(seq, restricted,&tree,args_info.dangles_arg,!args_info.pk_free_flag,args_info.pk_only_flag, !args_info.noGC_given,mark_candidates);

    cmdline_parser_free(&args_info);

    energy_t mfe = spark.fold();
    std::string structure = spark.trace_back();

    std::ostringstream smfe;
    smfe << std::setiosflags(std::ios::fixed) << std::setprecision(2) << mfe / 100.0;
    std::cout << seq << std::endl;
    std::cout << structure << " (" << smfe.str() << ")" << std::endl;

    if (verbose) {
        std::cout << std::endl;
        spark.print_trace_arrows();
        std::cout << std::endl;
        spark.print_candidates();
    }

    return 0;
}