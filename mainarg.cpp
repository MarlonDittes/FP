#include <iostream>
#include "argtable3.h"

//shai test

int main(int argc, char* argv[]) {


//    int target_node = -1;
    const char* graph_filename = NULL;

    // Define the argument structures
//    struct arg_int* arg_source = arg_int0("s", "source", "<s>", "Source node");
    struct arg_str* arg_graph = arg_str0("g", "graph", "<file>", "Graph filename");
    struct arg_end* end = arg_end(20);

    // Create an array of pointers to the argument structures
    void* argtable[] = {arg_graph, end }; // arg_target, 

    // Parse the command-line arguments
    int nerrors = arg_parse(argc, argv, argtable);

    // Check for parsing errors
    if (nerrors > 0) {
        arg_print_errors(stderr, end, "my_program");
        arg_print_syntax(stderr, argtable, "\n");
        arg_print_glossary(stderr, argtable, "  %-25s %s\n");
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        exit(1);
    }

    // Get the values of the parsed arguments
 //   target_node = arg_target->count > 0 ? arg_target->ival[0] : -1;
    graph_filename = arg_graph->count > 0 ? arg_graph->sval[0] : NULL;

    // Free the memory allocated by Argtable
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

    // Check if all required arguments are provided
    if (graph_filename == NULL) { //target_node == -1 ||
        printf("Error: Missing required arguments\n");
        exit(1);
    }

 //   printf("Source node: %d\n", source_node);
    printf("Graph filename: %s\n", graph_filename);  


}
