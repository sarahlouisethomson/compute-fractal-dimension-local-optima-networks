/**
 
I gratefully acknowledge that parts of this program were written by Lazaros Gallos in an implementation of the "MEMB" algorithm
and were obtained from: https://hmakse.ccny.cuny.edu/software-and-data/

The program has been edited and augmented by Sarah Thomson to compute "deterministic" multifractal dimension spectra for a local optima network 
(based on network edge distance and on fitness distances).

The algorithm is an updated variant of algorithm 1 in [1] and will be described in the upcoming publication 
"The Fractal Geometry of Fitness Landscapes at the Local Optima Level." 

It is a specialised implementation of the "Sandbox" algorithm, which was proposed for fractal analysis of complex networks in [2].

To compile: gcc mf1.c -o multifractal-analysis -lm

... or "make" with the makefile provided.

To run you will need a local optima network in Pajek format including fitness values for each node.

Nodes must be named as 0 - (n-1). 

Each line of the input network text file will be: NODENAME1 NODEFITNESS1 NODENAME2 NODEFITNESS2\n 

To run: ./multifractal-analysis INPUTNETWORK.txt N OUTPUTFILE.txt NETWORKDIAMETER NUMBERCENTRES DISTANCESTABLE.TXT

where N is the number of network nodes; NETWORKDIAMETER is the diameter; NUMBERCENTRES is the number of sandbox centres, and DISTANCESTABLE.TXT is a text file containing a matrix which is N*N (N being network size) containing the pairwise distances between nodes.

An example: ./multifractal-analysis had12-l0-p3.txt 133 had12-l0-p3fractal.txt 6 50 had12-l0-p3.distancetable

[1] Thomson, S. L., Verel, S., Ochoa, G., Veerapen, N., & Cairns, D. (2018, July). Multifractality and dimensional determinism in local optima networks. In Proceedings of the Genetic and Evolutionary Computation Conference (pp. 371-378).

[2] Liu, Jin-Long, Zu-Guo Yu, and Vo Anh. "Determination of multifractal dimensions of complex networks by means of the sandbox algorithm." Chaos: An Interdisciplinary Journal of Nonlinear Science 25.2 (2015): 023103.

**/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

typedef struct node
{
    int data;
    struct node *next;
    float fitness;
    
} *NODEPTR;

int *ivector(long nl, long nh);
float *fvector(int size);
void clear_lists(int cluster_size);
void free_ivector(int *v, long nl, long nh);
float mean(int m, int a[]);
int d = 0;
struct node ** neigh_list;
struct node ** full_vertex_set; 
int * sandbox_center_nodes;
int diameter_of_network; 
int * all_vertices;
int *neighbors;
int *box_sizes;
float *fitnesses;
float *boxed_fitnesses;
int * size_of_boxes;
int cycle_count = 0 ;

int **readmatrix(size_t *rows, size_t *cols, const char *filename)
{
    if(rows == NULL || cols == NULL || filename == NULL)
        return NULL;

    *rows = 0;
    *cols = 0;

    FILE *fp = fopen(filename, "r");
    int **matrix = NULL, **tmp;

    char line[1024];

    while(fgets(line, sizeof line, fp))
    {
        if(*cols == 0)
        {
            // determine the size of the columns based on
            // the first row
            char *scan = line;
            int dummy;
            int offset = 0;
            while(sscanf(scan, "%d%n", &dummy, &offset) == 1)
            {
                scan += offset;
                (*cols)++;
            }
        }

        tmp = realloc(matrix, (*rows + 1) * sizeof *matrix);

        if(tmp == NULL)
        {
            fclose(fp);
            return matrix; // return all you've parsed so far
        }

        matrix = tmp;

        matrix[*rows] = calloc(*cols, sizeof *matrix[*rows]);

        if(matrix[*rows] == NULL)
        {
            fclose(fp);
            if(*rows == 0) // failed in the first row, free everything
            {
                fclose(fp);
                free(matrix);
                return NULL;
            }

            return matrix; // return all you've parsed so far
        }

        int offset = 0;
        char *scan = line;
        for(size_t j = 0; j < *cols; ++j)
        {
            if(sscanf(scan, "%d%n", matrix[*rows] + j, &offset) == 1)
                scan += offset;
            else
                matrix[*rows][j] = 0; // could not read, set cell to 0
        }

        // incrementing rows
        (*rows)++;
    }

    fclose(fp);

    return matrix;
}

float mean(int m, int a[]) {
    int sum=0, i;
    for(i=0; i<m; i++)
        sum+=a[i];
    return((float)sum/m);
}

float *fvector(int size)
/* allocate an float vector of size n1 */
{
    float *v;
    v=(float *)malloc((size_t) ((size)*sizeof(float)));
    if (!v) {//printf("allocation failure in ivector()\n");exit(-1);
    }
    return v;
}
int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
    int *v;
    v=(int *)malloc((size_t) ((nh-nl+2)*sizeof(int)));
    if (!v) {//printf("allocation failure in ivector()\n");exit(-1);
    }
    return v-nl+1;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
    free((char*) (v+nl-1));
}

float SANDBOX(char *fileName, int cluster_size, int l_b, float q_frommain, int r_frommain, float e_frommain, int number_centres, int **matrix){
    
    cycle_count++;
    int edge_counter = 0;
    int source_node,destination_node;
    int source_fitness, destination_fitness;
    FILE *fp;
    fp = fopen(fileName,"r");
    neighbors=ivector(0,cluster_size-1);
    box_sizes=ivector(0, number_centres);
    int sizee = cluster_size-1;
    fitnesses=fvector(cluster_size);
    size_of_boxes = ivector(0,cluster_size);
    sandbox_center_nodes = ivector(0, number_centres);
    for(int i=0;i<number_centres; i++) sandbox_center_nodes[i]=0;
    for(int i=0;i<cluster_size;i++) neighbors[i]=0;
    for(int i=0;i<number_centres;i++) box_sizes[i]=0;
    for(int i=0;i<cluster_size;i++) size_of_boxes[i]=0;
    for(int i=0;i<cluster_size;i++) fitnesses[i]=0.0;
    
    if(cycle_count == 1){
        boxed_fitnesses=fvector(sizee);
        
        for(int i=0;i<sizee;i++) boxed_fitnesses[i]=0.0;
    }

    source_node = 0;
    source_fitness = 0;
    destination_fitness = 0;
    destination_node = 0;
    
    while(!feof(fp)){
      
        fscanf(fp,"%d %d %d %d\n",&source_node,&source_fitness,&destination_node,&destination_fitness);
        neighbors[source_node]++;
        fitnesses[source_node]++;
        neighbors[destination_node]++;
        fitnesses[destination_node]++;
        source_fitness = fabs(source_fitness);
        destination_fitness = fabs(destination_fitness);
        
        if(cycle_count == 1){
            boxed_fitnesses[source_node] = source_fitness;
            boxed_fitnesses[destination_node] = destination_fitness;
        }
    }
    fclose(fp); 

    for(int vert = 0; vert < number_centres; vert++){
    int random_number = rand() % (cluster_size - 1) + 1;
        sandbox_center_nodes[vert] = random_number;
      //printf("id of centre node: %d\n", random_number);
    }
  
    int current_sandbox_center;
 
   for(int iterator = 0; iterator < number_centres; iterator++){
      
       current_sandbox_center = sandbox_center_nodes[iterator];  
      
       size_of_boxes[current_sandbox_center] = 0;
       
        d = r_frommain;
        
        for(int node_possibly_within_sphere = 0; node_possibly_within_sphere < cluster_size; node_possibly_within_sphere ++){
        
           int current_node = 0; 
        
                current_node = node_possibly_within_sphere;
                int while_count = 0;
                    
                   
                    float fitness_one = boxed_fitnesses[node_possibly_within_sphere];
                   
                    fitness_one = fabs(fitness_one);
                    
                    float fitness_two = boxed_fitnesses[current_sandbox_center];
                    
                    fitness_two = fabs(fitness_two);
                    
                    float fitness_discrepancy = log(fitness_one/ fitness_two);
                    
                    fitness_discrepancy = fabs(fitness_discrepancy);
                    
                    if(matrix[current_sandbox_center][current_node] == 1){
                    
                        size_of_boxes[current_sandbox_center]++;
                        continue;
                      	 
                    }
            
                    else{
                    
                    if( matrix[current_sandbox_center][current_node] == d - 1 && fitness_discrepancy < e_frommain){
                        
                        size_of_boxes[current_sandbox_center]++;
           
                    	}
                        
                   }
                
            }
 
        }

   for(int i = 0; i < number_centres; i++){
	box_sizes[i] = size_of_boxes[sandbox_center_nodes[i]];
 	//printf("size of one box: %d\n", size_of_boxes[sandbox_center_nodes[i]]);
   }

    float mean_nodes_in_bubble = mean(number_centres, box_sizes);

    free_ivector(neighbors,0,cluster_size-1);
    free_ivector(size_of_boxes, 0,cluster_size);
    free_ivector(box_sizes, 0, number_centres);
    free_ivector(sandbox_center_nodes, 0, number_centres-1);
    
    //printf("mean box size: %f\n", mean_nodes_in_bubble);
    
    return mean_nodes_in_bubble;
}

int main(int argc, char *argv[])
{
    int size, lb;
    float mb;
    FILE *fp;
    FILE *distances;
    size = atoi(argv[2]);
    int number_centres = atoi(argv[5]);
    size_t cols, rows;
    cols = size;
    rows = size;
    int **matrix = readmatrix(&rows, &cols, argv[6]);    
    all_vertices = ivector(0, size);
  
    printf("%s", argv[1]);
    diameter_of_network = atoi(argv[4]);
 
    for(int k = 0; k < size; k++){
        all_vertices[k] = k;
    }

    int diameter_end = 10;
    if(diameter_end > diameter_of_network - 2){
        diameter_end = diameter_of_network - 2;
    }
    int combinations = 20*(diameter_end)*2;
    printf("%d combinations.\n", combinations);
    
    float q_passed = 3.00;
    for(int i = 0; i < 20; i++){
        int r_passed = 2;
	 for(int j = 0; j < diameter_end; j++){
         float e_passed = 0.05;
	 	 for(int k = 0; k < 2; k++){
                mb = (double)size;
                lb = 0;
                lb ++;
                mb = SANDBOX(argv[1], size, lb, q_passed, r_passed, e_passed, number_centres, matrix);               
                double dimension = 0.0;
                double proportional_observed_detail = (double)mb/size;
            	 //printf("proportion of network: %f\n", proportional_observed_detail);
                double non_logarithmic_top_fractional =  (double)pow(proportional_observed_detail, q_passed-1);
                double logarithmic_top_fractional = (double)log(non_logarithmic_top_fractional);
                double rad = (double)r_passed;
                double diam = (double)diameter_of_network;
                double denominator = (q_passed-1)*log(rad/diam);
                dimension = (double)(logarithmic_top_fractional)/denominator;
                fp = fopen(argv[3], "a");
                fprintf(fp, "%f\n",dimension);
                fclose(fp);
                printf("fractal dimension: %f\n", dimension);
                e_passed = e_passed + 0.05;
                        }
                r_passed = r_passed + 1;
                    }
            q_passed = q_passed + 0.3;
                }
    
    return 0;
}
