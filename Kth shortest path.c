//PDC FINAL PROJECT 
//MUHAMMAD BIN BASIT MASOOD 21I-2486
//ALI MAAZ CHEEMA 21I-2486
//SECTION E
                                                                                                                                                                                                                                                                                                                                                                                                                     #include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpi.h>
#include <omp.h>

#define MAX_EDGES 500000
#define MAX_NODES 400000
#define FILENAME "Email-Enron.txt"
#define NUM_PAIRS 10
#define K 3

void readGraphAndInitialize(int *numNodes, int edges[][3], int *numEdges)
{
    FILE *file = fopen(FILENAME, "r");
    if (file == NULL) 
    {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    //printing the first few lines to see if file is being read correctly
    printf("First 15 of the file:\n");
    for (int i = 0; i < 15; i++) 
    {
        char line[256];
        if (fgets(line, sizeof(line), file) != NULL) 
        {
            printf("%s", line);
        } 
        else 
        {
            break; //end of file reached
        }
    }

    //rewinding the file pointer to read from the beginning again
    rewind(file);

    //skipping unnecessary lines
    char line[256];
    while (fgets(line, sizeof(line), file) != NULL) 
    {
        if (sscanf(line, "# Nodes: %d Edges: %d", numNodes, numEdges) == 2) 
        {
            break; //line with node and edge info found
        }
    }

    //checkinf if the line with node and edge info was found
    if (*numNodes == 0 || *numEdges == 0) 
    {
        fprintf(stderr, "Error reading number of nodes and edges\n");
        exit(EXIT_FAILURE);
    }

    //printinf total number of nodes and edges
    printf("Total number of nodes: %d\n", *numNodes);
    printf("Total number of edges: %d\n", *numEdges);

     //reading edges from file
    int currentEdge = 0;
    printf("Edges from file:\n");
    while (fgets(line, sizeof(line), file) != NULL && currentEdge < MAX_EDGES)
    {
    if (sscanf(line, "%d %d", &edges[currentEdge][0], &edges[currentEdge][1]) == 2) 
    {
        edges[currentEdge][2] = 1; // Assuming all edges have weight 1 for simplicity
        printf("%d %d\n", edges[currentEdge][0], edges[currentEdge][1]);
        currentEdge++;
    } 
    else 
    {
        fprintf(stderr, "Error reading edge from file\n");
    }
    }

    if (currentEdge == MAX_EDGES && !feof(file)) 
    {
        fprintf(stderr, "Warning: Maximum number of edges reached, some edges may not be read\n");
    } 
    else if (currentEdge < MAX_EDGES && !feof(file)) 
    {
        fprintf(stderr, "Error: Premature end of file reached\n");
    }

    *numEdges = currentEdge; // Update the total number of edges

    fclose(file);
}

//func to manually enter pairs of nodes
void enterNodePairs(int pairs[][2]) 
{
    printf("Enter %d pairs of source and destination nodes:\n", NUM_PAIRS);
    for (int i = 0; i < NUM_PAIRS; i++)
    {
        printf("Pair %d: ", i + 1);
        scanf("%d %d", &pairs[i][0], &pairs[i][1]);
    }
}

//structure to represent edges
struct Edge 
{
    int dest;
    int cost;
    struct Edge *next;
};

//func to find K shortest path lengths
void findKShortestSerial(int edges[][3], int n, int m, int k, int source, int dest) 
{
    //initializing graph
    struct Edge **g = (struct Edge **)malloc((n + 1) * sizeof(struct Edge *));
    for (int i = 0; i <= n; i++) 
    {
        g[i] = NULL;
    }

    for (int i = 0; i < m; i++) 
    {
        int src = edges[i][0];
        int dst = edges[i][1];
        int cost = edges[i][2];
        
        struct Edge *newEdge = (struct Edge *)malloc(sizeof(struct Edge));
        newEdge->dest = dst;
        newEdge->cost = cost;
        newEdge->next = g[src];
        g[src] = newEdge;
    }

    //vector to store distances
    int **dis = (int **)malloc((n + 1) * sizeof(int *));
    for (int i = 0; i <= n; i++) 
    {
        dis[i] = (int *)malloc(k * sizeof(int));
        for (int j = 0; j < k; j++) 
        {
            dis[i][j] = INT_MAX;
        }
    }

    //initializing priority queue
    struct PriorityQueue 
    {
        int distance;
        int node;
        struct PriorityQueue *next;
    } *pq = NULL;

    //implementating priority queue functions
    void push(struct PriorityQueue **pq, int distance, int node) 
    {
        struct PriorityQueue *temp = (struct PriorityQueue *)malloc(sizeof(struct PriorityQueue));
        temp->distance = distance;
        temp->node = node;
        temp->next = *pq;
        *pq = temp;
    }

    void pop(struct PriorityQueue **pq) 
    {
        struct PriorityQueue *temp = *pq;
        *pq = (*pq)->next;
        free(temp);
    }

    int topDistance(struct PriorityQueue *pq) 
    {
        return pq->distance;
    }

    int topNode(struct PriorityQueue *pq) 
    {
        return pq->node;
    }

    int isEmpty(struct PriorityQueue *pq) 
    {
        return pq == NULL;
    }

    push(&pq, 0, source);
    dis[source][0] = 0;

    //while pq has elements
    while (!isEmpty(pq)) 
    {
        //storing the node value
        int u = topNode(pq);

        //storing the distance value
        int d = topDistance(pq);
        pop(&pq);
        if (dis[u][k - 1] < d)
            continue;
        struct Edge *current = g[u];

        //traversing the adjacency list
        while (current != NULL) 
        {
            int dest = current->dest;
            int cost = current->cost;

            //checking for the cost
            if (d + cost < dis[dest][k - 1]) 
            {
                dis[dest][k - 1] = d + cost;

                //sorting the distances
                for (int i = k - 1; i > 0; i--) 
                {
                    if (dis[dest][i] < dis[dest][i - 1]) 
                    {
                        int temp = dis[dest][i];
                        dis[dest][i] = dis[dest][i - 1];
                        dis[dest][i - 1] = temp;
                    } 
                    else
                    {
                        break;
                    }
                }

                //pushing elements to priority queue
                push(&pq, (d + cost), dest);
            }
            current = current->next;
        }
    }

printf("Source: %d, Destination: %d\n", source, dest);
    //printing K shortest paths
    for (int i = 0; i < k; i++)
    {
        printf("%d ", dis[dest][i]);
    }
printf("\n");

    //freeing memory
    for (int i = 0; i <= n; i++) 
    {
        struct Edge *temp = g[i];
        while (temp != NULL) 
        {
            struct Edge *toDelete = temp;
            temp = temp->next;
            free(toDelete);
        }
    }
    
    free(g);
    for (int i = 0; i <= n; i++) 
    {
        free(dis[i]);
    }
    free(dis);
}

//parallel implementation of K shortest paths within MPI processes
void findKShortestParallel(int edges[][3], int n, int m, int k, int pairs[][2]) 
{
    //initializing MPI
    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //calculating chunk size for each MPI process
    int chunk_size = NUM_PAIRS / world_size;
    int start = rank * chunk_size;
    int end = (rank == world_size - 1) ? NUM_PAIRS : start + chunk_size;

    //computing K shortest paths for each pair of nodes assigned to this MPI process
    #pragma omp parallel for
    for (int i = start; i < end; i++)
    {
        findKShortestSerial(edges, n, m, k, pairs[i][0], pairs[i][1]);
    }
}

int main(int argc, char *argv[]) 
{
    int numNodes, numEdges;
    int edges[MAX_EDGES][3];
    int pairs[NUM_PAIRS][2];

    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //only the root process reads the graph and initializes parameters
    if (rank == 0) 
    {
        //read graph and initialize parameters
        readGraphAndInitialize(&numNodes, edges, &numEdges);
       
        //broadcast numNodes, numEdges, and edges to all processes
        MPI_Bcast(&numNodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&numEdges, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(edges, numEdges * 3, MPI_INT, 0, MPI_COMM_WORLD);

        //prompting the user to enter pairs of nodes
        enterNodePairs(pairs);
    } 
    else 
    {
        //receive numNodes, numEdges, and edges from root process
        MPI_Bcast(&numNodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&numEdges, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(edges, numEdges * 3, MPI_INT, 0, MPI_COMM_WORLD);
    }

    //broadcasting pairs to all processes
    MPI_Bcast(pairs, NUM_PAIRS * 2, MPI_INT, 0, MPI_COMM_WORLD);

    //computing K shortest paths using serial implementation in process 0
    if (rank == 0)
    {
        printf("Serial Execution:\n");
        for (int i = 0; i < NUM_PAIRS; i++) 
        {
            findKShortestSerial(edges, numNodes, numEdges, K, pairs[i][0], pairs[i][1]);
        };

     printf("\n");
     printf("Parallel Execution:\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //compute K shortest paths using parallel implementation
    findKShortestParallel(edges, numNodes, numEdges, K, pairs);

    MPI_Finalize();

    return 0;
}

