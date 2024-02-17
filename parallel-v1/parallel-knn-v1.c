#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "mpi.h"

#define COORDINATOR 0

static inline void insertDistance(int *neighbours, double *min_distances, double distance, int neigh_idx, int num_neigh, int pos) {
    //Update min distances and neighbours
    int i;
    for (i = num_neigh - 1; i > pos; i--) {
        min_distances[i] = min_distances[i-1];
        neighbours[i] = neighbours[i-1];
    }
    min_distances[pos] = distance;
    neighbours[pos] = neigh_idx;
}

static inline int findPosition(double *min_distances, double distance, int num_neigh) {
    int left_limit, right_limit, middle;
    bool found;

    left_limit = 0;
    right_limit = num_neigh - 1;
    middle = (int) (left_limit + right_limit)/2;
    found = false;

    // Binary search for position
    while((left_limit < right_limit) && !found) {
        if (min_distances[middle] < distance) {
            left_limit = middle + 1;
            middle = (middle + right_limit) / 2;
        }
        else if (min_distances[middle] > distance) {
            right_limit = middle - 1;
            middle = (middle + left_limit) / 2;
        }
        else
            found = true;
    }

    return middle;
}

static inline void updateNeighbours(int *neighbours, double *min_distances, int num_neigh, double *distances, int num_points, int points_b_first_idx, bool against_self, bool isRoot){
    int i, j, distance_mat_offset, update_mat_offset, pos, steps;
    double distance;

    for (i=0; i < num_points; i++) {
        distance_mat_offset = i*num_points;
        update_mat_offset = i*num_neigh;

        if (against_self)
            steps = i;
        else
            steps = num_points;

        for (j=0; j < steps; j++) {
            distance = distances[distance_mat_offset + j];

            pos = findPosition(&min_distances[update_mat_offset], distance, num_neigh);

            while (pos < (num_neigh - 1) && min_distances[(update_mat_offset)+pos] <= distance)
                pos++;

            if (min_distances[(update_mat_offset)+pos] > distance) {
                insertDistance(&neighbours[update_mat_offset], &min_distances[update_mat_offset], distance, (points_b_first_idx + j), num_neigh, pos);
            }
        }

        if (against_self) {
            distance_mat_offset = i;
            update_mat_offset = i * num_neigh;

            for (j = i + 1; j < num_points; j++) {
                distance = distances[distance_mat_offset + (j * num_points)];

                pos = findPosition(&min_distances[update_mat_offset], distance, num_neigh);

                while (pos < (num_neigh - 1) && min_distances[(update_mat_offset) + pos] <= distance)
                    pos++;

                if (min_distances[(update_mat_offset) + pos] > distance) {
                    insertDistance(&neighbours[update_mat_offset], &min_distances[update_mat_offset], distance, (points_b_first_idx + j), num_neigh, pos);
                }
            }
        }

    }

}

static inline void computeDistances(double *points_a, double *points_b, double *distances, int num_points, bool against_self) {
    int i, j, steps;
    double x, y, z;

    for (i=0; i < num_points; i++) {
        x = points_a[3*i];
        y = points_a[1 + 3*i];
        z = points_a[2 + 3*i];

        if (against_self)
            steps = i;
        else
            steps = num_points;

        for (j=0; j < steps; j ++) {
            distances[i*num_points + j] = sqrt(pow((points_b[3*j] - x),2) + pow((points_b[1 + 3*j] - y),2) + pow((points_b[2 + 3*j] - z),2));
        }
    }

}

void generate_points(int numPoints, double *points, int limit) {
    int i;
    for (i = 0; i < 3*numPoints; i++) {
        int num = rand();
        if (num == 0) {
            points[i] = 0;
        }
        else if ((num % limit) == 0){
            points[i] = limit;
        }
        else
            points[i] = (rand() % limit);
    }
}

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    int N;          // Amount of points generated
    int K;          // Number of neighbours to account for
    int L;          // Boundary number for point generation
    int rank;       // Process rank
    int numProcs;   // Amount of processes
    int i,j,stripsize;

    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double local_start, local_finish, local_elapsed, elapsed;

    // Data structures for all processes
    double *my_points, *other_points, *min_dist, *distances;
    int *near_neigh;

    // CHECKS if parameters are ok
    if ((argc != 4) || ((N = atoi(argv[1])) <= 0) || ((K = atoi(argv[2])) <= 0) || (K >= N) ||  ((L = atoi(argv[3])) <= 0))
    {
        printf("\nUse: %s N K L \nN: Amount of points to generate\nK: Amount of neighbours to compute for each point. K must be less than N\nL: Boundary number for point generation\n", argv[0]);
        exit(1);
    }

    if ((N % numProcs) != 0) {
        printf("Number of points N must be divisible by the number of processes P");
        exit(1);
    }

    stripsize = N/numProcs;

    if (K > stripsize) {
        printf("Number of neighbours must be less than N/P.");
        exit(1);
    }

    // ================ POINTS GENERATION AND STRUCTURES INITIALIZATION ================
    if (rank == COORDINATOR) {
        my_points = (double *)malloc(sizeof(double) * 3 * N);
        srand(time(0)); // Randomize generation
        generate_points(N, my_points, L);
        min_dist = (double *) malloc(sizeof(double) * K * N);
        near_neigh = (int *) malloc(sizeof(int) * K * N);
    } else {
        my_points = (double *)malloc(sizeof(double) * 3 * stripsize);
        min_dist = (double *) malloc(sizeof(double) * K * stripsize);
        near_neigh = (int *) malloc(sizeof(int) * K * stripsize);
    }

    other_points = (double *)malloc(sizeof(double) * 3 * stripsize);
    distances = (double *)malloc(sizeof(double) * stripsize * stripsize);

    // Initializes minimum distances matrix with a big number grater than the biggest distance that could be computed within the boundaries.

    for (i=0; i < K*stripsize; i++) {
        min_dist[i] = sqrt(3 * pow(L,2)) + 1;
    }

    // Ensures every process has all structures initialized
    MPI_Barrier(MPI_COMM_WORLD);


    // ===================================================================================
    // ============================== C O M P U T A T I O N ==============================
    // ===================================================================================

    // =============== Data scattering ===============

    local_start = MPI_Wtime();

    MPI_Scatter(my_points, stripsize * 3, MPI_DOUBLE, my_points, stripsize * 3, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);


    int other_block_first_idx;

    for (i=0; i < numProcs; i++) {

        if (rank == COORDINATOR) {
            int offset = i*stripsize*3;
            for (j=0; j<stripsize*3; j++) {
                other_points[j] = my_points[offset + j];
            }
        }

    
        MPI_Bcast(other_points, stripsize*3, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
        
        computeDistances(my_points, other_points, distances, stripsize, i == rank);
        
        other_block_first_idx = i * stripsize;
        updateNeighbours(near_neigh, min_dist, K, distances, stripsize, other_block_first_idx, i == rank, rank == COORDINATOR);

    }

    MPI_Gather(min_dist, K*stripsize, MPI_DOUBLE, min_dist, K*stripsize, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
    MPI_Gather(near_neigh, K*stripsize, MPI_INT, near_neigh, K*stripsize, MPI_INT, COORDINATOR, MPI_COMM_WORLD);

    local_finish = MPI_Wtime();
    local_elapsed = local_finish - local_start;
   
    MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, COORDINATOR, MPI_COMM_WORLD);

    if (rank == COORDINATOR) {
        
        printf("(N=%d) (P=%d) (K=%d)\tTotal time=%lf\n", N, numProcs, K, elapsed);

        //UNCOMMENT IF WANT TO SAVE TO FILES
        /*FILE *fpt1, *fpt2, *fpt3;
        fpt1 = fopen("points.csv", "w+");
        fpt2 = fopen("neighbours.csv", "w+");
        fpt3 = fopen("min-distances.csv", "w+");
        for (i=0; i<N; i++) {
            fprintf(fpt1, "%d,%d,%d\n", (int) my_points[3*i], (int)my_points[3*i + 1], (int)my_points[3*i + 2]);
            for(j=0; j<K-1; j++) {
                fprintf(fpt2,"%d,", near_neigh[i*K + j]);
                fprintf(fpt3,"%f,", min_dist[i*K + j]);
            }
            fprintf(fpt2,"%d\n", near_neigh[(i+1)*K -1]);
            fprintf(fpt3,"%f\n", min_dist[(i+1)*K -1]);
        }
        fclose(fpt1);
        fclose(fpt2);
        fclose(fpt3);*/
    }

    free(my_points);
    free(near_neigh);
    free(min_dist);
    free(distances);
    free(other_points);

    MPI_Finalize();

    return 0;
}
