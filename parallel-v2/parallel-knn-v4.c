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

static inline void updateNeighbours(int *neighbours, double *min_distances, int num_neigh, double *distances, int num_points, int points_b_first_idx, bool against_self){
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

void transposeMatrix(double *mat, double *transpose, int dimension){
    int i,j, offset;
    for (i=0; i<dimension; i++) {
        offset = i*dimension;
        for (j=0; j<dimension; j++) {
            transpose[j*dimension + i] = mat[offset + j];
        }
    }
}

void saveDistancesToMat(double *main_mat, double *distances, int dimension){
    int i;
    for (i=0; i<dimension*dimension; i++) {
        main_mat[i] = distances[i];
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

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    int N;          // Amount of points generated
    int K;          // Number of neighbours to account for
    int L;          // Boundary number for point generation
    int world_rank;       // Process rank
    int numProcs;   // Amount of processes
    int i,j,stripsize;

    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Variables to track time
    double local_start, local_elapsed, elapsed;

    // Data structures for all processes
    double *my_points, *other_points, *min_dist, *my_distances, *dist_mat, *dist_trasp_mat;
    int *near_neigh;


    // Check if parameters are fine
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
    if (world_rank == COORDINATOR) {
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
    my_distances = (double *)malloc(sizeof(double) * stripsize * N);
    dist_mat = (double *)malloc(sizeof(double) * stripsize * stripsize);
    dist_trasp_mat = (double *)malloc(sizeof(double) * stripsize * stripsize);

    //Define communicators
    MPI_Comm communicators[numProcs];
    communicators[0] = MPI_COMM_WORLD;

    int color;
    for (i=1; i<numProcs; i++) {
        if (world_rank < i) {
            color = MPI_UNDEFINED;
        }
        else {
            color = 0;
        }
        MPI_Comm_split(MPI_COMM_WORLD,color,world_rank, &communicators[i]);
    }

    // Initializes minimum distances matrix with a big number grater than the biggest distance that could be computed within the boundaries.
    double max_val = sqrt(3 * pow(L,2)) + 1;
    for (i=0; i < K*stripsize; i++) {
        min_dist[i] = max_val;
    }

    // Ensures every process has all structures initialized
    MPI_Barrier(MPI_COMM_WORLD);


    // ===================================================================================
    // ============================== C O M P U T A T I O N ==============================
    // ===================================================================================

    local_start = MPI_Wtime();
    // ================ 1ST STAGE - Scatter generated points ================
    MPI_Scatter(my_points, 3*stripsize, MPI_DOUBLE, my_points, 3*stripsize, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);


    // ================ 2ND STAGE - Broadcast and Compute Distances ================
    i = 0;
    while (world_rank >= i) {

        if (world_rank == i) {
            for (j=0; j<stripsize*3; j++) {
                other_points[j] = my_points[j];
            }
        }

        if (i != numProcs - 1) {
            MPI_Bcast(other_points, 3*stripsize, MPI_DOUBLE, 0, communicators[i]);

            computeDistances(my_points, other_points, &my_distances[i*stripsize*stripsize], stripsize, world_rank == i);

            //saveDistancesToMat(&my_distances[i*stripsize*stripsize], dist_mat, stripsize);

            if (world_rank != i) {
                transposeMatrix(&my_distances[i*stripsize*stripsize], dist_trasp_mat, stripsize);
            }

            if (world_rank == i) {
                MPI_Gather(&my_distances[i*stripsize*stripsize], stripsize*stripsize, MPI_DOUBLE, &my_distances[i*stripsize*stripsize], stripsize*stripsize, MPI_DOUBLE, 0,  communicators[i]);
            } else {
                MPI_Gather(dist_trasp_mat, stripsize*stripsize, MPI_DOUBLE, &my_distances[i*stripsize*stripsize], stripsize*stripsize, MPI_DOUBLE, 0,  communicators[i]);
            }

        }

        else {
            computeDistances(my_points, other_points, &my_distances[i*stripsize*stripsize], stripsize, world_rank == i);
            //saveDistancesToMat(&my_distances[i*stripsize*stripsize], dist_mat, stripsize);
        }

        i++;
    }


    // ================ 3RD STAGE - Update Neighbours & Min Distances ================
    for (i=0; i<numProcs; i++) {
        updateNeighbours(near_neigh, min_dist, K, &my_distances[i*stripsize*stripsize], stripsize, i*stripsize, world_rank == i);
    }


    // ================ 4TH STAGE - Gather all data in COORDINATOR ================
    MPI_Gather(min_dist, K*stripsize, MPI_DOUBLE, min_dist, K*stripsize, MPI_DOUBLE, COORDINATOR, MPI_COMM_WORLD);
    MPI_Gather(near_neigh, K*stripsize, MPI_INT, near_neigh, K*stripsize, MPI_INT, COORDINATOR, MPI_COMM_WORLD);

    local_elapsed = MPI_Wtime() - local_start;

    MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, COORDINATOR, MPI_COMM_WORLD);

    if (world_rank == COORDINATOR) {
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
    free(my_distances);
    free(dist_mat);
    free(dist_trasp_mat);
    free(other_points);

    MPI_Finalize();

    return 0;
}
