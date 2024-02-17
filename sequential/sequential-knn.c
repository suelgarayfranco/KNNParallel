#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

double dwalltime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double sec = tv.tv_sec + tv.tv_usec / 1000000.0;

    return sec;
}

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

static inline void updateNeighboursRowAx(int *neighbours, double *min_distances, int num_neigh, double *distances, int num_points, int points_b_first_idx, bool against_self){
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

            while (pos < num_neigh && min_distances[(update_mat_offset)+pos] <= distance)
                pos++;

            if (min_distances[(update_mat_offset)+pos] > distance) {
                insertDistance(&neighbours[update_mat_offset], &min_distances[update_mat_offset], distance, (points_b_first_idx + j), num_neigh, pos);
            }

        }
    }
}

static inline void updateNeighboursColAx(int *neighbours, double *min_distances, int num_neigh, double *distances, int num_points, int points_b_first_idx, bool against_self) {
    int i, j, distance_mat_offset, update_mat_offset, pos, initial;
    double distance;

    for (i = 0; i < num_points; i++) {
        distance_mat_offset = i;
        update_mat_offset = i * num_neigh;

        if (against_self)
            initial = i + 1;
        else
            initial = 0;

        for (j = initial; j < num_points; j++) {
            distance = distances[distance_mat_offset + (j * num_points)];

            pos = findPosition(&min_distances[update_mat_offset], distance, num_neigh);

            while (pos < num_neigh && min_distances[(update_mat_offset) + pos] <= distance)
                pos++;

            if (min_distances[(update_mat_offset) + pos] > distance) {
                insertDistance(&neighbours[update_mat_offset], &min_distances[update_mat_offset], distance, (points_b_first_idx + j), num_neigh, pos);
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

    int N;          // Amount of points generated
    int K;          // Number of neighbours to account for
    int BS;
    int L;          // Boundary number for point generation

    int i,j,blocks;

    double tick;     // Arrays to track communication time

    // Data structures for all processes
    double *points, *min_dist, *distances;
    int *near_neigh;


    // CHECKS if parameters are ok
    if ((argc != 5) || ((N = atoi(argv[1])) <= 0) || ((K = atoi(argv[2])) <= 0) || (K >= N) || ((BS = atoi(argv[3])) <= 0) || ((L = atoi(argv[4])) <= 0))
    {
        printf("\nUse: %s N K L \nN: Amount of points to generate\nK: Amount of neighbours to compute for each point. K must be less than N\nBS: Block size for points splitting\nL: Boundary number for point generation\n", argv[0]);
        exit(1);
    }

    if ((N % BS) != 0) {
        printf("\nNumber of points must be divisible by block size");
        exit(1);
    }


    // ================ POINTS GENERATION AND STRUCTURES INITIALIZATION ================
    points = (double *)malloc(sizeof(double) * 3 * N);
    srand(time(0)); // Randomize generation
    generate_points(N, points, L);
    blocks = N/BS;
    min_dist = (double *) malloc(sizeof(double) * K * N);
    near_neigh = (int *) malloc(sizeof(int) * K * N);
    distances = (double *) malloc(sizeof(double) * BS * BS);

    // Initializes minimum distances matrix with a big number grater than the biggest distance that could be computed within the boundaries.
    for (i=0; i < K*N; i++) {
        min_dist[i] = sqrt(3 * pow(L,2)) + 1;
    }


    // ===================================================================================
    // ============================== C O M P U T A T I O N ==============================
    // ===================================================================================


    int column_block_first_idx, row_block_first_idx;

    tick = dwalltime();
    for (j = 0; j < blocks; j++) {
        computeDistances(&points[j*BS*3], &points[j*BS*3], distances, BS, true);
        column_block_first_idx = (j*BS);
        updateNeighboursRowAx(&near_neigh[column_block_first_idx*K], &min_dist[column_block_first_idx*K], K, distances, BS, column_block_first_idx, true);
        updateNeighboursColAx(&near_neigh[column_block_first_idx*K], &min_dist[column_block_first_idx*K], K, distances, BS, column_block_first_idx, true);
    }

    for (j = 0; j < blocks - 1; j++) {
        for (i = 0; i < blocks - j - 1; i++) {
            computeDistances(&points[j*BS*3], &points[(j+i+1)*BS*3], distances, BS, false);

            column_block_first_idx = (j*BS);
            row_block_first_idx = (j+i+1)*BS;

            updateNeighboursRowAx(&near_neigh[column_block_first_idx*K], &min_dist[column_block_first_idx*K], K, distances, BS, row_block_first_idx, false);
            updateNeighboursColAx(&near_neigh[row_block_first_idx*K], &min_dist[row_block_first_idx*K], K, distances, BS, column_block_first_idx, false);
        }
    }

    printf("(N=%d) (K=%d) (BS=%d) (Blocks:%d)\tTotal time=%lf\n", N, K, BS, blocks, dwalltime() - tick);

    // UNCOMMENT IF WANT TO SAVE TO FILES
    /*FILE *fpt1, *fpt2, *fpt3;
    fpt1 = fopen("points.csv", "w+");
    fpt2 = fopen("neighbours.csv", "w+");
    fpt3 = fopen("min-distances.csv", "w+");
    for (i=0; i<N; i++) {
        fprintf(fpt1, "%d,%d,%d\n", (int) points[3*i], (int)points[3*i + 1], (int)points[3*i + 2]);
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

    free(points);
    free(min_dist);
    free(distances);
    free(near_neigh);

    return 0;
}
