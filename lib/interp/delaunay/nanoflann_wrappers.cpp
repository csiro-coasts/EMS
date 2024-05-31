//nanoflann_wrappers.cpp
#include "ems_conf.h"
#ifdef USE_KDTREE

#include <cstdlib>
#include <ctime>
#include <algorithm> // For std::min
#include "nanoflann.hpp"
#include "utils.h" 

typedef struct {
  double x;
  double y;
  double z;
  double *v;
} point;

extern "C" {
class KDTreeContext {
public:
    using num_t = float;
    PointCloud<num_t> cloud;  // from utils.h
    // Defines KDTree as a type alias -> for nanoflann::KDTreeSingleIndexAdaptor template class from nanoflann namespace, configured for 2D points with L2 (Euclidean) distance, 
    using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<num_t, PointCloud<num_t>>,
        PointCloud<num_t>, 2 /* dim */
    >;
    // declare ptr:
    KDTree* tree;

    // Adjusted Constructor
    KDTreeContext(const point* points, size_t N, size_t max_leaf_size) {
        cloud.pts.resize(N);
        for (size_t i = 0; i < N; ++i) {
            cloud.pts[i].x = points[i].x; 
            cloud.pts[i].y = points[i].y;
            // printf("%i = %1.8f , %1.8f\n",i,cloud.pts[i].x,cloud.pts[i].y);
        }
        tree = new KDTree(2, cloud, {max_leaf_size});
    }

    // Destructor
    ~KDTreeContext() {
        delete tree;
    }
};

void* create_kdtree(const point* points, size_t N, size_t max_leaf_size) {
    return new KDTreeContext(points, N, max_leaf_size);
}

int find_nearest_vertex(void* context, const float query_pt[2]) {
    // convert void* to KDTreeContext*
    KDTreeContext* ctx = static_cast<KDTreeContext*>(context);

    // do we want this safety check or do we want a failure?
    if (!ctx || !ctx->tree) return -1;

    const size_t num_results = 1; // just one NN required
    size_t ret_index; // index of nearest vertex found
    float out_dist_sqr; // Local variable with the distance, still needed for KNNResultSet but not used outside

    // create new result set type, initialise
    nanoflann::KNNResultSet<float> resultSet(num_results);
    resultSet.init(&ret_index, &out_dist_sqr);
    float query_pt_2D[2] = {query_pt[0], query_pt[1]}; 
    // printf("NN %1.8f , %1.8f ",query_pt[0], query_pt[1]);
    // execute search
    ctx->tree->findNeighbors(resultSet, &query_pt_2D[0]);
    // printf(" %i\n", ret_index);
    return static_cast<int>(ret_index);
}

void delete_kdtree(void* context) {
    KDTreeContext* ctx = static_cast<KDTreeContext*>(context);
    delete ctx;
}
}
#endif