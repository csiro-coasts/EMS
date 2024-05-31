// nanoflann_wrappers.h

#include <stdio.h> // Include for printf


#ifdef __cplusplus
extern "C" {
#endif

// typedef struct {
//   double x;
//   double y;
//   double z;
//   double *v;
// } point;

typedef struct KDTreeContext KDTreeContext;
void* create_kdtree(const point* points, size_t N, size_t max_leaf_size);
int find_nearest_vertex(void* context, const float query_pt[3]); 
void delete_kdtree(void* context);

#ifdef __cplusplus
}
#endif