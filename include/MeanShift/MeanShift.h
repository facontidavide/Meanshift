#pragma once 

#include <vector>
#include <math.h>
#include "simd_math_prims.h"

template <typename PointType>
struct ClusterType {
    PointType center;
    std::vector<PointType> original_points;
    std::vector<PointType> shifted_points;
};

inline double GaussianKernel(double distance, double kernel_bandwidth)
{
    const double u = distance / kernel_bandwidth;
    const double val = -0.5 * u*u;
    if( val < -10 ) return 0.0;
    return expapprox(val);
}

inline double QuarticKernel(double distance, double kernel_bandwidth)
{
    const double u = std::abs( distance / kernel_bandwidth );
    if( u > 1.0 ) return 0.0;
    const double tmp = ( 1.0 - u*u);
    return (15.0/16.0)*tmp*tmp;
}

inline double ParabolicKernel(double distance, double kernel_bandwidth)
{
    const double u = std::abs( distance / kernel_bandwidth );
    if( u > 1.0 ) return 0.0;
    const double tmp = ( 1.0 - u*u);
    return 0.75*tmp;
}


template <typename PointType>
class MeanShift {
public:
    typedef PointType            Point;
    typedef ClusterType<Point>   Cluster;
    typedef std::vector<Point>   PointsVector;
    typedef std::vector<Cluster> ClustersVector;

    MeanShift():
        _kernel_func(GaussianKernel),
        _cluster_epsilon(0.5),
        _meanshift_epsilon(0.0000001)
    { }

    MeanShift(double (*_kernel_func)(double,double)) { setKernelFunction(QuarticKernel); }
    PointsVector meanshift(const PointsVector& points,
                           double kernel_bandwidth) const;

    ClustersVector cluster(const PointsVector&points, double kernel_bandwidth) const ;

    void setMeanshiftEps(double new_eps) { _meanshift_epsilon = new_eps; }
    void setClusterEps(double new_eps)   { _cluster_epsilon = new_eps; }

    double meanshiftEps()const { return _meanshift_epsilon; }
    double clusterEps()const   { return _cluster_epsilon; }
    void setKernelFunction(double (*_kernel_func)(double,double));

private:
    double (*_kernel_func)(double,double);
    void shift_point(const Point&, const PointsVector&, double, Point&) const;
    ClustersVector cluster(const PointsVector& points,
                           const PointsVector& shifted_points) const;

    double _cluster_epsilon;
    double _meanshift_epsilon;
};


//---------------------------------------------------------------

template <typename PointType> inline
double euclidean_distance_sqr(const PointType&point_a, const PointType&point_b){
    double total = 0;
    for(int j = 0; j<point_a.size(); j++){
        const double temp = (point_a[j] - point_b[j]);
        total += temp*temp;
    }
    return (total);
}

template <typename PointType> inline
double euclidean_distance(const PointType &point_a, const PointType &point_b){
    return sqrt(euclidean_distance_sqr(point_a, point_b));
}

template <typename PointType> inline
void reset_point(PointType& point)
{
    for(int j = 0; j<point.size(); j++){
        point[j] = 0;
    }
}

// equals to ( result += point*multiplier )
template <typename PointType> inline
void multiply_and_accumulate(const PointType& point, double multiplier,  PointType& result)
{
    for(int j = 0; j<point.size(); j++){
        result[j] += point[j] * multiplier;
    }
}

// equals to ( result = point*multiplier )
template <typename PointType> inline
void multiply(const PointType& point, double multiplier,  PointType& result)
{
    for(int j = 0; j<point.size(); j++){
        result[j] = point[j] * multiplier;
    }
}

template <typename PointType> inline
void add(const PointType& point_A, const PointType& point_B,  PointType& result)
{
    for(int j = 0; j<point_A.size(); j++){
        result[j] = point_A[j] + point_B[j];
    }
}

//------------------------------------------------------------

template <typename PointType> inline
void MeanShift<PointType>::setKernelFunction( double (*_kernel_func)(double,double) ) {
    if(!_kernel_func){
        _kernel_func = GaussianKernel;
    } else {
        _kernel_func = _kernel_func;
    }
}

template <typename PointType> inline
void MeanShift<PointType>::shift_point(const Point &point,
                                       const PointsVector &points,
                                       double kernel_bandwidth,
                                       Point &shifted_point) const {
    shifted_point = point;
    reset_point(shifted_point);

    double total_weight = 0;
    for(int i=0; i<points.size(); i++)
    {
        const Point& temp_point = points[i];
        double distance = euclidean_distance(point, temp_point);
        double weight = _kernel_func(distance, kernel_bandwidth);
        multiply_and_accumulate(temp_point, weight, shifted_point);
        total_weight += weight;
    }

    const double total_weight_inv = 1.0/total_weight;
    multiply(shifted_point, total_weight_inv, shifted_point);
}

template <typename PointType> inline
typename MeanShift<PointType>::PointsVector MeanShift<PointType>::meanshift(const PointsVector &points,
                                                                            double kernel_bandwidth) const{
    const double EPSILON_SQR = _meanshift_epsilon*_meanshift_epsilon;
    std::vector<bool> stop_moving(points.size(), false);
    PointsVector shifted_points = points;
    double max_shift_distance;
    Point point_new;
    do {
        max_shift_distance = 0;
        for(int i=0; i<shifted_points.size(); i++)
        {
            if (!stop_moving[i])
            {
                shift_point(shifted_points[i], points, kernel_bandwidth, point_new);
                double shift_distance_sqr = euclidean_distance_sqr(point_new, shifted_points[i]);
                if(shift_distance_sqr > max_shift_distance){
                    max_shift_distance = shift_distance_sqr;
                }
                if(shift_distance_sqr <= EPSILON_SQR) {
                    stop_moving[i] = true;
                }
                shifted_points[i] = point_new;
            }
        }
    } while (max_shift_distance > EPSILON_SQR);
    return shifted_points;
}

template <typename PointType> inline
typename MeanShift<PointType>::ClustersVector MeanShift<PointType>::cluster(const PointsVector& points,
                                                                            const PointsVector& shifted_points) const
{
    const double CLUSTER_EPSILON_SQR = _cluster_epsilon*_cluster_epsilon;
    ClustersVector clusters;

    for (int i = 0; i < shifted_points.size(); i++)
    {
        bool added_to_cluster = false;
        for (int c = 0; c < clusters.size(); c++) {
            if (euclidean_distance_sqr(shifted_points[i], clusters[c].mode) <= CLUSTER_EPSILON_SQR)
            {
                added_to_cluster = true;
                clusters[c].original_points.push_back(points[i]);
                clusters[c].shifted_points.push_back(shifted_points[i]);
                break;
            }
        }

        if ( !added_to_cluster ) {
            Cluster clus;
            clus.center = points[i];
            clus.original_points.push_back(points[i]);
            clus.shifted_points.push_back(shifted_points[i]);
            clusters.push_back(clus);
        }
    }

    for (Cluster& clus: clusters)
    {
        reset_point( clus.center );
        for (const Point& point: clus.shifted_points)
        {
            add( clus.center, point, clus.center );
        }
        double inv_total = 1.0 / static_cast<double>(clus.shifted_points.size());
        multiply( clus.center, inv_total, clus.center);
    }

    return clusters;
}

template <typename PointType> inline
typename MeanShift<PointType>::ClustersVector MeanShift<PointType>::cluster(const PointsVector&points, double kernel_bandwidth) const
{
    PointsVector shifted_points = meanshift(points, kernel_bandwidth);
    return cluster(points, shifted_points);
}
