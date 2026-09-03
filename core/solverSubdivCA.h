#pragma once

#include "patchSubdivisionCache.h"
#include "utils.h"

#include <openGJK/openGJK.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>

template <typename ParamObj1, typename ParamObj2, typename ParamBound1, typename ParamBound2>
class SolverSubdivCA {
    using ControlPoints1 = std::array<Vector3d, ParamObj1::cntCp>;
    using ControlPoints2 = std::array<Vector3d, ParamObj2::cntCp>;
    using PatchCache1 = PatchSubdivisionCache<ParamObj1, ParamBound1>;
    using PatchCache2 = PatchSubdivisionCache<ParamObj2, ParamBound2>;

  public:
    static double solveCCD(const ParamObj1 &position1, const ParamObj1 &velocity1,
                           const ParamObj2 &position2, const ParamObj2 &velocity2,
                           Array2d &uv1, Array2d &uv2, const double distanceTolerance,
                           const double upperTime = DeltaT) {
        if (distanceTolerance <= 0.0 || upperTime <= 0.0)
            return -1.0;

        struct PatchPairNode {
            std::size_t node1;
            std::size_t node2;
            double candidateTime = 0.0;
            bool needsAdvance = true;
            double spatialWidth = 1.0;

            bool operator<(const PatchPairNode &other) const {
                if (candidateTime != other.candidateTime)
                    return candidateTime > other.candidateTime;
                return spatialWidth > other.spatialWidth;
            }
        };

        PatchCache1 cache1(position1, velocity1);
        PatchCache2 cache2(position2, velocity2);
        std::priority_queue<PatchPairNode> heap;
        heap.push({0, 0, 0.0, true,
                   std::max(cache1[0].bound.width(), cache2[0].bound.width())});

        const double terminalWidth = 0.5 * std::sqrt(distanceTolerance);
        while (!heap.empty()) {
            PatchPairNode current = heap.top();
            heap.pop();

            if (current.needsAdvance) {
                const double competingTime =
                    heap.empty() ? std::numeric_limits<double>::infinity()
                                 : heap.top().candidateTime;
                const AdvanceResult advance = advanceSubpatchPair(
                    cache1[current.node1].position, cache1[current.node1].velocity,
                    cache2[current.node2].position, cache2[current.node2].velocity,
                    distanceTolerance, current.candidateTime, upperTime, competingTime);
                if (advance.time < 0.0)
                    continue;
                if (advance.time > current.candidateTime) {
                    current.candidateTime = advance.time;
                    current.needsAdvance = !advance.readyForSubdivision;
                    heap.push(current);
                    continue;
                }
                current.candidateTime = advance.time;
                current.needsAdvance = false;
            }

            if (cache1[current.node1].bound.width() < terminalWidth &&
                cache2[current.node2].bound.width() < terminalWidth) {
                refineContactParameters(cache1, cache2, current.candidateTime, current.node1,
                                        current.node2, distanceTolerance, uv1, uv2);
                return current.candidateTime;
            }

            const auto children1 = cache1.subdivide(current.node1);
            const auto children2 = cache2.subdivide(current.node2);
            for (const std::size_t child1 : children1) {
                for (const std::size_t child2 : children2) {
                    const double spatialWidth = std::max(cache1[child1].bound.width(),
                                                         cache2[child2].bound.width());
                    heap.push({child1, child2, current.candidateTime, true, spatialWidth});
                }
            }
        }
        return -1.0;
    }

  private:
    template <std::size_t PointCount> struct PolytopeStorage {
        std::array<std::array<gkFloat, 3>, PointCount> coordinates;
        std::array<gkFloat *, PointCount> pointers;
    };

    struct HullClosestPoints {
        bool valid = false;
        double distance = std::numeric_limits<double>::infinity();
        Vector3d normalFrom1To2 = Vector3d::Zero();
    };

    struct AdvanceResult {
        double time = -1.0;
        bool readyForSubdivision = false;
    };

    template <std::size_t PointCount>
    static gkPolytope buildPolytope(const std::array<Vector3d, PointCount> &position,
                                    const std::array<Vector3d, PointCount> &velocity,
                                    const double time,
                                    PolytopeStorage<PointCount> &storage) {
        for (std::size_t i = 0; i < PointCount; ++i) {
            const Vector3d point = position[i] + time * velocity[i];
            storage.coordinates[i] = {static_cast<gkFloat>(point.x()),
                                      static_cast<gkFloat>(point.y()),
                                      static_cast<gkFloat>(point.z())};
            storage.pointers[i] = storage.coordinates[i].data();
        }

        gkPolytope polytope{};
        polytope.numpoints = static_cast<int>(PointCount);
        polytope.coord = storage.pointers.data();
        return polytope;
    }

    static HullClosestPoints hullClosestPoints(const ControlPoints1 &position1,
                                                const ControlPoints1 &velocity1,
                                                const ControlPoints2 &position2,
                                                const ControlPoints2 &velocity2,
                                                const double time) {
        PolytopeStorage<ParamObj1::cntCp> storage1;
        PolytopeStorage<ParamObj2::cntCp> storage2;
        const gkPolytope polytope1 = buildPolytope(position1, velocity1, time, storage1);
        const gkPolytope polytope2 = buildPolytope(position2, velocity2, time, storage2);
        gkSimplex simplex{};
        const double distance = compute_minimum_distance(polytope1, polytope2, &simplex);

        HullClosestPoints result;
        if (!std::isfinite(distance))
            return result;
        result.valid = true;
        result.distance = distance;
        if (distance > 1e-12) {
            const Vector3d point1(simplex.witnesses[0][0], simplex.witnesses[0][1],
                                  simplex.witnesses[0][2]);
            const Vector3d point2(simplex.witnesses[1][0], simplex.witnesses[1][1],
                                  simplex.witnesses[1][2]);
            result.normalFrom1To2 = (point2 - point1) / distance;
        }
        return result;
    }

    static HullClosestPoints hullClosestPointsAt(const PatchCache1 &cache1,
                                                 const PatchCache2 &cache2,
                                                 const std::size_t node1,
                                                 const std::size_t node2,
                                                 const double time) {
        return hullClosestPoints(cache1[node1].position, cache1[node1].velocity,
                                 cache2[node2].position, cache2[node2].velocity, time);
    }

    static void refineContactParameters(
        PatchCache1 &cache1, PatchCache2 &cache2, const double time,
        const std::size_t initialNode1, const std::size_t initialNode2,
        const double parameterTolerance, Array2d &uv1, Array2d &uv2) {
        struct RefinementNode {
            std::size_t node1;
            std::size_t node2;
            double distanceLowerBound = 0.0;
            double spatialWidth = 1.0;

            bool operator<(const RefinementNode &other) const {
                if (distanceLowerBound != other.distanceLowerBound)
                    return distanceLowerBound > other.distanceLowerBound;
                return spatialWidth > other.spatialWidth;
            }
        };

        const HullClosestPoints initialClosest =
            hullClosestPointsAt(cache1, cache2, initialNode1, initialNode2, time);
        std::priority_queue<RefinementNode> heap;
        heap.push({initialNode1, initialNode2,
                   initialClosest.valid ? initialClosest.distance : 0.0,
                   std::max(cache1[initialNode1].bound.width(),
                            cache2[initialNode2].bound.width())});

        while (!heap.empty()) {
            const RefinementNode current = heap.top();
            heap.pop();

            if (cache1[current.node1].bound.width() < parameterTolerance &&
                cache2[current.node2].bound.width() < parameterTolerance) {
                uv1 = cache1[current.node1].bound.centerParam();
                uv2 = cache2[current.node2].bound.centerParam();
                return;
            }

            if (cache1[current.node1].bound.width() >=
                cache2[current.node2].bound.width()) {
                const auto children = cache1.subdivide(current.node1);
                for (const std::size_t child : children) {
                    const HullClosestPoints closest =
                        hullClosestPointsAt(cache1, cache2, child, current.node2, time);
                    const double lowerBound = closest.valid ? closest.distance : 0.0;
                    const double spatialWidth = std::max(
                        cache1[child].bound.width(), cache2[current.node2].bound.width());
                    heap.push({child, current.node2, lowerBound, spatialWidth});
                }
            } else {
                const auto children = cache2.subdivide(current.node2);
                for (const std::size_t child : children) {
                    const HullClosestPoints closest =
                        hullClosestPointsAt(cache1, cache2, current.node1, child, time);
                    const double lowerBound = closest.valid ? closest.distance : 0.0;
                    const double spatialWidth = std::max(
                        cache1[current.node1].bound.width(), cache2[child].bound.width());
                    heap.push({current.node1, child, lowerBound, spatialWidth});
                }
            }
        }

        uv1 = cache1[initialNode1].bound.centerParam();
        uv2 = cache2[initialNode2].bound.centerParam();
    }

    template <std::size_t PointCount>
    static double projectedVelocityUpperBound(
        const std::array<Vector3d, PointCount> &velocity, const Vector3d &direction) {
        double maximum = -std::numeric_limits<double>::infinity();
        for (const Vector3d &pointVelocity : velocity)
            maximum = std::max(maximum, pointVelocity.dot(direction));
        return maximum;
    }

    static AdvanceResult advanceSubpatchPair(
        const ControlPoints1 &position1, const ControlPoints1 &velocity1,
        const ControlPoints2 &position2, const ControlPoints2 &velocity2,
        const double distanceTolerance, const double lowerTime, const double upperTime,
        const double competingTime) {
        constexpr int MaxIterations = 1024;
        constexpr double VelocityEpsilon = 1e-12;

        double time = lowerTime;
        double lastEvaluatedTime = lowerTime;
        for (int iteration = 0; iteration < MaxIterations; ++iteration) {
            const HullClosestPoints closest =
                hullClosestPoints(position1, velocity1, position2, velocity2, time);
            if (!closest.valid)
                return {};

            lastEvaluatedTime = time;
            if (closest.distance <= 0.5 * distanceTolerance)
                return {time, true};

            const double closingVelocity =
                projectedVelocityUpperBound(velocity1, closest.normalFrom1To2) +
                projectedVelocityUpperBound(velocity2, -closest.normalFrom1To2);
            if (closingVelocity <= VelocityEpsilon)
                return {};

            time += closest.distance / closingVelocity;
            if (!std::isfinite(time) || time < 0.0 || time > upperTime)
                return {};
            if (time > competingTime)
                return {time, false};
        }

        return {lastEvaluatedTime, false};
    }
};
