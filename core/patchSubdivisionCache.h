#pragma once

#include "paramBound.h"

#include <array>
#include <cstddef>
#include <vector>

template <typename ParamObj, typename ParamBound> class PatchSubdivisionCache {
  public:
    using ControlPoints = std::array<Vector3d, ParamObj::cntCp>;

    struct Node {
        ParamBound bound;
        ControlPoints position;
        ControlPoints velocity;
        std::array<std::size_t, 4> children{};
        bool subdivided = false;
    };

    PatchSubdivisionCache(const ParamObj &position, const ParamObj &velocity)
        : position_(position), velocity_(velocity) {
        append(ParamBound{});
    }

    const Node &operator[](const std::size_t id) const { return nodes_[id]; }

    std::array<std::size_t, 4> subdivide(const std::size_t id) {
        if (nodes_[id].subdivided)
            return nodes_[id].children;

        const ParamBound parentBound = nodes_[id].bound;
        std::array<std::size_t, 4> children;
        for (int child = 0; child < 4; ++child)
            children[child] = append(parentBound.interpSubpatchParam(child));
        nodes_[id].children = children;
        nodes_[id].subdivided = true;
        return children;
    }

  private:
    const ParamObj &position_;
    const ParamObj &velocity_;
    std::vector<Node> nodes_;

    std::size_t append(const ParamBound &bound) {
        const std::size_t id = nodes_.size();
        nodes_.push_back(
            {bound, position_.divideBezierPatch(bound), velocity_.divideBezierPatch(bound)});
        return id;
    }
};
