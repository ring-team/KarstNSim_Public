#pragma once
#include <KarstNSim/basics.h>

namespace KarstNSim {

    struct ResultPoint {
        Vector3 p; //!< The 3D coordinates of the point.
        float cost; //!< The cost associated with the point.
        float equivalent_radius; //!< The equivalent radius at the point.
        int branch_id; //!< The branch ID associated with the point.
        bool vadose_flag; //!< A flag indicating if the point is in the vadose zone.
    };

    struct ResultSegment {
        ResultPoint start; //!< The starting point of the segment.
        ResultPoint end; //!< The ending point of the segment.
    };

    struct KarstNetworkResult {
        std::vector<ResultSegment> segments; //!< List of segments in the karst network.

        //!< Convert the result to a string representation, for easy saving to a file
        std::string to_string() const {
            auto header_line = "Index\tX\tY\tZ\tcost\tequivalent_radius\tbranch_id\tvadose_flag";
            std::string result = header_line;
            int index = 0;
            for (const auto& segment : segments) {
                // index is duped for both start and end points of the segment
                for (const auto& point : { segment.start, segment.end }) {
                    result += "\n" + std::to_string(index) + "\t" +
                        std::to_string(point.p.x) + "\t" +
                        std::to_string(point.p.y) + "\t" +
                        std::to_string(point.p.z) + "\t" +
                        std::to_string(point.cost) + "\t" +
                        std::to_string(point.equivalent_radius) + "\t" +
                        std::to_string(point.branch_id) + "\t" +
                        std::to_string(point.vadose_flag ? 1 : 0);
                }
                index++;
            }
            return result;
        }

        //!< Add a segment to the karst network result.
        void add_segment(const ResultPoint& start, const ResultPoint& end) {
            segments.push_back({ start, end });
        }
    };
}