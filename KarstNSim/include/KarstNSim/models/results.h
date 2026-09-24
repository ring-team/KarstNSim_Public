#pragma once

#include <KarstNSim/basics.h>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>

namespace KarstNSim {

    struct ResultPoint {
        Vector3 p; //!< The 3D coordinates of the point.
        float cost = 0.0f; //!< The cost associated with the point.
        float equivalent_radius = -99999.0f; //!< The equivalent radius at the point.
        int branch_id = -1; //!< The branch ID associated with the point. Intersections keep branch_id == -1.

        std::vector<float> vadose_flags; //!< Vadose/phreatic flags exported per physical water-table channel.

        float vadose_flag = -99999.0f; //!< Legacy first exported vadose flag, kept for backward compatibility.
        float external_drift = -99999.0f; //!< External drift value at the point.
        float kriging_weight = -99999.0f; //!< Kriging weight associated with the point.

        //! Index of this point's node in the native karst skeleton (KarsticSkeleton::nodes) that produced the
        //! result. Two points with the same node_id are the same skeleton node; use it, not coordinates, to
        //! rebuild topology. UINT32_MAX means unknown. Not part of the legacy text serialization.
        std::uint32_t node_id = UINT32_MAX;
    };

    struct ResultSegment {
        ResultPoint start; //!< The starting point of the segment.
        ResultPoint end; //!< The ending point of the segment.
    };

    struct KarstNetworkResult {
        std::vector<ResultSegment> segments; //!< List of segments in the karst network.

        std::vector<std::string> vadose_property_names; //!< Names of exported vadose flag properties.
        bool has_drift_properties = false; //!< Whether external drift and kriging weight must be written.

        //! Solved inlet/outlet matrix with the same content and layout as the legacy
        //! `<name>_connectivity_matrix.txt` export: rows and columns are laid out exactly like the input
        //! connectivity matrix (one row per sink, one column per spring). Filled only when
        //! ParamsSource::create_solved_connectivity_matrix is true. Not part of the legacy text serialization.
        std::vector<std::vector<int>> solved_connectivity_matrix;

        /*!
        \brief Convert the result to a string representation for ASCII export.
        \return Tab-separated string representation of the karst network.
        */
        std::string to_string() const {
            constexpr float NO_DATA_VALUE = -99999.0f;

            std::ostringstream out;

            out << "Index\tX\tY\tZ\tcost\tequivalent_radius\tbranch_id";

            for (const std::string& property_name : vadose_property_names) {
                out << "\t" << property_name;
            }

            if (has_drift_properties) {
                out << "\texternal_drift\tkriging_weight";
            }

            int index = 0;

            for (const auto& segment : segments) {
                // The segment index is duplicated for both start and end points.
                for (const auto& point : { segment.start, segment.end }) {
                    out << "\n"
                        << index << "\t"
                        << std::fixed << std::setprecision(5)
                        << point.p.x << "\t"
                        << point.p.y << "\t"
                        << point.p.z << "\t"
                        << point.cost << "\t"
                        << point.equivalent_radius << "\t"
                        << std::fixed << std::setprecision(0)
                        << point.branch_id;

                    for (std::size_t i = 0; i < vadose_property_names.size(); ++i) {
                        const float vadose_value =
                            i < point.vadose_flags.size()
                            ? point.vadose_flags[i]
                            : NO_DATA_VALUE;

                        out << "\t" << vadose_value;
                    }

                    if (has_drift_properties) {
                        out << "\t"
                            << std::fixed << std::setprecision(5)
                            << point.external_drift << "\t"
                            << point.kriging_weight;
                    }
                }

                index++;
            }

            return out.str();
        }

        /*!
        \brief Add a segment to the karst network result.
        \param start Starting point of the segment.
        \param end Ending point of the segment.
        */
        void add_segment(const ResultPoint& start, const ResultPoint& end) {
            segments.push_back({ start, end });
        }
    };
}