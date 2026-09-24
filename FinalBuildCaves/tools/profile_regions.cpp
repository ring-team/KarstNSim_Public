#include "fbs/caves.hpp"
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

// Sequential streaming workload. Retain only the current and previous region
// and validate their real seam. Use an external process RSS tool for memory.
int main(int argc,char** argv) {
    try {
        const std::string input=argc>1?argv[1]:"16";
        if(input.empty() || input.find_first_not_of("0123456789")!=std::string::npos)
            throw std::runtime_error("Expected a region count");
        const auto count=std::stoull(input);
        if(count==0 || count>4096) throw std::runtime_error("Count must be 1..4096");
        std::uint64_t caverns=0,tunnels=0,profiles=0,work=0,points=0;
        double length=0;
        fbs::caves::Region previous;
        const auto start=std::chrono::steady_clock::now();
        for(std::uint64_t i=0;i<count;++i) {
            fbs::caves::Request request;request.world_seed=42;
            request.key={static_cast<std::int64_t>(i)-static_cast<std::int64_t>(count/2),-3,0};
            auto region=fbs::caves::generate(request);
            caverns+=region.caverns.size();tunnels+=region.tunnels.size();
            work+=region.work_used;points+=region.support_points;
            if(i) tunnels+=fbs::caves::stitch(previous,region).size();
            for(const auto& tunnel:region.tunnels) {
                profiles+=tunnel.centerline.size();
                for(std::size_t j=1;j<tunnel.centerline.size();++j) {
                    const auto a=tunnel.centerline[j-1].center,b=tunnel.centerline[j].center;
                    const auto x=b.x-a.x,y=b.y-a.y,z=b.z-a.z;
                    length+=std::sqrt(x*x+y*y+z*z);
                }
            }
            previous=std::move(region);
        }
        const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout<<"{\"regions\":"<<count<<",\"caverns\":"<<caverns
                 <<",\"tunnels_including_stitches\":"<<tunnels<<",\"profiles\":"<<profiles
                 <<",\"local_tunnel_meters\":"<<length<<",\"work_used\":"<<work
                 <<",\"support_points_total\":"<<points<<",\"elapsed_seconds\":"<<seconds<<"}\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
