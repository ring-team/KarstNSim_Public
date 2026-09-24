#include <fbs/caves.hpp>
#include <iostream>
int main() {
    try {
        fbs::caves::Request request;request.world_seed=42;request.loop_count=0;
        const auto region=fbs::caves::generate(request);
        fbs::caves::validate(region);
        std::cout<<region.caverns.size()<<" caverns, "<<region.tunnels.size()<<" tunnels\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
