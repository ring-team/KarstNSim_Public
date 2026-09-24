#include "fbs/caves.hpp"
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
using namespace fbs::caves;
static void require(bool ok,const char* what) { if(!ok) throw std::runtime_error(what); }
static void rejected(const std::string& input,ErrorCode code=ErrorCode::invalid_input) {
    try { (void)request_from_json(input); }
    catch(const Error& e) {require(e.code==code,"wrong JSON rejection code");return;}
    throw std::runtime_error("invalid request accepted");
}
int main() {
    try {
        Request r;r.world_seed=std::numeric_limits<std::uint64_t>::max();r.key={-123,456,-789};
        r.style_id=std::numeric_limits<std::uint64_t>::max();
        const auto text=request_to_json(r);const auto copy=request_from_json(text);
        require(request_to_json(copy)==text,"request roundtrip changed bytes");
        require(copy.world_seed==r.world_seed && copy.key.x==-123,"integer precision lost");
        require(request_from_json("{}").chamber_count==8,"defaults changed");
        rejected("{\"world_seed\":-1}");rejected("{\"world_seed\":1.5}");
        rejected("{\"world_seed\":18446744073709551616}");
        rejected("{\"neighbors\":4294967296}");rejected("{\"shelves\":1}");
        rejected("{\"key\":[1.5,0,0]}");rejected("{\"key\":[0,0]}");
        rejected("{\"size\":[1,1,1e999]}");rejected("{\"neighbours\":20}");
        rejected("{\"world_seed\":1,\"world_seed\":2}");
        rejected("{\"version\":2}",ErrorCode::incompatible_version);
        rejected("{\"authored_ports\":[{\"state\":\"unknown\"}]}");
        rejected(std::string(65,'[')+"0"+std::string(65,']'));
        rejected(std::string(16u*1024u*1024u+1,' '));
        r.sampling_radius=std::numeric_limits<double>::quiet_NaN();
        bool nonfinite=false;
        try { (void)request_to_json(r); } catch(const Error& e) {nonfinite=e.code==ErrorCode::invalid_input;}
        require(nonfinite,"non-finite request serialized as null");
        r.sampling_radius=0.08;r.authored_caverns.emplace_back();
        r.authored_caverns.back().name=std::string(1,static_cast<char>(0xff));
        bool utf8=false;
        try { (void)request_to_json(r); } catch(const Error& e) {utf8=e.code==ErrorCode::invalid_input;}
        require(utf8,"invalid UTF-8 escaped the public error contract");
        std::cout << "Transport roundtrip, integer precision, bounded and strict parsing passed\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
