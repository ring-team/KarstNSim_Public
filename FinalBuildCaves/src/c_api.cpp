#include "fbs/caves.h"
#include "fbs/caves.hpp"
#include "transport_internal.hpp"
#include <algorithm>
#include <cstring>
#include <memory>
#include <new>
#include <string>

struct fbs_caves_result { fbs::caves::Region region; std::string json; };
struct fbs_caves_buffer { std::string text; };
namespace {
void message(char* target, size_t capacity, const char* text) noexcept {
    if (!target || !capacity) return;
    const auto size=std::min(capacity-1,std::strlen(text));
    std::memcpy(target,text,size);target[size]='\0';
}
fbs_caves_status status(fbs::caves::ErrorCode code) noexcept {
    using E=fbs::caves::ErrorCode;
    switch(code) {
    case E::invalid_input:return FBS_CAVES_INVALID_INPUT;
    case E::incompatible_version:return FBS_CAVES_VERSION;
    case E::cancelled:return FBS_CAVES_CANCELLED;
    case E::budget:return FBS_CAVES_BUDGET;
    case E::no_route:return FBS_CAVES_NO_ROUTE;
    case E::constraint:return FBS_CAVES_CONSTRAINT;
    case E::incompatible_boundary:return FBS_CAVES_BOUNDARY;
    case E::internal:return FBS_CAVES_INTERNAL;
    }
    return FBS_CAVES_INTERNAL;
}
template<class Function> fbs_caves_status guarded(Function fn,char* error,size_t cap) noexcept {
    try { fn();message(error,cap,"");return FBS_CAVES_OK; }
    catch(const fbs::caves::Error& e) { message(error,cap,e.what());return status(e.code); }
    catch(const std::bad_alloc&) { message(error,cap,"Allocation failed");return FBS_CAVES_MEMORY; }
    catch(const std::exception& e) { message(error,cap,e.what());return FBS_CAVES_INTERNAL; }
    catch(...) { message(error,cap,"Unexpected native exception");return FBS_CAVES_INTERNAL; }
}
}
extern "C" {
fbs_caves_options fbs_caves_options_default(void) {
    return {sizeof(fbs_caves_options),FBS_CAVES_ABI_VERSION,100000000,50000,2000000,nullptr,nullptr};
}
const char* fbs_caves_version(void) { return "1.0.0"; }
fbs_caves_status fbs_caves_generate_json(const char* request,size_t size,const fbs_caves_options* supplied,
    fbs_caves_result** out,char* error,size_t cap) {
    return guarded([&] {
        using namespace fbs::caves;
        if (!request || !size || size>16u*1024u*1024u || !out) throw Error(ErrorCode::invalid_input,"Invalid request buffer or output pointer");
        if(supplied && (supplied->struct_size!=sizeof(fbs_caves_options)||supplied->abi_version!=FBS_CAVES_ABI_VERSION)) throw Error(ErrorCode::incompatible_version,"Unsupported options ABI");
        const auto cfg=supplied?*supplied:fbs_caves_options_default();
        Options options;options.work_limit=cfg.work_limit;options.max_points=cfg.max_points;options.max_edges=cfg.max_edges;
        if(cfg.is_cancelled) options.cancelled=[cfg] { return cfg.is_cancelled(cfg.user)!=0; };
        auto result=std::make_unique<fbs_caves_result>();
        result->region=generate(request_from_json(std::string(request,size)),options);
        if(options.cancelled && options.cancelled()) throw Error(ErrorCode::cancelled,"Generation cancelled");
        result->json=detail::serialize_validated_region(result->region);
        if(options.cancelled && options.cancelled()) throw Error(ErrorCode::cancelled,"Generation cancelled");
        *out=result.release();
    },error,cap);
}
const char* fbs_caves_result_json(const fbs_caves_result* result,size_t* size) {
    if(size) *size=result?result->json.size():0;
    return result?result->json.c_str():nullptr;
}
void fbs_caves_result_destroy(fbs_caves_result* result) { const std::unique_ptr<fbs_caves_result> owner(result); }
fbs_caves_status fbs_caves_export_moocow_csv(const fbs_caves_result* const* results,size_t count,
    const char* map_id,const char* map_name,fbs_caves_buffer** out,char* error,size_t cap) {
    return guarded([&] {
        using namespace fbs::caves;
        if(!results||!count||count>1024||!map_id||!out) throw Error(ErrorCode::invalid_input,"Invalid export arguments");
        std::vector<Region> regions;regions.reserve(count);
        for(size_t i=0;i<count;++i) {
            if(!results[i]) throw Error(ErrorCode::invalid_input,"Null region handle");
            regions.push_back(results[i]->region);
        }
        auto result=std::make_unique<fbs_caves_buffer>();
        result->text=to_moocow_csv(regions,map_id,map_name?map_name:"Generated caves");
        *out=result.release();
    },error,cap);
}
const char* fbs_caves_buffer_data(const fbs_caves_buffer* buffer,size_t* size) {
    if(size) *size=buffer?buffer->text.size():0;
    return buffer?buffer->text.c_str():nullptr;
}
void fbs_caves_buffer_destroy(fbs_caves_buffer* buffer) { const std::unique_ptr<fbs_caves_buffer> owner(buffer); }
}
