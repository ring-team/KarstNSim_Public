#include "fbs/caves.h"
#include <cstring>
#include <iostream>
#include <stdexcept>
static void require(bool ok,const char* text) {if(!ok) throw std::runtime_error(text);}
static int cancel(void*) {return 1;}
int main() {
    try {
        auto cfg=fbs_caves_options_default();char error[256];
        fbs_caves_result* result=nullptr;
        require(fbs_caves_generate_json("{",1,&cfg,&result,error,sizeof(error))==FBS_CAVES_INVALID_INPUT,"malformed JSON accepted");
        require(result==nullptr && error[0],"error modified result or omitted message");
        cfg.abi_version=999;
        require(fbs_caves_generate_json("{}",2,&cfg,&result,error,sizeof(error))==FBS_CAVES_VERSION,"ABI accepted");
        cfg=fbs_caves_options_default();cfg.is_cancelled=cancel;
        require(fbs_caves_generate_json("{}",2,&cfg,&result,error,sizeof(error))==FBS_CAVES_CANCELLED,"cancellation failed");
        cfg=fbs_caves_options_default();
        const char* request="{\"world_seed\":17,\"loop_count\":0}";
        const auto status=fbs_caves_generate_json(request,std::strlen(request),&cfg,&result,error,sizeof(error));
        if(status!=FBS_CAVES_OK) throw std::runtime_error(error);
        size_t size=0;const auto* data=fbs_caves_result_json(result,&size);
        require(data && size && std::strstr(data,"caverns"),"missing result JSON");
        auto* previous=result;
        require(fbs_caves_generate_json("{",1,&cfg,&result,error,1)==FBS_CAVES_INVALID_INPUT,"second invalid call accepted");
        require(result==previous && error[0]=='\0',"failed generation destroyed accepted result");
        fbs_caves_buffer* csv=nullptr;const fbs_caves_result* inputs[]={result};
        require(fbs_caves_export_moocow_csv(inputs,1,"f0c148bd-517c-4c30-aa66-dceff189e87a","C ABI cave",&csv,error,sizeof(error))==FBS_CAVES_OK,error);
        data=fbs_caves_buffer_data(csv,&size);require(data && size && std::strstr(data,"=== TUNNELS ==="),"missing CSV");
        fbs_caves_buffer_destroy(csv);fbs_caves_result_destroy(result);
        fbs_caves_result_destroy(nullptr);fbs_caves_buffer_destroy(nullptr);
        require(fbs_caves_result_json(nullptr,&size)==nullptr && size==0,"null query");
        std::cout<<"C ABI generation, lifetime, export and transactional errors passed\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
