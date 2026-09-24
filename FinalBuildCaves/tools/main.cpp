#include "fbs/caves.hpp"
#include <nlohmann/json.hpp>
#include <csignal>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
namespace fs=std::filesystem;
using namespace fbs::caves;
namespace {
volatile std::sig_atomic_t cancelled=0;
void interrupt(int) {cancelled=1;}
void check_cancelled() {
    if(cancelled) throw Error(ErrorCode::cancelled,"Generation cancelled");
}
std::string load(const fs::path& path) {
    const auto size=fs::file_size(path);
    if(size>16u*1024u*1024u) throw Error(ErrorCode::invalid_input,"Request exceeds 16 MiB");
    std::ifstream input(path,std::ios::binary);
    if(!input) throw std::runtime_error("Cannot read "+path.string());
    std::string text(static_cast<std::size_t>(size),'\0');
    input.read(text.data(),static_cast<std::streamsize>(text.size()));
    if(!input && !text.empty()) throw std::runtime_error("Incomplete read: "+path.string());
    return text;
}
void save(const fs::path& path,const std::string& text) {
    std::ofstream out(path,std::ios::binary|std::ios::trunc);
    if(!out) throw std::runtime_error("Cannot write "+path.string());
    out.write(text.data(),static_cast<std::streamsize>(text.size()));out.close();
    if(!out) throw std::runtime_error("Incomplete write: "+path.string());
}
std::uint64_t number(const std::string& text) {
    if(text.empty() || text.find_first_not_of("0123456789")!=std::string::npos) throw std::runtime_error("Expected an unsigned integer");
    return std::stoull(text);
}
void usage() {
    std::cout<<"fbs-caves --request region.json [--request neighbor.json] --out NEW_DIRECTORY\n"
        "          --map-id UUID [--name NAME] [--work-limit N] [--max-points N] [--max-edges N]\n"
        "fbs-caves --default-request\n"
        "Generates region JSON, current-format MOOCoW CSV and a completion manifest.\n"
        "Work, point and edge limits apply independently to each region request.\n"
        "Existing output directories are rejected. No database is opened or changed.\n";
}
}
int main(int argc,char** argv) {
    fs::path output;bool owns_output=false;
    try {
        std::vector<fs::path> requests;
        std::string map_id,map_name="Generated caves";Options options;
        for(int i=1;i<argc;++i) {
            const std::string arg=argv[i];
            if(arg=="--help") {usage();return 0;}
            if(arg=="--default-request") {std::cout<<request_to_json(Request{});return 0;}
            if(i+1==argc) throw std::runtime_error("Missing value for "+arg);
            const std::string value=argv[++i];
            if(arg=="--request") requests.emplace_back(value);
            else if(arg=="--out") output=value;
            else if(arg=="--map-id") map_id=value;
            else if(arg=="--name") map_name=value;
            else if(arg=="--work-limit") options.work_limit=number(value);
            else if(arg=="--max-points") options.max_points=number(value);
            else if(arg=="--max-edges") options.max_edges=number(value);
            else throw std::runtime_error("Unknown argument: "+arg);
        }
        if(requests.empty()||requests.size()>1024||output.empty()||map_id.empty()) {usage();return 2;}
        if(fs::exists(output)) throw std::runtime_error("Output already exists: "+output.string());
        std::signal(SIGINT,interrupt);std::signal(SIGTERM,interrupt);
        options.cancelled=[] {return cancelled!=0;};
        std::vector<Region> regions;regions.reserve(requests.size());
        for(const auto& path:requests) regions.push_back(generate(request_from_json(load(path)),options));
        // Validate the complete export before creating any output directory.
        const auto csv=to_moocow_csv(regions,map_id,map_name);
        check_cancelled();
        if(!fs::create_directory(output)) throw std::runtime_error("Cannot reserve new output directory");
        owns_output=true;
        nlohmann::json manifest={{"schema_version",1},{"algorithm",algorithm_version},
            {"map_id",map_id},{"units","meters"},{"vertical_axis","Z"},
            {"csv_format","existing MOOCoW sectioned CSV"},{"database_migration_required",false},
            {"regions",nlohmann::json::array()}};
        for(std::size_t i=0;i<regions.size();++i) {
            check_cancelled();
            const auto name="region-"+std::to_string(i)+".json";
            save(output/name,to_json(regions[i]));
            manifest["regions"].push_back({{"file",name},{"id",regions[i].id},{"settings_id",regions[i].settings_id}});
        }
        save(output/"dataset.csv",csv);
        // Consumers accept a directory only after this last file is present.
        check_cancelled();
        manifest["complete"]=true;save(output/"manifest.json.tmp",manifest.dump(2)+"\n");
        check_cancelled();
        fs::rename(output/"manifest.json.tmp",output/"manifest.json");
        std::cout<<fs::absolute(output).string()<<'\n';return 0;
    } catch(const std::exception& e) {
        if(owns_output) {std::error_code ignored;fs::remove_all(output,ignored);}
        std::cerr<<"fbs-caves: "<<e.what()<<'\n';return 1;
    }
}
