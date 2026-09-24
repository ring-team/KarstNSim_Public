#include <fbs/caves.h>
#include <stdio.h>
#include <string.h>
int main(void) {
    const char* request="{\"world_seed\":42,\"loop_count\":0}";
    fbs_caves_options options=fbs_caves_options_default();
    fbs_caves_result* result=NULL;
    char error[512];size_t length=0;
    const fbs_caves_status status=fbs_caves_generate_json(request,strlen(request),&options,&result,error,sizeof(error));
    if(status!=FBS_CAVES_OK) {fprintf(stderr,"%s\n",error);return 1;}
    if(!fbs_caves_result_json(result,&length)||length==0) {fbs_caves_result_destroy(result);return 2;}
    printf("%zu bytes of owned generation output\n",length);
    fbs_caves_result_destroy(result);return 0;
}
