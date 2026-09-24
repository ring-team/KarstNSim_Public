#ifndef FBS_CAVES_H
#define FBS_CAVES_H
#include <stddef.h>
#include <stdint.h>
#include "caves_export.h"
#ifdef __cplusplus
extern "C" {
#endif
#define FBS_CAVES_ABI_VERSION 1u
typedef enum fbs_caves_status {
 FBS_CAVES_OK=0, FBS_CAVES_INVALID_INPUT=1, FBS_CAVES_VERSION=2,
 FBS_CAVES_CANCELLED=3, FBS_CAVES_BUDGET=4, FBS_CAVES_NO_ROUTE=5,
 FBS_CAVES_CONSTRAINT=6, FBS_CAVES_BOUNDARY=7, FBS_CAVES_INTERNAL=8,
 FBS_CAVES_MEMORY=9
} fbs_caves_status;
typedef struct fbs_caves_options {
 uint32_t struct_size, abi_version;
 uint64_t work_limit, max_points, max_edges;
 /* Synchronous callback. Return nonzero to cancel. Never throw across C ABI. */
 int (*is_cancelled)(void* user);
 void* user;
} fbs_caves_options;
typedef struct fbs_caves_result fbs_caves_result;
typedef struct fbs_caves_buffer fbs_caves_buffer;

FBS_CAVES_API fbs_caves_options fbs_caves_options_default(void);
FBS_CAVES_API const char* fbs_caves_version(void);
/* UTF-8 request JSON, bounded to 16 MiB. out remains unchanged on failure.
   Independent jobs may run concurrently. No filesystem or stdout side effects.
   Error is optional and truncated/NUL-terminated when error_capacity > 0. */
FBS_CAVES_API fbs_caves_status fbs_caves_generate_json(
 const char* request, size_t request_size, const fbs_caves_options* options,
 fbs_caves_result** out, char* error, size_t error_capacity);
/* Read-only borrowed UTF-8 JSON. Valid until result_destroy. */
FBS_CAVES_API const char* fbs_caves_result_json(const fbs_caves_result*, size_t* size);
FBS_CAVES_API void fbs_caves_result_destroy(fbs_caves_result*);
/* New CSV text in MOOCoW's EXISTING format; never writes a file or database.
   At most 1024 regions per export. map_id must be a nonzero UUID string. */
FBS_CAVES_API fbs_caves_status fbs_caves_export_moocow_csv(
 const fbs_caves_result* const* results, size_t count,
 const char* map_id, const char* map_name, fbs_caves_buffer** out,
 char* error, size_t error_capacity);
FBS_CAVES_API const char* fbs_caves_buffer_data(const fbs_caves_buffer*, size_t* size);
FBS_CAVES_API void fbs_caves_buffer_destroy(fbs_caves_buffer*);
#ifdef __cplusplus
}
#endif
#endif
