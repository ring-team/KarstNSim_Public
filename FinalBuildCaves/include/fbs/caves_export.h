#ifndef FBS_CAVES_EXPORT_H
#define FBS_CAVES_EXPORT_H
#if defined(_WIN32) && defined(FBS_CAVES_SHARED)
# if defined(FBS_CAVES_BUILD)
#  define FBS_CAVES_API __declspec(dllexport)
# else
#  define FBS_CAVES_API __declspec(dllimport)
# endif
#elif defined(__GNUC__) || defined(__clang__)
# define FBS_CAVES_API __attribute__((visibility("default")))
#else
# define FBS_CAVES_API
#endif
#endif
