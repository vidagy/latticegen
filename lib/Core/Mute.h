
#ifndef LATTICEGEN_MUTE_H
#define LATTICEGEN_MUTE_H

#if defined(__clang__)
#define LATTICEGEN_MUTE_BEGIN _Pragma("GCC diagnostic push")
#define LATTICEGEN_MUTE_EIGEN
#define LATTICEGEN_MUTE_UNDEFINED_VAR_TEMPLATE _Pragma("GCC diagnostic ignored \"-Wundefined-var-template\"")
#define LATTICEGEN_MUTE_UNUSED_VAR _Pragma("GCC diagnostic ignored \"-Wunused-parameter\"")
#define LATTICEGEN_MUTE_END _Pragma("GCC diagnostic pop")

#elif defined(__GNUG__)
#define LATTICEGEN_MUTE_BEGIN _Pragma("GCC diagnostic push")
#define LATTICEGEN_MUTE_EIGEN _Pragma("GCC diagnostic ignored \"-Wint-in-bool-context\"")
#define LATTICEGEN_MUTE_UNDEFINED_VAR_TEMPLATE
#define LATTICEGEN_MUTE_UNUSED_VAR _Pragma("GCC diagnostic ignored \"-Wunused-parameter\"")
#define LATTICEGEN_MUTE_END _Pragma("GCC diagnostic pop")
#endif


#endif //LATTICEGEN_MUTE_H
