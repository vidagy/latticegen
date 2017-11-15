
#ifndef LATTICEGEN_MUTE_H
#define LATTICEGEN_MUTE_H

#ifdef __GNUG__

#define LATTICEGEN_MUTE_BEGIN _Pragma("GCC diagnostic push")
#define LATTICEGEN_MUTE_EIGEN _Pragma("GCC diagnostic ignored \"-Wint-in-bool-context\"")
#define LATTICEGEN_MUTE_END _Pragma("GCC diagnostic pop")

#endif

#endif //LATTICEGEN_MUTE_H
