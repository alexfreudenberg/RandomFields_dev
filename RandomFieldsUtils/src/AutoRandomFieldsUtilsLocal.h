#ifndef auto_rfutils_local_h
#define auto_rfutils_local_h 1



// Reihenfolge nie aendern!!
typedef enum la_modes {LA_AUTO, LA_INTERN, LA_R, LA_GPU} la_modes;
#define LA_LAST LA_GPU
// Reihenfolge nie aendern!!
typedef enum pivot_modes {PIVOT_NONE, PIVOT_AUTO, PIVOT_DO, PIVOT_IDX,
			  PIVOT_UNDEFINED} pivot_modes;
#define PIVOT_LAST PIVOT_UNDEFINED

#define PIVOTSPARSE_MMD 1 // for spam
#define PIVOTSPARSE_RCM 2 // for spam


extern const char *LA_NAMES[LA_LAST + 1], *PIVOT_NAMES[PIVOT_LAST + 1];



#endif
  
