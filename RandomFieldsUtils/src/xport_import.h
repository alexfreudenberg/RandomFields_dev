
#ifndef RandomFieldsUtilsxport_H
#define RandomFieldsUtilsxport_H 1

extern int CORES;
extern int PL;

bool parallel();


extern utilsparam GLOBAL;

#define prefixN 2
extern const char * prefixlist[prefixN], **all[prefixN];
extern int allN[prefixN];

#endif
