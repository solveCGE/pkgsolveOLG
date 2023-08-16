#ifdef ___USEGLOBAL___
  // do nothing
#else
  #ifdef ___USEFIXED___
    #include "auto/classfixed_unpack.h" // unpack dataOLG using references to elements
  #else
    #include "auto/class_unpack.h" // unpack dataOLG using references to elements
  #endif
#endif
