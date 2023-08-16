#ifdef ___USEGLOBAL___
  #ifdef ___USEFIXED___
    //Rcout << "CASE: useglobal & usefixed" << endl;
  #else
    //Rcout << "CASE: useglobal & !usefixed" << endl;
  #endif
  import2cpp(dataOLGin);
#else
  #ifdef ___USEFIXED___
    //Rcout << "CASE: !useglobal & usefixed" << endl;
    DataOLG* dataOLG = new DataOLG(dataOLGin);
  #else
    //Rcout << "CASE: !useglobal & !usefixed" << endl;
    DataOLG dataOLG = DataOLG(dataOLGin);
  #endif
#endif
