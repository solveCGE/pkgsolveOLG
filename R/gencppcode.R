# takes a data list as input
gencppcode <- function(datalist,prefix="") {
  
  # parameters
  namein   = "datain";  # name of datainput variable name in C++
  instname = "dataOLG"; # name of the class instance
  
  # start of code
  number_of_list_elements <- length(datalist);
  names_of_list_elements  <- names(datalist);
  
  # get type of every variable
  vars = data.frame("varname"=names_of_list_elements, "type"=rep("",number_of_list_elements), "type_fixed"=rep("",number_of_list_elements));
  
  for (counter in 1:number_of_list_elements) {

    if (is.matrix(datalist[[counter]])) {
      
      if (dim(datalist[[counter]])[1] == 1) {
        typemat = "rowvec"; # row vector
      } else if (dim(datalist[[counter]])[2] == 1) {
        typemat = "colvec"; # row vector
      } else {
        typemat = "mat"; # matrix
      }
      
      dims    = dim(datalist[[counter]]);
      dimtext = paste0("::fixed<")
      if (typemat == "rowvec") dimtext = paste0(dimtext,dims[2]);
      if (typemat == "colvec") {dimtext = paste0(dimtext,dims[1]); colvectext = paste0(dimtext,">");}
      if (typemat == "mat")    {dimtext = paste0(dimtext,dims[1],",",dims[2]); tendtext = datalist[["tend"]]; nagtext = datalist[["nag"]];}
      dimtext = paste0(dimtext,">");

      typetext       = paste0("arma::",typemat);
      typetext_fixed = paste0("arma::",typemat,dimtext);

    } else if (is.double(datalist[[counter]])) {
      
      # double
      typemat  = NULL;
      typetext = "double";
      typetext_fixed = typetext;

    } else if (is.integer(datalist[[counter]])) {
      
      # int
      typemat  = NULL;
      typetext = "int";
      typetext_fixed = typetext;

    } else if (is.logical(datalist[[counter]])) {
      
      # bool
      typemat  = NULL;
      typetext = "bool";
      typetext_fixed = typetext;

    } else {
      stop("Unknown data type in data list!");
    }
    
    vars[counter,"type"]       = typetext;
    vars[counter,"type_fixed"] = typetext_fixed;
 
  }
  
  varnamemax    = max(nchar(vars[,"varname"]))+3;
  typemax       = max(nchar(vars[,"type"]))+3;
  type_fixedmax = max(nchar(vars[,"type_fixed"]))+3;
  
  varname_fill        = sapply(vars[,"varname"], function(x) {paste0(x, paste0(rep(" ",varnamemax-nchar(x)),collapse=""))});
  type_fill           = sapply(vars[,"type"], function(x) {paste0(x, paste0(rep(" ",typemax-nchar(x)),collapse=""))});
  type_fixed_fill     = sapply(vars[,"type_fixed"], function(x) {paste0(x, paste0(rep(" ",type_fixedmax-nchar(x)),collapse=""))});
  type_fill_ref       = sapply(vars[,"type"], function(x) {paste0(x,"&", paste0(rep(" ",typemax-nchar(x)+1),collapse=""))});
  type_fixed_fill_ref = sapply(vars[,"type_fixed"], function(x) {paste0(x,"&", paste0(rep(" ",type_fixedmax-nchar(x)+1),collapse=""))});
  
  txt_init           = paste0(type_fill,vars[,"varname"],";");
  txt_init_fixed     = paste0(type_fixed_fill,vars[,"varname"],";");
  txt_init_ext       = paste0("extern ",txt_init);
  txt_init_fixed_ext = paste0("extern ",txt_init_fixed);
  
  txt_import         = paste0(varname_fill,"= as<", vars[,"type"],">(",namein,"[\"",vars[,"varname"],"\"]);")
  
  txt_export             = paste0("out[",1:number_of_list_elements-1,"] = ",sapply(vars[,"varname"], function(x) {tmppos = which(vars[,"varname"]==x); if (substr(vars[tmppos,"type"],1,4)=="arma") paste0(vars[tmppos,"type"],"(",x,")") else x}),";\nnames[",1:number_of_list_elements-1,"] = \"", vars[,"varname"], "\";");
  txt_export_class       = paste0("out[",1:number_of_list_elements-1,"] = ",sapply(vars[,"varname"], function(x) {tmppos = which(vars[,"varname"]==x); if (substr(vars[tmppos,"type"],1,4)=="arma") paste0(vars[tmppos,"type"],"(",instname,".",x,")") else paste0(instname,".",x)}),";\nnames[",1:number_of_list_elements-1,"] = \"", vars[,"varname"], "\";");
  txt_export_class_fixed = paste0("out[",1:number_of_list_elements-1,"] = ",sapply(vars[,"varname"], function(x) {tmppos = which(vars[,"varname"]==x); if (substr(vars[tmppos,"type"],1,4)=="arma") paste0(vars[tmppos,"type"],"(",instname,"->",x,")") else paste0(instname,"->",x)}),";\nnames[",1:number_of_list_elements-1,"] = \"", vars[,"varname"], "\";");
  
  txt_unpack         = paste0(type_fill_ref,varname_fill,"= ",instname,".",vars[,"varname"],";");
  txt_unpack_fixed   = paste0(type_fixed_fill_ref,varname_fill,"= ",instname,"->",vars[,"varname"],";");
  
  # write h-files
  code_class_init             = paste0("// copy this code to auto/class_init.h and rebuild package.\n// Time stamp: ", Sys.time(), "\n");
  code_classfixed_init        = paste0("// copy this code to auto/classfixed_init.h and rebuild package.\n// Time stamp: ", Sys.time(), "\n");
  code_global_init            = paste0("// copy this code to auto/global_init.h and rebuild package.\n// Time stamp: ", Sys.time(), "\n");
  code_globalfixed_init       = paste0("// copy this code to auto/globalfixed_init.h and rebuild package.\n// Time stamp: ", Sys.time(), "\n");
  code_global_init_ext        = paste0("// copy this code to auto/global_init_ext.h and rebuild package.\n// Time stamp: ", Sys.time(), "\n");
  code_globalfixed_init_ext   = paste0("// copy this code to auto/globalfixed_init_ext.h and rebuild package.\n// Time stamp: ", Sys.time(), "\n");
  
  code_class_import           = paste0("// copy this code to auto/class_import.h and rebuild package.\n// Time stamp: ", Sys.time(), "\n");
  code_global_import          = paste0("// copy this code to auto/global_import.h and rebuild package.\n// Time stamp: ", Sys.time(), "\n");
  
  code_class_export           = paste0("// copy this code to auto/class_export.h and rebuild package.\n// Time stamp: ", Sys.time(), "\n");
  code_classfixed_export      = paste0("// copy this code to auto/classfixed_export.h and rebuild package.\n// Time stamp: ", Sys.time(), "\n");
  code_global_export          = paste0("// copy this code to auto/global_export.h and rebuild package.\n// Time stamp: ", Sys.time(), "\n");
  
  code_class_unpack           = paste0("// copy this code to auto/class_unpack.h and rebuild package.\n// Time stamp: ", Sys.time(), "\n");
  code_classfixed_unpack      = paste("// copy this code to auto/classfixed_unpack.h and rebuild package.\n// Time stamp: ", Sys.time(), "\n");

  ### init files
  temp_define                 = paste0("\n#define ___COLVECFIXED___ arma::colvec\n#define ___TEND___ -1\n#define ___NAG___ -1");
  temp_fixed_define           = paste0("\n#define ___COLVECFIXED___ arma::colvec",colvectext,"\n#define ___TEND___ ", tendtext,"\n#define ___NAG___ ",nagtext);
  
  code_class_init             = c(code_class_init,"class DataOLG {\n\npublic:\nDataOLG(const List&);\n",txt_init,"\n};\n",temp_define);
  code_classfixed_init        = c(code_classfixed_init,"class DataOLG {\n\nprivate:\n~DataOLG(){} // private destructor: only allows dynamically allocated instances of DataOLG\n\npublic:\nDataOLG(const List&);\nfriend void deleteDataOLG(DataOLG*);\n",txt_init_fixed,"\n};\n",paste0("void deleteDataOLG(DataOLG* ",instname,");\n"),temp_fixed_define);
  code_global_init            = c(code_global_init,txt_init,temp_define);
  code_global_init_ext        = c(code_global_init_ext,txt_init_ext,temp_define);
  code_globalfixed_init       = c(code_globalfixed_init,txt_init_fixed,temp_fixed_define);
  code_globalfixed_init_ext   = c(code_globalfixed_init_ext,txt_init_fixed_ext,temp_fixed_define);

  ### import files
  code_class_import           = c(code_class_import,paste0("DataOLG::DataOLG(const List& ",namein,") {\n"),txt_import,"\n}");
  code_global_import          = c(code_global_import,paste0("void import2cpp(const List& ",namein,") {\n"),txt_import,"\n}");
  
  ### export files
  code_class_export           = c(code_class_export,paste0("List export2R(DataOLG ",instname,") {\n"),paste0("List out(",number_of_list_elements,");\nCharacterVector names(",number_of_list_elements,");\n"), txt_export_class, "\nout.attr(\"names\") = names;\nreturn out;\n\n}");
  code_classfixed_export      = c(code_classfixed_export,paste0("List export2R(DataOLG* ",instname,") {\n"),paste0("List out(",number_of_list_elements,");\nCharacterVector names(",number_of_list_elements,");\n"), txt_export_class_fixed, "\nout.attr(\"names\") = names;\nreturn out;\n\n}");
  code_global_export          = c(code_global_export,"List export2R() {\n",paste0("List out(",number_of_list_elements,");\nCharacterVector names(",number_of_list_elements,");\n"), txt_export, "\nout.attr(\"names\") = names;\nreturn out;\n\n}");                                
  
  ### unpack files
  code_class_unpack           = c(code_class_unpack,txt_unpack);
  code_classfixed_unpack      = c(code_classfixed_unpack,txt_unpack_fixed);
  
  writeLines(code_class_init, paste0(prefix,"class_init.h"));
  writeLines(code_classfixed_init, paste0(prefix,"classfixed_init.h"));
  writeLines(code_global_init, paste0(prefix,"global_init.h"));
  writeLines(code_global_init_ext, paste0(prefix,"global_init_ext.h"));
  writeLines(code_globalfixed_init, paste0(prefix,"globalfixed_init.h"));
  writeLines(code_globalfixed_init_ext, paste0(prefix,"globalfixed_init_ext.h"));
  
  writeLines(code_class_import, paste0(prefix,"class_import.h"));
  writeLines(code_global_import, paste0(prefix,"global_import.h"));
  
  writeLines(code_class_export, paste0(prefix,"class_export.h"));
  writeLines(code_classfixed_export, paste0(prefix,"classfixed_export.h"));  
  writeLines(code_global_export, paste0(prefix,"global_export.h"));
  
  writeLines(code_class_unpack, paste0(prefix,"class_unpack.h"));
  writeLines(code_classfixed_unpack, paste0(prefix,"classfixed_unpack.h"));
  
}
