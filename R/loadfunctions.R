# formatted reporting
report <- function(reporttext,reportcalc) {
  cursorstart     = 45;
  
  countlinebreaks = 0;
  for (i in 1:nchar(reporttext)) {
    if (substring(reporttext,1,i) != "\n") break;
    countlinebreaks = countlinebreaks + 1;
  }
  
  cursorstart     = max(nchar(reporttext)-countlinebreaks,cursorstart)+3*(reportcalc>=0)+2*(reportcalc<0);
  
  charfill = paste(rep(" ",cursorstart-nchar(reporttext)+countlinebreaks),collapse="");
  
  cat(paste(reporttext,charfill,reportcalc,"\n",sep=""));
  
}

# initializes variables
initvars <- function(timedep,agedep,varnames,float=T) {
  
  varlist = character(0);
  
  for (varname in varnames) {
    
    if (timedep && agedep) {
      assign(paste0(varname,"v"), zerosmat(nag,tend), envir = .GlobalEnv);
      assign(paste0(varname,"z"), zerosmat(nag,ncoh), envir = .GlobalEnv);
      assign(paste0(varname,"v0"), zeroscol(nag), envir = .GlobalEnv);
      varlist = c(varlist,paste0(varname,c("v","z","v0")));
    }
    
    if (!timedep && agedep) {
      assign(varname, zeroscol(nag), envir = .GlobalEnv);
      varlist = c(varlist,varname);
    }
    
    if (timedep && !agedep) {
      assign(varname, zerosrow(tend), envir = .GlobalEnv);
      assign(paste0(varname,"0"), 0, envir = .GlobalEnv);
      varlist = c(varlist,paste0(varname,c("","0")));
    }
    
    if (!timedep && !agedep && float) {
      assign(varname, 0, envir = .GlobalEnv);
      varlist = c(varlist,varname);
    }
    
    if (!timedep && !agedep && !float) {
      assign(varname, 0L, envir = .GlobalEnv);
      varlist = c(varlist,varname);
    }
  }
  return(invisible(varlist));
}

# fill time-dependent variables with calibration values
fillvars <- function(agedep,varnames) {
  
  for (varname in varnames) {
    
    if (agedep) {
      if (!exists(paste0(varname,"z")) || !exists(paste0(varname,"v")) || !exists(paste0(varname,"v0"))) stop("Variable ", paste0(varname), " is not properly initialized!");
      assign(paste0(varname,"v"), kronecker(eval(parse(text=paste0(varname,"v0"))),onesrow(tend)), envir = .GlobalEnv);
      assign(paste0(varname,"z"), kronecker(eval(parse(text=paste0(varname,"v0"))),onesrow(ncoh)), envir = .GlobalEnv);
    }
    
    if (!agedep) {
      if (!exists(varname) || !exists(paste0(varname,"0"))) stop("Variable ", paste0(varname), " is not properly initialized!");
      assign(varname, eval(parse(text=paste0(varname,"0")))*onesrow(tend), envir = .GlobalEnv);
    }
  }
}

zeroscol <- function(dim1) {
  return(matrix(rep(0,dim1),nrow=dim1,ncol=1));
}

onescol <- function(dim1) {
  return(matrix(rep(1,dim1),nrow=dim1,ncol=1));
}

zerosrow <- function(dim1) {
  return(matrix(rep(0,dim1),ncol=dim1,nrow=1));
}

onesrow <- function(dim1) {
  return(matrix(rep(1,dim1),ncol=dim1,nrow=1));
}

zerosmat <- function(dim1,dim2) {
  return(matrix(rep(0,dim1*dim2),nrow=dim1,ncol=dim2));
}

onesmat <- function(dim1,dim2) {
  return(matrix(rep(1,dim1*dim2),nrow=dim1,ncol=dim2));
}

format_dec <- function(x, k) trimws(format(round(x, k), nsmall=k));
format_sci <- function(numb) formatC(numb, format = "e", digits = 2);

# R-wrapper for C++ solveOLG_()
solveOLG <- function(starttime = 1, maxiter = 200, tol = 1e-4, damping_budget = 1.0, damping_assets = 1.0, damping_ab = 1.0, damping_r = 0.5, damping_new_assets = 0.7, nthrdsin = nthrds, dataOLGnamesin = dataOLGnames) {
  
  cat("\nRunning Tatonnement Algorithm for Transition:\n\n");

  pkgdims = fixeddim();
  
  if ((pkgdims[1]>0) && (pkgdims[1]!=tend))  stop("'pkgsolveOLG' was compiled with fixed dimensions (tend=",pkgdims[1],") that are incompatible with chosen 'tend'. 'pkgsolveOLG' has to be recompiled.");
  if ((pkgdims[2]>0) && (pkgdims[2]!=nag))   stop("'pkgsolveOLG' was compiled with fixed dimensions (nag=",pkgdims[2],") that are incompatible with chosen 'nag'. 'pkgsolveOLG' has to be recompiled.");
  
  dataOLG = collectglobals(dataOLGnamesin);

  dataOLG = solveOLG_(starttime, maxiter, tol, damping_budget, damping_assets, damping_ab, damping_r, damping_new_assets, nthrdsin, dataOLG);
  
  makeelementsglobal(dataOLG, env=parent.frame());
  
  cat(paste0("CHECK SOLUTION:\t\t", max(abs(edy)+abs(edl)+abs(edg)+abs(eda)+abs(ediv)+abs(edab))),"\n");
  
} 

# collect global variables in a list
collectglobals <- function(globalnamelist) {
  number_of_list_elements <- length(globalnamelist);
  
  outlist = list();
  
  for (counter in 1:number_of_list_elements) {
    outlist[[globalnamelist[counter]]] = get(globalnamelist[counter]);
  }
  return(outlist);
}

# write list items back to global workspace
makeelementsglobal <- function(inputlist, env = environment()) {
  number_of_list_elements <- length(inputlist);
  names_of_list_elements  <- names(inputlist);
  
  for (counter in 1:number_of_list_elements) {
    assign(names_of_list_elements[counter],inputlist[[counter]], envir = env);
  }
}

