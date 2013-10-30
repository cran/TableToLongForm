##----------------------------------------------------------------------
## The code in this .R file is machine generated from the literate
##  program, TableToLongForm.Rnw
## Documentation can be found in the literate description for this
##  program, TableToLongForm.pdf
##----------------------------------------------------------------------
TableToLongForm =
  function(Table, IdentResult = NULL,
           fulloutput = FALSE, diagnostics = FALSE){
    if(is.data.frame(Table)){
      warning("Table supplied is a data.frame. TableToLongForm ",
              "is designed for a character matrix. The data.frame ",
              "is being coerced to a matrix but this may lead to ",
              "unexpected results.")
      as.matrix(Table)
    }
    if(!is.matrix(Table))
      stop("Table argument must be a matrix or a data.frame")
    if(diagnostics != FALSE){
      if(!is.character(diagnostics))
        diagnostics = deparse(substitute(Table))
      assign("TCRunout", envir = TTLFBaseEnv,
             file(paste0(diagnostics, ".TCRunout"), "w"))
      on.exit({
        with(TTLFBaseEnv, close(TCRunout))
        rm("TCRunout", envir = TTLFBaseEnv)
      })
    }
    fullout = ReconsMain(matFull = Table, IdentResult)
    if(fulloutput) fullout else fullout$datafr
  }
IdentMain =
  function(matFull){
    rowNonempty = (1:nrow(matFull))[IdentNonEmpty(matFull, 1)]
    colNonempty = (1:ncol(matFull))[IdentNonEmpty(matFull, 2)]
    rowData = IdentMostCommonBoundary(matFull, 2)
    colData = IdentMostCommonBoundary(matFull, 1)
    TCRsink("CIMCB", rowData, colData)
    rowslist = list(label = rowNonempty[rowNonempty < rowData[1]],
                    data = rowNonempty[(rowNonempty >= rowData[1]) &
                                       (rowNonempty <= rowData[2])])
    colslist = list(label = colNonempty[colNonempty < colData[1]],
                    data = colNonempty[(colNonempty >= colData[1]) &
                                       (colNonempty <= colData[2])])
    TCRsink("CRAC", rowslist, colslist)
    matRowLabel = matFull[rowslist$data, colslist$label, drop = FALSE]
    if(!all(is.na(matRowLabel)) && ncol(matRowLabel) > 1){
      RowLabelNonempty = IdentNonEmpty(matRowLabel, 2)
      if(max(RowLabelNonempty) < ncol(matRowLabel)){
        toshift = (max(RowLabelNonempty) + 1):ncol(matRowLabel)
        colslist$data = c(colslist$label[toshift], colslist$data)
        colslist$label = colslist$label[-toshift]
      }
    }
    IdentResult = list(rows = rowslist, cols = colslist)
    SeqResult = IdentSequence(matFull, IdentResult)
    if(!all(is.na(SeqResult)))
      IdentResult = SeqResult
    IdentResult
  }
IdentNonEmpty =
  function(mat, margin, emptyident = is.na){
    isnonempty = apply(mat, margin, function(x) !all(emptyident(x)))
    which(isnonempty)
  }
IdentPattern =
  function(vec){
    matchvec = match(vec, unique(vec))
    for(i in 1:length(unique(matchvec))){
      repind = unique(diff(which(matchvec == i)))
      if(length(repind) == 0)
        repind = length(vec)
      if(length(repind) == 1)
        break
    }
    curseg = paste0("^(", paste(vec[1:repind], collapse = ""), ")+$")
    if(length(grep(curseg, paste(vec, collapse = ""))) > 0)
      repind else length(vec)
  }
IdentMostCommonBoundary =
  function(matFull, margin){
    isnumber = suppressWarnings(apply(matFull, margin,
      function(x) which(!is.na(as.numeric(x)))))
    nstarts = table(sapply(isnumber,
      function(x) if(length(x) > 0) min(x) else NA))
    nends = table(sapply(isnumber,
      function(x) if(length(x) > 0) max(x) else NA))
    as.numeric(names(c(which.max(nstarts), which.max(rev(nends)))))
  }
IdentSequence =
  function(matFull, IdentResult)
  with(IdentResult, {
    matRowLabel = matFull[rows$data, cols$label]
    if(all(is.na(matRowLabel))){
      cols$label = cols$data[1]
      cols$data = cols$data[-1]
      IdentSequence(matFull, list(rows = rows, cols = cols))
    }
    else{
      matRowLabel = suppressWarnings(as.numeric(matRowLabel))
      if(length(unique(matRowLabel)) > 1 &&
         length(unique(diff(matRowLabel))) == 1)
        list(rows = rows, cols = cols)
      else NA
    }
  })
PareFront =
  function(matRowLabel)
  PareMain(matSub = matRowLabel, plist =
           list(rows = 1:nrow(matRowLabel), cols = 1:ncol(matRowLabel)))
PareCol =
  function(matFull, IdentResult){
    matColLabel = with(IdentResult,
      matFull[rows$label, cols$data, drop = FALSE])
    matData = with(IdentResult,
      matFull[rows$data, cols$data, drop = FALSE])
    colsData = IdentNonEmpty(matData, 2)
    colsLabels = IdentNonEmpty(matColLabel, 2)
    if(length(colsData) == length(colsLabels))
      if(ncol(matData) != length(colsData)){
        IdentResult$cols$data = IdentResult$cols$data[colsData]
        matColLabel = matColLabel[,colsLabels,drop = FALSE]
        matData = matData[,colsData,drop = FALSE]
      }
    for(i in 1:nrow(matColLabel)){
      currow = matColLabel[i,]
      curPattern =
        if(all(is.na(currow))) NA
        else if(any(is.na(currow))) IdentPattern(is.na(currow))
        else IdentPattern(currow)
      if(!is.na(curPattern)){
        nParents = length(currow)/curPattern
        for(j in 1:nParents){
          curcols = 1:curPattern + curPattern * (j - 1)
          cursub = currow[curcols]
          TCRsink("ACP", cursub)
          currow[curcols] = c(cursub[!is.na(cursub)], cursub[is.na(cursub)])
          TCRsink("ACP", currow[curcols])
        }
        matColLabel[i,] = currow
      }
    }
    fullrows = apply(matColLabel, 1, function(x) all(!is.na(x)))
    if(any(diff(fullrows) > 1))
      warning("full rows followed by not full rows!")
    pastestring = ""
    pasterows = which(fullrows)
    for(i in 1:length(pasterows))
      pastestring[i] = paste0("matColLabel[", pasterows[i],
                   ",,drop = FALSE]")
    collapsedlabels =
      eval(parse(text = paste0("paste(",
                   paste(pastestring, collapse = ", "), ")")))
    
    matColLabel = rbind(matColLabel[!fullrows,, drop = FALSE],
      collapsedlabels)
    list(colplist = PareFront(t(matColLabel)), IdentResult = IdentResult)
}
PareMain =
  function(matSub, plist){
    if(length(plist$cols) == 1){
      res = structure(plist$rows, .Names = matSub[plist$rows, plist$cols])
      res = attrLoc(res, cols = plist$col)
      TCRsink("IOOC", plist, res)
    }
    else if(all(is.na(matSub[plist$rows, plist$cols[1]]))){
      plist$cols = plist$cols[-1]
      res = PareMain(matSub, plist)
    }
    else if(length(plist$rows) == 1){
      res = structure(plist$rows,
        .Names = matSub[plist$rows, plist$cols[length(plist$cols)]])
      res = attrLoc(res, cols = plist$cols[length(plist$cols)])
      for(i in (length(plist$cols) - 1):1){
        res = list(res)
        names(res) = matSub[plist$rows, plist$cols[i]]
        res = attrLoc(res, rows = plist$rows, cols = plist$cols[i])
      }
      TCRsink("IOOR", plist, res)
    }
    else if(is.na(matSub[plist$rows[1], plist$cols[1]])){
      warning("cell[1, 1] is empty")
      print(plist)
      res = NA
    }
    else{
      res = PareByEmptyRight(matSub, plist)
      if(any(is.na(res)))
        res = PareByEmptyBelow(matSub, plist)
      for(i in 1:length(res))
        res[[i]] = PareMain(matSub, res[[i]])
      res
    }
    class(res) = "plist"
    res
  }
PareByEmptyRight =
  function(matSub, plist)
  with(plist,
       if(all(is.na(matSub[rows[1], cols[-1]]))){
         emptyrights = apply(matSub[rows, cols[-1], drop = FALSE], 1,
           function(x) all(is.na(x)))
         rowemptyright = rows[emptyrights]
         if(length(rowemptyright) == 1){
           res = list(list(rows = rows[-1], cols = cols))
           names(res) = matSub[rows[1], cols[1]]
           res = attrLoc(res, rows = rows[1], cols = cols[1])
           TCRsink("CSER", res)
         }
         else{
           rowdiff = diff(rowemptyright)
           if(any(rowdiff == 1))
             rowemptyright = rowemptyright[c(rowdiff == 1, FALSE)]
           
           rowstart = pmin(rowemptyright + 1, max(rows))
           rowend = c(pmax(rowemptyright[-1] - 1, min(rows)), max(rows))
           
           res = list()
           for(i in 1:length(rowstart))
             res[i] = list(list(rows = rowstart[i]:rowend[i], cols = cols))
           names(res) = matSub[rowemptyright, cols[1]]
           res = attrLoc(res, rows = rowemptyright, cols = cols[1])
           TCRsink("CMER", res)
         }
         res
       } else NA)
PareByEmptyBelow =
  function(matSub, plist)
  with(plist, {
    emptybelow = is.na(matSub[rows, cols[1]])
    rowstart = rows[!emptybelow]
    rowend = c(rowstart[-1] - 1, max(rows))
    res = list()
    for(i in 1:length(rowstart))
      res[i] = list(list(rows = rowstart[i]:rowend[i], cols = cols[-1]))
    names(res) = matSub[rowstart, cols[1]]
    res = attrLoc(res, rows = rowstart, cols = cols[1])
    TCRsink("PBEB", res)
    res
  })
ReconsMain =
  function(matFull, IdentResult){
    if(is.null(IdentResult))
      IdentResult = IdentMain(matFull)
    matRowLabel = with(IdentResult,
      matFull[rows$data, cols$label, drop = FALSE])
    matRowLabel = matRowLabel[,
      IdentNonEmpty(matRowLabel, 2), drop = FALSE]
    rowplist = PareFront(matRowLabel)
    rowvecs = ReconsRowLabels(rowplist)
    TCRsink("RRL", rowplist, rowvecs[1:4,])
    PareColres = PareCol(matFull, IdentResult)
    colplist = PareColres$colplist
    IdentResult = PareColres$IdentResult
    matData = with(IdentResult,
      matFull[rows$data[unlist(rowplist)], cols$data])
    res = ReconsColLabels(colplist, matData, rowvecs)
    TCRsink("RCL", colplist, res[1:4,])
    list(datafr = res, oriTable = matFull, IdentResult = IdentResult,
         rowplist = rowplist, colplist = colplist)
  }
ReconsRowLabels =
  function(plist)
  if(is.list(plist)){
    rowvecs = as.list(names(plist))
    for(i in 1:length(rowvecs))
      rowvecs[[i]] = cbind(rowvecs[[i]], ReconsRowLabels(plist[[i]]))
    do.call(rbind, rowvecs)
  } else as.matrix(names(plist))
ReconsColLabels =
  function(plist, matData, rowvecs){
    if(is.list(plist)){
      colvecs = as.list(names(plist))
      for(i in 1:length(colvecs)){
        colvecs[[i]] = cbind(colvecs[[i]],
                 ReconsColLabels(plist[[i]], matData, rowvecs))
        colnames(colvecs[[i]])[1] = "UNKNOWN"
      }
      datfr = do.call(rbind, colvecs)
    }
    else{
      datbit = matData[,plist]
      datlist = NULL
      for(j in 1:ncol(datbit)){
        asnumer = suppressWarnings(as.numeric(datbit[,j]))
        if(all(is.na(datbit[,j])) || !all(is.na(asnumer)))
          datlist[[j]] = asnumer
        else
          datlist[[j]] = datbit[,j]
      }
      datbit = do.call(cbind, datlist)
      ## Specify row.names to avoid annoying warnings
      datfr =
        cbind(as.data.frame(rowvecs, row.names = 1:nrow(rowvecs)), datbit)
      colnames(datfr) =
        c(rep("UNKNOWN", length = ncol(rowvecs)), names(plist))
    }
    datfr
  }
print.plist = function(x, ...){
  plistC = function(plist){
    pLoc = attr(plist, "Loc")
    if(is.list(plist)){
      namevec = names(plist)
      if(!is.null(pLoc))
        namevec = paste0(names(plist),
          " (", pLoc[,"rows"], ", ", pLoc[,"cols"], ")")
      namelist = as.list(namevec)
      for(i in 1:length(namelist))
        namelist[[i]] =
          c(paste("+", namelist[[i]]),
            paste("-", plistC(plist[[i]])))
      do.call(c, namelist)
    } else{
      if(!is.null(names(plist))){
        namevec = names(plist)
        if(!is.null(pLoc))
          namevec = paste0(names(plist),
            " (", plist, ", ", pLoc[,"cols"], ")")
        paste("+", namevec)
      } else paste(plist, collapse = " ")
    }
  }
  
  cat(plistC(x), sep = "\n")
}
attrLoc =
  function(plist, rows = NULL, cols = NULL){
    attr(plist, "Loc") = cbind(rows, cols)
    class(plist) = "plist"
    plist
  }
TCRsink =
  function(ID, ...)
  if(exists("TCRunout", envir = TTLFBaseEnv)){
    varlist = list(...)
    names(varlist) = gsub(" ", "", as.character(match.call()[-(1:2)]))
    with(TTLFBaseEnv, sink(TCRunout))
    for(i in 1:length(varlist)){
      cat("###TCR", ID, names(varlist)[i], "\n")
      print(varlist[[i]])
    }
    sink()
  }
TTLFBaseEnv = new.env()
