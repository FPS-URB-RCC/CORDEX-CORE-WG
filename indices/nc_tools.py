# netcdf tools

errormsg = 'ERROR -- error -- ERROR -- error'
warnmsg = 'WARNING -- warning -- WARNING -- warning'

####### Content
# concatenate_HMTvariable: Function to concatenate a single variable from a series of 
#   similar files located in a folder

def concatenate_HMTvariable(infold, Hn, Mn, Tn, varn, dimn, sliceS):
    """ Function to concatenate a single variable from a series of similar files 
        located in a folder
      infold: folder with the files
      Hn: header section of the names of the files
      Mn: middle section of the names of the files
      Tn: tail section of the names of the files
      varn: name of the variable
      dimn: name of the dimension along which concatenate the variable
      sliceS = [dimn1]:[slice1],...,[dimnN]:[sliceN] slice configuration for the 
          variable in the other dimensions (None, no slice)
        * [integer]: which value of the dimension
        * -1: all along the dimension
        * -9: last value of the dimension
        * [beg]@[end]@[freq] slice from [beg] to [end] every [freq]
    """
    fname = 'concatenate_HMTvariable'
    
    files = gen.files_folder_HMT(infold, Hn, Mn, Tn)
    
    # Getting the length of the final dimension for the concatenation
    Nfiles = len(files)

    if sliceS is not None:
        slicedic = gen.stringS_dictvar(sliceS, ',', ':')
        # Fixing slicedic
        for kn in slicedic.keys():
            dicv = slicedic[kn]
            if dicv.count('@') != 0: slicedic[kn] = gen.str_list_k(dicv, '@', 'I')
            else: slicedic[kn] = int(dicv)
    else:
        slicedic = None
        

    iff = 0
    Lcondim = 0
    fillv = None
    oncs = []
    slicevar = []
    for filen in files:
        onc = NetCDFFile(str(filen), 'r')
        oncs.append(onc)
        
        check_varInfile(fname, filen, onc, varn)
        
        ovar = onc.variables[varn]
        check_dimInvar (fname, filen, onc, varn, dimn, 'concatenating')

        dimvar = list(ovar.dimensions)
        idim = dimvar.index(dimn)
    
        Lcondim = Lcondim + len(onc.dimensions[dimn])
    
        if iff == 0:
            varattrs = ovar.ncattrs()
            varshape = ovar.shape
            vartype = ovar.dtype
            if gen.searchInlist(varattrs, '_FillValue'): 
                fillv = ovar.getncattr('_fillValue')
            if slicedic is not None:
                slicevar, slcd = SliceVarDict_keep(ovar, slicedic)
                varshape = ovar[tuple(slicevar)].shape
            else:
                for dn in dimvar:
                    dsize = len(onc.dimensions[dn])
                    if dn != dimn: 
                        slicevar.append(slice(0,dsize))
                    else:
                        slicevar.append(slice(0,dsize))
            
        iff = iff + 1
    
    # Defining the new variable
    concatshp = list(varshape)
    concatshp[idim] = Lcondim
    
    if fillv is not None:
        concatvar = np.full(tuple(concatshp), fillv)
    else:    
        concatvar = np.full(tuple(concatshp), gen.typemod(0,vartype))
    
    # Filling the variable
    newslicevar = list(slicevar)
    icon = 0
    for iff in range(Nfiles):
        onc = oncs[iff]
        ovar = onc.variables[varn]

        dimcon = ovar.shape[idim]
        slicevar[idim] = slice(0,dimcon)
        newslicevar[idim] = slice(icon,icon+dimcon)

        concatvar[tuple(newslicevar)] = ovar[tuple(slicevar)]

        icon = icon + dimcon

        oncs[iff].close()

    return concatvar

