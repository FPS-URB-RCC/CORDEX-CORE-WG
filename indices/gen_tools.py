# Generic tools

errormsg = 'ERROR -- error -- ERROR -- error'
warnmsg = 'WARNING -- warning -- WARNING -- warning'

####### Content
# check_arguments: Function to check the number of arguments if they are coincident
# files_folder_HMT: Function to retrieve a list of files from a folder [fold] and 
#   files named [head]*[middle]*[tail]
# Str_Bool: Function to transform from a String value to a boolean one

def files_folder_HMT(folder='.',head='',middle='',tail='', rmfolder=True):
    """ Function to retrieve a list of files from a folder [fold] and files named 
         [head]*[middle]*[tail]
      [folder]: fname of the folder with files ('.', default)
      [head]: header of the name of the files ('', default)
      [middle]: middle part of the name of the files ('', default)
      [tail]: ending part of the name of the files ('', default)
      [rmfolder]: whether folder section should be removed from the found files 
        (True, default)
    >>> files_folder_HMT('.','t','wrf','mean.nc')
    ['clt_wrfout_tmean.nc', 'tas_wrfout_tmean.nc', 'ta_wrfout_xmean.nc']
    >>> files_folder_HMT(folder='/WRF/current', head='hfls', tail='.nc')
    ['hfls_wrfout.nc', 'hfls_wrfout_tmean.nc']
    """
    fname = 'files_folder_HMT'

    if folder == 'h':
        print ('  ' + fname + '_____________________________________________________')
        print (files_folder_HMT.__doc__)
        quit()

    # FROM: http://stackoverflow.com/questions/9997048/python-subprocess-wildcard-usage
    ins = folder + "/" + head + "*" + middle + '*' + tail
    files = sub.Popen("/bin/ls -1R " + ins, shell=True, stdout=sub.PIPE)
    fileslist = files.communicate()
    listfiles0 = str(fileslist).replace("'",'').replace('(','').replace(')','').split('\\n')


    # Filtering output
    listfiles = []
    if rmfolder:
        iif = 0
        for ifile in listfiles0:
            # Removing folder section
            filen = ifile[len(folder)+1:len(ifile)+1]
            if ifile.find(head) != -1 and  ifile.find(middle) != -1 and              \
              ifile.find(tail) != -1:
                indhead = filen.index(head)
                indmiddle = filen.index(middle)
                if indmiddle == 0: indmiddle = indhead
                indtail = filen.index(tail)
                if indtail == 0: indtail = indhead + indmiddle
                if indhead <= indmiddle <= indtail: 
                    if iif == 0:
                        # Removing 'b' at the beginning
                        Lfile = len(ifile)
                        listfiles.append(ifile[1:Lfile])
                    else:
                        listfiles.append(ifile)
            iif = iif + 1
    else:
        listfiles = list(listfiles0)

    return listfiles

def Str_Bool(val):
    """ Function to transform from a String value to a boolean one
    >>> Str_Bool('True')
    True
    >>> Str_Bool('0')
    False
    >>> Str_Bool('no')
    False
    """

    fname = 'Str_Bool'

    if val == 'True' or val == 'true' or val == '1' or val == 'yes': 
        boolv = True
    elif val == 'False' or val == 'false' or val == '0' or val== 'no':
        boolv = False
    else:
        print (errormsg)
        print ('  ' + fname + ": value '" + val + "' not ready!!")
        quit(-1)

    return boolv

def check_arguments(funcname, args, expectargs, char):
    """ Function to check the number of arguments if they are coincident
    check_arguments(funcname,Nargs,args,char)
      funcname= name of the function/program to check
      args= passed arguments
      expectargs= expected arguments
      char= character used to split the arguments
    """
    fname = 'check_arguments'

    Nvals = len(args.split(char))
    Nargs = len(expectargs.split(char))

    if Nvals != Nargs:
        print (errormsg)
        print ('  '+fname + ': wrong number of arguments:', Nvals, " passed to  '",  \
          funcname, "' which requires:", Nargs, '!!')
        print ("    passed string: '" + args + "'")
        print ("    expected string: '" + expectargs + "'")
        print ('     # given expected _______')
        Nmax = np.max([Nvals, Nargs])
        for ia in range(Nmax):
            if ia > Nvals-1:
                aval = ' '
            else:
                aval = args.split(char)[ia]
            if ia > Nargs-1:
                aexp = ' '
            else:
                aexp = expectargs.split(char)[ia]

            print ('    ', ia, aval, aexp)
        quit(-1)

    return


