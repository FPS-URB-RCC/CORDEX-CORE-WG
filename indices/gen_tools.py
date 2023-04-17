# Generic tools

errormsg = 'ERROR -- error -- ERROR -- error'
warnmsg = 'WARNING -- warning -- WARNING -- warning'

####### Content
# Text
# check_arguments: Function to check the number of arguments if they are coincident
# files_folder_HMT: Function to retrieve a list of files from a folder [fold] and 
#   files named [head]*[middle]*[tail]
# Str_Bool: Function to transform from a String value to a boolean one
# str_list: Function to obtain a list from a string giving a split character
# str_list_kinds: Function to obtain a list of types of values from a string giving a 
#   split character and a list of kinds of values
# stringS_dictvar: Function to provide a dictionary from a String which list of DVs 
#   separated [key] [value] couples

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

# Section: Text

def str_list(string, cdiv, empty=False):
    """ Function to obtain a list from a string giving a split character
      string= String from which to obtain a list ('None' for None)
      cdiv= character to use to split the string
      empty= whether empty values should be retained (False, default)
        ....)
    >>> str_list('d01 1995-01-10_00:00:00  alloc_space_field: domain  ' + 
      '          2 ,               4281600  bytes allocated', ' ')
    ['d01', '1995-01-10_00:00:00', 'alloc_space_field:', 'domain', '2', ',', 
     '4281600', 'bytes', 'allocated']
    """
    fname = 'str_list_k'

    if string.find(cdiv) != -1:
        listv = string.split(cdiv)
    else:
        if string == 'None':
            listv = None
        else:
            listv = [string]

    if listv is not None:
        finalist = []
        if not empty:
            for lv in listv:
                if len(lv) > 0: finalist.append(lv)
        else:
            finalist = list(listv)
    else:
        finalist = None

    return finalist

def str_list_kinds(string, cdiv, kinds):
    """ Function to obtain a list of types of values from a string giving a split 
        character and a list of kinds of values
      string= String from which to obtain a list ('None' for None)
      cdiv= character to use to split the string
      kinds= list of kind of desired type of values (as string like: 'np.float', 
        'int', 'np.float64', ....)
    >>> str_list_kinds('1:@#:$:56', ':', ['I', 'S', 'S', 'F'])
    [1, '@#', '$', 56.0]
    >>> str_list_kinds('1:3.4:12.3', ':', ['I', 'np.float64', 'F'])
    [1, 3.4, 12.3]
    """
    fname = 'str_list_kinds'

    if string.count(cdiv) != 0:
        listv = string.split(cdiv)
    else:
        if string == 'None':
            listv = None
        else:
            listv = [string]

    Nv = len(listv)
    newlist = []
    # Checking for 0-length values
    for il in range(Nv):
        if len(listv[il]) > 0: newlist.append(listv[il])
    Nv = len(newlist)
    listv = list(newlist)
   
    Nk = len(kinds)
    if Nv != Nk:
        print (errormsg)
        print ('  ' + fname + ": the amount of values", Nv, "and the amount of " +   \
          "kinds", Nk, 'do not coincide !!')
        print ('  String of values:', string)
        print ('    i value kind _______')
        if Nv > Nk:
            for iv in range(Nk):
                print ('   ', iv, listv[iv], kinds[iv])
            print ('  exceess of values:', listv[Nk:Nv])
        else:
            for iv in range(Nv):
                print ('   ', iv, listv[iv], kinds[iv])
            print ('  exceess of kinds:', kinds[Nv:Nk])
        quit(-1)
            
    if listv is not None:
        finalist = []
        iv = 0
        for lv in listv:
            finalist.append(typemod(lv, kinds[iv]))
            iv = iv + 1
    else:
        finalist = None

    return finalist

def stringS_dictvar(stringv, Dc=',', DVs=':'):
    """ Function to provide a dictionary from a String which list of DVs separated 
          [key] [value] couples
      stringv: String with a dictionary content as [key1][DVs][val1][Dc][key2][DVs][Val2][Dc]...
      Sc: character to separate key values
      DVs: character to separate key-value couples
    >>> stringS_dictvar('i:1,iv:4,vii:7',',',':')
    {'i': '1', 'vii': '7', 'iv': '4'}
    """
    fname = 'stringS_dictvar'

    if stringv.count(Dc) != 0:
        strvv = stringv.split(Dc)
    else:
        strvv = [stringv]

    dictv = {}
    for Svals in strvv:
        keyn = Svals.split(DVs)[0]
        valn = Svals.split(DVs)[1]
        dictv[keyn] = valn

    return dictv

