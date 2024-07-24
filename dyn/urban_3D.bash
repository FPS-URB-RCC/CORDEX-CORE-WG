#!/bin/bash

pyBIN=/usr/bin/python3
export HDF5_USE_FILE_LOCKING=FALSE

fscratch=false
gscratch=false

# List of cities to process
#   [dom]|[CityName]|[lon],[lat]
#cities='EUR-11|Paris|2.35,48.85:EUR-11|Berlin|13.4039,52.4683:'
#cities=${cities}'EUR-11|London|-0.13,51.50:EUR-11|Madrid|-3.70,40.42:'
#cities=${cities}'SAM-22|Buenos!Aires|-58.416,-34.559:'
#cities=${cities}'EAS-22|Beijing|116.41,39.90'
cityfilen='CORDEX-CORE-WG/city_info.csv'
#cityfilen='city_info_mod.csv'
cities=`cat ${cityfilen} | tr ' ' '!' | tr ',' ' ' | awk '{print $10"|"$1"|"$7","$8}' | tr '\n' ':'`

## Single city
#cityval='Mumbai,India,1074.0,21755882.364945933,20256.873710377964,1.79,73.03655296648097,19.13765841208688,Am,WAS'
#cities=`echo ${cityval} | tr ' ' '!' | tr ',' ' ' | awk '{print $10"|"$1"|"$7","$8}' | tr '\n' ':'`

infold=/datos/MOD/CORDEX/CMIP5

# Variable to process
varn=ta

# Diagnostics for a given variable
#   [varn]|[diagn]
diagn='ta|tha:ua|ua:va|va'

# +/- points from city center 
Npoints=10

#######    ######
## MAIN
    #######
errmsg="ERROR -- error -- ERROR -- error"
cits=`echo ${cities}  | tr ':' ' '`
diags=`echo ${diagn}  | tr ':' ' '`
Lvarn=`expr length ${varn}`
for city in ${cits}; do
  if test ${city} = 'domain|city|lon,lat'; then continue; fi
  dom=`echo ${city} | tr '|' ' ' | awk '{print $1}'`
  if test ${dom} = 'EUR'; then
    dom=${dom}-11
  else
    dom=${dom}-22
  fi
  citn=`echo ${city} | tr '|' ' ' | awk '{print $2}'`
  if test ${citn} = 'UNNAMED'; then continue; fi
  lonlat=`echo ${city} | tr '|' ' ' | awk '{print $3}'`

  echo "city: "${domn}" "${citn}" "${lonlat}
  cityS=`echo ${citn} | tr '!' '_'`

  idir=${infold}/${dom}
  gcms=`ls -1 ${idir}`

  for gcm in ${gcms}; do
    idirg=${infold}/${dom}/${gcm}/historical/r1i1p1/*/*/mon
    echo "idirg: "${idirg}
    vars=`ls -1 ${idirg}`

    for vn in ${vars}; do
      #echo "     vn: "${vn}
      if test ${vn:0:${Lvarn}} = ${varn}; then
        idirv=${idirg}/${vn}/*
        files=`ls -1 ${idirv}/${vn}_*nc`

        Nfiles=`echo ${files} | wc -w | awk '{print $1}'`
        if test ${Nfiles} -ne 0; then
          echo " ${gcm} ${vn} files found !!"

          ifile=`echo ${files} | awk '{print $1}'`

	  vals='lon:lat:time|-9'
          ijcity=`$pyBIN $pyHOME3/nc_var.py -o get_point -S ${vals} -f ${ifile} -v ${lonlat}`
	  if test $? -ne 0; then
            echo ${errmsg}
	    echo "  python failed !!"
	    echo nc_var.py -o get_point -S ${vals} -f ${ifile} -v ${lonlat}
	    exit -1
	  fi
	  word1=`echo ${ijcity} | awk '{print $1}'`
	  if test $word1 = 'INFORMATION'; then
	    ijcity=`echo ${ijcity} | awk '{print $33}'`
	  fi
	  ipt=`echo ${ijcity} | tr ',' ' ' | awk '{print $1}'` 
	  jpt=`echo ${ijcity} | tr ',' ' ' | awk '{print $2}'` 

	  ipti=`expr ${ipt} - ${Npoints}`
	  ipte=`expr ${ipt} + ${Npoints}`
	  jpti=`expr ${jpt} - ${Npoints}`
	  jpte=`expr ${jpt} + ${Npoints}`

	  odir=${dom}/${cityS}
	  mkdir -p ${odir}

	  ifileS=`basename ${ifile}`
	  ifiledir=`dirname ${ifile}`

	  rcm=`echo ${ifileS} | tr '_' ' ' | awk '{print $6}'`
	  echo "  gcm: "${gcm}" rcm: "${rcm}

	  # x,y dims are center dependent! Not fully CF compilant !!
          # FROM: https://stackoverflow.com/questions/21688553/bash-expr-index-command
	  Lrcm=`expr length ${rcm}`
	  prefix=${rcm%%ICTP*}
	  iposICTP=${#prefix}
	  prefix=${rcm%%GERICS*}
	  iposGERICS=${#prefix}
	  echo "ICTP: "${iposICTP}" GERICS: "${iposGERICS}
	  if test $iposICTP != ${Lrcm}; then
            xdimn='x'
	    ydimn='y'
	  elif test $iposGERICS != ${Lrcm}; then
            xdimn='rlon'
            ydimn='rlat'
          else
            echo "avoiding non CORDEX-CORE ..."
	    continue
          fi
	  # Domain extremes
	  values='dim:'${xdimn}':None:None'
          xdS=`$pyBIN $pyHOME3/nc_var.py -o get_file_inf -S ${values} -f ${ifile}`
	  if test $? -ne 0; then
            echo ${errmsg}
	    echo "  python failed !!"
	    echo nc_var.py -o get_file_inf -S ${values} -f ${ifile}
	    exit -1
	  fi
	  values='dim:'${ydimn}':None:None'
          ydS=`$pyBIN $pyHOME3/nc_var.py -o get_file_inf -S ${values} -f ${ifile}`
	  if test $? -ne 0; then
            echo ${errmsg}
	    echo "  python failed !!"
	    echo nc_var.py -o get_file_inf -S ${values} -f ${ifile}
	    exit -1
	  fi
          domdimx=`echo ${xdS} | tr '|' ' ' | awk '{print $2}'`
          domdimy=`echo ${ydS} | tr '|' ' ' | awk '{print $2}'`

	  echo "    ij SW: "${ipti}", "${jpti}" NE: "${ipte}", "${jpte} 
	  echo "    "${dom}" domain size: "${domdimx}", "${domdimy}
	  if test ${ipti} -lt 1 || test ${jpti} -lt 1 || test ${ipte} -gt ${domdimx} \
	    || test ${jpte} -gt ${domdimy}; then
            echo ${errmsg}
	    echo "  city '"${cityn}"' too close to borders of CORDEX domain '"${dom}"' !!"
	    echo "    select another domain for the city"
	    echo "    ij SW: "${ipti}", "${jpti}" NE: "${ipte}", "${jpte} 
	    echo "    "${dom}" domain size: "${domdimx}", "${domdimy}
	    exit -1
	  fi

	  ofile=${odir}/${vn}_${gcm}_${rcm}.nc
	  ofile2=${odir}/${vn}_${gcm}_${rcm}_ymonmean.nc
	  vnS=${vn//[[:digit:]]}
	  if ${fscratch}; then rm ${ofile}; fi
          if test ! -f ${ofile}; then
	    values=${ydimn}'|'${jpti}'@'${jpte}'@1:'${xdimn}'|'${ipti}'@'${ipte}'@1,'
	    values=${values}${ifiledir}',time,time'
	    ${pyBIN} ${pyHOME3}/nc_var.py -o netcdf_fold_slice_concatenation_HMT -S ${values}\
              -v all -f ${vn},${rcm},nc
  	    if test $? -ne 0; then
              echo ${errmsg}
	      echo "  python failed !!"
	      echo nc_var.py -o netcdf_fold_slice_concatenation_HMT -S ${values} -v all \
		-f ${vn},${rcm},nc
	      exit -1
	    fi
	    presv=`echo ${vn} | tr -d [a-z]`
            #${pyBIN} ${pyHOME3}/nc_var.py -o addDim -S 1 -f netcdf_fold_slice_concatenated_HMT.nc -v plev
  	    #if test $? -ne 0; then
            #  echo ${errmsg}
	    #  echo "  python failed !!"
	    #  echo nc_var.py -o addDim -S 1 -f netcdf_fold_slice_concatenated_HMT.nc -v plev
	    #  exit -1
            #fi
	    #vals='plev|pressure@Pressure@Pa|f'
            #${pyBIN} ${pyHOME3}/nc_var.py -o addVar -S ${vals}                       \
	    #  -f netcdf_fold_slice_concatenated_HMT.nc -v pressure
  	    #if test $? -ne 0; then
            #  echo ${errmsg}
	    #  echo "  python failed !!"
	    #  echo nc_var.py -o addVar -S ${vals} -f netcdf_fold_slice_concatenated_HMT.nc  \
#		-v pressure
	    #  exit -1
            #fi
	    #${pyBIN} ${pyHOME3}/nc_var.py -o fill_varNCfile                          \
	    #  -f netcdf_fold_slice_concatenated_HMT.nc -S constant,${presv}00 -f pressure
  	    #if test $? -ne 0; then
            #  echo ${errmsg}
	    #  echo "  python failed !!"
	    #  echo nc_var.py -o fill_varNCfile -f netcdf_fold_slice_concatenated_HMT.nc  \
#		-S constant,${presv}00 -f pressure
#	      exit -1
#	    fi
            mv netcdf_fold_slice_concatenated_HMT.nc ${ofile}

	    # Statistics
	    values=${xdimn},${ydimn}:cycle,time,0,12,1,12:yes:mean:yes
	    ${pyBIN} ${pyHOME3}/nc_var.py -o stats_varfile -S ${values} -f ${ofile} -v ${vn}
	    echo nc_var.py -o stats_varfile -S ${values} -f ${ofile} -v ${vn}
  	    if test $? -ne 0; then
              echo ${errmsg}
	      echo "  python failed !!"
	      echo nc_var.py -o stats_varfile -S ${values} -f ${ofile} -v ${vn}
	      exit -1
	    fi

	    ${pyBIN} ${pyHOME3}/nc_var.py -o fvaradd -S ${ofile} -f stats_varfile.nc -v lon
	    ${pyBIN} ${pyHOME3}/nc_var.py -o fvaradd -S ${ofile} -f stats_varfile.nc -v lat
	    mv stats_varfile.nc ${ofile2}

	  fi

        fi
      fi

    # end of variables
    done

    # Plotting
    Lofile2=`expr length ${ofile2}0`
    if test ${Lofile2} -lt 2; then continue; fi
    echo "ofile2: "${ofile2}
    if test -f ${ofile2}; then
      din=`echo ${diagn} | tr ':' '\n' | grep ${vnS} | tr '|' ' ' | awk '{print $2}'`
      echo "    "${vnS}" - "${din}

      ofig=${odir}/${din}ymon-anom_${gcm}_${rcm}_SN.png
      if ${gscratch}; then rm ${ofig}; fi
      if test ! -f ${ofig}; then
	echo "  plotting '"${ofig}"' ..."
        python3 tamon_pressure_ver.py -d ${dom} -c ${cityS} -g ${gcm} -r ${rcm} -v no
        if test $? -ne 0; then
          echo ${errmsg}
          echo "  python failed !!"
          echo tamon_pressure_ver.py -d ${dom} -c ${cityS} -g ${gcm} -r ${rcm} -v no
        fi
	#exit
      fi
    fi

    #exit
  # end of gcms
  done

# end of cities
done
