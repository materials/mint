#!/bin/bash



# Name of directory to copy old files
backupDir="backup"



# Function to exit the script if something went wrong
function exitFunction
{
	echo ""
	echo "====================================="
	echo ""
	echo "Update failed because "$1"."
	echo ""
	echo "Your old files are safe though and located in the \""$backupDir"\" directory that has been created."
	echo ""
	echo "You can download the source code manually in either of the following formats:"
	echo "    Tarball: https://github.com/materials/mint/tarball/master"
	echo "    Archive: https://github.com/materials/mint/archive/master.zip"
	echo ""
	exit
}



# Function to make a backup copy of files
function backupFiles
{
	echo "Creating \""$backupDir"\" directory and moving all files to it."
	origFiles=`ls`
	mkdir $backupDir
	for cur in $origFiles; do
		if [ $cur == "mint" ]; then
			cp $cur $backupDir
			echo "    A copy of the \"mint\" executable was left in the current directory."
		else
			mv $cur $backupDir
		fi
	done
}



# Function to get tarball from github
function getNewCode
{
	echo "Downloading source code and other files from GitHub."
	wget --no-check-certificate -q https://github.com/materials/mint/tarball/master 2> /dev/null || \
		curl -LOks https://github.com/materials/mint/tarball/master > /dev/null || \
		exitFunction "tarball could not be downloaded from GitHub. Make sure that either \"wget\" or \"curl\" is available"
	if [ ! -e master ]; then
		exitFunction "source code was not downloaded correctly"
	fi
}



# Extract code from tarball
function extractNewCode
{
	echo "Extracting new files from the tarball that was downloaded."
	tar -xzf master || exitFunction "tarball extraction was not successful"
	if [ -e master ]; then
		rm -f master
	fi
	newDir="materials-mint"*
	if [ ${#newDir[@]} -eq 0 ]; then
		exitFunction " extracted directory was not found"
	fi
	mv $newDir/* .
	rm -fr $newDir
}



# Copy parameters from old makefile into new one
function updateMakefile
{
	
	# Makefile parameters and defaults
	makeParams=("CC" "g++" "OPT" "\-O3" "COMP" "\-c \-g" "LINK" "\-g" "BLAS" "\-lblas" "LAPACK" "\-llapack" \
			 "DEFINE" "" "MPIRUN" "mpirun \-np")
	
	# Get current values if previous makefile exists
	if [ -e $backupDir/Makefile ]; then
		echo "Saving previous Makefile settings."
		i=0; j=1
		while [ $i -lt ${#makeParams[@]} ]; do
			vars=(`grep ^${makeParams[$i]} $backupDir/Makefile | sed 's/:/ /' | sed 's/=/ /'`)
			if [ ${#vars[@]} -gt 0 ]; then
				makeParams[$j]=`echo ${vars[@]:1} | sed 's/\-/\\\-/g'`
			fi
			i=`expr $i + 2`; j=`expr $j + 2`
		done
	
	# No Makefile existed
	else
		echo "Using default Makefile settings. Check these before compiling!"
	fi
	
	# Set values in Makefile
	i=0; j=1
	while [ $i -lt ${#makeParams[@]} ]; do
		sed 's/#'"${makeDef[$i]}"'#/'"${makeDef[$j]}"'/' < Makefile > temp
		mv temp Makefile
		i=`expr $i + 2`; j=`expr $j + 2`
	done
}



# Run update
backupFiles
getNewCode
extractNewCode
updateMakefile


exit








# Get external version
externVersion=(`wget -qO- palestrina.northwestern.edu/mint/version`)
if [[ ${#externVersion[@]} -lt 4 ]]; then
	echo "Could not connect to server"
	exit
fi

# Check if version file exists
update=0
if [[ ! -e version ]]; then
	update=1
	localVersion=(0 0 0 0)

# Only compare version if file exists
else

	# Get local version
	localVersion=(`cat version`)

	# Check if versions are different
	if [[ ${localVersion[0]} -ne ${externVersion[0]} ]]; then
		update=1
	elif [[ 10#${localVersion[1]} -ne 10#${externVersion[1]} ]]; then
		update=1
	elif [[ 10#${localVersion[2]} -ne 10#${externVersion[2]} ]]; then
		update=1
	elif [[ ${localVersion[3]} -ne ${externVersion[3]} ]]; then
		update=1
	fi
fi

# Code is up to date
if [[ $update -eq 0 ]]; then
	echo "Code is up to date so no changes will be made"

# Need to update code
else

	# Set versions
	oldVersion=${localVersion[0]}"."${localVersion[1]}"."${localVersion[2]}"."${localVersion[3]}
	newVersion=${externVersion[0]}"."${externVersion[1]}"."${externVersion[2]}"."${externVersion[3]}

	# Backup files
	echo "Moving files in current directory to "$oldVersion
	mkdir $oldVersion
	for cur in `ls`; do
		if [[ $cur == "mint" ]]; then
			cp $cur $oldVersion
		elif [[ $cur != $oldVersion && $cur != ".." && $cur != "." ]]; then
			mv $cur $oldVersion
		fi
	done

	# Output
	echo "Updating code to version "$newVersion

	# Copy files
	wget -qr --no-parent --reject "index.html*" palestrina.northwestern.edu/mint/
	mv palestrina.northwestern.edu/mint/* .
        chmod u+x update.sh
	rm -fr palestrina.northwestern.edu    
	
	# Makefile parameters and defaults
	makeDef=("CC" "g++" "OPT" "\-O3" "CFLAGS" "\-c \-g" "LFLAGS" "\-g" "BLAS" "\-lblas" "LAPACK" "\-llapack" \
			 "DEFINE" "" "MPIRUN" "mpirun \-np")
	
	# Loop over parameters and set their values in the makefile
	i=0
	j=1
	while [[ $i -lt ${#makeDef[@]} ]]; do
	
		# Makefile does not exist so set value
		if [[ ! -e $oldVersion"/Makefile" ]]; then
			sed 's/#'"${makeDef[$i]}"'#/'"${makeDef[$j]}"'/' < Makefile > temp
		
		# Makefile exists so get current value
		else
			
			# Get current value
			vals=()
			while read line; do
				vars=(`echo $line | sed 's/:/ /' | sed 's/=/ /'`)
				if [[ ${#vars[@]} -gt 0 ]]; then
					if [[ ${vars[0]} == ${makeDef[$i]} ]]; then
						vals=(${vars[@]})
						break
					fi
				fi
			done < $oldVersion"/Makefile"
			
			# Value was not found
			if [[ ${#vals[@]} -eq 0 ]]; then
				sed 's/#'"${makeDef[$i]}"'#/'"${makeDef[$j]}"'/' < Makefile > temp
			
			# Value was found
			else
				line=`echo ${vals[@]:1} | sed 's/\-/\\\-/g'`
				sed 's/#'"${makeDef[$i]}"'#/'"$line"'/' < Makefile > temp
			fi
		fi
		
		# Replace makefile
		if [[ -e temp ]]; then
			mv temp Makefile
		fi
		
		# Go to next value
		i=`expr $i + 2`
		j=`expr $j + 2`
	done
fi
