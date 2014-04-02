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
	echo "Your old files are safe though and located in the directory \""$backupDir"\"."
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



# Get dlib package
function getDLIB
{
	if [ ! -d "dlib" ]; then
		
		# Download the package
		echo "Downloading dlib package."
		wget -O dlib.tar.bz2 --no-check-certificate -q http://sourceforge.net/projects/dclib/files/latest/download 2> /dev/null || \
			curl -Lks http://sourceforge.net/projects/dclib/files/latest/download 1> dlib.tar.bz2 2> /dev/null || \
			exitFunction "dlib package could not be downloaded"
		
		# Unzip the directory and delete extra files
		mkdir dlib dlibTemp
		tar -xjf dlib.tar.bz2 -C dlibTemp
		mv dlibTemp/*/dlib/* dlib
		rm -fr dlibTemp dlib.tar.bz2
	fi
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
				if [ ${vars[1]:0:1} != '#' ]; then
					makeParams[$j]=`echo ${vars[@]:1} | sed 's/\-/\\\-/g'`
				fi
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
		sed 's/#'"${makeParams[$i]}"'#/'"${makeParams[$j]}"'/' < Makefile > temp
		mv temp Makefile
		i=`expr $i + 2`; j=`expr $j + 2`
	done
}



# Run update
backupFiles
getNewCode
extractNewCode
getDLIB
updateMakefile
echo "Update finished successfully. It should be safe to delete the directory \""$backupDir"\"."
echo "Check \"Makefile\" and run \"make\" once you're ready to compile."
