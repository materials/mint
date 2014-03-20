#====================================================================#
#                                                                    #
#   This is the Makefile for the Mint program. It should be edited   #
#   to match your system.                                            #
#                                                                    #
#   Once edited, simply type "make" from the command line to build   #
#   the mint executable.                                             #
#                                                                    #
#====================================================================#
#                                                                    #
#   Contact: Kyle Michel (kylemichel@gmail.com)                      #
#            Logan Ward (LoganWard2012@u.northwestern.edu)           #
#                                                                    #
#====================================================================#
#                                                                    #
#   Copyright 2011-2014 Kyle Michel, Logan Ward                      #
#                                                                    #
#                                                                    #
#   This file is part of Mint.                                       #
#                                                                    #
#   Mint is free software: you can redistribute it and/or modify     #
#   it under the terms of the GNU Lesser General Public License as   #
#   published by the Free Software Foundation, either version 3 of   #
#   the License, or (at your option) any later version.              #
#                                                                    #
#   Mint is distributed in the hope that it will be useful, but      #
#   WITHOUT ANY WARRANTY; without even the implied warranty of       #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    #
#   GNU Lesser General Public License for more details.              #
#                                                                    #
#   You should have received a copy of the GNU Lesser General        #
#   Public License along with Mint.  If not, see                     #
#   <http://www.gnu.org/licenses/>.                                  #
#                                                                    #
#====================================================================#



# C++ compiler and optimization level
# If you are compiling with MPI support then the compiler must be MPI compatible!
CC  := #CC#
OPT := #OPT#



# Any other compiler and linker flags
# Mint requires linking to blas/lapack libraries
COMP   := #COMP#
LINK   := #LINK#
BLAS   := #BLAS#
LAPACK := #LAPACK#



# Definitions (each separated by a space)
# Possible values:
#    MPI - enable mpi support (compiler must be MPI compatible!)
#    MKL - using MKL libraries (do not add this definition if MKL is not used)
DEFINE := #DEFINE#



# String used to launch an MPI program for current system (mpirun, aprun, etc)
# If number of processors is not determined automatically, then use e.g. "mpirun -np"
# This only matters if you will be launching external mpi programs from within program
MPIRUN := #MPIRUN#



#====================================================================#
#                                                                    #
#   You should not have to edit anything past this point.            #
#                                                                    #
#====================================================================#

# Version
VERSION := "2014.03.19.1"

# Directories
SRCD := src
OBJD := obj

# Setup
EXE   := mint
FDEF  := $(addprefix -DMINT_,$(DEFINE))
FMPI  := -DMPIRUN=\"$(MPIRUN)\"
FVERS := -DVERSION=\"$(VERSION)\"
FTIME := -DCOMP=\""`date +'%B %d, %Y at %r %Z'`"\"
FALL  := $(COMP) $(OPT) $(FDEF)

# Build everything
all : touchAbout $(EXE)
touchAbout :
	@touch $(SRCD)/about.cpp

# Builds for object files
$(OBJD)/about.o : $(SRCD)/about.cpp $(SRCD)/about.h $(SRCD)/output.h $(SRCD)/num.h $(SRCD)/list.h $(SRCD)/text.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(FVERS) $(FTIME) $(SRCD)/about.cpp -o $@
$(OBJD)/bonds.o : $(SRCD)/bonds.cpp $(SRCD)/bonds.h $(SRCD)/iso.h $(SRCD)/elements.h $(SRCD)/num.h $(SRCD)/text.h $(SRCD)/list.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/bonds.cpp -o $@
$(OBJD)/cif.o : $(SRCD)/cif.cpp $(SRCD)/cif.h $(SRCD)/elements.h $(SRCD)/symmetry.h $(SRCD)/spaceGroup.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/num.h $(SRCD)/iso.h $(SRCD)/text.h $(SRCD)/list.h $(SRCD)/fileSystem.h $(SRCD)/pointGroup.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/cif.cpp -o $@
$(OBJD)/constants.o : $(SRCD)/constants.cpp $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/constants.cpp -o $@
$(OBJD)/crystalMaker.o : $(SRCD)/crystalMaker.cpp $(SRCD)/crystalMaker.h $(SRCD)/bonds.h $(SRCD)/language.h $(SRCD)/num.h $(SRCD)/output.h $(SRCD)/iso.h $(SRCD)/text.h $(SRCD)/elements.h $(SRCD)/list.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/crystalMaker.cpp -o $@
$(OBJD)/diffraction.o : $(SRCD)/diffraction.cpp $(SRCD)/multi.h $(SRCD)/diffraction.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/iso.h $(SRCD)/elements.h $(SRCD)/symmetry.h $(SRCD)/fileSystem.h $(SRCD)/list.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/diffraction.cpp -o $@
$(OBJD)/elements.o : $(SRCD)/elements.cpp $(SRCD)/elements.h $(SRCD)/output.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/list.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/elements.cpp -o $@
$(OBJD)/espresso.o : $(SRCD)/espresso.cpp $(SRCD)/multi.h $(SRCD)/espresso.h $(SRCD)/structureIO.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/iso.h $(SRCD)/elements.h $(SRCD)/kpoints.h $(SRCD)/list.h $(SRCD)/fileSystem.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/espresso.cpp -o $@
$(OBJD)/ewald.o : $(SRCD)/ewald.cpp $(SRCD)/multi.h $(SRCD)/num.h $(SRCD)/ewald.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/text.h $(SRCD)/list.h $(SRCD)/constants.h $(SRCD)/locPotential.h $(SRCD)/iso.h $(SRCD)/elements.h $(SRCD)/symmetry.h $(SRCD)/potential.h $(SRCD)/fileSystem.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/ewald.cpp -o $@
$(OBJD)/extPotential.o : $(SRCD)/extPotential.cpp $(SRCD)/extPotential.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/vasp.h $(SRCD)/espresso.h $(SRCD)/potential.h $(SRCD)/iso.h $(SRCD)/elements.h $(SRCD)/symmetry.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/list.h $(SRCD)/fileSystem.h $(SRCD)/kpoints.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/extPotential.cpp -o $@
$(OBJD)/fileSystem.o : $(SRCD)/fileSystem.cpp $(SRCD)/multi.h $(SRCD)/fileSystem.h $(SRCD)/output.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/list.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/fileSystem.cpp -o $@
$(OBJD)/findsym.o : $(SRCD)/findsym.cpp $(SRCD)/findsym.h $(SRCD)/num.h $(SRCD)/output.h $(SRCD)/iso.h $(SRCD)/text.h $(SRCD)/list.h $(SRCD)/constants.h $(SRCD)/elements.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/findsym.cpp -o $@
$(OBJD)/gaPredict.o : $(SRCD)/gaPredict.cpp $(SRCD)/gaPredict.h $(SRCD)/structureIO.h $(SRCD)/randomStructure.h $(SRCD)/fileSystem.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/num.h $(SRCD)/ga.h $(SRCD)/iso.h $(SRCD)/symmetry.h $(SRCD)/potential.h $(SRCD)/diffraction.h $(SRCD)/random.h $(SRCD)/text.h $(SRCD)/list.h $(SRCD)/spaceGroup.h $(SRCD)/constants.h $(SRCD)/mtwist.h $(SRCD)/randistrs.h $(SRCD)/elements.h $(SRCD)/pointGroup.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/gaPredict.cpp -o $@
$(OBJD)/help.o : $(SRCD)/help.cpp $(SRCD)/help.h $(SRCD)/output.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/list.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/help.cpp -o $@
$(OBJD)/interstitial.o : $(SRCD)/interstitial.cpp $(SRCD)/multi.h $(SRCD)/interstitial.h $(SRCD)/constants.h $(SRCD)/output.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/iso.h $(SRCD)/symmetry.h $(SRCD)/list.h $(SRCD)/elements.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/interstitial.cpp -o $@
$(OBJD)/iso.o : $(SRCD)/iso.cpp $(SRCD)/multi.h $(SRCD)/iso.h $(SRCD)/output.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/elements.h $(SRCD)/list.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/iso.cpp -o $@
$(OBJD)/json.o : $(SRCD)/json.cpp $(SRCD)/json.h $(SRCD)/bonds.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/list.h $(SRCD)/iso.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/elements.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/json.cpp -o $@
$(OBJD)/kmc.o : $(SRCD)/kmc.cpp $(SRCD)/kmc.h $(SRCD)/symmetry.h $(SRCD)/unique.h $(SRCD)/language.h $(SRCD)/num.h $(SRCD)/iso.h $(SRCD)/elements.h $(SRCD)/structureIO.h $(SRCD)/text.h $(SRCD)/random.h $(SRCD)/fileSystem.h $(SRCD)/constants.h $(SRCD)/output.h $(SRCD)/list.h $(SRCD)/mtwist.h $(SRCD)/randistrs.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/kmc.cpp -o $@
$(OBJD)/kpoints.o : $(SRCD)/kpoints.cpp $(SRCD)/kpoints.h $(SRCD)/output.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/list.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/kpoints.cpp -o $@
$(OBJD)/language.o : $(SRCD)/language.cpp $(SRCD)/language.h $(SRCD)/text.h $(SRCD)/list.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/language.cpp -o $@
$(OBJD)/launcher.o : $(SRCD)/launcher.cpp $(SRCD)/multi.h $(SRCD)/num.h $(SRCD)/launcher.h $(SRCD)/settings.h $(SRCD)/about.h $(SRCD)/help.h $(SRCD)/randomStructure.h $(SRCD)/unique.h $(SRCD)/pointGroup.h $(SRCD)/spaceGroup.h $(SRCD)/interstitial.h $(SRCD)/gaPredict.h $(SRCD)/pdf.h $(SRCD)/fileSystem.h $(SRCD)/language.h $(SRCD)/timer.h $(SRCD)/output.h $(SRCD)/text.h $(SRCD)/list.h $(SRCD)/constants.h $(SRCD)/iso.h $(SRCD)/structureIO.h $(SRCD)/symmetry.h $(SRCD)/potential.h $(SRCD)/phonons.h $(SRCD)/kmc.h $(SRCD)/diffraction.h $(SRCD)/random.h $(SRCD)/ga.h $(SRCD)/elements.h $(SRCD)/mtwist.h $(SRCD)/randistrs.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/launcher.cpp -o $@
$(OBJD)/locPotential.o : $(SRCD)/locPotential.cpp $(SRCD)/locPotential.h $(SRCD)/pairPotential.h $(SRCD)/ewald.h $(SRCD)/relax.h $(SRCD)/output.h $(SRCD)/potential.h $(SRCD)/iso.h $(SRCD)/symmetry.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/list.h $(SRCD)/elements.h $(SRCD)/constants.h $(SRCD)/fileSystem.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/locPotential.cpp -o $@
$(OBJD)/mint.o : $(SRCD)/mint.cpp $(SRCD)/multi.h $(SRCD)/output.h $(SRCD)/launcher.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/list.h $(SRCD)/iso.h $(SRCD)/structureIO.h $(SRCD)/symmetry.h $(SRCD)/potential.h $(SRCD)/phonons.h $(SRCD)/kmc.h $(SRCD)/diffraction.h $(SRCD)/random.h $(SRCD)/elements.h $(SRCD)/constants.h $(SRCD)/fileSystem.h $(SRCD)/mtwist.h $(SRCD)/randistrs.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/mint.cpp -o $@
$(OBJD)/mintStructure.o : $(SRCD)/mintStructure.cpp $(SRCD)/num.h $(SRCD)/mintStructure.h $(SRCD)/elements.h $(SRCD)/symmetry.h $(SRCD)/spaceGroup.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/list.h $(SRCD)/constants.h $(SRCD)/iso.h $(SRCD)/text.h $(SRCD)/fileSystem.h $(SRCD)/pointGroup.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/mintStructure.cpp -o $@
$(OBJD)/mtwist.o : $(SRCD)/mtwist.c $(SRCD)/mtwist.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/mtwist.c -o $@
$(OBJD)/multi.o : $(SRCD)/multi.cpp $(SRCD)/multi.h $(SRCD)/language.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/list.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(FMPI) $(SRCD)/multi.cpp -o $@
$(OBJD)/output.o : $(SRCD)/output.cpp $(SRCD)/multi.h $(SRCD)/output.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/list.h $(SRCD)/constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/output.cpp -o $@
$(OBJD)/pairPotential.o : $(SRCD)/pairPotential.cpp $(SRCD)/multi.h $(SRCD)/pairPotential.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/list.h $(SRCD)/constants.h $(SRCD)/locPotential.h $(SRCD)/iso.h $(SRCD)/elements.h $(SRCD)/symmetry.h $(SRCD)/potential.h $(SRCD)/fileSystem.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/pairPotential.cpp -o $@
$(OBJD)/pdf.o : $(SRCD)/pdf.cpp $(SRCD)/pdf.h $(SRCD)/elements.h $(SRCD)/num.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/iso.h $(SRCD)/diffraction.h $(SRCD)/text.h $(SRCD)/fileSystem.h $(SRCD)/list.h $(SRCD)/constants.h $(SRCD)/symmetry.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/pdf.cpp -o $@
$(OBJD)/phonons.o : $(SRCD)/phonons.cpp $(SRCD)/phonons.h $(SRCD)/language.h $(SRCD)/constants.h $(SRCD)/output.h $(SRCD)/iso.h $(SRCD)/symmetry.h $(SRCD)/potential.h $(SRCD)/text.h $(SRCD)/fileSystem.h $(SRCD)/num.h $(SRCD)/list.h $(SRCD)/elements.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/phonons.cpp -o $@
$(OBJD)/pointGroup.o : $(SRCD)/pointGroup.cpp $(SRCD)/pointGroup.h $(SRCD)/output.h $(SRCD)/num.h $(SRCD)/iso.h $(SRCD)/symmetry.h $(SRCD)/text.h $(SRCD)/list.h $(SRCD)/constants.h $(SRCD)/elements.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/pointGroup.cpp -o $@
$(OBJD)/potential.o : $(SRCD)/potential.cpp $(SRCD)/potential.h $(SRCD)/locPotential.h $(SRCD)/extPotential.h $(SRCD)/language.h $(SRCD)/num.h $(SRCD)/iso.h $(SRCD)/elements.h $(SRCD)/symmetry.h $(SRCD)/fileSystem.h $(SRCD)/list.h $(SRCD)/text.h $(SRCD)/output.h $(SRCD)/constants.h $(SRCD)/vasp.h $(SRCD)/espresso.h $(SRCD)/kpoints.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/potential.cpp -o $@
$(OBJD)/randistrs.o : $(SRCD)/randistrs.c $(SRCD)/mtwist.h $(SRCD)/randistrs.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/randistrs.c -o $@
$(OBJD)/random.o : $(SRCD)/random.cpp $(SRCD)/num.h $(SRCD)/multi.h $(SRCD)/random.h $(SRCD)/output.h $(SRCD)/constants.h $(SRCD)/list.h $(SRCD)/text.h $(SRCD)/mtwist.h $(SRCD)/randistrs.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/random.cpp -o $@
$(OBJD)/randomStructure.o : $(SRCD)/randomStructure.cpp $(SRCD)/randomStructure.h $(SRCD)/output.h $(SRCD)/constants.h $(SRCD)/num.h $(SRCD)/iso.h $(SRCD)/spaceGroup.h $(SRCD)/symmetry.h $(SRCD)/random.h $(SRCD)/list.h $(SRCD)/text.h $(SRCD)/elements.h $(SRCD)/pointGroup.h $(SRCD)/mtwist.h $(SRCD)/randistrs.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/randomStructure.cpp -o $@
$(OBJD)/relax.o : $(SRCD)/relax.cpp $(SRCD)/relax.h $(SRCD)/output.h $(SRCD)/iso.h $(SRCD)/symmetry.h $(SRCD)/locPotential.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/list.h $(SRCD)/constants.h $(SRCD)/elements.h $(SRCD)/potential.h $(SRCD)/fileSystem.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/relax.cpp -o $@
$(OBJD)/settings.o : $(SRCD)/settings.cpp $(SRCD)/settings.h $(SRCD)/language.h $(SRCD)/fileSystem.h $(SRCD)/output.h $(SRCD)/iso.h $(SRCD)/structureIO.h $(SRCD)/ga.h $(SRCD)/gaPredict.h $(SRCD)/text.h $(SRCD)/list.h $(SRCD)/num.h $(SRCD)/constants.h $(SRCD)/elements.h $(SRCD)/random.h $(SRCD)/mtwist.h $(SRCD)/randistrs.h $(SRCD)/symmetry.h $(SRCD)/potential.h $(SRCD)/diffraction.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/settings.cpp -o $@
$(OBJD)/spaceGroup.o : $(SRCD)/spaceGroup.cpp $(SRCD)/spaceGroup.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/num.h $(SRCD)/iso.h $(SRCD)/symmetry.h $(SRCD)/pointGroup.h $(SRCD)/text.h $(SRCD)/list.h $(SRCD)/constants.h $(SRCD)/elements.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/spaceGroup.cpp -o $@
$(OBJD)/structureIO.o : $(SRCD)/structureIO.cpp $(SRCD)/structureIO.h $(SRCD)/mintStructure.h $(SRCD)/crystalMaker.h $(SRCD)/vasp.h $(SRCD)/findsym.h $(SRCD)/espresso.h $(SRCD)/json.h $(SRCD)/cif.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/iso.h $(SRCD)/text.h $(SRCD)/fileSystem.h $(SRCD)/num.h $(SRCD)/elements.h $(SRCD)/list.h $(SRCD)/constants.h $(SRCD)/kpoints.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/structureIO.cpp -o $@
$(OBJD)/symmetry.o : $(SRCD)/symmetry.cpp $(SRCD)/multi.h $(SRCD)/symmetry.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/constants.h $(SRCD)/text.h $(SRCD)/num.h $(SRCD)/list.h $(SRCD)/iso.h $(SRCD)/elements.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/symmetry.cpp -o $@
$(OBJD)/text.o : $(SRCD)/text.cpp $(SRCD)/text.h $(SRCD)/list.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/text.cpp -o $@
$(OBJD)/timer.o : $(SRCD)/timer.cpp $(SRCD)/timer.h $(SRCD)/text.h $(SRCD)/list.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/timer.cpp -o $@
$(OBJD)/unique.o : $(SRCD)/unique.cpp $(SRCD)/unique.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/iso.h $(SRCD)/symmetry.h $(SRCD)/text.h $(SRCD)/list.h $(SRCD)/num.h $(SRCD)/constants.h $(SRCD)/elements.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/unique.cpp -o $@
$(OBJD)/vasp.o : $(SRCD)/vasp.cpp $(SRCD)/multi.h $(SRCD)/num.h $(SRCD)/vasp.h $(SRCD)/kpoints.h $(SRCD)/structureIO.h $(SRCD)/language.h $(SRCD)/output.h $(SRCD)/text.h $(SRCD)/list.h $(SRCD)/constants.h $(SRCD)/iso.h $(SRCD)/elements.h $(SRCD)/fileSystem.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(SRCD)/vasp.cpp -o $@

# Linker
$(EXE) : $(OBJD)/mint.o $(OBJD)/multi.o $(OBJD)/output.o $(OBJD)/launcher.o $(OBJD)/text.o $(OBJD)/iso.o $(OBJD)/structureIO.o $(OBJD)/symmetry.o $(OBJD)/potential.o $(OBJD)/phonons.o $(OBJD)/kmc.o $(OBJD)/diffraction.o $(OBJD)/random.o $(OBJD)/elements.o $(OBJD)/constants.o $(OBJD)/fileSystem.o $(OBJD)/mtwist.o $(OBJD)/randistrs.o $(OBJD)/language.o $(OBJD)/settings.o $(OBJD)/about.o $(OBJD)/help.o $(OBJD)/randomStructure.o $(OBJD)/unique.o $(OBJD)/pointGroup.o $(OBJD)/spaceGroup.o $(OBJD)/interstitial.o $(OBJD)/gaPredict.o $(OBJD)/pdf.o $(OBJD)/timer.o $(OBJD)/mintStructure.o $(OBJD)/crystalMaker.o $(OBJD)/vasp.o $(OBJD)/findsym.o $(OBJD)/espresso.o $(OBJD)/json.o $(OBJD)/cif.o $(OBJD)/kpoints.o $(OBJD)/locPotential.o $(OBJD)/extPotential.o $(OBJD)/bonds.o $(OBJD)/pairPotential.o $(OBJD)/ewald.o $(OBJD)/relax.o
	$(CC) $(LINK) $(BLAS) $(LAPACK) $(OPT) $(OBJD)/mint.o $(OBJD)/multi.o $(OBJD)/output.o $(OBJD)/launcher.o $(OBJD)/text.o $(OBJD)/iso.o $(OBJD)/structureIO.o $(OBJD)/symmetry.o $(OBJD)/potential.o $(OBJD)/phonons.o $(OBJD)/kmc.o $(OBJD)/diffraction.o $(OBJD)/random.o $(OBJD)/elements.o $(OBJD)/constants.o $(OBJD)/fileSystem.o $(OBJD)/mtwist.o $(OBJD)/randistrs.o $(OBJD)/language.o $(OBJD)/settings.o $(OBJD)/about.o $(OBJD)/help.o $(OBJD)/randomStructure.o $(OBJD)/unique.o $(OBJD)/pointGroup.o $(OBJD)/spaceGroup.o $(OBJD)/interstitial.o $(OBJD)/gaPredict.o $(OBJD)/pdf.o $(OBJD)/timer.o $(OBJD)/mintStructure.o $(OBJD)/crystalMaker.o $(OBJD)/vasp.o $(OBJD)/findsym.o $(OBJD)/espresso.o $(OBJD)/json.o $(OBJD)/cif.o $(OBJD)/kpoints.o $(OBJD)/locPotential.o $(OBJD)/extPotential.o $(OBJD)/bonds.o $(OBJD)/pairPotential.o $(OBJD)/ewald.o $(OBJD)/relax.o -o $@

