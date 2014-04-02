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
VPATH := $(SRCD):dlib

# Build everything
all : touchAbout $(EXE)
touchAbout :
	@touch about.cpp

# Builds for object files
$(OBJD)/about.o : about.cpp about.h output.h num.h list.h text.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(FVERS) $(FTIME) about.cpp -o $@
$(OBJD)/bonds.o : bonds.cpp bonds.h iso.h elements.h num.h text.h list.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) bonds.cpp -o $@
$(OBJD)/cif.o : cif.cpp cif.h elements.h symmetry.h spaceGroup.h language.h output.h num.h iso.h text.h list.h fileSystem.h pointGroup.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) cif.cpp -o $@
$(OBJD)/constants.o : constants.cpp constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) constants.cpp -o $@
$(OBJD)/crystalMaker.o : crystalMaker.cpp crystalMaker.h bonds.h language.h num.h output.h iso.h text.h elements.h list.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) crystalMaker.cpp -o $@
$(OBJD)/diffraction.o : diffraction.cpp multi.h diffraction.h language.h output.h text.h num.h iso.h elements.h symmetry.h fileSystem.h list.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) diffraction.cpp -o $@
$(OBJD)/elements.o : elements.cpp elements.h output.h text.h num.h list.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) elements.cpp -o $@
$(OBJD)/espresso.o : espresso.cpp multi.h espresso.h structureIO.h language.h output.h text.h num.h iso.h elements.h kpoints.h list.h fileSystem.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) espresso.cpp -o $@
$(OBJD)/ewald.o : ewald.cpp multi.h num.h ewald.h language.h output.h text.h list.h constants.h locPotential.h iso.h elements.h symmetry.h potential.h fileSystem.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) ewald.cpp -o $@
$(OBJD)/extPotential.o : extPotential.cpp extPotential.h language.h output.h vasp.h espresso.h potential.h iso.h elements.h symmetry.h text.h num.h list.h fileSystem.h kpoints.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) extPotential.cpp -o $@
$(OBJD)/fileSystem.o : fileSystem.cpp multi.h fileSystem.h output.h text.h num.h list.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) fileSystem.cpp -o $@
$(OBJD)/findsym.o : findsym.cpp findsym.h num.h output.h iso.h text.h list.h constants.h elements.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) findsym.cpp -o $@
$(OBJD)/gaPredict.o : gaPredict.cpp gaPredict.h structureIO.h randomStructure.h fileSystem.h language.h output.h num.h ga.h iso.h symmetry.h potential.h diffraction.h random.h text.h list.h spaceGroup.h constants.h mtwist.h randistrs.h elements.h pointGroup.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) gaPredict.cpp -o $@
$(OBJD)/help.o : help.cpp help.h output.h text.h num.h list.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) help.cpp -o $@
$(OBJD)/interstitial.o : interstitial.cpp multi.h interstitial.h constants.h output.h text.h num.h iso.h symmetry.h list.h elements.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) interstitial.cpp -o $@
$(OBJD)/iso.o : iso.cpp multi.h iso.h output.h text.h num.h elements.h list.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) iso.cpp -o $@
$(OBJD)/json.o : json.cpp json.h bonds.h language.h output.h list.h iso.h text.h num.h elements.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) json.cpp -o $@
$(OBJD)/kmc.o : kmc.cpp kmc.h symmetry.h unique.h language.h num.h iso.h elements.h structureIO.h text.h random.h fileSystem.h constants.h output.h list.h mtwist.h randistrs.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) kmc.cpp -o $@
$(OBJD)/kpoints.o : kpoints.cpp kpoints.h output.h text.h num.h list.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) kpoints.cpp -o $@
$(OBJD)/language.o : language.cpp language.h text.h list.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) language.cpp -o $@
$(OBJD)/launcher.o : launcher.cpp multi.h num.h launcher.h settings.h about.h help.h randomStructure.h unique.h pointGroup.h spaceGroup.h interstitial.h gaPredict.h pdf.h fileSystem.h language.h timer.h output.h text.h list.h constants.h iso.h structureIO.h symmetry.h potential.h phonons.h kmc.h diffraction.h random.h ga.h elements.h mtwist.h randistrs.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) launcher.cpp -o $@
$(OBJD)/locPotential.o : locPotential.cpp locPotential.h pairPotential.h ewald.h relax.h output.h potential.h iso.h symmetry.h text.h num.h list.h elements.h constants.h fileSystem.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) locPotential.cpp -o $@
$(OBJD)/mint.o : mint.cpp multi.h output.h launcher.h text.h num.h list.h iso.h structureIO.h symmetry.h potential.h phonons.h kmc.h diffraction.h random.h elements.h constants.h fileSystem.h mtwist.h randistrs.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) mint.cpp -o $@
$(OBJD)/mintStructure.o : mintStructure.cpp num.h mintStructure.h elements.h symmetry.h spaceGroup.h language.h output.h list.h constants.h iso.h text.h fileSystem.h pointGroup.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) mintStructure.cpp -o $@
$(OBJD)/mtwist.o : mtwist.c mtwist.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) mtwist.c -o $@
$(OBJD)/multi.o : multi.cpp multi.h language.h text.h num.h list.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) $(FMPI) multi.cpp -o $@
$(OBJD)/output.o : output.cpp multi.h output.h text.h num.h list.h constants.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) output.cpp -o $@
$(OBJD)/pairPotential.o : pairPotential.cpp multi.h pairPotential.h language.h output.h text.h num.h list.h constants.h locPotential.h iso.h elements.h symmetry.h potential.h fileSystem.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) pairPotential.cpp -o $@
$(OBJD)/pdf.o : pdf.cpp pdf.h elements.h num.h language.h output.h iso.h diffraction.h text.h fileSystem.h list.h constants.h symmetry.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) pdf.cpp -o $@
$(OBJD)/phonons.o : phonons.cpp phonons.h language.h constants.h output.h iso.h symmetry.h potential.h text.h fileSystem.h num.h list.h elements.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) phonons.cpp -o $@
$(OBJD)/pointGroup.o : pointGroup.cpp pointGroup.h output.h num.h iso.h symmetry.h text.h list.h constants.h elements.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) pointGroup.cpp -o $@
$(OBJD)/potential.o : potential.cpp potential.h locPotential.h extPotential.h language.h num.h iso.h elements.h symmetry.h fileSystem.h list.h text.h output.h constants.h vasp.h espresso.h kpoints.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) potential.cpp -o $@
$(OBJD)/randistrs.o : randistrs.c mtwist.h randistrs.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) randistrs.c -o $@
$(OBJD)/random.o : random.cpp num.h multi.h random.h output.h constants.h list.h text.h mtwist.h randistrs.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) random.cpp -o $@
$(OBJD)/randomStructure.o : randomStructure.cpp randomStructure.h output.h constants.h num.h iso.h spaceGroup.h symmetry.h random.h list.h text.h elements.h pointGroup.h mtwist.h randistrs.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) randomStructure.cpp -o $@
$(OBJD)/relax.o : relax.cpp relax.h output.h iso.h symmetry.h locPotential.h text.h num.h list.h constants.h elements.h potential.h fileSystem.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) relax.cpp -o $@
$(OBJD)/settings.o : settings.cpp settings.h language.h fileSystem.h output.h iso.h structureIO.h ga.h gaPredict.h text.h list.h num.h constants.h elements.h random.h mtwist.h randistrs.h symmetry.h potential.h diffraction.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) settings.cpp -o $@
$(OBJD)/spaceGroup.o : spaceGroup.cpp spaceGroup.h language.h output.h num.h iso.h symmetry.h pointGroup.h text.h list.h constants.h elements.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) spaceGroup.cpp -o $@
$(OBJD)/structureIO.o : structureIO.cpp structureIO.h mintStructure.h crystalMaker.h vasp.h findsym.h espresso.h json.h cif.h language.h output.h iso.h text.h fileSystem.h num.h elements.h list.h constants.h kpoints.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) structureIO.cpp -o $@
$(OBJD)/symmetry.o : symmetry.cpp multi.h symmetry.h language.h output.h constants.h text.h num.h list.h iso.h elements.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) symmetry.cpp -o $@
$(OBJD)/text.o : text.cpp text.h list.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) text.cpp -o $@
$(OBJD)/timer.o : timer.cpp timer.h text.h list.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) timer.cpp -o $@
$(OBJD)/unique.o : unique.cpp unique.h language.h output.h iso.h symmetry.h text.h list.h num.h constants.h elements.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) unique.cpp -o $@
$(OBJD)/vasp.o : vasp.cpp multi.h num.h vasp.h kpoints.h structureIO.h language.h output.h text.h list.h constants.h iso.h elements.h fileSystem.h 
	@mkdir -p $(@D)
	$(CC) $(FALL) vasp.cpp -o $@

# Linker
$(EXE) : $(OBJD)/mint.o $(OBJD)/multi.o $(OBJD)/output.o $(OBJD)/launcher.o $(OBJD)/text.o $(OBJD)/iso.o $(OBJD)/structureIO.o $(OBJD)/symmetry.o $(OBJD)/potential.o $(OBJD)/phonons.o $(OBJD)/kmc.o $(OBJD)/diffraction.o $(OBJD)/random.o $(OBJD)/elements.o $(OBJD)/constants.o $(OBJD)/fileSystem.o $(OBJD)/mtwist.o $(OBJD)/randistrs.o $(OBJD)/language.o $(OBJD)/settings.o $(OBJD)/about.o $(OBJD)/help.o $(OBJD)/randomStructure.o $(OBJD)/unique.o $(OBJD)/pointGroup.o $(OBJD)/spaceGroup.o $(OBJD)/interstitial.o $(OBJD)/gaPredict.o $(OBJD)/pdf.o $(OBJD)/timer.o $(OBJD)/mintStructure.o $(OBJD)/crystalMaker.o $(OBJD)/vasp.o $(OBJD)/findsym.o $(OBJD)/espresso.o $(OBJD)/json.o $(OBJD)/cif.o $(OBJD)/kpoints.o $(OBJD)/locPotential.o $(OBJD)/extPotential.o $(OBJD)/bonds.o $(OBJD)/pairPotential.o $(OBJD)/ewald.o $(OBJD)/relax.o
	$(CC) $(LINK) $(BLAS) $(LAPACK) $(OPT) $(OBJD)/mint.o $(OBJD)/multi.o $(OBJD)/output.o $(OBJD)/launcher.o $(OBJD)/text.o $(OBJD)/iso.o $(OBJD)/structureIO.o $(OBJD)/symmetry.o $(OBJD)/potential.o $(OBJD)/phonons.o $(OBJD)/kmc.o $(OBJD)/diffraction.o $(OBJD)/random.o $(OBJD)/elements.o $(OBJD)/constants.o $(OBJD)/fileSystem.o $(OBJD)/mtwist.o $(OBJD)/randistrs.o $(OBJD)/language.o $(OBJD)/settings.o $(OBJD)/about.o $(OBJD)/help.o $(OBJD)/randomStructure.o $(OBJD)/unique.o $(OBJD)/pointGroup.o $(OBJD)/spaceGroup.o $(OBJD)/interstitial.o $(OBJD)/gaPredict.o $(OBJD)/pdf.o $(OBJD)/timer.o $(OBJD)/mintStructure.o $(OBJD)/crystalMaker.o $(OBJD)/vasp.o $(OBJD)/findsym.o $(OBJD)/espresso.o $(OBJD)/json.o $(OBJD)/cif.o $(OBJD)/kpoints.o $(OBJD)/locPotential.o $(OBJD)/extPotential.o $(OBJD)/bonds.o $(OBJD)/pairPotential.o $(OBJD)/ewald.o $(OBJD)/relax.o -o $@

# Clean
clean:
	rm -fr $(OBJD) $(EXE)
