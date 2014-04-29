# Mint (Materials Interface)


## Installation


1. Get a copy of "update.sh" onto your local machine
  - Method 1
    * Copy code from update.sh at https://raw.githubusercontent.com/materials/mint/master/update.sh to your local machine
    * Run on your local machine "chmod u+x update.sh"
  - Method 2
    * Download all files in the repository using the "Download ZIP" link at https://github.com/materials/mint
    * Unzip the files on your local machine - one of these files will be update.sh
2. Run on your local machine "./update.sh"
3. Edit the Makefile with your settings and run "make"
4. In the future, just run "./update.sh" and then "make" anytime that you want to update to a new version


## The Interface to Mint


### General

A single call to mint should have the form: 

    mint _files -_function1 _arguments1 -_function2 _arguments2

Here, _files is a list of files that will be read on startup. All files should appear before the first function call. Each function is preceded by a - sign. In general, only the first few letters are needed for a function name. A list of arguments needed for the function should appear after the call. Any number of functions can be called at a given time; with the exception of settings, they will be executed in the order in which they appear. Those functions that change settings (e.g. -tolerance) will be executed before all others.


Any number of files can be passed in a single call to mint. Their order does not matter and their formats will be determined automatically. The following file formats are supported as input:
- Structure (vasp 4, vasp 5, mint, cif, quantum espresso)
- Settings
- Potential
- Xray pattern (mint, pdf xml)

As a practical example, the call to

    mint structureFile -conventional -symmetry -print screen

would read the structure in structureFile, convert it to a conventional cell, print the symmetry of the conventional cell, and print the structure to stdout.


### Settings


Settings can be changed using a settings file supplied to mint using e.g.:

    mint _settingsFile _files -_function1 ...

where _settingsFile is the name of the settings file, _files are any other files, and -_function1 is the first function to perform. Multiple settings files can be supplied in a single call. For example, in the call to

    mint _settings1 _settings2 _files -_function1 ...

settings in _settings1 would be applied, followed by _settings2 (overwriting values in _settings1 if defined also in _settings2). Any settings applied on the command line (e.g. -tolerance) overwrite anything supplied by file.


A global settings file can also be created to change default values. This file should be named ".mint_settings" and appear in your home directory. Mint will automatically read this file if it is present. These settings have the lowest weight and are overwritten by any files or settings passed in on the command line.



### Output


General output will be directed to stdout. However, for some functions (e.g. optimizations) more detailed log files will be created in order to organize output more efficiently. General output is divided into two catagories: runtime and results.


Runtime output shows details of each function call as they occur. This can be useful be degbugging purposes and may give more detailed information than what is shown in the results. Each line from runtime output is preceeded by o:, w:, or e:, which coorespond to ordinary output, warnings, and errors, respectively. Warnings and errors are always shown, but ordinary output can be supressed since it often gives far more information than is needed (see the setting "displaylevel" or the function "-display").


The other output category, results, is always printed. Results are not preceded by any identifier (o, w, or e) and will be printed to stdout (with the exception of structures, which may be printed to a file - see the setting "usestdout" or the function "-print").


## Functions Available in Mint


List of functions in this document and brief descriptions:

            -help   Print help file
               -n   Set the number of processors in external mpi jobs
         -display   Set the level of runtime output to show
            -time   Print the time to complete a single call to mint
       -tolerance   Set the maximum cartesian distance for values to be equal
           -print   Control format and location of writing structures
            -name   Change the file name to which a structure is written
             -fix   Fix cell parameters and positions that are defined
          -remove   Remove atoms from the structure
       -neighbors   Print nearest neighbors lists for atoms in the structure
           -shell   Print nearest neighbor shells for atoms in the structure
       -transform   Apply a transformation matrix to convert between unit cells
          -rotate   Rotate all positions of atoms with fixed unit cell vectors
           -shift   Shift all positions of atoms by a set amount
         -reduced   Convert the current cell to its Niggli-reduced form
       -primitive   Convert the current cell to a primitive one
    -conventional   Convert the current cell to its conventional form
           -ideal   Convert to an "ideal" form of the current cell
        -symmetry   Print the symmetry of the current unit cell
           -point   Get point group of structure or print about point groups
           -space   Get space group of structure or print about space groups
           -about   Print listed information about the structure or atoms in it
          -refine   Adjust cell to satisfy symmetry to numerical precision
         -perturb   Randomly adjust unit cell parameters or positions
          -unique   Print unique atoms, groups of atoms, or jumps between atoms
      -equivalent   Same as -unique but prints all atoms/groups/jumps
    -interstitial   Generate interstitial sites in the structure
          -energy   Calculate energy of the structure under supplied potential
          -forces   Calculate forces on all atoms under supplied potential
     -diffraction   Calculate diffraction patterns and R factors
        -optimize   Global optimization of the structure to minimize energy
         -compare   Compare two structures to determine if they are similar
             -kmc   Run kinetic Monte Carlo (KMC) diffusion simulation



### -help

General: Print help information.



### -n / -np

######General:
Set the number of processors to use when launching an external mpi job. This is only used if the launching function (MPIRUN in the makefile) contains -n or -np. Otherwise it is ignored since the launching function determines the number of available processors on its own.

######Arguments: 
- Any integer number to set the number of processors to use

######Default: 
1

######Examples:
    "-np 2" to launch external mpi programs on 2 processors



### -display / -output

######General: 
Set the level of runtime display to show. Level 1 is the most broad
    while the detail increases with level. This controls ordinary output only; 
    warnings and errors are always printed.

######Arguments: 
- Any integer number to set the max level that is displayed (0 is off)
- "all" or "everything" to show all output
- "none", "nothing", or "off" to show no output

######Default: 
All runtime display is shown

######Examples:
    "-display 3"   show the first three levels of runtime output
    "-display off" show no output (except warnings/errors)



### -time

######General: 
Show the total time to complete all functions.

######Arguments: 
- Any integer number sets the number of places to show after decimal
- "seconds" to print as seconds only, not hours/minutes/seconds

###### Default: 
- four places after the decimal
- Time is formatted into hours/minutes/seconds

######Examples:
    "-time"     print the time to finish the run
    "-time 2"   print the time with two places after decimal
    "-time sec" print, for example, 90 seconds, not 1 minute, 30 seconds



### -tolerance

######General: 
Set the absolute tolerance used during comparisons between positions
    and distances. For a tolerance _tol, two positions are taken as equal if
    the distance between them, _dist, satisfies _dist <= _tol. This tolerance
    is used in many functions, such as those for determining symmetries.

######Arguments: 
- Any number (units: angstroms)

######Default: 
1e-4

######Examples:
    "-tolerance 1e-2" set the tolerance to 1e-2 angstroms



### -print / -write

######General: 
Controls where structures are printed and their format. If some
    positions or the lattice vectors are not defined, then they will be
    assigned random values; this can be turned off when printing to the mint
    format by passing "free" as an argument.

######Arguments:
- "fractional" or "direct" to print using fractional coordinates
- "cartesian" to print to cartesian coordinates if available in set format
- "screen" or "stdout" to print to standard out instead of a file
- "mint" or ".mint" for mint format
- "vasp", "vasp5", "vasp 5", ".vasp", or ".vasp5" for vasp 5 format
- "vasp 4", "vasp4", or ".vasp4" for vasp 4 format
- "cif" or ".cif" for cif format
- "crystalmaker", "cm", "cmtx", or ".cmtx" for crystal maker format
- "findsym", ".findsym", "fs", or ".fs" for findsym format
- "quantum", "espresso", "qe", or ".qe" for quantum espresso format
- "free" prevent atoms/basis from being randomly set (only for mint format)

######Default: 
structures are printed to a file in the same format as the original

######Examples:
    "-print screen"         print to standard out
    "-print vasp"           print to vasp format
    "-print vasp cartesian" print to vasp format in cartesian coordinates



### -name

######General: 
Change the file name for writing a structure.

######Arguments: 
- Any valid file name

######Default: 
- Same file name as the original. If the setting "overwrite" is set to false, then a number may be appended to the file name (i.e. name_1) in  order to avoid writing over an existing file. If this setting is true, then overwrites will be allowed.

######Examples:
    "-name newname" print structure to a file named newname



### -remove

######General: 
Remove atoms from the structure by element type, index, or position.

######Arguments:
- Integer to remove atom with that index
- Element to remove all atoms of set element
- Element and integer to remove instance of set element
- Three decimal values to remove atom at set position

######Examples:
    "-rem 1 2"         remove atoms 1 and 2
    "-rem Al"          remove all Al atoms
    "-rem Al 2 H 3"    remove second Al and third H atoms
    "-rem 0.5 0.0 0.0" remove atom at 0.5, 0.0, 0.0



### -fix

######General: 
Fix atomic coordinates and lattice parameters that have been
    explicitly defined for the structure and are not already fixed. Any
    positions or parameters that were not set will still be allowed to change.
    This is useful when importing a structure that is to be optimized but for
    which some properties are known; in this case, the properties that are
    fixed by this command will not be changed.

######Arguments:
- "atoms" or "positions" to fix positions
- "basis" or "lattice" to fix lattice parameters

######Default: 
Without arguments, known positions and lattice parameters are fixed

######Examples:
    "-fix"       fix all positions and lattice parameters that are defined
    "-fix atoms" fix all positions that are defined, basis is still free
    "-fix basis" fix lattice parameters, all positions are still free



### -neighbors

######General: 
Calculate nearest neighbors list. For each atom, the distance in
    angstroms is printed as well as three integers that give the cell in which
    the nearest image of the atom occurs.

######Arguments:
- No arguments to print neighbors list for all atoms
- Integer to print neighbors to atom with that index
- Element to print neighbors to all atoms of set element
- Element and integer to print neighbors to instance of set element
- Three decimal values to print neighbors to atom at set position

######Examples:
    "-neigh"             print neighbors list for all atoms
    "-neigh 1 2"         print neighbors for atoms 1 and 2
    "-neigh Al"          print neighbors for all Al atoms
    "-neigh Al 2 H 3"    print neighbors for second Al and third H atoms
    "-neigh 0.5 0.0 0.0" print neighbors for atom at 0.5, 0.0, 0.0



### -shell

######General: 
Calculate nearest neighbor shells. Given a particular atom, A0, in the
    structure, all other other atoms (including periodic images) out to a set
    maximum are grouped by distance from A0. When grouping atoms into shells,
    two distances are considered the same if the difference between them is
    less than some tolerance; this can be changed using the -tol command or the
    "tolerance" tag in a settings file.

######Arguments:
- "details" to print information about each atom in the shell
- No arguments to print shells for all atoms
- Integer to print shells of atom with that index
- Element to print shells of all atoms of set element
- Element and integer to print shells of instance of set element
- Decimal value to set the max distance to calculate shells

######Default: 
Default max distance is 3 Ang

######Examples:
    "-shell"             print shells for all atoms
    "-shell 1 2"         print shells for atoms 1 and 2
    "-shell Al"          print shells for all Al atoms
    "-shell Al 2 H 3"    print shells for second Al and third H atoms
    "-shell 5.0"         print shells out to 5 Ang
    "-shell det"         print shells and information about each member atom



### -transform

######General: 
Transform a cell. The new lattice vectors are equal to Lnew = M*Lorig
    where M is the transformation matrix and Lorig is the original matrix of
    lattice vectors (each vector along a row of the matrix).

######Arguments:
- Three numbers (n1 to n3) to set the transformation matrix  
			| n1  0  0 |  
		M = | 0  n2  0 |  
			| 0   0 n3 |  
- Nine numbers (n1 to n9) to set the transformation matrix  
			| n1 n2 n3 |  
		M = | n4 n5 n6 |  
			| n7 n8 n9 |

######Examples:
    "-transform 2 1 2"             transform with 2, 1, 2 on matrix diagonal
    "-transform 1 1 1 0 1 1 0 0 1" transform with corresponding matrix



### -rotate

General: Rotate fractional positions of atoms in the cell. For a rotation
    matrix R, the position of an atom after rotation is R*x0 where x0 is its
    original position.

Arguments:
    Nine numbers (r1 to r9) to set the rotation matrix
            | r1 r2 r3 |
        R = | r4 r5 r6 |
            | r7 r8 r9 |

Examples:
    "-rotate 1 1 0 1 0 0 0 0 1" rotate positions with corresponding matrix



### -shift

General: Shift atom positions by a set amount.

Arguments: (all numbers are given as fractions of lattice vectors)
    One number (n) to shift by the vector n, n, n
    Three numbers (n1 to n3) to shift by the vector n1, n2, n3

Examples:
    "-shift 0.5"     shift positions by the vector 0.5, 0.5, 0.5
    "-shift 0.5 0 0" shift positions by the vector 0.5, 0.0, 0.0



### -reduced

General: Perform Niggli reduction of current cell. Note that this DOES NOT
    convert to the primitive cell first; for a reduced primitive cell, run
    -primitive -reduced.

Arguments: none



### -primitive

General: Convert to the primitive form of the current cell. Note that the
    primitive cell is not unique; a unique primitive cell can be obtained by
    reducing (-reduced) the primitive cell.

Arguments: none



### -conventional

General: Convert current structure to conventional cell. The space group of the
    structure is first determined and the appropriate transformation and origin
    shift are applied to convert the structure so that the symmetry operations
    agree with those in the International Tables of Crystallography.

Arguments: none



### -ideal

General: Convert current structure to its most ideal form. The most ideal form
    is defined as the unit cell that maximizes the minimum image distance under
    periodic boundary conditions while preserving a set number of atoms or,
    equivalently, a set volume. The transformation that is found also ensures
    that all rotational symmetry operation elements are integer. This last
    requirement can significantly restrict the cell shapes that are allowed
    for high-symmetry structures, which may lead to cell sizes that are very
    different than what is targeted.

Arguments:
    None to convert current cell to its most ideal form without changing size
    Integer number to set the minimum number of atoms in the ideal cell
    "max" or "most" to change to setting the maximum number of atoms
    "dis" and a number (Ang) to set ideal cell by the minimum image distance

Examples:
    "-ideal"         to convert to ideal cell without changing size
    "-ideal 100"     to convert to ideal cell with at least 100 atoms
    "-ideal max 100" to convert to ideal cell with at most 100 atoms
    "-ideal dis 12"  to convert to cell with minimum image distance of 12 Ang



### -symmetry

General: Print the symmetry of the current unit cell. A list of unique symmetry
    operations are printed first; there should be as many operations in this
    list as there are operations in the primitive cell of the structure. A
    second list is printed of centering vectors. For a primitive cell, this
    will simply contain 0, 0, 0, but can contain more entries for non-primitive
    cells. Each centering vector can be added to any of the unique symmetry
    operations so that in total there are Nsym*Ncent symmetry operations in the
    unit cell where Nsym is the number of unique operations and Ncent is the
    number of centering vectors. Note that these symmetries apply the current
    cell, which may not necessarily be the conventional one. To obtain the
    conventional cell symmetry, run -conventional -symmetry.

Arguments:
    "matrix" to print operations in matrix form

Defaults: Without arguments operations are printed in Jones-Faithful notation

Expample:
    "-symmetry"        to print operations in Jones-Faithful notation
    "-symmetry matrix" to print operations in matrix form



### -point

General: Print the point group of the structure, a list of all point groups,
    or information about a point group.

Arguments:
    No arguments and no supplied structure to print a list of all point groups
    No arguments and a structure to print the point group of the structure
    Point group name (no spaces) for information about that point group

Examples:
    "mint -point"               print a list of all point groups
    "mint structureFile -point" get point group of structure in structureFile
    "mint -point mmm 6"         get information about point groups mmm and 6



### -space

General: Print the space group of the structure, a list of all space groups,
    or information about a space group.

Arguments:
    No arguments and no supplied structure to print a list of all space groups
    No arguments and a structure to print the space group of the structure
    Space group name (no spaces) for information about that space group

Examples:
    "mint -space"               print a list of all space groups
    "mint structureFile -space" get space group of structure in structureFile
    "mint -space P3 C2"         get information about space groups P3 and C2



### -about

General: Print the following information about the structure:
      * Point group
      * Space group
      * Unit cell volume
      * Unit cell basis vectors
      * Unit cell lengths and angles
      * Symmetry operations
      * All atomic positions grouped by equivalence under structure symmetry

    If the keyword "atom" is passed, then the following information for each
    selected atom is printed:
      * Atom number and element
      * Fractional and cartesian coordinates
      * Site symmetry operations
      * List of atoms in the same orbit

Arguments:
    No arguments to print information about the structure
    "atom" to print information about selected atoms
    Integer to print information about atom with that index
    Element to print information about all atoms of set element
    Element and integer to print information about instance of set element
    Three decimal values to print information about atom at set position

Examples:
    "-about"                  print information about the structure
    "-about atom 1 2"         print information about atoms 1 and 2
    "-about atom Al"          print information about all Al atoms
    "-about atom Al 2 H 3"    print information about 2nd Al and 3rd H atoms
    "-about atom 0.5 0.0 0.0" print information about atom at 0.5, 0.0, 0.0



### -refine

General: Refine the lattice vectors and atomic positions so that a set of
    symmetry operations that are satisfied only to some tolerance are made
    exact. A wide tolerance (e.g. 0.1 angstroms or higher, set with -tolerance)
    may be needed for structures that are slightly disordered. Once the
    symmetry has been determined, the lattice vectors are refined to preserve
    the lattice symmetry that is found (e.g. cubic, tetragonal, etc.). Atomic
    positions are refined essentially by averaging over sites that are
    equivalent to one another under these symmetries.

    If a diffraction pattern has been supplied as input, the refinement
    has a different purpose. In this case, the atomic coordinates are relaxed
    so as to minimize the difference between the calculated and supplied
    patterns.

Arguments:
    "atoms" or "positions" to refine positions only
    "basis" or "lattice" to refine lattice parameters only

Default: Without arguments, positions and lattice parameters are both refined

Examples:
    "-refine"       refine atomic positions and lattice vectors
    "-refine atoms" refine atomic positions only
    "-refine basis" refine lattice vectors only



### -perturb

General: Randomly adjust the atomic positions and/or lattice vectors.

Arguments:
    One number to perturb by this amount in random direction (units: angstroms)
    Two numbers for random amount within their range (units: angstroms)
    "basis" or "lattice" to perturb lattice vectors only
    "atoms", "coordinates", or "positions" to perturb positions only
    "symmetry" to preserve any symmetries that exist in the structure

Defaults: displacement amount of 0.25 angstroms

Examples:
    "-perturb"          perturb atoms and basis vectors by default
    "-perturb atoms"    perturb atoms only
    "-perturb basis"    perturb lattice vectors only
    "-perturb 0.1"      perturb atoms and basis by 0.1 angstroms
    "-perturb 0.1 0.2"  perturb atoms and basis randomly 0.1 to 0.2 angstroms
    "-perturb symmetry" make perturbations while maintaining cell symmetry



### -unique / -equivalent

General: Search for unique atoms, groups of atoms, or transitions between atomic
    sites in a structure. -unique prints only the unique groups that are found.
    -equivalent prints the unique groups and all equivalent groups for each.
    If searching for unique atoms or groups (not transitions), a copy is made
    of the structure for each unique group and the corresponding atoms are
    removed to make a vacancy. For transitions, no distance is needed but a
    max distance can be supplied. The program will use several criteria to
    determine those jumps that are likely the most important, but none of them
    will be larger than this max distance.

Arguments:
    List of elements in the group
    Number to set the distance from the first atom when searching for groups
    "Jumps" to search for jumps between sites
    "Interstitial" to search for atoms that are marked as interstitials

Defaults:
    Distance from first atom for a group is 2 Ang
    Max distance for a jump is 7 Ang

Examples:
    "-unique Al"       find all unique Al atoms in a structure
    "-equiv Al"        find unique Al atoms all equivalent atoms to each
    "-unique Al H 2"   find all Al-H pairs within 2 ang
    "-unique Al H H 2" find Al-H-H with neither H more than 2 ang from Al
    "-unique Al jumps" find all jumps between Al sites
    "-unique Al int"   find all unique Al atoms marked as interstitials



### -interstitial

General: Search for interstitial sites in a structure. This algorithm works by
    placing an exponential decay (Exp[-r/a]) at each atomic site and searching
    for minima in the resulting function. For each unique site that is found, a
    copy of the structure is made and an atom of the set element is placed at
    one of the sites. These structures can be printed by specifying the -print
    command. These sites are purely geometric and should be refined using a
    more accurate metric (e.g. pair potential or DFT).

Arguments:
    Element (symbol or name) to set the element to place at unique sites
    Integer to set number of start points (per atom) when searching for minima
    Float to set decay rate for exponentials (a in Exp[-r/a])
    "Expand" to generate all symmetrically equivalent sites

Defaults:
    Hydrogen is placed at unique sites
    25 starting points per atom
    Decay rate of 0.25

Examples:
    "-int"           search for interstitial sites
    "-int expand"    search for interstitial sites and populate all equivalent
    "-int Li"        place Li at interstitial sites
    "-int 50"        search for interstitial sites with 50 starting points
    "-int 0.2"       search for interstitial sites with decay rate of 0.2
    "-int 0.2 Li 50" previous three examples combined



### -energy

General: Calculate the energy of a structure under a supplied potential. Both
    the total energy of the unit cell and the energy per atom are returned.

Arguments: none



### -forces

General: Calculate the forces on atoms in a structure under a supplied
    potential.

Arguments: none



### -diffraction

General: Calculate the powder xray diffraction pattern for a structure. If a
    diffraction pattern was supplied as input, the r-factor will be calculated
    by optimizing the scaling between calculated and reference patterns, as well
    as the atomic displacement factors for the atoms in the structure (see
    -refine to optimize atomic positions).

Arguments: additional xray patterns to compare to reference
    "broaden" to print a broadened diffraction pattern
    "wavelength" and number to set the xray wavelength
    "fwhm" and number to set the FWHM for broadening (or use variance)
    "variance" and number to set the variance for broadending (or use fwhm)
    "resolution" and number to set resolution for printing broadened points
    "minimum" and number to set the minimum two-theta value to calculate
    "maximum" and number to set the maximum two-theta value to calculate

Examples:
    "mint str -xray"        get diffraction pattern for structure in str
    "mint x.in str -xray"   get diffraction pattern and compare to x.in
    "mint x.in -xray x2.in" compare patterns in x.in and x2.in
    "-xray broaden"         print a broadened pattern
    "-xray wave 1.5"        set the xray wavelength to 1.5 Ang
    "-xray fwhm 0.5"        set the fwhm to 0.5 deg for broadending
    "-xray res 0.05"        print a point every 0.05 deg when broadending
    "-xray min 10"          calculate intensities starting at 2theta of 10
    "-xray max 100"         calculate intensities up to 2theta of 100



### -optimize

General: Run global optimization for a structure. Any free parameters will be
    changed while any that are fixed will be preserved. Multiple optimization
    metrics can be supplied, and all will be used to rank structures during
    crossover operations, but only one of these metrics is used to determine
    which structure is the "best" one. For example, if a potential energy
    function and reference x-ray pattern are both supplied, but the energy is
    set as the metric to optimize, then both metrics are used to rank structures
    during crossover, but the structure with the overall lowest energy is taken
    as the "best" structure.

Arguments:
    "energy" or "potential" to optimize the energy
    "xray" or "diffraction" to optimize the r-factor

Examples:
    "mint str xray.in -opt"             optimize r-factor
    "mint str pot.in -opt"              optimize energy
    "mint str xray.in pot.in -opt xray" optimize r-factor, bias with energy
    "mint str xray.in pot.in -opt pot"  optimize energy, bias with r-factor



### -compare

General: Compare structures to determine if they are cells derived from the
    same lattice. Comparisons are made in the following ways:

        Test 1: Are elements the same in both structures?
        Action: Convert both cells to primitive form
        Test 2: Are the number of atoms of each element the same in both
                primitive cells?
        Test 3: (optional) Are the primitive cell volumes the same?
        Action: Convert both primitive cells to reduced form
        Test 4: Are the internal angles of both reduced cells the same? If not,
                convert reduced cell two from all acute angles to all obtuse
                angles (or vice versa) and compare angles again.
        Test 5: Are all ratios of reduced cell vector lengths the same? eg Is
                it true that a1/b1 = a2/b2 where a1 is the length of the a
                lattice vector in reduced cell 1?
        Test 6: Is there a rotation and translation that when applied to
                reduced cell 2, maps the positions onto those of reduced
                cell 1? Possible rotations are determined using a modified
                version of the Le Page algorithm to determine lattice symmetry.

    Two structures are considered the same only if they pass all tests. It may
    be necessary to increase the tolerance (-tol) to identify similar structures
    that are slightly different.

Arguments:
    "volume" to compare the volumes of the primitive cells

Default: Volumes are not compared

Examples:
    "mint str1 str2 -compare"     compare structures in str1 and str2
    "mint str1 str2 -compare vol" compare structures + volume in str1 and str2



### -kmc

General: Run lattice-based kinetic Monte Carlo (KMC) simulation of atomic
    diffusion. These simulations are initialized over two steps and run in
    a third. In the first step, the sub-lattice on which diffusion occurs is
    defined and the user performs relaxations of defects at symmetrically
    unique sites. In the second step, a network of jumps between sites in the
    structure is generated and the user calculates the energy along each path.
    In the final step, the KMC simulation is performed and the diffusivity
    prefactor and activation energy are calculated assuming Arrhenius behavior
    of the diffusivity.

    In the first step, the user supplies a structure and declares on which
    sub-lattice diffusion takes place. The program will then print a file called
    "kmc.setup" and a series of directories. The first several lines of this
    file will begin with the word "SITE". Each site corresponds to one of the
    directories than was created. The user should fully relax the structures in
    each of these directories. Then edit the kmc.setup file so that each site
    line points to these unrelaxed and relaxed structures.

    In the second step, the user supplies their modified kmc.setup file and the
    program generates a second file, "kmc.in", and another set of directories.
    Following the "SITE" lines in kmc.in are a series of lines that begin with
    "JUMP". Each jump line corresponds to one of the directories that was
    printed in this step. Each of the jump directories contains two files:
    "initial" and "final". The user should find the energetic barrier
    between these two structures (for example using NEB calculations). The user
    should then modify the kmc.in file. On each "SITE" line, the user enters
    the formation energy of the defect at that particular site, the total energy
    of the structure, and, optionally, the name of a file containing phonon
    frequencies for this structure. On each "JUMP" line, the user should enter
    the total energy of the structure in its transition state along this jump
    and, optionally, the name of a file containing phonon frequencies of the
    structure in its transition state.

    In the final step, the user supplies their modified kmc.in file and the
    program runs the KMC simulation. Simulations are performed at a series of
    temperatures and a fit of the Arrhenius equation (D = D0*exp[-Q/kt]) is
    applied to the resulting diffusivities; the values of the pre-exponential
    factor, D0, and the activation energy, Q, are printed to standard out.
    Diffusivities from each simulation are printed to files.

    In order to run diffusion for interstitial atoms, the initial structure in
    step 1 must be given in mint structure format, with the interstitial sites
    labeled as such.

Arguments:
    Step 1: Element to specify the sub-lattice on which diffusion is calculated
            "interstitial" to specify that diffusion is on intersitial sites
            One of the formats in "-print" to set the output file format
    Step 2: None
    Step 3: None

Default:
    Step 1: Diffusion is on parent structure sites (ie not interstitials)
            Structures are printed in the same format as the input one
    Step 2: N/A
    Step 3: N/A

Examples:
    Step 1: "mint str -kmc Al"      Set up self-diffusion on Al sites
            "mint str -kmc Al int"  Set up diffusion on interstitial Al sites
            "mint str -kmc Al vasp" Print all structure files in vasp format
    Step 2: "mint kmc.setup -kmc"   Generate kmc.in file from kmc.setup
    Step 3: "mint kmc.in -kmc"      Perform kmc simulation using kmc.in


## Settings

List of settings and brief descriptions (detailed descriptions follow):

             numprocs   Number of processors when launching external MPI programs
         displaylevel   Level of runtime output to display (0 is no output)
           displaytab   Indentation for each nested level of runtime output
             timeshow   Set whether the time to complete a run is printed
             timeprec   Number of decimal places to show when printing run time
           formattime   Control whether time is converted to minutes/seconds
            tolerance   Tolerance used to test floating point values for equality
           clustertol   Equality tolerance when expanding positions by symmetry
		    usestdout   Structures are printed to stdout, instead of to a file
          coordinates   Print positions in fractional or cartesian coordinates
            strformat   Control the file format for printing structures
            overwrite   Allow structure files to be overwritten
         addextension   Add a file extension when printing a structure to file
	  randstrmaxloops   Maximum number of loops when generating a random structure
       randstrminbond   Minimum random bond length
          gaoptnumsim   Number of unique runs during GA optimization
         gaoptpopsize   Number of structure in GA optimization
     gaoptcellmutprob   Probability of cell mutation in GA optimization
      gaoptposmutprob   Probability of position mutation in GA optimization
     gaoptwyckmutprob   Probability of Wyckoff site mutation in GA optimization
          gaoptmetric   Metric to optimize in GA-based optimization
        gaoptconverge   Number of GA generations without change for convergence
         gaoptmaxgens   Maximum number of GA generations allowed during optimization
       gaoptnumtokeep   Number of structures to keep between GA generations
       gaoptselection   Selection method to use during GA optimization
       gaoptenergytol   Tolerance for an energy to be a new best during GA
         gaoptdifftol   Tolerance for an r-factor to be a new best during GA
    gaoptscreenmethod   Screening method used during GA optimization
       gaoptscreennum   Number of trial structures to screen during GA
          wyckoffbias   Biasing level for choosing random Wyckoff positions
          minimagedis   Minimum image distance when generating supercells
      maxjumpdistance   Maximum jump distance when generating jumps between sites
      kmcjumpsperatom   Number of jumps in a KMC simulation per atom
          kmcconverge   Convergence for a KMC simulation to be complete



### numprocs

General: Mint has the ability to launch external programs such as VASP and
    Quantum Espresso, which may be MPI-enabled. This setting controls the number
    of processors that are allocated to the external job. In this way, Mint can
    be launched on a single core, but still call external programs running on
    many processors. This setting can be controlled from the command line
    using the -n/-np function.

Values: Any integer

Default: 1



### displaylevel

General: By default, no runtime output is shown, only results. However, there
    is a large amount of information that can be printed by increasing the
    value of this setting. Level 1 output is very general where higher levels
    become more specific. This setting can be controlled from the command line
    with the -display function.

Values: Any integer number

Default: 0



### displaytab

General: Each successive level of runtime output is indented by a set tab value.
    This setting controls the number of spaces that define this tab. This is a
    purely aesthetic setting.

Values: Any integer number

Default: 4



### timeshow

General: The time to complete a single call to Mint can be displayed at the
    end if a run with this setting. This can also be controlled from the command
    line using the -time function.

Values: True (to show the time)
        False (not to show the time)

Default: False (time is not shown)



### timeprec

General: Precision to use when printing the time to complete a run. This setting
    controls the number of places that are shown after the decimal. This can
    also be controlled from the command line using the -time function.

Values: Any integer number

Default: 4



### formattime

General: Control whether time is printed as the total number of seconds to
    complete a single run, or formatted into days/hours/minutes/seconds.

Values: True (to convert to days/hours/minutes/seconds)
        False (to print as the total number of seconds)

Default: True (time is converted to days/hours/minutes/seconds)



### tolerance

General: Global tolerance used to test floating point values for equality. Two
    numbers are taken as equal if the absolute difference between them is less
    than this value. Vectors are taken as equal if the length of the vector that
    connects their tips is less than this value. All values are converted to the
    cartesian coordinate system before testing for equality. This setting can
    be controlled from the command line using the -tol function.

Values: Any floating point number (in units of Angstroms where applicable)

Default: 1e-4



### clustertol

General: In several functions, Mint expands positions of atoms from the
    assymetric unit using the symmetries of the crystal. For atoms not lying on
    a general position, a single coordinate can be generated more than one time.
    This setting controls the tolerance for atoms generated using symmetry to be
    taken as equal. That is, if the distance between two atoms generated by
    symmetry is less than this value, then they are considered the same atom
    and their positions averaged.

Values: Any floating point number (in units of Angstroms)

Default: 0.2



### usestdout

General: Control whether structures are printed to stdout or to a file. This
    can be controlled from the command line using the -print function.

Values: True (print structures to stdout)
        False (print structures to file)

Default: False (structures are printed to file)



### coordinates

General: Control whether positions in a crystal are printed using fractional
    (direct) or cartesian coordinates. This can be controlled from the command
    line using the -print function. Some structure formats only allow
    fractional coordinates (but not cartesian), while others the opposite. For
    such formats, this setting is ignored and the correct coordinates are used.

Values: "Fractional" or "cartesian"

Default: Fractional



### strformat

General: Set the format used when printing a structure. This setting can be
    controlled from the command line using the -print function.

Values:
    "mint" or ".mint" for mint
    "vasp", "vasp5", "vasp 5", ".vasp", or ".vasp5" for vasp 5
    "vasp 4", "vasp4", or ".vasp4" for vasp 4
    "cif" or ".cif" for cif
    "crystalmaker", "cm", "cmtx", or ".cmtx" for crystal maker
    "findsym", ".findsym", "fs", or ".fs" for findsym
    "quantum", "espresso", "qe", or ".qe" for quantum espresso

Default: Same as input file



### overwrite

General: Control whether structure files can be overwritten. Mint names output
    files based on the name of the supplied file, but can alter the name so
    that it does not overwrite the original when printed. For example, if a
    structure is supplied in the file STR, then when printed, Mint will name
    the file STR_1 (or the first STR_N that is not taken, with N as an integer).

Values: True (to allow original files to be overwritten)
        False (to create a unique name each time that a file is printed)

Default: False (a unique name is generated each time that a file is printed)



### addextension

General: Add an extension to structure file names when printing them. The
    extension used will depened on the format of the file. For example, files
    printed to the VASP format will use the .vasp extension.

Values: True (add an extension to structure file names)
        False (do not add an extension to structure file names)

Default: True (extensions are added to structure file names)



### randstrmaxloops

General: Control the number of attempts to create a "good" structure when
    generated randomly. A structure is considered "good" if no two atoms
    are separated by a distance less than the minimum desired bond length (see
    the randstrminbond setting). If randstrmaxloops loops are made and no
    "good" structure has been generated, then the best structure up to that
    point is taken.

Values: Any integer number

Default: 100



### randstrminbond

General: Control the minimum bond length when generating random structures.
    The ideal bond length between two atoms, d_ideal, is taken as the sum of
    their average covalent radii (http://en.wikipedia.org/wiki/Covalent_radius).
    The minimum bond length between two atoms is then defined as
    d_min = randstrminbond * d_ideal. In other words, randstrminbond sets the
    minimum fraction of the ideal bond length that is allowed. This setting is
    tied to the randstrmaxloops setting.

Values: Any floating point number

Default: 0.5 (minimum distance is half the ideal distance)



### gaoptnumsim

General: Control the number of unique simulations to perform during GA-based
    structure prediction/solution. For each unique simulation, a new random
    initial population is generated. The best structure is taken as the one
    that optimizes the set metric (see gaoptmetric setting) best over all
    simulations. In other words, this controls the number of restarts that are
    performed during GA-based structure prediction.

Values: Any integer number

Default: 1 (only one unique simulation is performed)



### gaoptpopsize

General: Set the number of structures in the population during GA-based
    structure prediction/solution. This number is kept constant throughout the
    simulation.

Values: Any integer number

Default: 10



### gaoptcellmutprob

General: Probability of making a change to the unit cell parameters during
    GA-based structure prediction/solution. The unit cell is modified with this
    probability for each child generated during crossover at every generation
    of the GA. Modifications include changes the cell lengths and internal
    angles (where allowed by symmetry). This setting is ignored if the unit
    cell parameters are fixed for the structure.

Values: Floating point number between 0 (no mutations) and 1 (always mutate)

Default: 0.1 (probability of a mutation is 1/10)



### gaoptposmutprob

General: Probability of making a change to the positions of atoms in the unit
    cell during GA-based structure prediction/solution. Positions in a cell are
    modified with this probability for each child generated during crossover at
    every generation of the GA. All positions in the cell are changed along a
    different random vector by a distance randomly selected on the interval
    from 0.1 to 0.5 Angstroms. Symmetry is always preserved during mutation
    and atoms that are fixed are not moved.

Values: Floating point number between 0 (no mutations) and 1 (always mutate)

Default: 0.1 (probability of a mutation is 1/10)



### gaoptwyckmutprob

General: Probability of changing which Wyckoff sites are occupied in the
    structure during GA-based structure prediction/solution. Wyckoff site
    occupations are modified with this probability for each child generated
    during crossover at every generation of the GA. When an atom changes Wyckoff
    sites, the position of the atom is randomly generated so that it obeys its
    new site symmetry.

Values: Floating point number between 0 (no mutations) and 1 (always mutate)

Default: 0.1 (probability of a mutation is 1/10)



### gaoptmetric

General: Metric to optimize during GA-based structure prediction/solution. While
    multiple metrics can be used to guide the search, this setting controls
    which of them is used to determine the "best" structure. For example, an
    optimization supplied with a potential function and diffraction pattern to
    match will use both of these during optimization, but the "best" structure
    will be the one that has the lowest energy or lowest R-factor (depending on
    which metric is passed to this setting). This can also be controlled from
    the command line using the -optimize function.

Values:
    "energy" or "potential" to optimize the energy
    "xray" or "diffraction" to optimize the r-factor

Default: Energy if it is supplied, otherwise r-factor



### gaoptconverge

General: Number of generations without a change for a GA-based structure
    prediction/solution calculation to be considered converged. In other words,
    if no new "best" structure is found after gaoptconverge generations, the
    simulation terminates.

Values: Any integer number

Default: 10



### gaoptmaxgens

General: Maximum number of generations allowed during GA-based structure
    prediction/solution calculation. If the simulation has not converged after
    this many generations, then it terminates and the best structure to that
    point is returned.

Values: Any integer number

Default: 1000



### gaoptnumtokeep

General: Number of structures to keep from one generation to the next during
    a GA-based structure prediction/solution calculation. These structures are
    unchanged going to the next generation (i.e. no mutations are applied).
    Only the best gaoptnumtokeep structures are kept from the generation.

Values: Any integer number

Default: 0



### gaoptselection

General: Selection method for choosing structures to mate during GA-based
    structure prediction/solution calculations.

Values: "Roulette" or "tournament"

Default: Tournament



### gaoptenergytol

General: Tolerance used to decide whether an energy is lower than another during
    a GA-based structure prediction/solution calculation. For example, if the
    "best" structure has an energy of 4.00001 eV/atom, it is rather
    insignificant if a second structure has an energy of exactly 4.0 eV/atom -
    for all practical purposes, the energies are the same. This speeds up
    convergence of simulations by preventing "best" structure updates when
    energy differences are negligibly small.

Values: Any floating point number (in units of eV/atom)

Default: 0.001



### gaoptdifftol

General: Tolerance used to decide whether an R-factor is lower than another
    during a GA-based structure prediction/solution calculation. For example,
    if the "best" structure has an R-factor for 0.020001, it is rather
    insignificant if a second structure has an R-factor of exactly 0.02 - for
    all practical purposes, the R-factors are the same. This speeds up
    convergence of simulations by preventing "best" structure updates when
    R-factor differences are negligibly small.

Values: Any floating point number

Default: 1e-4



### gaoptscreenmethod

General: Set the screening method used when generating children during a
    GA-based structure prediction/solution calculation. For example, if the
    screening method is set to "energy" and gaoptscreennum to 10, then for
    each mating operation, ten candidates are generated for each child and the
    one with the lowest energy (without relaxation) is taken. Similarly, if
    the screening method is set to "xray", then of the ten candidates, the
    one with the lowest R-factor (again without relaxing positions) is
    accepted. Screens are used to improve simulation times by selecting
    children nearest to a minimum so that local relaxations occur more rapidly.

    This is distinct from the number of children in a simulation. For example,
    in a simulation with 10 structures in the population and gaoptscreennum
    set to 5, for each of the 10 children in a new generation, 5 candidates are
    produced and the one that best optimizes the metric defined by
    gaoptscreenmethod is accepted as the child. The selected children are
    relaxed locally using whichever metric is being optimized in the simulation.

    In order for a screen to be applied, gaoptscreenmethod must be set and
    gaoptscreennum must have a non-zero value. If one or neither of these are
    true, then a screen is not used.

Values: 
    "energy" or "potential" to screen using the energy
    "xray" or "diffraction" to screen using the r-factor

Default: None (no screen is applied)



### gaoptscreennum

General: Set the number of candidates to generate and screen for each child
    in each generation of a GA-based structure prediction/solution calculation.
    See gaoptscreenmethod for more details.

Values: Any integer number

Default: 0 (no screen is applied)



### wyckoffbias

General: From analysis of structures in the Inorganic Crystal Structure Database
    it was found that atoms tend to occupy Wyckoff site combinations so as to
    minimize the number of symmetrically-unique atoms in the structure. For
    example, assume that a space group has three Wyckoff sites, two with
    multiplicity 1 and one with multiplicity 2. If two atoms were placed in the
    structure, then based on this analysis, they would overwhelmingly prefer
    to occupy the site with multiplicity 2, rather than some combination of
    two sites with multiplicity 1.

    The exact statistics from this analysis are hard-coded into Mint and used
    to determine which Wyckoff sites are selected when randomly placing atoms
    into a structure with a set space group. This setting controls the extent
    to which these statistics are used. If set to 0, then no biasing is applied
    and Wyckoff site combinations are chosen completely at random. If set to 1,
    then biasing is maximized. Typically, biasing should be somewhere between
    these two extremes. If it is too low, then the simulation may take much
    longer to converge. It it is too high, then it can be difficult or even
    impossible to find structures that differ significantly from those in the
    ICSD, at least with respect to the Wyckoff site occupancies.

Values: Any floating point number between 0 (no biasing) and 1 (full biasing)

Default: 0.5



### minimagedis

General: Several functions within Mint generate supercells of structures. This
    setting controls the size of the supercell by controlling the minimum image
    distance that is allowed in the structure under periodic boundary
    conditions. A larger minimum image distance will always lead to a cell
    that is at least as large to, if not larger than, a cell with a smaller
    minimum image distance.

Values: Any floating point number (in units of Angstroms)

Default: 8.0



### maxjumpdistance

General: Maximum jumps distance allowed when generating jumps between atomic
    sites in a structure.

Values: Any floating point number (in units of Angstroms)

Default: 7.0



### kmcjumpsperatom

General: Number of jumps to perform in a KMC simulation. This setting in given
    as a function of the number of atoms in a structure so the actual number
    of jumps will depend on the simulation cell size. A larger number will
    lead to smaller errors in calculated diffusivities.

Values: Any integer number (in units of jumps/atom)

Default: 100



### kmcconverge

General: Convergence criterion for a KMC simulation to be considered complete.
    The simulation is terminated once the standard error is below this value.

Values: Any floating point number (in units of percent)

Default: 0.5
