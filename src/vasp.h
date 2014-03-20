/* vasp.h -- VASP functionality
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef VASP_H
#define VASP_H



#include "iso.h"
#include "elements.h"
#include "text.h"
#include "list.h"
#include "fileSystem.h"



// Namespace for VASP functionality
namespace Vasp
{
	
	// Classes
	class Job;
	class Structure;
	class Settings;
	class Potential;
	
	// Enumerated types
	enum Calculation {VC_NONE, VC_SINGLE, VC_RELAX_ALL, VC_RELAX_POSITIONS, VC_MD, VC_NEB, VC_CINEB};
	enum Relaxation {VR_UNKNOWN, VR_QUASI_NEWTON, VR_CONJUGATE_GRADIENT, VR_DAMPED_MD};
	enum Accuracy {VA_UNKNOWN, VA_HIGH, VA_NORMAL, VA_LOW};
	enum Material {VM_UNKNOWN, VM_METAL, VM_INSULATOR};
	enum IncarGroup {IG_UNKNOWN, IG_GENERAL, IG_ELECTRONIC, IG_ELECTRONIC_SMEARING, IG_ELECTRONIC_PROJECTORS, \
		IG_ELECTRONIC_STARTUP, IG_ELECTRONIC_SPIN, IG_ELECTRONIC_EXCHANGE_CORRELATION, IG_ELECTRONIC_CONVERGENCE, \
		IG_ELECTRONIC_CONVERGENCE_DETAIL, IG_ELECTRONIC_MIXER, IG_ELECTRONIC_MIXER_DETAILS, \
		IG_ELECTRONIC_DIPOLE_CORRECTION, IG_GRIDS, IG_IONIC, IG_IONIC_MD, IG_SYMMETRY, IG_DOS, IG_WRITING, \
		IG_PERFORMANCE, IG_MISCELLANEOUS, IG_OTHER, IG_VARIATIONAL_TRANSITION_STATE};
	enum IncarTag {IT_UNKNOWN, IT_SYSTEM, IT_LCOMPAT, IT_PREC, IT_ENMAX, IT_ENAUG, IT_EDIFF, IT_ALGO, IT_IALGO, \
		IT_IWAVPR, IT_NBANDS, IT_NELECT, IT_ISMEAR, IT_SIGMA, IT_LREA, IT_ROPT, IT_LMAXPAW, IT_LMAXMIX, IT_ISTART, \
		IT_ICHARG, IT_INIWAV, IT_ISPIN, IT_LNONCOLLINEAR, IT_MAGMOM, IT_NUPDOWN, IT_LSORBIT, IT_SAXIS, IT_LSPIRAL, \
		IT_QSPIRAL, IT_LZEROZ, IT_GGA, IT_VOSKOWN, IT_LASPH, IT_LMETAGGA, IT_GGA2, IT_NELM, IT_NELMDL, IT_NELMIN, \
		IT_ENINI, IT_LDIAG, IT_WEIMIN, IT_EBREAK, IT_DEPER, IT_NRMM, IT_TIME, IT_AMIX, IT_BMIX, IT_AMIN, IT_AMIX_MAG, \
		IT_BMIX_MAG, IT_IMIX, IT_MAXMIX, IT_WC, IT_INIMIX, IT_MIXPRE, IT_MREMOVE, IT_IDIPOL, IT_LDIPOL, IT_DIPOL, \
		IT_EFIELD, IT_NGX, IT_NGY, IT_NGZ, IT_NGXF, IT_NGYF, IT_NGZF, IT_ADDGRID, IT_NSW, IT_IBRION, IT_ISIF, \
		IT_PSTRESS, IT_EDIFFG, IT_NFREE, IT_POTIM, IT_SMASS, IT_TEBEG, IT_TEEND, IT_NBLOCK, IT_KBLOCK, IT_NPACO, \
		IT_APACO, IT_ISYM, IT_SYMPREC, IT_LORBIT, IT_RWIGS, IT_NEDOS, IT_EMIN, IT_EMAX, IT_NWRITE, IT_LWAVE, \
		IT_LCHARG, IT_LPARD, IT_LVTOT, IT_LELF, IT_LOPTICS, IT_STM, IT_NPAR, IT_NSIM, IT_NBLK, IT_LPLANE, \
		IT_LCRITICAL_MEM, IT_LSCALAPACK, IT_LSCALU, IT_LASYNC, IT_IDIOT, IT_LMUSIC, IT_POMASS, IT_LCORR, IT_LREAL, \
		IT_LREAL_COMPAT, IT_GGA_COMPAT, IT_LBERRY, IT_ICORELEVEL, IT_LDAU, IT_I_CONSTRAINED_M, IT_ICHAIN, IT_IMAGES, \
		IT_SPRING, IT_LCLIMB, IT_LTANGENTOLD, IT_LDNEB, IT_NEBCELL, IT_IOPT};
	
	// Functions
	Accuracy accuracy(const Word& word);
	Word accuracy(Accuracy accuracy);
}



// Vasp job
class Vasp::Job
{
	
	// Control variables
	bool _saveFiles;
	Calculation _calculation;
	Relaxation _relaxation;
	Accuracy _accuracy;
	Material _material;
	Word _executable;
	OList<Element> _elements;
	OList<Word> _potentials;
	
	// Variables to store results
	List<double> _energy;
	List<double>::D2 _convergence;
	OList<Vector3D >::D2 _forces;
	ISO _structure;
	OList<ISO> _structures;
	OList<Word>::D2 _tags;
	
	// Static member variables
	static bool _completed;
	static Word _dirName;
	static Word _prevDir;
	
	// Functions
	void writeFiles(const ISO& iso, bool restart);
	void writeNEBFiles(const OList<ISO>& iso, bool restart);
	void gatherResults();
	void gatherNEBResults(int numDirs);

public:
	
	// Constructor and destructor
	Job();
	
	// Setup functions
	void saveFiles(bool saveFiles)					{ _saveFiles = saveFiles; }
	void calculation(Calculation calculation)		{ _calculation = calculation; }
	void relaxation(Relaxation relaxation)			{ _relaxation = relaxation; }
	void accuracy(Accuracy accuracy)				{ _accuracy = accuracy; }
	void material(Material material)				{ _material = material; }
	void executable(const Word& executable)			{ _executable = executable; }
	void elements(const OList<Element>& elements)	{ _elements = elements; }
	void potentials(const OList<Word>& potentials)	{ _potentials = potentials; }
	void setTags(const OList<Word>::D2& tags)		{ _tags = tags; }
	
	// Functions
	bool submit(const ISO& iso, bool restart, bool quitIfError = true);
	bool submitNEB(const OList<ISO>& iso, bool restart, bool quitIfError = true);
	
	// Access functions
	bool completed() const							{ return _completed; }
	const List<double>& energy() const				{ return _energy; }
	const List<double>::D2& convergence() const		{ return _convergence; }
	const OList<Vector3D >::D2& forces() const		{ return _forces; }
	void updateStructure(ISO& iso) const;
	void updateStructures(OList<ISO>& isos) const;
};



// Vasp structure functions
class Vasp::Structure
{
public:

	// Read structure
	static ISO read(const Text& content, bool useVasp5);
	static ISO read(const Word& file, bool useVasp5)	{ return read(Read::text(file), useVasp5); }
	
	// Write structure
    static void write(const Word& file, const ISO& iso, CoordinateType coordinates, bool useVasp5);

	// Check if file is in correct format
	static bool isVersion4(const Word& file)		{ return isVersion4(Read::text(file)); }
    static bool isVersion4(const Text& content);
	static bool isVersion5(const Word& file)		{ return isVersion5(Read::text(file)); }
	static bool isVersion5(const Text& content);
};



// Vasp settings functions
class Vasp::Settings
{
	
	// Single setting
	struct Setting
	{
		bool forced;
		IncarGroup group;
		IncarTag tag;
		Word name;
		Word value;
	};
	
	// Variables
	OList<Setting> _settings;
	
	// Functions
	static bool match(const Word& name, IncarGroup* group, IncarTag* tag, const char* comp, IncarGroup compGroup, \
		IncarTag compTag);
	void write(IncarGroup group, const char* header) const;
	
	// Default functions
	void generalDefaults(bool replaceExisting);
	void calculationDefaults(Calculation calculation, bool replaceExisting);
	void materialDefaults(Material material, bool replaceExisting);
	void accuracyDefaults(Accuracy accuracy, bool replaceExisting);
	
public:
	
	// Constructor
	Settings()	{ defaults(VC_SINGLE, VR_QUASI_NEWTON, VM_METAL, VA_NORMAL); }
	
	// Setup functions
	void clear();
	void defaults(Calculation calculation, Relaxation relaxation, Material material, Accuracy accuracy, \
		const List<double>* encut = 0);
	bool add(const Word& name, const Word& value, bool replaceExisting, bool forced = false);
	void remove(const Word& name);
	
	// Read and write functions
	void read(const Text& content);
	void read(const Word& file)			{ read(Read::text(file)); }
	bool isFormat(const Text& content);
	bool isFormat(const Word& file)		{ return isFormat(Read::text(file)); }
	void write(const Word& file) const;
	static bool isSetting(const Word& name, IncarGroup* group = 0, IncarTag* tag = 0);
};



// Vasp potential
class Vasp::Potential
{
public:
	static void write(const Word& file, const OList<Element>& elements, const OList<Word>& files, const ISO& iso, \
		List<double>* encut = 0);
};



// =====================================================================================================================
// Job
// =====================================================================================================================

/* inline Vasp::Job::Job()
 *
 * Constructor for vasp job
 */

inline Vasp::Job::Job()
{
	_saveFiles = false;
	_calculation = VC_SINGLE;
	_relaxation = VR_CONJUGATE_GRADIENT;
	_accuracy = VA_NORMAL;
	_material = VM_METAL;
}



// =====================================================================================================================
// Settings
// =====================================================================================================================

/* inline void Vasp::Settings::clear()
 *
 * Clear settings
 */

inline void Vasp::Settings::clear()
{
	_settings.clear();
}



#endif
