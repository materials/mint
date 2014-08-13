/* Copyright 2011-2014 Kyle Michel, Logan Ward, Christopher Wolverton
 *
 * Contact: Kyle Michel (kylemichel@gmail.com)
 *			Logan Ward (LoganWard2012@u.northwestern.edu)
 *
 *
 * This file is part of Mint.
 *
 * Mint is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * Mint is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with Mint.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef DIFFRACTION_H
#define DIFFRACTION_H

#include "num.h"
#include "iso.h"
#include "elements.h"
#include "symmetry.h"
#include "fileSystem.h"
#include "text.h"
#include "list.h"
#include "constants.h"
#include <cmath>
#include <vector>
#include <set>
#include <valarray>
#include "dlib/optimization.h"
#include "output.h"

// Set == 1 to print out diffraction pattern at each stage
#define LW_EXCESSIVE_PRINTING 1

// Diffraction pattern determination methods
enum Method {DM_XRAY, DM_NEUTRON, DM_SIMPLE, DM_NONE};

/**
 * Utility class designed to store peak positions and intensities. Also contains 
 *  the ability to calculate integrated peak intensities for patterns calculated from
 *  a known structure.
 */
class DiffractionPeak {
public:
	
    /**
     * Create a peak found in a measured diffraction pattern (only requires the 
     *  measured intensity and angle)
     * @param twoThetaDegrees Bragg angle of peak in degrees
     * @param intensity Intensity of peak (units don't matter)
     */
    DiffractionPeak(double twoThetaDegrees, double intensity) {
        this->TwoThetaDeg = twoThetaDegrees;
        this->_twoThetaRad = Num<double>::toRadians(twoThetaDegrees);
        this->peakIntensity = intensity;
        this->patternIndex = -1;
    }
    
    // --> Operations to access data about this peak
    /** 
     * Get Bragg angle (twoTheta) for this peak in degrees.
     */
    double getAngle() const { return TwoThetaDeg; }
    
    // Get integrated intensity for this peak
    double getIntensity() const { return this->peakIntensity; }
    
    // Index of matching peak in reference pattern
    int patternIndex;
	
    // ---> Utility operations (allow this class to be used efficiently)
    DiffractionPeak& operator= (const DiffractionPeak& rhs) {
            if (this != &rhs) {
                    patternIndex = rhs.patternIndex;
                    TwoThetaDeg = rhs.TwoThetaDeg;
                    _twoThetaRad = rhs._twoThetaRad;
                    peakIntensity = rhs.peakIntensity;
            }
            return *this;
    }
    
    bool operator< (const DiffractionPeak rhs) const {
        return TwoThetaDeg < rhs.TwoThetaDeg;
    }
    
protected:
        
    // Angle of peak in radians
    double _twoThetaRad;
    // Angle of peak in degrees
    double TwoThetaDeg;
    // Intensity of peak 
    double peakIntensity;
    
};

/**
 * Create a new peak corresponding to a peak in a calculated diffraction pattern. 
 * For computational efficiency, this constructor assumes many of the parameters
 * related to symmetry have already been calculated. 
 * 
 * Developer's note: Since the Bragg angle and reciprocal lattice vectors are supplied,
 *  peaks diffraction angles will not change if the lattice parameters
 *  of the structure used to generate this object are changed.
 * 
 * @param wavelength Wavelength of incident radiation
 * @param twoThetaDegrees Bragg angle of peak in degrees
 * @param multiplicity Multiplicity factor for this peak
 * @param hkl A plane index that represents this peak
 * @param equivHKL Reciprocal lattice vectors of all planes
 *   contributing to this peak
 */
class CalculatedPeak : public DiffractionPeak {
	
public:
	/**
     * Create a new peak corresponding to a peak in a calculated diffraction pattern. 
     * For computational efficiency, this constructor assumes many of the parameters
     * related to symmetry have already been calculated. 
     * 
     * Developer's note: Since the Bragg angle and reciprocal lattice vectors are supplied,
     *  peaks diffraction angles will not change if the lattice parameters
     *  of the structure used to generate this object are changed.
	 * 
	 * @param method Method used to compute diffraction peaks
	 * @param structure ISO representing the structure from which this peak is calculated
	 * @param symmetry Symmetry object representing symmetry of structure
     * @param wavelength Wavelength of incident radiation
     * @param multiplicity Multiplicity factor for this peak
     * @param hkl A plane index that represents this peak
     * @param reciprocalLatticeVectors Reciprocal lattice vectors of all planes
     *   contributing to this peak
     */
    CalculatedPeak(Method method, const ISO* structure, const Symmetry* symmetry, double wavelength, 
			Vector3D hkl, std::vector<Vector3D> equivHKL) :	DiffractionPeak(-1, -1) {
		this->method = method;
		this->structure = structure;
		this->symmetry = symmetry;
		this->_twoThetaRad = CalculatedPeak::getDiffractionAngle(structure->basis(),hkl,wavelength);
		this->TwoThetaDeg = Num<double>::toDegrees(this->_twoThetaRad);
        this->wavelength = wavelength;
		this->multiplicity = equivHKL.size();
        this->hkl = hkl;
        this->equivHKL = equivHKL;
        this->lpFactor = getLPFactor(_twoThetaRad / 2.0);
    }

    void updateCalculatedIntensity(vector<double> BFactors, List<double>::D2 atfParams);
	
	static double getDiffractionAngle(const Basis& basis, const Vector3D& hkl, double wavelength);	
	
	CalculatedPeak& operator= (const CalculatedPeak& rhs) {
            if (this != &rhs) {
                    patternIndex = rhs.patternIndex;
                    TwoThetaDeg = rhs.TwoThetaDeg;
                    _twoThetaRad = rhs._twoThetaRad;
                    peakIntensity = rhs.peakIntensity;
                    lpFactor = rhs.lpFactor;
                    multiplicity = rhs.multiplicity;
                    hkl = rhs.hkl;
                    equivHKL = rhs.equivHKL;
            }
            return *this;
    }
	
private:
	
	// Method used to compute diffraction intensity
	Method method;
	// Pointer to the structure 
	const ISO* structure;
	// Pointer the symmetry description of structure
	const Symmetry* symmetry;
    // Wavelength of incident radiation
    double wavelength;
    // Index of family of planes which result in this peak
    Vector3D hkl;
    // Reciprocal lattice vectors of each member of the family of planes
    std::vector<Vector3D> equivHKL;
    // Lorentz polarization factor for this peak
    double lpFactor;
    // Multiplicity for this peak
    double multiplicity;
	
	// ---> Operations used to calculate diffraction peak intensity
    static double getLPFactor(double angle);
    static double thermalFactor(double angle, double wavelength, double Bfactor = 1);
    static double getAbsorptionFactor(double angle, double uEff);
    static double getTexturingFactor(Vector3D preferredOrientation, double tau, Linked<Vector3D> recipLatticeVectors );
    static double structureFactorSquared(Method method, double wavelength, const Symmetry& symmetry, double angle, \
			const Vector3D& hkl, vector<double> BFactors, List<double>::D2 atfParams);
    static double atomicScatteringFactor(List<double>& atfParams, double angle, double wavelength);
};

/**
 * Object used to calculate, store, and compare diffraction peaks. 
 */
class Diffraction {
	
public:
    
    // Pattern type
    enum PatternType {PT_EXP_RAW, // Raw pattern from experiment, raw pattern
        PT_EXP_INT, // Raw pattern from experiment, integrated intensities
        PT_CALCULATED}; // Pattern calculated from a crystal structure
		
	// Methods for calculating R factors
    enum Rmethod {DR_ABS, /** Method often used to report matches in lit : <br>
                           * R = sum(|I_ref - scale * I_calc|) / sum(I_ref) */
    DR_SQUARED /** Method used in refinement with integrated intensities 
                * (because it is differentiable):<br>
                * R = sum[ (I_ref - I_calc) ^ 2 ] / sum[ I_ref^2 ] */
    };
        
protected:
      	
    // ==============================
    // General variables and methods
    // ==============================
    
    // Type of diffraction pattern
    PatternType _type;
	
	// Experimental technique used to generate a pattern
    Method _method;
        
    // Wavelength of diffracted radiation
    double _wavelength;
    // Minimum angle at which diffracted intensities were measured
    double _minTwoTheta;
    // Maximum angle at which diffracted intensities were measured
    double _maxTwoTheta;
    // Minimum distance between any two features in a diffraction pattern (degrees)
    static double _resolution;
	
	// Index of peaks in this pattern which match a given peak in the reference
    vector<vector<int> > _matchingPeaks;
    // Index of peaks that do not match something in the calculated pattern
    vector<int> _unmatchedPeaks;
	
	// Used when calculating match between two problems
	
	virtual void matchPeaksToReference(const Diffraction& referencePattern);
	
	double getCurrentRFactor(const Diffraction& referencePattern, Rmethod rMethod = DR_ABS);
                
public:
	
	/**
	 * Given diffraction angle, return diffracted intensity
     * @param twoTheta [in/out] Diffraction angle will be sorted
     * @return Diffracted intensity at each specified angle
     */
	virtual vector<double> getDiffractedIntensity(vector<double>& twoTheta) const = 0;
	
	/**
	 * Get a list of diffraction peaks
	 * @return List of diffraction peaks
     */
	virtual vector<DiffractionPeak> getDiffractedPeaks() const = 0;
	
	// Constructors
	Diffraction();
	
	// Set settings
	void setMethod(Method input)		{ _method = input; }
	void setWavelength(double input)	{ _wavelength = input; }
	void setMinTwoTheta(double input)	{ _minTwoTheta = input; }
	void setMaxTwoTheta(double input)	{ _maxTwoTheta = input; }
	
	// Setup functions
	virtual void clear() {
		_matchingPeaks.clear();
		_unmatchedPeaks.clear();
	}
	
	// Get R factor for two patterns
	double rFactor(const Diffraction& reference);
	
	// Print functions
	void print(const Word& file, bool broaden = false) const;
	
	// Access functions
	PatternType patternType() { return _type; }
	bool isSet() const { return (_method != DM_NONE); }
	double wavelength() const { return _wavelength; }
	Method method() const { return _method; }
	double minTwoTheta() const { return _minTwoTheta; }
	double maxTwoTheta() const { return _maxTwoTheta; }         
};

/**
 * Represents a powder diffraction pattern from an experimental source. 
 * 
 * <p>Key functionality:
 * <ul>
 * <li>Read in diffraction patterns from file
 * </ul>
 */
class ExperimentalPattern : public Diffraction {
private: 
	// =====================================
	// Data describing the diffraction pattern
	// =====================================
	// Continuous x-ray pattern: Only filled when a raw pattern is imported
    valarray<double> _continuousTwoTheta; // LW 12Aug14: Rename this after porting
    // Continuous x-ray pattern: Only filled when a raw pattern is imported
    valarray<double> _continuousIntensity;
	
	/** 
	 * Stores each diffraction peak in the pattern. These values are either 
	 * stored directly from input data if integrated intensities are provided,
	 * or calculated directly from a raw pattern.
     */
    vector<DiffractionPeak> _diffractionPeaks;
	
	// =========================================
	// Stuff related to processing raw patterns
	// =========================================

	// Analyze a raw pattern
	void set(vector<double>& twoTheta, vector<double>& intensity);
	void smoothData(const vector<double>& rawTwoTheta, vector<double>& rawIntensity, \
			const int numPerSide = 2, const double power = 0.25);
	void removeBackground(vector<double>& rawTwoTheta, vector<double>& rawIntensity);
	void locatePeaks(vector<vector<double> >& peakTwoTheta, vector<vector<double> >& peakIntensity, \
			const vector<double>& rawTwoTheta, const vector<double>& rawIntensity);
	void getPeakIntensities(const vector<vector<double> >& peakTwoTheta, \
			const vector<vector<double> >& peakIntensity);
	vector<double> getFirstDerivative(const vector<double>& x, const vector<double>& y);
	vector<double> getSecondDerivative(const vector<double>& x, const vector<double>& y);

	// Functions for fitting Gaussian function
	// Params order H, 2*theta_k, I0
	double gaussian(const Vector& params, double twoTheta);
	Vector gaussianDerivs(const Vector& params, double twoTheta);
	double compositeGaussian(const Vector& params, double twoTheta);
	Vector compositeGaussianDerivs(const Vector& params, double twoTheta);

	// Functions for fitting pseudo-Voigt function
	// Params order: eta0, eta1, eta2, 2*theta_k, u, v, w, I0
	double PV(const Vector& params, double twoTheta);
	Vector PVderivs(const Vector& params, double twoTheta);
	double PVderiv(const Vector& params, double twoTheta);
	double compositePV(const Vector& params, double twoTheta);
	Vector compositePVDerivs(const Vector& params, double twoTheta);

	// Debugging functions
	void savePattern(const Word& filename, const vector<double>& twoTheta,
			const vector<double>& Intensity, const vector<double>& otherIntensity = vector<double>(0) );
	
public:

	void clear() {
		_continuousIntensity.resize(0);
		_continuousTwoTheta.resize(0);
		_diffractionPeaks.clear();
	}
	

	virtual vector<DiffractionPeak> getDiffractedPeaks() const;

	vector<double> getDiffractedIntensity(vector<double>& twoTheta) const;

	// Set from file or other data
	void set(const Text& text);
	void set(const Word& file)	{ set(Read::text(file)); }
	void set(const Linked<double>& twoTheta, const Linked<double>& intensity);
	
	// Check if file is correct format
	static bool isFormat(const Text& text);
	static bool isFormat(const Word& file)	{ return isFormat(Read::text(file)); }
	
	// Friendships
	friend class PVPeakFunction;
};

/**
 * Represents the calculated powder diffraction pattern from a structure.
 * 
 * Key functionality:
 * <ul>
 * <li>Generate peak locations and intensities
 * <li>Generate intensity as a continuous function of angle
 * <li>Refine crystal structure and other parameters to better match a reference pattern
 * </ul>
 */
class CalculatedPattern : public Diffraction {
public:
	// Used by Dlib
    typedef dlib::matrix<double, 0, 1> column_vector;
	
private:
	
	/**
	 * Holds each reflection from the structure.
	 */
	vector<CalculatedPeak> _reflections;
	
	// =========================================
	// Parameters used when generating a pattern
	// =========================================
	
	// Holds the internal degrees of freedom (atomic positions)
    const Symmetry* _symmetry;
    // Holds information about basis vectors and such
    const ISO* _structure;
    // B factors of each symmetrically-unique set of atoms
    vector<double> _BFactors;
	// Atomic form factor parameters for each symmetrically-unique set of atoms
    List<double>::D2 _atfParams;
	// Minimum allowed B factor
    double _minBFactor;
    // Maximum allowed B factor
    double _maxBFactor;
	
	// ==============================================================
	// Parameters used when refining structural model against pattern
	// ==============================================================
	
	// Parameters that can be refined
	enum RefinementParameters { RF_BFACTORS, RF_POSITIONS };
	// Parameters that are currently being refined
    std::set<RefinementParameters> _currentlyRefining;
	
	// =================================================================
	// Operations that ready this object to calculate peaks / be refined
	// =================================================================
	
	// --> Operations used to define refinement problem
	void defineReferencePattern(const Diffraction& reference);
    bool referencePatternIsDefined() const;
    void defineStructure(const ISO& structure, const Symmetry& symmetry);
    bool structureIsDefined();
    void initializeRefinementParameters();
    void setATFParams();
    void calculatePeakLocations();
	void matchPeaksToReference(Diffraction* referencePattern);
	
	// --> Operation employed by user to refine a calculated pattern
    void refineParameters(const Diffraction* referencePattern, std::set<RefinementParameters> toRefine);
	
	// --> Operations are used when calculating peak intensities
    void calculatePeakIntensities();
    void scaleIntensity(double max);
	
	// --> Operations used directly by refineStructure
    bool willRefine(RefinementParameters parameter, std::set<RefinementParameters> toRefine);
    double runRefinement(const Diffraction* ref);
    column_vector getRefinementParameters();
    column_vector getRefinementParameterDerivatives(column_vector x);
    column_vector getRefinementParameterLowerBoundary();
    column_vector getRefinementParameterUpperBoundary();
    void setAccordingToParameters(column_vector params);
	
	// --> Utility operations used during refinement
    void setPositions(const Symmetry& symmetry, const Vector& positions);
    void symPositions(const Symmetry& Symmetry, Vector& positions);
	
public:
	
	vector<DiffractionPeak> getDiffractedPeaks() const;
	
	vector<double> getDiffractedIntensity(vector<double>& twoTheta) const;
	
	virtual void clear() {
		_structure = 0;
	}
	
	CalculatedPattern() {
		_symmetry = 0;
        _structure = 0;
		_minBFactor = 0.1;
		_maxBFactor = 4.0;
	}
	
	// Friendships
	friend class RFactorFunctionModel;
	
	// LW 12Aug14: Eventually, make this the constructor
	double set(const ISO& iso, const Symmetry& symmetry, const Diffraction* ref = 0, bool fitBfactors = false);
	
	// Call refinement
	double refine(ISO& iso, Symmetry& symmetry, const Diffraction& reference, bool showWarnings = true);
};

/**
 * This class provides an interface to the Dlib optimization libraries for the purposes
 *  of evaluating the R factor as a function of some parameters.
 * 
 * Before calling Dlib, create an instance of this model by calling the constructor
 *  with a reference to "this" (the class describing the calculated diffraction pattern).
 * 
 * @param toRefine Diffraction pattern that will be refined.
 * @param reference Reference diffraction pattern
 */
class RFactorFunctionModel {
public:
    RFactorFunctionModel (CalculatedPattern* toRefine, const Diffraction* reference) {
        _toRefine = toRefine;
		_referencePattern = reference;
    }
    
    double operator() (const CalculatedPattern::column_vector& arg) const {
        // Change parameters of diffraction pattern
        _toRefine->setAccordingToParameters(arg);
        // Recalculate peak intensities
        _toRefine->calculatePeakIntensities();
        // Get the current R factor
        return _toRefine->getCurrentRFactor(*_referencePattern, Diffraction::DR_SQUARED);
    }
    
private:
    CalculatedPattern* _toRefine;
	const Diffraction* _referencePattern;
};

/**
 * This class is created for the sole purpose of interfacing with dlib's numerical
 *  integration routines. 
 * @param pattern
 * @param params
 */
class PVPeakFunction {
public:
    PVPeakFunction(ExperimentalPattern* pattern, Vector& params) {
        _pattern = pattern;
        _params = params;
    }
    
    double operator() (double x) const {
        return _pattern->PV(_params, x);
    }
    
private:
    ExperimentalPattern* _pattern;
    Vector _params;
    
};

// =====================================================================================================================
// Diffraction
// =====================================================================================================================

/**
 * Constructor for Diffraction object
 */
inline Diffraction::Diffraction()
{
	_method = DM_NONE;
	_wavelength = 1.5418;
	_minTwoTheta = 10;
	_maxTwoTheta = 100;
}

/**
 * Calculate Lorentz polarization factor
 * @param angle [in] Angle of peak in diffraction pattern
 * @return Lorentz polarization factor for that peak
 */
inline double CalculatedPeak::getLPFactor(double angle)
{
	return (1 + pow(cos(2 * angle), 2)) / (cos(angle) * pow(sin(angle), 2));
}

/**
 * Calculate absorption factor for a certain peak. Uses the following relationship:
 * <center>A = 1 - exp(-2*u<sub>eff</sub>/sin(theta))</center>
 * Note that thickness of the sample has been absorbed into u<sub>eff</sub>.
 * @param angle Angle of diffraction peak.
 * @param uEff Effective absorption coefficient.
 * @return Absorption factor for this peak
 */
inline double CalculatedPeak::getAbsorptionFactor(double angle, double uEff) {
    return 1 - exp(-2 * uEff / sin(angle));
}

/**
 * Calculate the texturing factor for a peak. 
 * 
 * Currently uses the March-Dollase function: <br> 
 * <code>T_hkl = 1 / N * sum_i{ ( tau^2 * cos(phi^i_hkl)^2 + 1 / tau * sin(phi^i_hkl)^2 )^-(3/2) }</code>
 * @param preferredOrientation Reciprocal lattice vector of preferred orientation
 * @param tau Texturing parameter (generally, 1 means no texturing effect)
 * @param recipLatticeVectors Reciprocal lattice vectors of all planes contributing
 *  to this diffraction peak
 * @return Texturing factor
 */
inline double CalculatedPeak::getTexturingFactor(Vector3D preferredOrientation, double tau, Linked<Vector3D> recipLatticeVectors) {
    double output = 0.0;
    double preNorm = preferredOrientation.magnitude();
    for (Linked<Vector3D>::iterator iter = recipLatticeVectors.begin(); 
            iter != recipLatticeVectors.end(); iter++) {
        double cosphi = preferredOrientation * *iter / preNorm / (*iter).magnitude();
        cosphi *= cosphi;
        output += pow(1 + (tau * tau - 1) * cosphi, -0.5);
        // output += pow(tau * tau * cosphi + (1 - cosphi) / tau, -1.5);
    }
    output /= (double) recipLatticeVectors.length();
    return output;
}

/**
 * Calculate a thermal factor
 * @param angle [in] Angle at which a certain peak diffractions
 * @param wavelength [in] Wavelength of incident radiation
 * @param Bfactor [in] Thermal factor parameter for a certain atom
 * @return Thermal factor
 */
inline double CalculatedPeak::thermalFactor(double angle, double wavelength, double Bfactor)
{
	return exp(-Bfactor * pow(sin(angle) / wavelength, 2.0));
}

/**
 * Given a set of basis vector and a certain plane index, calculate the angle at which
 *  the corresponding diffraction peak will appear. Wavelength value is stored 
 *  internally.
 * @param basis [in] Basis vectors of a particular unit cell
 * @param hkl [in] Index of a particular plane
 * @param wavelength [in] Wavelength of incident radiation
 * @return Angle at which that plane will appear in radians
 */
inline double CalculatedPeak::getDiffractionAngle(const Basis& basis, const Vector3D& hkl, double wavelength)
{
	double arg = (basis.inverse() * hkl).magnitude() * wavelength / 2;
	if ((arg >= -1) && (arg <= 1))
		return asin(arg);
	else if (arg < -1)
		return -Constants::pi / 2;
	return Constants::pi / 2;
}


/**
 * Evaluate a Gaussian function
 * @param params [in] Parameters of Gaussian (Cg, mu, sigma)
 * @param twoTheta [in] Value at which function is evaluated
 * @return Value of Gaussian at that angle
 */
inline double ExperimentalPattern::gaussian(const Vector& params, double twoTheta) {
	// Set variables
	double pi = Constants::pi;
	double Cg  = 4*log(2);
	double dif = twoTheta - params[1];
	double ex  = exp(-Cg * dif*dif / params[0]);
	
	// Return result
	return params[2] * sqrt(Cg) * ex / sqrt(pi * params[0]);
}



/**
 * Calculate Derivatives of a Gaussian function 
 * 
 * @param params [in] Parameters of Gaussian (Cg, mu, sigma)
 * @param twoTheta [in] Value at which derivative is evaluated
 * @return Derivatives wrt each parameter
 */
inline Vector ExperimentalPattern::gaussianDerivs(const Vector& params, double twoTheta) {
	
	// Set variables
	double pi = Constants::pi;
	double Cg  = 4*log(2);
	double dif = twoTheta - params[1];
	double ex  = exp(-Cg * dif*dif / params[0]);

        // Calculate derivatives
	Vector res(3);
	res[0] = params[2] * sqrt(Cg) * (Cg*dif*dif - params[0]) * ex / (2 * sqrt(pi) * pow(params[0], 2.5));
	res[1] = 2 * pow(Cg, 1.5) * params[2] * dif * ex / (sqrt(pi * params[0])*params[0]);
	res[2] = sqrt(Cg) * ex / sqrt(pi * params[0]);
	
	// Return derivatives
	return res;
}

/**
 * Calculate the composite of multiple Gaussian functions
 * <p>The first 3 parameters are for the first Gaussian, the second 3 are for the 
 * second, and so on.
 * @param params [in] Parameters of Gaussian functions (Cg_1, mu_1, sigma_1, Cg_2, ...)
 * @param twoTheta [in] Angle at which function is evaluated 
 * @return Value of multiple Gaussians at that angle
 */
inline double ExperimentalPattern::compositeGaussian(const Vector& params, double twoTheta) {
    double output = 0;
    for (int f = 0; f < params.length()/3; f++) {
        Vector subParams(3);
        for (int i=0; i<3; i++) subParams[i] = params[f * 3 + i];
        output += gaussian(subParams, twoTheta);
    }
    return output;
}

/**
 * Calculate derivatives of the composite of multiple Gaussian functions
 * <p>The first 3 parameters are for the first Gaussian, the second 3 are for the 
 * second, and so on.
 * @param params [in] Parameters of Gaussian functions (Cg_1, mu_1, sigma_1, Cg_2, ...)
 * @param twoTheta [in] Angle at which derivatives are evaluated 
 * @return Derivatives with respect to each parameter
 */ 
inline Vector ExperimentalPattern::compositeGaussianDerivs(const Vector& params, double twoTheta) {
    Vector output(params.length(), 0);
    for (int f = 0; f < params.length()/3; f++) {
        Vector subParams(3);
        for (int i=0; i<3; i++) subParams[i] = params[f * 3 + i];
        Vector subDerivs = gaussianDerivs(subParams, twoTheta);
        for (int i=0; i<3; i++) output[f * 3 + i] = subDerivs[i];
    }
    return output;
}

/**
 * Evaluate a Psuedo-Voigt function
 * @param params [in] Parameters of Psuedo-Voight function. (eta0, eta1, eta2, 2*Theta_k, 
 *  u, V, W, I)
 * @param twoTheta [in] Angle at which function is evaluated
 * @return Value of PS at that angle
 */
inline double ExperimentalPattern::PV(const Vector& params, double twoTheta)
{
	
	// Set variables
	double pi = Constants::pi;
	double Cg  = 4*log(2);
	double dif = twoTheta - params[3];
	double tTT = tan(Num<double>::toRadians(twoTheta / 2));
	double sfw = params[4] + params[5]*tTT + params[6]*tTT*tTT;
	double ex  = exp(-Cg * dif*dif / sfw);
	double eta = params[0] + params[1]*twoTheta + params[2]*twoTheta*twoTheta;
	double den = 1 + 4*dif*dif / sfw;
	
	// Return result
 	return params[7] * (sqrt(Cg) * ex * eta / sqrt(pi * sfw) + \
                2 * (1 - eta) / (pi * sqrt(sfw) * den));
}



/**
 * Calculate Derivatives of Pseudo-Voigt function
 * @param params [in] Parameters of Psuedo-Voight function. (eta0, eta1, eta2, 2*Theta_k, 
 *  u, V, W, I)
 * @param twoTheta [in] Angle at which derivatives is evaluated
 * @return Derivatives wrt each parameter
 */
inline Vector ExperimentalPattern::PVderivs(const Vector& params, double twoTheta)
{
	
	// Set variables
	double pi = Constants::pi;
	double Cg  = 4*log(2);
	double dif = twoTheta - params[3];
	double tTT = tan(Num<double>::toRadians(twoTheta / 2));
	double sfw = params[4] + params[5]*tTT + params[6]*tTT*tTT;
	double ex  = exp(-Cg * dif*dif / sfw);
	double eta = params[0] + params[1]*twoTheta + params[2]*twoTheta*twoTheta;
	double den = 1 + 4*dif*dif / sfw;
	
	// Derivatives
	Vector res(8);
	
	// Derivatives wrt eta0, eta1, eta2
	res[0] = params[7] * (sqrt(Cg) * ex / sqrt(pi * sfw) - 2 / (pi * sqrt(sfw) * den));
	res[1] = twoTheta * res[0];
	res[2] = twoTheta * res[1];
	
	// Derivative wrt 2*Theta_k
	double sfw32 = pow(sfw, 1.5);
	double den2 = den*den;
	res[3] = params[7] * (2*pow(Cg,1.5) * ex * eta * dif / (sqrt(pi) * sfw32) + 16*(1 - eta)*dif / (pi * sfw32 * den2));
	
	// Derivatives wrt u, v, w
	double sfw52 = pow(sfw, 2.5);
	double term1 = pow(Cg,1.5) * ex * eta * dif * dif / (sqrt(pi) * sfw52);
	double term2 = sqrt(Cg) * ex * eta / (2 * sqrt(pi) * sfw32);
	double term3 = 8 * (1 - eta) * dif * dif / (pi * sfw52 * den2);
	double term4 = (1 - eta) / (pi * sfw32 * den);
	res[4] = params[7] * (term1 - term2 + term3 - term4);
	res[5] = tTT * res[4];
	res[6] = tTT * res[5];
	
	// Derivative wrt I
	res[7] = sqrt(Cg) * ex * eta / sqrt(pi * sfw) + 2*(1 - eta) / (pi * sqrt(sfw) * den);
	
	// Return result
	return res;
}



/**
 * Calculate derivative of Psuedo-Voight wrt twoTheta
 * @param params [in] Parameters of Psuedo-Voight function. (eta0, eta1, eta2, 2*Theta_k, 
 *  u, V, W, I)
 * @param twoTheta [in] Angle at which derivative is evaluated 
 * @return Value of derivative wrt twoTheta
 */
inline double ExperimentalPattern::PVderiv(const Vector& params, double twoTheta)
{
	
	// Set variables
	double pi = Constants::pi;
	double Cg  = 4*log(2);
	double TT2 = Num<double>::toRadians(twoTheta / 2);
	double dif = twoTheta - params[3];
	double tTT = tan(TT2);
	double cTT = 1.0/cos(TT2);
	double sfw = params[4] + params[5]*tTT + params[6]*tTT*tTT;
	double sct = pi*cTT*cTT/180 * (params[5]/2 + params[6]*tTT);
	double ex  = exp(-Cg * dif*dif / sfw);
	double eta = params[0] + params[1]*twoTheta + params[2]*twoTheta*twoTheta;
	double den = 1 + 4*dif*dif / sfw;
	
	// Return result
	return params[7] * \
		(-sqrt(Cg) * ex * pi * eta * sct / (2 * pow(pi * sfw, 1.5)) + \
		sqrt(Cg) * ex * (params[1] + 2*params[2]*twoTheta) / sqrt(pi * sfw) + \
		pow(Cg, 1.5) * dif * ex * eta * (dif*sct / sfw - 2) / (sfw * sqrt(pi*sfw)) - \
		8 * (1 - eta) * dif * (2 - dif*sct / sfw) / (pi * pow(sfw, 1.5) * den*den) - \
		(1 - eta) * sct / (pi * pow(sfw, 1.5) * den) + \
		2 * (-params[1] - 2*params[2]*twoTheta) / (pi * sqrt(sfw) * den));
}


/**
 * Calculate the composite of multiple pseudo-Voight functions
 * <p>The first 8 parameters are for the first function, the second 8 are for the 
 * second, and so on.
 * @param params [in] Parameters of PV functions (eto0_1, ..., I_1, eta0_1)
 * @param twoTheta [in] Angle at which function is evaluated 
 * @return Value of multiple PVs at that angle
 */
inline double ExperimentalPattern::compositePV(const Vector& params, double twoTheta) {
    double output = 0;
    for (int f = 0; f < params.length()/8; f++) {
        Vector subParams(8);
        for (int i=0; i<8; i++) subParams[i] = params[f * 8 + i];
        output += PV(subParams, twoTheta);
    }
    return output;
}

/**
 * Calculate derivatives of the composite of multiple pseudo-Voight functions
 * <p>The first 8 parameters are for the first function, the second 8 are for the 
 * second, and so on.
 * @param params [in] Parameters of PV functions (eto0_1, ..., I_1, eta0_1)
 * @param twoTheta [in] Angle at which function is evaluated 
 * @return Value of derivatives of each parameter at that angle
 */ 
inline Vector ExperimentalPattern::compositePVDerivs(const Vector& params, double twoTheta) {
    Vector output(params.length(), 0);
    for (int f = 0; f < params.length()/8; f++) {
        Vector subParams(8);
        for (int i=0; i<8; i++) subParams[i] = params[f * 8 + i];
        Vector subDerivs = PVderivs(subParams, twoTheta);
        for (int i=0; i<8; i++) output[f * 8 + i] = subDerivs[i];
    }
    return output;
}

#endif
