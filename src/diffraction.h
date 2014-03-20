/* Copyright 2011-2014 Kyle Michel, Logan Ward
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



// Class to generate diffraction peaks
class Diffraction
{
	
public:
	
	// Evaluation methods
	enum Method {DM_XRAY, DM_NEUTRON, DM_SIMPLE};
	
	// Structure to store a single peak
	struct Peak;
	
private:
	
	// Methods for calculating R factors
	enum Rmethod {DR_ABS, DR_SQUARED};
	
	// ==============================
	// General variables and methods
	// ==============================
	
		// Variables
		List<double> _twoTheta;
		List<double> _intensity;
	
		// Settings
		Method _method;
		double _wavelength;
		double _minTwoTheta;
		double _maxTwoTheta;
		double _minRelativeIntensity;
		double _fwhm;
		double _variance;
		double _fitAccuracyRFactor;
		double _fitAccuracyDerivs;
		double _minBFactor;
		double _maxBFactor;
		static double _resolution;
		
		// Functions
		void sortByTwoTheta(int left, int right);
	
	
	// ====================================
	// Calculate a pattern for a structure
	// ====================================
	
		// Variable to store working peaks
		OList<Peak>::D2 _peaks;
		
		// Variable to store atomic form factor paramaters
		List<double>::D2 _atfParams;
		
		// Calculate peaks
		void setATFParams(const Symmetry& symmetry);
		void setPeakLocations(const ISO& iso, const Symmetry& symmetry);
		void setPeakIntensities(const Symmetry& symmetry, const Vector& Bfactors, bool getBderivs, \
			bool getPosDerivs, double scale = 1);
		void scaleIntensity(double max);
		void extractIntensities();
		void matchPeaksToReference(const Diffraction& reference);
		double optIntensity(const Diffraction& reference, double* maxIntensity = 0);
		
		// Helper functions
		double thermalFactor(double angle, double Bfactor = 1);
		double thermalFactorDeriv(double angle, double Bfactor = 1);
		double lpFactor(double angle);
		double angle(const Basis& basis, const Vector3D& hkl);
		double structureFactorSquared(const Symmetry& symmetry, double angle, const Vector3D& hkl, \
				const Vector& Bfactors, Vector& Bderivs, Vector& posDerivs);
		double atomicScatteringFactor(int index, double angle);
	
	
	// ===================================
	// Optimize pattern against reference
	// ===================================
	
		// Compare patterns
		double optStructure(const Diffraction& reference, const Symmetry& symmetry, bool optBfactors, \
			bool optPositions, bool showWarnings, const Vector3D* lengths = 0);
		bool isConverged(double value)	{ return (Num<double>::abs(value) < _fitAccuracyDerivs); }
		bool isConverged(Vector& values, Vector* Bfactors = 0);
		void constrainBFactors(Vector& values);
		void setPositions(const Symmetry& symmetry, const Vector& positions);
		void symDerivatives(const Symmetry& symmetry, Vector& derivs);
		double curRFactor(const Diffraction& reference, Rmethod method, double* scale = 0, double* scaleDeriv = 0, \
			Vector* Bderivs = 0, Vector* posDerivs = 0);
		void setNewIntensities(double oldMax, double newMax);
		
	
	// =========================================
	// Stuff related to processing raw patterns
	// =========================================
	
		// Fitting variables
		Vector _PSparams;
	
		// Analyze a raw pattern
		void smoothData(Linked<double>& rawTwoTheta, Linked<double>& rawIntensity);
		void removeBackground(Linked<double>& rawTwoTheta, Linked<double>& rawIntensity);
		void getPeaks(List<double>::D2& peakTwoTheta, List<double>::D2& peakIntensity, \
				const Linked<double>& rawTwoTheta, const Linked<double>& rawIntensity);
		void fitPeaks(const List<double>::D2& peakTwoTheta, const List<double>::D2& peakIntensity);
		static double smoothPoint(const double* twoTheta, const double* intensity, double scale, int numPoints);
	
		// Functions for fitting Gaussian function
		// Params order H, 2*theta_k, I0
		double gaussian(const Vector& params, double twoTheta);
		Vector gaussianDerivs(const Vector& params, double twoTheta);
	
		// Functions for fitting pseudo-Voigt function
		// Params order: eta0, eta1, eta2, 2*theta_k, u, v, w, I0
		double PStwoTheta(double twoTheta)		{ return PS(_PSparams, twoTheta); }
		double PS(const Vector& params, double twoTheta);
		Vector PSderivs(const Vector& params, double twoTheta);
		double PSderiv(const Vector& params, double twoTheta);
		double PSderiv(double twoTheta)			{ return PSderiv(_PSparams, twoTheta); }
	
public:
	
	// Constructors
	Diffraction();
	
	// Set settings
	void method(Method input)				{ _method = input; }
	void wavelength(double input)			{ _wavelength = input; }
	void minTwoTheta(double input)			{ _minTwoTheta = input; }
	void maxTwoTheta(double input)			{ _maxTwoTheta = input; }
	void minRelativeIntensity(double input)	{ _minRelativeIntensity = input; }
	void fwhm(double input)					{ _fwhm = input; _variance = pow(input / (2 * sqrt(2 * log(2))), 2); }
	void variance(double input)				{ _variance = input; _fwhm = 2 * sqrt(2 * log(2) * input); }
	void fitAccuracyRFactor(double input)	{ _fitAccuracyRFactor = input; }
	void fitAccuracyDerivs(double input)	{ _fitAccuracyDerivs = input; }
	
	// Setup functions
	void clear();
	
	// Check if file is correct format
	static bool isFormat(const Text& text);
	static bool isFormat(const Word& file)	{ return isFormat(Read::text(file)); }
	
	// Set from file, existing data, or a structure
	void set(const Text& text);
	void set(const Word& file)	{ set(Read::text(file)); }
	void set(const Linked<double>& twoTheta, const Linked<double>& intensity);
	double set(const ISO& iso, const Symmetry& symmetry, const Diffraction* ref = 0, bool fitBfactors = false);
	
	// Refine a structure
	double refine(ISO& iso, Symmetry& symmetry, const Diffraction& reference, bool showWarnings = true);
	
	// Get R factor for two patterns
	double rFactor(const Diffraction& reference);
	
	// Print functions
	void print(const Word& file, bool broaden = false) const;
	
	// Access functions
	bool isSet() const						{ return (_twoTheta.length() > 0); }
	double wavelength() const				{ return _wavelength; }
	const List<double>& twoTheta() const	{ return _twoTheta; }
	const List<double>& intensity() const	{ return _intensity; }
};



// Structure to store a diffraction pattern peak
struct Diffraction::Peak
{
	
	// General variables
	int patternIndex;
	double twoThetaDeg;
	double twoThetaRad;
	double intensity;
	
	// Variables to speed up calculations
	double lpFactor;
	double multiplicity;
	Vector3D hkl;
	
	// Derivatives
	Vector derivBfactors;
	Vector derivPositions;
	
	// Functions
	Peak& operator= (const Peak& rhs)
	{
		if (this != &rhs)
		{
			patternIndex = rhs.patternIndex;
			twoThetaDeg = rhs.twoThetaDeg;
			twoThetaRad = rhs.twoThetaRad;
			intensity = rhs.intensity;
			lpFactor = rhs.lpFactor;
			multiplicity = rhs.multiplicity;
			hkl = rhs.hkl;
			derivBfactors = rhs.derivBfactors;
			derivPositions = rhs.derivPositions;
		}
		return *this;
	}
};



// =====================================================================================================================
// Diffraction
// =====================================================================================================================

/* inline Diffraction::Diffraction()
 *
 * Constructor for Diffraction object
 */

inline Diffraction::Diffraction()
{
	_method = DM_XRAY;
	_wavelength = 1.5418;
	_minTwoTheta = 10;
	_maxTwoTheta = 100;
	_minRelativeIntensity = 1e-3;
	fwhm(0.1);
	_fitAccuracyRFactor = 1e-4;
	_fitAccuracyDerivs = 1e-5;
	_minBFactor = 0.1;
	_maxBFactor = 4.0;
}



/* inline void Diffraction::clear()
 *
 * Clear data in Diffraction object
 */

inline void Diffraction::clear()
{
	_twoTheta.clear();
	_intensity.clear();
	_peaks.clear();
	_atfParams.clear();
}



/* inline double Diffraction::lpFactor(double angle)
 *
 * Return the Lorentz polarization factor
 */

inline double Diffraction::lpFactor(double angle)
{
	return (1 + pow(cos(2 * angle), 2)) / (cos(angle) * pow(sin(angle), 2));
}



/* inline double Diffraction::thermalFactor(double Bfactor)
 *
 * Return the thermal factor
 */

inline double Diffraction::thermalFactor(double angle, double Bfactor)
{
	return exp(-Bfactor * pow(sin(angle) / _wavelength, 2.0));
}



/* inline double Diffraction::thermalFactorDeriv(double angle, double Bfactor)
 * 
 * Return the derivative of the thermal factor wrt Bfactor
 */

inline double Diffraction::thermalFactorDeriv(double angle, double Bfactor)
{
	double temp = sin(angle) / _wavelength;
	temp *= temp;
	return -temp * exp(-Bfactor * temp);
}



/* inline double Diffraction::angle(const Basis& basis, const Vector3D& hkl)
 *
 * Return the angle of a vector
 */

inline double Diffraction::angle(const Basis& basis, const Vector3D& hkl)
{
	double arg = (basis.inverse() * hkl).magnitude() * _wavelength / 2;
	if ((arg >= -1) && (arg <= 1))
		return asin(arg);
	else if (arg < -1)
		return -Constants::pi / 2;
	return Constants::pi / 2;
}



/* inline bool Diffraction::isConverged(Vector& values, Vector* Bfactors)
 *
 * Return whether all values are within tolerance
 */

inline bool Diffraction::isConverged(Vector& values, Vector* Bfactors)
{
	for (int i = 0; i < values.length(); ++i)
	{
		if (Num<double>::abs(values[i]) > _fitAccuracyDerivs)
		{
			if (Bfactors)
			{
				if ((Num<double>::neq((*Bfactors)[i], _minBFactor, 1e-6)) && \
					(Num<double>::neq((*Bfactors)[i], _maxBFactor, 1e-6)))
					return false;
			}
			else
				return false;
		}
	}
	return true;
}



/* inline void Diffraction::constrainBFactors(Vector& values)
 *
 * Make sure that all thermal factors (B) are within set range
 */

inline void Diffraction::constrainBFactors(Vector& values)
{
	int i;
	for (i = 0; i < values.length(); ++i)
	{
		values[i] = values[i] < _minBFactor ? _minBFactor : values[i];
		values[i] = values[i] > _maxBFactor ? _maxBFactor : values[i];
	}
}



/* inline double Diffraction::gaussian(const Vector& params, double twoTheta)
 *
 * Gaussian function
 */

inline double Diffraction::gaussian(const Vector& params, double twoTheta)
{
	
	// Set variables
	double pi = Constants::pi;
	double Cg  = 4*log(2);
	double dif = twoTheta - params[1];
	double ex  = exp(-Cg * dif*dif / params[0]);
	
	// Return result
	return params[2] * sqrt(Cg) * ex / (sqrt(pi) * params[0]);
}



/* inline Vector Diffraction::gaussianDerivs(const Vector& params, double twoTheta)
 *
 * Derivatives of gaussian function
 */

inline Vector Diffraction::gaussianDerivs(const Vector& params, double twoTheta)
{
	
	// Set variables
	double pi = Constants::pi;
	double Cg  = 4*log(2);
	double dif = twoTheta - params[1];
	double ex  = exp(-Cg * dif*dif / params[0]);
	
	// Get derivatives
	Vector res(3);
	res[0] = params[2] * sqrt(Cg) * (Cg*dif*dif - params[0]) * ex / (sqrt(pi) * pow(params[0], 3.0));
	res[1] = 2 * pow(Cg, 1.5) * params[2] * dif * ex / (sqrt(pi) * params[0]*params[0]);
	res[2] = sqrt(Cg) * ex / (sqrt(pi) * params[0]);
	
	// Return derivatives
	return res;
}



/* inline double Diffraction::PS(const Vector& params, double twoTheta)
 *
 * Psuedo-Voigt function
 */

inline double Diffraction::PS(const Vector& params, double twoTheta)
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
 	return params[7] * (sqrt(Cg) * ex * eta / sqrt(pi * sfw) + 2*(1 - eta) / (pi * sqrt(sfw) * den));
}



/* inline Vector Diffraction::PSderivs(const Vector& params, double twoTheta)
 *
 * Derivatives of Pseudo-Voigt function
 */

inline Vector Diffraction::PSderivs(const Vector& params, double twoTheta)
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



/* inline double Diffraction::PSderiv(const Vector& params, double twoTheta)
 *
 * Return value of derivative wrt two theta
 */

inline double Diffraction::PSderiv(const Vector& params, double twoTheta)
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



#endif
