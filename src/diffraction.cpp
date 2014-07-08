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

#include "multi.h"
#include "diffraction.h"
#include "language.h"
#include "output.h"
#include "random.h"
#include "launcher.h"
#include "timer.h"
#include "dlib/numerical_integration.h"
#include <cstdlib>
#include <vector>
#include <deque>
#include <algorithm>
#include <numeric>


double Diffraction::_resolution = 1e-2;

/**
 * Define the structure used to generate a diffraction pattern. Also populates 
 * a list of atoms whose positions can be refined, and calculates locations at which
 * diffraction peaks will appear. 
 * 
 * Note: This operation must be run before any sort of diffraction pattern calculation
 * is run.
 * @param structure [in] Structure to be used in calculations / refinements
 * @param symmetry [in] Symmetry information of that structure
 */
void Diffraction::defineStructure(const ISO& structure, const Symmetry& symmetry) {
    // Clear out any old data
    _symmetry = &symmetry;
    _structure = &structure;

    // Get the atomic parameters 
    setATFParams();
    calculatePeakLocations();

    // Initialize guesses for fitting parameters
    initializeRefinementParameters();
    
    // If a reference pattern is set, match peas to it
    if (referencePatternIsDefined())
        matchPeaksToReference();
}

/**
 * Check whether user has defined the structure from which to calculate the 
 *  diffraction pattern.
 * @return Whether defineStructure has been run
 */
bool Diffraction::structureIsDefined() {
    return _symmetry != 0;
}

void Diffraction::defineReferencePattern(const Diffraction& reference) {
    _referencePattern = &reference;
    // Retrieve parameters of reference pattern
    _method = reference._method;
    _wavelength = reference._wavelength;
    _minTwoTheta = reference._minTwoTheta;
    _maxTwoTheta = reference._maxTwoTheta;
    // If a structure is defined, redefine peak locations
    if (structureIsDefined())
        calculatePeakLocations();
    // Match peaks to reference pattern
    matchPeaksToReference();
}

/**
 * Whether the user has defined the reference pattern against which the structure 
 *  is refined.
 */
bool Diffraction::referencePatternIsDefined() {
    return _referencePattern != 0;
}

/**
 * Make initial guesses for refinement parameters based on the provided structure. 
 * 
 * Currently, this retrieves the internal degrees of freedom from the structure.
 * @param structure [in] Structure used to calculate diffraction pattern
 * @param symmetry [in] Symmetry information of structure
 */
void Diffraction::initializeRefinementParameters() {
    // Set initial guess for B factors to 1.0
    _BFactors.resize(_symmetry->orbits().length());
    fill(_BFactors.begin(), _BFactors.end(), 0.5);
    // Allocate position / B factor derivatives
    _BFactorDerivs.resize(_BFactors.size(), 0.0);
    _PositionDerivs.resize(_BFactorDerivs.size() * 3, 0.0);
    // Note: Internal parameters are kept in _symmetry
}

/**
 * Determine whether a specific parameter is in the set of parameters to be 
 *  refined.
 * @param parameter Parameter in question
 * @param toRefine Set of parameters that will be refined
 * @return Whether this parameter will be refined
 */
bool Diffraction::willRefine(RefinementParameters parameter, std::set<RefinementParameters> toRefine) {
    return toRefine.find(parameter) != toRefine.end();
}

/**
 * Calculate the diffraction pattern of a structure and store in this object. If provided with
 * a reference pattern and <code>fitBfactors</code> is true, will fit B factors to best
 * match the provided pattern.
 * 
 * @param iso [in] Structure from which to calculate diffraction pattern
 * @param symmetry [in] Describes symmetry of structure 
 * @param ref [in] Reference pattern to fit intensities and (if desired) BFactors against. Equals 0 if no pattern supplied
 * @param fitBFactors [in] Whether to fit B factors (if pattern provided)
 * @return R factor (if pattern provided)
 */
double Diffraction::set(const ISO& iso, const Symmetry& symmetry, const Diffraction* ref, bool fitBfactors) {
    // Clear space
    clear();

    // Set this structure as the structure to be refined
    defineStructure(iso, symmetry);

    // Define the reference pattern
    if (ref != 0)
        defineReferencePattern(*ref);

    // Output
    Output::newline();
    Output::print("Calculating peak intensities for the structure");
    Output::increase();

    // If a reference pattern was passed then optimize max intensity and B factors
    double rFactor = 0;
    if (ref) {
        // Output
        Output::newline();
        Output::print("Optimizing against reference pattern");
        Output::increase();

        // Refine the scale factor and B factors
        std::set<RefinementParameters> toRefine;
        if (fitBfactors)
            toRefine.insert(RF_BFACTORS);
        refineParameters(toRefine);

        rFactor = getCurrentRFactor(DR_ABS);

        // Output
        Output::newline();
        Output::print("Optimal R factor: ");
        Output::print(rFactor);
        Output::decrease();
    }
        // Calculate the peak intensities using initial guesses
    else {
        calculatePeakIntensities();
    }

    // Extract the intensities
    extractIntensities();
    
    // Scale the intensities to 1000
    _integratedIntensity *= 1000 / _integratedIntensity.max();

    // Print intensities
    Output::newline();
    Output::print("Generated ");
    Output::print(_peakTwoTheta.size());
    Output::print(" peak");
    if (_peakTwoTheta.size() != 1)
        Output::print("s");
    Output::increase();
    for (int i = 0; i < _peakTwoTheta.size(); ++i) {
        if (_integratedIntensity[i] < 50) 
            continue;
        Output::newline();
        Output::print("Two-theta and intensity of ");
        Output::print(_peakTwoTheta[i]);
        Output::print(" ");
        Output::print(_integratedIntensity[i]);
    }
    Output::decrease();

    // Output
    Output::decrease();

    // Return the R factor
    return rFactor;
}

/**
 * Refine a structure (which has already been stored internally) against a reference
 *  diffraction pattern
 * @param toRefine [in] What parameters of the diffracted intensity should be refined
 */
void Diffraction::refineParameters(std::set<RefinementParameters> toRefine) {
    if (!referencePatternIsDefined()) {
        Output::newline(ERROR);
        Output::print("Internal Error: Reference pattern not yet defined.");
    }
    if (!structureIsDefined()) {
        Output::newline(ERROR);
        Output::print("Internal Error: Structure not yet defined.");
    }
    Output::increase();

    // Clear what parameters are currently being refined
    _currentlyRefining.clear();

    // Optimize atomic positions, if desired
    if (willRefine(RF_POSITIONS, toRefine)) {
        Output::newline();
        Output::print("Refining atomic positions. Current R Factor: ");
        _currentlyRefining.insert(RF_POSITIONS);
        double curRFactor = runRefinement();
        Output::print(curRFactor, 3);
    }

    // Next, refine both B factors and atomic positions
    if (willRefine(RF_BFACTORS, toRefine)) {
        Output::newline();
        Output::print("Also refining isotropic thermal factors. Current R Factor: ");
        _currentlyRefining.insert(RF_BFACTORS);
        double curRFactor = runRefinement();
        Output::print(curRFactor, 3);
    }
    Output::decrease();

    // Output results
    if (willRefine(RF_BFACTORS, toRefine)) {
        for (int i = 0; i < _BFactors.size(); ++i) {
            Output::newline();
            Output::print("Optimized B factor for atom ");
            Output::print(_symmetry->orbits()[i].atoms()[0]->atomNumber() + 1);
            Output::print(" (");
            Output::print(_symmetry->orbits()[i].atoms()[0]->element().symbol());
            Output::print("): ");
            Output::print(_BFactors[i]);
        }
    }
    if (willRefine(RF_POSITIONS, toRefine)) {
        for (int i = 0; i < _symmetry->orbits().length(); ++i) {
            Output::newline();
            Output::print("Optimized position for atom ");
            Output::print(_symmetry->orbits()[i].atoms()[0]->atomNumber() + 1);
            Output::print(" (");
            Output::print(_symmetry->orbits()[i].atoms()[0]->element().symbol());
            Output::print("): ");
            for (int j = 0; j < 3; ++j) {
                Output::print(_symmetry->orbits()[i].atoms()[0]->fractional()[j]);
                if (j != 2)
                    Output::print(", ");
            }
        }
    }
}

/**
 * Called from refineParameters. Refine any parameters currently defined in 
 *  _currentlyRefining.
 * @return Minimal R factor (using DR_ABS)
 */
double Diffraction::runRefinement() {
    column_vector params;
    column_vector x_low = getRefinementParameterLowerBoundary();
    column_vector x_high = getRefinementParameterUpperBoundary();
    params = getRefinementParameters();
    RFactorFunctionModel f(this);
    // Technical issue (as of 12Mar14): 
    //  Derivatives calculated during getCurrentRFactor are wrong. That function 
    //   sets the "optimal scale factor," but assumes it is constant as we adjust 
    //   B factors and another parameters. We need the derivative of the scale factor
    //   with respect to each refinement parameter, which would be tedious. Instead,
    //   approximate derivatives are used, which seem to work well. 
    // RFactorDerivativeFunctionalModel der(this);
    dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
            dlib::objective_delta_stop_strategy(1e-8),
            f, dlib::derivative(f), params, x_low, x_high);
    setAccordingToParameters(params);
    calculatePeakIntensities();
    return getCurrentRFactor(DR_ABS);
}

/**
 * Get a vector representing the parameters which are being refined. Always arranged
 *  in the following order:
 * <ol>
 * <li>Atomic positions</li>
 * <li>Thermal factors</li>
 * </ol>
 * @return Current values of parameters to be optimized
 */
Diffraction::column_vector Diffraction::getRefinementParameters() {
    queue<double> params;
    if (willRefine(RF_POSITIONS, _currentlyRefining)) {
        for (int orbit = 0; orbit < _symmetry->orbits().length(); orbit++)
            for (int dir = 0; dir < 3; dir++) {
                double pos = _symmetry->orbits()[orbit].atoms()[0]->fractional()[dir];
                params.push(pos);
            }
    }
    if (willRefine(RF_BFACTORS, _currentlyRefining)) {
        for (int i = 0; i < _BFactors.size(); i++)
            params.push(_BFactors[i]);
    }

    // Copy parameters to an appropriate container
    int nParams = params.size();
    column_vector output(nParams);
    for (int i = 0; i < nParams; i++) {
        output(i) = params.front();
        params.pop();
    }
    return output;
}

/**
 * Generate the first derivatives of R factor with respect to each current 
 * refinement parameter, in a format that dlib's optimization algorithms can use.
 * 
 * Note: dlib operation algorithms ensure that a call to this algorithm is always preceeded by
 *  a call to currentRFactor, which updates the derivatives for each parameter.
 * 
 * @return A dlib "column_vector"
 */
Diffraction::column_vector Diffraction::getRefinementParameterDerivatives(column_vector x) {
    queue<double> params;
    if (willRefine(RF_POSITIONS, _currentlyRefining))
        for (int i = 0; i < _PositionDerivs.size(); i++)
            params.push(_PositionDerivs[i]);
    if (willRefine(RF_BFACTORS, _currentlyRefining))
        for (int i = 0; i < _BFactors.size(); i++)
            params.push(_BFactorDerivs[i]);
    // Copy parameters to an appropriate container
    int nParams = params.size();
    column_vector output(nParams);
    for (int i = 0; i < nParams; i++) {
        output(i) = params.front();
        params.pop();
    }
    return output;
}

/**
 * Get the lower boundary of each refinement parameter. 
 * @return Column vector detailing lower boundary for each parameter, same order
 *  as getRefinementParameters.
 */
Diffraction::column_vector Diffraction::getRefinementParameterLowerBoundary() {
    queue<double> params;
    if (willRefine(RF_POSITIONS, _currentlyRefining))
        for (int i = 0; i < _PositionDerivs.size(); i++)
            params.push(-1);
    if (willRefine(RF_BFACTORS, _currentlyRefining))
        for (int i = 0; i < _BFactors.size(); i++)
            params.push(_minBFactor);
    // Copy parameters to an appropriate container
    int nParams = params.size();
    column_vector output(nParams);
    for (int i = 0; i < nParams; i++) {
        output(i) = params.front();
        params.pop();
    }
    return output;
}

/**
 * Get the upper boundary of each refinement parameter. 
 * @return Column vector detailing upper boundary for each parameter, same order
 *  as getRefinementParameters.
 */
Diffraction::column_vector Diffraction::getRefinementParameterUpperBoundary() {
    queue<double> params;
    if (willRefine(RF_POSITIONS, _currentlyRefining))
        for (int i = 0; i < _PositionDerivs.size(); i++)
            params.push(2);
    if (willRefine(RF_BFACTORS, _currentlyRefining))
        for (int i = 0; i < _BFactors.size(); i++)
            params.push(_maxBFactor);
    // Copy parameters to an appropriate container
    int nParams = params.size();
    column_vector output(nParams);
    for (int i = 0; i < nParams; i++) {
        output(i) = params.front();
        params.pop();
    }
    return output;
}

/**
 * Set all refinement parameters according to those contained in an input vector.
 *  The values in this vector depending on which parameters are _currentlyRefining. 
 *  Their order should be the same as defined in getCurrentParameters
 *  
 * @param params New values of parameters being refined
 */
void Diffraction::setAccordingToParameters(column_vector params) {
    int position = 0;
    if (willRefine(RF_POSITIONS, _currentlyRefining)) {
        Vector newPositions(_symmetry->orbits().length()*3);
        for (int i = 0; i < newPositions.length(); i++) {
            newPositions[i] = params(position++);
        }
        symPositions(*_symmetry, newPositions);
        setPositions(*_symmetry, newPositions);
    }
    if (willRefine(RF_BFACTORS, _currentlyRefining)) {
        for (int i = 0; i < _BFactors.size(); i++)
            _BFactors[i] = params(position++);
    }
}

/**
 * Refine a structure against reference pattern
 * 
 * @param iso [in,out] Structure to be refined. Returns refined coordinates
 * @param symmetry [in,out] Symmetry information about structure. Returns refined coordinates
 * @param reference [in] Pattern to refine against
 * @param showWarnings [in] Whether to print warnings
 * @return Optimized R Factor
 */
double Diffraction::refine(ISO& iso, Symmetry& symmetry, const Diffraction& reference, bool showWarnings) {
    // Clear out any old information
    clear();
    
    // Store the pattern information  
    defineReferencePattern(reference);
    
    // Store structure and initialize guesses
    defineStructure(iso, symmetry);

    // Output
    Output::newline();
    Output::print("Refining structure against reference pattern");
    Output::increase();

    // Refine B factors, intensity scale, and positions
    std::set<RefinementParameters> toRefine;
    toRefine.insert(RF_BFACTORS);
    toRefine.insert(RF_POSITIONS);

    // Run refinement
    refineParameters(toRefine);
    double rFactor = getCurrentRFactor(DR_ABS);

    // Output
    Output::newline();
    Output::print("Optimal R factor: ");
    Output::print(rFactor);

    // Output
    Output::decrease();

    // Return R factor
    return rFactor;
}

/** 
 * Calculate peaks that will appear in diffraction pattern of a 
 * particular structure. Store them internally in _integratedPeaks.
 */
void Diffraction::calculatePeakLocations() {
    Output::increase();
    // Clear any previously-calculated peaks
    _diffractionPeaks.clear();

    // Calculate the hkl range
    int i, j;
    double range[3];
    double maxMag = 2 * sin(Num<double>::toRadians(_maxTwoTheta / 2)) / _wavelength;
    Vector3D vec;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j)
            vec[j] = _structure->basis().reducedInverse()(j, i);
        range[i] = Num<double>::abs(Num<double>::ceil(maxMag / vec.magnitude()));
    }

    // Conversion matrix to take reduced cell reciprocal lattice vector to unit cell reciprocal lattice vector
    Matrix3D convHKL = _structure->basis().unitPointToReduced().transpose();

    // Generate symmetry operations for reduced cell
    Matrix3D P = _structure->basis().unitToReduced().transpose();
    Matrix3D Q = P.inverse();
    OList<Matrix3D > operations(_symmetry->operations().length());
    for (i = 0; i < _symmetry->operations().length(); ++i) {
        operations[i] = P;
        operations[i] *= _symmetry->operations()[i].rotation();
        operations[i] *= Q;
        operations[i] = operations[i].transpose();
    }

    // Remove identity operation
    Matrix3D identity;
    identity.makeIdentity();
    for (i = 0; i < operations.length(); ++i) {
        if (operations[i] == identity) {
            operations.remove(i);
            break;
        }
    }

    // Get the intrinsic part of all symmetry operations
    OList<Vector3D >::D2 translations(_symmetry->operations().length());
    for (i = 0; i < _symmetry->operations().length(); ++i) {
        translations[i].length(_symmetry->operations()[i].translations().length());
        for (j = 0; j < _symmetry->operations()[i].translations().length(); ++j)
            translations[i][j] = Symmetry::intrinsicTranslation(_symmetry->operations()[i].rotation(), \
				_symmetry->operations()[i].translations()[j]);
    }

    // Get the peak intensities
    bool found;
    int mult;
    double product;
    double twoTheta;
    Vector3D hkl;
    Vector3D redHKL;
    Vector3D symHKL;
    Linked<Vector3D > equivPoints;
    Linked<Vector3D >::iterator itEquiv;
    for (redHKL[0] = -range[0]; redHKL[0] <= range[0]; ++redHKL[0]) {
        for (redHKL[1] = -range[1]; redHKL[1] <= range[1]; ++redHKL[1]) {
            for (redHKL[2] = -range[2]; redHKL[2] <= range[2]; ++redHKL[2]) {
                // Loop over operations to generate equivalent points
                mult = 1;
                equivPoints.clear();
                equivPoints += redHKL;
                for (i = 0; i < operations.length(); ++i) {
                    // Get new point
                    symHKL = operations[i] * redHKL;
                    for (j = 0; j < 3; ++j)
                        symHKL[j] = Num<double>::round(symHKL[j], 1);

                    // Check if hkl is equilvalent of something that has been generated earlier
                    // If so, set mult = 0
                    if (symHKL[0] < redHKL[0] - 1e-4)
                        mult = 0;
                    else if (Num<double>::abs(symHKL[0] - redHKL[0]) < 1e-4) {
                        if (symHKL[1] < redHKL[1] - 1e-4)
                            mult = 0;
                        else if (Num<double>::abs(symHKL[1] - redHKL[1]) < 1e-4) {
                            if (symHKL[2] < redHKL[2] - 1e-4)
                                mult = 0;
                        }
                    }
                    if (mult == 0)
                        break;

                    // Check if point is already known as an equivalent point
                    found = false;
                    for (itEquiv = equivPoints.begin(); itEquiv != equivPoints.end(); ++itEquiv) {
                        if ((Num<double>::abs((*itEquiv)[0] - symHKL[0]) < 1e-4) && \
                                (Num<double>::abs((*itEquiv)[1] - symHKL[1]) < 1e-4) && \
                                (Num<double>::abs((*itEquiv)[2] - symHKL[2]) < 1e-4)) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        ++mult;
                        equivPoints += symHKL;
                    }
                }

                // Multiplier is zero so skip
                if (mult == 0)
                    continue;

                // Convert current reduced basis hkl to unit cell
                hkl = convHKL * redHKL;
                vector<Vector3D> recipLatVectors;
                recipLatVectors.reserve(mult);
                for (itEquiv = equivPoints.begin(); itEquiv != equivPoints.end(); ++itEquiv) {
                    *itEquiv = convHKL * *itEquiv;
                    recipLatVectors.push_back(_structure->basis().inverse() * *itEquiv);
                }

                // Check if direction will be a systematic absence
                found = false;
                for (i = 0; i < _symmetry->operations().length(); ++i) {

                    // Check if R*hkl = hkl
                    symHKL = _symmetry->operations()[i].rotation() * hkl;
                    if ((Num<double>::abs(symHKL[0] - hkl[0]) > 1e-4) || \
                            (Num<double>::abs(symHKL[1] - hkl[1]) > 1e-4) || \
                            (Num<double>::abs(symHKL[2] - hkl[2]) > 1e-4))
                        continue;

                    // Loop over intrinsic translations and check if ti*hkl = integer for all
                    for (j = 0; j < translations[i].length(); ++j) {
                        product = translations[i][j] * hkl;
                        if (Num<double>::abs(Num<double>::round(product, 1) - product) > 1e-4) {
                            found = true;
                            break;
                        }
                    }

                    // Break if hkl is a system absence
                    if (found)
                        break;
                }

                // Current hkl will be a systematic absence (keep just in case)
                // if (found)
                   //continue;

                // Get the current angle
                twoTheta = 2 * Num<double>::toDegrees(getDiffractionAngle(_structure->basis(), hkl));

                // Skip if diffraction angle is too small or too large
                if ((twoTheta < _minTwoTheta) || (twoTheta > _maxTwoTheta))
                    continue;
                
                // Add peak to list of known peaks
                Peak newPeak(this, twoTheta, mult, hkl, recipLatVectors);
                this->_diffractionPeaks.push_back(newPeak);
            }
        }
    }

    // Ensure peaks are in proper order
    std::sort(_diffractionPeaks.begin(), _diffractionPeaks.end());
    
    Output::newline();
    Output::print("Total number of peaks: ");
    Output::print(_diffractionPeaks.size());
	Output::decrease();
}

/**
 * Calculate the peak intensities, given complete information about a structure. Stores
 * resulting peak heights internally. 
 */
void Diffraction::calculatePeakIntensities() {
    for (int i = 0; i < _diffractionPeaks.size(); i++) {
        _diffractionPeaks[i].updateCalculatedIntensity();
    }
}

/**
 * Update the calculated integrated intensity for a single peak. Also, calculate derivatives
 *  of peak integrated intensity with respect to each parameter
 */
void Diffraction::Peak::updateCalculatedIntensity() {
    // Set derivatives to zero
    derivBfactors.resize(_sourcePattern->_BFactors.size());
    derivBfactors = 0;
    derivPositions.resize(3 * _sourcePattern->_BFactors.size());
    derivPositions = 0;

    // Calculate integrated intensity. (Note abscence of scale factor, which 
    //  is always optimized when calculating R factor)
    peakIntensity = structureFactorSquared(*(_sourcePattern->_symmetry), _twoThetaRad / 2,
            hkl, _sourcePattern->_BFactors, derivBfactors, derivPositions);
    // Everything but the structure factor
    double otherFactors = multiplicity * lpFactor;
    peakIntensity *= otherFactors;
    // otherFactors *= getTexturingFactor(preferredOrientation, 1.666, _peaks[i][j].recipLatticeVectors);

    // Update B and position factor derivatives (which are currently derivatives of structure
    //  factor) to those of the peak intensity. Convert simply by multiply all other factors
    derivBfactors *= otherFactors;
    derivPositions *= otherFactors;
}

/** 
 * Extract calculated peak intensities. Peaks that are closer than a certain 
 *  tolerance or match the same peak in a reference pattern are combined.
 * 
 * Results are stored internally in _peakTwoTheta and _integratedIntensity
 */
void Diffraction::extractIntensities() {
    // Build a list of intensities
    vector<double> tempTwoTheta; tempTwoTheta.reserve(_diffractionPeaks.size());
    vector<double> tempIntensity; tempIntensity.reserve(_diffractionPeaks.size());
    tempTwoTheta.push_back(_diffractionPeaks[0].TwoThetaDeg);
    tempIntensity.push_back(_diffractionPeaks[0].peakIntensity);
    if (referencePatternIsDefined()) {
        // Combine peaks with the same pattern index
        int lastPatternIndex = _diffractionPeaks[0].patternIndex;
        for (int i=1; i<_diffractionPeaks.size(); i++) {
            if (_diffractionPeaks[i].patternIndex == -1 || 
                    _diffractionPeaks[i].patternIndex != lastPatternIndex) {
                tempTwoTheta.push_back(_diffractionPeaks[i].TwoThetaDeg);
                tempIntensity.push_back(_diffractionPeaks[i].peakIntensity);
            } else 
                tempIntensity.back() += _diffractionPeaks[i].peakIntensity;
        }
    } else {
        // Combine peaks closer than 0.15 degrees
        int lastAngle = -100;
        for (int i=1; i<_diffractionPeaks.size(); i++) {
            if (_diffractionPeaks[i].TwoThetaDeg - lastAngle > 0.15) {
                tempTwoTheta.push_back(_diffractionPeaks[i].TwoThetaDeg);
                tempIntensity.push_back(_diffractionPeaks[i].peakIntensity);
                lastAngle = tempTwoTheta.back();
            } else {
                tempIntensity.back() += _diffractionPeaks[i].peakIntensity;
            }
        }
    }
    
    // Save peaks
    _peakTwoTheta.resize(tempTwoTheta.size());
    _integratedIntensity.resize(tempTwoTheta.size());
    for (int i = 0; i < tempTwoTheta.size(); ++i) {
        _peakTwoTheta[i] = tempTwoTheta[i];
        _integratedIntensity[i] = tempIntensity[i];
    }
}

/**
 * Given a diffraction pattern, identifies which calculated peaks match up to 
 *  each peak in the stored reference pattern. Stores matching peak in _diffractionPeak data structure.
 * 
 * This function also creates a "matching peaks" object, which contains the indexes of
 *  peaks in this pattern that match each peak in the reference pattern. While that may
 *  be a mouthful to explain, it will be very useful when calculating R factors, 
 *  which rely on the fact that several peaks from this pattern may match to a single
 *  peak in the reference.
 */
void Diffraction::matchPeaksToReference() {
    // Return if there are no reference peaks
    if (_referencePattern == 0) {
        Output::newline(ERROR);
        Output::print("Internal Error: Reference pattern not set");
        return;
    }
    
    // Reset lookup table for which peaks in this pattern match to a each 
    //  peak in the reference
    _matchingPeaks.clear(); 
    _matchingPeaks.reserve(_referencePattern->_diffractionPeaks.size());
    for (int i=0; i<_referencePattern->_diffractionPeaks.size(); i++) {
        vector<int> empty;
        _matchingPeaks.push_back(empty);
    }
    
    // Clear list of peaks that don't match anything
    _unmatchedPeaks.clear();

    // Tolerance for peaks to be aligned
    double tol = 0.15;

    // Loop over diffraction peaks in this pattern
    for (int thisPeak = 0; thisPeak < _diffractionPeaks.size(); ++thisPeak) {
        // Find the nearest diffraction peak in reference pattern
        int nearIndex = 0;
        double nearDif = abs(_diffractionPeaks[thisPeak].TwoThetaDeg - \
                _referencePattern->_diffractionPeaks[thisPeak].TwoThetaDeg);
        for (int refPeak = 1; refPeak < _referencePattern->_diffractionPeaks.size(); ++refPeak) {
            double curDif = abs(_diffractionPeaks[thisPeak].TwoThetaDeg \
                - _referencePattern->_diffractionPeaks[refPeak].TwoThetaDeg);
            if (curDif < nearDif) {
                nearIndex = refPeak;
                nearDif = curDif;
            }
        }

        // If no peak within tolerance of this peak mark its index as -1 (not matched)
        //  and move on to the next peak
        if (nearDif > tol) {
            _diffractionPeaks[thisPeak].patternIndex = -1;
            _unmatchedPeaks.push_back(thisPeak);
            continue;
        }

        // Otherwise, save index of peak in reference as the pattern index for this peak
        _diffractionPeaks[thisPeak].patternIndex = nearIndex;
        
        // Add this peak to the list of those that match nearIndex
        _matchingPeaks[nearIndex].push_back(thisPeak);
    }
}

/**
 * Calculate the squared structure factor for all atoms in a crystal for certain 
 *  plane at a specific anlge.
 * @param symmetry [in] Symmetrical information about a structure
 * @param angle [in] Angle at which radiation is reflected
 * @param hkl [in] Crystallographic plane being considered
 * @param BFactors [in] Thermal factors for each orbit of atoms
 * @param Bderivs [out] Derivatives of squared structure factor wrt to each B factor
 * @param posDerivates [out] Derivatives of squared structure factor wrt to each B factor
 * @return Squared structure factor for this condition
 */
double Diffraction::Peak::structureFactorSquared(const Symmetry& symmetry, double angle, \
        const Vector3D& hkl, vector<double> BFactors, valarray<double>& Bderivs, valarray<double>& posDerivs) {

    // Variables to store real and imaginary parts of derivatives
    valarray<double> realBderivs(0.0, Bderivs.size());
    valarray<double> imagBderivs(0.0, Bderivs.size());
    valarray<double> realPosDerivs(0.0, posDerivs.size());
    valarray<double> imagPosDerivs(0.0, posDerivs.size());

    // Loop over all atoms and calculate magnitude squared
    int i, j, k;
    double dot;
    double pre;
    double real = 0;
    double imag = 0;
    double sinTerm;
    double cosTerm;
    double preBderiv;
    double prePosDeriv;
    double thermFactor;
    double scatteringFactor;
    double thermFactorDeriv = 0;
    Atom* curAtom;
	for (i = 0; i < symmetry.orbits().length(); ++i) {

		// Get scaling factors for current atom
		curAtom = symmetry.orbits()[i].atoms()[0];
		scatteringFactor = atomicScatteringFactor(i, angle);
		thermFactor = (_sourcePattern->_method == DM_SIMPLE) ? 1 : thermalFactor(angle, BFactors[i]);
		if (Bderivs.size() > 0)
			thermFactorDeriv = (_sourcePattern->_method == DM_SIMPLE) ? 0 : thermalFactorDeriv(angle, BFactors[i]);

		// Loop over equivalent atoms
		for (j = 0; j < symmetry.orbits()[i].atoms().length(); ++j) {

			// Calculate contribution from current atom
			curAtom = symmetry.orbits()[i].atoms()[j];
			dot = 2 * Constants::pi * (hkl * curAtom->fractional());
			sinTerm = sin(dot);
			cosTerm = cos(dot);

			// Add contribution
			pre = scatteringFactor * thermFactor * curAtom->occupancy();
			real += pre * cosTerm;
			imag += pre * sinTerm;

			// Add the B factor derivative if passed
			if (Bderivs.size() > 0) {
				preBderiv = scatteringFactor * thermFactorDeriv * curAtom->occupancy();
				realBderivs[i] += preBderiv * cosTerm;
				imagBderivs[i] += preBderiv * sinTerm;
			}

			// Add position derivatives if passed
			if ((j == 0) && (posDerivs.size() > 0)) {
				for (k = 0; k < 3; ++k) {
					prePosDeriv = pre * 2 * Constants::pi * hkl[k] * symmetry.orbits()[i].atoms().length();
					realPosDerivs[3 * i + k] += prePosDeriv * sinTerm;
					imagPosDerivs[3 * i + k] += prePosDeriv * cosTerm;
				}
			}
		}
	}

	// Save B factor derivatives
	for (i = 0; i < Bderivs.size(); ++i)
		Bderivs[i] = 2 * (real * realBderivs[i] + imag * imagBderivs[i]);

	// Save position derivatives
	if (posDerivs.size() > 0) {
		Vector3D tempDeriv;
		for (i = 0; i < posDerivs.size(); ++i)
			posDerivs[i] = 2 * (-real * realPosDerivs[i] + imag * imagPosDerivs[i]);
		for (i = 0; i < symmetry.orbits().length(); ++i) {
			for (j = 0; j < 3; ++j)
				tempDeriv[j] = posDerivs[3 * i + j];
			tempDeriv *= symmetry.orbits()[i].specialPositions()[0].rotation();
			if (symmetry.orbits()[i].anyAtomsFixed())
				tempDeriv = 0.0;
			for (j = 0; j < 3; ++j)
				posDerivs[3 * i + j] = tempDeriv[j];
		}
	}

    // Return the square of the magnitude
    return real * real + imag*imag;
}

/**
 * Given a list of angles and intensities, determine locations and intensities of 
 *  diffraction peaks. Store peak information internally.
 * 
 * If the intensities are uniformly spaced, this function assumes that it was provided
 *  with a raw diffraction pattern. It will subsequently remove the noise and background
 *  signal, and then detect and locate peaks. 
 * 
 * If the spacings are random (as far as the code is concerned), it will assume that
 *  the peak intensities have already been measured.
 * 
 * In either case, this object will contain the angles and integrated intensities of 
 *  each provided diffraction peak at the end of of the operation 
 * 
 * @param twoTheta [in] List of angles at which diffracted intensity is measured
 * @param intensity [in] Intensity measured at each angle
 */
void Diffraction::set(const Linked<double>& twoTheta, const Linked<double>& intensity) {
	// Copy over two theta
	vector<double> twoThetaCopy(twoTheta.length());
	Linked<double>::iterator iter = twoTheta.begin();
	for (int i = 0; i < twoTheta.length(); i++, iter++) twoThetaCopy[i] = *iter;

	// Copy over intensity
	vector<double> intensityCopy(intensity.length());
	iter = intensity.begin();
	for (int i = 0; i < intensity.length(); i++, iter++) intensityCopy[i] = *iter;

	// Call the real set function
	set(twoThetaCopy, intensityCopy);
}

/**
 * Given a list of angles and intensities, determine locations and intensities of 
 *  diffraction peaks. Store peak information internally.
 * 
 * If the intensities are uniformly spaced, this function assumes that it was provided
 *  with a raw diffraction pattern. It will subsequently remove the noise and background
 *  signal, and then detect and locate peaks. 
 * 
 * If the spacings are random (as far as the code is concerned), it will assume that
 *  the peak intensities have already been measured.
 * 
 * In either case, this object will contain the angles and integrated intensities of 
 *  each provided diffraction peak at the end of of the operation 
 * 
 * @param twoTheta [in, out] List of angles at which diffracted intensity is measured. Returned in sorted order
 * @param intensity [in] Intensity measured at each angle. Sorted in same order as twoTheta
 */
void Diffraction::set(vector<double>& twoTheta, vector<double>& intensity) {
    // Clear space
    clear();

    // Ensure that arrays are sorted in the ascending order.
    //  LW 8Apr14: Requires copying to pair vector and back, ick! It might be reasonable to 
    //             treat twoTheta and Intensity as a pair always. Consider for later
    vector<pair<double,double> > twoThetaIntensityPairs(twoTheta.size());
    for (int i=0; i<twoTheta.size(); i++) {
	pair<double,double> newPoint(twoTheta[i], intensity[i]);
	twoThetaIntensityPairs[i] = newPoint;
    }
    sort(twoThetaIntensityPairs.begin(), twoThetaIntensityPairs.end());
    for (int i=0; i<twoTheta.size(); i++) {
         twoTheta[i] = twoThetaIntensityPairs[i].first;
         intensity[i] = twoThetaIntensityPairs[i].second;
    }
    twoThetaIntensityPairs.clear();


	// Loop over two theta values and get max and min distances between them
	double curDif;
	double minDif = 0;
	double maxDif = 0;
	if (twoTheta.size() >= 2)
		minDif = maxDif = twoTheta[1] - twoTheta[0];
	for (int i = 1; i < twoTheta.size(); i++) {
		curDif = twoTheta[i] - twoTheta[i - 1];
		if (curDif < minDif)
			minDif = curDif;
		else if (curDif > maxDif)
			maxDif = curDif;
	}

	// Save peaks if data is already processed
	if ((maxDif > 1.1 * minDif) || (maxDif == 0)) {
		// Talk about what we are doing here
		Output::newline();
		Output::print("Importing an already-processed pattern");

		// Mark the type of pattern 
		_type = PT_EXP_INT;
		_diffractionPeaks.reserve(twoTheta.size());

		for (int i = 0; i < twoTheta.size(); ++i) {
			Diffraction::Peak newPeak(twoTheta[i], intensity[i]);
			_diffractionPeaks.push_back(newPeak);
		}

		// Ensure peaks are in sorted order
		std::sort(_diffractionPeaks.begin(), _diffractionPeaks.end());

		// Save min and max two theta
		_minTwoTheta = _diffractionPeaks[0].getAngle() - _resolution;
		_maxTwoTheta = _diffractionPeaks.back().getAngle() + _resolution / 2;
	} else { // Save peaks after processing

		// Output
		Output::newline();
		Output::print("Processing raw diffraction pattern");
		Output::increase();

		// Store raw pattern
		_continuousTwoTheta.resize(twoTheta.size());
		_continuousIntensity.resize(twoTheta.size());
		for (int i = 0; i < _continuousTwoTheta.size(); i++) {
			_continuousTwoTheta[i] = twoTheta[i];
			_continuousIntensity[i] = intensity[i];
		}

		// Make a copy of the data to process
		vector<double> twoThetaCopy(twoTheta);
		vector<double> intensityCopy(intensity);
		_minTwoTheta = *std::min_element(twoTheta.begin(), twoTheta.end());
		_maxTwoTheta = *std::max_element(twoTheta.begin(), twoTheta.end());

		// Prepare data for peak fitting
		smoothData(twoThetaCopy, intensityCopy);

		if (LW_EXCESSIVE_PRINTING == 1)
			savePattern("xray-smoothed.out", twoThetaCopy, intensityCopy);
		removeBackground(twoThetaCopy, intensityCopy);

		// Get peaks
		vector<vector<double> > peakTwoTheta;
		vector<vector<double> > peakIntensity;
		locatePeaks(peakTwoTheta, peakIntensity, twoThetaCopy, intensityCopy);
		getPeakIntensities(peakTwoTheta, peakIntensity);

		// Output
		Output::decrease();
	}

	extractIntensities();
}

/**
 * Given a diffraction pattern, report its match to the pattern stored in this object.
 * @param reference [in] Pattern to match against
 * @return Match between patterns, expressed as an R factor
 */
double Diffraction::rFactor(const Diffraction& reference) {
    // Output
    Output::newline();
    Output::print("Calculating R factor compared to reference pattern");
    Output::increase();

    // Match peaks in this pattern to reference pattern
    defineReferencePattern(reference);

    // Get the R factor (refining nothing)
    double optMaxIntensity;
    _currentlyRefining.clear();
    double rFactor = getCurrentRFactor(DR_ABS);

    // Output
    Output::newline();
    Output::print("Optimal R factor of ");
    Output::print(rFactor);

    // Output
    Output::decrease();

    // Return result
    return rFactor;
}

/**
 * Used during structural refinement. Sets positions of atoms.
 *
 * @param symmetry [in,out] Symmetry object describing structure (contains atomic positions)
 * @param positions [in] Vector containing all relevant position data
 */
void Diffraction::setPositions(const Symmetry& symmetry, const Vector& positions) {
    int i, j, k;
            Vector3D newPos;
    for (i = 0; i < symmetry.orbits().length(); ++i) {
        for (j = 0; j < symmetry.orbits()[i].atoms().length(); ++j) {

            for (k = 0; k < 3; ++k)
                    newPos[k] = positions[3 * i + k];
                    newPos *= symmetry.orbits()[i].generators()[j].rotation();
                    newPos += symmetry.orbits()[i].generators()[j].translations()[0];
                    ISO::moveIntoCell(newPos);
                    symmetry.orbits()[i].atoms()[j]->fractional(newPos);
            }
    }
}

/**
 * Ensure that derivatives of a property wrt position obey symmetry of a crystal.
 * 
 * @param symmetry [in] Object describing the symmetry of the structure
 * @param derivs [in,out] Derivatives with respect to each position parameter. Will be adjusted
 */
void Diffraction::symDerivatives(const Symmetry& symmetry, Vector& derivs) {
    int i, j;
            Vector3D tempDeriv;
    for (i = 0; i < symmetry.orbits().length(); ++i) {
        for (j = 0; j < 3; ++j)
                tempDeriv[j] = derivs[3 * i + j];
                tempDeriv *= symmetry.orbits()[i].specialPositions()[0].rotation();
            if (symmetry.orbits()[i].anyAtomsFixed())
                    tempDeriv = 0.0;

                for (j = 0; j < 3; ++j)
                        derivs[3 * i + j] = tempDeriv[j];
                }
}

/**
 * Ensure that positions obey symmetry of a crystal.
 * 
 * @param symmetry [in] Object describing the positions
 * @param postions [in,out] Derivatives with respect to each position parameter. Will be adjusted
 */
void Diffraction::symPositions(const Symmetry& symmetry, Vector& position) {
    int i, j;
            Vector3D tempPos;
    for (i = 0; i < symmetry.orbits().length(); ++i) {
        for (j = 0; j < 3; ++j)
                tempPos[j] = position[3 * i + j];
                tempPos -= symmetry.orbits()[i].specialPositions()[0].translation();
                tempPos *= symmetry.orbits()[i].specialPositions()[0].rotation();
                tempPos += symmetry.orbits()[i].specialPositions()[0].translation();

            for (j = 0; j < 3; ++j)
                    position[3 * i + j] = tempPos[j];
            }
}

/**
 * Calculate the current R factor based on the intensities stored in _peaks. This operation
 *  also calculates derivatives of R factor with respect to every variable that can 
 *  be refined.
 * 
 * User Notices:
 * <ol>
 * <li>You should have already have defined the reference pattern and structure
 * <li>Current intensity values in peaks are used. Make sure you run calculatePeakIntensities
 * with your new parameters before running this operation.
 * </ol>
 * 
 * @param method [in] Method used to calculate R factor
 * @return R factor describing match between this pattern and reference
 */
double Diffraction::getCurrentRFactor(Rmethod method) {
    if (!referencePatternIsDefined()) {
        Output::newline(ERROR);
        Output::print("Internal Error: Reference pattern not yet defined.");
        return -1;
    }
    
    // ---> Part #1: Store the peak intensities for both the reference pattern and 
    //   peaks from the calculated pattern that match those peaks
    
    // Intensities of peaks in the reference pattern that are matched by peaks 
    //  in the calculated pattern.
    vector<double> referenceIntensity(_referencePattern->_diffractionPeaks.size(), 0.0);
    // For each peak in the reference, total intensity of all matching peaks from this pattern
    vector<double> matchedIntensity(referenceIntensity.size(), 0.0);
    // Intensity of peaks in this pattern that do not match a reference peak
    vector<double> unmatchedIntensity(_unmatchedPeaks.size());
    
    for (int i = 0; i < referenceIntensity.size(); i++) {
        referenceIntensity[i] = _referencePattern->_diffractionPeaks[i].peakIntensity;
        for (int j=0; j < _matchingPeaks[i].size(); j++)
            matchedIntensity[i] += _diffractionPeaks[_matchingPeaks[i][j]].peakIntensity;
    }
    
    for (int i=0; i<unmatchedIntensity.size(); i++)
        unmatchedIntensity[i] = _diffractionPeaks[_unmatchedPeaks[i]].peakIntensity;

    // ---> Part #2: Calculate the normalization factor, which is what the total
    //  error is divided by to make it an R factor. This is generally the sum of
    //  the integrated intensities from the reference pattern
    double norm;
    switch (method) {
        case DR_SQUARED:
            norm = 0;
            for (int i = 0; i < referenceIntensity.size(); i++)
                norm += referenceIntensity[i] * referenceIntensity[i];
            break;
        case DR_ABS:
            norm = std::accumulate(referenceIntensity.begin(), referenceIntensity.end(), 0.0);
            break;
        default:
            Output::newline(ERROR);
            Output::print("Internal Error: No method set to determine normalization factor with this Rmethod.");
            return -1;
    }

    // ---> Part #3: Determine optimal scale factor
    double optimalScale = 1.0;
    if (method == DR_SQUARED) {
        // Very simple for this case. Since R is quadratic wrt scale factor it is 
        //  possible to solve where the first derivative is equal to zero:
        //  s_optimal = sum[ I_ref * I_calc ] / sum[ I_calc * I_calc ]
        // Note: We do not need to worry about where peaks are not matched , because
        //       I_ref * I_calc == 0 for those cases
        optimalScale = std::inner_product(matchedIntensity.begin(), matchedIntensity.end(), \
                referenceIntensity.begin(), 0.0);
        double denom = 0;
        for (int i=0; i<_diffractionPeaks.size(); i++) 
            denom += _diffractionPeaks[i].peakIntensity * _diffractionPeaks[i].peakIntensity;
        optimalScale /= denom;
    } else if (method == DR_ABS) {
        // In this case, the minimum occurs when at least one calculated peak intensity
        //  is exactly equal to reference peak intensity. So, we are going to loop through
        //  each of these conditions
        double minimumError = 1e100;
        for (int i = 0; i < matchedIntensity.size(); i++) {
            if (matchedIntensity[i] == 0) continue;
                double curScale = referenceIntensity[i] / matchedIntensity[i];
                double curError = 0;
                for (int j = 0; j < matchedIntensity.size(); j++)
                    curError += abs(referenceIntensity[j] - curScale * matchedIntensity[j]);
                for (int j = 0; j < _unmatchedPeaks.size(); j++) 
                    curError += abs(curScale * _diffractionPeaks[_unmatchedPeaks[j]].peakIntensity);
                if (curError < minimumError) {
                    minimumError = curError; optimalScale = curScale; 
                }
        }
    } else {
        Output::newline(ERROR);
        Output::print("Internal Error: No method set to determine optimal scale with this Rmethod.");
        return -1;
    }
    
    // ---> Part #5: Calculate R factor
    double rFactor = 0;
    switch (method) {
        case DR_SQUARED:
            for (int i=0; i < matchedIntensity.size(); i++)
                rFactor += pow(referenceIntensity[i] - optimalScale * matchedIntensity[i],2);
            for (int i=0; i<unmatchedIntensity.size(); i++)
                rFactor += pow(optimalScale * unmatchedIntensity[i], 2);
            rFactor /= norm;
            break;
        case DR_ABS:
            for (int i=0; i < matchedIntensity.size(); i++)
                rFactor += abs(referenceIntensity[i] - optimalScale * matchedIntensity[i]);
            for (int i=0; i<unmatchedIntensity.size(); i++)
                rFactor += abs(optimalScale * unmatchedIntensity[i]);
            rFactor /= norm;
            break;
        default:
            Output::newline(ERROR);
            Output::print("Internal Error: No method set to determine R factor with this Rmethod.");
            return -1;
    }
    
    // ---> Part #5: Calculate derivative of R factor with respect to each refinement parameter
    // This is only necessary for calculated patterns
    // LW 8Apr14: MARKED FOR DELETION. I don't think we are going to use analytical derivatives, 
    //            numerical ones seem to work well enough. See discussion in runRefinement().
    if (_type == PT_CALCULATED) {
        // Zero all derivatives (assume they have already been allocated)
        _BFactorDerivs =  0.0;
        _PositionDerivs = 0.0;
        if (method == DR_SQUARED) {
            // Accrue total derivatives from each peaks contributing at this angle
            for (int p=0; p < _diffractionPeaks.size(); p++) {
                // Get the error for this peak
                int patternIndex = _diffractionPeaks[p].patternIndex;
                double error;
                if (patternIndex == -1) error = _diffractionPeaks[p].peakIntensity;
                else error = referenceIntensity[patternIndex] 
                        - optimalScale * matchedIntensity[patternIndex];
                // Get the contribution to the derivatives from this peak
                valarray<double> curBderivs = _diffractionPeaks[p].derivBfactors;
                curBderivs *= -2 * optimalScale * error / norm;
                valarray<double> curPosDerivs = _diffractionPeaks[p].derivPositions;
                curPosDerivs *= -2 * optimalScale * error / norm;
                // Add the contribution to the total derivatives
                _BFactorDerivs += curBderivs;
                _PositionDerivs += curBderivs;
            }
        } else if (method == DR_ABS) {
            for (int p=0; p < _diffractionPeaks.size(); p++) {
                // Get the error for this peak
                int patternIndex = _diffractionPeaks[p].patternIndex;
                double error;
                if (patternIndex == -1) error = _diffractionPeaks[p].peakIntensity;
                else error = referenceIntensity[patternIndex] 
                        - optimalScale * matchedIntensity[patternIndex];
                // Get the contribution to the derivatives from this peak
                valarray<double> curBderivs = _diffractionPeaks[p].derivBfactors;
                curBderivs *= optimalScale / norm;
                if (error < 0) curBderivs *= -1;
                valarray<double> curPosDerivs = _diffractionPeaks[p].derivPositions;
                curPosDerivs *= optimalScale / norm;
                if (error < 0) curPosDerivs *= -1;
                // Add the contribution to the total derivatives
                _BFactorDerivs += curBderivs;
                _PositionDerivs += curBderivs;
            }
        }
    }
    
    return rFactor;
}

/**
 * @param text [in] Text object containing contents of file
 * @return Whether file contains diffraction data
 */
bool Diffraction::isFormat(const Text& text) {
    int pairCount = 0;
            int lineCount = 0;
    for (int i = 0; i < text.length(); ++i) {
        if (!text[i].length())
            continue;
            if (Language::isComment(text[i][0]))
                continue;
                ++lineCount;
                if (text[i].length() < 2)
                    continue;
                    if ((Language::isNumber(text[i][0])) && (Language::isNumber(text[i][1])))
                            ++pairCount;
                    }
    if (!lineCount)
        return false;
        if ((double) pairCount / lineCount < 0.5)

            return false;
            return true;
        }

/**
 * Extract diffraction data from file, store in this object. File must be in the 
 *  following format:
 * <p>wavelength [wavelength]<br>
 * fwhm [fwhm]<br>
 * variance [variance]<br>
 * [twoTheta1] [intensity1]<br>
 * {...}<br>
 * 
 * Wavelength and other text words are optional. As long as the file contains more
 *  at least two columns of numerical data, this functional will read it as a diffraction pattern.
 * 
 * @param text [in] Text object containing contents of a file
 */
void Diffraction::set(const Text& text) {
    // Output
    Output::newline();
    Output::print("Reading diffraction data from file");
    Output::increase();

    // Clear space
    clear();

    // Loop over lines in file
    int i;
    vector<double> rawTwoTheta;
    vector<double> rawIntensity;

    rawTwoTheta.reserve(text.length());
    rawIntensity.reserve(text.length());
    for (i = 0; i < text.length(); ++i) {

        // Skip if line is blank
        if (!text[i].length())
            continue;

        // Skip if line is too short
        if (text[i].length() < 2)
            continue;

        // Found wavelength
        if (text[i][0].equal("wavelength", false, 4)) {
            if (Language::isNumber(text[i][1]))
                    _wavelength = atof(text[i][1].array());
            else {
                Output::newline(ERROR);
                        Output::print("Did not recognize wavelength value in diffraction file (");
                        Output::print(text[i][1]);
                        Output::print(")");
                        Output::quit();
            }
        }
            // Found fwhm
        else if (text[i][0].equal("fwhm", false, 4)) {
            if (Language::isNumber(text[i][1]))
                    fwhm(atof(text[i][1].array()));
            else {
                Output::newline(ERROR);
                        Output::print("Did not recognize FWHM value in diffraction file (");
                        Output::print(text[i][1]);
                        Output::print(")");
                        Output::quit();
            }
        }
            // Found variance
        else if (text[i][0].equal("variance", false, 3)) {
            if (Language::isNumber(text[i][1]))
                    variance(atof(text[i][1].array()));
            else {
                Output::newline(ERROR);
                        Output::print("Did not recognize variance value in diffraction file (");
                        Output::print(text[i][1]);
                        Output::print(")");
                        Output::quit();
            }
        }
            // Found a data line
        else if ((Language::isNumber(text[i][0])) && (Language::isNumber(text[i][1]))) {
            rawTwoTheta.push_back(atof(text[i][0].array()));
            rawIntensity.push_back(atof(text[i][1].array()));
        }
    }

    // Reduce data arrays to minimal size
    rawTwoTheta.resize(rawTwoTheta.size());
    rawIntensity.resize(rawIntensity.size());

    // Process data
    set(rawTwoTheta, rawIntensity);

    // Output
    Output::newline();
    Output::print("Found ");
    Output::print(_peakTwoTheta.size());
    Output::print(" peak");
    if (_peakTwoTheta.size() != 1)
            Output::print("s");
            Output::increase();
    for (i = 0; i < _peakTwoTheta.size(); ++i) {

        Output::newline();
        Output::print("Two-theta and intensity of ");
        Output::print(_peakTwoTheta[i]);
        Output::print(" ");
        Output::print(_integratedIntensity[i]);
    }
    Output::decrease();

    // Output
    Output::decrease();
}

/** 
 * Apply smoothing function to intensities. Works by taking a weighted average of the
 *  value of a certain point and its closest neighbors. 
 * 
 * The weight of the center point is always equal to 1. It is up to the user
 * to determine how many neighbors to include in the smoothing, and the weight 
 * of the point farthest from the center. The weights of the other points is
 * determined through linear interpolation.
 * 
 * @param rawTwoTheta [in] Angles at which intensity is recorded
 * @param rawIntensity [in,out] Intensity measured at each angle. This function
 *  should remove noise from this data.
 * @param numPerSide [in] Number of points to use on either side for smoothing.
 * @param power [in] How much weight to apply to farthest points. Should range between 0 and 1.
 */
void Diffraction::smoothData(const vector<double>& rawTwoTheta, vector<double>& rawIntensity,
        const int numPerSide, const double power) {

    // Weight for point at max distance away
    double farWeight = power;

            // Calculate weights (linear scaling)
            int numSmoothPoints = numPerSide * 2 + 1;
            double weight[numSmoothPoints];
            weight[numPerSide] = 1.0;
            double totalWeight = 1.0;
    for (int i = 1; i <= numPerSide; i++) {
        double temp = 1.0 + (farWeight - 1.0) * (double) i / (double) numPerSide;
                totalWeight += 2 * temp;
                weight[numPerSide - i ] = temp;
                weight[numPerSide + i ] = temp;
    }
    for (int i = 0; i < numSmoothPoints; i++)
            weight[i] /= totalWeight;

            // Original intensity (before filtering)
            vector<double> initialInt(rawIntensity);

            // Calculate new intensities
        for (int i = numPerSide; i < rawIntensity.size() - numPerSide; i++) {
            double newValue = 0.0;
                    int startPoint = i - numPerSide;

            for (int j = 0; j < numSmoothPoints; j++)
                    newValue += weight[j] * initialInt[startPoint + j];
                    rawIntensity[i] = newValue;
            }
}

/**
 * Given a raw powder diffraction pattern, remove the background signal
 *
 * @param rawTwoTheta [in] Angles at which diffracted intensity is measured
 * @param rawItensity [in,out] Intensity measured at each angle. Background will be 
 *                      removed from these measurements.
 */
void Diffraction::removeBackground(vector<double>& rawTwoTheta, vector<double>& rawIntensity) {
    // Determine how many points to include in smoothing
    double boxSize = 4.0;
            int nPoints = (int) (boxSize / (rawTwoTheta[1] - rawTwoTheta[0]));
            int pointsPerSide = nPoints / 2;

            // Weight points based on the squared inverse of their intensities
            vector<double> fitWeight(rawIntensity.size());
    for (int i = 0; i < fitWeight.size(); i++) {
        fitWeight[i] = rawIntensity[i] > 0 ? 1.0 / rawIntensity[i] : 10;
                fitWeight[i] *= fitWeight[i]; fitWeight[i] *= fitWeight[i];
    }


    // Determine background signal as weighted average of background 
    vector<double> backgroundSignal(rawIntensity.size(), 0.0);
    for (int point = 0; point < backgroundSignal.size(); point++) {
        double totalWeight = 0.0;
                // Determine how many points to use in average
                int toAverage = min(point, pointsPerSide);
                toAverage = min(toAverage, (int) backgroundSignal.size() - 1 - point);
        for (int neigh = -toAverage; neigh <= toAverage; neigh++) {
            backgroundSignal[point] += fitWeight[point + neigh] * rawIntensity[point + neigh];
                    totalWeight += fitWeight[point + neigh];
        }
        backgroundSignal[point] /= totalWeight;
    }

    // Subtract background
    for (int i = 0; i < rawIntensity.size(); i++)
            rawIntensity[i] -= backgroundSignal[i];


    // If desired, print out background-less signal
    if (LW_EXCESSIVE_PRINTING == 1)
            savePattern("xray-nobackground.out", rawTwoTheta, rawIntensity, backgroundSignal);
}



/**
 * Detects diffraction peaks from raw diffraction pattern data.
 * 
 * @param peakTwoTheta [out] 2D list. For each detected peak, contains list of angles corresponding to measured peak
 * @param peakIntensity [out] 2D list. For each detected peak, contains intensities measured at each angle
 * @param rawTwoTheta [in] Raw diffraction pattern: Angles at which intensities were measured
 * @param rawIntensity [in] Raw diffraction pattern: Intensities measured at each angle
 */
void Diffraction::locatePeaks(vector<vector<double> >& peakTwoTheta,
        vector<vector<double> >& peakIntensity,
        const vector<double>& rawTwoTheta, const vector<double>& rawIntensity) {

    // Tolerance for point being on a peak. Fraction of maximum intensity
    double peakTol = 0.01;

    // Clear space for results
    peakTwoTheta.clear();
    peakIntensity.clear();

    // Loop over points to get maximum
    double maxHeight = *std::max_element(rawIntensity.begin(), rawIntensity.end());

    // Set absolute peak tolerance
    peakTol *= maxHeight;

    // First derivative of intensity wrt twoTheta
    vector<double> firstDerivative = getFirstDerivative(rawTwoTheta, rawIntensity);
    smoothData(rawTwoTheta, firstDerivative, 3, 1.0);
    if (LW_EXCESSIVE_PRINTING == 1)
    savePattern("xray-firstDerivative.out", rawTwoTheta, firstDerivative);

    // Second derivative of intensity wrt twoTheta
    vector<double> secondDerivative = getSecondDerivative(rawTwoTheta, rawIntensity);
    smoothData(rawTwoTheta, secondDerivative, 3, 1.0);
    if (LW_EXCESSIVE_PRINTING == 1)
        savePattern("xray-secondDerivative.out", rawTwoTheta, secondDerivative);
            
    // Position of center of peaks
    int position = 0;
    // Peak positions (stored as index in each array)
    vector<int> peakPosition;
    // Loop through entire pattern
    while (position < rawTwoTheta.size()) {
        // Loop until we find intensities that are above peak tolerance and a positive 
        //  second derivative
        while (rawIntensity[position] < peakTol || secondDerivative[position] < 0) {
            position++; if (position == rawTwoTheta.size()) break;
            }
        // See if we are done searching
        if (position == rawTwoTheta.size()) break;

                // Loop until second derivatives crosses zero
            while (secondDerivative[position] > 0) {
                position++; if (position == rawTwoTheta.size()) break;
                }
        if (position == rawTwoTheta.size()) break;

                // Loop until first derivative crosses zero (marks center)
            while (firstDerivative[position] > 0) {
                position++; if (position == rawTwoTheta.size()) break;
                }
        if (position == rawTwoTheta.size()) break;
                peakPosition.push_back(position);

                // Loop until end of peak (second derivative goes positive)
            while (secondDerivative[position] < 0) {
                position++; if (position == rawTwoTheta.size()) break;
                }
        if (position == rawTwoTheta.size()) { // Peak did not finish
            peakPosition.pop_back(); break;
        }
    }

    // Part 2: Store edges of peaks
    peakIntensity.clear();
    peakIntensity.reserve(peakPosition.size());
    peakTwoTheta.clear();
    peakTwoTheta.reserve(peakPosition.size());

    // Lower edge of peak
    int leftMinimum = 0;
    double temp = rawIntensity[0];
    for (int i = 1; i < peakPosition[0]; i++) {
        if (rawIntensity[i] < temp) {
            temp = rawIntensity[i]; leftMinimum = i;
        }
    }
    // Upper edge of peak
    int rightMinimum = 0;
    for (int i = 0; i < peakPosition.size(); i++) {
        // Step 1: Find the minimum between next peak and current peak
        temp = rawIntensity[peakPosition[i]];
                rightMinimum = peakPosition[i];
                double rightMaximum = i == peakPosition.size() - 1 ? rawIntensity.size()
                : peakPosition[i + 1];
        for (position = peakPosition[i]; position < rightMaximum; position++) {
            if (rawIntensity[position] < temp) {
                temp = rawIntensity[position]; rightMinimum = position; }
        }
        // Step 2: Define edges of peak (either the minimum between peaks, 
        //         or where intensity crosses zero). Store in a deque temporarily
        deque<double> dequePeakTwoTheta, dequePeakIntensity;
                position = peakPosition[i];
        while (position >= leftMinimum && rawIntensity[position] > 0) {
            dequePeakTwoTheta.push_front(rawTwoTheta[position]);
                    dequePeakIntensity.push_front(rawIntensity[position]);
                    position--;
        }
        position = peakPosition[i] + 1;
        while (position <= rightMinimum && rawIntensity[position] > 0) {
            dequePeakTwoTheta.push_back(rawTwoTheta[position]);
                    dequePeakIntensity.push_back(rawIntensity[position]);
                    position++;
        }
        // Step 3: Prepare to move to next peak
        // Store current peak
        vector<double> tempVector(dequePeakTwoTheta.begin(),
                dequePeakTwoTheta.end());
        if (tempVector.size() > 0) {
            peakTwoTheta.push_back(tempVector);
                    tempVector.assign(dequePeakIntensity.begin(), dequePeakIntensity.end());
                    peakIntensity.push_back(tempVector);
        }
        // Update boundary
        leftMinimum = rightMinimum;
    }

    // Part 3: Combine peaks that are smaller than certain criteria
    position = 0;
    while (position < peakTwoTheta.size()) {
        bool toRemove = false;
                // Step 1: Check maximum intensity of peak is above 2% of the maximum
                double peakHeight = *max_element(peakIntensity[position].begin(), \
                    peakIntensity[position].end());
                toRemove = peakHeight < 0.02 * maxHeight;
                // Step 2: If tall enough, check if peak width is greater than 0.1 degrees
                double peakWidth = peakTwoTheta[position].back() - peakTwoTheta[position][0];
        if (!toRemove) toRemove = peakWidth < 0.05;
                // Step 3: If it fails either test, remove this peak.
            if (toRemove) {
                // First, check if peak on right is connected
                if (position != peakTwoTheta.size() - 1 && \
                        peakTwoTheta[position].back() == peakTwoTheta[position + 1][0]) {
                    // If so, add this peak to that one
                    peakTwoTheta[position + 1].insert(peakTwoTheta[position + 1].begin(), \
                            peakTwoTheta[position].begin(), peakTwoTheta[position].end());
                            peakIntensity[position + 1].insert(peakIntensity[position + 1].begin(), \
                            peakIntensity[position].begin(), peakIntensity[position].end());
                } else if (position != 0 && \
                        peakTwoTheta[position][0] == peakTwoTheta[position - 1].back()) {
                    peakTwoTheta[position - 1].insert(peakTwoTheta[position - 1].end(), \
                            peakTwoTheta[position].begin(), peakTwoTheta[position].end());
                            peakIntensity[position - 1].insert(peakIntensity[position - 1].end(), \
                            peakIntensity[position].begin(), peakIntensity[position].end());
                }
                // Now delete, the peak
                peakTwoTheta.erase(peakTwoTheta.begin() + position);
                        peakIntensity.erase(peakIntensity.begin() + position);
            }

              else position++;
            }
    Output::decrease();
}

/** 
 * Given a list of peaks extracted from a raw diffraction pattern, fit interpolation functions to 
 *  each peak and use them to calculate area under peak. Locations and integrated intensities
 *  of each peak are stored internally in _twoTheta and _intensity.
 * 
 * The two inputs to this function are created by running a background-subtracted, 
 *  smoothed raw diffraction pattern through getPeaks().
 * 
 * @param peakTwoTheta [in] 2D array containing angles corresponding to each peak. Peaks must 
 *  be arranged in ascending order with two theta.
 * @param peakIntensity [in] 2D array containing intensities at each measured angle for each peak
 */
void Diffraction::getPeakIntensities(const vector<vector<double> >& peakTwoTheta, \
        const vector<vector<double> >& peakIntensity) {
    // Clear the old peak data out
    _diffractionPeaks.clear();
    _diffractionPeaks.reserve(peakIntensity.size());
    
    // Functors used in fitting
    Functor<Diffraction> gaussFun(this, &Diffraction::gaussian);
    Functor<Diffraction> compositeGaussFun(this, &Diffraction::compositeGaussian);
    Functor<Diffraction> psFun(this, &Diffraction::PV);
    Functor<Diffraction> compositePVFun(this, &Diffraction::compositePV);
    VectorFunctor<Diffraction> gaussDeriv(this, &Diffraction::gaussianDerivs);
    VectorFunctor<Diffraction> compositeGaussDeriv(this, &Diffraction::compositeGaussianDerivs);
    VectorFunctor<Diffraction> psDeriv(this, &Diffraction::PVderivs);
    VectorFunctor<Diffraction> compositePVDeriv(this, &Diffraction::compositePVDerivs);


    // Part #0: Save two-theta/intensity pairs in a format compatable with 
    //  Mint's fitting programs
    List<double>::D3 singlePeakPoints(peakTwoTheta.size());
    for (int curPeak = 0; curPeak < peakTwoTheta.size(); curPeak++) {
        singlePeakPoints[curPeak].length(peakTwoTheta[curPeak].size());
                vector<double> curPeakTwoTheta = peakTwoTheta[curPeak];
                vector<double> curPeakIntensity = peakIntensity[curPeak];
        for (int j = 0; j < peakTwoTheta[curPeak].size(); ++j) {
            singlePeakPoints[curPeak][j].length(2);
                    singlePeakPoints[curPeak][j][0] = curPeakTwoTheta[j];
                    singlePeakPoints[curPeak][j][1] = curPeakIntensity[j];
        }
    }

    // Part #1: Fit individual peaks with Gaussian functions
    OList<Vector> gaussianParams(peakTwoTheta.size());
    for (int curPeak = 0; curPeak < peakTwoTheta.size(); ++curPeak) {
        // Initialize fitting parameters for Gaussian function
        Vector initialGaussian(3);
                initialGaussian[0] = 0.25;
                initialGaussian[1] = singlePeakPoints[curPeak][0][0];
                initialGaussian[2] = singlePeakPoints[curPeak][0][1];
        for (int j = 1; j < singlePeakPoints[curPeak].length(); ++j) {
            if (singlePeakPoints[curPeak][j][1] > initialGaussian[2]) {
                initialGaussian[1] = singlePeakPoints[curPeak][j][0];
                        initialGaussian[2] = singlePeakPoints[curPeak][j][1];
            }
        }

        // Fit Gaussian function
        gaussianParams[curPeak] = Fit::LM<Diffraction>(singlePeakPoints[curPeak], gaussFun, gaussDeriv, initialGaussian, 1e-5);
    }

    // Part #2: Decide which peaks need to be fitted together
    // Detect which peaks are in contact
    vector<vector<int> > peakGroup; peakGroup.reserve(peakTwoTheta.size());
    {
        vector<int> newVec; newVec.push_back(0); peakGroup.push_back(newVec); }
    for (int peak = 1; peak < peakTwoTheta.size(); peak++) {
        double peakStart = peakTwoTheta[peak].front();
                double lastGroupEnd = peakTwoTheta[peakGroup.back().back()].back();
        if (peakStart - lastGroupEnd < 0.1)
                peakGroup.back().push_back(peak);
        else {
            vector<int> newGroup(1, peak);
            peakGroup.push_back(newGroup);
        }
    }
    // Combine data from peak groups
	List<double>::D3 peakGroupPoints(peakGroup.size());
	for (int group = 0; group < peakGroup.size(); group++) {
		int totalPeakSize = 0;
		for (int subPeak = 0; subPeak < peakGroup[group].size(); subPeak++)
			totalPeakSize += singlePeakPoints[peakGroup[group][subPeak]].length();
		// Copying points from each singlePeak (in reverse order)
		peakGroupPoints[group].length(totalPeakSize);
		for (int subPeak = peakGroup[group].size() - 1; subPeak >= 0; subPeak--) {
			int curPeak = peakGroup[group][subPeak];
			for (int point = singlePeakPoints[curPeak].length() - 1; point >= 0; point--)
				peakGroupPoints[group][--totalPeakSize] = singlePeakPoints[curPeak][point];
		}
	}

    // Part #3: Fit groups of peaks with multiple Gaussian functions. Note
    //  that each Gaussian still corresponds to a single peak.
	for (int group = 0; group < peakGroup.size(); group++) {
		// If only a single peak in group, do nothing
		if (peakGroup[group].size() == 1) continue;
		// Step #1: Extract parameters from individual Gaussians
		Vector initialCompositeParams(3 * peakGroup[group].size());
		for (int peak = 0; peak < peakGroup[group].size(); peak++)
			for (int i = 0; i < 3; i++)
				initialCompositeParams[peak * 3 + i] = gaussianParams[peakGroup[group][peak]][i];
		// Step #2: Fit new parameters 
		Vector compositeParams = Fit::LM<Diffraction>(peakGroupPoints[group],
				compositeGaussFun, compositeGaussDeriv, initialCompositeParams, 1e-5);
		// Step #3: Copy parameters back
		for (int peak = 0; peak < peakGroup[group].size(); peak++)
			for (int i = 0; i < 3; i++)
				gaussianParams[peakGroup[group][peak]][i] = compositeParams[peak * 3 + i];
	}

    // Part #4: Fit groups of peaks with multiple pseudo-Voight functions
    OList<Vector> psParams(peakTwoTheta.size());
            // Step #1: Convert Guassian results into inital pseudo-Voight guesses
    for (int curPeak = 0; curPeak < peakTwoTheta.size(); curPeak++) {
        Vector initialPS(8);
        // Initialize pseudo-voigt fitting parameters using solution from Gaussian
        initialPS[0] = 1.0; // Purely weight on the Guassian
        initialPS[1] = initialPS[2] = 0.0;
        initialPS[3] = gaussianParams[curPeak][1];
        // double temp = tan(Num<double>::toRadians(initialPS[3]/2));
        // initialPS[4] = initialPS[5] = initialPS[6] = gaussianParams[curPeak][0] / (1 + temp + temp*temp);
        initialPS[4] = gaussianParams[curPeak][0];
        initialPS[5] = initialPS[6] = 0.0;
        initialPS[7] = gaussianParams[curPeak][2];

        // Store result
        psParams[curPeak] = initialPS;
	}

	for (int group = 0; group < peakGroup.size(); group++) {
		// Step #2: Extract parameters from individual peaks
		Vector initialCompositeParams(8 * peakGroup[group].size());
		for (int peak = 0; peak < peakGroup[group].size(); peak++)
			for (int i = 0; i < 8; i++)
				initialCompositeParams[peak * 8 + i] = psParams[peakGroup[group][peak]][i];
		// Step #3: Fit new parameters 
		Vector compositeParams = Fit::LM<Diffraction>(peakGroupPoints[group],
				compositePVFun, compositePVDeriv, initialCompositeParams, 1e-5);
		// Step #4: Copy parameters back
		for (int peak = 0; peak < peakGroup[group].size(); peak++)
			for (int i = 0; i < 8; i++)
				psParams[peakGroup[group][peak]][i] = compositeParams[peak * 8 + i];
	}

    // Part #5: Extract peak intensities, store in _diffractionPeaks
    Output::increase();
    for (int group = 0; group < peakGroup.size(); group++) {
        double groupMin = peakTwoTheta[peakGroup[group].front()].front();
        double groupMax = peakTwoTheta[peakGroup[group].back()].back();
        for (int subPeak = 0; subPeak < peakGroup[group].size(); subPeak++) {
            int curPeak = peakGroup[group][subPeak];
            double initialTwoTheta = psParams[curPeak][3];
            double twoThetaStep = 1e-3;
            Functor<Diffraction> psTT(this, &Diffraction::PV, &psParams[curPeak]);
            double location, intensity;
            
            // Find the maximum (location of peak)
            intensity = Solve<Diffraction>::maximize(psTT, 1e-8, initialTwoTheta, twoThetaStep, location);
            
            // Integrate the peak
            PVPeakFunction peakFunc(this, psParams[curPeak]);
            intensity = dlib::integrate_function_adapt_simp(peakFunc, groupMin, groupMax, 1e-8);
			
			// Check that results make sense
			if (intensity < 0.0) {
				Output::newline(WARNING);
				Output::print("Failure during peak integration - Negative intensity found near: ");
				Output::print(location, 3);
				Output::quit();
			}
			
			// Check that the maximum is within bounds of the measurement
			if (location < this->_minTwoTheta || location > this->_maxTwoTheta) {
				Output::newline(WARNING);
				Output::print("Failure during peak integration - Peak maximum outside of measured range: ");
				Output::print(location, 3);
				Output::quit();
			}
            
            // Make a new peak
            Peak newPeak(location, intensity);
            _diffractionPeaks.push_back(newPeak);
        }
    }
    Output::decrease();

    // Now that we are done, print out all the data
    if (LW_EXCESSIVE_PRINTING == 1) {
        for (int peak = 0; peak < peakTwoTheta.size(); peak++) {
            vector<double> fittedIntensity(peakTwoTheta[peak].size());
            for (int pos = 0; pos < peakTwoTheta[peak].size(); pos++)
                    fittedIntensity[pos] = PV(psParams[peak], peakTwoTheta[peak][pos]);
                    //fittedIntensity[pos] = gaussian(gaussianParams[peak], peakTwoTheta[peak][pos]);
                    Word filename = "peaks/peak";
                    filename += Language::numberToWord(peak) + Word(".out");
                    savePattern(filename, peakTwoTheta[peak], peakIntensity[peak], fittedIntensity);
            }
    }
}



/**
 * Print diffraction data stored in this object to a file.
 * 
 * @param file [in] Name of file to receive output. Can be "stdout"
 * @param broaden [in] Whether to print the diffraction pattern as a continuous function (true),
 *    or only the peak centers and integrated intensities (false).
 */
void Diffraction::print(const Word& file, bool broaden) const {

    // Open file for writing if needed
    int origStream = Output::streamID();
            PrintMethod origMethod = Output::method();
    if (file != "stdout")
            Output::setStream(Output::addStream(file));

    // Set output method
    Output::method(STANDARD);

    // If printing to file, then print settings
    Output message;
    if (file != "stdout") {
        Output::newline();
        Output::print("Wavelength ");
        Output::print(_wavelength);
        Output::newline();
        Output::print("FWHM ");
        Output::printSci(_fwhm);
        Output::newline();
        Output::print("Resolution ");
        Output::print(_resolution);
    }
    // If printing to screen then add header
    else {
        message.addLine();
        message.add("Two-theta");
        message.add("Intensity");
        message.addLine();
        message.add("---------");
        message.add("---------");
    }

    // Add peaks
    if (!broaden) {
        message.addLines(_peakTwoTheta.size());
        for (int i = 0; i < _peakTwoTheta.size(); ++i) {
            if (_integratedIntensity[i] < 1) 
                continue;
            message.addLine();
            message.add(_peakTwoTheta[i], 10);
            message.add(_integratedIntensity[i], 10);
        }
    }
    // Add broadened data
    else {
        int i;
        double intensity;
        double twoTheta = _minTwoTheta - 3 * _fwhm;
        message.addLines((int) (_maxTwoTheta - _minTwoTheta) / _resolution + 1);
        for (twoTheta = (twoTheta < 3) ? 3 : twoTheta; twoTheta <= _maxTwoTheta + 3 * _fwhm; twoTheta += _resolution) {
            intensity = 0;
            for (i = 0; i < _peakTwoTheta.size(); ++i)
                intensity += _integratedIntensity[i] * exp(-pow(twoTheta - _peakTwoTheta[i], 2.0) / (2 * _variance));
            message.addLine();
            message.add(twoTheta, 10);
            message.add(intensity, 10);
        }
    }

    // Print peaks
    Output::newline();
    Output::print(message, RIGHT);

    // Reset output
    if (file != "stdout")
        Output::removeStream(Output::streamID());
    Output::setStream(origStream);
    Output::method(origMethod);
}

/**
 * Internally store parameters used for atomic form factor (used to calculated atomic scattering
 * factor) . This function simply gets the atomic form factors for each symmetrically
 * unique set of atoms in the structure (should be already defined)
 */
void Diffraction::setATFParams() {
    if (!structureIsDefined()) {
        Output::newline(ERROR);
                Output::print("Structure has not yet been defined. Cannot get ATF parameters");
    }
    // Allocate space
    _atfParams.length(_symmetry->orbits().length());

            // Get parameters
            double a1, a2, a3, a4, b1, b2, b3, b4, c;
    for (int i = 0; i < _symmetry->orbits().length(); ++i) {
        const Element& element = _symmetry->orbits()[i].atoms()[0]->element();
        if (element.number() == 1) {
            a1 = 0.489918; b1 = 20.659300; a2 = 0.262003; b2 = 7.740390; a3 = 0.196767; b3 = 49.551899;
                    a4 = 0.049879; b4 = 2.201590; c = 0.001305;
        } else if (element.number() == 2) {
            a1 = 0.873400; b1 = 9.103700; a2 = 0.630900; b2 = 3.356800; a3 = 0.311200; b3 = 22.927601;
                    a4 = 0.178000; b4 = 0.982100; c = 0.006400;
        } else if (element.number() == 3) {
            a1 = 1.128200; b1 = 3.954600; a2 = 0.750800; b2 = 1.052400; a3 = 0.617500; b3 = 85.390503;
                    a4 = 0.465300; b4 = 168.261002; c = 0.037700;
        } else if (element.number() == 4) {
            a1 = 1.591900; b1 = 43.642700; a2 = 1.127800; b2 = 1.862300; a3 = 0.539100; b3 = 103.483002;
                    a4 = 0.702900; b4 = 0.542000; c = 0.038500;
        } else if (element.number() == 5) {
            a1 = 2.054500; b1 = 23.218500; a2 = 1.332600; b2 = 1.021000; a3 = 1.097900; b3 = 60.349800;
                    a4 = 0.706800; b4 = 0.140300; c = -0.193200;
        } else if (element.number() == 6) {
            a1 = 2.310000; b1 = 20.843901; a2 = 1.020000; b2 = 10.207500; a3 = 1.588600; b3 = 0.568700;
                    a4 = 0.865000; b4 = 51.651199; c = 0.215600;
        } else if (element.number() == 7) {
            a1 = 12.212600; b1 = 0.005700; a2 = 3.132200; b2 = 9.893300; a3 = 2.012500; b3 = 28.997499;
                    a4 = 1.166300; b4 = 0.582600; c = -11.529000;
        } else if (element.number() == 8) {
            a1 = 3.048500; b1 = 13.277100; a2 = 2.286800; b2 = 5.701100; a3 = 1.546300; b3 = 0.323900;
                    a4 = 0.867000; b4 = 32.908901; c = 0.250800;
        } else if (element.number() == 9) {
            a1 = 3.539200; b1 = 10.282500; a2 = 2.641200; b2 = 4.294400; a3 = 1.517000; b3 = 0.261500;
                    a4 = 1.024300; b4 = 26.147600; c = 0.277600;
        } else if (element.number() == 10) {
            a1 = 3.955300; b1 = 8.404200; a2 = 3.112500; b2 = 3.426200; a3 = 1.454600; b3 = 0.230600;
                    a4 = 1.125100; b4 = 21.718399; c = 0.351500;
        } else if (element.number() == 11) {
            a1 = 4.762600; b1 = 3.285000; a2 = 3.173600; b2 = 8.842200; a3 = 1.267400; b3 = 0.313600;
                    a4 = 1.112800; b4 = 129.423996; c = 0.676000;
        } else if (element.number() == 12) {
            a1 = 5.420400; b1 = 2.827500; a2 = 2.173500; b2 = 79.261101; a3 = 1.226900; b3 = 0.380800;
                    a4 = 2.307300; b4 = 7.193700; c = 0.858400;
        } else if (element.number() == 13) {
            a1 = 6.420200; b1 = 3.038700; a2 = 1.900200; b2 = 0.742600; a3 = 1.593600; b3 = 31.547199;
                    a4 = 1.964600; b4 = 85.088600; c = 1.115100;
        } else if (element.number() == 14) {
            a1 = 6.291500; b1 = 2.438600; a2 = 3.035300; b2 = 32.333698; a3 = 1.989100; b3 = 0.678500;
                    a4 = 1.541000; b4 = 81.693703; c = 1.140700;
        } else if (element.number() == 15) {
            a1 = 6.434500; b1 = 1.906700; a2 = 4.179100; b2 = 27.157000; a3 = 1.780000; b3 = 0.526000;
                    a4 = 1.490800; b4 = 68.164497; c = 1.114900;
        } else if (element.number() == 16) {
            a1 = 6.905300; b1 = 1.467900; a2 = 5.203400; b2 = 22.215099; a3 = 1.437900; b3 = 0.253600;
                    a4 = 1.586300; b4 = 56.172001; c = 0.866900;
        } else if (element.number() == 17) {
            a1 = 11.460400; b1 = 0.010400; a2 = 7.196400; b2 = 1.166200; a3 = 6.255600; b3 = 18.519400;
                    a4 = 1.645500; b4 = 47.778400; c = -9.557400;
        } else if (element.number() == 18) {
            a1 = 7.484500; b1 = 0.907200; a2 = 6.772300; b2 = 14.840700; a3 = 0.653900; b3 = 43.898300;
                    a4 = 1.644200; b4 = 33.392899; c = 1.444500;
        } else if (element.number() == 19) {
            a1 = 8.218600; b1 = 12.794900; a2 = 7.439800; b2 = 0.774800; a3 = 1.051900; b3 = 213.186996;
                    a4 = 0.865900; b4 = 41.684101; c = 1.422800;
        } else if (element.number() == 20) {
            a1 = 8.626600; b1 = 10.442100; a2 = 7.387300; b2 = 0.659900; a3 = 1.589900; b3 = 85.748398;
                    a4 = 1.021100; b4 = 178.436996; c = 1.375100;
        } else if (element.number() == 21) {
            a1 = 9.189000; b1 = 9.021300; a2 = 7.367900; b2 = 0.572900; a3 = 1.640900; b3 = 136.108002;
                    a4 = 1.468000; b4 = 51.353100; c = 1.332900;
        } else if (element.number() == 22) {
            a1 = 9.759500; b1 = 7.850800; a2 = 7.355800; b2 = 0.500000; a3 = 1.699100; b3 = 35.633801;
                    a4 = 1.902100; b4 = 116.105003; c = 1.280700;
        } else if (element.number() == 23) {
            a1 = 10.297100; b1 = 6.865700; a2 = 7.351100; b2 = 0.438500; a3 = 2.070300; b3 = 26.893801;
                    a4 = 2.057100; b4 = 102.477997; c = 1.219900;
        } else if (element.number() == 24) {
            a1 = 10.640600; b1 = 6.103800; a2 = 7.353700; b2 = 0.392000; a3 = 3.324000; b3 = 20.262600;
                    a4 = 1.492200; b4 = 98.739899; c = 1.183200;
        } else if (element.number() == 25) {
            a1 = 11.281900; b1 = 5.340900; a2 = 7.357300; b2 = 0.343200; a3 = 3.019300; b3 = 17.867399;
                    a4 = 2.244100; b4 = 83.754303; c = 1.089600;
        } else if (element.number() == 26) {
            a1 = 11.769500; b1 = 4.761100; a2 = 7.357300; b2 = 0.307200; a3 = 3.522200; b3 = 15.353500;
                    a4 = 2.304500; b4 = 76.880501; c = 1.036900;
        } else if (element.number() == 27) {
            a1 = 12.284100; b1 = 4.279100; a2 = 7.340900; b2 = 0.278400; a3 = 4.003400; b3 = 13.535900;
                    a4 = 2.348800; b4 = 71.169197; c = 1.011800;
        } else if (element.number() == 28) {
            a1 = 12.837600; b1 = 3.878500; a2 = 7.292000; b2 = 0.256500; a3 = 4.443800; b3 = 12.176300;
                    a4 = 2.380000; b4 = 66.342102; c = 1.034100;
        } else if (element.number() == 29) {
            a1 = 13.338000; b1 = 3.582800; a2 = 7.167600; b2 = 0.247000; a3 = 5.615800; b3 = 11.396600;
                    a4 = 1.673500; b4 = 64.812599; c = 1.191000;
        } else if (element.number() == 30) {
            a1 = 14.074300; b1 = 3.265500; a2 = 7.031800; b2 = 0.233300; a3 = 5.165200; b3 = 10.316300;
                    a4 = 2.410000; b4 = 58.709702; c = 1.304100;
        } else if (element.number() == 31) {
            a1 = 15.235400; b1 = 3.066900; a2 = 6.700600; b2 = 0.241200; a3 = 4.359100; b3 = 10.780500;
                    a4 = 2.962300; b4 = 61.413502; c = 1.718900;
        } else if (element.number() == 32) {
            a1 = 16.081600; b1 = 2.850900; a2 = 6.374700; b2 = 0.251600; a3 = 3.706800; b3 = 11.446800;
                    a4 = 3.683000; b4 = 54.762501; c = 2.131300;
        } else if (element.number() == 33) {
            a1 = 10.672300; b1 = 2.634500; a2 = 6.070100; b2 = 0.264700; a3 = 3.431300; b3 = 12.947900;
                    a4 = 4.277900; b4 = 47.797199; c = 2.531000;
        } else if (element.number() == 34) {
            a1 = 17.000601; b1 = 2.409800; a2 = 5.819600; b2 = 0.272600; a3 = 3.973100; b3 = 15.237200;
                    a4 = 4.354300; b4 = 43.816299; c = 2.840900;
        } else if (element.number() == 35) {
            a1 = 17.178900; b1 = 2.172300; a2 = 5.235800; b2 = 16.579599; a3 = 5.637700; b3 = 0.260900;
                    a4 = 3.985100; b4 = 41.432800; c = 2.955700;
        } else if (element.number() == 36) {
            a1 = 17.355499; b1 = 1.938400; a2 = 6.728600; b2 = 16.562300; a3 = 5.549300; b3 = 0.226100;
                    a4 = 3.537500; b4 = 39.397202; c = 2.825000;
        } else if (element.number() == 37) {
            a1 = 17.178400; b1 = 1.788800; a2 = 9.643500; b2 = 17.315100; a3 = 5.139900; b3 = 0.274800;
                    a4 = 1.529200; b4 = 164.934006; c = 3.487300;
        } else if (element.number() == 38) {
            a1 = 17.566299; b1 = 1.556400; a2 = 9.818400; b2 = 14.098800; a3 = 5.422000; b3 = 0.166400;
                    a4 = 2.669400; b4 = 132.376007; c = 2.506400;
        } else if (element.number() == 39) {
            a1 = 17.775999; b1 = 1.402900; a2 = 10.294600; b2 = 12.800600; a3 = 5.726290; b3 = 0.125599;
                    a4 = 3.265880; b4 = 104.353996; c = 1.912130;
        } else if (element.number() == 40) {
            a1 = 17.876499; b1 = 1.276180; a2 = 10.948000; b2 = 11.916000; a3 = 5.417320; b3 = 0.117622;
                    a4 = 3.657210; b4 = 87.662697; c = 2.069290;
        } else if (element.number() == 41) {
            a1 = 17.614201; b1 = 1.188650; a2 = 12.014400; b2 = 11.766000; a3 = 4.041830; b3 = 0.204785;
                    a4 = 3.533460; b4 = 69.795700; c = 3.755910;
        } else if (element.number() == 42) {
            a1 = 3.702500; b1 = 0.277200; a2 = 17.235600; b2 = 1.095800; a3 = 12.887600; b3 = 11.004000;
                    a4 = 3.742900; b4 = 61.658401; c = 4.387500;
        } else if (element.number() == 43) {
            a1 = 19.130100; b1 = 0.864132; a2 = 11.094800; b2 = 8.144870; a3 = 4.649010; b3 = 21.570700;
                    a4 = 2.712630; b4 = 86.847198; c = 5.404280;
        } else if (element.number() == 44) {
            a1 = 19.267401; b1 = 0.808520; a2 = 12.918200; b2 = 8.434670; a3 = 4.863370; b3 = 24.799700;
                    a4 = 1.567560; b4 = 94.292801; c = 5.378740;
        } else if (element.number() == 45) {
            a1 = 19.295700; b1 = 0.751536; a2 = 14.350100; b2 = 8.217580; a3 = 4.734250; b3 = 25.874901;
                    a4 = 1.289180; b4 = 98.606201; c = 5.328000;
        } else if (element.number() == 46) {
            a1 = 19.331900; b1 = 0.698655; a2 = 15.501700; b2 = 7.989290; a3 = 5.295370; b3 = 25.205200;
                    a4 = 0.605844; b4 = 76.898598; c = 5.265930;
        } else if (element.number() == 47) {
            a1 = 19.280800; b1 = 0.644600; a2 = 16.688499; b2 = 7.472600; a3 = 4.804500; b3 = 24.660500;
                    a4 = 1.046300; b4 = 99.815598; c = 5.179000;
        } else if (element.number() == 48) {
            a1 = 19.221399; b1 = 0.594600; a2 = 17.644400; b2 = 6.908900; a3 = 4.461000; b3 = 24.700800;
                    a4 = 1.602900; b4 = 87.482498; c = 5.069400;
        } else if (element.number() == 49) {
            a1 = 19.162399; b1 = 0.547600; a2 = 18.559601; b2 = 6.377600; a3 = 4.294800; b3 = 25.849899;
                    a4 = 2.039600; b4 = 92.802902; c = 4.939100;
        } else if (element.number() == 50) {
            a1 = 19.188900; b1 = 5.830300; a2 = 19.100500; b2 = 0.503100; a3 = 4.458500; b3 = 26.890900;
                    a4 = 2.466300; b4 = 83.957100; c = 4.782100;
        } else if (element.number() == 51) {
            a1 = 19.641800; b1 = 5.303400; a2 = 19.045500; b2 = 0.460700; a3 = 5.037100; b3 = 27.907400;
                    a4 = 2.682700; b4 = 75.282501; c = 4.590900;
        } else if (element.number() == 52) {
            a1 = 19.964399; b1 = 4.817420; a2 = 19.013800; b2 = 0.420885; a3 = 6.144870; b3 = 28.528400;
                    a4 = 2.523900; b4 = 70.840302; c = 4.352000;
        } else if (element.number() == 53) {
            a1 = 20.147200; b1 = 4.347000; a2 = 18.994900; b2 = 0.381400; a3 = 7.513800; b3 = 27.766001;
                    a4 = 2.273500; b4 = 66.877602; c = 4.071200;
        } else if (element.number() == 54) {
            a1 = 20.293301; b1 = 3.928200; a2 = 19.029800; b2 = 0.344000; a3 = 8.976700; b3 = 26.465900;
                    a4 = 1.990000; b4 = 64.265800; c = 3.711800;
        } else if (element.number() == 55) {
            a1 = 20.389200; b1 = 3.569000; a2 = 19.106199; b2 = 0.310700; a3 = 10.662000; b3 = 24.387899;
                    a4 = 1.495300; b4 = 213.904007; c = 3.335200;
        } else if (element.number() == 56) {
            a1 = 20.336100; b1 = 3.216000; a2 = 19.297001; b2 = 0.275600; a3 = 10.888000; b3 = 20.207300;
                    a4 = 2.695900; b4 = 167.201996; c = 2.773100;
        } else if (element.number() == 57) {
            a1 = 20.577999; b1 = 2.948170; a2 = 19.599001; b2 = 0.244475; a3 = 11.372700; b3 = 18.772600;
                    a4 = 3.287190; b4 = 133.123993; c = 2.146780;
        } else if (element.number() == 58) {
            a1 = 21.167101; b1 = 2.812190; a2 = 19.769501; b2 = 0.226836; a3 = 11.851300; b3 = 17.608299;
                    a4 = 3.330490; b4 = 127.112999; c = 1.862640;
        } else if (element.number() == 59) {
            a1 = 22.044001; b1 = 2.773930; a2 = 19.669701; b2 = 0.222087; a3 = 12.385600; b3 = 16.766899;
                    a4 = 2.824280; b4 = 143.643997; c = 2.058300;
        } else if (element.number() == 60) {
            a1 = 22.684500; b1 = 2.662480; a2 = 19.684700; b2 = 0.210628; a3 = 12.774000; b3 = 15.885000;
                    a4 = 2.851370; b4 = 137.903000; c = 1.984860;
        } else if (element.number() == 61) {
            a1 = 23.340500; b1 = 2.562700; a2 = 19.609501; b2 = 0.202088; a3 = 13.123500; b3 = 15.100900;
                    a4 = 2.875160; b4 = 132.720993; c = 2.028760;
        } else if (element.number() == 62) {
            a1 = 24.004200; b1 = 2.472740; a2 = 19.425800; b2 = 0.196451; a3 = 13.439600; b3 = 14.399600;
                    a4 = 2.896040; b4 = 128.007004; c = 2.209630;
        } else if (element.number() == 63) {
            a1 = 24.627399; b1 = 2.387900; a2 = 19.088600; b2 = 0.194200; a3 = 13.760300; b3 = 13.754600;
                    a4 = 2.922700; b4 = 123.174004; c = 2.574500;
        } else if (element.number() == 64) {
            a1 = 25.070900; b1 = 2.253410; a2 = 19.079800; b2 = 0.181951; a3 = 13.851800; b3 = 12.933100;
                    a4 = 3.545450; b4 = 101.398003; c = 2.419600;
        } else if (element.number() == 65) {
            a1 = 25.897600; b1 = 2.242560; a2 = 18.218500; b2 = 0.196143; a3 = 14.316700; b3 = 12.664800;
                    a4 = 2.953540; b4 = 115.362000; c = 3.589240;
        } else if (element.number() == 66) {
            a1 = 26.507000; b1 = 2.180200; a2 = 17.638300; b2 = 0.202172; a3 = 14.559600; b3 = 12.189900;
                    a4 = 2.965770; b4 = 111.874001; c = 4.297280;
        } else if (element.number() == 67) {
            a1 = 26.904900; b1 = 2.070510; a2 = 17.294001; b2 = 0.197940; a3 = 14.558300; b3 = 11.440700;
                    a4 = 3.638370; b4 = 92.656601; c = 4.567960;
        } else if (element.number() == 68) {
            a1 = 27.656300; b1 = 2.073560; a2 = 16.428499; b2 = 0.223545; a3 = 14.977900; b3 = 11.360400;
                    a4 = 2.982330; b4 = 105.703003; c = 5.920460;
        } else if (element.number() == 69) {
            a1 = 28.181900; b1 = 2.028590; a2 = 15.885100; b2 = 0.238849; a3 = 15.154200; b3 = 10.997500;
                    a4 = 2.987060; b4 = 102.960999; c = 6.756210;
        } else if (element.number() == 70) {
            a1 = 28.664101; b1 = 1.988900; a2 = 15.434500; b2 = 0.257119; a3 = 15.308700; b3 = 10.664700;
                    a4 = 2.989630; b4 = 100.417000; c = 7.566720;
        } else if (element.number() == 71) {
            a1 = 28.947599; b1 = 1.901820; a2 = 15.220800; b2 = 9.985190; a3 = 15.100000; b3 = 0.261033;
                    a4 = 3.716010; b4 = 84.329803; c = 7.976280;
        } else if (element.number() == 72) {
            a1 = 29.143999; b1 = 1.832620; a2 = 15.172600; b2 = 9.599900; a3 = 14.758600; b3 = 0.275116;
                    a4 = 4.300130; b4 = 72.028999; c = 8.581540;
        } else if (element.number() == 73) {
            a1 = 29.202400; b1 = 1.773330; a2 = 15.229300; b2 = 9.370460; a3 = 14.513500; b3 = 0.295977;
                    a4 = 4.764920; b4 = 63.364399; c = 9.243540;
        } else if (element.number() == 74) {
            a1 = 29.081800; b1 = 1.720290; a2 = 15.430000; b2 = 9.225900; a3 = 14.432700; b3 = 0.321703;
                    a4 = 5.119820; b4 = 57.056000; c = 9.887500;
        } else if (element.number() == 75) {
            a1 = 28.762100; b1 = 1.671910; a2 = 15.718900; b2 = 9.092270; a3 = 14.556400; b3 = 0.350500;
                    a4 = 5.441740; b4 = 52.086102; c = 10.472000;
        } else if (element.number() == 76) {
            a1 = 28.189400; b1 = 1.629030; a2 = 16.155001; b2 = 8.979480; a3 = 14.930500; b3 = 0.382661;
                    a4 = 5.675890; b4 = 48.164700; c = 11.000500;
        } else if (element.number() == 77) {
            a1 = 27.304899; b1 = 1.592790; a2 = 16.729601; b2 = 8.865530; a3 = 15.611500; b3 = 0.417916;
                    a4 = 5.833770; b4 = 45.001099; c = 11.472200;
        } else if (element.number() == 78) {
            a1 = 27.005899; b1 = 1.512930; a2 = 17.763901; b2 = 8.811740; a3 = 15.713100; b3 = 0.424593;
                    a4 = 5.783700; b4 = 38.610298; c = 11.688300;
        } else if (element.number() == 79) {
            a1 = 16.881901; b1 = 0.461100; a2 = 18.591299; b2 = 8.621600; a3 = 25.558201; b3 = 1.482600;
                    a4 = 5.860000; b4 = 36.395599; c = 12.065800;
        } else if (element.number() == 80) {
            a1 = 20.680901; b1 = 0.545000; a2 = 19.041700; b2 = 8.448400; a3 = 21.657499; b3 = 1.572900;
                    a4 = 5.967600; b4 = 38.324600; c = 12.608900;
        } else if (element.number() == 81) {
            a1 = 27.544600; b1 = 0.655150; a2 = 19.158400; b2 = 8.707510; a3 = 15.538000; b3 = 1.963470;
                    a4 = 5.525930; b4 = 45.814899; c = 13.174600;
        } else if (element.number() == 82) {
            a1 = 31.061701; b1 = 0.690200; a2 = 13.063700; b2 = 2.357600; a3 = 18.441999; b3 = 8.618000;
                    a4 = 5.969600; b4 = 47.257900; c = 13.411800;
        } else if (element.number() == 83) {
            a1 = 33.368900; b1 = 0.704000; a2 = 12.951000; b2 = 2.923800; a3 = 16.587700; b3 = 8.793700;
                    a4 = 6.469200; b4 = 48.009300; c = 13.578200;
        } else if (element.number() == 84) {
            a1 = 34.672600; b1 = 0.700999; a2 = 15.473300; b2 = 3.550780; a3 = 13.113800; b3 = 9.556420;
                    a4 = 7.025800; b4 = 47.004501; c = 13.677000;
        } else if (element.number() == 85) {
            a1 = 35.316299; b1 = 0.685870; a2 = 19.021099; b2 = 3.974580; a3 = 9.498870; b3 = 11.382400;
                    a4 = 7.425180; b4 = 45.471500; c = 13.710800;
        } else if (element.number() == 86) {
            a1 = 35.563099; b1 = 0.663100; a2 = 21.281601; b2 = 4.069100; a3 = 8.003700; b3 = 14.042200;
                    a4 = 7.443300; b4 = 44.247299; c = 13.690500;
        } else if (element.number() == 87) {
            a1 = 35.929901; b1 = 0.646453; a2 = 23.054701; b2 = 4.176190; a3 = 12.143900; b3 = 23.105200;
                    a4 = 2.112530; b4 = 150.645004; c = 13.724700;
        } else if (element.number() == 88) {
            a1 = 35.763000; b1 = 0.616341; a2 = 22.906401; b2 = 3.871350; a3 = 12.473900; b3 = 19.988701;
                    a4 = 3.210970; b4 = 142.324997; c = 13.621100;
        } else if (element.number() == 89) {
            a1 = 35.659698; b1 = 0.589092; a2 = 23.103201; b2 = 3.651550; a3 = 12.597700; b3 = 18.599001;
                    a4 = 4.086550; b4 = 117.019997; c = 13.526600;
        } else if (element.number() == 90) {
            a1 = 35.564499; b1 = 0.563359; a2 = 23.421900; b2 = 3.462040; a3 = 12.747300; b3 = 17.830900;
                    a4 = 4.807030; b4 = 99.172203; c = 13.431400;
        } else if (element.number() == 91) {
            a1 = 35.884701; b1 = 0.547751; a2 = 23.294800; b2 = 3.415190; a3 = 14.189100; b3 = 16.923500;
                    a4 = 4.172870; b4 = 105.250999; c = 13.428700;
        } else if (element.number() == 92) {
            a1 = 36.022800; b1 = 0.529300; a2 = 23.412800; b2 = 3.325300; a3 = 14.949100; b3 = 16.092699;
                    a4 = 4.188000; b4 = 100.612999; c = 13.396600;
        } else if (element.number() == 93) {
            a1 = 36.187401; b1 = 0.511929; a2 = 23.596399; b2 = 3.253960; a3 = 15.640200; b3 = 15.362200;
                    a4 = 4.185500; b4 = 97.490799; c = 13.357300;
        } else if (element.number() == 94) {
            a1 = 36.525398; b1 = 0.499384; a2 = 23.808300; b2 = 3.263710; a3 = 16.770700; b3 = 14.945500;
                    a4 = 3.479470; b4 = 105.980003; c = 13.381200;
        } else if (element.number() == 95) {
            a1 = 36.670601; b1 = 0.483629; a2 = 24.099199; b2 = 3.206470; a3 = 17.341499; b3 = 14.313600;
                    a4 = 3.493310; b4 = 102.273003; c = 13.359200;
        } else if (element.number() == 96) {
            a1 = 36.648800; b1 = 0.465154; a2 = 24.409599; b2 = 3.089970; a3 = 17.399000; b3 = 13.434600;
                    a4 = 4.216650; b4 = 88.483398; c = 13.288700;
        } else if (element.number() == 97) {
            a1 = 36.788101; b1 = 0.451018; a2 = 24.773600; b2 = 3.046190; a3 = 17.891899; b3 = 12.894600;
                    a4 = 4.232840; b4 = 86.002998; c = 13.275400;
        } else if (element.number() == 98) {
            a1 = 36.918499; b1 = 0.437533; a2 = 25.199499; b2 = 3.007750; a3 = 18.331699; b3 = 12.404400;
                    a4 = 4.243910; b4 = 83.788101; c = 13.267400;
        } else {

            Output::newline(ERROR);
                    Output::print("Atomic scattering factor is not defined for ");
                    Output::print(element.symbol());
                    Output::quit();
        }
        _atfParams[i].length(9);
                _atfParams[i][0] = a1;
                _atfParams[i][1] = a2;
                _atfParams[i][2] = a3;
                _atfParams[i][3] = a4;
                _atfParams[i][4] = b1;
                _atfParams[i][5] = b2;
                _atfParams[i][6] = b3;
                _atfParams[i][7] = b4;
                _atfParams[i][8] = c;
    }
}

/**
 * Calculate atomic scattering factor for a group of symmetrically-identical atoms
 * @param index [in] Which group of atoms (same order as in symmetry.orbits())
 * @param angle [in] Angle radiation is scattered into
 * @return Atomic scattering factor for those atoms
 */
double Diffraction::Peak::atomicScatteringFactor(int index, double angle) {
    // Get s value
    double s = sin(angle) / _sourcePattern->_wavelength;
            double s2 = s * s;
    if (s > 2) {

        Output::newline(WARNING);
                Output::print("Atomic scattering factor is not optimized for s greater than 2");
    }
            

    // Return result
    return _sourcePattern->_atfParams[index][0] * exp(-_sourcePattern->_atfParams[index][4] * s2) + \
           _sourcePattern->_atfParams[index][1] * exp(-_sourcePattern->_atfParams[index][5] * s2) + \
           _sourcePattern->_atfParams[index][2] * exp(-_sourcePattern->_atfParams[index][6] * s2) + \
           _sourcePattern->_atfParams[index][3] * exp(-_sourcePattern->_atfParams[index][7] * s2) + _sourcePattern->_atfParams[index][8];
}

/**
 * Prints out a diffraction pattern so that it can be visualized in another program.
 * 
 * Note: twoTheta, Intensity, and fittedIntensity must all be the same length.
 * 
 * @param filename [in] Name of file in which to save raw pattern
 * @param twoTheta [in] Angles at which intensity is measured
 * @param Intensity [in] Intensity measured at each angle
 * @param fittedIntensity [in] Optional: Some other function that should be expressed as a function of angle
 */
void Diffraction::savePattern(const Word& filename, const vector<double>& twoTheta, \
        const vector<double>& Intensity, const vector<double>& otherIntensity) {
    int origStream = Output::streamID();
            Output::setStream(Output::addStream(filename));
            Output::newline();
    for (int i = 0; i < twoTheta.size(); i++) {
        Output::printPadded(twoTheta[i], 10, RIGHT, 3);
                Output::printPaddedSci(Intensity[i], 15, RIGHT, 5);

        if (otherIntensity.size() != 0)
                Output::printPaddedSci(otherIntensity[i], 15, RIGHT, 5);
                Output::newline();
        }
    Output::removeStream(Output::streamID());
            Output::setStream(origStream);
}

/**
 * Calculate the first derivative of a curve. x and y must be the same length. 
 * 
 * <p>For now, assumes that the spacing of points along x is uniform.
 * @param x Array of values of independent variable
 * @param y f(x)
 * @return f'(x) for each point. Same length as x and y
 */
vector<double> Diffraction::getFirstDerivative(const vector<double>& x, const vector<double>& y) {
    // Initialize output
    vector<double> d(x.size(), 0.0);
            // Calculate derivatives
            double h = 2 * (x[1] - x[0]);
            d[0] = (d[1] - d[0]) / h * 2;
    for (int i = 1; i < x.size() - 1; i++)
            d[i] = (y[i + 1] - y[i - 1]) / h;
            d.back() = (d.back() - d[d.size() - 2]) / h * 2;

    return d;
}

/**
 * Calculate the second derivative of a curve. x and y must be the same length. 
 * <p>For now, assumes that the spacing of points along x is uniform.
 * @param x Array of values of independent variable
 * @param y f(x)
 * @return f''(x) for each point. Same length as x and y
 */
vector<double> Diffraction::getSecondDerivative(const vector<double>& x, const vector<double>& y) {
    // Initialize output
    vector<double> d(x.size(), 0.0);
            // Calculate derivatives
            double h2 = (x[1] - x[0]); h2 *= h2;
    for (int i = 1; i < x.size() - 1; i++)
            d[i] = (y[i + 1] - 2 * y[i] + y[i - 1]) / h2;
            // Just store the 2nd and 2nd to last values as the first and last, respectively
            d[0] = d[1];
            d.back() = d[d.size() - 2];
    return d;
}
