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
#include "iso.h"
#include "output.h"



/* Basis& Basis::operator= (const Basis& rhs)
 *
 * Copy Basis object
 */

Basis& Basis::operator= (const Basis& rhs)
{
	if (this != &rhs)
	{
		_volume = rhs._volume;
		_vectors = rhs._vectors;
		_inverse = rhs._inverse;
		_metric = rhs._metric;
		_vectorsTranspose = rhs._vectorsTranspose;
		_inverseTranspose = rhs._inverseTranspose;
		_lengths = rhs._lengths;
		_angles = rhs._angles;
		_reduced = rhs._reduced;
		_reducedTranspose = rhs._reducedTranspose;
		_reducedInverse = rhs._reducedInverse;
		_reducedMetric = rhs._reducedMetric;
		_unitToReduced = rhs._unitToReduced;
		_reducedPointToUnit = rhs._reducedPointToUnit;
		_unitPointToReduced = rhs._unitPointToReduced;
		_latticeSystem = rhs._latticeSystem;
		for (int i = 0; i < 3; ++i)
		{
			_lengthFixed[i] = rhs._lengthFixed[i];
			_angleFixed[i] = rhs._angleFixed[i];
		}
	}
	return *this;
}



/* void Basis::set(const Matrix3D& vectors, bool showOutput)
 *
 * Set basis by vectors
 */

void Basis::set(const Matrix3D& vectors, bool showOutput)
{
	
	// Save vectors
	_vectors = vectors;
	
	// Set lengths and angles
	_lengths = lengths(vectors);
	_angles = angles(vectors);
	
	// Finish setup
	finishSetup(showOutput);
}



/* void Basis::set(const Vector3D& lengths, const Vector3D& angles, bool showOutput)
 *
 * Set basis by lengths and angles
 */

void Basis::set(const Vector3D& lengths, const Vector3D& angles, bool showOutput)
{
	
	// Save lengths and angles
	_lengths = lengths;
	_angles = angles;
	
	// Set vectors
	_vectors = vectors(lengths, angles);
	
	// Finish setup
	finishSetup(showOutput);
}



/* double Basis::getAngle(const Vector3D& vec1, const Vector3D& vec2)
 *
 * Get angle between two vectors
 */

double Basis::getAngle(const Vector3D& vec1, const Vector3D& vec2)
{
	Vector3D cross = vec1.cross(vec2);
	double denom = vec1 * vec2;
	if ((denom >= 0) && (denom < 1e-8))
		denom = 1e-8;
	else if ((denom <= 0) && (denom > -1e-8))
		denom = -1e-8;
	return atan(cross.magnitude() / denom);
}



/* Vector3D Basis::rotationAxis(const Matrix3D& matrix)
 *
 * Return the axis about which the rotation occurs
 */

Vector3D Basis::rotationAxis(const Matrix3D& matrix)
{
	
	// Subtract identity from matrix 
	Matrix3D tempMat = matrix;
	if (tempMat.determinant() < 0)
		tempMat *= -1;
	tempMat -= Matrix3D::identity();
	
	// Return rotation axis
	return backSolve(tempMat);
}



/* Vector3D Basis::backSolve(const Matrix3D& rotation, const Vector3D* fill)
 *
 * Use back substitution to solve for eigenvector
 */

Vector3D Basis::backSolve(const Matrix3D& rotation, const Vector3D* fill)
{
	
	// Turn matrix into row echelon form
	Matrix3D newMat = rotation.rowEchelon(0, 0, true);
	
	// Get the column of the pivot on each row if one occurs
	int i, j;
	int pivCol[3] = {-1, -1, -1};
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			if (Num<double>::abs(newMat(i, j)) > 1e-4)
			{
				pivCol[i] = j;
				break;
			}
		}
	}
	
	// Arrange rows so that pivots occur along diagonal
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			if (pivCol[j] == i)
			{
				newMat.swapRows(i, j);
				swap(pivCol[i], pivCol[j]);
				break;
			}
		}
	}
	
	// Vector to store solution
	Vector3D res;
	
	// Solve for answers
	int curFill = 0;
	double total;
	for (i = 2; i >= 0; --i)
	{
		if (Num<double>::abs(newMat(i, i)) > 1e-4)
		{
			total = 0;
			for (j = i + 1; j < 3; ++j)
				total += newMat(i, j) * res[j];
			res[i] = -total / newMat(i, i);
		}
		else if (!fill)
			res[i] = 1;
		else
		{
			res[i] = (*fill)[curFill];
			curFill++;
		}
	}
	
	// Convert all values to integers
	int temp[2];
	double denom[3];
	for (i = 0; i < 3; i++)
	{
		Num<double>::decimalTofraction(temp, res[i], 1e-4);
		denom[i] = (double) temp[1];
	}
	
	// Convert all values to integers
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
			res[i] *= denom[j];
	}
	
	// Divide values by the greatest common factor
	int intRes[3] = {(int) Num<double>::round(res[0], 1), (int) Num<double>::round(res[1], 1), \
		(int) Num<double>::round(res[2], 1)};
	int gcf = Num<int>::gcf(3, intRes);
	for (i = 0; i < 3; ++i)
		res[i] /= gcf;
	
	// Convert to standard orientation
	for (i = 2; i >= 0; i--)
	{
		if (Num<double>::abs(res[i]) > 1e-4)
		{
			if (res[i] < 0)
			{
				for (j = 0; j < 3; j++)
					res[j] *= -1;
			}
			break;
		}
	}
	
	// Return result
	return res;
}



/* Matrix3D Basis::vectors(const Vector3D& lengths, const Vector3D& angles)
 *
 * Return the matrix of basis vectors given a set of lengths and angles
 */

Matrix3D Basis::vectors(const Vector3D& lengths, const Vector3D& angles)
{
	Matrix3D res;
	res(0, 0) = lengths[0];
	res(0, 1) = 0;
	res(0, 2) = 0;
	res(1, 0) = lengths[1]*cos(angles[2]);
	res(1, 1) = lengths[1]*sin(angles[2]);
	res(1, 2) = 0;
	res(2, 0) = lengths[2]*cos(angles[1]);
	res(2, 1) = (lengths[1]*lengths[2]*cos(angles[0]) - res(1, 0)*res(2, 0)) / res(1, 1);
	res(2, 2) = sqrt(lengths[2]*lengths[2] - res(2, 0)*res(2, 0) - res(2, 1)*res(2, 1));
	return res;
}



/* Vector3D Basis::lengths(const Matrix3D& vectors)
 * 
 * Return lengths of basis vectors
 */

Vector3D Basis::lengths(const Matrix3D& vectors)
{
	Vector3D res;
	res[0] = sqrt(vectors(0, 0)*vectors(0, 0) + vectors(0, 1)*vectors(0, 1) + vectors(0, 2)*vectors(0, 2));
	res[1] = sqrt(vectors(1, 0)*vectors(1, 0) + vectors(1, 1)*vectors(1, 1) + vectors(1, 2)*vectors(1, 2));
	res[2] = sqrt(vectors(2, 0)*vectors(2, 0) + vectors(2, 1)*vectors(2, 1) + vectors(2, 2)*vectors(2, 2));
	return res;
}



/* Vector3D Basis::angles(const Matrix3D& vectors)
 *
 * Return angles between basis vectors
 */

Vector3D Basis::angles(const Matrix3D& vectors)
{
	Vector3D vec0(vectors(0, 0), vectors(0, 1), vectors(0, 2));
	Vector3D vec1(vectors(1, 0), vectors(1, 1), vectors(1, 2));
	Vector3D vec2(vectors(2, 0), vectors(2, 1), vectors(2, 2));
	return Vector3D(vec1.angle(vec2), vec0.angle(vec2), vec0.angle(vec1));
}



/* void Basis::finishSetup(bool showOutput)
 *
 * Set basis properties
 */

void Basis::finishSetup(bool showOutput)
{
	
	// Turn off output if needed
	if (!showOutput)
		Output::quietOn();
	
	// Output
	Output::newline();
	Output::print("Setting the basis");
	Output::increase();
	
	// Calculate volume
	_volume = _vectors.volume();
	
	// Get the shortest basis vector length
	double minLength = _lengths[0];
	if (_lengths[1] < minLength)
		minLength = _lengths[1];
	if (_lengths[2] < minLength)
		minLength = _lengths[2];
	
	// Volume is negative
	double tol = minLength > 0.01 ? minLength / 100 : 0.01;
	if (_volume < -tol)
	{
		Output::newline(ERROR);
		Output::print("Setting basis with negative volume");
		Output::quit();
	}
	
	// Volume is zero
	if (_volume < tol)
	{
		Output::newline(ERROR);
		Output::print("Setting basis with zero volume");
		Output::quit();
	}
	
	// Calculate complementary matrices
	_inverse = _vectors.inverse();
	_vectorsTranspose = _vectors.transpose();
	_inverseTranspose = _inverse.transpose();
	_metric = _vectors * _vectorsTranspose;
	
	// Calculate the reduced cell
	_reduced = reduced(this->vectors(), &_unitToReduced);
	_reducedTranspose = _reduced.transpose();
	_reducedInverse = _reduced.inverse();
	_reducedMetric = _reduced * _reducedTranspose;
	
	// Save transformations of positions between unit and reduced cell
	_reducedPointToUnit = _unitToReduced.transpose();
	_unitPointToReduced = _reducedPointToUnit.inverse();
	
	// Print volume
	Output::newline();
	Output::print("Volume: ");
	Output::print(_volume);
	Output::print(" Ang^3");
	
	// Print basis vectors
	int i, j;
	Output::newline();
	Output::print("Basis vectors");
	for (i = 0; i < 3; ++i)
	{
		Output::newline();
		Output::tab();
		for (j = 0; j < 3; ++j)
		{
			Output::print(_vectors(i, j), 8);
			Output::print(" ");
		}
	}
	
	// Print basis metrices
	Output::newline();
	Output::print("Basis metrics");
	Output::newline();
	Output::tab();
	Output::print("A: ");
	Output::print(_lengths[0], 8);
	Output::print(" Ang");
	Output::newline();
	Output::tab();
	Output::print("B: ");
	Output::print(_lengths[1], 8);
	Output::print(" Ang");
	Output::newline();
	Output::tab();
	Output::print("C: ");
	Output::print(_lengths[2], 8);
	Output::print(" Ang");
	Output::newline();
	Output::tab();
	Output::print("Alpha: ");
	Output::print(Num<double>::toDegrees(_angles[0]), 4);
	Output::print(" deg");
	Output::newline();
	Output::tab();
	Output::print("Beta: ");
	Output::print(Num<double>::toDegrees(_angles[1]), 4);
	Output::print(" deg");
	Output::newline();
	Output::tab();
	Output::print("Gamma: ");
	Output::print(Num<double>::toDegrees(_angles[2]), 4);
	Output::print(" deg");
	
	// Print inverse basis vectors
	Output::newline();
	Output::print("Inverse of basis vectors");
	for (i = 0; i < 3; ++i)
	{
		Output::newline();
		Output::tab();
		for (j = 0; j < 3; ++j)
		{
			Output::print(_inverse(i, j), 8);
			Output::print(" ");
		}
	}
	
	// Print reduced basis vectors
	Output::newline();
	Output::print("Reduced basis vectors");
	for (i = 0; i < 3; ++i)
	{
		Output::newline();
		Output::tab();
		for (j = 0; j < 3; ++j)
		{
			Output::print(_reduced(i, j), 8);
			Output::print(" ");
		}
	}
	
	// Output
	Output::decrease();
	
	// Reset output
	if (!showOutput)
		Output::quietOff();
}



/* double Basis::distance(const Vector3D& pos1, CoordinateType type1, const Vector3D& pos2, \
 *		CoordinateType type2, Vector3D* cell) const
 *
 * Return the distance between two positions
 */

double Basis::distance(const Vector3D& pos1, CoordinateType type1, const Vector3D& pos2, \
	CoordinateType type2, Vector3D* cell) const
{
	
	// Get fractional coordinates
	_frac1 = (type1 == FRACTIONAL) ? pos1 : getFractional(pos1);
	_frac2 = (type2 == FRACTIONAL) ? pos2 : getFractional(pos2);
	
    // Get the difference between points in fractional coordinates
	_redDif = _unitPointToReduced * (_frac2 - _frac1);

	// Set which cells should be searched
	for (_i = 0; _i < 3; ++_i)
	{
		_numCellsToSearch[_i] = 1;
		_cellsToSearch[_i][0] = -Num<double>::floor(_redDif[_i]);
		_curDis = _cellsToSearch[_i][0] + _redDif[_i];
		if (_curDis < -1e-4)
		{
			_numCellsToSearch[_i] = 2;
			_cellsToSearch[_i][1] = _cellsToSearch[_i][0] + 1;
		}
		else if (_curDis > 1e-4)
		{
			_numCellsToSearch[_i] = 2;
			_cellsToSearch[_i][1] = _cellsToSearch[_i][0] - 1;
		}
	}
	
	// Loop over cells to search to find minimum distance
	_first = true;
	for (_i = 0; _i < _numCellsToSearch[0]; ++_i)
	{
		for (_j = 0; _j < _numCellsToSearch[1]; ++_j)
		{
			for (_k = 0; _k < _numCellsToSearch[2]; ++_k)
			{
				_curVec[0]  = _redDif[0] + _cellsToSearch[0][_i];
				_curVec[1]  = _redDif[1] + _cellsToSearch[1][_j];
				_curVec[2]  = _redDif[2] + _cellsToSearch[2][_k];
				_curDis  = _curVec * (_reducedMetric * _curVec);
				
				// Found a new minimum
				if ((_first) || (_curDis < _minDis))
				{
					_first = false;
					_minDis = _curDis;
					_minCell[0] = _cellsToSearch[0][_i];
					_minCell[1] = _cellsToSearch[1][_j];
					_minCell[2] = _cellsToSearch[2][_k];
				}
			}
		}
	}
	
	// Save cell if needed
	if (cell)
		*cell = _reducedPointToUnit * _minCell;
	
    // Return distance
	return sqrt(_minDis);
}



/* double Basis::secondDistance(const Vector3D& pos1, CoordinateType type1, const Vector3D& pos2, \
 *		CoordinateType type2, Vector3D* cell) const
 *
 * Return the second nearest distance between two positions
 */

double Basis::secondDistance(const Vector3D& pos1, CoordinateType type1, const Vector3D& pos2, \
	CoordinateType type2, Vector3D* cell) const
{
	
	// Get fractional coordinates
	_frac1 = (type1 == FRACTIONAL) ? pos1 : getFractional(pos1);
	_frac2 = (type2 == FRACTIONAL) ? pos2 : getFractional(pos2);
    
    // Get the difference between points in fractional coordinates
	_frac2[0] -= _frac1[0];
	_frac2[1] -= _frac1[1];
	_frac2[2] -= _frac1[2];
	_redDif = _unitPointToReduced * _frac2;

	// Loop over cells to search
	_first = true;
	_second = false;
	for (_i = -1; _i <= 1; ++_i)
	{
		for (_j = -1; _j <= 1; ++_j)
		{
			for (_k = -1; _k <= 1; ++_k)
			{
				_curVec[0]  = _redDif[0] + _i;
				_curVec[1]  = _redDif[1] + _j;
				_curVec[2]  = _redDif[2] + _k;
				_curDis  = _curVec * (_reducedMetric * _curVec);
				
				// This is the first cell
				if (_first)
				{
					_minDis = _curDis;
					_minCell[0] = _i;
					_minCell[1] = _j;
					_minCell[2] = _k;
					_first = false;
					_second = true;
				}
				
				// Found a new minimum distance
				else if (_curDis < _minDis)
				{
					_secondDis = _minDis;
					_minDis = _curDis;
					_secondCell[0] = _minCell[0];
					_secondCell[1] = _minCell[1];
					_secondCell[2] = _minCell[2];
					_minCell[0] = _i;
					_minCell[1] = _j;
					_minCell[2] = _k;
				}
				
				// Found a new second distance
				else if ((_second) || (_curDis < _secondDis))
				{
					_secondDis = _curDis;
					_secondCell[0] = _i;
					_secondCell[1] = _j;
					_secondCell[2] = _k;
				}
				
				// Save that no longer on second cell searched
				if (_second)
					_second = false;
			}
		}
	}
	
	// Save cell if needed
	if (cell)
		*cell = _reducedPointToUnit * _secondCell;

	// Return distance
	return sqrt(_secondDis);
}



/* Matrix3D Basis::reduced(const Matrix3D& vectors, Matrix3D* transformation)
 *
 * Get the reduced basis vectors
 */

Matrix3D Basis::reduced(const Matrix3D& vectors, Matrix3D* transformation)
{
	
	// Get the transformation to reduced basis
	Matrix3D trans = reducedTransformation(vectors);

	// Save the transformation if needed
	if (transformation)
		*transformation = trans;

    // Set the new basis
	return trans * vectors;
}



/* Matrix3D Basis::reducedTransformation(const Matrix3D& vectors)
 *
 * Get the transformation matrix from unit to reduced cell
 * Acta. Cryst. (1976) A32, 297
 * Acta. Cryst. (2003) A60, 1
 */

Matrix3D Basis::reducedTransformation(const Matrix3D& vectors)
{
	
	// Set the initial parameters
	double A;
	double B;
	double C;
	double ksi;
	double eta;
	double zeta;
	NiggliParams(vectors, A, B, C, ksi, eta, zeta);

	// Minimum reduction variables
	bool firstMinRedTest = true;
	bool minRedOn = false;
	double minRedMult = 10;
	double lastA;
	double lastB;
	double lastC;

    // Set initial tolerance
	double eps = (1e-5) * pow(vectors.volume(), 1.0/3.0);

    // Initialize transformation matrix to identity matrix
	Matrix3D transformation;
	transformation.makeIdentity();
	
	// Variables to control reduction
	int maxLoops = 1000;
	
	// Loop until basis is in reduced form
	int count;
	int val[3];
	int plusMin[3];
	for (count = 1; count <= maxLoops; count++)
	{

		// Test 1
        if ((Num<double>::gt(A, B, eps)) || 
			((Num<double>::eq(A, B, eps)) && (Num<double>::gt(Num<double>::abs(ksi), Num<double>::abs(eta), eps))))
        {	
            updateTransformation(transformation, 0, -1, 0, -1, 0, 0, 0, 0, -1);
			NiggliParams(transformation.transpose() * vectors, A, B, C, ksi, eta, zeta);
        }

		// Test 2
        if ((Num<double>::gt(B, C, eps)) || \
			((Num<double>::eq(B, C, eps)) && (Num<double>::gt(Num<double>::abs(eta), Num<double>::abs(zeta), eps))))
        {
            updateTransformation(transformation, -1, 0, 0, 0, 0, -1, 0, -1, 0);
			NiggliParams(transformation.transpose() * vectors, A, B, C, ksi, eta, zeta);
            continue;
        }
		
		// Test 3
		plusMin[0] = plusMin[1] = plusMin[2] = 0;
		if (Num<double>::lt(ksi, 0, eps))
			plusMin[0] = -1;
		else if (Num<double>::gt(ksi, 0, eps))
			plusMin[0] = 1;
		if (Num<double>::lt(eta, 0, eps))
			plusMin[1] = -1;
		else if (Num<double>::gt(eta, 0, eps))
			plusMin[1] = 1;
		if (Num<double>::lt(zeta, 0, eps))
			plusMin[2] = -1;
		else if (Num<double>::gt(zeta, 0, eps))
			plusMin[2] = 1;
		if (plusMin[0]*plusMin[1]*plusMin[2] == 1)
		{
            updateTransformation(transformation, plusMin[0], 0, 0, 0, plusMin[1], 0, 0, 0, plusMin[2]);
			NiggliParams(transformation.transpose() * vectors, A, B, C, ksi, eta, zeta);
        }
		
		// Test 4
		plusMin[0] = plusMin[1] = plusMin[2] = 0;
		if (Num<double>::lt(ksi, 0, eps))
			plusMin[0] = -1;
		else if (Num<double>::gt(ksi, 0, eps))
			plusMin[0] = 1;
		if (Num<double>::lt(eta, 0, eps))
			plusMin[1] = -1;
		else if (Num<double>::gt(eta, 0, eps))
			plusMin[1] = 1;
		if (Num<double>::lt(zeta, 0, eps))
			plusMin[2] = -1;
		else if (Num<double>::gt(zeta, 0, eps))
			plusMin[2] = 1;
		if ((plusMin[0]*plusMin[1]*plusMin[2] == 0) || (plusMin[0]*plusMin[1]*plusMin[2] == -1))
        {
			
			// Set values
			val[0] = (plusMin[0] == 1) ? -1 : 1;
			val[1] = (plusMin[1] == 1) ? -1 : 1;
			val[2] = (plusMin[2] == 1) ? -1 : 1;
			if (val[0]*val[1]*val[2] == -1)
			{
				if (plusMin[0] == 0)
					val[0] = -1;
				else if (plusMin[1] == 0)
					val[1] = -1;
				else if (plusMin[2] == 0)
					val[2] = -1;
			}
            updateTransformation(transformation, val[0], 0, 0, 0, val[1], 0, 0, 0, val[2]);
			NiggliParams(transformation.transpose() * vectors, A, B, C, ksi, eta, zeta);

			// Check minimum reduction
			if (!firstMinRedTest)
			{
				if (((minRedMult * A + (A - lastA)) - (A * minRedMult) == 0) && \
					((minRedMult * B + (B - lastB)) - (B * minRedMult) == 0) && \
					((minRedMult * C + (C - lastC)) - (C * minRedMult) == 0))
				{
					if ((minRedOn) && (count > 100))
						break;
					minRedOn = true;
				}
				else
					minRedOn = false;
			}
			lastA = A;
			lastB = B;
			lastC = C;
			firstMinRedTest = false;
        }
		
		// Test 5
        if ((Num<double>::gt(Num<double>::abs(ksi), B, eps)) || \
            ((Num<double>::eq(ksi, B, eps)) && (Num<double>::lt(2*eta, zeta, eps))) || \
            ((Num<double>::eq(ksi, -B, eps)) && (Num<double>::lt(zeta, 0, eps))))
        {
            updateTransformation(transformation, 1, 0, 0, 0, 1, -Num<double>::sign(ksi), 0, 0, 1);
			NiggliParams(transformation.transpose() * vectors, A, B, C, ksi, eta, zeta);
            continue;
        }

		// Test 6
        if ((Num<double>::gt(Num<double>::abs(eta), A, eps)) || \
            ((Num<double>::eq(eta, A, eps)) && (Num<double>::lt(2*ksi, zeta, eps))) || \
            ((Num<double>::eq(eta, -A, eps)) && (Num<double>::lt(zeta, 0, eps))))
		{
            updateTransformation(transformation, 1, 0, -Num<double>::sign(eta), 0, 1, 0, 0, 0, 1);
			NiggliParams(transformation.transpose() * vectors, A, B, C, ksi, eta, zeta);
            continue;
        }

		// Test 7
        if ((Num<double>::gt(Num<double>::abs(zeta), A, eps)) || \
            ((Num<double>::eq(zeta, A, eps)) && (Num<double>::lt(2*ksi, eta, eps))) || \
            ((Num<double>::eq(zeta, -A, eps)) && (Num<double>::lt(eta, 0, eps))))
        {
            updateTransformation(transformation, 1, -Num<double>::sign(zeta), 0, 0, 1, 0, 0, 0, 1);
			NiggliParams(transformation.transpose() * vectors, A, B, C, ksi, eta, zeta);
            continue;
        }

		// Test 8
        if ((Num<double>::lt(ksi + eta + zeta + A + B, 0, eps)) || \
			((Num<double>::eq(ksi + eta + zeta + A + B, 0, eps)) && (Num<double>::gt(2*(A + eta) + zeta, 0, eps))))
        {
            updateTransformation(transformation, 1, 0, 1, 0, 1, 1, 0, 0, 1);
			NiggliParams(transformation.transpose() * vectors, A, B, C, ksi, eta, zeta);
            continue;
        }

		// Finished if at this point
        break;
	}
	
	// Could not find reduced cell
	if (count == maxLoops)
	{
		Output::newline(ERROR);
        Output::print("Niggli reduction reached limit of ");
		Output::print(maxLoops);
		Output::print(" loop");
		if (maxLoops != 1)
			Output::print("s");
        Output::quit();
	}
	
    // Save result as transpose of transformation matrix
	return transformation.transpose();
}



/* void Basis::NiggliParams(const Matrix3D& vectors, double& A, double& B, double& C, double& ksi,
 *		double& eta, double& zeta)
 *
 * Set the metrics used in Niggli reduction
 */

void Basis::NiggliParams(const Matrix3D& vectors, double& A, double& B, double& C, double& ksi, \
	double& eta, double& zeta)
{
	
	// Basis lengths
	double a = sqrt(vectors(0, 0)*vectors(0, 0) + vectors(0, 1)*vectors(0, 1) + vectors(0, 2)*vectors(0, 2));
    double b = sqrt(vectors(1, 0)*vectors(1, 0) + vectors(1, 1)*vectors(1, 1) + vectors(1, 2)*vectors(1, 2));
    double c = sqrt(vectors(2, 0)*vectors(2, 0) + vectors(2, 1)*vectors(2, 1) + vectors(2, 2)*vectors(2, 2));

	// Basis angles
	Vector3D vec0(vectors(0, 0), vectors(0, 1), vectors(0, 2));
	Vector3D vec1(vectors(1, 0), vectors(1, 1), vectors(1, 2));
	Vector3D vec2(vectors(2, 0), vectors(2, 1), vectors(2, 2));
	double alpha = vec1.angle(vec2);
    double beta  = vec0.angle(vec2);
    double gamma = vec0.angle(vec1);
    
    // Set the main variables
    A = a * a;
    B = b * b;
    C = c * c;
    ksi =  2 * b * c * cos(alpha);
    eta =  2 * a * c * cos(beta);
    zeta = 2 * a * b * cos(gamma);
}



/* void Basis::updateTransformation(Matrix3D& transformation, int m00, int m01, int m02, int m10, int m11,
 *		int m12, int m20, int m21, int m22)
 *
 * Update the transformation matrix during Niggli reduction
 */

void Basis::updateTransformation(Matrix3D& transformation, int m00, int m01, int m02, int m10, int m11, \
	int m12, int m20, int m21, int m22)
{
	
	// Form the transformation matrix
	Matrix3D updateMatrix(m00, m01, m02, m10, m11, m12, m20, m21, m22);
	if (updateMatrix.determinant() < 0)
	{
		Output::newline(ERROR);
		Output::print("Negative transformation during Niggli reduction");
		Output::quit();
	}
	
	// Set the new transformation matrix
	transformation = transformation * updateMatrix;
}



/* void Basis::getPossibleRotations(Linked<Matrix3D>& rotations, const Matrix3D& vectors, double tol)
 *
 * Get list of allowed rotational symmetry operations
 */

void Basis::getPossibleRotations(Linked<Matrix3D>& rotations, const Matrix3D& vectors, double tol)
{
	
	// Clear space for result
	rotations.clear();
	
	// Set tolerance
	double angTol = atan(tol / pow(vectors.volume(), 1.0/3.0));
	
	// Set transformations
	Matrix3D realFracToCart = vectors.transpose();
	Matrix3D recipFracToCart = vectors.inverse();
	
	// Make list of possible two-fold operations
	Linked<Matrix3D> twoFoldRotations;
	twoFoldRotations.add(Matrix3D(-1, -1, -1, 0, 0, 1, 0, 1, 0));
	twoFoldRotations.add(Matrix3D(-1, -1, 0, 0, 1, 0, 0, -1, -1));
	twoFoldRotations.add(Matrix3D(-1, -1, 0, 0, 1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(-1, -1, 0, 0, 1, 0, 0, 1, -1));
	twoFoldRotations.add(Matrix3D(-1, -1, 1, 0, 0, -1, 0, -1, 0));
	twoFoldRotations.add(Matrix3D(-1, 0, -1, 0, -1, -1, 0, 0, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, -1, 0, -1, 0, 0, 0, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, -1, 0, -1, 1, 0, 0, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, -1, 0, -1, 1, -1, 0));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, -1, 0, 1, -1, 1, 0));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, -1, 1, -1, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, -1, 1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, -1, 1, 1, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, -1, -1, 0, 0, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, -1, 0, -1, -1, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, -1, 0, -1, 0, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, -1, 0, -1, 1, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, -1, 0, 0, -1, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, -1, 0, 0, 0, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, -1, 0, 0, 1, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, -1, 0, 1, -1, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, -1, 0, 1, 0, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, -1, 0, 1, 1, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, -1, 1, 0, 0, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, 0, -1, 0, -1, 0));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, 0, 1, 0, 1, 0));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, 1, -1, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, 1, 0, 0, -1, -1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, 1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, 1, 0, 0, 1, -1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 0, 1, 1, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 1, 0, -1, -1, -1, 0));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 1, 0, 1, 1, 1, 0));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 1, 1, -1, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 1, 1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(-1, 0, 0, 1, 1, 1, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(-1, 0, 1, 0, -1, -1, 0, 0, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, 1, 0, -1, 0, 0, 0, 1));
	twoFoldRotations.add(Matrix3D(-1, 0, 1, 0, -1, 1, 0, 0, 1));
	twoFoldRotations.add(Matrix3D(-1, 1, -1, 0, 0, -1, 0, -1, 0));
	twoFoldRotations.add(Matrix3D(-1, 1, 0, 0, 1, 0, 0, -1, -1));
	twoFoldRotations.add(Matrix3D(-1, 1, 0, 0, 1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(-1, 1, 0, 0, 1, 0, 0, 1, -1));
	twoFoldRotations.add(Matrix3D(-1, 1, 1, 0, 0, 1, 0, 1, 0));
	twoFoldRotations.add(Matrix3D(0, -1, -1, -1, 0, 1, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(0, -1, -1, 0, -1, 0, -1, 1, 0));
	twoFoldRotations.add(Matrix3D(0, -1, 0, -1, 0, 0, -1, 1, -1));
	twoFoldRotations.add(Matrix3D(0, -1, 0, -1, 0, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(0, -1, 0, -1, 0, 0, 1, -1, -1));
	twoFoldRotations.add(Matrix3D(0, -1, 1, -1, 0, -1, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(0, -1, 1, 0, -1, 0, 1, -1, 0));
	twoFoldRotations.add(Matrix3D(0, 0, -1, -1, -1, 1, -1, 0, 0));
	twoFoldRotations.add(Matrix3D(0, 0, -1, 0, -1, 0, -1, 0, 0));
	twoFoldRotations.add(Matrix3D(0, 0, -1, 1, -1, -1, -1, 0, 0));
	twoFoldRotations.add(Matrix3D(0, 0, 1, -1, -1, -1, 1, 0, 0));
	twoFoldRotations.add(Matrix3D(0, 0, 1, 0, -1, 0, 1, 0, 0));
	twoFoldRotations.add(Matrix3D(0, 0, 1, 1, -1, 1, 1, 0, 0));
	twoFoldRotations.add(Matrix3D(0, 1, -1, 0, -1, 0, -1, -1, 0));
	twoFoldRotations.add(Matrix3D(0, 1, -1, 1, 0, -1, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(0, 1, 0, 1, 0, 0, -1, -1, -1));
	twoFoldRotations.add(Matrix3D(0, 1, 0, 1, 0, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(0, 1, 0, 1, 0, 0, 1, 1, -1));
	twoFoldRotations.add(Matrix3D(0, 1, 1, 0, -1, 0, 1, 1, 0));
	twoFoldRotations.add(Matrix3D(0, 1, 1, 1, 0, 1, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(1, -1, -1, 0, -1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(1, -1, 0, 0, -1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(1, -1, 1, 0, -1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(1, 0, -1, 0, -1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(1, 0, 0, -1, -1, 0, -1, 0, -1));
	twoFoldRotations.add(Matrix3D(1, 0, 0, -1, -1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(1, 0, 0, -1, -1, 0, 1, 0, -1));
	twoFoldRotations.add(Matrix3D(1, 0, 0, 0, -1, 0, -1, 0, -1));
	twoFoldRotations.add(Matrix3D(1, 0, 0, 0, -1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(1, 0, 0, 0, -1, 0, 1, 0, -1));
	twoFoldRotations.add(Matrix3D(1, 0, 0, 1, -1, 0, -1, 0, -1));
	twoFoldRotations.add(Matrix3D(1, 0, 0, 1, -1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(1, 0, 0, 1, -1, 0, 1, 0, -1));
	twoFoldRotations.add(Matrix3D(1, 0, 1, 0, -1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(1, 1, -1, 0, -1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(1, 1, 0, 0, -1, 0, 0, 0, -1));
	twoFoldRotations.add(Matrix3D(1, 1, 1, 0, -1, 0, 0, 0, -1));
	
	// Loop over possible two-fold operations
	bool found;
	double delta;
	Vector3D realAxis;
	Vector3D recipAxis;
	Linked<double> deltaList;
	Linked<Matrix3D> matrixList;
	Linked<double>::iterator itDelta;
	Linked<Matrix3D>::iterator itMatrix;
	for (Linked<Matrix3D>::iterator mat = twoFoldRotations.begin(); mat != twoFoldRotations.end(); ++mat)
	{
		
		// Get real and reciprocal space lattice directions
		realAxis  = realFracToCart  * rotationAxis(*mat);
		recipAxis = recipFracToCart * rotationAxis((*mat).transpose());
		
		// Get angle between real and reciprocal vectors
		delta = Num<double>::abs(getAngle(realAxis, recipAxis));
		
		// Save if within tolerance
		if (delta < angTol)
		{
			
			// Check if delta is lower than known matrices
			found = false;
			itDelta = deltaList.begin();
			itMatrix = matrixList.begin();
			for (; itDelta != deltaList.end(); ++itDelta, ++itMatrix)
			{
				if (delta < *itDelta)
				{
					deltaList.addBefore(itDelta, delta);
					matrixList.addBefore(itMatrix, *mat);
					found = true;
					break;
				}
			}
			
			// Delta is larger than all known
			if (!found)
			{
				deltaList += delta;
				matrixList += *mat;
			}
		}
	}

	// Recursively create list of possible symmetry operations
	Linked<Matrix3D>::iterator prevLast;
	for (itMatrix = matrixList.begin(); itMatrix != matrixList.end(); ++itMatrix)
	{
		
		// Try adding rotation
		prevLast = rotations.last();
		if (!addPossibleRotation(rotations, *itMatrix))
		{

			// Remove operations that were added
			while (prevLast != rotations.last())
				rotations.remove(prevLast + 1);
		}
	}
	
	// Add negative of each operation
	int i;
	int length = rotations.length();
	for (i = 0, itMatrix = rotations.begin(); i < length; ++i, ++itMatrix)
		rotations += *itMatrix * -1;
}



/* bool Basis::addPossibleRotation(Linked<Matrix3D>& rotations, const Matrix3D& rotationToAdd)
 *
 * Add symmetry operation to set
 */

bool Basis::addPossibleRotation(Linked<Matrix3D>& rotations, const Matrix3D& rotationToAdd)
{
	
	// Loop over rotations
	Linked<Matrix3D>::iterator itRot;
	for (itRot = rotations.begin(); itRot != rotations.end(); ++itRot)
	{
		
		// Rotation already exists
		if (*itRot == rotationToAdd)
			return true;
	}
	
	// Check if operation is allowed
	int i, j;
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			if (Num<double>::abs(rotationToAdd(i, j)) > 1.5)
				return false;
		}
	}
	
	// Save if operation is not identity
	if (rotationToAdd != Matrix3D::identity())
		rotations += rotationToAdd;
	
	// Generate new operations if rotation was not found
	Matrix3D newRotation;
	for (itRot = rotations.begin(); itRot != rotations.end(); ++itRot)
	{
		
		// Generate forward operation
		newRotation = *itRot * rotationToAdd;
		if (!addPossibleRotation(rotations, newRotation))
			return false;
		
		// Generate reverse operation
		newRotation = rotationToAdd * *itRot;
		if (!addPossibleRotation(rotations, newRotation))
			return false;
	}
	
	// Return true if at this point
	return true;
}



/* void ISO:: fillRotations(Linked<Matrix3D>& rotations, Linked<Matrix3D>& origRotations)
 *
 * Recursively build list of possible rotations
 */

void ISO::fillRotations(Linked<Matrix3D>& rotations, Linked<Matrix3D>& origRotations)
{
	
	// Build list of rotations
	rotations.clear();
	for (Linked<Matrix3D>::iterator itOrig = origRotations.begin(); itOrig != origRotations.end(); ++itOrig)
		addRotation(rotations, *itOrig);
	addRotation(rotations, Matrix3D::identity());
	addRotation(rotations, Matrix3D::identity()*-1);
	
	// Move identity to first position
	Matrix3D identity = Matrix3D::identity();
	for (Linked<Matrix3D>::iterator itOrig = origRotations.begin(); itOrig != origRotations.end(); ++itOrig)
	{
		if (*itOrig == identity)
		{
			if (itOrig != origRotations.begin())
			{
				*itOrig = *origRotations.begin();
				*origRotations.begin() = identity;
			}
			break;
		}
	}
}



/* void ISO::addRotation(Linked<Matrix3D>& rotations, const Matrix3D& curRotation)
 *
 * Add operation to list
 */

void ISO::addRotation(Linked<Matrix3D>& rotations, const Matrix3D& curRotation)
{
	
	// Check if current rotation was already saved
	Linked<Matrix3D>::iterator it = rotations.begin();
	for (; it != rotations.end(); ++it)
	{
		if (*it == curRotation)
			break;
	}
	if (it != rotations.end())
		return;
	
	// Save rotation
	rotations += curRotation;
	
	// Save inverse of operation
	addRotation(rotations, curRotation.inverse());
	
	// Save multiples of rotations
	for (it = rotations.begin(); it != rotations.end(); ++it)
	{
		addRotation(rotations, *it * curRotation);
		addRotation(rotations, curRotation * *it);
	}
}



/* Atom::Atom()
 *
 * Constructor for Atom object
 */

Atom::Atom()
{
	_assigned = false;
	_fixed[0] = _fixed[1] = _fixed[2] = false;
	_interstitial = false;
	_cartesian = 0.0;
	_fractional = 0.0;
	_magneticMoment = 0.0;
	_occupancy = 1;
	_atomNumber = -1;
	_basis = 0;
}



/* void Atom::clear()
 *
 * Clear data in Atom object
 */

void Atom::clear()
{
	_tags.clear();
	_element.clear();
	_assigned = false;
	_fixed[0] = _fixed[1] = _fixed[2] = false;
	_interstitial = false;
	_cartesian = 0.0;
	_fractional = 0.0;
	_magneticMoment = 0.0;
	_occupancy = 1;
	_atomNumber = -1;
	_basis = 0;
	
}



/* Atom& Atom::operator= (const Atom& rhs)
 *
 * Assignment operator for Atom object
 */

Atom& Atom::operator= (const Atom& rhs)
{
	if (this != &rhs)
	{
		_assigned = rhs._assigned;
		_fixed[0] = rhs._fixed[0];
		_fixed[1] = rhs._fixed[1];
		_fixed[2] = rhs._fixed[2];
		_interstitial = rhs._interstitial;
		_occupancy = rhs._occupancy;
		_cartesian = rhs._cartesian;
		_fractional = rhs._fractional;
		_magneticMoment = rhs._magneticMoment;
		_element = rhs._element;
		_tags = rhs._tags;
		_atomNumber = rhs._atomNumber;
		_basis = rhs._basis;
	}
	return *this;
}



/* bool Atom::equal(const Atom& rhs, double tol, double* distance, Vector3D* cell) const
 *
 * Return whether two atoms are equal
 */

bool Atom::equal(const Atom& rhs, double tol, double* distance, Vector3D* cell) const
{
	
	// Check the element
	if (_element != rhs._element)
		return false;
	
	// Check if both are or are not interstitials
	if (_interstitial != rhs._interstitial)
		return false;
	
	// Check the magnetic moment
	if ((_magneticMoment[0] - rhs._magneticMoment[0])*(_magneticMoment[0] - rhs._magneticMoment[0]) + \
		(_magneticMoment[1] - rhs._magneticMoment[1])*(_magneticMoment[1] - rhs._magneticMoment[1]) + \
		(_magneticMoment[2] - rhs._magneticMoment[2])*(_magneticMoment[2] - rhs._magneticMoment[2]) > 2*tol)
		return false;
	
	// Check the positions
	double tempDis = _basis->distance(_fractional, FRACTIONAL, rhs._fractional, FRACTIONAL, cell);
	if (distance)
		*distance = tempDis;
	if (tempDis > tol)
		return false;
	
	// Return that atoms are the same if at this point
	return true;
}



/* void ISO::clear()
 *
 * Clear data is ISO object
 */

void ISO::clear()
{
	_comment.clear();
	_atoms.clear();
	_numAtoms = 0;
	_spaceGroup.clear();
}



/* ISO& ISO::operator= (const ISO& rhs)
 *
 * Assignment operator for ISO object
 */

ISO& ISO::operator= (const ISO& rhs)
{
	
	// Make sure not copying self
	if (this != &rhs)
	{
		
		// Clear space
		clear();
		
		// Save properties
		_comment = rhs._comment;
		_basis = rhs._basis;
		_spaceGroup = rhs._spaceGroup;
		
		// Save atoms
		int i, j;
		_numAtoms = rhs._numAtoms;
		_atoms.length(rhs._atoms.length());
		for (i = 0; i < rhs._atoms.length(); ++i)
		{
			_atoms[i].length(rhs._atoms[i].length());
			for (j = 0; j < rhs._atoms[i].length(); ++j)
			{
				_atoms[i][j] = rhs._atoms[i][j];
				_atoms[i][j].basis(&_basis);
			}
		}
	}
	
	// Return result
	return *this;
}



/* void ISO::basis(const Basis& basis, bool showOutput)
 *
 * Set the basis for a structure - if atoms are present, then just keep fractional coordinates the same
 */

void ISO::basis(const Basis& basis, bool showOutput)
{
	
	// Set the basis
	_basis.set(basis.vectors(), showOutput);
	
	// Set the atom positions
	int i, j;
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			if (_atoms[i][j].assigned())
				_atoms[i][j].fractional(_atoms[i][j].fractional());
		}
	}
}



/* void ISO::basis(const Matrix3D& vectors, bool showOutput)
 *
 * Set the basis for a structure
 */

void ISO::basis(const Matrix3D& vectors, bool showOutput)
{
	
	// Set the basis
	_basis.set(vectors, showOutput);
	
	// Set the atom positions
	int i, j;
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			if (_atoms[i][j].assigned())
				_atoms[i][j].fractional(_atoms[i][j].fractional());
		}
	}
}



/* void ISO::basis(const Vector3D& lengths, const Vector3D& angles, bool showOutput)
 *
 * Set the basis for a structure
 */

void ISO::basis(const Vector3D& lengths, const Vector3D& angles, bool showOutput)
{
	
	// Set the basis
	_basis.set(lengths, angles, showOutput);
	
	// Set the atom positions
	int i, j;
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			if (_atoms[i][j].assigned())
				_atoms[i][j].fractional(_atoms[i][j].fractional());
		}
	}
}



/* Matrix3D ISO::makeReduced(bool showOutput)
 *
 * Convert atoms to reduced positions
 */

Matrix3D ISO::makeReduced(bool showOutput)
{
	
	// Output
	if (!showOutput)
		Output::quietOn();
	
	// Output
	Output::newline();
	Output::print("Setting the cell to reduced form");
	Output::increase();
	
	// Save transformation for atoms to reduced cell
	Matrix3D transformation = _basis.unitPointToReduced();
	
	// Convert basis to reduced form
	Matrix3D origUnitToReduced = _basis.unitToReduced();
	_basis.set(_basis.unitToReduced() * _basis.vectors(), showOutput);
	
	// Output
	Output::newline();
	Output::print("Setting atom positions");
	Output::increase();
	
	// Save new atom positions
	int i, j;
	Vector3D tempPos;
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			
			// Output
			Output::newline();
			Output::print("Setting atom ");
			Output::print(_atoms[i][j].atomNumber() + 1);
			Output::increase();
			
			// Set position
			tempPos = transformation * _atoms[i][j].fractional();
			_atoms[i][j].fractional(tempPos);
			
			// Output
			Output::newline();
			Output::print("Fractional: ");
			Output::print(_atoms[i][j].fractional(), 8, true);
			Output::newline();
			Output::print("Cartesian: ");
			Output::print(_atoms[i][j].cartesian(), 8, true);
			Output::decrease();
		}
	}
	
	// Output
	Output::decrease();
	Output::decrease();
	if (!showOutput)
		Output::quietOff();
	
	// Return transformation
	return origUnitToReduced;
}



/* Matrix3D ISO::primitiveTransformation(double tol, bool showOutput) const
 *
 * Get the transformation from unit to primitive cell
 */

Matrix3D ISO::primitiveTransformation(double tol, bool showOutput) const
{
	
	// Output
	if (!showOutput)
		Output::quietOn();
	
	// Output
	Output::newline();
	Output::print("Searching for the transformation from unit cell to primitive cell");
	Output::increase();
	
	// Figure out which element has the fewest occurances
	int i;
	int minElem = 0;
	for (i = 1; i < _atoms.length(); ++i)
	{
		if (_atoms[i].length() < _atoms[minElem].length())
			minElem = i;
	}
	
	// Get the possible vectors
	int j;
	OList<Vector3D > vectors (_atoms[minElem].length() - 1);
	for (i = 1; i < _atoms[minElem].length(); ++i)
	{
		vectors[i-1] = _atoms[minElem][i].fractional() - _atoms[minElem][0].fractional();
		ISO::moveIntoCell(vectors[i-1]);
		for (j = 0; j < 3; ++j)
		{
			if (vectors[i-1][j] > 0.5)
				vectors[i-1][j] -= 1;
		}
	}
	
	// Get the distance of all atoms from the origin
	Vector3D origin = 0.0;
	List<double>::D2 distances(_atoms.length());
	for (i = 0; i < _atoms.length(); ++i)
	{
		distances[i].length(_atoms[i].length());
		for (j = 0; j < _atoms[i].length(); ++j)
			distances[i][j] = _basis.distance(_atoms[i][j].fractional(), FRACTIONAL, origin, FRACTIONAL);
	}
	
	// List to store translated atoms
	Atoms transAtoms(_atoms);
	
	// Initialize current cell as primitive
	double volume = 1;
	Matrix3D transformation;
	transformation.makeIdentity();
	
	// Number of loops to make
	int numLoops = (int) Num<double>::ceil((double) vectors.length() / Multi::worldSize());
	
	// Loop over vectors
	int k, m;
	bool areEqual;
	bool areEqualRec;
	int index;
	double newVolume;
	double volTol = 1e-3;
	Vector3D vectorRec;
	Matrix3D newTransformation;
	for (i = 0; i < numLoops; ++i)
	{
		
		// Figure out current translation
		index = i * Multi::worldSize() + Multi::rank();
		
		// Check if sites map on current processor
		areEqual = false;
		if (index < vectors.length())
		{
					
			// Translate all atoms by current vector
			for (j = 0; j < _atoms.length(); ++j)
			{
				for (k = 0; k < _atoms[j].length(); ++k)
					transAtoms[j][k].fractional(_atoms[j][k].fractional() + vectors[index]);
			}

			// Check if sites map
			areEqual = areSitesEqual(_basis, _atoms, transAtoms, tol, &(vectors[index]), &distances);
		}
		
		// Loop over processors
		for (j = 0; j < Multi::worldSize(); ++j)
		{	
			
			// Broadcast result
			areEqualRec = areEqual;
			Multi::broadcast(areEqualRec, j);
			
			// Sites did not map
			if (!areEqualRec)
				continue;
			
			// Broadcast vector
			if (j == Multi::rank())
				vectorRec = vectors[index];
			Multi::broadcast(vectorRec, j);
			
			// Loop over basis vectors
			for (k = 0; k < 3; ++k)
			{
				
				// Save new basis
				newTransformation = transformation;
				for (m = 0; m < 3; ++m)
					newTransformation(k, m) = vectorRec[m];
				
				// Get the volume of the transformation
				newVolume = newTransformation.volume();
				
				// Swap last two rows if a negative volume
				if (newVolume < volTol)
				{
					newTransformation.swapRows(1, 2);
					newVolume *= -1;
				}
				
				// Skip if the volume is zero
				if (newVolume < volTol)
					continue;
				
				// Skip if volume is larger than current
				if (newVolume - volume > volTol)
					continue;
				
				// Save new basis
				volume = newVolume;
				transformation = newTransformation;
			}
		}
	}
	
	// Output
	Output::newline();
    Output::print("Relative to the unit cell, the primitive basis vectors are:");
	for (i = 0; i < 3; ++i)
	{
		Output::newline();
		Output::tab();
		for (j = 0; j < 3; ++j)
		{
			Output::print(transformation(i, j), 8);
			Output::print(" ");
		}
	}
	
	// Output
	Output::decrease();
	
	// Output
	if (!showOutput)
		Output::quietOff();
	
	// Return the transformation
	return transformation;
}



/* Matrix3D ISO::makePrimitive(double tol, bool showOutput)
 *
 * Convert to primitive cell
 */

Matrix3D ISO::makePrimitive(double tol, bool showOutput)
{
	
	// Output
	if (!showOutput)
		Output::quietOn();
	
	// Output
	Output::newline();
	Output::print("Converting cell to primitive form");
	Output::increase();
	
	// Get the transformation to the primitive cell
	Matrix3D transformation = primitiveTransformation(tol, showOutput);
	
	// Set the new structure
	transform(transformation, tol, showOutput);
	
	// Output
	Output::decrease();
	
	// Output
	if (!showOutput)
		Output::quietOff();
	
	// Return the transformation
	return transformation;
}



/* void ISO::transform(const Matrix3D& transformation, double tol, bool showOutput)
 *
 * Transform a structure
 */

void ISO::transform(const Matrix3D& transformation, double tol, bool showOutput)
{
	
	// Set output
	if (!showOutput)
		Output::quietOn();
	
	// Output
	int i, j;
	Output::newline();
	Output::print("Applying transformation: [");
	for (i = 0; i < 3; ++i)
	{
		Output::print("(");
		for (j = 0; j < 3; ++j)
		{
			Output::print(transformation(i, j));
			if (j != 2)
				Output::print(", ");
		}
		Output::print(")");
		if (i != 2)
			Output::print(", ");
	}
	Output::print("]");
	Output::increase();
	
	// Check if the identity matrix
	if (transformation == Matrix3D::identity())
	{
		Output::newline();
		Output::print("Transforming using the identity matrix so nothing needs to be done");
		Output::decrease();
		if (!showOutput)
			Output::quietOff();
		return;
	}
	
	// Figure out the number of cells created under transformation
	double numCells = transformation.determinant();
	
	// Set the new basis
	_basis.set(transformation * _basis.vectors());
	
	// Get the atoms in the new cell
	Atoms newAtoms = numCells < 0.999999 ? contract(transformation, tol) : expand(transformation);
	
	// Set the new atoms
	_numAtoms = 0;
	_atoms.clear();
	_atoms.length(newAtoms.length());
	for (i = 0; i < newAtoms.length(); ++i)
	{
		_atoms[i].length(newAtoms[i].length());
		for (j = 0; j < newAtoms[i].length(); ++j)
		{
			
			// Output
			Output::newline();
			Output::print("Adding atom ");
			Output::print(_numAtoms + 1);
			Output::print(" (");
			Output::print(newAtoms[i][j].element().symbol());
			Output::print(") at ");

			// Save atom properties
			_atoms[i][j] = newAtoms[i][j];
			_atoms[i][j].atomNumber(_numAtoms++);

			// Print information about the new atom
			Output::print(_atoms[i][j].fractional(), 8, true);
		}
	}
	
	// Reset output
	Output::decrease();
	if (!showOutput)
		Output::quietOff();
}



/* Atoms ISO::expand(const Matrix3D& transformation)
 * 
 * Get atom positions in expanded cell
 */

Atoms ISO::expand(const Matrix3D& transformation)
{
	
	// Get the lattice points of the original cell in the new cell
	LatticePoints points = getLatticePoints(transformation);
	
	// Output
	int i, j;
	Output::newline();
	Output::print(points.length());
	Output::print(" lattice point");
	if (points.length() != 1)
		Output::print("s");
	Output::print(" of original cell in new cell");
	for (i = 0; i < points.length(); ++i)
	{
		Output::newline();
		Output::tab();
		for (j = 0; j < 3; ++j)
		{
			Output::print(points[i][j], 0);
			Output::print(" ");
		}
	}
	
	// Get the transformation for atoms in original cell to new cell
	Matrix3D conversion = transformation.inverse().transpose();
	
	// Variable to store result
	Atoms newAtoms(_atoms.length());
	
	// Loop over atoms in the orginal cell
	int k, m;
	bool found;
	int pos;
	Vector3D newPos;
	for (i = 0; i < _atoms.length(); ++i)
	{
		pos = 0;
		newAtoms[i].length(points.length() * _atoms[i].length());
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			
			// Loop over lattice points of original cell in new cell
			for (k = 0; k < points.length(); ++k)
			{

				// Set current atom
				newAtoms[i][pos] = _atoms[i][j];
				newPos = conversion * (points[k] + _atoms[i][j].fractional());
				newAtoms[i][pos].fractional(newPos);
				
				// Check if atom is already known
				found = false;
				for (m = 0; m < pos; ++m)
				{
					if (_basis.distance(newAtoms[i][pos].fractional(), FRACTIONAL, newAtoms[i][m].fractional(), \
						FRACTIONAL) < 1e-4)
					{
						found = true;
						break;
					}
				}
				if (found)
					continue;
				
				// Go to next atom
				++pos;
			}
		}
		
		// Set length of atoms
		newAtoms[i].length(pos);
	}
	
	// Return the new atoms
	return newAtoms;
}



/* Atoms ISO::contract(const Matrix3D& transformation, double tol)
 *
 * Get atom positions in contracted positions
 */

Atoms ISO::contract(const Matrix3D& transformation, double tol)
{
	
	// Get the transformation for atoms in original cell to new cell
	Matrix3D conversion = transformation.inverse().transpose();
	
	// Get the number of cells of original in new
	double numCells = transformation.determinant();
	
	// Output
	Output::newline();
	Output::print("Contracting cell and removing atoms");
	
	// Variable to store results
	Atoms newAtoms (_atoms.length());
	List<int>::D2 timesFound (_atoms.length());
	List<double>::D2 newDis (_atoms.length());
	
	// Loop over atoms in original cell
	int i, j, k, m;
	bool isNew;
	int curAtom;
	int numAtoms;
	double curDis;
	Vector3D nearCell;
	Vector3D newPos;
	Vector3D newMagMom;
	Vector3D origin(0.0);
	for (i = 0; i < _atoms.length(); ++i)
	{
		
		// Loop over atoms of current element
		curAtom = 0;
		numAtoms = (int) Num<double>::round(_atoms[i].length() * numCells, 1);
		newAtoms[i].length(numAtoms);
		newDis[i].length(numAtoms);
		timesFound[i].length(numAtoms);
		timesFound[i].fill(0);
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			
			// Get new position
			newPos = conversion * _atoms[i][j].fractional();
			ISO::moveIntoCell(newPos);
			curDis = _basis.distance(newPos, FRACTIONAL, origin, FRACTIONAL);

			// Loop over atoms of current element
			isNew = true;
			for (k = 0; k < curAtom; ++k)
			{
				
				// Check if positions are the same
				if (Num<double>::abs(curDis - newDis[i][k]) > tol)
					continue;
				if (_basis.distance(newAtoms[i][k].fractional(), FRACTIONAL, newPos, FRACTIONAL, &nearCell) > tol)
					continue;
				
				// Positions are the same but magnetic moments are different
				if (_basis.absoluteDistance(_atoms[i][j].magneticMoment(), CARTESIAN, \
					newAtoms[i][k].magneticMoment(), CARTESIAN) > tol)
				{
					Output::newline(ERROR);
					Output::print("Atoms generated at same position but with different magnetic moments");
					Output::quit();
				}
				
				// Average values
				for (m = 0; m < 3; ++m)
				{
					newPos[m] = (newPos[m] + nearCell[m] + timesFound[i][k] * newAtoms[i][k].fractional()[m]) /\
						(timesFound[i][k] + 1);
					newMagMom[m] = (_atoms[i][j].magneticMoment()[m] + \
						timesFound[i][k] * newAtoms[i][k].magneticMoment()[m]) / (timesFound[i][k] + 1);
				}
				
				// Update atom properties
				ISO::moveIntoCell(newPos);
				newAtoms[i][k].fractional(newPos);
				newAtoms[i][k].magneticMoment(newMagMom);
				timesFound[i][k]++;
				newDis[i][k] = _basis.distance(newAtoms[i][k].fractional(), FRACTIONAL, origin, FRACTIONAL);
				for (m = 0; m < 3; ++m)
				{
					if (_atoms[i][j].fixed()[m])
						newAtoms[i][k].fixed(m, true);
				}

				// Break since atom was found
				isNew = false;
				break;
			}

			// Atom was already known
			if (!isNew)
				continue;
				
			// An error occured if at this point
			if (curAtom >= numAtoms)
			{
				Output::newline(ERROR);
				Output::print("Transformation is not commensurate with the lattice");
				Output::quit();
			}
			
			// Save new atom
			newAtoms[i][curAtom] = _atoms[i][j];
			newAtoms[i][curAtom].fractional(newPos);
			timesFound[i][curAtom] = 1;
			newDis[i][curAtom] = curDis;
			curAtom++;
		}
		
		// Make sure that the number of atoms is not too small
		if (curAtom < numAtoms)
		{
			int n, p;
			for (j = -1; j <= 1; ++j)
			{
				for (k = -1; k <= 1; ++k)
				{
					for (m = -1; m <= 1; ++m)
					{
						
						// Skip if generating the original cell
						if ((!j) && (!k) && (!m))
							continue;
						
						// Loop over atoms in the structure
						for (n = 0; n < _atoms[i].length(); ++n)
						{
							
							// Generate new position
							newPos = _atoms[i][n].fractional();
							newPos[0] += j;
							newPos[1] += k;
							newPos[2] += m;
							newPos *= conversion;
							ISO::moveIntoCell(newPos);
							curDis = _basis.distance(newPos, FRACTIONAL, origin, FRACTIONAL);
						
							// Check if this is a new position
							for (p = 0; p < curAtom; ++p)
							{
								if (Num<double>::abs(curDis - newDis[i][p]) < tol)
								{
									if (_basis.distance(newAtoms[i][p].fractional(), FRACTIONAL, newPos, \
										FRACTIONAL) < tol)
										break;
								}
							}
						
							// Found a new atom
							if (p >= curAtom)
							{
								newAtoms[i][curAtom] = _atoms[i][n];
								newAtoms[i][curAtom].fractional(newPos);
								timesFound[i][curAtom] = 1;
								newDis[i][curAtom] = curDis;
								if (++curAtom == numAtoms)
									break;
							}
						}
						if (curAtom == numAtoms)
							break;
					}
					if (curAtom == numAtoms)
						break;
				}
				if (curAtom == numAtoms)
					break;
			}
		}
		
		// Did not find enough atoms
		if (curAtom < numAtoms)
		{
			Output::newline(ERROR);
			Output::print("Did not generate enough atoms during transformation");
			Output::quit();
		}
	}
	
	// Return result
	return newAtoms;
}



/* void ISO::rotatePositions(const Matrix3D& rotation, bool showOutput)
 *
 * Rotate positions of atoms keeping the basis fixed
 */

void ISO::rotatePositions(const Matrix3D& rotation, bool showOutput)
{
	// Set output
	if (!showOutput)
		Output::quietOn();
	
	// Output
	int i, j;
	Output::newline();
	Output::print("Applying rotation: ");
	for (i = 0; i < 3; ++i)
	{
		Output::print("(");
		for (j = 0; j < 3; ++j)
		{
			Output::print(rotation(i, j));
			if (j != 2)
				Output::print(", ");
		}
		Output::print(")");
		if (i != 2)
			Output::print(", ");
	}
	Output::print("]");
	Output::increase();
	
	// Rotate all positions
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
			_atoms[i][j].fractional(rotation * _atoms[i][j].fractional());
	}
	
	// Reset output
	Output::decrease();
	if (!showOutput)
		Output::quietOff();
}



/* void ISO::shift(const Vector3D& vector, bool showOutput)
 *
 * Shift atom positions
 */

void ISO::shift(const Vector3D& vector, bool showOutput)
{
	
	// Set output
	if (!showOutput)
		Output::quietOn();
	
	// Output
	Output::newline();
	Output::print("Shifting positions by ");
	Output::print(vector, 8, true);
	Output::increase();
	
	// Shift all atom positions
	int i, j;
	Vector3D tempPos;
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			tempPos = _atoms[i][j].fractional();
			tempPos += vector;
			_atoms[i][j].fractional(tempPos);
		}
	}
	
	// Reset output
	Output::decrease();
	if (!showOutput)
		Output::quietOff();
}



/* Atom* ISO::addAtom(const Element& element)
 *
 * Add atom to the structure
 */

Atom* ISO::addAtom(const Element& element)
{
	
	// Pointer to added atom
	Atom* newAtom = 0;
	
	// Loop over known elements
	int i;
	bool isNewElement = true;
	for (i = 0; i < _atoms.length(); ++i)
	{
		
		// Found current element
		if (element == _atoms[i][0].element())
		{
			
			// Add atom
			_atoms[i].add();
			_atoms[i].last().basis(&_basis);
			_atoms[i].last().element(element);
			_atoms[i].last().atomNumber(_numAtoms);
			newAtom = &(_atoms[i].last());
			
			// Break since found
			isNewElement = false;
			break;
		}
	}
	
	// Found a new element
	if (isNewElement)
	{
		_atoms.add();
		_atoms.last().add();
		_atoms.last().last().basis(&_basis);
		_atoms.last().last().element(element);
		_atoms.last().last().atomNumber(_numAtoms);
		newAtom = &(_atoms.last().last());
	}
	
	// Save that number of atoms increased
	++_numAtoms;
	
	// Return pointer to added atom
	return newAtom;
}



/* Atom* ISO::addAtom(const Atom& atom)
 *
 * Add atom to structure
 */

Atom* ISO::addAtom(const Atom& atom)
{
	
	// Add by element
	Atom* newAtom = addAtom(atom.element());
	
	// Save atom properties
	newAtom->fractional(atom.fractional());
	newAtom->fixed(atom.fixed());
	newAtom->occupancy(atom.occupancy());
	newAtom->isInterstitial(atom.isInterstitial());
	newAtom->magneticMoment(atom.magneticMoment());
	for (int i = 0; i < atom.tags().length(); ++i)
		newAtom->addTag(atom.tags()[i]);
	
	// Return pointer to atom
	return newAtom;
}



/* void ISO::removeAtom(int index, bool showOutput)
 *
 * Remove atom from the structure
 */

void ISO::removeAtom(int index, bool showOutput)
{
	
	// Look for atom in structure
	int i, j;
	bool found = false;
	int elemIndex = 0;
	int atomIndex = 0;
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			if (_atoms[i][j].atomNumber() == index)
			{
				elemIndex = i;
				atomIndex = j;
				found = true;
				break;
			}
		}
		if (found)
			break;
	}
	
	// Atom is out of range
	if (!found)
	{
		Output::newline(WARNING);
		Output::print("Atom ");
		Output::print(index + 1);
		Output::print(" is out of range so cannot be removed");
		return;
	}
	
	// Output
	if (!showOutput)
		Output::quietOn();
	
	// Output
	Output::newline();
	Output::print("Removing atom ");
	Output::print(index + 1);
	Output::print(" from the structure");
	
	// Remove atom
	_atoms[elemIndex].remove(atomIndex);
	--_numAtoms;
	
	// Reset atom numbers
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			if (_atoms[i][j].atomNumber() > index)
				_atoms[i][j].atomNumber(_atoms[i][j].atomNumber() - 1);
		}
	}
	
	// No elements left of type
	if (!_atoms[elemIndex].length())
		_atoms.remove(elemIndex);
	
	// Output
	if (!showOutput)
		Output::quietOff();
}



/* void ISO::setElement(int index, const Element& element, bool showOutput)
 *
 * Set the element of an atom
 */

void ISO::setElement(int index, const Element& element, bool showOutput)
{
	
	// Look for atom in structure
	int i, j;
	bool found = false;
	int elemIndex = 0;
	int atomIndex = 0;
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			if (_atoms[i][j].atomNumber() == index)
			{
				elemIndex = i;
				atomIndex = j;
				found = true;
				break;
			}
		}
		if (found)
			break;
	}
	
	// Atom is out of range
	if (!found)
	{
		Output::newline(WARNING);
		Output::print("Atom ");
		Output::print(index + 1);
		Output::print(" is out of range so element cannot be set");
		return;
	}
	
	// Output
	if (!showOutput)
		Output::quietOn();
	
	// Output
	Output::newline();
	Output::print("Setting element of atom ");
	Output::print(index + 1);
	Output::print(" to ");
	Output::print(element.symbol());
	
	// Atom is already of current element
	if (element == _atoms[elemIndex][atomIndex].element())
		return;
	
	// Check if new element is already in the system
	int newElementNumber = -1;
	for (i = 0; i < _atoms.length(); ++i)
	{
		if (element == _atoms[i][0].element())
		{
			newElementNumber = i;
			break;
		}
	}
	
	// Adding a new element
	if (newElementNumber == -1)
	{
		newElementNumber = _atoms.length();
		_atoms.add();
	}
	
	// Save element
	_atoms[newElementNumber].add();
	_atoms[newElementNumber].last() = _atoms[elemIndex][atomIndex];
	_atoms[newElementNumber].last().element(element);
	_atoms[elemIndex].remove(atomIndex);
	
	// No elements left of type
	if (!_atoms[elemIndex].length())
		_atoms.remove(elemIndex);
	
	// Output
	if (!showOutput)
		Output::quietOff();
}



/* List<Atom*> ISO::coordination(const Atom* atom, double fracTol) const
 *
 * Get the coordination for the atoms
 */

List<Atom*> ISO::coordination(const Atom* atom, double fracTol) const
{
	
	// Variable to store result
	List<Atom*> res;
	
	// Return if there are no atoms
	if (!_atoms.length())
		return res;
	
	// Loop over atoms to figure out the nearest to current
	int i, j;
	double curDis;
	double minDis = (atom == &_atoms[0][0]) ? \
		_basis.secondDistance(_atoms[0][0].fractional(), FRACTIONAL, atom->fractional(), FRACTIONAL) :
		_basis.distance(_atoms[0][0].fractional(), FRACTIONAL, atom->fractional(), FRACTIONAL);
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			curDis = (atom == &_atoms[i][j]) ? \
				_basis.secondDistance(_atoms[i][j].fractional(), FRACTIONAL, atom->fractional(), FRACTIONAL) :
				_basis.distance(_atoms[i][j].fractional(), FRACTIONAL, atom->fractional(), FRACTIONAL);
			if (curDis < minDis)
				minDis = curDis;
		}
	}
	
	// Initialize image interator
	ImageIterator images;
	images.setCell(_basis, minDis*(1 + fracTol));
	
	// Loop over atoms in the structure
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			
			// Set iterator for current pair
			images.reset(atom->fractional(), _atoms[i][j].fractional());
			
			// Loop over images
			while (!images.finished())
			{
				
				// Save if distance is not zero
				if (++images > 1e-8)
					res += &_atoms[i][j];
			}
		}
	}
	
	// Return result
	return res;
}



/* OList<Atom>::D2 ISO::shells(const Atom* atom, double maxDistance, double tol) const
 *
 * Get a list of atoms within set range and group into shells
 */

OList<Atom>::D2 ISO::shells(const Atom* atom, double maxDistance, double tol) const
{
	
	// Output
	Output::newline();
	Output::print("Calculating nearest neighbor shells out to ");
	Output::print(maxDistance);
	Output::print(" Ang for atom ");
	Output::print(atom->atomNumber() + 1);
	Output::increase();
	
	// Variable to store result
	OList<Atom>::D2 res;
	List<double>::D2 distances;
	
	// Loop over atoms in the structure
	int i, j, k, m;
	double dis;
	ImageIterator images(_basis, maxDistance);
	for (i = 0; i < _atoms.length(); ++i)
	{
		for (j = 0; j < _atoms[i].length(); ++j)
		{
			
			// Loop over periodic images of current atom
			images.reset(atom->fractional(), _atoms[i][j].fractional());
			while (!images.finished())
			{
				
				// Get current distance
				dis = ++images;
				if (dis < 1e-8)
					continue;
				
				// Loop over shells that are already saved and check if atom belongs in existing shell
				for (k = 0; k < distances.length(); ++k)
				{
					for (m = 0; m < distances[k].length(); ++m)
					{
						
						// Atom belongs in current shell
						if (Num<double>::eq(dis, distances[k][m], tol))
						{
							res[k] += _atoms[i][j];
							res[k].last().fractional(images.fracVector() + atom->fractional(), false);
							distances[k] += dis;
							break;
						}
					}
					if (m < distances[k].length())
						break;
				}
				
				// Atom belongs in a new shell
				if (k >= distances.length())
				{
					
					// Find inseration point
					for (k = 0; k < distances.length(); ++k)
					{
						if (dis < distances[k][0])
							break;
					}
					res.add();
					distances.add();
					
					// Move new shell to correct position
					for (m = distances.length() - 1; m > k; --m)
					{
						res.swap(m, m-1);
						distances.swap(m, m-1);
					}
					
					// Save atom and properties
					res[k] += _atoms[i][j];
					res[k].last().fractional(images.fracVector() + atom->fractional(), false);
					distances[k] += dis;
				}
			}
		}
	}
	
	// Output
	Output::decrease();
	
	// Return result
	return res;
}



/* bool ISO::equivalent(const ISO& compISO, double tol, bool matchVolume, bool matchCellParams) const
 *
 * Return whether two structures are the same
 */

bool ISO::equivalent(const ISO& compISO, double tol, bool matchVolume, bool matchCellParams) const
{
	
	// Output
	Output::newline();
	Output::print("Comparing structures to determine whether they are the same");
	Output::increase();
	
	// Output
	Output::newline();
	Output::print("Checking equality of the elements in the structure");
	Output::increase();
	
	// Make sure that number of elements is the same in both structures
	bool res = true;
	if (_atoms.length() != compISO.atoms().length())
	{
		Output::newline();
		Output::print("Number of elements in the structures are not the same");
		res = false;
	}
	
	// Only check structures if similar to this point
	int i, j;
	if (res)
	{
		
		// Make sure that the elements are the same in both structures
		for (i = 0; i < _atoms.length(); ++i)
		{
			for (j = 0; j < compISO.atoms().length(); ++j)
			{
				if (_atoms[i][0].element() == compISO.atoms()[j][0].element())
					break;
			}
			if (j >= compISO.atoms().length())
				break;
		}
		if (i < _atoms.length())
		{
			Output::newline();
			Output::print("Elements are not the same in the structures");
			res = false;
		}
	}
	
	// Output
	Output::decrease();
	
	// Only check structures if similar to this point
	ISO thisPrim;
	ISO compPrim;
	Matrix3D thisToPrim = 0.0;
	Matrix3D compToPrim = 0.0;
	if (res)
	{
		
		// Output
		Output::newline();
		Output::print("Converting structures to primitive form for further comparisons");
		Output::increase();
		
		// Get current structure in primitive form
		thisPrim = *this;
		thisToPrim = thisPrim.primitiveTransformation(tol);
		thisPrim.transform(thisToPrim, tol);
	
		// Get comparison structure in primitive form
		compPrim = compISO;
		compToPrim = compPrim.primitiveTransformation(tol);
		compPrim.transform(compToPrim, tol);
		
		// Make sure that the number of atoms is the same
		if (thisPrim.numAtoms() != compPrim.numAtoms())
		{
			Output::newline();
			Output::print("Number of atoms in primitive cells are not the same");
			res = false;
		}
		
		// Output
		Output::decrease();
	}
	
	// Check that the number of atoms of each element is the same
	if (res)
	{
		
		// Output
		Output::newline();
		Output::print("Checking the number of atoms of each element");
		Output::increase();
		
		// Check that the number of atoms are the same
		bool found = true;
		for (i = 0; i < thisPrim.atoms().length(); ++i)
		{
			found = false;
			for (j = 0; j < compPrim.atoms().length(); ++j)
			{
				if (thisPrim.atoms()[i][0].element() == compPrim.atoms()[j][0].element())
				{
					if (thisPrim.atoms()[i].length() == compPrim.atoms()[j].length())
						found = true;
					break;
				}
			}
			if (!found)
				break;
		}
		if (!found)
		{
			Output::newline();
			Output::print("Number of atoms of each element in the primitive cells do not match");
			res = false;
		}
		
		// Output
		Output::decrease();
	}
	
	// Check volume if needed
	if ((res) && (matchVolume))
	{
		
		// Output
		Output::newline();
		Output::print("Comparing the volumes of the primitive cells");
		Output::increase();
		
		// Get the volume tolerance for current structure
		double volTol = 0;
		for (i = 0; i < 3; ++i)
			volTol += pow(tol / thisPrim.basis().lengths()[0], 2);
		volTol = thisPrim.basis().volume() * sqrt(volTol);
		
		// Check if volumes are the same
		if (Num<double>::abs(thisPrim.basis().volume() - compPrim.basis().volume()) > volTol)
		{
			Output::newline();
			Output::print("Volumes are not the same");
			res = false;
		}
		
		// Output
		Output::decrease();
	}
	
	// Check the cell parameters
	Matrix3D thisPrimToRed = 0.0;
	Matrix3D compPrimToRed = 0.0;
	if (res)
	{
		
		// Output
		if (matchCellParams)
		{
			Output::newline();
			Output::print("Comparing basis metrics");
			Output::increase();
		}
		
		// Output
		Output::newline();
		Output::print("Converting primitive cells to reduced form");
		Output::increase();
		
		// Set the tolerances
		double lenTol = tol / pow(thisPrim.basis().volume(), 1.0/3.0);
		double angTol = atan(lenTol);
		
		// Make cells reduced
		thisPrimToRed = Basis::reducedTransformation(thisPrim.basis().vectors());
		thisPrim.transform(thisPrimToRed, tol);
		compPrimToRed = Basis::reducedTransformation(compPrim.basis().vectors());
		compPrim.transform(compPrimToRed, tol);
		
		// Output
		Output::decrease();
		
		// Only run this section if comparing the cell parameters
		if (matchCellParams)
		{
			
			// Check the current angles
			const Vector3D& thisAngles = thisPrim.basis().angles();
			const Vector3D& compAngles = compPrim.basis().angles();
			if ((Num<double>::abs(thisAngles[0] - compAngles[0]) > angTol) || \
				(Num<double>::abs(thisAngles[1] - compAngles[1]) > angTol) || \
				(Num<double>::abs(thisAngles[2] - compAngles[2]) > angTol))
			{
			
				// Output
				Output::newline();
				Output::print("Attempting transformation between reduced cell types I and II");
				Output::increase();
			
				// Try converting cell type for comparison cell
				Matrix3D changeType (-1, 0, 0, 0, -1, 0, 0, 0, 1);
				compPrim.transform(changeType, tol);
				compPrimToRed *= changeType;
			
				// Output
				Output::decrease();
			
				// Check new angles
				if ((Num<double>::abs(thisAngles[0] - compAngles[0]) > angTol) || \
					(Num<double>::abs(thisAngles[1] - compAngles[1]) > angTol) || \
					(Num<double>::abs(thisAngles[2] - compAngles[2]) > angTol))
				{
					Output::newline();
					Output::print("Angles in primitive cell do not match");
					res = false;
				}
			}
		
			// Compare lattice parameter ratios
			if (res)
			{
				const Vector3D& thisLengths = thisPrim.basis().lengths();
				const Vector3D& compLengths = compPrim.basis().lengths();
				if ((Num<double>::abs(thisLengths[0]/thisLengths[1] - compLengths[0]/compLengths[1]) > lenTol) || \
					(Num<double>::abs(thisLengths[0]/thisLengths[2] - compLengths[0]/compLengths[2]) > lenTol) || \
					(Num<double>::abs(thisLengths[1]/thisLengths[2] - compLengths[1]/compLengths[2]) > lenTol))
				{
					Output::newline();
					Output::print("Ratios of primitive cell lattice vector lengths do not match");
					res = false;
				}
			}
		
			// Output
			Output::decrease();
		}
	}
	
	// Check atom positions
	Vector3D curVec = 0.0;
	Matrix3D compRedToThisRed = 0.0;
	if (res)
	{
		
		// Output
		Output::newline();
		Output::print("Comparing positions of atoms by searching for a rotation and vector that overlap structures");
		Output::increase();
		
		// Make copy of comparison atoms and order the same as current atoms
		Atoms compAtoms = compPrim.atoms();
		for (i = 0; i < thisPrim.atoms().length(); ++i)
		{
			if (thisPrim.atoms()[i][0].element() != compAtoms[i][0].element())
			{
				for (j = i + 1; j < compAtoms.length(); ++j)
				{
					if (thisPrim.atoms()[i][0].element() == compAtoms[j][0].element())
					{
						compAtoms.swap(i, j);
						break;
					}
				}
			}
		}
		
		// Figure out which position permutations to test to account for lattice symmetry
		Linked<Matrix3D> rotations1;
		Linked<Matrix3D> rotations2;
		Basis::getPossibleRotations(rotations1, thisPrim.basis().vectors(), tol);
		Basis::getPossibleRotations(rotations2, compPrim.basis().vectors(), tol);
		
		// Build full list of rotations
		Linked<Matrix3D> rotations;
		fillRotations(rotations, (rotations1.length() > rotations2.length()) ? rotations1 : rotations2);
		
		// Output
		Output::newline();
		Output::print("Testing ");
		Output::print(rotations.length());
		Output::print(" possible rotation");
		if (rotations.length() != 1)
			Output::print("s");
		Output::print(" based on the lattice symmetry");
		
		// Make a list of distances of all atoms from the origin
		Vector3D origin = 0.0;
		List<double>::D2 distances;
		distances.length(thisPrim.atoms().length());
		for (i = 0; i < thisPrim.atoms().length(); ++i)
		{
			distances[i].length(thisPrim.atoms()[i].length());
			for (j = 0; j < thisPrim.atoms()[i].length(); ++j)
				distances[i][j] = thisPrim.basis().distance(thisPrim.atoms()[i][j].fractional(), FRACTIONAL, \
					origin, FRACTIONAL);
		}
		
		// Figure out which element occurs least in the structure
		int minElem = 0;
		for (i = 1; i < thisPrim.atoms().length(); ++i)
		{
			if (thisPrim.atoms()[i].length() < thisPrim.atoms()[minElem].length())
				minElem = i;
		}
		
		// Loop over all conversions
		int k;
		Vector3D newPos = 0.0;
		Atoms rotAtoms = compAtoms;
		Atoms transAtoms = compAtoms;
		Linked<Matrix3D>::iterator itRot = rotations.begin();
		for (; itRot != rotations.end(); ++itRot)
		{
			
			// Apply conversion to all positions
			for (i = 0; i < compAtoms.length(); ++i)
			{
				for (j = 0; j < compAtoms[i].length(); ++j)
				{
					newPos  = compAtoms[i][j].fractional();
					newPos *= *itRot;
					ISO::moveIntoCell(newPos);
					rotAtoms[i][j].fractional(newPos);
				}
			}
			
			// Loop over atoms and try translations
			for (i = 0; i < rotAtoms[minElem].length(); ++i)
			{
			
				// Get vector connecting atoms and translate all atoms
				curVec = thisPrim.atoms()[minElem][0].fractional() - rotAtoms[minElem][i].fractional();
				for (j = 0; j < rotAtoms.length(); ++j)
				{
					for (k = 0; k < rotAtoms[j].length(); ++k)
					{
						newPos  = rotAtoms[j][k].fractional();
						newPos += curVec;
						ISO::moveIntoCell(newPos);
						transAtoms[j][k].fractional(newPos);
					}
				}
			
				// Check if atoms overlap
				if (areSitesEqual(thisPrim.basis(), thisPrim.atoms(), transAtoms, tol, &curVec, &distances))
					break;
			}
		
			// Match was made
			if (i < compAtoms[minElem].length())
			{
				compRedToThisRed = *itRot;
				break;
			}
		}
		if (itRot == rotations.end())
		{
			Output::newline();
			Output::print("Could not find a rotation and vector to overlap atomic positions");
			res = false;
		}
		
		// Output
		Output::decrease();
	}
	
	// Output
	Output::newline();
	Output::print("Structures are");
	if (!res)
		Output::print(" not");
	Output::print(" the same");
	
	// Print conversion if needed
	if (res)
	{
		
		// Create transformation
		Matrix3D compToOrig = Matrix3D::identity();
		compToOrig *= compToPrim;
		compToOrig *= compPrimToRed;
		compToOrig *= thisPrimToRed.inverse();
		compToOrig *= thisToPrim.inverse();
		
		// Create rotation
		Matrix3D rotation = compRedToThisRed;
		rotation = thisPrimToRed.transpose() * rotation * thisPrimToRed.transpose().inverse();
		rotation = thisToPrim.transpose() * rotation * thisToPrim.transpose().inverse();
		
		// Create translation
		curVec *= thisPrimToRed.transpose();
		curVec *= thisToPrim.transpose();
		moveIntoCell(curVec);
		
		// Print transformation
		Output::increase();
		Output::newline();
		Output::print("Transformation to take second cell to first");
		Output::newline();
		Output::print("    Transform: [");
		for (i = 0; i < 3; ++i)
		{
			if (i)
				Output::print(", ");
			Output::print("(");
			for (j = 0; j < 3; ++j)
			{
				if (j)
					Output::print(", ");
				Output::print(compToOrig(i, j));
			}
			Output::print(")");
		}
		Output::print("]");
		Output::newline();
		Output::print("       Rotate: [");
		for (i = 0; i < 3; ++i)
		{
			if (i)
				Output::print(", ");
			Output::print("(");
			for (j = 0; j < 3; ++j)
			{
				if (j)
					Output::print(", ");
				Output::print(rotation(i, j));
			}
			Output::print(")");
		}
		Output::print("]");
		Output::newline();
		Output::print("        Shift: (");
		for (i = 0; i < 3; ++i)
		{
			if (i)
				Output::print(", ");
			Output::print(curVec[i]);
		}
		Output::print(")");
		Output::decrease();
	}
	
	// Output
	Output::decrease();
	
	// Return the result
	return res;
}



/* bool ISO::areSitesEqual(const Basis& basis, const Atoms& origAtoms, const Atoms& newAtoms, double tol,
 *		Vector3D* vector, List<double>::D2* origDistances)
 *
 * Return whether two sets of sites are the same
 */

bool ISO::areSitesEqual(const Basis& basis, const Atoms& origAtoms, const Atoms& newAtoms, double tol, \
	Vector3D* vector, List<double>::D2* origDistances)
{
	
	// Set break tolerance
    double breakTol = 4 * tol;
	
	// Sort elements in original structure by how often they appear, fewest first
	int i, j;
	List<int> order;
	order.length(origAtoms.length());
	for (i = 0; i < origAtoms.length(); ++i)
	{
		order[i] = i;
		for (j = i; j > 0; --j)
		{
			if (origAtoms[order[j]].length() < origAtoms[order[j-1]].length())
				order.swap(j, j-1);
			else
				break;
		}
	}
	
	// Variable to store nearest pairs of atoms
	Linked<Atom*> pairOrig;
	Linked<Atom*> pairNew;
	Linked<Vector3D > error;
	
	// Loop over elements
	int elem;
	double curDis;
	double nearDis;
	double origDis;
	Linked<Atom*>::iterator itNew;
	Linked<Atom*>::iterator itNewNear;
	Linked<double>::iterator itDis;
	Linked<double>::iterator itDisNear;
	Vector3D curCell;
	Vector3D nearCell;
	Vector3D origin = 0.0;
	for (i = 0; i < order.length(); ++i)
	{
		
		// Make linked lists of new atoms and their distances from origin
		elem = order[i];
		Linked<Atom*> newList;
		Linked<double> newDis;
		for (j = 0; j < newAtoms[elem].length(); ++j)
		{
			newList.add(&(newAtoms[elem][j]));
			newDis.add(basis.distance(newAtoms[elem][j].fractional(), FRACTIONAL, origin, FRACTIONAL));
		}
		
		// Loop over original atoms of current element
		for (j = 0; j < origAtoms[elem].length(); ++j)
		{
			
			// Get distance from origin for current original atom
			if (origDistances)
				origDis = (*origDistances)[elem][j];
			else
				origDis = basis.distance(origAtoms[elem][j].fractional(), FRACTIONAL, origin, FRACTIONAL);
			
			// Loop over current atoms to find one that could be nearest
			for (itNew = newList.begin(), itDis = newDis.begin(); itNew != newList.end(); ++itNew, ++itDis)
			{
				if (Num<double>::abs(origDis - *itDis) <= breakTol)
				{
					if (origAtoms[elem][j].equal(**itNew, breakTol, &nearDis, &nearCell))
					{
						itNewNear = itNew;
						itDisNear = itDis;
						break;
					}
				}
			}
			
			// Did not find a possible match
			if (itNew == newList.end())
				return false;
			
			// Loop over atoms
			for (; itNew != newList.end(); ++itNew, ++itDis)
			{
				
				// Skip if atoms cannot be the same based on distance from the origin
				if (Num<double>::abs(origDis - *itDis) > breakTol)
					continue;
				
				// Check if atoms are the same
				if (!origAtoms[elem][j].equal(*(*itNew), breakTol, &curDis, &curCell))
					continue;
				
				// Found a new minimum
				if (curDis < nearDis)
				{
					
					// Save properties
					itNewNear = itNew;
					itDisNear = itDis;
					nearDis = curDis;
					nearCell = curCell;
					
					// Found a match so break
					if (nearDis < tol/10)
						break;
				}
			}
			
			// Save the pair
			pairOrig.add(&(origAtoms[elem][j]));
			pairNew.add(*itNewNear);
			
			// Save the error
			error.add(origAtoms[elem][j].fractional() - ((*itNewNear)->fractional() + nearCell));
			
			// Remove the matched atom
			newList.remove(itNewNear);
			newDis.remove(itDisNear);
		}
	}
	
	// Get the average error
	Vector3D avgError = 0.0;
	Linked<Vector3D >::iterator it;
	for (it = error.begin(); it != error.end(); ++it)
		avgError += *it;
	avgError /= error.length();

	// Loop over all pairs of atoms and check if modified distance is within error
	Vector3D tempPos;
	Linked<Atom*>::iterator itPairNew = pairNew.begin();
	Linked<Atom*>::iterator itPairOrig = pairOrig.begin();
	for (; itPairOrig != pairOrig.end(); ++itPairOrig, ++itPairNew)
	{
		
		// Get the position of the error modified atom
		tempPos = (*itPairNew)->fractional() + avgError;
		ISO::moveIntoCell(tempPos);
		
		// Atoms are not within tolerance
		if (basis.distance((*itPairOrig)->fractional(), FRACTIONAL, tempPos, FRACTIONAL) > tol)
			return false;
	}
	
	// Update vector
	if (vector)
		*vector += avgError;
	
	// Return that lists are the same
	return true;
}



/* LatticePoints ISO::getLatticePoints(const Matrix3D& transformationToNewCell)
 *
 * Get the lattice points of an original basis in a new basis obtained from the transformation
 */

LatticePoints ISO::getLatticePoints(const Matrix3D& transformationToNewCell)
{
	
	// Tolerance used in this section
	double tol = 1e-6;
	
	// Get the inverse transpose of the transformation
	Matrix3D transpose (transformationToNewCell.transpose());
	Matrix3D invTrans (transpose.inverse());
	
	// Loop over lattice directions
	int i, j, k;
	bool isNew;
	Vector3D curVec;
	Vector3D transVec;
	LatticePoints latVecs[3];
	LatticePoints transVecs;
	for (i = 0; i < 3; ++i)
	{
		
		// Loop until points are repeated
		isNew = true;
		curVec = 0.0;
		transVecs = curVec;
		latVecs[i] = curVec;
		while (isNew)
		{
			
			// Generate new position
			++curVec[i];
			transVec = invTrans * curVec;
			ISO::moveIntoCell(transVec);
			
			// Check if position is known
			for (j = 0; j < transVecs.length(); ++j)
			{
				isNew = false;
				for (k = 0; k < 3; ++k)
				{
					if (!Num<double>::modEq(transVec[k], transVecs[j][k], tol))
					{
						isNew = true;
						break;
					}
				}
				if (!isNew)
					break;
			}
			
			// Position is new
			if (isNew)
			{
				latVecs[i] += curVec;
				transVecs += transVec;
			}
		}
	}
	
	// Loop over permutations of vectors
	int m, n;
	transVecs.clear();
	for (i = 0; i < latVecs[0].length(); ++i)
	{
		for (j = 0; j < latVecs[1].length(); ++j)
		{
			for (k = 0; k < latVecs[2].length(); ++k)
			{
				
				// Save current vector
				curVec = latVecs[0][i];
				curVec += latVecs[1][j];
				curVec += latVecs[2][k];
				
				// Get vector in cell
				transVec  = curVec;
				transVec *= invTrans;
				ISO::moveIntoCell(transVec);
				
				// Check if position is known
				isNew = true;
				for (m = 0; m < transVecs.length(); ++m)
				{
					isNew = false;
					for (n = 0; n < 3; ++n)
					{
						if (!Num<double>::modEq(transVec[n], transVecs[m][n], tol))
						{
							isNew = true;
							break;
						}
					}
					if (!isNew)
						break;
				}
				
				// Found a new position
				if (isNew)
					transVecs += transVec;
			}
		}
	}
	
	// Transform points back to original cell
	LatticePoints res(transVecs);
	for (i = 0; i < res.length(); ++i)
		res[i] *= transpose;
	
	// Return result
	return res;
}



/* Word ISO::system(LatticeSystem input)
 *
 * Return the name of a lattice system
 */

Word ISO::system(LatticeSystem input)
{
	switch (input)
	{
		case LS_TRICLINIC:
			return "Triclinic";
		case LS_MONOCLINIC:
			return "Monoclinic";
		case LS_ORTHORHOMBIC:
			return "Orthorhombic";
		case LS_TETRAGONAL:
			return "Tetragonal";
		case LS_RHOMBOHEDRAL:
			return "Rhombohedral";
		case LS_HEXAGONAL:
			return "Hexagonal";
		case LS_CUBIC:
			return "Cubic";
		default:
			return "Unknown";
	}
}



/* Word ISO::centering(LatticeCentering input)
 *
 * Return the name of a lattice centering
 */

Word ISO::centering(LatticeCentering input)
{
	switch (input)
	{
		case LC_P:
			return "Primitive";
		case LC_H:
			return "Primitive";
		case LC_C:
			return "C face";
		case LC_A:
			return "A face";
		case LC_B:
			return "B face";
		case LC_I:
			return "Body";
		case LC_R:
			return "Rhombohedral";
		case LC_F:
			return "All faces";
		default:
			return "Unknown";
	}
}



/* LatticeSystem ISO::system(const Word& input)
 *
 * Return the lattice system given by a word
 */

LatticeSystem ISO::system(const Word& input)
{
	if (input.equal("triclinic", false, 3))
		return LS_TRICLINIC;
	if (input.equal("monoclinic", false, 3))
		return LS_MONOCLINIC;
	if (input.equal("orthorhombic", false, 4))
		return LS_ORTHORHOMBIC;
	if (input.equal("tetragonal", false, 3))
		return LS_TETRAGONAL;
	if (input.equal("rhombohedral", false, 4))
		return LS_RHOMBOHEDRAL;
	if (input.equal("hexagonal", false, 3))
		return LS_HEXAGONAL;
	if (input.equal("cubic", false, 3))
		return LS_CUBIC;
	return LS_UNKNOWN;
}
