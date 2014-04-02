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



#ifndef NUM_H
#define NUM_H



#include "list.h"
#include "constants.h"
#include <cmath>



// Numerical integration methods
enum IntegrationMethod {GAUSSLEGENDRE};

// Numerical root finding methods
enum RootMethod {BRENT, BISECTION};



// Complex numbers
struct Complex
{
	double real;
	double imag;
};


// Rename blas and lapack functions if not using mkl
#ifndef MINT_MKL
#define DCOPY  dcopy_
#define ZCOPY  zcopy_
#define DSCAL  dscal_
#define DAXPY  daxpy_
#define DDOT   ddot_
#define DNRM2  dnrm2_
#define DGEMV  dgemv_
#define DGEMM  dgemm_
#define DGETRF dgetrf_
#define DGESV  dgesv_
#define DGESVD dgesvd_
#define DGEEV  dgeev_
#define DSYEV  dsyev_
#define DSYEVD dsyevd_
#define ZGEEV  zgeev_
#define ZHEEV  zheev_
#define ZHEEVD zheevd_
#endif // MINT_MKL


// BLAS and LAPACK
extern "C"
{
	
	// Copy array to another
	void   DCOPY  (const int* length, double* origVec, const int* incOrig, double* newVec, const int* incNew);
	
	// Copy complex array to another
	void   ZCOPY  (const int* length, Complex* origVec, const int* incOrig, Complex* newVec, const int* incNew);
	
	// Compute x = ax where a is a constant
	void   DSCAL  (const int* length, double* a, double* vector, const int* inc);
	
	// Compute y = ax + y where a is a constant
	void   DAXPY  (const int* length, double* a, double* vector, const int* incVec, double* res, const int* incRes);
	
	// Dot product
	double DDOT   (const int* length, double* vec1, const int* inc1, double* vec2, const int* inc2);
	
	// Euclidean norm
	double DNRM2  (const int* length, double* vector, const int* inc);
	
	// Matrix multiplied onto vector
	void   DGEMV  (char* trans, const int* numRows, const int* numCols, double* alpha, double* matrix, \
				   const int* ldMat, const double* vector, int* vecInc, double* beta, double* res, int* resInc);

	// Matrix multiplied onto matrix
	void   DGEMM  (char* transA, char* transB, const int* numRowsA, const int* numColsB, const int* numColsA, \
				   double* alpha, double* matA, const int* ldA, double* matB, const int* ldB, double* beta, \
				   double* res, const int* ldRes);

	// LU decomposition
	void   DGETRF (const int* numRows, const int* numCols, double* matrix, const int* ldA, int* pivots, int* info);
	
	// Solve Mx = y
	void   DGESV  (const int *numMRows, const int* numYCols, double* matrix, const int* ldA, int* pivots, \
				   double* rhs, const int* ldrhs, int* info);
		
	// Singular value decomposition
	void   DGESVD (const char* jobU, const char* jobVT, const int* numRowsA, const int* numColsA, double* A, \
				   const int* ldA, double* singVals, double* U, const int* ldU, double* VT, const int* ldVT, \
				   double* work, const int* lenWork, int* info);
	
	// Eigenvalues/eigenvectors for real general matrix
	void   DGEEV  (const char* jobVL, const char* jobVR, const int* sizeA, double* A, const int* ldA, \
				   double* realEVals, double* complexEVals, double* leftEVecs, const int* ldLEV, \
				   double* rightEVecs, const int* ldREV, double* work, const int* ldWork, int* info);
	
	// Eigenvalues/eigenvectors of real symmetric matrix
	void   DSYEV  (const char* jobZ, const char* upOrLow, const int* sizeA, double* A, const int* ldA, \
				   double* evals, double* work, const int* ldWork, int* info);
	
	// Eigenvalues/eigenvectors of real symmetric matrix using divide and conquer
	void   DSYEVD (const char* jobZ, const char* upOrLow, const int* sizeA, double* A, const int* ldA, \
				   double* evals, double* work, const int* ldWork, int* iWork, const int* ldIWork, int* info);
	
	// Complex eigenvalues/eigenvectors for general matrix
	void   ZGEEV  (const char* jobVL, const char* jobVR, const int* sizeA, Complex* A, const int* ldA, \
				   Complex* evals, Complex* leftEVecs, const int* ldLEV, Complex* rightEVecs, const int* ldREV, \
				   Complex* lWork, const int* ldLWork, double* rWork, int* info);
	
	// Complex eigenvalues/eigenvectors for Hermitian matrix
	void   ZHEEV  (const char* jobZ, const char* upOrLow, const int* sizeA, Complex* A, const int* ldA, \
				   double* evals, Complex* work, const int* lenWork, double* rWork, int* info);
	
	// Complex eigenvalues/eigenvectors for Hermitian matrix using divide and conquer
	void   ZHEEVD (const char* jobZ, const char* upOrLow, const int* sizeA, Complex* A, const int* ldA, \
				   double* evals, Complex* work, const int* lenWork, double* rWork, const int* lenRWork, \
				   int* iWork, const int* lenIWork, int* info);
};



// Declare that matrix class will be defined later
class Matrix;
class Matrix3D;



// Class to store a vector
class Vector
{
	
	// Variables
	int _length;
	int _actualLength;
 	double* _vector;
	
	// Functions
	void initialize();
	
public:
	
	// Constructors
	Vector();
	Vector(const Vector& copy);
	Vector(int inLength);
	Vector(int inLength, double value);
	
	// Destructor
	~Vector();
	
	// Delete all data
	void clear();
	
	// Setup functions
	void length (int inLength);
	void fill   (double input);
	void append (double input);
	
	// General functions
	bool operator== (const Vector& rhs) const;
	
	// Assignment functions
	Vector& operator= (double rhs);
	Vector& operator= (const Vector& rhs);
	
	// Math operators
	double  magnitude  ()                  const;
	Vector  operator*  (double rhs)        const;
	Vector  operator/  (double rhs)        const;
	double  operator*  (const Vector& rhs) const;
	Vector  operator+  (const Vector& rhs) const;
	Vector  operator-  (const Vector& rhs) const;
	Vector& operator+= (const Vector& rhs);
	Vector& operator-= (const Vector& rhs);
	Vector& operator*= (const Matrix& rhs);
	Vector& operator*= (double rhs);
	Vector& operator/= (double rhs);
	
	// Access functions
	int length     				()          const	{ return _length; }
	double*		  	operator() 	()					{ return _vector; }
	const double*	operator() 	()			const	{ return _vector; }
	double&      	operator[]	(int index)			{ return _vector[index]; }
	const double&	operator[]	(int index) const	{ return _vector[index]; }
	
	// Friends
	friend class Matrix;
	friend class Multi;
};



// Complex vector
class CVector
{

	// Variables
	int _length;
	int _actualLength;
 	Complex* _vector;
	
	// Functions
	void initialize();

public:
	
	// Constructors
	CVector();
	CVector(const CVector& copy);
	CVector(int inLength);
	
	// Destructor
	~CVector();
	
	// Delete all data
	void clear();
	
	// Setup functions
	void length (int inLength);
	
	// Assignment functions
	CVector& operator= (const CVector& rhs);
	
	// Access functions
	int length     				()          const	{ return _length; }
	Complex*	    operator() 	()					{ return _vector; }
	const Complex*	operator() 	()			const	{ return _vector; }
	Complex&      	operator[]	(int index)			{ return _vector[index]; }
	const Complex&	operator[]	(int index) const	{ return _vector[index]; }
	
	// Friends
	friend class CMatrix;
};



// Class to store a vector of fixed length 3
class Vector3D
{
	
	// Variables
 	double _vector[3];
	
public:
	
	// Constructors
	Vector3D();
	Vector3D(double value);
	Vector3D(const double* vector);
	Vector3D(double val1, double val2, double val3);
	Vector3D(const Vector3D& copy);
	
	// Setup functions
	void fill (double input);
	void set  (double val1, double val2, double val3);
	void set  (const double* vector);
	
	// General functions
	bool operator== (const Vector3D& rhs) const;
	bool operator!= (const Vector3D& rhs) const;
	
	// Assignment functions
	Vector3D& operator= (double value);
	Vector3D& operator= (const double* vector);
	Vector3D& operator= (const Vector3D& rhs);
	
	// Math operators
	double    magnitude  ()                      const;
	double    operator*  (const Vector3D& rhs)   const;
	double    angle      (const Vector3D& vec2)  const;
	Vector3D  operator+  (const Vector3D& rhs)   const;
	Vector3D  operator-  (const Vector3D& rhs)   const;
	Vector3D  operator*  (double rhs)            const;
	Vector3D  operator/  (double rhs)            const;
	Vector3D& operator+= (const Vector3D& rhs);
	Vector3D& operator+= (const double* rhs);
	Vector3D& operator-= (const Vector3D& rhs);
	Vector3D& operator*= (const Matrix3D& rhs);
	Vector3D& operator*= (double rhs);
	Vector3D& operator/= (double rhs);
	Vector3D  cross      (const Vector3D& rhs)   const;
	
	// Geometry functions
	double distanceToPoint(const Vector3D& rhs) const;
	double distanceToLineSegment(const Vector3D& begin, const Vector3D& end) const;
	double nearestPointOnPlane(const Vector3D& pointOnPlane, const Vector3D& norm, Vector3D* nearPoint = 0) const;
	double nearestPointOnPlane(const Vector3D& direction, const Vector3D& pointOnPlane, const Vector3D& norm, \
		Vector3D* nearPoint = 0) const;
	void voronoi(const OList<Vector3D>& points, double tol, OList<Vector3D>* vertices = 0, \
		OList<Vector3D>* centers = 0, OList<Vector3D>* norms = 0, List<int>* indices = 0) const;
	void voronoi(const OList<Vector3D>& points, const List<double>& weights, double tol, \
		OList<Vector3D>* vertices = 0, OList<Vector3D>* centers = 0, OList<Vector3D>* norms = 0, \
		List<int>* indices = 0) const;
	
	// Access functions
	double*		  operator() ()					{ return _vector; }
	const double* operator() ()			 const	{ return _vector; }
	double& 	  operator[] (int index)		{ return _vector[index]; }
	const double& operator[] (int index) const	{ return _vector[index]; }
	
	// Friends
	friend class Matrix3D;
	friend class Multi;
};



// Class to store a matrix
class Matrix
{
	
	// Variables
	int _numRows;
	int _numCols;
	int _length;
	double* _matrix;
	
	// Functions
	void initialize();
	
public:
	
	// Constructors
	Matrix();
	Matrix(const Matrix& copy);
	Matrix(int inSize);
	Matrix(int inRows, int inCols);
	
	// Destructor
	~Matrix();
	
	// Delete all data
	void clear();
	
	// Setup functions
	void size         (int inSize, bool preserve = false);
	void size         (int inRows, int inCols, bool preserve = false);
	void fill         (double input);
	void swapRows     (int row1, int row2);
	void swapColumns  (int col1, int col2);
	void makeIdentity ();
	
	// General functions
	bool operator== (const Matrix& rhs) const;
	
	// Assignment functions
	Matrix& operator= (const Matrix& rhs);
	
	// Math operators
	bool    isDiagonal    ()                  const;
	double  determinant   ()                  const;
	Vector  operator*     (const Vector& rhs) const;
	Matrix  operator*     (const Matrix& rhs) const;
	Matrix  transpose     ()                  const;
	Matrix& makeTranspose ();
	Matrix  inverse       ()                  const;
	Vector  solve         (const Vector& rhs) const;
	Vector  eigenvalues   (Matrix* eigenvectors = 0, bool isSymmetric = false) const;
	
	// Reductions
	Matrix rowEchelon(int* swaps = 0, Matrix* operations = 0, bool integer = false) const;
	
	// Access functions
	int           numRows ()                    const	{ return _numRows; }
	int           numCols ()                    const	{ return _numCols; }
	const double& operator() (int row, int col) const	{ return _matrix[_numRows*col + row]; }
	double&       operator() (int row, int col)			{ return _matrix[_numRows*col + row]; }
	
	// Static member functions
	static Matrix identity(int size);
	static Matrix identity(int numRows, int numCols);
	
	// Friends
	friend class Vector;
};



// Class to store a complex matrix
class CMatrix
{
	
	// Variables
	int _numRows;
	int _numCols;
	int _length;
	Complex* _matrix;
	
	// Functions
	void initialize();
	
public:
	
	// Constructors
	CMatrix();
	CMatrix(const CMatrix& copy);
	CMatrix(int inSize);
	CMatrix(int inRows, int inCols);
	
	// Destructor
	~CMatrix();
	
	// Delete all data
	void clear();
	
	// Setup functions
	void size        (int inSize, bool preserve = false);
	void size        (int inRows, int inCols, bool preserve = false);
	void swapRows    (int row1, int row2);
	void swapColumns (int col1, int col2);
	
	// Assignment functions
	CMatrix& operator= (const CMatrix& rhs);
	
	// Functions
	CVector eigenvalues(CMatrix* eigenvectors = 0, bool isHermitian = false) const;
	
	// Access functions
	int            numRows ()                    const	{ return _numRows; }
	int            numCols ()                    const	{ return _numCols; }
	const Complex& operator() (int row, int col) const	{ return _matrix[_numRows*col + row]; }
	Complex&       operator() (int row, int col)		{ return _matrix[_numRows*col + row]; }
};



// Class to store a 3D matrix
class Matrix3D
{
	
	// Variables
	double _matrix[9];
	
public:
	
	// Constructors
	Matrix3D();
	Matrix3D(const Matrix3D& copy);
	Matrix3D(double value);
	Matrix3D(double val11, double val12, double val13, double val21, double val22, double val23, double val31, \
			 double val32, double val33);
	
	// Setup functions
	void fill         (double input);
	void set          (double val11, double val12, double val13, double val21, double val22, double val23, \
					   double val31, double val32, double val33);
	void setRow       (int row, const double* vector);
	void setRow       (int row, const Vector3D& vector);
	void swapRows     (int row1, int row2);
	void makeIdentity ();
	
	// General functions
	bool operator== (const Matrix3D& rhs) const;
	bool operator!= (const Matrix3D& rhs) const;
	bool isInteger(double tol) const;
	
	// Assignment functions
	Matrix3D& operator= (const Matrix3D& rhs);
	Matrix3D& operator= (double rhs);
	
	// Math operators
	int		  rank          ()                           const;
	double    trace         ()                           const;
	double    determinant   ()                           const;
	double    volume        ()                           const;
	Matrix3D& operator+=    (const Matrix3D& rhs);
	Matrix3D& operator-=    (const Matrix3D& rhs);
	Vector3D  operator*     (const double* rhs)          const;
	Vector3D  operator*     (const Vector3D& rhs)        const;
	Matrix3D  operator*     (double rhs)                 const;
	Matrix3D  operator*     (const Matrix3D& rhs)        const;
	Matrix3D& operator*=    (const Matrix3D& rhs);
	Matrix3D& operator*=    (double rhs);
	Matrix3D& operator/=    (double rhs);
	Matrix3D  transpose     ()                           const;
	Matrix3D& makeTranspose ();
	Matrix3D  inverse       ()                           const;
	Vector3D  eigenvalues   (Matrix3D* eigenvectors = 0) const;
	
	// Reductions
	Matrix3D rowEchelon(int* swaps = 0, Matrix3D* operations = 0, bool integer = false) const;
	
	// Access functions
	const double& operator() (int row, int col) const	{ return  _matrix[3*row + col]; }
	double&       operator() (int row, int col)			{ return  _matrix[3*row + col]; }
	const double* operator[] (int row)          const   { return &_matrix[3*row]; }
	double*       operator[] (int row)                  { return &_matrix[3*row]; }
	
	// Static member functions
	static Matrix3D identity();
	
	// Friends
	friend class Vector3D;
};



template <class Tclass>
class Functor
{

    const Vector* _params;
    
	// Variables
	Tclass* _objPtr;
	double (*_funPtr)(double);
	double (Tclass::*_objFunPtr)(double);
	double (*_funPtrVec)(const Vector&);
	double (Tclass::*_objFunPtrVec)(const Vector&);
	double (*_funPtrVecArg)(const Vector&, double);
	double (Tclass::*_objFunPtrVecArg)(const Vector&, double);
	
	// Functions
	void initialize();
	
public:
	
	// Constructor
	Functor ();
	Functor (double (*funPtr)(double));
	Functor (Tclass* objPtr, double (Tclass::*objFunPtr)(double));
	Functor (double (*funPtr)(const Vector&));
	Functor (Tclass* objPtr, double (Tclass::*funPtr)(const Vector&));
	Functor (double (*funPtr)(const Vector&, double), const Vector* params = 0);
	Functor (Tclass* objPtr, double (Tclass::*funPtr)(const Vector&, double), const Vector* params = 0);
	
	// Access functions
	double operator() (double arg);
	double operator() (const Vector& params);
	double operator() (const Vector& params, double arg);
};




// VectorFunctor definition
template <class Tclass>
class VectorFunctor
{
	
	// Variables
	Tclass* _objPtr;
	Vector (*_funPtr)            (const Vector&, double);
	Vector (Tclass::*_objFunPtr) (const Vector&, double);
	
	// Functions
	void initialize();
	
public:
	
	// Constructor
	VectorFunctor();
	VectorFunctor(Vector (*funPtr)(const Vector&, double));
	VectorFunctor(Tclass* objPtr, Vector (Tclass::*objFunPtr)(const Vector&, double));
	
	// Access functions
	Vector operator() (const Vector& params, double arg);
};



// Math functions
template <class T>
class Num
{
	
	// Integration point methods
	static void GaussLegendrePoints(T min, T max, Matrix& points);
	static T Legendre(const Matrix& coeffs, int order, T point);
	
public:
	
	// Simple functions
	static T    ceil  (T value);
	static T    floor (T value);
	static T    abs   (T value);
	static T    min   (T val1, T val2);
	static T    max   (T val1, T val2);
	static T    sign  (T value);
	static T    sign  (T value, double tol);
	static T    round (T value, double roundTo);
	static void swap  (T& val1, T& val2);
	static T    next  (T value, T max);
	static T    prev  (T value, T max);
	static T    delta (T x, T y);
	static T    mod   (T value, T max);
	
	// Comparisons
	static bool lt    (T x, T y, double tol)	{ return (x < y - tol); }
    static bool gt    (T x, T y, double tol)	{ return (x > y + tol); }
    static bool le    (T x, T y, double tol)	{ return !(x > y + tol); }
    static bool ge    (T x, T y, double tol)	{ return !(x < y - tol); }
    static bool eq    (T x, T y, double tol)	{ return !((x < y - tol) || (x > y + tol)); }
	static bool neq   (T x, T y, double tol)	{ return !eq(x, y, tol); }
	static bool modEq (T x, T y, double tol)	{ return abs(mod(x - y, 1-tol)) < tol; }
	
	// Convert between radians and degrees
	static T toRadians (T input)				{ return input * Constants::pi / 180; }
	static T toDegrees (T input)				{ return input * 180 / Constants::pi; }
	
	// Convert between eV and Rydberg
	static T energyFromAU (T input)				{ return input * Constants::Ry; }
	static T energyToAU   (T input)				{ return input / Constants::Ry; }
	
	// Convert between bohr and angstrom
	static T distanceFromAU (T input)			{ return input * Constants::bohr; }
	static T distanceToAU   (T input)			{ return input / Constants::bohr; }
	
	// Convert between eV/ang and Ry/bohr
	static T forceFromAU (T input)				{ return input * Constants::Ry / Constants::bohr; }
	static T forceToAU   (T input)				{ return input / Constants::Ry * Constants::bohr; }
	
	// Convert frequency to Hz
	static T frequencyToHz   (T input)			{ return input*Constants::meter*sqrt(Constants::kg/Constants::joule); }
	static T frequencyFromHz (T input)			{ return input/Constants::meter/sqrt(Constants::kg/Constants::joule); }
	
	// Convert frequency to cm^-1
	static T frequencyToCm   (T input)			{ return frequencyToHz(input) / Constants::c / 100; }
	static T frequencyFromCm (T input)			{ return frequencyFromHz(input) * Constants::c * 100; }
	
	// Factors
	static T               gcf           (int length, int* values);
	static List<T>         factors       (int value);
	static OList<List<T> > tripleFactors (int value);
	
	// Fraction functions
	static void decimalTofraction(int* res, T number, double tol);
	
	// Vector functions
	static T dot       (int length, const T* vec1, const T* vec2);
	static T magnitude (int length, const T* vec);
	static T angle     (const T* vec1, const T* vec2);
	static T volume    (const T* vec1, const T* vec2, const T* vec3);
	
	// Integration
	template<class Tclass>
	static T integrate(Functor<Tclass>& fun, T min, T max, IntegrationMethod method = GAUSSLEGENDRE);
};



// Solving functions
template <class Tclass>
class Solve
{
	
public:
	
	// Simplex methods
	static double NelderMead (Functor<Tclass>& fun, double convergeRes, const Vector& initial, const Vector& steps, \
							  Vector& argRes, int mult);
	static void sortVertices (OList<Vector>& vertices, List<double>& values, int left, int right);
	
private:
	
	// Root finding methods
	static double brent     (Functor<Tclass>& fun, double target, double convergeArg, double convergeRes, \
							 double min, double max, double& argRes);
	static double bisection (Functor<Tclass>& fun, double target, double convergeArg, double convergeRes, \
							 double min, double max, double& argRes);
	
	// Static helpers
	static Functor<Tclass>* _functorFun;
	static double NelderHelper(const Vector& args)	{ return (*_functorFun)(args[0]); }
	
public:
	
	// Find roots of a function
	static double findRoot (Functor<Tclass>& fun, double target, double convergeArg, double convergeRes, \
							double min, double max);
	static double findRoot (Functor<Tclass>& fun, double target, double convergeArg, double convergeRes, \
							double min, double max, RootMethod method);
	static double findRoot (Functor<Tclass>& fun, double target, double convergeArg, double convergeRes, \
							double min, double max, double& argRes);
	static double findRoot (Functor<Tclass>& fun, double target, double convergeArg, double convergeRes, \
							double min, double max, RootMethod method, double& argRes);
	
	// Minimize a function of one variable
	static double minimize (Functor<Tclass>& fun, double convergeRes, double initial, double step);
	static double minimize (Functor<Tclass>& fun, double convergeRes, double initial, double step, double& argRes);
	
	// Maximize a function of one variable
	static double maximize (Functor<Tclass>& fun, double convergeRes, double initial, double step);
	static double maximize (Functor<Tclass>& fun, double convergeRes, double initial, double step, double& argRes);
	
	// Minimize a function of multiple variable
	static double minimize (Functor<Tclass>& fun, double convergeRes, const Vector& initial, const Vector& step);
	static double minimize (Functor<Tclass>& fun, double convergeRes, const Vector& initial, const Vector& step, \
							Vector& argRes);
	
	// Maximize a function of multiple variable
	static double maximize (Functor<Tclass>& fun, double convergeRes, const Vector& initial, const Vector& step);
	static double maximize (Functor<Tclass>& fun, double convergeRes, const Vector& initial, const Vector& step, \
							Vector& argRes);
};



// Fitting functions
class Fit
{
	
	// General helper functions
	template <class Tclass>
	static double getRSquared   (const OList<List<double> >& data, Functor<Tclass>& fun, Vector& params);
	
	// Helper functions for Levenberg–Marquardt algorithm
	template <class Tclass>
	static void LMBuildJacobian (const OList<List<double> >& data, VectorFunctor<Tclass>& deriv, Vector& params, \
								 Matrix& jacobian, Matrix& jacobianTranspose);
	static void LMSetMatrix     (Matrix& jacobian, Matrix& jacobianTranspose, double damper, Matrix& matrix);
	template <class Tclass>
	static void LMGetResidual   (const OList<List<double> >& data, Functor<Tclass>& fun, Vector& params, \
								 Vector& residuals);
	
	// Helpers for Nelder-Mead minimization
	template <class Tclass>
	struct NMHelp
	{
		Functor<Tclass>* _functor;
		const OList<List<double> >* _NMData;
		double NMRSquared(const Vector& params);
	};
	
public:
	
	// Polynomial expansion
	static Vector polynomial(const OList<List<double> >& data, int order0, int numTerms);
	
	// Levenberg–Marquardt for general fit
	template <class Tclass>
	static Vector LM (const OList<List<double> >& data, Functor<Tclass>& fun, VectorFunctor<Tclass>& deriv, \
					  const Vector& initial, double tol = 1e-3);
	
	// General fit without derivatives using Nelder-Mead minimization method
	template <class Tclass>
	static Vector NM (const OList<List<double> >& data, Functor<Tclass>& fun, const Vector& initial, \
					  const Vector& steps, double tol = 1e-3);
};



// =====================================================================================================================
// Complex numbers
// =====================================================================================================================

/* inline Complex operator+ (const Complex& lhs, const Complex& rhs)
 *
 * Addition
 */

inline Complex operator+ (const Complex& lhs, const Complex& rhs)
{
	Complex res;
	res.real = lhs.real + rhs.real;
	res.imag = lhs.imag + rhs.imag;
	return res;
}

inline Complex& operator+= (Complex& lhs, const Complex& rhs)
{
	lhs.real += rhs.real;
	lhs.imag += rhs.imag;
	return lhs;
}

inline Complex operator+ (double lhs, const Complex& rhs)
{
	Complex res;
	res.real = lhs + rhs.real;
	res.imag = rhs.imag;
	return res;
}

inline Complex operator+ (const Complex& lhs, double rhs)
{
	Complex res;
	res.real = lhs.real + rhs;
	res.imag = lhs.imag;
	return res;
}

inline Complex& operator+= (Complex& lhs, double rhs)
{
	lhs.real += rhs;
	return lhs;
}



/* inline Complex operator- (const Complex& lhs, const Complex& rhs)
 *
 * Subtraction
 */

inline Complex operator- (const Complex& lhs, const Complex& rhs)
{
	Complex res;
	res.real = lhs.real - rhs.real;
	res.imag = lhs.imag - rhs.imag;
	return res;
}

inline Complex& operator-= (Complex& lhs, const Complex& rhs)
{
	lhs.real -= rhs.real;
	lhs.imag -= rhs.imag;
	return lhs;
}

inline Complex operator- (double lhs, const Complex& rhs)
{
	Complex res;
	res.real = lhs - rhs.real;
	res.imag = rhs.imag;
	return res;
}

inline Complex operator- (const Complex& lhs, double rhs)
{
	Complex res;
	res.real = lhs.real - rhs;
	res.imag = lhs.imag;
	return res;
}

inline Complex& operator-= (Complex& lhs, double rhs)
{
	lhs.real -= rhs;
	return lhs;
}



/* inline Complex operator* (const Complex& lhs, const Complex& rhs)
 *
 * Multiplication
 */

inline Complex operator* (const Complex& lhs, const Complex& rhs)
{
	Complex res;
	res.real = lhs.real*rhs.real - lhs.imag*rhs.imag;
	res.imag = lhs.real*rhs.imag + lhs.imag*rhs.real;
	return res;
}

inline Complex& operator*= (Complex& lhs, const Complex& rhs)
{
	double r = lhs.real;
	double i = lhs.imag;
	lhs.real = r*rhs.real - i*rhs.imag;
	lhs.imag = r*rhs.imag + i*rhs.real;
	return lhs;
}

inline Complex operator* (double lhs, const Complex& rhs)
{
	Complex res;
	res.real = lhs*rhs.real;
	res.imag = lhs*rhs.imag;
	return res;
}

inline Complex operator* (const Complex& lhs, double rhs)
{
	Complex res;
	res.real = lhs.real*rhs;
	res.imag = lhs.imag*rhs;
	return res;
}

inline Complex& operator*= (Complex& lhs, double rhs)
{
	lhs.real *= rhs;
	lhs.imag *= rhs;
	return lhs;
}



/* inline Complex operator/ (const Complex& lhs, const Complex& rhs)
 *
 * Division
 */

inline Complex operator/ (const Complex& lhs, const Complex& rhs)
{
	Complex res;
	double denom = rhs.real*rhs.real + rhs.imag*rhs.imag;
	res.real = (lhs.real*rhs.real + lhs.imag*rhs.imag) / denom;
	res.imag = (lhs.imag*rhs.real - lhs.real*rhs.imag) / denom;
	return res;
}

inline Complex& operator/= (Complex& lhs, const Complex& rhs)
{
	double r = lhs.real;
	double i = lhs.imag;
	double denom = rhs.real*rhs.real + rhs.imag*rhs.imag;
	lhs.real = (r*rhs.real + i*rhs.imag) / denom;
	lhs.imag = (i*rhs.real - r*rhs.imag) / denom;
	return lhs;
}

inline Complex operator/ (double lhs, const Complex& rhs)
{
	Complex res;
	double denom = rhs.real*rhs.real + rhs.imag*rhs.imag;
	res.real =  lhs*rhs.real / denom;
	res.imag = -lhs*rhs.imag / denom;
	return res;
}

inline Complex operator/ (const Complex& lhs, double rhs)
{
	Complex res;
	res.real = lhs.real / rhs;
	res.imag = lhs.imag / rhs;
	return res;
}

inline Complex& operator/= (Complex& lhs, double rhs)
{
	lhs.real /= rhs;
	lhs.imag /= rhs;
	return lhs;
}



// =====================================================================================================================
// Vector
// =====================================================================================================================

/* inline Vector::Vector()
 *
 * Constructors
 */

inline Vector::Vector()
{
	initialize();
}

inline Vector::Vector(const Vector& copy)
{
	initialize();
	*this = copy;
}

inline Vector::Vector(int inLength)
{
	initialize();
	length(inLength);
}
	
inline Vector::Vector(int inLength, double value)
{
	initialize();
	length(inLength);
	fill(value);
}
	


/* inline Vector::~Vector()
 *
 * Destructor for vector
 */

inline Vector::~Vector()
{
	clear();
}



/* inline void Vector::initialize()
 *
 * Initialize vector
 */

inline void Vector::initialize()
{
	_length = 0;
	_vector = 0;
	_actualLength = 0;
}



/* inline void Vector::clear()
 *
 * Clear data in vector
 */

inline void Vector::clear()
{
	if (_vector)
		delete [] _vector;
	initialize();
}



/* inline void Vector::length(int inLength)
 *
 * Set length of vector
 */

inline void Vector::length(int inLength)
{
	
	// Length is already long enough
	if (inLength <= _actualLength)
	{
		_length = inLength;
		return;
	}
	
	// Create new list
	double* temp = new double [inLength];
	
	// Copy old values
	int inc = 1;
	DCOPY(&_length, _vector, &inc, temp, &inc);
	
	// Delete old list
	if (_vector)
		delete [] _vector;
	
	// Save list
	_vector = temp;
	_length = inLength;
	_actualLength = inLength;
}



/* inline void Vector::fill(double input)
 *
 * Fill vector with constant value
 */

inline void Vector::fill(double input)
{
	int i;
	for (i = 0; i < _length-5;)
	{
		_vector[i++] = input;
		_vector[i++] = input;
		_vector[i++] = input;
		_vector[i++] = input;
		_vector[i++] = input;
	}
	for ( ; i < _length; )
		_vector[i++] = input;
}



/* inline bool Vector::operator== (const Vector& rhs) const
 *
 * Test if two vectors are the same
 */

inline bool Vector::operator== (const Vector& rhs) const
{
	
	// Lengths are not the same
	if (_length != rhs._length)
		return false;
	
	// Get the max value in both lists
	int i;
	double cur;
	double scale = Num<double>::max(Num<double>::abs(_vector[0]), Num<double>::abs(rhs._vector[0]));
	for (i = 1; i < _length; ++i)
	{
		cur = Num<double>::max(Num<double>::abs(_vector[i]), Num<double>::abs(rhs._vector[i]));
		scale = (cur > scale) ? cur : scale;
	}
	scale = scale < 1e-14 ? 1 : scale;
	
	// Check values in each list
	for (int i = 0; i < _length; ++i)
	{
		if (Num<double>::abs(_vector[i] - rhs._vector[i]) / scale > 1e-6)
			return false;
	}
	
	// Return that lists are the same if at this point
	return true;
}



/* inline void Vector::append(double input)
 *
 * Append value to vector
 */

inline void Vector::append(double input)
{
	length(_length + 1);
	_vector[_length - 1] = input;
}



/* inline Vector& Vector::operator= (double rhs)
 *
 * Assign vector to value
 */

inline Vector& Vector::operator= (double rhs)
{
	
	// Allocate space
	length(1);
	
	// Save value
	_vector[0] = rhs;
	
	// Return result
	return *this;
}



/* inline Vector& Vector::operator= (const Vector& rhs)
 *
 * Copy vector
 */

inline Vector& Vector::operator= (const Vector& rhs)
{
	
	// Make sure not copying self
	if (this != &rhs)
	{

		// Allocate space
		length(rhs._length);

		// Copy list
		int inc = 1;
		DCOPY(&rhs._length, rhs._vector, &inc, _vector, &inc);
	}
	
	// Return result
	return *this;
}



/* inline double Vector::magnitude() const
 *
 * Return the magnitude of the vector
 */

inline double Vector::magnitude() const
{
	int one = 1;
	return DNRM2(&_length, _vector, &one);
}



/* inline Vector Vector::operator* (double rhs) const
 *
 * Multiply vector by constant
 */

inline Vector Vector::operator* (double rhs) const
{
	int one = 1;
	Vector res(*this);
	DSCAL(&res._length, &rhs, res._vector, &one);	
	return res;
}



/* inline Vector Vector::operator/ (double rhs) const
 *
 * Divide vector by constant
 */

inline Vector Vector::operator/ (double rhs) const
{
	return *this * (1.0/rhs);
}



/* inline double Vector::operator* (const Vector& rhs) const
 *
 * Dot product
 */

inline double Vector::operator* (const Vector& rhs) const
{
	int inc = 1;
	return DDOT(&_length, _vector, &inc, rhs._vector, &inc);
}



/* inline Vector Vector::operator+ (const Vector& rhs) const
 *
 * Add two vectors
 */

inline Vector Vector::operator+ (const Vector& rhs) const
{
	Vector res(*this);
	res += rhs;
	return res;
}



/* inline Vector Vector::operator- (const Vector& rhs) const
 *
 * Subtract one vector from another
 */

inline Vector Vector::operator- (const Vector& rhs) const
{
	Vector res(*this);
	res -= rhs;
	return res;
}



/* inline Vector& Vector::operator+= (const Vector& rhs)
 *
 * Add vector to current
 */

inline Vector& Vector::operator+= (const Vector& rhs)
{
	
	// Add new values
	int one = 1;
	double scale = 1;
	DAXPY(&_length, &scale, rhs._vector, &one, _vector, &one);
	
	// Return result
	return *this;
}



/* inline Vector& Vector::operator-= (const Vector& rhs)
 *
 * Subtract vector from current
 */

inline Vector& Vector::operator-= (const Vector& rhs)
{
	
	// Add new values
	int one = 1;
	double scale = -1;
	DAXPY(&_length, &scale, rhs._vector, &one, _vector, &one);
	
	// Return result
	return *this;
}



/* inline Vector& Vector::operator*= (const Matrix& rhs)
 *
 * Multiply vector by a matrix
 */

inline Vector& Vector::operator*= (const Matrix& rhs)
{
	
	// Save current values
	Vector orig(*this);
	
	// Set the new length
	length(rhs.numRows());
	
	// Multiply by matrix
	int inc = 1;
	char trans = 'n';
	double beta = 0;
	double alpha = 1;
	DGEMV(&trans, &rhs._numRows, &rhs._numCols, &alpha, rhs._matrix, &rhs._numRows, orig._vector, &inc, &beta, \
		   _vector, &inc);
	
	// Return result
	return *this;
}



/* inline Vector& Vector::operator*= (double rhs)
 *
 * Multiply vector by constant value
 */

inline Vector& Vector::operator*= (double rhs)
{
	int inc = 1;
	DSCAL(&_length, &rhs, _vector, &inc);
	return *this;
}



/* inline Vector& Vector::operator/= (double rhs)
 *
 * Divide vector by constant value
 */

inline Vector& Vector::operator/= (double rhs)
{
	return this->operator*=(1.0/rhs);
}



// =====================================================================================================================
// Complex vector
// =====================================================================================================================

/* inline void CVector::initialize()
 *
 * Initialize vector
 */

inline void CVector::initialize()
{
	_length = 0;
	_vector = 0;
	_actualLength = 0;
}



/* inline CVector::CVector()
 *
 * Constructors for CVector
 */

inline CVector::CVector()
{
	initialize();
}

inline CVector::CVector(const CVector& copy)
{
	initialize();
	*this = copy;
}

inline CVector::CVector(int inLength)
{
	initialize();
	length(inLength);
}



/* inline CVector::~CVector()
 * 
 * Destructor for CVector
 */

inline CVector::~CVector()
{
	clear();
}



/* inline void CVector::clear()
 *
 * Clear data in CVector
 */

inline void CVector::clear()
{
	if (_vector)
		delete [] _vector;
	initialize();
}



/* inline void CVector::length(int inLength)
 *
 * Set length of CVector
 */

inline void CVector::length(int inLength)
{
	
	// Length is already long enough
	if (inLength <= _actualLength)
	{
		_length = inLength;
		return;
	}
	
	// Create new list
	Complex* temp = new Complex [inLength];
	
	// Copy old values
	int inc = 1;
	ZCOPY(&_length, _vector, &inc, temp, &inc);
	
	// Delete old list
	if (_vector)
		delete [] _vector;
	
	// Save list
	_vector = temp;
	_length = inLength;
	_actualLength = inLength;
}



/* inline CVector::CVector& operator= (const CVector& rhs)
 *
 * Assignment operator for CVector
 */

inline CVector& CVector::operator= (const CVector& rhs)
{
	
	// Make sure not copying self
	if (this != &rhs)
	{

		// Allocate space
		length(rhs._length);

		// Copy list
		int inc = 1;
		ZCOPY(&rhs._length, rhs._vector, &inc, _vector, &inc);
	}
	
	// Return result
	return *this;
}



// =====================================================================================================================
// 3D vector
// =====================================================================================================================

/* inline Vector3D::Vector3D()
 *
 * Constructors for Vector3D object
 */

inline Vector3D::Vector3D()
{ }

inline Vector3D::Vector3D(double value)
{
	fill(value);
}

inline Vector3D::Vector3D(const double* vector)
{
	set(vector);
}

inline Vector3D::Vector3D(double val1, double val2, double val3)
{
	set(val1, val2, val3);
}

inline Vector3D::Vector3D(const Vector3D& copy)
{
	*this = copy;
}



/* inline void Vector3D::fill(double input)
 *
 * Fill vector with constant value
 */

inline void Vector3D::fill(double input)
{
	_vector[0] = input;
	_vector[1] = input;
	_vector[2] = input; 
}



/* inline void Vector3D::set(double val1, double val2, double val3)
 *
 * Set values
 */

inline void Vector3D::set(double val1, double val2, double val3)
{
	_vector[0] = val1;
	_vector[1] = val2;
	_vector[2] = val3;
}



/* inline void Vector3D::set(const double* vector)
 *
 * Set values
 */

inline void Vector3D::set(const double* vector)
{
	_vector[0] = vector[0];
	_vector[1] = vector[1];
	_vector[2] = vector[2];
}



/* inline double Vector3D::magnitude() const
 *
 * Get vector magnitude
 */

inline double Vector3D::magnitude() const
{
	return sqrt(_vector[0]*_vector[0] + _vector[1]*_vector[1] + _vector[2]*_vector[2]);
}



/* inline bool Vector3D::operator== (const Vector3D& rhs) const
 *
 * Return whether two vectors are the same
 */

inline bool Vector3D::operator== (const Vector3D& rhs) const
{
	
	// Get the max value in both lists
	int i;
	double cur;
	double scale = Num<double>::max(Num<double>::abs(_vector[0]), Num<double>::abs(rhs._vector[0]));
	for (i = 1; i < 3; ++i)
	{
		cur = Num<double>::max(Num<double>::abs(_vector[i]), Num<double>::abs(rhs._vector[i]));
		scale = (cur > scale) ? cur : scale;
	}
	scale = scale < 1e-14 ? 1 : scale;
	
	// Check values in each list
	for (int i = 0; i < 3; ++i)
	{
		if (Num<double>::abs(_vector[i] - rhs._vector[i]) / scale > 1e-6)
			return false;
	}
	
	// Return that lists are the same if at this point
	return true;
}



/* inline bool Vector3D::operator!= (const Vector3D& rhs) const
 *
 * Return whether two vectors are not the same
 */

inline bool Vector3D::operator!= (const Vector3D& rhs) const
{
	return !(*this == rhs);
}



/* inline Vector3D& Vector3D::operator= (double value)
 *
 * Assign constant value to vector
 */

inline Vector3D& Vector3D::operator= (double value)
{
	_vector[0] = value;
	_vector[1] = value;
	_vector[2] = value;
	return *this;
}



/* inline Vector3D& Vector3D::operator= (const double* vector)
 *
 * Assign array to vector
 */

inline Vector3D& Vector3D::operator= (const double* vector)
{
	_vector[0] = vector[0];
	_vector[1] = vector[1];
	_vector[2] = vector[2];
	return *this;
}



/* inline Vector3D& Vector3D::operator= (const Vector3D& rhs)
 *
 * Assignment operator for Vector3D object
 */

inline Vector3D& Vector3D::operator= (const Vector3D& rhs)
{
	if (this != &rhs)
	{
		_vector[0] = rhs._vector[0];
		_vector[1] = rhs._vector[1];
		_vector[2] = rhs._vector[2];
	}
	return *this;
}



/* inline double Vector3D::operator* (const Vecto3D& rhs) const
 * 
 * Dot product
 */

inline double Vector3D::operator* (const Vector3D& rhs) const
{
	return _vector[0]*rhs._vector[0] + _vector[1]*rhs._vector[1] + _vector[2]*rhs._vector[2];
}



/* inline double Vector3D::angle(const Vector3D& vec2) const
 *
 * Return the angle between vector and another
 */

inline double Vector3D::angle(const Vector3D& vec2) const
{
	double temp = (*this * vec2) / (magnitude() * vec2.magnitude());
	if ((temp > 1) || (temp < -1))
        return 0;
    return acos(temp);
}



/* inline Vector3D Vector3D::operator+ (const Vector3D& rhs) const
 *
 * Add two vectors
 */

inline Vector3D Vector3D::operator+ (const Vector3D& rhs) const
{
	return Vector3D(_vector[0]+rhs._vector[0], _vector[1]+rhs._vector[1], _vector[2]+rhs._vector[2]);
}



/* inline Vector3D Vector3D::operator- (const Vector3D& rhs) const
 *
 * Get difference between two vectors
 */

inline Vector3D Vector3D::operator- (const Vector3D& rhs) const
{
	return Vector3D(_vector[0]-rhs._vector[0], _vector[1]-rhs._vector[1], _vector[2]-rhs._vector[2]);
}



/* inline Vector3D Vector3D::operator* (double rhs) const
 *
 * Return vector multiplied by constant
 */

inline Vector3D Vector3D::operator* (double rhs) const
{
	return Vector3D(_vector[0]*rhs, _vector[1]*rhs, _vector[2]*rhs);
}



/* inline Vector3D Vector3D::operator/ (double rhs) const
 *
 * Return vector divided by constant
 */

inline Vector3D Vector3D::operator/ (double rhs) const
{
	return Vector3D(_vector[0]/rhs, _vector[1]/rhs, _vector[2]/rhs);
}



/* inline Vector3D& Vector3D::operator+= (const Vector3D& rhs)
 *
 * Add vector to current
 */

inline Vector3D& Vector3D::operator+= (const Vector3D& rhs)
{
	_vector[0] += rhs._vector[0];
	_vector[1] += rhs._vector[1];
	_vector[2] += rhs._vector[2];
	return *this;
}

inline Vector3D& Vector3D::operator+= (const double* rhs)
{
	_vector[0] += rhs[0];
	_vector[1] += rhs[1];
	_vector[2] += rhs[2];
	return *this;
}



/* inline Vector3D& Vector3D::operator-= (const Vector3D& rhs)
 *
 * Add vector to current
 */

inline Vector3D& Vector3D::operator-= (const Vector3D& rhs)
{
	_vector[0] -= rhs._vector[0];
	_vector[1] -= rhs._vector[1];
	_vector[2] -= rhs._vector[2];
	return *this;
}



/* inline Vector3D& Vector3D::operator*= (const Matrix3D& rhs)
 *
 * Get product of vector with matrix
 */

inline Vector3D& Vector3D::operator*= (const Matrix3D& rhs)
{
	
	// Save current values
	double orig[3] = {_vector[0], _vector[1], _vector[2]};
	
	// Multiply by matrix
	_vector[0] = rhs._matrix[0]*orig[0] + rhs._matrix[1]*orig[1] + rhs._matrix[2]*orig[2];
	_vector[1] = rhs._matrix[3]*orig[0] + rhs._matrix[4]*orig[1] + rhs._matrix[5]*orig[2];
	_vector[2] = rhs._matrix[6]*orig[0] + rhs._matrix[7]*orig[1] + rhs._matrix[8]*orig[2];
	
	// Return result
	return *this;
}



/* inline Vector3D& Vector3D::operator*= (double rhs)
 *
 * Multiple vector by constant value
 */

inline Vector3D& Vector3D::operator*= (double rhs)
{
	_vector[0] *= rhs;
	_vector[1] *= rhs;
	_vector[2] *= rhs;
	return *this;
}



/* inline Vector3D& Vector3D::operator/= (double rhs)
 *
 * Divide vector by constant value
 */

inline Vector3D& Vector3D::operator/= (double rhs)
{
	_vector[0] /= rhs;
	_vector[1] /= rhs;
	_vector[2] /= rhs;
	return *this;
}



/* inline Vector3D Vector3D::cross(const Vector3D& rhs) const
 *
 * Return cross product between vectors
 */

inline Vector3D Vector3D::cross(const Vector3D& rhs) const
{
	return Vector3D(_vector[1] * rhs[2] - _vector[2] * rhs[1], \
					_vector[2] * rhs[0] - _vector[0] * rhs[2], \
					_vector[0] * rhs[1] - _vector[1] * rhs[0]);
}



/* inline double Vector3D::distanceToPoint(const Vector3D& rhs) const
 *
 * Return the distance between vectors 
 */

inline double Vector3D::distanceToPoint(const Vector3D& rhs) const
{
	double res = 0;
	double temp = _vector[0] - rhs[0];
	res += temp*temp;
	temp = _vector[1] - rhs[1];
	res += temp*temp;
	temp = _vector[2] - rhs[2];
	res += temp*temp;
	return sqrt(res);
}



/* inline double Vector3D::distanceToLineSegment(const Vector3D& begin, const Vector3D& end) const
 *
 * Return the nearest distance of point to line segment determined by begin and end
 */

inline double Vector3D::distanceToLineSegment(const Vector3D& begin, const Vector3D& end) const
{

	// Line segment length is zero so return distance from point to start/end
	Vector3D vec = end;
	vec -= begin;
	double lenSquared = vec*vec;
	if (lenSquared < 1e-10)
		return distanceToPoint(begin);
	
	// Get the t value in parametric equation for the line
	double t = 0;
	for (int i = 0; i < 3; i++)
		t += (_vector[i] - begin[i])*vec[i];
	t /= lenSquared;
	
	// Nearest distance is before start
	if (t < 0)
		return distanceToPoint(begin);
	
	// Nearest distance is after end
	if (t > 1)
		return distanceToPoint(end);
	
	// Get the point on the line segment where the nearest distance occurs
	vec *= t;
	vec += begin;
	
	// Return distance to point
	return distanceToPoint(vec);
}



/* inline double Vector3D::nearestPointOnPlane(const Vector3D& pointOnPlane, const Vector3D& norm,
 *		Vector3D* nearPoint) const
 *
 * Calculate nearest point on set plane
 */

inline double Vector3D::nearestPointOnPlane(const Vector3D& pointOnPlane, const Vector3D& norm, \
	Vector3D* nearPoint) const
{
	
	// Calculate vector from point to point on plane
	Vector3D vec = pointOnPlane;
	vec -= *this;
	
	// Save normalized norm
	Vector3D normNorm = norm;
	normNorm /= normNorm.magnitude();
	
	// Calculate distance from point to plane
	double distance = vec * normNorm;
	
	// Save point if needed
	if (nearPoint)
	{
		*nearPoint  = normNorm;
		*nearPoint *= distance;
		*nearPoint += *this;
	}
	
	// Return near distance
	return distance;
}



/* inline double Vector3D::nearestPointOnPlane(const Vector3D& direction, const Vector3D& pointOnPlane, \
 *		const Vector3D& norm, Vector3D* nearPoint) const
 *
 * Calculate nearest point on set plane along set direction
 */

inline double Vector3D::nearestPointOnPlane(const Vector3D& direction, const Vector3D& pointOnPlane, \
	const Vector3D& norm, Vector3D* nearPoint) const
{
	
	// Calculate vector from point to point on plane
	Vector3D vec = pointOnPlane;
	vec -= *this;
	
	// Calculate distance
	double distance = direction * norm;
	distance = distance == 0 ? 1e10 : (vec * norm) / distance;
	
	// Save point
	if (nearPoint)
	{
		*nearPoint  = direction;
		*nearPoint *= distance;
		*nearPoint += *this;
	}
	
	// Return near distance
	return distance * direction.magnitude();
}



/* inline void Vector3D::voronoi(const OList<Vector3D>& points, double tol, OList<Vector3D>* vertices, \
 *		OList<Vector3D>* centers, OList<Vector3D>* norms, List<int>* indices) const
 *
 * Calculate Voronoi volume
 */

inline void Vector3D::voronoi(const OList<Vector3D>& points, double tol, OList<Vector3D>* vertices, \
	OList<Vector3D>* centers, OList<Vector3D>* norms, List<int>* indices) const
{
	List<double> weights(points.length(), 0.5);
	voronoi(points, weights, tol, vertices, centers, norms, indices);
}



/* inline void Vector3D::voronoi(const OList<Vector3D>& points, const List<double>& weights, double tol,
 *		OList<Vector3D>* vertices, OList<Vector3D>* centers, OList<Vector3D>* norms, List<int>* indices) const
 *
 * Calculate Voronoi volume
 */

inline void Vector3D::voronoi(const OList<Vector3D>& points, const List<double>& weights, double tol, \
	OList<Vector3D>* vertices, OList<Vector3D>* centers, OList<Vector3D>* norms, List<int>* indices) const
{
	
	// Variables to store results
	Linked<int> linkedIndices;
	Linked<double> linkedDistances;
	Linked<Vector3D> linkedCenters;
	Linked<Vector3D> linkedNorms;
	
	// Loop over points
	int i;
	bool saveCurrent;
	double curDis;
	double intercept;
	Vector3D curVec;
	Vector3D curPoint;
	Linked<int>::iterator itI;
	Linked<double>::iterator itD;
	Linked<Vector3D>::iterator itC;
	Linked<Vector3D>::iterator itN;
	for (i = 0; i < points.length(); ++i)
	{
		
		// Save length to plane along current path
		curVec  = points[i];
		curVec -= *this;
		curVec /= curVec.magnitude();
		curDis  = weights[i] * distanceToPoint(points[i]);
		curPoint  = curVec;
		curPoint *= curDis;
		curPoint += *this;
		
		// Loop over saved planes and check if current passes through any before correct distance
		saveCurrent = true;
		for (itC = linkedCenters.begin(), itN = linkedNorms.begin(); itC != linkedCenters.end(); ++itC, ++itN)
		{
			
			// Current plane intersection is within distance
			intercept = nearestPointOnPlane(curVec, *itC, *itN);
			if ((intercept > 0) && (Num<double>::abs(intercept) < curDis+tol))
			{
				saveCurrent = false;
				break;
			}
		}
		
		// Skip if not saving current plane
		if (!saveCurrent)
			continue;
		
		// Loop over saved planes and check if any are now removed
		itI = linkedIndices.begin();
		itD = linkedDistances.begin();
		itC = linkedCenters.begin();
		itN = linkedNorms.begin();
		Linked<Linked<int>::iterator> itIToRemove;
		Linked<Linked<double>::iterator> itDToRemove;
		Linked<Linked<Vector3D>::iterator> itCToRemove;
		Linked<Linked<Vector3D>::iterator> itNToRemove;
		for(; itI != linkedIndices.end(); ++itI, ++itD, ++itC, ++itN)
		{
			
			// Intersection with new plane breaks old
			intercept = nearestPointOnPlane(*itN, curPoint, curVec);
			if ((intercept > 0) && (Num<double>::abs(intercept) < *itD+tol))
			{
				itIToRemove += itI;
				itDToRemove += itD;
				itCToRemove += itC;
				itNToRemove += itN;
			}
		}
		
		// Remove planes that were replaced
		Linked<Linked<int>::iterator>::iterator itITR = itIToRemove.begin();
		Linked<Linked<double>::iterator>::iterator itDTR = itDToRemove.begin();
		Linked<Linked<Vector3D>::iterator>::iterator itCTR = itCToRemove.begin();
		Linked<Linked<Vector3D>::iterator>::iterator itNTR = itNToRemove.begin();
		for (; itITR != itIToRemove.end(); ++itITR, ++itDTR, ++itCTR, ++itNTR)
		{
			linkedIndices.remove(*itITR);
			linkedDistances.remove(*itDTR);
			linkedCenters.remove(*itCTR);
			linkedNorms.remove(*itNTR);
		}
		
		// Save new plane
		linkedIndices += i;
		linkedDistances += curDis;
		linkedNorms += curVec;
		linkedCenters += curPoint;
	}
	
	// Save any values that are needed
	if (centers)
	{
		centers->length(linkedCenters.length());
		for (i = 0, itC = linkedCenters.begin(); itC != linkedCenters.end(); ++i, ++itC)
			(*centers)[i] = *itC;
	}
	if (norms)
	{
		norms->length(linkedNorms.length());
		for (i = 0, itN = linkedNorms.begin(); itN != linkedNorms.end(); ++i, ++itN)
			(*norms)[i] = *itN;
	}
	if (indices)
	{
		indices->length(linkedIndices.length());
		for (i = 0, itI = linkedIndices.begin(); itI != linkedIndices.end(); ++i, ++itI)
			(*indices)[i] = *itI;
	}
	
	// Save vertices if needed
	if (vertices)
	{
		
		// Clear space
		vertices->length(0);
		
		// Loop over all triplets of planes
		double det;
		Matrix3D normMat;
		Linked<Vector3D>::iterator itC2, itC3;
		Linked<Vector3D>::iterator itN2, itN3;
		for (itC = linkedCenters.begin(), itN = linkedNorms.begin(); itC != linkedCenters.end(); ++itC, ++itN)
		{
			for (itC2 = itC + 1, itN2 = itN + 1; itC2 != linkedCenters.end(); ++itC2, ++itN2)
			{
				for (itC3 = itC2 + 1, itN3 = itN2 + 1; itC3 != linkedCenters.end(); ++itC3, ++itN3)
				{
					
					// Save norms as a matrix
					for (i = 0; i < 3; ++i)
					{
						normMat(i, 0) = (*itN )[i];
						normMat(i, 1) = (*itN2)[i];
						normMat(i, 2) = (*itN3)[i];
					}
					
					// There is not intersection if determinant is zero
					det = normMat.determinant();
					if (Num<double>::abs(det) < 1e-8)
						continue;
					
					// Save new vertex
					vertices->add();
					vertices->last()  = (*itN2).cross(*itN3) * (*itC  * *itN ) / det;
					vertices->last() += (*itN3).cross(*itN ) * (*itC2 * *itN2) / det;
					vertices->last() += (*itN ).cross(*itN2) * (*itC3 * *itN3) / det;
				}
			}
		}
	}
}



// =====================================================================================================================
// Matrix
// =====================================================================================================================

/* inline Matrix::Matrix()
 *
 * Matrix constructors
 */

inline Matrix::Matrix()
{
	initialize();
}

inline Matrix::Matrix(const Matrix& copy)
{
	initialize();
	*this = copy;
}

inline Matrix::Matrix(int inSize)
{
	initialize();
	size(inSize, inSize);
}

inline Matrix::Matrix(int inRows, int inCols)
{
	initialize();
	size(inRows, inCols); 
}
	


/* inline Matrix::~Matrix()
 *
 * Destructor for matrix
 */

inline Matrix::~Matrix()
{
	clear();
}



/* inline void Matrix::initialize()
 *
 * Initialize matrix
 */

inline void Matrix::initialize()
{
	_numRows = 0;
	_numCols = 0;
	_length = 0;
	_matrix = 0;
}



/* inline void Matrix::clear()
 *
 * Clear data in matrix
 */

inline void Matrix::clear()
{
	if (_matrix)
		delete [] _matrix;
	initialize();
}



/* inline void Matrix::size(int inSize, bool preserve)
 *
 * Set size of matrix
 */

inline void Matrix::size(int inSize, bool preserve)
{
	size(inSize, inSize, preserve);
}



/* inline void Matrix::size(int inRows, int inCols, bool preserve)
 *
 * Set the size of a matrix
 */

inline void Matrix::size(int inRows, int inCols, bool preserve)
{
	
	// Return if no change
	if ((inRows == _numRows) && (inCols == _numCols))
		return;
	
	// Not worried about keeping original data
	int newLength = inRows*inCols;
	if (!preserve)
	{
		
		// Increase size of matrix if needed
		if (newLength > _length)
		{
			clear();
			_matrix = new double[newLength];
			_length = newLength;
		}
	}
	
	// Need to keep original data
	else
	{
		
		// Make a copy of the current data
		int inc = 1;
		int origSize = _numRows*_numCols;
		double origData[origSize];
		DCOPY(&origSize, _matrix, &inc, origData, &inc);
		
		// Increase the size of the matrix if needed
		if (newLength > _length)
		{
			clear();
			_matrix = new double[newLength];
			_length = newLength;
		}
		
		// Copy old data
		int i, j, k, m;
		for (i = 0, j = 0; (i < inCols) && (j < _numCols); ++i)
		{
			for (j = 0, k = inRows*i, m = _numRows*i; (j < inRows) && (j < _numRows); ++j, ++k, ++m)
				_matrix[k] = origData[m];
		}
	}
	
	// Save new size
	_numRows = inRows;
	_numCols = inCols;
}



/* inline void Matrix::fill(double input)
 *
 * Fill matrix with constant value
 */

inline void Matrix::fill(double input)
{
	int i;
	int size = _numRows*_numCols;
	for (i = 0; i < size-5;)
	{
		_matrix[i++] = input;
		_matrix[i++] = input;
		_matrix[i++] = input;
		_matrix[i++] = input;
		_matrix[i++] = input;
	}
	for ( ; i < size; )
		_matrix[i++] = input;
}



/* inline void Matrix::swapRows(int row1, int row2)
 *
 * Swap rows of matrix
 */

inline void Matrix::swapRows(int row1, int row2)
{
	int index1;
	int index2;
	double temp;
	for (int i = 0; i < _numCols; ++i)
	{
		index1 = _numRows*i + row1;
		index2 = _numRows*i + row2;
		temp = _matrix[index1];
		_matrix[index1] = _matrix[index2];
		_matrix[index2] = temp;
	}
}



/* inline void Matrix::swapColumns(int col1, int col2)
 *
 * Swap columns of matrix
 */

inline void Matrix::swapColumns(int col1, int col2)
{
	int index1 = col1 * _numRows;
	int index2 = col2 * _numRows;
	double temp;
	for (int i = 0; i < _numRows; ++i, ++index1, ++index2)
	{
		temp = _matrix[index1];
		_matrix[index1] = _matrix[index2];
		_matrix[index2] = temp;
	}
}



/* inline void Matrix::makeIdentity()
 *
 * Make identity matrix
 */

inline void Matrix::makeIdentity()
{
	fill(0);
	for (int i = 0; (i < _numRows) && (i < _numCols); ++i)
		_matrix[i*_numRows + i] = 1;
}



/* inline Matrix Matrix::identity(int size)
 *
 * Return identity matrix
 */

inline Matrix Matrix::identity(int size)
{
	Matrix res(size);
	res.makeIdentity();
	return res;
}



/* inline Matrix Matrix::identity(int numRows, int numCols)
 *
 * Return identity matrix
 */

inline Matrix Matrix::identity(int numRows, int numCols)
{
	Matrix res(numRows, numCols);
	res.makeIdentity();
	return res;
}



/* inline bool void Matrix::operator== (const Matrix& rhs) const
 *
 * Return whether two matrices are the same
 */

inline bool Matrix::operator== (const Matrix& rhs) const
{
	
	// Sizes are not the same
	if (_numRows != rhs._numRows)
		return false;
	if (_numCols != rhs._numCols)
		return false;
	
	// Get the max value in both lists
	int i;
	double cur;
	double scale = Num<double>::max(Num<double>::abs(_matrix[0]), Num<double>::abs(rhs._matrix[0]));
	for (i = 1; i < _numRows*_numCols; ++i)
	{
		cur = Num<double>::max(Num<double>::abs(_matrix[i]), Num<double>::abs(rhs._matrix[i]));
		scale = (cur > scale) ? cur : scale;
	}
	scale = scale < 1e-14 ? 1 : scale;
	
	// Check values in each list
	for (int i = 0; i < _numRows*_numCols; ++i)
	{
		if (Num<double>::abs(_matrix[i] - rhs._matrix[i]) / scale > 1e-6)
			return false;
	}
	
	// Return that matrices are the same if at this point
	return true;
}



/* inline Matrix& Matrix::operator= (const Matrix& rhs)
 *
 * Assignment operator for Matrix object
 */

inline Matrix& Matrix::operator= (const Matrix& rhs)
{
	
	// Make sure not copying self
	if (this != &rhs)
	{

		// Allocate space
		size(rhs._numRows, rhs._numCols, false);

		// Copy matrix
		int inc = 1;
		int size = _numRows * _numCols;
		DCOPY(&size, rhs._matrix, &inc, _matrix, &inc);
	}
	
	// Return result
	return *this;
}



/* inline bool Matrix::isDiagonal() const
 *
 * Return whether matrix is in diagonal form
 */

inline bool Matrix::isDiagonal() const
{
	
	// Finished if no rows or no columns
	if ((!_numRows) || (!_numCols))
		return true;
	
	// Look for max value on diagonal
	int i;
	double temp;
	double tol = Num<double>::abs(_matrix[0]);
	for (i = 0; i < Num<int>::min(_numRows, _numCols); ++i)
	{
		temp = Num<double>::abs(_matrix[i*_numCols + i]);
		if (temp > tol)
			tol = temp;
	}
	
	// Set tolerance
	tol = (tol == 0) ? 1e-10 : tol / 1e10;
	
	// Check for non zero entries off diagonal
	int j, k;
	for (i = 0; i < _numCols; ++i)
	{
		for (j = 0, k = i*_numRows; j < _numRows; ++j, ++k)
		{
			if (i == j)
				continue;
			if (Num<double>::abs(_matrix[k]) > tol)
				return false;
		}
	}
	
	// Return that matrix is diagonal if at this point
	return true;
}



/* inline double Matrix::determinant() const
 *
 * Return the determinant of the matrix
 */

inline double Matrix::determinant() const
{
	
	// 2x2 matrix
	if ((_numRows == 2) && (_numCols == 2))
		return _matrix[0] * _matrix[3] - _matrix[1] * _matrix[2];
	
	// 3x3 matrix
	if ((_numRows == 3) && (_numCols == 3))
		return _matrix[0] * (_matrix[4] * _matrix[8] - _matrix[5] * _matrix[7]) - \
			   _matrix[3] * (_matrix[1] * _matrix[8] - _matrix[7] * _matrix[2]) + \
			   _matrix[6] * (_matrix[1] * _matrix[5] - _matrix[4] * _matrix[2]);
	
	// Other size use LU decomposition
	int i;
	int info;
	int pivots[Num<int>::min(_numRows, _numCols)];
	double work[_numRows*_numCols];
	for (i = 0; i < _numRows*_numCols; ++i)
		work[i] = _matrix[i];
	
	// lapack routine for LU decomposition
	DGETRF(&_numRows, &_numCols, work, &_numRows, pivots, &info);
	
	// Determinant is product of diagonal elements
	double det = 1;
	for (i = 0; (i < _numRows) && (i < _numCols); ++i)
		det *= work[i*_numRows + i];
	return det;
}



/* inline Vector Matrix::operator* (const Vector& rhs) const
 *
 * Multiply vector by matrix
 */

inline Vector Matrix::operator* (const Vector& rhs) const
{
	
	// Variable to store result
	Vector res(_numRows);
	
	// Calculate product
	int inc = 1;
	char trans = 'n';
	double beta = 0;
	double alpha = 1;
	DGEMV(&trans, &_numRows, &_numCols, &alpha, _matrix, &_numRows, rhs(), &inc, &beta, res(), &inc);
	
	// Return result
	return res;
}



/* inline Matrix Matrix::operator* (const Matrix& rhs) const
 *
 * Multiply two matrices
 */

inline Matrix Matrix::operator* (const Matrix& rhs) const
{
	
	// Variable to store result
	Matrix res(_numRows, rhs._numCols);
	
	// Calculate product
	char trans = 'n';
	double beta = 0;
	double alpha = 1;
	DGEMM(&trans, &trans, &_numRows, &rhs._numCols, &_numCols, &alpha, _matrix, &_numRows, rhs._matrix, \
		   &rhs._numRows, &beta, res._matrix, &_numRows);
	
	// Return result
	return res;
}



/* inline Matrix Matrix::transpose() const
 *
 * Return the transpose of a matrix
 */

inline Matrix Matrix::transpose() const
{
	int i, j;
	Matrix res(_numCols, _numRows);
	for (i = 0; i < _numCols; ++i)
	{
		for (j = 0; j < _numRows; ++j)
			res._matrix[_numCols*j + i] = _matrix[_numRows*i + j];
	}
	return res;
}



/* inline Matrix& Matrix::makeTranspose()
 *
 * Return the transpose of a matrix
 */

inline Matrix& Matrix::makeTranspose()
{
	*this = transpose();
	return *this;
}



/* inline Matrix Matrix::inverse() const
 *
 * Return the inverse of a matrix
 */

inline Matrix Matrix::inverse() const
{
	
	// Matrix is not 3x3
	if ((_numRows != 3) || (_numCols != 3))
		return Matrix(0);
	
	// Variable to store result
	Matrix res(3);
	
	// Calculate inverse
	double det = determinant();
	res._matrix[0] = (_matrix[4] * _matrix[8] - _matrix[5] * _matrix[7]) / det;
	res._matrix[1] = (_matrix[7] * _matrix[2] - _matrix[8] * _matrix[1]) / det;
	res._matrix[2] = (_matrix[1] * _matrix[5] - _matrix[2] * _matrix[4]) / det;
	res._matrix[3] = (_matrix[6] * _matrix[5] - _matrix[8] * _matrix[3]) / det;
	res._matrix[4] = (_matrix[0] * _matrix[8] - _matrix[2] * _matrix[6]) / det;
	res._matrix[5] = (_matrix[3] * _matrix[2] - _matrix[5] * _matrix[0]) / det;
	res._matrix[6] = (_matrix[3] * _matrix[7] - _matrix[4] * _matrix[6]) / det;
	res._matrix[7] = (_matrix[6] * _matrix[1] - _matrix[7] * _matrix[0]) / det;
	res._matrix[8] = (_matrix[0] * _matrix[4] - _matrix[1] * _matrix[3]) / det;
	
	// Return result
	return res;
}



/* inline Vector Matrix::solve(const Vector& rhs) const
 *
 * Solve for x in Mat.x = rhs
 */

inline Vector Matrix::solve(const Vector& rhs) const
{
	
	// Storage used in solution
	int pivots[_numCols];
	Matrix workMat(*this);
	
	// Solve system
	int one = 1;
	int info;
	Vector res(rhs);
	DGESV(&_numRows, &one, workMat._matrix, &_numCols, pivots, res(), &_numRows, &info);
	
	// Return result
	return res;
}



/* inline Vector Matrix::eigenvalues(Matrix* eigenvectors, bool isSymmetric) const
 *
 * Get eigenvalues for real matrix
 */

inline Vector Matrix::eigenvalues(Matrix* eigenvectors, bool isSymmetric) const
{
	
	// Matrix is not square
	if (_numRows != _numCols)
		return Vector(0);
	
	// Common variables
	int info;
	Matrix matCopy(*this);
	
	// Variable to store eigenvalues
	Vector res(_numRows);
	
	// Matrix is symmetric
	if (isSymmetric)
	{
		
		// Work space
		int lenIWork = eigenvectors == 0 ? 2                  : 2 * (3 + 5*_numRows);
		int lenWork  = eigenvectors == 0 ? 2 * (_numRows + 1) : 2 * (1 + 6*_numRows + 2*_numRows*_numRows);
		int    iWork[lenIWork];
		double  work[lenWork];
		
		// Set job type
		char jobZ = eigenvectors == 0 ? 'N' : 'V';
		char upOrLow = 'U';
		
		// Get eigenvalues
		DSYEVD(&jobZ, &upOrLow, &_numRows, matCopy._matrix, &_numRows, res._vector, work, &lenWork, iWork, \
			   &lenIWork, &info);
		
		// Save eigenvectors
		if (eigenvectors)
			*eigenvectors = matCopy;
	}
	
	// Matrix is not symmetric
	else
	{
		
		// Work space
		int one = 1;
		int workSize = 8*_numRows;
		double work[workSize];
		double imagEVals[_numRows];

		// Set job type
		char jobVL = 'N';
		char jobVR = eigenvectors == 0 ? 'N' : 'V';
		
		// Make sure there is room for eigenvectors if needed
		if (eigenvectors)
			eigenvectors->size(_numRows);
		
		// Dummy variable or pointer to eigenvectors
		double* evecs = eigenvectors == 0 ? 0 : eigenvectors->_matrix;
		
		// Get eigenvalues
		DGEEV(&jobVL, &jobVR, &_numRows, matCopy._matrix, &_numRows, res._vector, imagEVals, 0, &one, \
			  evecs, &_numRows, work, &workSize, &info);
	}
	
	// Return eigenvectors
	return res;
}



/* inline Matrix Matrix::rowEchelon(int* swaps, Matrix* operations, bool integer) const
 *
 * Return row-echelon form for matrix
 */

inline Matrix Matrix::rowEchelon(int* swaps, Matrix* operations, bool integer) const
{
	
	// Finished if no rows or no columns
	if ((!_numRows) || (!_numCols))
		return *this;
	
	// Make copy of matrix
	Matrix res(*this);
	
	// Set the tolerance
	int i;
	double cur;
	double tol = Num<double>::abs(res._matrix[0]);
	for (i = 1; i < res._numRows*res._numCols; ++i)
	{
		cur = Num<double>::abs(res._matrix[i]);
		if (cur > tol)
			tol = cur;
	}
	tol /= 1e8;

	// Loop through all rows
	bool found;
	int rowIndex;
	int colIndex;
	int pivRow = 0;
	int pivCol = 0;
	double tempVal;
	while ((pivRow < res._numRows) && (pivCol < res._numCols))
	{

		// Check if column has any non-zero entries at or below current row
		found = false;
		for (rowIndex = pivRow; rowIndex < res._numRows; ++rowIndex)
		{

			// Current value is nonzero
			if (Num<double>::abs(res(rowIndex, pivCol)) > tol)
			{

				// Swap values on rows
				found = true;
				if (rowIndex != pivRow)
				{
					res.swapRows(pivRow, rowIndex);
					if (swaps)
						Num<int>::swap(swaps[pivRow], swaps[rowIndex]);
					if (operations)
						operations->swapRows(pivRow, rowIndex);
				}
				break;
			}
		}

		// Increase column if needed
		if (!found)
		{
			++pivCol;
			continue;
		}

		// Loop over rest of rows and move smallest absolute value to top
		for (rowIndex = pivRow + 1; rowIndex < res._numRows; ++rowIndex)
		{

			// Found a new minimum
			if ((Num<double>::abs(res(rowIndex, pivCol)) > tol) && \
				(Num<double>::abs(res(rowIndex, pivCol)) < Num<double>::abs(res(pivRow, pivCol))))
			{

				// Swap values on rows
				res.swapRows(pivRow, rowIndex);
				if (swaps)
					Num<int>::swap(swaps[pivRow], swaps[rowIndex]);
				if (operations)
					operations->swapRows(pivRow, rowIndex);
			}
		}

		// Set pivot to a positive value
		if (res(pivRow, pivCol) < 0)
		{
			for (colIndex = 0; colIndex < res._numCols; ++colIndex)
				res(pivRow, colIndex) *= -1;
			if (operations)
			{
				for (colIndex = 0; colIndex < operations->_numCols; ++colIndex)
					(*operations)(pivRow, colIndex) *= -1;
			}
		}

		// Set other rows
		for (rowIndex = pivRow + 1; rowIndex < res._numRows; ++rowIndex)
		{

			// Subtract multiple of pivot row
			if (integer)
				tempVal = Num<double>::floor(res(rowIndex, pivCol) / res(pivRow, pivCol));
			else
				tempVal = res(rowIndex, pivCol) / res(pivRow, pivCol);
			if (Num<double>::abs(tempVal) < tol)
				continue;
			for (colIndex = 0; colIndex < res._numCols; ++colIndex)
				res(rowIndex, colIndex) -= tempVal * res(pivRow, colIndex);
			if (operations)
			{
				for (colIndex = 0; colIndex < operations->_numCols; ++colIndex)
					(*operations)(rowIndex, colIndex) -= tempVal * (*operations)(pivRow, colIndex);
			}
		}

		// Check if all zero below pivot
		found = false;
		for (rowIndex = pivRow + 1; rowIndex < res._numRows; ++rowIndex)
		{
			if (Num<double>::abs(res(rowIndex, pivCol)) > tol)
			{
				found = true;
				break;
			}
		}

		// Go to next row and column
		if (!found)
		{
			++pivRow;
			++pivCol;
		}
	}
	
	// Return matrix
	return res;
}



// =====================================================================================================================
// Complex matrix
// =====================================================================================================================

/* inline void CMatrix::CMatrix::initialize()
 *
 * Initialize CMatrix
 */

inline void CMatrix::CMatrix::initialize()
{
	_numRows = 0;
	_numCols = 0;
	_length = 0;
	_matrix = 0;
}



/* inline CMatrix::CMatrix()
 *
 * Constructors for CMatrix
 */

inline CMatrix::CMatrix()
{
	initialize();
}

inline CMatrix::CMatrix(const CMatrix& copy)
{
	initialize();
	*this = copy;
}

inline CMatrix::CMatrix(int inSize)
{
	initialize();
	size(inSize, inSize);
}

inline CMatrix::CMatrix(int inRows, int inCols)
{
	initialize();
	size(inRows, inCols);
}



/* inline CMatrix::~CMatrix()
 *
 * Destructor for CMatrix
 */

inline CMatrix::~CMatrix()
{
	clear();
}



/* inline void CMatrix::clear()
 *
 * Clear data in CMatrix
 */

inline void CMatrix::clear()
{
	if (_matrix)
		delete [] _matrix;
	initialize();
}



/* inline void CMatrix::size(int inSize, bool preserve)
 *
 * Set size of CMatrix
 */

inline void CMatrix::size(int inSize, bool preserve)
{
	size(inSize, inSize, preserve);
}

inline void CMatrix::size(int inRows, int inCols, bool preserve)
{

	// Return if no change
	if ((inRows == _numRows) && (inCols == _numCols))
		return;
	
	// Not worried about keeping original data
	int newLength = inRows*inCols;
	if (!preserve)
	{
		
		// Increase size of matrix if needed
		if (newLength > _length)
		{
			clear();
			_matrix = new Complex [newLength];
			_length = newLength;
		}
	}
	
	// Need to keep original data
	else
	{
		
		// Make a copy of the current data
		int inc = 1;
		int origSize = _numRows*_numCols;
		Complex origData[origSize];
		ZCOPY(&origSize, _matrix, &inc, origData, &inc);
		
		// Increase the size of the matrix if needed
		if (newLength > _length)
		{
			clear();
			_matrix = new Complex [newLength];
			_length = newLength;
		}
		
		// Copy old data
		int i, j, k, m;
		for (i = 0, j = 0; (i < inCols) && (j < _numCols); ++i)
		{
			for (j = 0, k = inRows*i, m = _numRows*i; (j < inRows) && (j < _numRows); ++j, ++k, ++m)
				_matrix[k] = origData[m];
		}
	}
	
	// Save new size
	_numRows = inRows;
	_numCols = inCols;
}



/* inline void CMatrix::swapRows(int row1, int row2)
 *
 * Swap rows of matrix
 */

inline void CMatrix::swapRows(int row1, int row2)
{
	int index1;
	int index2;
	Complex temp;
	for (int i = 0; i < _numCols; ++i)
	{
		index1 = _numRows*i + row1;
		index2 = _numRows*i + row2;
		temp = _matrix[index1];
		_matrix[index1] = _matrix[index2];
		_matrix[index2] = temp;
	}
}



/* inline void CMatrix::swapColumns(int col1, int col2)
 *
 * Swap columns of matrix
 */

inline void CMatrix::swapColumns(int col1, int col2)
{
	int index1 = col1 * _numRows;
	int index2 = col2 * _numRows;
	Complex temp;
	for (int i = 0; i < _numRows; ++i, ++index1, ++index2)
	{
		temp = _matrix[index1];
		_matrix[index1] = _matrix[index2];
		_matrix[index2] = temp;
	}
}



/* inline CMatrix& CMatrix::operator= (const CMatrix& rhs)
 *
 * Assignment operator for CMatrix
 */

inline CMatrix& CMatrix::operator= (const CMatrix& rhs)
{
	
	// Make sure not copying self
	if (this != &rhs)
	{

		// Allocate space
		size(rhs._numRows, rhs._numCols, false);

		// Copy matrix
		int inc = 1;
		int size = _numRows * _numCols;
		ZCOPY(&size, rhs._matrix, &inc, _matrix, &inc);
	}
	
	// Return result
	return *this;
}



/* inline CVector CMatrix::eigenvalues(CMatrix* eigenvectors, bool isHermitian) const
 *
 * Get eigenvalues of matrix
 */

inline CVector CMatrix::eigenvalues(CMatrix* eigenvectors, bool isHermitian) const
{
	
	// Matrix is not square
	if (_numRows != _numCols)
		return CVector(0);
	
	// Common variables
	int info;
	CMatrix matCopy(*this);
	
	// Variable to store eigenvalues
	CVector res(_numRows);
	
	// Matrix is Hermitian
	if (isHermitian)
	{
		
		// Work space
		int lenIWork = eigenvectors == 0 ? 2                  : 2 * (3 + 5*_numRows);
		int lenRWork = eigenvectors == 0 ? 2 * _numRows       : 2 * (1 + 5*_numRows + 2*_numRows*_numRows);
		int lenWork  = eigenvectors == 0 ? 2 * (_numRows + 1) : 2 * (2*_numRows + _numRows*_numRows);
		int     iWork[lenIWork];
		double  rWork[lenRWork];
		Complex work [lenWork];
		
		// Set job type
		char jobZ = eigenvectors == 0 ? 'N' : 'V';
		char upOrLow = 'U';
		
		// Get eigenvalues
		double evals[_numRows];
		ZHEEVD(&jobZ, &upOrLow, &_numRows, matCopy._matrix, &_numRows, evals, work, &lenWork, rWork, \
			   &lenRWork, iWork, &lenIWork, &info);
		
		// Save eigenvalues
		for (int i = 0; i < _numRows; ++i)
		{
			res._vector[i].real = evals[i];
			res._vector[i].imag = 0;
		}
		
		// Save eigenvectors
		if (eigenvectors)
			*eigenvectors = matCopy;
	}
	
	// Matrix is not hermitian
	else
	{
		
		// Work space
		int one = 1;
		int lWorkSize = 4*_numRows;
		Complex lWork[lWorkSize];
		double  rWork[2*_numRows];

		// Set job type
		char jobVL = 'N';
		char jobVR = eigenvectors == 0 ? 'N' : 'V';
		
		// Make sure there is room for eigenvectors if needed
		if (eigenvectors)
			eigenvectors->size(_numRows);
		
		// Dummy variable or pointer to eigenvectors
		Complex* evecs = eigenvectors == 0 ? 0 : eigenvectors->_matrix;
		
		// Get eigenvalues
		ZGEEV(&jobVL, &jobVR, &_numRows, matCopy._matrix, &_numRows, res._vector, 0, &one, evecs, &_numRows, lWork, \
			  &lWorkSize, rWork, &info);
	}
	
	// Return eigenvectors
	return res;
}



// =====================================================================================================================
// 3D matrix
// =====================================================================================================================

/* inline Matrix3D::Matrix3D()
 * 
 * Matrix3D constructors
 */

inline Matrix3D::Matrix3D()
{  }

inline Matrix3D::Matrix3D(const Matrix3D& copy)
{
	*this = copy;
}

inline Matrix3D::Matrix3D(double value)
{
	fill(value);
}

inline Matrix3D::Matrix3D(double val11, double val12, double val13, double val21, double val22, double val23, \
						  double val31, double val32, double val33)
{
	set(val11, val12, val13, val21, val22, val23, val31, val32, val33);
}



/* inline void Matrix3D::fill(double input)
 *
 * Set all values in matrix
 */

inline void Matrix3D::fill(double input)
{
	_matrix[0] = _matrix[1] = _matrix[2] = _matrix[3] = _matrix[4] = _matrix[5] = _matrix[6] = \
	_matrix[7] = _matrix[8] = input;
}



/* inline void Matrix3D::set(double val11, double val12, double val13, double val21, double val22, double val23,
 *							 double val31, double val32, double val33)
 *
 * Set all values in a matrix
 */

inline void Matrix3D::set(double val11, double val12, double val13, double val21, double val22, double val23, \
						  double val31, double val32, double val33)
{
	_matrix[0] = val11;
	_matrix[1] = val12;
	_matrix[2] = val13;
	_matrix[3] = val21;
	_matrix[4] = val22;
	_matrix[5] = val23;
	_matrix[6] = val31;
	_matrix[7] = val32;
	_matrix[8] = val33;
}



/* inline void Matrix3D::setRow(int row, const double* vector)
 *
 * Set row in matrix
 */

inline void Matrix3D::setRow(int row, const double* vector)
{
	row *= 3;
	for (int i = 0; i < 3; ++i, ++row)
		_matrix[row] = vector[i];
}

inline void Matrix3D::setRow(int row, const Vector3D& vector)
{
	row *= 3;
	for (int i = 0; i < 3; ++i, ++row)
		_matrix[row] = vector[i];
}



/* inline void Matrix3D::swapRows(int row1, int row2)
 *
 * Swap two rows in a matrix
 */

inline void Matrix3D::swapRows(int row1, int row2)
{
	double temp;
	row1 *= 3;
	row2 *= 3;
	temp = _matrix[row1];
	_matrix[row1] = _matrix[row2];
	_matrix[row2] = temp;
	++row1;
	++row2;
	temp = _matrix[row1];
	_matrix[row1] = _matrix[row2];
	_matrix[row2] = temp;
	++row1;
	++row2;
	temp = _matrix[row1];
	_matrix[row1] = _matrix[row2];
	_matrix[row2] = temp;
}



/* inline void Matrix3D::makeIdentity()
 *
 * Make identity matrix
 */

inline void Matrix3D::makeIdentity()
{
	_matrix[0] = _matrix[4] = _matrix[8] = 1;
	_matrix[1] = _matrix[2] = _matrix[3] = _matrix[5] = _matrix[6] = _matrix[7] = 0;
}



/* inline Matrix3D Matrix3D::identity()
 *
 * Return identity matrix
 */

inline Matrix3D Matrix3D::identity()
{
	Matrix3D res;
	res.makeIdentity();
	return res;
}



/* inline bool Matrix3D::operator== (const Matrix3D& rhs) const
 *
 * Return whether two matrices are the same
 */

inline bool Matrix3D::operator== (const Matrix3D& rhs) const
{
	
	// Get the max value in both lists
	int i;
	double cur;
	double scale = Num<double>::max(Num<double>::abs(_matrix[0]), Num<double>::abs(rhs._matrix[0]));
	for (i = 1; i < 9; ++i)
	{
		cur = Num<double>::max(Num<double>::abs(_matrix[i]), Num<double>::abs(rhs._matrix[i]));
		scale = (cur > scale) ? cur : scale;
	}
	scale = scale < 1e-14 ? 1 : scale;
	
	// Check values in each list
	for (int i = 0; i < 9; ++i)
	{
		if (Num<double>::abs(_matrix[i] - rhs._matrix[i]) / scale > 1e-6)
			return false;
	}
	
	// Return that matrices are the same if at this point
	return true;
}



/* inline bool Matrix3D::operator!= (const Matrix3D& rhs) const
 *
 * Return whether two matrices are not the same
 */

inline bool Matrix3D::operator!= (const Matrix3D& rhs) const
{
	return !(*this == rhs);
}



/* inline bool Matrix3D::isInteger(double tol) const
 *
 * Return whether matrix contains all integer values up to tol
 */
inline bool Matrix3D::isInteger(double tol) const
{
        for (int i = 0; i < 9; ++i)
        {
                if (Num<double>::abs(_matrix[i] - Num<double>::round(_matrix[i], 1)) > tol)
                        return false;
        }
        return true;
}



/* inline Matrix3D& Matrix3D::operator= (const Matrix3D& rhs)
 *
 * Assignment operator for Matrix3D object
 */

inline Matrix3D& Matrix3D::operator= (const Matrix3D& rhs)
{
	
	// Make sure not copying self
	if (this != &rhs)
	{
		_matrix[0] = rhs._matrix[0];
		_matrix[1] = rhs._matrix[1];
		_matrix[2] = rhs._matrix[2];
		_matrix[3] = rhs._matrix[3];
		_matrix[4] = rhs._matrix[4];
		_matrix[5] = rhs._matrix[5];
		_matrix[6] = rhs._matrix[6];
		_matrix[7] = rhs._matrix[7];
		_matrix[8] = rhs._matrix[8];
	}
	
	// Return result
	return *this;
}

inline Matrix3D& Matrix3D::operator= (double rhs)
{
	fill(rhs);
	return *this;
}



/* inline int Matrix3D::rank() const
 *
 * 
 */

inline int Matrix3D::rank() const
{
	
	// Save transpose of matrix
	int size = 3;
	double trans[9] = {_matrix[0], _matrix[3], _matrix[6], _matrix[1], _matrix[4], _matrix[7], _matrix[2], \
					   _matrix[5], _matrix[8]};
	
	// Set tolerance
	double cur;
	double tol = Num<double>::abs(_matrix[0]);
	for (int i = 1; i < 9; ++i)
	{
		cur = Num<double>::abs(_matrix[i]);
		if (cur > tol)
			tol = cur;
	}
	tol /= 1e6;
	
	// Perform singular value decomposition
	int info;
	int one = 1;
	int workSize = 25;
	char job = 'N';
	double work[workSize];
	double singularValues[3];
	DGESVD(&job, &job, &size, &size, trans, &size, singularValues, 0, &one, 0, &one, work, &workSize, &info);
	
	// Count number of non-zero values
	int res = 0;
	if (Num<double>::abs(singularValues[0]) > tol)
		++res;
	if (Num<double>::abs(singularValues[1]) > tol)
		++res;
	if (Num<double>::abs(singularValues[2]) > tol)
		++res;
	
	// Return rank
	return res;
}



/* inline double Matrix3D::trace() const
 *
 * Return trace of the matrix
 */

inline double Matrix3D::trace() const
{
	return _matrix[0] + _matrix[4] + _matrix[8];
}



/* inline double Matrix3D::determinant() const
 *
 * Return the determinant of a matrix
 */

inline double Matrix3D::determinant() const
{
	return _matrix[0] * (_matrix[4] * _matrix[8] - _matrix[7] * _matrix[5]) - \
		   _matrix[1] * (_matrix[3] * _matrix[8] - _matrix[5] * _matrix[6]) + \
		   _matrix[2] * (_matrix[3] * _matrix[7] - _matrix[4] * _matrix[6]);
}



/* inline double Matrix3D::volume() const
 *
 * Return volume of matrix
 */

inline double Matrix3D::volume() const
{
	return determinant();
}



/* inline Matrix3D& Matrix3D::operator+= (const Matrix3D& rhs)
 *
 * Add matrix to current
 */

inline Matrix3D& Matrix3D::operator+= (const Matrix3D& rhs)
{
	_matrix[0] += rhs._matrix[0];
	_matrix[1] += rhs._matrix[1];
	_matrix[2] += rhs._matrix[2];
	_matrix[3] += rhs._matrix[3];
	_matrix[4] += rhs._matrix[4];
	_matrix[5] += rhs._matrix[5];
	_matrix[6] += rhs._matrix[6];
	_matrix[7] += rhs._matrix[7];
	_matrix[8] += rhs._matrix[8];
	return *this;
}



/* inline Matrix3D& Matrix3D::operator-= (const Matrix3D& rhs)
 *
 * Subtract matrix from current
 */

inline Matrix3D& Matrix3D::operator-= (const Matrix3D& rhs)
{
	_matrix[0] -= rhs._matrix[0];
	_matrix[1] -= rhs._matrix[1];
	_matrix[2] -= rhs._matrix[2];
	_matrix[3] -= rhs._matrix[3];
	_matrix[4] -= rhs._matrix[4];
	_matrix[5] -= rhs._matrix[5];
	_matrix[6] -= rhs._matrix[6];
	_matrix[7] -= rhs._matrix[7];
	_matrix[8] -= rhs._matrix[8];
	return *this;
}



/* inline Matrix3D Matrix3D::operator* (double rhs) const
 * 
 * Multiply matrix by a constant value
 */
 
inline Matrix3D Matrix3D::operator* (double rhs) const
{
	return Matrix3D(_matrix[0]*rhs, _matrix[1]*rhs, _matrix[2]*rhs, \
					_matrix[3]*rhs, _matrix[4]*rhs, _matrix[5]*rhs, \
					_matrix[6]*rhs, _matrix[7]*rhs, _matrix[8]*rhs);
}



/* inline Vector3D Matrix3D::operator* (const double* rhs) const
 *
 * Return product of matrix with a vector
 */

inline Vector3D Matrix3D::operator* (const double* rhs) const
{
	return Vector3D(_matrix[0]*rhs[0] + _matrix[1]*rhs[1] + _matrix[2]*rhs[2], \
					_matrix[3]*rhs[0] + _matrix[4]*rhs[1] + _matrix[5]*rhs[2], \
					_matrix[6]*rhs[0] + _matrix[7]*rhs[1] + _matrix[8]*rhs[2]);
}



/* inline Vector3D Matrix3D::operator* (const Vector3D& rhs) const
 *
 * Return the product of matrix with a vector
 */

inline Vector3D Matrix3D::operator* (const Vector3D& rhs) const
{
	return Vector3D(_matrix[0]*rhs[0] + _matrix[1]*rhs[1] + _matrix[2]*rhs[2], \
					_matrix[3]*rhs[0] + _matrix[4]*rhs[1] + _matrix[5]*rhs[2], \
					_matrix[6]*rhs[0] + _matrix[7]*rhs[1] + _matrix[8]*rhs[2]);
}



/* inline Matrix3D Matrix3D::operator* (const Matrix3D& rhs) const
 *
 * Multiply two matrices
 */

inline Matrix3D Matrix3D::operator* (const Matrix3D& rhs) const
{
	double temp[9] = {rhs._matrix[0], rhs._matrix[3], rhs._matrix[6], rhs._matrix[1], rhs._matrix[4], \
					  rhs._matrix[7], rhs._matrix[2], rhs._matrix[5], rhs._matrix[8]};
	return Matrix3D \
		(_matrix[0]*temp[0] + _matrix[1]*temp[1] + _matrix[2]*temp[2], \
		 _matrix[0]*temp[3] + _matrix[1]*temp[4] + _matrix[2]*temp[5], \
		 _matrix[0]*temp[6] + _matrix[1]*temp[7] + _matrix[2]*temp[8], \
		 _matrix[3]*temp[0] + _matrix[4]*temp[1] + _matrix[5]*temp[2], \
		 _matrix[3]*temp[3] + _matrix[4]*temp[4] + _matrix[5]*temp[5], \
		 _matrix[3]*temp[6] + _matrix[4]*temp[7] + _matrix[5]*temp[8], \
		 _matrix[6]*temp[0] + _matrix[7]*temp[1] + _matrix[8]*temp[2], \
		 _matrix[6]*temp[3] + _matrix[7]*temp[4] + _matrix[8]*temp[5], \
		 _matrix[6]*temp[6] + _matrix[7]*temp[7] + _matrix[8]*temp[8]);
}



/* inline Matrix3D& Matrix3D::operator*= (const Matrix3D& rhs)
 *
 * Multiply current matrix by another
 */

inline Matrix3D& Matrix3D::operator*= (const Matrix3D& rhs)
{
	double temp[9] = {_matrix[0], _matrix[3], _matrix[6], _matrix[1], _matrix[4], _matrix[7], _matrix[2], \
					  _matrix[5], _matrix[8]};
	_matrix[0] = rhs._matrix[0]*temp[0] + rhs._matrix[1]*temp[1] + rhs._matrix[2]*temp[2];
	_matrix[1] = rhs._matrix[0]*temp[3] + rhs._matrix[1]*temp[4] + rhs._matrix[2]*temp[5];
	_matrix[2] = rhs._matrix[0]*temp[6] + rhs._matrix[1]*temp[7] + rhs._matrix[2]*temp[8];
	_matrix[3] = rhs._matrix[3]*temp[0] + rhs._matrix[4]*temp[1] + rhs._matrix[5]*temp[2];
	_matrix[4] = rhs._matrix[3]*temp[3] + rhs._matrix[4]*temp[4] + rhs._matrix[5]*temp[5];
	_matrix[5] = rhs._matrix[3]*temp[6] + rhs._matrix[4]*temp[7] + rhs._matrix[5]*temp[8];
	_matrix[6] = rhs._matrix[6]*temp[0] + rhs._matrix[7]*temp[1] + rhs._matrix[8]*temp[2];
	_matrix[7] = rhs._matrix[6]*temp[3] + rhs._matrix[7]*temp[4] + rhs._matrix[8]*temp[5];
	_matrix[8] = rhs._matrix[6]*temp[6] + rhs._matrix[7]*temp[7] + rhs._matrix[8]*temp[8];
	return *this;
}



/* inline Matrix3D& Matrix3D::operator*= (double rhs)
 *
 * Multiply matrix by a value
 */

inline Matrix3D& Matrix3D::operator*= (double rhs)
{
	_matrix[0] *= rhs;
	_matrix[1] *= rhs;
	_matrix[2] *= rhs;
	_matrix[3] *= rhs;
	_matrix[4] *= rhs;
	_matrix[5] *= rhs;
	_matrix[6] *= rhs;
	_matrix[7] *= rhs;
	_matrix[8] *= rhs;
	return *this;
}



/* inline Matrix3D& Matrix3D::operator/= (double rhs)
 *
 * Divide matrix by a value
 */

inline Matrix3D& Matrix3D::operator/= (double rhs)
{
	_matrix[0] /= rhs;
	_matrix[1] /= rhs;
	_matrix[2] /= rhs;
	_matrix[3] /= rhs;
	_matrix[4] /= rhs;
	_matrix[5] /= rhs;
	_matrix[6] /= rhs;
	_matrix[7] /= rhs;
	_matrix[8] /= rhs;
	return *this;
}



/* inline Matrix3D Matrix3D::transpose() const
 *
 * Return the transpose of a matrix
 */

inline Matrix3D Matrix3D::transpose() const
{
	return Matrix3D(_matrix[0], _matrix[3], _matrix[6], _matrix[1], _matrix[4], _matrix[7], _matrix[2], \
					_matrix[5], _matrix[8]);
}



/* inline Matrix3D& Matrix3D::makeTranspose()
 *
 * Return the transpose of a matrix
 */

inline Matrix3D& Matrix3D::makeTranspose()
{
	double temp = _matrix[1];
	_matrix[1] = _matrix[3];
	_matrix[3] = temp;
	temp = _matrix[2];
	_matrix[2] = _matrix[6];
	_matrix[6] = temp;
	temp = _matrix[5];
	_matrix[5] = _matrix[7];
	_matrix[7] = temp;
	return *this;
}



/* inline Matrix3D Matrix3D::inverse() const
 *
 * Calculate the inverse of the matrix
 */

inline Matrix3D Matrix3D::inverse() const
{
	double det = determinant();
	return Matrix3D((_matrix[4] * _matrix[8] - _matrix[7] * _matrix[5]) / det, \
					(_matrix[2] * _matrix[7] - _matrix[8] * _matrix[1]) / det, \
					(_matrix[1] * _matrix[5] - _matrix[4] * _matrix[2]) / det, \
    				(_matrix[5] * _matrix[6] - _matrix[8] * _matrix[3]) / det, \
					(_matrix[0] * _matrix[8] - _matrix[6] * _matrix[2]) / det, \
					(_matrix[2] * _matrix[3] - _matrix[5] * _matrix[0]) / det, \
					(_matrix[3] * _matrix[7] - _matrix[6] * _matrix[4]) / det, \
					(_matrix[1] * _matrix[6] - _matrix[7] * _matrix[0]) / det, \
					(_matrix[0] * _matrix[4] - _matrix[3] * _matrix[1]) / det);
}



/* inline Vector3D Matrix3D::eigenvalues(Matrix3D* eigenvectors) const
 *
 * Calculate eigenvalues for matrix
 */

inline Vector3D Matrix3D::eigenvalues(Matrix3D* eigenvectors) const
{
	
	// Get max value in matrix
	int i;
	double cur;
	double tol = Num<double>::abs(_matrix[0]);
	for (i = 1; i < 9; ++i)
	{
		cur = Num<double>::abs(_matrix[i]);
		if (cur > tol)
			tol = cur;
	}

	// Setup
	tol = (tol == 0) ? 1 : tol / 1e10;
	double thisTrace = this->trace();
	Matrix3D squared = *this * *this;
	double squaredTrace = squared.trace();

	// Get coefficients of characteristic polynomial
	double a2 = -thisTrace;
	double a1 = (thisTrace*thisTrace - squaredTrace)/2.0;
	double a0 = -determinant();

	// Define values
	double Q = (3*a1 - a2*a2) / 9;
	double R = (9*a1*a2 - 27*a0 - 2*pow(a2, 3)) / 54;
	double theta;
	double temp = (Q >= 0) ? 0 : pow(-Q, 1.5);
	if ((Num<double>::abs(R) < tol) && (Num<double>::abs(temp) < tol))
		theta = 0;
	else
	{
		temp = R / temp;
		if (temp >= 1)
			theta = 0;
		else if (temp <= -1)
			theta = Constants::pi;
		else
			theta = acos(temp);
	}

	// Save eigenvalues
	double sqrtQ = (Q >= 0) ? 0 : sqrt(-Q);
	Vector3D eigenvalues(2 * sqrtQ * cos(theta / 3) - a2 / 3, \
						 2 * sqrtQ * cos((theta + 2*Constants::pi) / 3) - a2 / 3, \
						 2 * sqrtQ * cos((theta + 4*Constants::pi) / 3) - a2 / 3);

	// Get eigenvectors if needed
	if (eigenvectors)
	{

		// Loop over eigenvalues
		int j, k;
		double mag;
		Matrix3D eigMat;
		for (i = 0; i < 3; ++i)
		{
			
			// Create matrix
			eigMat = *this;
			eigMat._matrix[0] -= eigenvalues[i];
			eigMat._matrix[4] -= eigenvalues[i];
			eigMat._matrix[8] -= eigenvalues[i];
			
			// Convert matrix to row echelon form
			eigMat = eigMat.rowEchelon();

			// Backward subsitution
			tol = Num<double>::abs(eigenvalues[i]) / 1e8;
			for (j = 2; j >= 0; --j)
			{
				(*eigenvectors)(j, i) = 0;
				for (k = j + 1; k < 3; ++k)
					(*eigenvectors)(j, i) -= eigMat(j, k) * (*eigenvectors)(k, i);
				if (Num<double>::abs(eigMat(j, j)) <= tol)
					(*eigenvectors)(j, i) = 1;
				else
					(*eigenvectors)(j, i) /= eigMat(j, j);
			}

			// Save normalized result
			mag = sqrt((*eigenvectors)(0, i)*(*eigenvectors)(0, i) + (*eigenvectors)(1, i)*(*eigenvectors)(1, i) + \
					   (*eigenvectors)(2, i)*(*eigenvectors)(2, i));
			if (mag != 0)
			{
				(*eigenvectors)(0, i) /= mag;
				(*eigenvectors)(1, i) /= mag;
				(*eigenvectors)(2, i) /= mag;
			}
		}
	}

	// Return eigenvalues
	return eigenvalues;
}



/* inline Matrix3D Matrix3D::rowEchelon(int* swaps, bool integer) const
 *
 * Return row-echelon form for matrix
 */

inline Matrix3D Matrix3D::rowEchelon(int* swaps, Matrix3D* operations, bool integer) const
{
	
	// Make copy of matrix
	Matrix3D res(*this);
	
	// Set the tolerance
	int i;
	double cur;
	double tol = Num<double>::abs(res._matrix[0]);
	for (i = 1; i < 9; ++i)
	{
		cur = Num<double>::abs(res._matrix[i]);
		if (cur > tol)
			tol = cur;
	}
	tol /= 1e8;

	// Loop through all rows
	bool found;
	int rowIndex;
	int colIndex;
	int pivRow = 0;
	int pivCol = 0;
	double tempVal;
	while ((pivRow < 3) && (pivCol < 3))
	{

		// Check if column has any non-zero entries at or below current row
		found = false;
		for (rowIndex = pivRow; rowIndex < 3; ++rowIndex)
		{

			// Current value is nonzero
			if (Num<double>::abs(res(rowIndex, pivCol)) > tol)
			{

				// Swap values on rows
				found = true;
				if (rowIndex != pivRow)
				{
					res.swapRows(pivRow, rowIndex);
					if (swaps)
						Num<int>::swap(swaps[pivRow], swaps[rowIndex]);
					if (operations)
						operations->swapRows(pivRow, rowIndex);
				}
				break;
			}
		}

		// Increase column if needed
		if (!found)
		{
			++pivCol;
			continue;
		}

		// Loop over rest of rows and move smallest absolute value to top
		for (rowIndex = pivRow + 1; rowIndex < 3; ++rowIndex)
		{

			// Found a new minimum
			if ((Num<double>::abs(res(rowIndex, pivCol)) > tol) && \
				(Num<double>::abs(res(rowIndex, pivCol)) < Num<double>::abs(res(pivRow, pivCol))))
			{

				// Swap values on rows
				res.swapRows(pivRow, rowIndex);
				if (swaps)
					Num<int>::swap(swaps[pivRow], swaps[rowIndex]);
				if (operations)
					operations->swapRows(pivRow, rowIndex);
			}
		}

		// Set pivot to a positive value
		if (res(pivRow, pivCol) < 0)
		{
			for (colIndex = 0; colIndex < 3; ++colIndex)
				res(pivRow, colIndex) *= -1;
			if (operations)
			{
				for (colIndex = 0; colIndex < 3; ++colIndex)
					(*operations)(pivRow, colIndex) *= -1;
			}
		}

		// Set other rows
		for (rowIndex = pivRow + 1; rowIndex < 3; ++rowIndex)
		{

			// Subtract multiple of pivot row
			if (integer)
				tempVal = Num<double>::floor(res(rowIndex, pivCol) / res(pivRow, pivCol));
			else
				tempVal = res(rowIndex, pivCol) / res(pivRow, pivCol);
			if (Num<double>::abs(tempVal) < tol)
				continue;
			for (colIndex = 0; colIndex < 3; ++colIndex)
				res(rowIndex, colIndex) -= tempVal * res(pivRow, colIndex);
			if (operations)
			{
				for (colIndex = 0; colIndex < 3; ++colIndex)
					(*operations)(rowIndex, colIndex) -= tempVal * (*operations)(pivRow, colIndex);
			}
		}

		// Check if all zero below pivot
		found = false;
		for (rowIndex = pivRow + 1; rowIndex < 3; ++rowIndex)
		{
			if (Num<double>::abs(res(rowIndex, pivCol)) > tol)
			{
				found = true;
				break;
			}
		}

		// Go to next row and column
		if (!found)
		{
			++pivRow;
			++pivCol;
		}
	}
	
	// Return matrix
	return res;
}




// =====================================================================================================================
// Functors
// =====================================================================================================================

/* inline void Functor<Tclass>::initialize()
 *
 * Initialize functor
 */

template <class Tclass>
inline void Functor<Tclass>::initialize()
{
    _params = 0;
	_objPtr = 0;
	_funPtr = 0;
	_objFunPtr = 0;
	_funPtrVec = 0;
	_objFunPtrVec = 0; 
	_funPtrVecArg = 0;
	_objFunPtrVecArg = 0;
}



/* inline Functor<Tclass>::Functor()
 *
 * Functor constructors
 */

template <class Tclass>
inline Functor<Tclass>::Functor()
{
	initialize();
}

template <class Tclass>
inline Functor<Tclass>::Functor(double (*funPtr)(double))
{
	initialize();
	_funPtr = funPtr;
}

template <class Tclass>
inline Functor<Tclass>::Functor(Tclass* objPtr, double (Tclass::*objFunPtr)(double))
{
	initialize();
	_objPtr = objPtr;
	_objFunPtr = objFunPtr;
}

template <class Tclass>
inline Functor<Tclass>::Functor(double (*funPtr)(const Vector&))
{
	initialize();
	_funPtrVec = funPtr;
}

template <class Tclass>
inline Functor<Tclass>::Functor(Tclass* objPtr, double (Tclass::*funPtr)(const Vector&))
{
	initialize();
	_objPtr = objPtr;
	_objFunPtrVec = funPtr;
}

template <class Tclass>
inline Functor<Tclass>::Functor(double (*funPtr)(const Vector&, double), const Vector* params)
{
	initialize();
        if (params != 0) _params = params;
	_funPtrVecArg = funPtr;
}

template <class Tclass>
inline Functor<Tclass>::Functor(Tclass* objPtr, double (Tclass::*funPtr)(const Vector&, double), const Vector* params)
{
	initialize();
        if (params != 0) _params = params;
	_objPtr = objPtr;
	_objFunPtrVecArg = funPtr;
}



/* inline double Functor<Tclass>>::operator() (double arg)
 *
 * Return value of functor
 */

template <class Tclass>
inline double Functor<Tclass>::operator() (double arg)
{
    if (_params != 0)
	return (_funPtr == 0) ? (*_objPtr.*_objFunPtrVecArg)(*_params, arg) : (*_funPtrVecArg)(*_params, arg);
    else
        return (_funPtr == 0) ? (*_objPtr.*_objFunPtr)(arg) : (*_funPtr)(arg);
}

template <class Tclass>
inline double Functor<Tclass>::operator()(const Vector& params)
{
	return (_funPtrVec == 0) ? (*_objPtr.*_objFunPtrVec)(params) : (*_funPtrVec)(params);
}

template <class Tclass>
inline double Functor<Tclass>::operator()(const Vector& params, double arg)
{
	return (_funPtrVecArg == 0) ? (*_objPtr.*_objFunPtrVecArg)(params, arg) : (*_funPtrVecArg)(params, arg);
}



// =====================================================================================================================
// Vector functors
// =====================================================================================================================

/* inline void VectorFunctor<Tclass>::initialize()
 *
 * Initialize vector functor
 */

template <class Tclass>
inline void VectorFunctor<Tclass>::initialize()
{
	_objPtr = 0;
	_funPtr = 0;
	_objFunPtr = 0;
}



/* inline VectorFunctor<Tclass>::VectorFunctor()
 *
 * Constructors for vector functor
 */

template <class Tclass>
inline VectorFunctor<Tclass>::VectorFunctor()
{
	initialize();
}

template <class Tclass>
inline VectorFunctor<Tclass>::VectorFunctor(Vector (*funPtr)(const Vector&, double))
{
	initialize();
	_funPtr = funPtr;
}

template <class Tclass>
inline VectorFunctor<Tclass>::VectorFunctor(Tclass* objPtr, Vector (Tclass::*objFunPtr)(const Vector&, double))
{
	initialize();
	_objPtr = objPtr;
	_objFunPtr = objFunPtr;
}



/* inline Vector VectorFunctor<Tclass>::operator() (const Vector& params, double arg)
 *
 * Evaluate vector functor
 */

template <class Tclass>
inline Vector VectorFunctor<Tclass>::operator() (const Vector& params, double arg)
{
	return (_funPtr == 0) ? (*_objPtr.*_objFunPtr)(params, arg) : (*_funPtr)(params, arg);
}



// =====================================================================================================================
// General math functions
// =====================================================================================================================

/* inline T Num<T>::ceil(T value)
 *
 * Return the ceil of a value
 */

template <class T>
inline T Num<T>::ceil(T value)
{
	return std::ceil((double)value);
}



/* inline T Num<T>::floor(T value)
 *
 * Return the floor of a value
 */

template <class T>
inline T Num<T>::floor(T value)
{
	return std::floor((double)value);
}



/* inline T Num<T>::abs(T value)
 *
 * Return the absolute value of a number
 */

template <class T>
inline T Num<T>::abs(T value)
{
	return (value < 0) ? value*-1 : value;
}



/* inline T Num<T>::min(T val1, T val2)
 *
 * Return smaller of two values
 */

template <class T>
inline T Num<T>::min(T val1, T val2)
{
	return (val1 < val2) ? val1 : val2;
}



/* inline T Num<T>::max(T val1, T val2)
 *
 * Return larger of two values 
 */

template <class T>
inline T Num<T>::max(T val1, T val2)
{
	return (val1 > val2) ? val1 : val2;
}



/* inline T Num<T>::sign(T value)
 *
 * Return the sign of a value
 */

template <class T>
inline T Num<T>::sign(T value)
{
	return value > 0 ? 1 : (value < 0 ? -1 : 0);
}



/* inline T Num<T>::sign(T value, double tol)
 *
 * Return the sign of a value
 */

template <class T>
inline T Num<T>::sign(T value, double tol)
{
	return lt(value, 0, tol) ? -1 : (gt(value, 0, tol) ? 1 : 0);
}



/* inline T Num<T>::round(T value, double roundTo)
 *
 * Round a number
 */

template <class T>
inline T Num<T>::round(T value, double roundTo)
{
	T temp = value / roundTo;
	if (temp >= 0)
		return roundTo*((temp - Num<double>::floor(temp) >= 0.5) ? Num<double>::ceil(temp) : Num<double>::floor(temp));
	return roundTo*((temp - Num<double>::ceil(temp) <= -0.5) ? Num<double>::floor(temp) : Num<double>::ceil(temp));
}



/* inline void Num<T>::swap(T& val1, T& val2)
 *
 * Swap two numbers 
 */

template <class T>
inline void Num<T>::swap(T& val1, T& val2)
{
	T temp = val1;
	val1 = val2;
	val2 = temp;
}



/* inline T Num<T>::next(T value, T max)
 *
 * Return next value in series
 */

template <class T>
inline T Num<T>::next(T value, T max)
{
	return (value < max) ? value + 1 : 0;
}



/* inline T Num<T>::prev(T value, T max)
 *
 * Return previous value in series
 */

template <class T>
inline T Num<T>::prev(T value, T max)
{
	return (value > 0) ? value - 1: max;
}



/* inline T Num<T>::delta(T x, T y)
 *
 * Delta function
 */

template <class T>
inline T Num<T>::delta(T x, T y)
{
	return (x == y) ? 1 : 0;
}



/* inline T Num<T>::mod(T value, T max)
 *
 * Return decimal value of value % 1
 */

template <class T>
inline T Num<T>::mod(T value, T max)
{
	T temp = value - Num<T>::floor(value);
	return (temp > max) ? temp - 1 : temp;
}



/* inline T Num<T>::gcf(int length, int* values)
 *
 * Get the greate common factor of a list of integers
 */

template <class T>
inline T Num<T>::gcf(int length, int* values)
{
	
	// Get absolute values of numbers
	int i;
	int absVals[length];
	for (i = 0; i < length; ++i)
		absVals[i] = Num<int>::abs(values[i]);
	
	// Initialize smallest value to first that is not zero
	int minIndex;
	for (minIndex = 0; minIndex < length; ++minIndex)
	{
		if (absVals[minIndex] > 0)
			break;
	}
	
	// No nonzero values were found
	if (minIndex == length)
		return 1;
	
	// Find the smallest value
	for (i = minIndex + 1; i < length; ++i)
	{
		if ((absVals[i] > 0) && (absVals[i] < absVals[minIndex]))
			minIndex = i;
	}
	
	// Get the factors of the minimum index
	List<T> factors = Num<T>::factors(absVals[minIndex]);
	
	// Loop over factors and find the largest that factors all numbers
	int j;
	bool found;
	double temp;
	T res = 1;
	for (i = 0; i < factors.length(); ++i)
	{
		if (factors[i] <= res)
			continue;
		found = true;
		for (j = 0; j < length; ++j)
		{
			temp = absVals[j] / (double) factors[i];
			if (Num<double>::abs(temp - Num<double>::floor(temp)) > 1e-10)
			{
				found = false;
				break;
			}
		}
		if (found)
			res = factors[i];
	}
	
	// Return result
	return res;
}



/* inline List<T> Num<T>::factors(int value)
 *
 * Return the factors of a number
 */

template <class T>
inline List<T> Num<T>::factors(int value)
{
	
	// Variable to store result
	List<T> res;
	
	// Value is 1
	if (value == 1)
	{
		res += 1;
		return res;
	}
	
	// Loop over possible values
	for (int i = 1; i <= value / 2.0; ++i)
	{
		
		// Found a factor
		if (value % i == 0)
		{
			res += i;
			if (value / i > value / 2.0)
				res += value / i;
		}
	}
	
	// Return result
	return res;
}



/* inline OList<List<T> > Num<T>::tripleFactors(int value)
 *
 * Return the triple factors of a number
 */

template <class T>
inline OList<List<T> > Num<T>::tripleFactors(int value)
{
	
	// Variable to store result
	OList<List<T> > res;
	
	// Value is 1
	if (value == 1)
	{
		res.add();
		res[0].length(3);
		res[0][0] = res[0][1] = res[0][2] = 1;
		return res;
	}
	
	// Value is 2
	if (value == 2)
	{
		res.add();
		res[0].length(3);
		res[0][0] = res[0][1] = 1;
		res[0][2] = 2;
		return res;
	}
	
	// Find factors
	int i, j;
	int temp;
	for (i = 1; i <= value / 3.0; ++i)
	{
		
		// Current number is a factor
		if (value % i == 0)
		{
			
			// Find second factor
			temp = value / i;
			for (j = i; j <= temp / 2.0; ++j)
			{
				
				// Found a factor
				if (temp % j == 0)
				{
					
					// Last factor is too small
					if (temp / j < j)
						continue;
					
					// Save factors
					res.add();
					res.last().length(3);
					res.last()[0] = i;
					res.last()[1] = j;
					res.last()[2] = temp / j;
				}
			}
		}
	}
	
	// Return result
	return res;
}



/* inline void Num<T>::decimalTofraction(int* res, T number, double tol)
 *
 * Convert number to fraction
 */

template <class T>
inline void Num<T>::decimalTofraction(int* res, T number, double tol)
{
	
	// Convert to fraction
	T mult = number < 0 ? -1 : 1;
	number = Num<T>::abs(number);
	T Z = number;
    T denom = 1;
    T prevDenom = 0;
    T num = Z;
    T temp;
	if (fabs(Z - (int)Z) > tol/10)
	{
		do
	    {
			Z = 1.0/(Z - (int)Z);
	        temp = denom;
	        denom = denom * (int)Z + prevDenom;
	        prevDenom = temp;
	        num = (int)(number * denom + 0.5);
	    }
	    while ((fabs(number - num/denom) > tol) && (Z != (int)Z));
	}

	// Save result
	res[0] = (int) (mult * Num<T>::round(num, 1));
	res[1] = (int) Num<T>::round(denom, 1);
}



/* inline T Num<T>::dot(int length, const T* vec1, const T* vec2)
 *
 * Calculate dot product between two vectors
 */

template <class T>
inline T Num<T>::dot(int length, const T* vec1, const T* vec2)
{
	T res = 0;
	for (int i = 0; i < length; ++i)
		res += vec1[i] * vec2[i];
	return res;
}



/* inline T Num<T>::magnitude(int length, const T* vec)
 *
 * Calculate magnitude of a vector
 */

template <class T>
inline T Num<T>::magnitude(int length, const T* vec)
{
	return (T) sqrt(dot(length, vec, vec));
}



/* inline T Num<T>::angle(const T* vec1, const T* vec2)
 *
 * Return the angle between two vectors
 */

template <class T>
inline T Num<T>::angle(const T* vec1, const T* vec2)
{
	T temp = dot(3, vec1, vec2) / (magnitude(3, vec1) * magnitude(3, vec2));
	if ((temp > 1) || (temp < -1))
        return 0;
    return acos(temp);
}



/* inline T Num<T>::volume(const T* vec1, const T* vec2, const T* vec3)
 *
 * Calculate volume enclosed by three vectors
 */

template <class T>
inline T Num<T>::volume(const T* vec1, const T* vec2, const T* vec3)
{
	return vec1[0] * (vec2[1] * vec3[2] - vec2[2] * vec3[1]) + \
		   vec1[1] * (vec2[2] * vec3[0] - vec2[0] * vec3[2]) + \
		   vec1[2] * (vec2[0] * vec3[1] - vec2[1] * vec3[0]);
}



/* inline T Num<T>::integrate(T (*fun)(T), T min, T max, IntegrationMethod method)
 *
 * Numerical integration
 */

template <class T>
template <class Tclass>
inline T Num<T>::integrate(Functor<Tclass>& fun, T min, T max, IntegrationMethod method)
{
	
	// Loop until converged
	int i;
	T curInt;
	T prevInt;
	bool first = true;
	double tol = 1e-4;
	Matrix points;
        int order;
	for (order = 6; order <= 48; order += 3)
	{
		
		// Save previous result
		prevInt = curInt;
		
		// Generate points
		points.size(order, 2);
		if (method == GAUSSLEGENDRE)
			GaussLegendrePoints(min, max, points);
		
		// Evaluate integral
		curInt = 0;
		for (i = 0; i < order; ++i)
			curInt += points(i, 1) * fun(points(i, 0));
		curInt *= (max - min) / 2;
		
		// Break if converged
		if (first)
			first = false;
		else if (Num<T>::abs(curInt - prevInt) / curInt < tol)
			break;
	}
        
	// Return integral
	return curInt;
}



/* inline void Num<T>::GaussLegendrePoints(T min, T max, Matrix& points)
 *
 * Generate points for Gauss-Legendre quadrature integration
 */

template <class T>
inline void Num<T>::GaussLegendrePoints(T min, T max, Matrix& points)
{
	
	// Generate coefficients
	int i, j;
	int order = points.numRows();
	Matrix coeffs(order + 1);
	coeffs.fill(0);
	coeffs(0, 0) = 1;
	coeffs(1, 1) = 1;
	for (i = 2; i <= order; ++i)
	{
		coeffs(i, 0) = -(i - 1) * coeffs(i-2, 0) / i;
		for (j = 1; j <= order; ++j)
			coeffs(i, j) = ((2*i - 1) * coeffs(i-1, j-1) - (i - 1)*coeffs(i-2, j)) / i;
	}
	
	// Generate points and weights
	T x;
	T x1;
	for (i = 1; i <= order; ++i)
	{
		x = cos(Constants::pi * (i - 0.25) / (order + 0.5));
		do
		{
			x1 = x;
			x -= Legendre(coeffs, order, x) * (x*x - 1) / order /\
			 	(x * Legendre(coeffs, order, x) - Legendre(coeffs, order - 1, x));
		} while (abs(x - x1) > 1e-8);
		points(i-1, 0) = x;
		x1 = order*(x*Legendre(coeffs, order, x) - Legendre(coeffs, order - 1, x)) / (x*x - 1);
		points(i-1, 1) = 2 / ((1 - x*x) * x1*x1);
	}
	
	// Convert points from -1:1 to min:max
	for (i = 0; i < order; ++i)
		points(i, 0) = points(i, 0) * (max - min) / 2 + (max + min) / 2;
}



/* inline T Num<T>::Legendre(const Matrix<T>& coeffs, int order, T point)
 *
 * Generate value of Legendre polynomial at point
 */

template <class T>
inline T Num<T>::Legendre(const Matrix& coeffs, int order, T point)
{
	T res = coeffs(order, order);
	for (int i = order; i > 0; --i)
		res = res * point + coeffs(order, i-1);
	return res;
}



// =====================================================================================================================
// Solve
// =====================================================================================================================

// Static member variables
template <class Tclass>
Functor<Tclass>* Solve<Tclass>::_functorFun;



/* inline double Solve<Tclass>::findRoot(Functor<Tclass>& fun, double target, double convergeArg, double convergeRes,
 *		double min, double max, RootMethod method, double& argRes)
 *
 * Find root of a function
 */

template <class Tclass>
inline double Solve<Tclass>::findRoot(Functor<Tclass>& fun, double target, double convergeArg, double convergeRes, \
	double min, double max)
{
	double argRes;
	return findRoot(fun, target, convergeArg, convergeRes, min, max, BRENT, argRes);
}

template <class Tclass>
inline double Solve<Tclass>::findRoot(Functor<Tclass>& fun, double target, double convergeArg, double convergeRes, \
	double min, double max, RootMethod method)
{
	double argRes;
	return findRoot(fun, target, convergeArg, convergeRes, min, max, method, argRes);
}

template <class Tclass>
inline double Solve<Tclass>::findRoot(Functor<Tclass>& fun, double target, double convergeArg, double convergeRes, \
	double min, double max, double& argRes)
{
	return findRoot(fun, target, convergeArg, convergeRes, min, max, BRENT, argRes);
}

template <class Tclass>
inline double Solve<Tclass>::findRoot(Functor<Tclass>& fun, double target, double convergeArg, double convergeRes, \
	double min, double max, RootMethod method, double& argRes)
{
	switch (method)
	{
		case BRENT:
			return brent(fun, target, convergeArg, convergeRes, min, max, argRes);
		case BISECTION:
			return bisection(fun, target, convergeArg, convergeRes, min, max, argRes);
		default:
			return min;
	}
}



/* inline double Solve<Tclass>::brent(Functor<Tclass>& fun, double target, double convergeArg, double convergeRes, 
 *		double min, double max, double& argRes)
 *
 * Brent's method to find function root
 */

template <class Tclass>
inline double Solve<Tclass>::brent(Functor<Tclass>& fun, double target, double convergeArg, double convergeRes, \
	double min, double max, double& argRes)
{
	
	// Set a and b
    double a = min;
    double b = max;
    double fa = fun(a) - target;
    double fb = fun(b) - target;
    if (Num<double>::abs(fa) < Num<double>::abs(fb))
		Num<double>::swap(a, b);
    
    // Loop until max iterations is reach
    bool mflag = true;
    double c = a;
    double fc;
    double s;
    double fs;
    double d = 0;
	double res = 0;
	double lastRes;
	int count = 0;
	int maxCount = 10000;
    for (count = 0; count < maxCount; ++count)
    {
        
        // Get function values
        fa = fun(a) - target;
        fb = fun(b) - target;
        fc = fun(c) - target;

        // Inverse quadratic interpolation
		if ((Num<double>::abs(fa - fc) > convergeRes) && (Num<double>::abs(fb - fc) > convergeRes))
			s = (a * fb * fc) / ((fa - fb) * (fa - fc)) + (b * fa * fc) / ((fb - fa) * (fb - fc)) + \
				(c * fa * fc) / ((fc - fa) * (fc - fb));
        
        // Secant rule
		else
			s = b - fb * (b - a) / (fb - fa);
		
        // Test if bisection should be used
        if (((s < (3*a + b) / 4) || (s > b)) || \
            ((mflag)  && (Num<double>::abs(s - b) >= Num<double>::abs(b - c)/2)) || \
            ((!mflag) && (Num<double>::abs(s - b) >= Num<double>::abs(c - d)/2)) || \
            ((mflag)  && (Num<double>::abs(b - c) < convergeArg)) || \
            ((!mflag) && (Num<double>::abs(c - d) < convergeArg)))
        {
            s = (a + b)/2;
            mflag = false;
        }
        else
            mflag = false;

        // Set values
		lastRes = res;
		argRes = s;
		fs = fun(s) - target;
		res = fs + target;
        d = c;
        c = b;
        if (fa*fs < 0)
            b = s;
        else
            a = s;
        if (Num<double>::abs(fa) < Num<double>::abs(fb))
			Num<double>::swap(a, b);
        
        // Break if converged
		//if ((count) && (Num<double>::abs(lastRes - res) <= convergeRes))
		//	break;
		if (Num<double>::abs(fb) <= convergeRes)
        {
            res = fb + target;
			argRes = b;
            break;
        }
        if (Num<double>::abs(fs) <= convergeRes)
            break;
    }

	// Return result
	return res;
}



/* inline double Solve<Tclass>::bisection(Functor<Tclass>& fun, double target, double convergeArg, double convergeRes,
 *		double min, double max, double& argRes)
 *
 * Bisection method to find function root
 */

template <class Tclass>
inline double Solve<Tclass>::bisection(Functor<Tclass>& fun, double target, double convergeArg, double convergeRes, \
	double min, double max, double& argRes)
{
	
	// Variables to store values
	double mid;
	double fmid;
	double fmin = fun(min) - target;
	double fmax = fun(max) - target;
	
	// Points do not bound root
	if (((fmin > 0) && (fmax > 0)) || ((fmin < 0) && (fmax < 0)))
		return min;
	
	// Swap min and max if needed
	if (fmin > fmax)
	{
		Num<double>::swap(min, max);
		Num<double>::swap(fmin, fmax);
	}
	
	// Loop until solved
	int maxCount = 10000;
	for (int count = 0; count < maxCount; ++count)
	{
		
		// Create new point
		mid = (min + max) / 2;
		fmid = fun(mid) - target;
		
		// Check if converged
		if (Num<double>::abs(fmid) <= convergeRes)
			break;
		if (mid - min <= convergeArg)
			break;
		
		// Replace min or max
		if (((fmid > 0) && (fmin > 0)) || ((fmid < 0) && (fmin < 0)))
		{
			min = mid;
			fmin = fmid;
		}
		else
		{
			max = mid;
			fmax = fmid;
		}
	}
	
	// Return result
	argRes = mid;
	return fmid + target;
}



/* inline double Solve<Tclass>::minimize(Functor<Tclass>& fun, double convergeRes, double initial, double step,
 *		double& argRes)
 *
 * Minimize a function of one variable
 */

template <class Tclass>
inline double Solve<Tclass>::minimize(Functor<Tclass>& fun, double convergeRes, double initial, double step)
{
	double argRes;
	return minimize(fun, convergeRes, initial, step, argRes);
}

template <class Tclass>
inline double Solve<Tclass>::minimize(Functor<Tclass>& fun, double convergeRes, double initial, double step, \
	double& argRes)
{
	
	// Vector variables
	Vector initVec(1);
	Vector stepVec(1);
	Vector argResVec(1);
	initVec[0] = initial;
	stepVec[0] = step;
	
	// Functor with vector
	_functorFun = &fun;
	Functor<Solve<Tclass> > funVec(&Solve<Tclass>::NelderHelper);
	
	// Evaluate function
	double res = Solve<Solve<Tclass> >::NelderMead(funVec, convergeRes, initVec, stepVec, argResVec, 1);
	
	// Return result
	argRes = argResVec[0];
	return res;
}



/* inline double Solve<Tclass>::maximize(Functor<Tclass>& fun, double convergeRes, double initial, double step,
 *		double& argRes)
 *
 * Maximize a function of one variable
 */

template <class Tclass>
inline double Solve<Tclass>::maximize(Functor<Tclass>& fun, double convergeRes, double initial, double step)
{
	double argRes;
	return maximize(fun, convergeRes, initial, step, argRes);
}

template <class Tclass>
inline double Solve<Tclass>::maximize(Functor<Tclass>& fun, double convergeRes, double initial, double step, \
	double& argRes)
{
	
	// Vector variables
	Vector initVec(1);
	Vector stepVec(1);
	Vector argResVec(1);
	initVec[0] = initial;
	stepVec[0] = step;
	
	// Functor with vector
	_functorFun = &fun;
	Functor<Solve<Tclass> > funVec(&Solve<Tclass>::NelderHelper);
	
	// Evaluate function
	double res = Solve<Solve<Tclass> >::NelderMead(funVec, convergeRes, initVec, stepVec, argResVec, -1);
	
	// Return result
	argRes = argResVec[0];
	return res;
}



/* inline double Solve<Tclass>::minimize(Functor<Tclass>& fun, double convergeRes, const Vector& initial,
 *		const Vector& step, Vector& argRes)
 *
 * Minimize a function of multiple variable
 */

template <class Tclass>
inline double Solve<Tclass>::minimize(Functor<Tclass>& fun, double convergeRes, const Vector& initial, \
	const Vector& step)
{
	Vector argRes;
	return NelderMead(fun, convergeRes, initial, step, argRes, 1);
}

template <class Tclass>
inline double Solve<Tclass>::minimize(Functor<Tclass>& fun, double convergeRes, const Vector& initial, \
	const Vector& step, Vector& argRes)
{
	return NelderMead(fun, convergeRes, initial, step, argRes, 1);
}


/* inline double Solve<Tclass>::maximize(Functor<Tclass>& fun, double convergeRes, const Vector& initial,
 *		const Vector& step, Vector& argRes)
 *
 * Maximize a function of multiple variable
 */

template <class Tclass>
inline double Solve<Tclass>::maximize(Functor<Tclass>& fun, double convergeRes, const Vector& initial, \
	const Vector& step)
{
	Vector argRes;
	return NelderMead(fun, convergeRes, initial, step, argRes, -1);
}

template <class Tclass>
inline double Solve<Tclass>::maximize(Functor<Tclass>& fun, double convergeRes, const Vector& initial, \
	const Vector& step, Vector& argRes)
{
	return NelderMead(fun, convergeRes, initial, step, argRes, -1);
}



/* inline double Solve<Tclass>::NelderMead(Functor<Tclass>& fun, double convergeRes, const Vector& initial,
 *		const Vector& steps, Vector& argRes, int mult)
 *
 * Nelder-Mead simplex method
 */

template <class Tclass>
inline double Solve<Tclass>::NelderMead(Functor<Tclass>& fun, double convergeRes, const Vector& initial, \
	const Vector& steps, Vector& argRes, int mult)
{
	
	// Return if number of parameters does not match steps
	if (initial.length() != steps.length())
	{
		argRes = Vector(0);
		return 0;
	}
	
	// Variables to store vertices
	OList<Vector> vertices(initial.length() + 1);
	List<double>  values (vertices.length());
	
	// Initialize vertices
	int i;
	for (i = 0; i < vertices.length(); ++i)
	{
		vertices[i] = initial;
		if (i > 0)
			vertices[i][i-1] += steps[i-1];
		values[i] = mult * fun(vertices[i]);
	}
	
	// Scaling factors
	double cReflect = 1.0;
	double cExpand = 1.0 + 2.0 / initial.length();
	double cContract = 0.75 - 1.0 / (2*initial.length());
	double cReduce = 1.0 - 1.0/initial.length();
	
	// Loop until converged
	int j;
	int maxCount = 10000;
	double valueReflected;
	double valueExpanded;
	double valueContracted;
	Vector center(initial);
	Vector vertexReflected;
	Vector vertexExpanded;
	Vector vertexContracted;
	for (int count = 0; count < maxCount; ++count)
	{
		
		// Sort values
		sortVertices(vertices, values, 0, vertices.length() - 1);
		
		// Test for convergence
		if ((values.last() - values[0]) / (1 + Num<double>::abs(values.last())) < convergeRes)
			break;
		
		// Get centroid
		for (i = 0; i < vertices[0].length(); ++i)
		{
			center[i] = 0;
			for (j = 0; j < vertices.length() - 1; ++j)
				center[i] += vertices[j][i];
			center[i] /= vertices.length() - 1;
		}
		
		// Reflection
		vertexReflected = center;
		vertexReflected += (center - vertices.last()) * cReflect;
		valueReflected = mult * fun(vertexReflected);
		if ((valueReflected >= values[0]) && (valueReflected < values[values.length() - 2]))
		{
			values.last() = valueReflected;
			vertices.last() = vertexReflected;
			continue;
		}
		
		// Expansion
		if (valueReflected < values[0])
		{
			vertexExpanded = center;
			vertexExpanded += (vertexReflected - center) * cExpand;
			valueExpanded = mult * fun(vertexExpanded);
			if (valueExpanded < valueReflected)
			{
				values.last() = valueExpanded;
				vertices.last() = vertexExpanded;
			}
			else
			{
				values.last() = valueReflected;
				vertices.last() = vertexReflected;
			}
			continue;
		}
		
		// Outside contraction
		if (valueReflected < values.last())
		{
			vertexContracted = center;
			vertexContracted += (vertexReflected - center) * cContract;
			valueContracted = mult * fun(vertexContracted);
			if (valueContracted <= valueReflected)
			{
				values.last() = valueContracted;
				vertices.last() = vertexContracted;
			}
		}
		
		// Inside contraction
		else
		{
			vertexContracted = center;
			vertexContracted -= (vertexReflected - center) * cContract;
			valueContracted = mult * fun(vertexContracted);
			if (valueContracted < values.last())
			{
				values.last() = valueContracted;
				vertices.last() = vertexContracted;
			}
		}
		
		// Reduction
		for (i = 1; i < vertices.length(); ++i)
		{
			vertices[i] = vertices[0] + (vertices[i] - vertices[0]) * cReduce;
			values[i] = mult * fun(vertices[i]);
		}
	}
	
	// Return result
	argRes = vertices[0];
	return values[0] / mult;
}



/* inline void Solve<Tclass>::sortVertices(OList<Vector>& vertices, List<double>& values, int left, int right)
 *
 * Sort values during Nelder-Mead method using quicksort
 */

template <class Tclass>
inline void Solve<Tclass>::sortVertices(OList<Vector>& vertices, List<double>& values, int left, int right)
{
	
	// Current partition has no width
	if (left >= right)
		return;
	
	// Choose pivot index
	int pivotIndex = (left + right) / 2;
	double pivot = values[pivotIndex];
	
	// Move pivot to end
	values.swap(pivotIndex, right);
	vertices.swap(pivotIndex, right);
	
	// Iterate through current partition
	int newPivotIndex = left;
	for (int i = left; i < right; ++i)
	{
		if (values[i] < pivot)
		{
			values.swap(i, newPivotIndex);
			vertices.swap(i, newPivotIndex);
			newPivotIndex++;
		}
	}
	
	// Move pivot to final position
	values.swap(newPivotIndex, right);
	vertices.swap(newPivotIndex, right);
	
	// Recursive calls to next sorts
	sortVertices(vertices, values, left, newPivotIndex - 1);
	sortVertices(vertices, values, newPivotIndex + 1, right);
}



// =====================================================================================================================
// Fit
// =====================================================================================================================

/* inline Vector Fit::polynomial(const OList<List<double> >& data, int order0, int numTerms)
 *
 * Fit to polynomial expansion
 */

inline Vector Fit::polynomial(const OList<List<double> >& data, int order0, int numTerms)
{
	
	// Build M and y for Mx = y
	int i, j, k;
	List<double> powers(2*numTerms);
	Vector y(numTerms);
	Matrix M(numTerms);
	y.fill(0);
	M.fill(0);
	for (i = 0; i < data.length(); ++i)
	{
		if (numTerms != 0)
			powers[0] = pow(data[i][0], order0);
		for (j = 1; j < 2*numTerms; ++j)
			powers[j] = powers[j-1] * data[i][0];
		for (j = 0; j < numTerms; ++j)
			y[j] += powers[j] * data[i][1];
		for (j = 0; j < numTerms; ++j)
		{
			for (k = 0; k < numTerms; ++k)
				M(j, k) += powers[j + k];
		}
	}
	
	// Get fitting parameters
	return M.solve(y);
}



/* inline Vector Fit::LM(const OList<List<double> >& data, Functor<Tclass>& fun, VectorFunctor<Tclass>& deriv,
 *		const Vector& initial, double tol)
 *
 * Levenberg–Marquardt algorithm to fit general functions
 * At each step, solve for delta in
 *		[Transpose(J)*J + lambda*diag(Transpose(J)*J)]*delta = Transpose(J)*[y - f(beta)]
 *		J is the Jacobian matrix where there are as many rows as data points
 *			Each column gives the derivative of the fitting function with respect to a fitting parameter
 *		lambda is the damping factor
 *		y is the list of all data points
 *		f(beta) is the value of the fitting function with a fit set beta evaluated at all data points
 */

template <class Tclass>
inline Vector Fit::LM(const OList<List<double> >& data, Functor<Tclass>& fun, VectorFunctor<Tclass>& deriv, \
	const Vector& initial, double tol)
{
	
	// Variable to store result
	Vector params = initial;
	Vector paramsNew = initial;
	
	// Storage variables
	Matrix jacobian(data.length(), params.length());
	Matrix jacobianTranspose(params.length(), data.length());
	Matrix matrix(data.length(), data.length());
	
	// Save initial residual and magnitude
	Vector residual(data.length());
	LMGetResidual(data, fun, params, residual);
	double residualMag = residual.magnitude();
	
	// Loop until converged
	double residualMagNew;
	double damper = 100;
	double damperUp = 10;
	double damperDown = 10;
	double cutoff = 1 - tol;
	Vector residualNew(data.length());
	for (int i = 0; i < 1000; ++i)
	{
		
		// Get the new parameters
		LMBuildJacobian(data, deriv, params, jacobian, jacobianTranspose);
		LMSetMatrix(jacobian, jacobianTranspose, damper, matrix);
		paramsNew = params + matrix.solve(jacobianTranspose*residual);
		
		// Get the new residual
		LMGetResidual(data, fun, paramsNew, residualNew);
		residualMagNew = residualNew.magnitude();
		
		// Found a new best solution
		if (residualMagNew < residualMag)
		{
			
			// Save parameters
			params = paramsNew;
			
			// Check if converged
			if (residualMagNew / residualMag > cutoff)
				break;
			
			// Save values
			residual = residualNew;
			residualMag = residualMagNew;
			damper /= damperDown;
		}
		
		// Solution is worse than previous
		else
			damper *= damperUp;
	}
	
	// Return result
	return params;
}



/* inline void Fit::LMBuildJacobian(const OList<List<double> >& data, VectorFunctor<Tclass>& deriv, Vector& params,
 *		Matrix& jacobian, Matrix& jacobianTranspose)
 *
 * Build the Jacobian matrix and its transpose for Levenberg–Marquardt algorithm
 */

template <class Tclass>
inline void Fit::LMBuildJacobian(const OList<List<double> >& data, VectorFunctor<Tclass>& deriv, Vector& params, \
	Matrix& jacobian, Matrix& jacobianTranspose)
{
	int i, j;
	Vector derivs;
	for (i = 0; i < data.length(); ++i)
	{
		derivs = deriv(params, data[i][0]);
		for (j = 0; j < derivs.length(); ++j)
			jacobian(i, j) = jacobianTranspose(j, i) = derivs[j];
	}
}



/* inline void Fit::LMSetMatrix(Matrix& jacobian, Matrix& jacobianTranspose, double damper, Matrix& matrix)
 *
 * Set the matrix values for Levenberg–Marquardt algorithm
 */

inline void Fit::LMSetMatrix(Matrix& jacobian, Matrix& jacobianTranspose, double damper, Matrix& matrix)
{
	matrix = jacobianTranspose * jacobian;
	for (int i = 0; i < matrix.numRows(); ++i)
		matrix(i, i) += matrix(i, i) * damper;
}



/* inline void Fit::LMGetResidual(const OList<List<double> >& data, Functor<Tclass>& fun, Vector& params, Vector& residuals)
 *
 * Get residual in Levenberg-Marquardt algorithm
 */

template <class Tclass>
inline void Fit::LMGetResidual(const OList<List<double> >& data, Functor<Tclass>& fun, Vector& params, Vector& residuals)
{
	for (int i = 0; i < data.length(); ++i)
		residuals[i] = data[i][1] - fun(params, data[i][0]);
}



/* inline Vector Fit::NM(const OList<List<double> >& data, Functor<Tclass>& fun, const Vector& initial, 
 *		const Vector& steps, double tol)
 *
 * Fit funtion to data points using Nelder-Meads minimization of R^2
 */

template <class Tclass>
inline Vector Fit::NM(const OList<List<double> >& data, Functor<Tclass>& fun, const Vector& initial, \
	const Vector& steps, double tol)
{
	
	// Create functor object
	NMHelp<Tclass> help;
	help._functor = &fun;
	help._NMData = &data;
	Functor<NMHelp<Tclass> > functor(&help, &help.NMRSquared);
	
	// Return fitting parameters
	return Solve<Tclass>::minimize(functor, tol, initial, steps);
}



/* inline double Fit::NMHelp<Tclass>::NMRSquared(const Vector& params)
 *
 * R squared value for Nelder-Meads method 
 */

template <class Tclass>
inline double Fit::NMHelp<Tclass>::NMRSquared(const Vector& params)
{
	double temp;
	double res = 0;
	for (int i = 0; i < _NMData->length(); ++i)
	{
		temp = (*_NMData)[i][1] - (*_functor)(params, (*_NMData)[i][0]);
		res += temp*temp;
	}
	return res;
}



/* inline double Fit::getRSquared(const OList<List<double> >& data, Functor<Tclass>& fun, Vector& params)
 *
 * Get R squared value for fit
 */

template <class Tclass>
inline double Fit::getRSquared(const OList<List<double> >& data, Functor<Tclass>& fun, Vector& params)
{
	double temp;
	double res = 0;
	for (int i = 0; i < data.length(); ++i)
	{
		temp = data[i][1] - fun(params, data[i][0]);
		res += temp*temp;
	}
	return res;
}



#endif
