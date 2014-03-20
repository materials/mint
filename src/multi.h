/* multi.h -- MPI functions
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef MULTI_H
#define MULTI_H



#ifdef MINT_MPI
	#include <mpi.h>
#endif
#include "text.h"
#include "num.h"



// Class for MPI functions
class Multi
{
	
	// Variables
	static int _rank;
	static int _worldSize;
	static bool _mpiOn;
	
	// Variables to control job
	static bool _setExitFun;
	static bool _runningFun;
	static int _jobSize;
	
	// Exit function to control jobs
	static void exitfun()	{ if (_runningFun) throw 1; }
	
public:
	
	// Setup functions
	static void initialize(int argc, char** argv);
	static void finalize();
	static void jobSize(int input)	{ _jobSize = input; }
	
	// General functions
	static void barrier();
	static bool safeCall(void (&function)(), const char* stdoutFile = 0, const char* stderrFile = 0);
	static void external(const Word& exe, const char* flags = 0, const char* stdoutFile = 0);
	
	// Send single values
	static void send(             bool& value, int root, int dest)	{ send(&value, 1, root, dest); }
	static void send(              int& value, int root, int dest)	{ send(&value, 1, root, dest); }
	static void send(unsigned long int& value, int root, int dest)	{ send(&value, 1, root, dest); }
	static void send(           double& value, int root, int dest)	{ send(&value, 1, root, dest); }
	static void send(             char& value, int root, int dest)	{ send(&value, 1, root, dest); }
	
	// Send arrays
	static void send(             bool* array, int length, int root, int dest);
	static void send(              int* array, int length, int root, int dest);
	static void send(unsigned long int* array, int length, int root, int dest);
	static void send(           double* array, int length, int root, int dest);
	static void send(             char* array, int length, int root, int dest);
	
	// Send special values
	static void send(Vector3D& vector, int root, int dest);
	
	// Broadcast single values
	static void broadcast(             bool& value, int root)	{ broadcast(&value, 1, root); }
	static void broadcast(              int& value, int root)	{ broadcast(&value, 1, root); }
	static void broadcast(unsigned long int& value, int root)	{ broadcast(&value, 1, root); }
	static void broadcast(           double& value, int root)	{ broadcast(&value, 1, root); }
	static void broadcast(             char& value, int root)	{ broadcast(&value, 1, root); }
	
	// Broadcast arrays
	static void broadcast(             bool* array, int length, int root);
	static void broadcast(              int* array, int length, int root);
	static void broadcast(unsigned long int* array, int length, int root);
	static void broadcast(           double* array, int length, int root);
	static void broadcast(             char* array, int length, int root);
	
	// Broadcast special values
	static void broadcast(Vector3D& vector, int root);
	static void broadcast(  Vector& vector, int root);
	
	// Access functions
	static int rank()		{ return _rank; }
	static int worldSize()	{ return _worldSize; }
	static bool mpiOn()		{ return _mpiOn; }
};



/* inline void Multi::barrier()
 *
 * Add barrier
 */

inline void Multi::barrier()
{
	#ifdef MINT_MPI
		MPI_Barrier(MPI_COMM_WORLD);
	#endif
}



/* inline void Multi::send(TYPE* array, int length, int root, int dest)
 *
 * Send an array of variables
 */

inline void Multi::send(bool* array, int length, int root, int dest)
{
	#ifdef MINT_MPI
		if (root == dest)
			return;
		if (_rank == root)
			MPI_Send(array, length, MPI_BYTE, dest, dest, MPI_COMM_WORLD);
		if (_rank == dest)
			MPI_Recv(array, length, MPI_BYTE, root, dest, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	#endif
}

inline void Multi::send(int* array, int length, int root, int dest)
{
	#ifdef MINT_MPI
		if (root == dest)
			return;
		if (_rank == root)
			MPI_Send(array, length, MPI_INT, dest, dest, MPI_COMM_WORLD);
		if (_rank == dest)
			MPI_Recv(array, length, MPI_INT, root, dest, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	#endif
}

inline void Multi::send(unsigned long int* array, int length, int root, int dest)
{
	#ifdef MINT_MPI
		if (root == dest)
			return;
		if (_rank == root)
			MPI_Send(array, length, MPI_UNSIGNED_LONG, dest, dest, MPI_COMM_WORLD);
		if (_rank == dest)
			MPI_Recv(array, length, MPI_UNSIGNED_LONG, root, dest, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	#endif
}

inline void Multi::send(double* array, int length, int root, int dest)
{
	#ifdef MINT_MPI
		if (root == dest)
			return;
		if (_rank == root)
			MPI_Send(array, length, MPI_DOUBLE, dest, dest, MPI_COMM_WORLD);
		if (_rank == dest)
			MPI_Recv(array, length, MPI_DOUBLE, root, dest, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	#endif
}

inline void Multi::send(char* array, int length, int root, int dest)
{
	#ifdef MINT_MPI
		if (root == dest)
			return;
		if (_rank == root)
			MPI_Send(array, length, MPI_CHAR, dest, dest, MPI_COMM_WORLD);
		if (_rank == dest)
			MPI_Recv(array, length, MPI_CHAR, root, dest, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	#endif
}

inline void Multi::send(Vector3D& vector, int root, int dest)
{
	#ifdef MINT_MPI
		if (root == dest)
			return;
		if (_rank == root)
			MPI_Send(vector._vector, 3, MPI_DOUBLE, dest, dest, MPI_COMM_WORLD);
		if (_rank == dest)
			MPI_Recv(vector._vector, 3, MPI_DOUBLE, root, dest, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	#endif
}



/* inline void Multi::broadcast(TYPE* array, int length, int root)
 *
 * Broadcast an array of variables
 */

inline void Multi::broadcast(bool* array, int length, int root)
{
	#ifdef MINT_MPI
		MPI_Bcast(array, length, MPI_BYTE, root, MPI_COMM_WORLD);
	#endif
}

inline void Multi::broadcast(int* array, int length, int root)
{
	#ifdef MINT_MPI
		MPI_Bcast(array, length, MPI_INT, root, MPI_COMM_WORLD);
	#endif
}

inline void Multi::broadcast(unsigned long int* array, int length, int root)
{
	#ifdef MINT_MPI
		MPI_Bcast(array, length, MPI_UNSIGNED_LONG, root, MPI_COMM_WORLD);
	#endif
}

inline void Multi::broadcast(double* array, int length, int root)
{
	#ifdef MINT_MPI
		MPI_Bcast(array, length, MPI_DOUBLE, root, MPI_COMM_WORLD);
	#endif
}

inline void Multi::broadcast(char* array, int length, int root)
{
	#ifdef MINT_MPI
		MPI_Bcast(array, length, MPI_CHAR, root, MPI_COMM_WORLD);
	#endif
}

inline void Multi::broadcast(Vector3D& vector, int root)
{
	#ifdef MINT_MPI
		MPI_Bcast(vector._vector, 3, MPI_DOUBLE, root, MPI_COMM_WORLD);
	#endif
}

inline void Multi::broadcast(Vector& vector, int root)
{
	#ifdef MINT_MPI
		int length = vector._length;
		MPI_Bcast(&length, 1, MPI_INT, root, MPI_COMM_WORLD);
		vector.length(length);
		MPI_Bcast(vector._vector, length, MPI_DOUBLE, root, MPI_COMM_WORLD);
	#endif
}



#endif
