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



#ifndef LIST_H
#define LIST_H



#include <cmath>



// Sorting methods
enum SortMethod {QUICKSORT};



// Class to store a list of objects
template <class T>
class OList
{
	
	// Buffer size
	static int _buffer;

	// Variables
	int _actualLength;
	int _length;
	T** _list;
	
	// Functions
	void initialize();
	
public:
	
	// Constructors
	OList() 							{ initialize(); }
	OList(const T& copy)				{ initialize(); *this = copy; }
	OList(const OList<T>& copy)			{ initialize(); *this = copy; }
	OList(int inLength)					{ initialize(); length(inLength); }
	OList(int inLength, const T& copy)	{ initialize(); length(inLength); fill(copy); }
	
	// Destructor
	~OList() { clear(); }
	
	// Delete data
	void clear();
	void remove(int index);
	
	// Setup functions
	void length(int inLength);
	void add()			{ length(_length + 1); }
	void add(int num)	{ length(_length + num); }
	void fill(const T& input);
	
	// General functions
	bool operator== (const OList<T>& rhs) const;
	void swap(int index1, int index2);
	
	// Assignment functions
	OList<T>& operator=  (const T& rhs);
	OList<T>& operator=  (const OList<T>& rhs);
	OList<T>  operator+  (const T& rhs) const;
	OList<T>  operator+  (const OList<T>& rhs) const;
	OList<T>& operator+= (const T& rhs);
	OList<T>& operator+= (const OList<T>& rhs);
	
	// Access functions
	int length() const				{ return _length; }
	T& operator[] (int index) const { return *(_list[index]); }
	T& last() const					{ return *(_list[_length - 1]); }
	
	// Higher dimension lists
	typedef OList<T> D1;
	typedef OList<OList<T> > D2;
	typedef OList<OList<OList<T> > > D3;
	typedef OList<OList<OList<OList<T> > > > D4;
};



// Class to store a simple list
template <class T>
class List
{
	
private:
	
	// Buffer size
	static int _buffer;

	// Variables
	int _actualLength;
	int _length;
	T* _list;
	
	// Functions
	void initialize();
	void quicksort(int left, int right);
	
public:
	
	// Constructors
	List() 								{ initialize(); }
	List(const List<T>& copy)			{ initialize(); *this = copy; }
	List(int inLength)					{ initialize(); length(inLength); }
	List(int inLength, const T* array)	{ initialize(); set(inLength, array); }
	List(int inLength, T value)			{ initialize(); length(inLength); fill(value); }
	
	// Destructor
	~List() { clear(); }
	
	// Delete data
	void clear();
	void remove(int index);
	
	// Setup functions
	void length(int inLength);
	void add()							{ length(_length + 1); }
	void add(int num)					{ length(_length + num); }
	void fill(const T& input);
	
	// General functions
	bool operator== (const List<T>& rhs) const;
	void swap(int index1, int index2);
	
	// Assignment functions
	void set(int num, const T* array);
	List<T>& operator=  (const T& rhs);
	List<T>& operator=  (const List<T>& rhs);
	List<T>  operator+  (const T& rhs) const;
	List<T>  operator+  (const List<T>& rhs) const;
	List<T>& operator+= (const T& rhs);
	List<T>& operator+= (const List<T>& rhs);
	
	// Complex functions
	void sort(SortMethod method = QUICKSORT);
	
	// Access functions
	int length() const			{ return _length; }
	T& operator[] (int index)	{ return _list[index]; }
	T& last()					{ return _list[_length - 1]; }
	const T* array() const		{ return _list; }

	// Constant access functions
	const T& operator[] (int index) const	{ return _list[index]; }
	const T& last()	const					{ return _list[_length - 1]; }
	
	// Higher dimension lists
	typedef List<T> D1;
	typedef OList<List<T> > D2;
	typedef OList<OList<List<T> > > D3;
	typedef OList<OList<OList<List<T> > > > D4;
};



// Class to store a node in a linked list
template <class T>
class Link
{
	
	// Variables
	mutable T _value;
	Link<T>* _prev;
	Link<T>* _next;
	
	// Functions
	void initialize()	{ _prev = 0; _next = 0; }
	
public:
	
	// Constructors
	Link()														{ initialize(); }
	Link(const T& input)										{ initialize(); _value = input; }
	Link(Link<T>* prevPtr, Link<T>* nextPtr)					{ _prev = prevPtr; _next = nextPtr; }
	Link(const T& input, Link<T>* prevPtr, Link<T>* nextPtr)	{ _value = input; _prev = prevPtr; _next = nextPtr; }	
	
	// Access functions
	T& value() const	{ return _value; }
	Link<T>*& prev()	{ return _prev; }
	Link<T>*& next() 	{ return _next; }
	
	// Friends
	template <class TLinked>
	friend class Linked;
};



// Class to store a linked list
template <class T>
class Linked
{
	
public:
	
	// Class to store iterator
	class iterator
	{
		Link<T>* _ptr;
	public:
		iterator()											{}
		iterator(const iterator& copy) : _ptr(copy._ptr)	{}
		iterator(Link<T>* input) : _ptr(input)				{}
		iterator& operator= (const iterator& rhs)			{ _ptr = rhs._ptr; return *this; }
		void operator++ ()									{ _ptr = _ptr->next(); }
		void operator++ (int)								{ _ptr = _ptr->next(); }
		void operator-- ()									{ _ptr = _ptr->prev(); }
		void operator-- (int)								{ _ptr = _ptr->prev(); }
		Link<T>* operator-> ()								{ return _ptr; }
		Link<T>* pointer()									{ return _ptr; }
		T& operator* ()										{ return _ptr->value(); }
		bool operator== (const iterator& rhs) const			{ return (_ptr == rhs._ptr); }
		bool operator== (List<T>* rhs) const				{ return (_ptr == rhs); }
		bool operator!= (const iterator& rhs) const			{ return (_ptr != rhs._ptr); }
		bool operator!= (List<T>* rhs) const				{ return (_ptr != rhs); }
		void operator+= (int num)							{ for (int i = 0; i < num; i++) _ptr = _ptr->next(); }
		void operator-= (int num)							{ for (int i = 0; i < num; i++) _ptr = _ptr->prev(); }
		iterator operator+ (int num)
		{
			iterator res(*this);
			for (int i = 0; i < num; i++)
				++res;
			return res;
		}
		iterator operator- (int num)
		{
			iterator res(*this);
			for (int i = 0; i < num; i++)
				--res;
			return res;
		}
	};
	
private:

	// Variables
	Link<T>* _begin;
	Link<T>* _last;
	iterator _end;
	int _length;
	
	// Functions
	void initialize() { _begin = 0; _last = 0; _length = 0; _end = iterator(0); }
	
public:
	
	// Constructor
	Linked() { initialize(); }
	Linked(const Linked<T>& rhs);
	
	// Destructor
	~Linked()	{ clear(); }
	
	// Clear all data
	void clear();
	
	// Setup functions
	Linked<T>& operator= (const Linked<T>& rhs);
	void add();
	void add(const T& value);
	void remove(iterator link);
	void operator+= (const T& value)	{ add(value); }
	void addBefore(iterator& it, const T& value);
	void addAfter(iterator& it, const T& value);
	
	// Acess functions
	int length() const			{ return _length; }
	T& operator[] (int index)	{ return *(iterator(begin) + index); }
	
	// Iterator
	iterator begin() const	{ return iterator(_begin); }
	iterator last() const	{ return iterator(_last); }
	iterator end() const	{ return _end; }
};



// =====================================================================================================================
// Object list
// =====================================================================================================================

// Buffer size for OList object
template <class T>
int OList<T>::_buffer = 4;



/* inline void OList<T>::initialize()
 *
 * Initialize OList object
 */

template <class T>
inline void OList<T>::initialize()
{
	_actualLength = 0;
	_length = 0;
	_list = 0;
}



/* inline void OList<T>::clear()
 *
 * Clear data in OList object
 */

template <class T>
inline void OList<T>::clear()
{
	for (int i = 0; i < _actualLength; ++i)
	{
		if (_list[i])
			delete _list[i];
	}
	if (_list)
		delete [] _list;
	initialize();
}



/* void OList<T>::remove(int index)
 *
 * Remove value from list at index
 */

template <class T>
void OList<T>::remove(int index)
{
	
	// Delete value
	if (_list[index])
		delete _list[index];
	
	// Copy over value to remove
	for (int i = index; i < _length-1; ++i)
		_list[i] = _list[i+1];
	
	// Set last to null pointer and reset length
	if (_length)
		_list[--_length] = 0;
}



/* void OList<T>::length(int inLength)
 *
 * Set the length of the list
 */

template <class T>
void OList<T>::length(int inLength)
{
	
	// Requested length is greater than or equal to actual length
	if (inLength > _actualLength)
	{
		
		// Figure out new buffered size
		int newLength = _buffer * (int) ceil((double) inLength / _buffer);
		
		// Create new list
		T** temp = new T* [newLength];
		
		// Copy old values
		int i;
		for (i = 0; i < _length; ++i)
			temp[i] = _list[i];
		
		// Delete old list
		if (_list)
			delete [] _list;
		
		// Create new values
		for (i = _length; i < newLength; ++i)
			temp[i] = new T;
		
		// Save list
		_list = temp;
		
		// Save new actual length
		_actualLength = newLength;
	}
	
	// Set length
	int i = _length;
	_length  = inLength;
	
	// Make sure that all values are allocated
	for (; i < _length; ++i)
	{
		if (!_list[i])
			_list[i] = new T;
	}
}



/* inline void OList<T>::fill(const T& input)
 *
 * Fill list with constant value
 */

template <class T>
inline void OList<T>::fill(const T& input)
{
	for (int i = 0; i < _length; ++i)
		_list[i]->operator=(input);
}



/* inline bool OList<T>::operator== (const OList<T>& rhs) const
 *
 * Test whether two lists are the same
 */

template <class T>
inline bool OList<T>::operator== (const OList<T>& rhs) const
{
	if (_length != rhs._length)
		return false;
	for (int i = 0; i < _length; ++i)
	{
		if (!(*(_list[i]) == *(rhs._list[i])))
			return false;
	}	
	return true;
}



/* inline void OList<T>::swap(int index1, int index2)
 *
 * Swap two values in list
 */

template <class T>
inline void OList<T>::swap(int index1, int index2)
{
	T* temp = _list[index1];
	_list[index1] = _list[index2];
	_list[index2] = temp;
}



/* inline OList<T>& OList<T>::operator= (const T& rhs)
 *
 * Assignment of object to list
 */

template <class T>
inline OList<T>& OList<T>::operator= (const T& rhs)
{
	
	// Remove all data from current list
	clear();
	
	// Allocate space
	length(1);
	
	// Save value
	_list[0]->operator=(rhs);
	
	// Return result
	return *this;
}



/* inline OList<T>& OList<T>::operator= (const OList& rhs)
 *
 * Assignment operator for OList object
 */

template <class T>
inline OList<T>& OList<T>::operator= (const OList& rhs)
{
	
	// Make sure not copying self
	if (this != &rhs)
	{
		
		// Remove all data from current list
		clear();

		// Allocate space
		length(rhs._length);

		// Copy list
		for (int i = 0; i < rhs._length; ++i)
			*(_list[i]) = *(rhs._list[i]);
	}
	
	// Return result
	return *this;
}



/* OList<T> OList<T>::operator+ (const T& rhs) const
 *
 * Add value to list
 */

template <class T>
OList<T> OList<T>::operator+ (const T& rhs) const
{
	
	// Variable to store new list
	OList<T> res;
	res.length(_length + 1);
	
	// Save current list values
	for (int i = 0; i < _length; ++i)
		*(res._list[i]) = *(_list[i]);
	
	// Save new value
	res._list[_length]->operator=(rhs);
	
	// Return result
	return res;
}



/* OList<T> OList<T>::operator+ (const OList<T>& rhs) const
 *
 * Add two lists
 */

template <class T>
OList<T> OList<T>::operator+ (const OList<T>& rhs) const
{
	
	// Variable to store new list
	OList<T> res;
	res.length(_length + rhs._length);
	
	// Save current list values
	int i;
	for (i = 0; i < _length; ++i)
		*(res._list[i]) = *(_list[i]);

	// Save rhs list values
	for (i = 0; i < rhs._length; ++i)
		*(res._list[i + _length]) = *(rhs._list[i]);
	
	// Return new list
	return res;
}



/* inline OList<T>& OList<T>::operator+= (const T& rhs)
 *
 * Add value to current list
 */

template <class T>
inline OList<T>& OList<T>::operator+= (const T& rhs)
{
	
	// Allocate space
	length(_length + 1);
	
	// Save new value
	_list[_length - 1]->operator=(rhs);
	
	// Return result
	return *this;
}



/* inline OList<T>& OList<T>::operator+= (const OList<T>& rhs)
 *
 * Add list to current
 */

template <class T>
inline OList<T>& OList<T>::operator+= (const OList<T>& rhs)
{
	
	// Allocate space
	length(_length + rhs._length);
	
	// Save new list values
	for (int i = 0; i < rhs._length; ++i)
		*(_list[_length - rhs._length + i]) = *(rhs._list[i]);
	
	// Return result
	return *this;
}



// =====================================================================================================================
// Simple list
// =====================================================================================================================

// Buffer size for List object
template <class T>
int List<T>::_buffer = 4;



/* inline void List<T>::initialize()
 *
 * Initialize List object
 */

template <class T>
inline void List<T>::initialize()
{
	_actualLength = 0;
	_length = 0;
	_list = 0;
}



/* inline void List<T>::clear()
 *
 * Clear data in List object
 */

template <class T>
inline void List<T>::clear()
{
	if (_list)
		delete [] _list;
	initialize();
}



/* void List<T>::remove(int index)
 *
 * Remove value from list at index
 */

template <class T>
void List<T>::remove(int index)
{
	
	// Copy over value to remove
	for (int i = index; i < _length-1; ++i)
		_list[i] = _list[i+1];
	
	// Reset length
	_length--;
}



/* void List<T>::length(int inLength)
 *
 * Set the length of the list
 */

template <class T>
void List<T>::length(int inLength)
{
	
	// Requested length is greater than or equal to actual length
	if (inLength > _actualLength)
	{
		
		// Figure out new buffered size
		int newLength = _buffer * (int) ceil((double) inLength / _buffer);
		
		// Create new list
		T* temp = new T [newLength];
		
		// Copy old values
		for (int i = 0; i < _length; ++i)
			temp[i] = _list[i];
		
		// Delete old list
		if (_list)
			delete [] _list;
		
		// Save list
		_list = temp;
		
		// Save new actual length
		_actualLength = newLength;
	}
	
	// Set length
	_length  = inLength;
}



/* inline void List<T>::fill(const T& input)
 *
 * Fill list with constant value
 */

template <class T>
inline void List<T>::fill(const T& input)
{
	for (int i = 0; i < _length; ++i)
		_list[i] = input;
}



/* inline bool List<T>::operator== (const List<T>& rhs) const
 *
 * Test whether two lists are the same
 */

template <class T>
inline bool List<T>::operator== (const List<T>& rhs) const
{
	if (_length != rhs._length)
		return false;
	for (int i = 0; i < _length; ++i)
	{
		if (_list[i] != rhs._list[i])
			return false;
	}
	return true;
}



/* inline void List<T>::swap(int index1, int index2)
 *
 * Swap two values in list
 */

template <class T>
inline void List<T>::swap(int index1, int index2)
{
	T temp = _list[index1];
	_list[index1] = _list[index2];
	_list[index2] = temp;
}



/* inline void List<T>::set(int num, const T* array)
 *
 * Set array
 */

template <class T>
inline void List<T>::set(int num, const T* array)
{
	
	// Remove all data from current list
	clear();
	
	// Allocate space
	length(num);
	
	// Save data
	for (int i = 0; i < num; ++i)
		_list[i] = array[i];
}



/* inline List<T>& List<T>::operator= (const T& rhs)
 *
 * Assignment of value to list
 */

template <class T>
inline List<T>& List<T>::operator= (const T& rhs)
{
	
	// Remove all data from current list
	clear();
	
	// Allocate space
	length(1);
	
	// Save value
	_list[0] = rhs;
	
	// Return result
	return *this;
}



/* inline List<T>& List<T>::operator= (const List& rhs)
 *
 * Assignment operator for List object
 */

template <class T>
inline List<T>& List<T>::operator= (const List& rhs)
{
	
	// Make sure not copying self
	if (this != &rhs)
	{
		
		// Remove all data from current list
		clear();

		// Allocate space
		length(rhs._length);

		// Copy list
		for (int i = 0; i < rhs._length; ++i)
			_list[i] = rhs._list[i];
	}
	
	// Return result
	return *this;
}



/* List<T> List<T>::operator+ (const T& rhs) const
 *
 * Add value to list
 */

template <class T>
List<T> List<T>::operator+ (const T& rhs) const
{
	
	// Variable to store new list
	List<T> res;
	res.length(_length + 1);
	
	// Save current list values
	for (int i = 0; i < _length; ++i)
		res._list[i] = _list[i];
	
	// Save new value
	res._list[_length] = rhs;
	
	// Return result
	return res;
}



/* List<T> List<T>::operator+ (const List<T>& rhs) const
 *
 * Add two lists
 */

template <class T>
List<T> List<T>::operator+ (const List<T>& rhs) const
{
	
	// Variable to store new list
	List<T> res;
	res.length(_length + rhs._length);
	
	// Save current list values
	int i;
	for (i = 0; i < _length; ++i)
		res._list[i] = _list[i];

	// Save rhs list values
	for (i = 0; i < rhs._length; ++i)
		res.list[i + _length] = rhs._list[i];
	
	// Return new list
	return res;
}



/* inline List<T>& List<T>::operator+= (const T& rhs)
 *
 * Add value to current list
 */

template <class T>
inline List<T>& List<T>::operator+= (const T& rhs)
{
	
	// Allocate space
	length(_length + 1);
	
	// Save new value
	_list[_length - 1] = rhs;
	
	// Return result
	return *this;
}



/* inline List<T>& List<T>::operator+= (const List<T>& rhs)
 *
 * Add list to current
 */

template <class T>
inline List<T>& List<T>::operator+= (const List<T>& rhs)
{
	
	// Allocate space
	length(_length + rhs._length);
	
	// Save new list values
	for (int i = 0; i < rhs._length; ++i)
		_list[_length - rhs._length + i] = rhs._list[i];
	
	// Return result
	return *this;
}



/* inline void List<T>::sort(SortMethod method)
 *
 * Sort values in list
 */

template <class T>
inline void List<T>::sort(SortMethod method)
{
	if (method == QUICKSORT)
		quicksort(0, _length - 1);
}



/* void List<T>::quicksort(int left, int right)
 *
 * Quicksort algorithm
 */

template <class T>
void List<T>::quicksort(int left, int right)
{
	
	// Current partition has no width
	if (left >= right)
		return;
	
	// Choose pivot index
	int pivotIndex = (left + right) / 2;
	T pivot = _list[pivotIndex];
	
	// Move pivot to end
	swap(pivotIndex, right);
	
	// Iterate through current partition
	int newPivotIndex = left;
	for (int i = left; i < right; ++i)
	{
		if (_list[i] < pivot)
		{
			swap(i, newPivotIndex);
			newPivotIndex++;
		}
	}
	
	// Move pivot to final position
	swap(newPivotIndex, right);
	
	// Recursive calls to next sorts
	quicksort(left, newPivotIndex - 1);
	quicksort(newPivotIndex + 1, right);
}



// =====================================================================================================================
// Linked list
// =====================================================================================================================

/* inline Linked<T>::Linked(const Linked<T>& rhs)
 *
 * Copy constructor for linked list
 */

template <class T>
inline Linked<T>::Linked(const Linked<T>& rhs)
{
	initialize();
	for (iterator it = rhs.begin(); it != rhs.end(); ++it)
		add(*it);
}



/* void Linked<T>::clear()
 *
 * Clear data in linked list
 */

template <class T>
void Linked<T>::clear()
{
	Link<T>* next;
	Link<T>* current = _begin;
	while (current)
	{
		next = current->next();
		delete current;
		current = next;
	}
	_begin = 0;
	_last = 0;
	_length = 0;
}



/* Linked<T>& Linked<T>::operator= (const Linked<T>& rhs)
 *
 * Assignment operator for Linked object
 */

template <class T>
Linked<T>& Linked<T>::operator= (const Linked<T>& rhs)
{
	clear();
	for (iterator it = rhs.begin(); it != rhs.end(); ++it)
		add(*it);
}



/* void Linked<T>::add()
 *
 * Add value to linked list
 */

template <class T>
void Linked<T>::add()
{
	
	// Save that length is increased
	++_length;
	
	// List is empty
	if (!_begin)
	{
		_begin = new Link<T>(0, 0);
		_last = _begin;
		return;
	}
	
	// Add value to end of list
	_last->next() = new Link<T> (_last, 0);
	_last = _last->next();
}



/* void Linked<T>::add(const T& value)
 *
 * Add value to linked list
 */

template <class T>
void Linked<T>::add(const T& value)
{
	
	// Save that length is increased
	++_length;
	
	// List is empty
	if (!_begin)
	{
		_begin = new Link<T>(value, 0, 0);
		_last = _begin;
		return;
	}
	
	// Add value to end of list
	_last->next() = new Link<T> (value, _last, 0);
	_last = _last->next();
}



/* void Linked<T>::addBefore(iterator& it, const T& value)
 *
 * Add value within linked list
 */

template <class T>
void Linked<T>::addBefore(iterator& it, const T& value)
{
	
	// Save that length is increased
	++_length;
	
	// Adding to end
	if (it == iterator(0))
	{
		Link<T>* newLink = new Link<T>(value, _last, 0);
		if (_last)
			_last->next() = newLink;
		_last = newLink;
		if (!_begin)
			_begin = _last;
	}
	
	// Adding to beginning
	else if (it == iterator(_begin))
	{
		Link<T>* newLink = new Link<T>(value, 0, _begin);
		_begin->prev() = newLink;
		_begin = newLink;
	}
	
	// Adding in middle
	else
	{
		Link<T>* newLink = new Link<T>(value, it.pointer()->prev(), it.pointer());
		newLink->prev()->next() = newLink;
		newLink->next()->prev() = newLink;
	}
}



/* void Linked<T>::addAfter(iterator& it, const T& value)
 *
 * Add value within linked list
 */

template <class T>
void Linked<T>::addAfter(iterator& it, const T& value)
{
	
	// Save that length is increased
	++_length;
	
	// Adding to beginning
	if (it == iterator(0))
	{
		Link<T>* newLink = new Link<T>(value, 0, _begin);
		if (_begin)
			_begin->prev() = newLink;
		_begin = newLink;
		if (!_last)
			_last = _begin;
	}
	
	// Adding to end
	else if (it == iterator(_last))
	{
		Link<T>* newLink = new Link<T>(value, _last, 0);
		_last->next() = newLink;
		_last = newLink;
	}
	
	// Adding in middle
	else
	{
		Link<T>* newLink = new Link<T>(value, it.pointer(), it.pointer()->next());
		newLink->prev()->next() = newLink;
		newLink->next()->prev() = newLink;
	}
}



/* void Linked<T>::remove(iterator link)
 *
 * Remove link for list
 */

template <class T>
void Linked<T>::remove(iterator link)
{
	
	// Point previous link to next
	if (link.pointer()->prev() != 0)
		link.pointer()->prev()->next() = link.pointer()->next();
	else
		_begin = link.pointer()->next();
	
	// Point next link to previous
	if (link.pointer()->next() != 0)
		link.pointer()->next()->prev() = link.pointer()->prev();
	else
		_last = link.pointer()->prev();
	
	// Delete link
	delete link.pointer();
	
	// Save that list is shorter
	--_length;
}



#endif
