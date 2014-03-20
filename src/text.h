/* text.h -- Storage for a word and text
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#ifndef TEXT_H
#define TEXT_H



#include "list.h"
#include <iostream>
#include <cmath>
using namespace std;



// Class to store a word
class Word
{
	
	// Buffer size
	static int _bufferSize;
	
	// Variables
	int _actualLength;
	int _length;
	char* _word;
	
	// Helper variables
	static char _buffer[100];
	
	// Functions
	void setBlank();
	void empty();
	void initialize() { length(0); }
	
	// Helper functions
	static int charLength(const char* array);
	
public:
	
	// Constructors
	Word()										{ setBlank(); initialize(); }
	Word(const char copy)						{ setBlank(); *this = copy; }
	Word(const char* copy)						{ setBlank(); *this = copy; }
	Word(const Word& copy)						{ setBlank(); *this = copy; }
	Word(int size, const char* copy)			{ setBlank(); set(size, copy); }
	Word(int size, const Word& copy)			{ setBlank(); set(size, copy); }
	Word(int start, int size, const char* copy)	{ setBlank(); set(start, size, copy); }
	Word(int start, int size, const Word& copy)	{ setBlank(); set(start, size, copy); }
	
	// Destructor
	~Word() { empty(); }
	
	// Delete all data
	void clear() { empty(); initialize(); }
	
	// Comparison functions
	bool equal(const Word& rhs, bool caseMatters, int minToMatch) const;
	bool equal(const Word& rhs, bool caseMatters) const					{ return equal(rhs, caseMatters, 0); }
	bool equal(const Word& rhs, int minToMatch) const					{ return equal(rhs, true, minToMatch); }
	bool equal(const Word& rhs) const									{ return equal(rhs, true, 0); }
	bool operator== (const Word& rhs) const								{ return equal(rhs, true, 0); }
	bool operator!= (const Word& rhs) const								{ return (!equal(rhs, true, 0)); }
	bool equal(const char* rhs, bool caseMatters, int minToMatch) const;
	bool equal(const char* rhs, bool caseMatters) const					{ return equal(rhs, caseMatters, 0); }
	bool equal(const char* rhs, int minToMatch) const					{ return equal(rhs, true, minToMatch); }
	bool equal(const char* rhs) const									{ return equal(rhs, true, 0); }
	bool operator== (const char* rhs) const								{ return equal(rhs, true, 0); }
	bool equal(const Word& rhs, bool caseMatters, int minToMatch, char ignore) const;
	bool equal(const Word& rhs, bool caseMatters, char ignore) const	{ return equal(rhs, caseMatters, 0, ignore); }
	bool equal(const Word& rhs, int minToMatch, char ignore) const		{ return equal(rhs, true, minToMatch, ignore); }
	bool equal(const Word& rhs, char ignore) const						{ return equal(rhs, true, 0, ignore); }
	bool equal(const char* rhs, bool caseMatters, int minToMatch, char ignore) const;
	bool equal(const char* rhs, bool caseMatters, char ignore) const	{ return equal(rhs, caseMatters, 0, ignore); }
	bool equal(const char* rhs, int minToMatch, char ignore) const		{ return equal(rhs, true, minToMatch, ignore); }
	bool equal(const char* rhs, char ignore) const						{ return equal(rhs, true, 0, ignore); }
	bool contains(const Word& rhs, bool caseMatters = true) const;
	bool contains(const char* rhs, bool caseMatters = true) const;
	
	// Assignment functions
	void length(int size);
	void set(int size, const char* rhs);
	void set(int size, const Word& rhs);
	void set(int start, int size, const char* rhs);
	void set(int start, int size, const Word& rhs);
	Word& operator=  (const char  rhs);
	Word& operator=  (const char* rhs);
	Word& operator=  (const Word& rhs);
	Word  operator+  (const char  rhs) const;
	Word  operator+  (const char* rhs) const;
	Word  operator+  (const Word& rhs) const;
	Word& operator+= (const char  rhs);
	Word& operator+= (const char* rhs);
	Word& operator+= (const Word& rhs);
	void strip(char value);
	void cutAt(int indexOfNewLast);
	Word tolower() const;
	Word toupper() const;
	Word& insert(const char* rhs, int startIndex);
	Word& insert(const Word& rhs, int startIndex);
	
	// Insertion and extraction
	friend ostream& operator<< (ostream& output, const Word& word);
	friend istream& operator>> (istream& input, Word& word);
	
	// Access functions
	int length() const					{ return _length; }
	char& operator[] (int index) const	{ return _word[index]; }
	char& last() const					{ return _word[_length - 1]; }
	char* array() const					{ return _word; }
};



// List of words
class Words
{
	
	// Variables
	OList<Word> _words;
	
public:
	
	// Constructors
	Words()							{}
	Words(const Word& copy)			{ _words = copy; }
	Words(const Words& copy)		{ _words = copy._words; }
	Words(const OList<Word>& copy)	{ _words = copy; }
	Words(int inLength)				{ _words.length(inLength); }
	
	// Setup functions
	void clear()				{ _words.clear(); }
	void add()					{ _words.add(); }
	void add(int num)			{ _words.add(num); }
	void remove(int index)		{ _words.remove(index); }
	void length(int inLength)	{ _words.length(inLength); }
	
	// Assignment functions
	Words& operator=  (const Word& rhs)			{ _words = rhs; return *this; }
	Words& operator=  (const Words& rhs)		{ if (this != &rhs) _words = rhs._words; return *this; }
	Words  operator+  (const Word& rhs) const	{ Words res; res._words = _words + rhs; return res; }
	Words  operator+  (const Words& rhs) const	{ Words res; res._words = _words + rhs._words; return res; }
	Words& operator+= (const Word& rhs)			{ _words += rhs; return *this; }
	Words& operator+= (const Words& rhs)		{ _words += rhs._words; return *this; }
	
	// Inseration
	friend ostream& operator<< (ostream& output, const Word& word);
	
	// Access functions
	int length() const							{ return _words.length(); }
	const Word& operator[] (int index) const	{ return _words[index]; }
	const Word& last() const					{ return _words.last(); }
	Word& operator[] (int index)				{ return _words[index]; }
	Word& last()								{ return _words.last(); }
};



// Class to store a block of text
class Text
{
	
	// Variables
	int _curLine;
	int _curWord;
	OList<Word>::D2 _text;
	
	// Functions
	void initialize() { _curLine = -1; _curWord = -1; }
	bool getRange(int& start, int& end, int searchStart, const char* line);
	
public:
	
	// Constructor
	Text() { initialize(); }
	
	// Clear data
	void clear() { _text.clear(); initialize(); }
	
	// Setup functions
	void addLine();
	void addLine(const char* line);
	void addWord(const Word& word);
	void addWord(const char* word);
	void addLines(int numToAdd)				{ _text.length(_text.length() + numToAdd); }
	void addWords(int numToAdd)				{ _text[_curLine].length(_text[_curLine].length() + numToAdd); }
	void addWords(int line, int numToAdd)	{ _text[line].length(_text[line].length() + numToAdd); }
	void addWords(const Words& words)
		{ addWords(words.length()); for (int i = 0; i < words.length(); ++i) addWord(words[i]); }
	void addWords(const OList<Word>& words)
		{ addWords(words.length()); for (int i = 0; i < words.length(); ++i) addWord(words[i]); }
	
	// Modifier functions
	void split(char value);
	
	// Access functions
	int length() const							{ return _text.length(); }
	OList<Word>& operator[] (int line) const	{ return _text[line]; }
	
	// Static member functions
	static bool equal(const char* lhs, const char* rhs, bool caseMatters, int minToMatch);
	static bool equal(const char* lhs, const char* rhs, bool caseMatters) { return equal(lhs, rhs, caseMatters, 0); }
	static bool equal(const char* lhs, const char* rhs, int minToMatch)	  { return equal(lhs, rhs, true, minToMatch); }
	static bool equal(const char* lhs, const char* rhs)					  { return equal(lhs, rhs, true, 0); }
	static bool equal(const char* lhs, const char* rhs, bool caseMatters, int minToMatch, char ignore);
	static bool equal(const char* lhs, const char* rhs, bool caseMatters, char ignore)
		{ return equal(lhs, rhs, caseMatters, 0); }
	static bool equal(const char* lhs, const char* rhs, int minToMatch, char ignore)
		{ return equal(lhs, rhs, true, minToMatch); }
	static bool equal(const char* lhs, const char* rhs, char ignore)
		{ return equal(lhs, rhs, true, 0); }
	static bool contains(const char* text, const char* item, bool caseMatters = true);
	static int length(const char* rhs);
};



// =====================================================================================================================
// Word
// =====================================================================================================================

/* inline void Word::setBlank()
 *
 * Set word to null
 */

inline void Word::setBlank()
{
	_word = 0;
	_length = 0;
	_actualLength = 0;
}



/* inline void Word::empty()
 *
 * Clear data in Word object
 */

inline void Word::empty()
{
	if (_word)
		delete [] _word;
	setBlank();
}



/* inline void Word::length(int size)
 *
 * Add space to word
 */

inline void Word::length(int size)
{
	
	// Requested length is greater than or equal to actual length
	if (size > _actualLength - 1)
	{
		
		// Figure out new buffered size
		int bufLength = _bufferSize * (int) ceil((double) size / _bufferSize) + 1;
		
		// Create new list
		char* temp = new char [bufLength];
		temp[size] = '\0';
		
		// Copy old values
		for (int i = 0; i < _length; ++i)
			temp[i] = _word[i];
		
		// Delete old list
		if (_word)
			delete [] _word;
		
		// Save list
		_word = temp;
		
		// Save new actual length
		_actualLength = bufLength - 1;
	}
	
	// Set null character
	_length = size;
	_word[size] = '\0';
}



/* inline bool Word::equal(const Word& rhs, bool caseMatters, int minToMatch) const
 *
 * Return whether two words are equal
 */

inline bool Word::equal(const Word& rhs, bool caseMatters, int minToMatch) const
{
	
	// Lengths are not the same
	if (_length != rhs._length)
	{
		
		// Must match all characters
		if (!minToMatch)
			return false;
		
		// Lengths are not the same and at least one is below the number to match
		else if ((_length < minToMatch) || (rhs._length < minToMatch))
			return false;
	}
	
	// Return whether words are the same
	return Text::equal(_word, rhs._word, caseMatters, minToMatch);
}



/* inline bool Word::equal(const char* rhs, bool caseMatters, int minToMatch) const
 *
 * Return whether two words are equal
 */

inline bool Word::equal(const char* rhs, bool caseMatters, int minToMatch) const
{
	return Text::equal(_word, rhs, caseMatters, minToMatch);
}



/* inline bool Word::equal(const Word& rhs, bool caseMatters, int minToMatch, char ignore) const
 *
 * Return whether two words are the same, ignoring set character
 */

inline bool Word::equal(const Word& rhs, bool caseMatters, int minToMatch, char ignore) const
{
	return Text::equal(_word, rhs._word, caseMatters, minToMatch, ignore);
}



/* inline bool Word::equal(const char* rhs, bool caseMatters, int minToMatch, char ignore) const
 *
 * Return whether two words are the same, ignoring set character
 */

inline bool Word::equal(const char* rhs, bool caseMatters, int minToMatch, char ignore) const
{
	return Text::equal(_word, rhs, caseMatters, minToMatch, ignore);
}



/* inline bool Word::contains(const Word& rhs, bool caseMatters) const
 *
 * Return whether a word contains a set of characters
 */

inline bool Word::contains(const Word& rhs, bool caseMatters) const
{
	if (rhs.length() == 0)
		return true;
	int i, j, k;
	for (i = 0; i < _length; ++i)
	{
		if (((caseMatters) && (_word[i] == rhs[0])) || \
			((!caseMatters) && (std::tolower(_word[i]) == std::tolower(rhs[0]))))
		{
			for (j = 0, k = i; (j < rhs._length) && (k < _length); ++j, ++k)
			{
				if (caseMatters)
				{
					if (rhs._word[j] != _word[k])
						break;
				}
				else
				{
					if (std::tolower(rhs._word[j]) != std::tolower(_word[k]))
						break;
				}
			}
			if (j == rhs._length)
				return true;
		}
	}
	return false;
}



/* inline bool Word::contains(const char* rhs, bool caseMatters) const
 *
 * Return whether a word contains a set of characters
 */

inline bool Word::contains(const char* rhs, bool caseMatters) const
{
	if (rhs == 0)
		return true;
	if (rhs[0] == '\0')
		return true;
	int i, j, k;
	for (i = 0; i < _length; ++i)
	{
		if (((caseMatters) && (_word[i] == rhs[0])) || \
			((!caseMatters) && (std::tolower(_word[i]) == std::tolower(rhs[0]))))
		{
			for (j = 0, k = i; (rhs[j] != '\0') && (k < _length); ++j, ++k)
			{
				if (caseMatters)
				{
					if (rhs[j] != _word[k])
						break;
				}
				else
				{
					if (std::tolower(rhs[j]) != std::tolower(_word[k]))
						break;
				}
			}
			if (rhs[j] == '\0')
				return true;
		}
	}
	return false;
}



/* inline void Word::set(int size, const char* rhs)
 *
 * Set word
 */

inline void Word::set(int size, const char* rhs)
{
	clear();
	length(size);
	for (int i = 0; i < size; ++i)
		_word[i] = rhs[i];
}



/* inline void Word::set(int size, const Word& rhs)
 *
 * Set word
 */

inline void Word::set(int size, const Word& rhs)
{
	if (this != &rhs)
	{
		clear();
		length(size);
		for (int i = 0; i < size; ++i)
			_word[i] = rhs[i];
	}
}



/* inline void Word::set(int start, int size, const char* rhs)
 *
 * Set word
 */

inline void Word::set(int start, int size, const char* rhs)
{
	clear();
	length(size);
	for (int i = start; i < start+size; ++i)
		_word[i-start] = rhs[i];
}



/* inline void Word::set(int start, int size, const Word& rhs)
 *
 * Set word
 */

inline void Word::set(int start, int size, const Word& rhs)
{
	if (this != &rhs)
	{
		clear();
		length(size);
		for (int i = start; i < start+size; ++i)
			_word[i-start] = rhs[i];
	}
}



/* inline Word& Word::operator= (const char rhs)
 *
 * Assign character to word
 */

inline Word& Word::operator= (const char rhs)
{
	
	// Clear space
	empty();
	
	// Allocate space
	length(1);
	
	// Save character
	_word[0] = rhs;
	
	// Return result
	return *this;
}



/* inline Word& Word::operator= (const char* rhs)
 *
 * Assign char array to word
 */

inline Word& Word::operator= (const char* rhs)
{
	
	// Clear space
	empty();
	
	// Allocate space
	length(charLength(rhs));

	// Save characters
	for (int i = 0; i < _length; ++i)
		_word[i] = rhs[i];
	
	// Return result
	return *this;
}



/* inline Word& Word::operator= (const Word& rhs)
 *
 * Copy Word object
 */

inline Word& Word::operator= (const Word& rhs)
{
	
	// Make sure not copying self
	if (this != &rhs)
	{
		
		// Clear space
		empty();
		
		// Allocate space
		length(rhs._length);

		// Save characters
		for (int i = 0; i < rhs._length; ++i)
			_word[i] = rhs._word[i];
	}
	
	// Return result
	return *this;
}



/* inline Word Word::operator+ (const char rhs) const
 *
 * Add character to word
 */

inline Word Word::operator+ (const char rhs) const
{
	
	// Variable to store result
	Word res;
	res.length(_length + 1);
	
	// Save current word
	for (int i = 0; i < _length; ++i)
		res[i] = _word[i];
	
	// Save character
	res[_length] = rhs;
	
	// Return result
	return res;
}



/* inline Word Word::operator+ (const char* rhs) const
 *
 * Add char array to word
 */

inline Word Word::operator+ (const char* rhs) const
{
	
	// Variable to store result
	Word res;
	int rhsLength = charLength(rhs);
	res.length(_length + rhsLength);
	
	// Save current word
	int i;
	for (i = 0; i < _length; ++i)
		res[i] = _word[i];
	
	// Save character array
	for (i = 0; i < rhsLength; ++i)
		res[i + _length] = rhs[i];
	
	// Return result
	return res;
}



/* inline Word Word::operator+ (const Word& rhs) const
 *
 * Add word to word
 */

inline Word Word::operator+ (const Word& rhs) const
{
	
	// Variable to store result
	Word res;
	res.length(_length + rhs._length);
	
	// Save current word
	int i;
	for (i = 0; i < _length; ++i)
		res[i] = _word[i];
	
	// Save character array
	for (i = 0; i < rhs._length; ++i)
		res[i + _length] = rhs._word[i];
	
	// Return result
	return res;
}



/* inline Word& Word::operator+= (const char rhs)
 *
 * Add character to word
 */

inline Word& Word::operator+= (const char rhs)
{
	
	// Allocate space
	length(_length + 1);
	
	// Save new character
	_word[_length - 1] = rhs;
	
	// Return result
	return *this;
}



/* inline Word& Word::operator+= (const char* rhs)
 *
 * Add char array to word
 */

inline Word& Word::operator+= (const char* rhs)
{
	
	// Save current length
	int origLen = _length;
	
	// Allocate space
	int rhsLength = charLength(rhs);
	length(_length + rhsLength);
	
	// Save new char array
	for (int i = 0; i < rhsLength; ++i)
		_word[i + origLen] = rhs[i];
	
	// Return result
	return *this;
}



/* inline Word& Word::operator+= (const Word& rhs)
 *
 * Add word to word
 */

inline Word& Word::operator+= (const Word& rhs)
{
	
	// Save current length
	int origLen = _length;
	
	// Allocate space
	length(_length + rhs._length);
	
	// Save new char array
	for (int i = 0; i < rhs._length; ++i)
		_word[i + origLen] = rhs._word[i];
	
	// Return result
	return *this;
}



/* inline void Word::strip(char value)
 *
 * Remove character from word
 */

inline void Word::strip(char value)
{
	int i, j;
	for (i = _length - 1; i >= 0; --i)
	{
		if (_word[i] == value)
		{
			for (j = i; j < _length; ++j)
				_word[j] = _word[j+1];
			length(_length - 1);
		}
	}
}



/* inline void Word::cutAt(int indexOfNewLast)
 *
 * Set a new end to word
 */

inline void Word::cutAt(int indexOfNewLast)
{
	_length = indexOfNewLast + 1;
	_word[_length] = '\0';
}



/* inline Word Word::tolower() const
 *
 * Convert all characters to lowercase
 */

inline Word Word::tolower() const
{
	Word res;
	res.length(_length);
	for (int i = 0; i < _length; ++i)
		res[i] = std::tolower(_word[i]);
	return res;
}



/* inline Word Word::toupper() const
 *
 * Convert all characters to uppercase
 */

inline Word Word::toupper() const
{
	Word res;
	res.length(_length);
	for (int i = 0; i < _length; ++i)
		res[i] = std::toupper(_word[i]);
	return res;
}



/* inline Word& Word::insert(const char* rhs, int startIndex)
 *
 * Insert characters into word
 */

inline Word& Word::insert(const char* rhs, int startIndex)
{
	
	// Allocate space
	int rhsLength = charLength(rhs);
	length(_length + rhsLength);
	
	// Move values ahead
	int i;
	for (i = _length; i >= startIndex + rhsLength; --i)
		_word[i] = _word[i - rhsLength];
	
	// Save new char array
	for (i = 0; i < rhsLength; ++i)
		_word[i + startIndex] = rhs[i];
	
	// Return result
	return *this;
}



/* inline Word& Word::insert(const Word& rhs, int startIndex)
 *
 * Insert characters into word
 */

inline Word& Word::insert(const Word& rhs, int startIndex)
{
	
	// Allocate space
	length(_length + rhs._length);
	
	// Move values ahead
	int i;
	for (i = _length; i >= startIndex + rhs._length; --i)
		_word[i] = _word[i - rhs._length];
	
	// Save new char array
	for (i = 0; i < rhs._length; ++i)
		_word[i + startIndex] = rhs[i];
	
	// Return result
	return *this;
}



/* inline int Word::charLength(const char* array)
 *
 * Return length of char array
 */

inline int Word::charLength(const char* array)
{
	int count = 0;
	while (array[count] != '\0')
		count++;
	return count;
}



/* inline ostream& operator<< (ostream& output, const Word& word)
 *
 * Inseration operator
 */

inline ostream& operator<< (ostream& output, const Word& word)
{
	output << word._word;
	return output;
}



/* inline istream& operator>> (istream& input, Word& word)
 *
 * Extraction operator
 */

inline istream& operator>> (istream& input, Word& word)
{
	input.get(Word::_buffer, 100);
	word = Word::_buffer;
	return input;
}



// =====================================================================================================================
// Words
// =====================================================================================================================

/* inline ostream& operator<< (ostream& output, const Words& words)
 *
 * Insertion operator
 */

inline ostream& operator<< (ostream& output, const Words& words)
{
	for (int i = 0; i < words.length(); ++i)
	{
		output << words[i];
		if (i != words.length() - 1)
			output << " ";
	}
	return output;
}



// =====================================================================================================================
// Text
// =====================================================================================================================

/* inline void Text::addLine()
 *
 * Add line to Text oject
 */

inline void Text::addLine()
{
	_text.length(++_curLine + 1);
	_curWord = -1;
}



/* inline void Text::addWord(const Word& word)
 *
 * Add word to Text object
 */

inline void Text::addWord(const Word& word)
{
	_text[_curLine].length(++_curWord + 1);
	_text[_curLine][_curWord] = word;
}


/* inline void Text::addWord(const char* word)
 *
 * Add word to text object
 */

inline void Text::addWord(const char* word)
{
	_text[_curLine].length(++_curWord + 1);
	_text[_curLine][_curWord] = word;
}



/* inline int Text::length(const char* rhs)
 *
 * Return the length of a string of text
 */

inline int Text::length(const char* rhs)
{
	int res = 0;
	while (rhs[res] != '\0')
		++res;
	return res;
}



/* inline bool Text::getRange(int& start, int& end, int searchStart, const char* line)
 *
 * Get the range of the next word in a line, return false if no word is found
 */

inline bool Text::getRange(int& start, int& end, int searchStart, const char* line)
{
	
	// Loop over contents and skip white space
	start = searchStart;
	while (((line[start] == ' ') || (line[start] == '\t')) && (line[start] != '\0') && (line[start] != '\r'))
		++start;
	
	// No word was found
	if ((line[start] == '\0') || (line[start] == '\r'))
		return false;
	
	// Loop through word
	end = start;
	while ((line[end] != ' ') && (line[end] != '\t') && (line[end] != '\0') && (line[end] != '\r'))
		++end;
	
	// Return that a word was found
	return true;
}



#endif
