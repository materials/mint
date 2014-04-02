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



#include "text.h"



// Buffer size for Word object
int Word::_bufferSize = 8;

// Stream buffer
char Word::_buffer[100];



/* void Text::addLine(const char* line)
 *
 * Add line of text, breaking by spaces
 */

void Text::addLine(const char* line)
{
	
	// Add new line
	addLine();
	
	// Count the number of words in line
	int start = 0;
	int end = 0;
	int numWords = 0;
	while (getRange(start, end, end, line))
		++numWords;
	
	// Allocate space for words
	_text[_curLine].length(numWords);
	
	// Save words
	int i;
	start = 0;
	end = 0;
	while (getRange(start, end, end, line))
	{
		
		// Allocate space in word
		_text[_curLine][++_curWord].length(end-start);
		
		// Save characters
		for (i = start; i < end; ++i)
			_text[_curLine][_curWord][i-start] = line[i];
	}
}



/* void Text::split(char value)
 *
 * Split text at a given character
 */

void Text::split(char value)
{
	
	// Loop over lines
	int i, j, k, m;
	Word temp;
	for (i = 0; i < _text.length(); ++i)
	{
		
		// Loop over words on line
		for (j =  _text[i].length() - 1; j >= 0; --j)
		{
			
			// Loop over characters
			for (k = _text[i][j].length() - 1; k >= 0; --k)
			{
				
				// Found split character
				if (_text[i][j][k] == value)
				{
					
					// Value is at end of word
					if (k == _text[i][j].length() - 1)
					{
						
						// Value is only character in word
						if (k == 0)
							_text[i].remove(j);
						
						// Characters exist before
						else
						{
							_text[i][j].length(k);
							++j;
						}
					}
					
					// Value is at beginning of word
					else if (k == 0)
					{
						
						// Save word
						temp = _text[i][j];
						_text[i][j].set(1, temp.length() - 1, temp);
					}
					
					// Value occurs in middle of word
					else
					{
						
						// Save word
						temp = _text[i][j];
						
						// Add word to line and move words back
						_text[i].add();
						for (m = _text[i].length() - 1; m > j; --m)
							_text[i].swap(m, m - 1);
						
						// Set words
						_text[i][j].set(0, k, temp);
						_text[i][j+1].set(k+1, temp.length() - k - 1, temp);
						++j;
					}
					
					// Break since found
					break;
				}
			}
		}
	}
}



/* bool Text::equal(const char* lhs, const char* rhs, bool caseMatters, int minToMatch)
 *
 * Return whether two char arrays are the same
 */

bool Text::equal(const char* lhs, const char* rhs, bool caseMatters, int minToMatch)
{
	
	// Must match all characters
	if (!minToMatch)
	{
		
		// Check if characters are the same, including case
		if (caseMatters)
		{
			const char* itLHS;
			const char* itRHS;
			for (itLHS = lhs, itRHS = rhs; 1; itLHS++, itRHS++)
			{
				if (*itLHS != *itRHS)
					return false;
				if (*itLHS == '\0')
					break;
			}
		}
		
		// Check if characters are the same, ignoring case
		else
		{
			const char* itLHS;
			const char* itRHS;
			for (itLHS = lhs, itRHS = rhs; 1; itLHS++, itRHS++)
			{
				if (tolower(*itLHS) != tolower(*itRHS))
					return false;
				if (*itLHS == '\0')
					break;
			}
		}
	}
	
	// Only the first set number of chacters need to match
	else
	{
		
		// Check if characters are the same, including case
		if (caseMatters)
		{
			int count = 0;
			const char* itLHS;
			const char* itRHS;
			for (itLHS = lhs, itRHS = rhs; 1; itLHS++, itRHS++)
			{
				if (*itLHS != *itRHS)
					return false;
				if ((++count == minToMatch) || (*itLHS == '\0') || (*itRHS == '\0'))
					break;
			}
		}
		
		// Check if characters are the same, ignoring case
		else
		{
			int count = 0;
			const char* itLHS;
			const char* itRHS;
			for (itLHS = lhs, itRHS = rhs; 1; itLHS++, itRHS++)
			{
				if (tolower(*itLHS) != tolower(*itRHS))
					return false;
				if ((++count == minToMatch) || (*itLHS == '\0') || (*itRHS == '\0'))
					break;
			}
		}
	}
	
	// Return that words are the same if at this point
	return true;
}



/* bool Text::equal(const char* lhs, const char* rhs, bool caseMatters, int minToMatch, char ignore)
 *
 * Check if two words are equal, ignoring set character
 */

bool Text::equal(const char* lhs, const char* rhs, bool caseMatters, int minToMatch, char ignore)
{
	
	// Must match all characters
	if (!minToMatch)
	{
		
		// Check if characters are the same, including case
		if (caseMatters)
		{
			const char* itLHS;
			const char* itRHS;
			for (itLHS = lhs, itRHS = rhs; 1; ++itLHS, ++itRHS)
			{
				while ((*itLHS == ignore) && (*itLHS != '\0'))
					++itLHS;
				while ((*itRHS == ignore) && (*itRHS != '\0'))
					++itRHS;
				if (*itLHS != *itRHS)
					return false;
				if (*itLHS == '\0')
					break;
			}
		}
		
		// Check if characters are the same, ignoring case
		else
		{
			const char* itLHS;
			const char* itRHS;
			for (itLHS = lhs, itRHS = rhs; 1; ++itLHS, ++itRHS)
			{
				while ((tolower(*itLHS) == ignore) && (*itLHS != '\0'))
					++itLHS;
				while ((tolower(*itRHS) == ignore) && (*itRHS != '\0'))
					++itRHS;
				if (tolower(*itLHS) != tolower(*itRHS))
					return false;
				if (*itLHS == '\0')
					break;
			}
		}
	}
	
	// Only the first set number of chacters need to match
	else
	{
		
		// Check if characters are the same, including case
		if (caseMatters)
		{
			int count = 0;
			const char* itLHS;
			const char* itRHS;
			for (itLHS = lhs, itRHS = rhs; 1; ++itLHS, ++itRHS)
			{
				while ((*itLHS == ignore) && (*itLHS != '\0'))
					++itLHS;
				while ((*itRHS == ignore) && (*itRHS != '\0'))
					++itRHS;
				if (*itLHS != *itRHS)
					return false;
				if ((++count == minToMatch) || (*itLHS == '\0') || (*itRHS == '\0'))
					break;
			}
		}
		
		// Check if characters are the same, ignoring case
		else
		{
			int count = 0;
			const char* itLHS;
			const char* itRHS;
			for (itLHS = lhs, itRHS = rhs; 1; ++itLHS, ++itRHS)
			{
				while ((tolower(*itLHS) == ignore) && (*itLHS != '\0'))
					++itLHS;
				while ((tolower(*itRHS) == ignore) && (*itRHS != '\0'))
					++itRHS;
				if (tolower(*itLHS) != tolower(*itRHS))
					return false;
				if ((++count == minToMatch) || (*itLHS == '\0') || (*itRHS == '\0'))
					break;
			}
		}
	}
	
	// Return that words are the same if at this point
	return true;
}



/* bool Text::contains(const char* text, const char* item, bool caseMatters)
 *
 * Check whether char array contains item
 */

bool Text::contains(const char* text, const char* item, bool caseMatters)
{
	if (item == 0)
		return true;
	if (item[0] == '\0')
		return true;
	int i, j, k;
	for (i = 0; text[i] != '\0'; ++i)
	{
		if (((caseMatters) && (text[i] == item[0])) || \
			((!caseMatters) && (std::tolower(text[i]) == std::tolower(item[0]))))
		{
			for (j = 0, k = i; (item[j] != '\0') && (text[k] != '\0'); ++j, ++k)
			{
				if (caseMatters)
				{
					if (item[j] != text[k])
						break;
				}
				else
				{
					if (std::tolower(item[j]) != std::tolower(text[k]))
						break;
				}
			}
			if (item[j] == '\0')
				return true;
		}
	}
	return false;
}
