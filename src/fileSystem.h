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



#ifndef FILESYSTEM_H
#define FILESYSTEM_H



#include "text.h"



// Namespace for functions to read from file
namespace Read
{
	Text text(const Word& file);
}



// Namespace for file functions
namespace File
{
    void copy(const Word& from, const Word& to);
    bool exists(const Word& file, bool exitIfNotFound = false);
    void create(const Word& file, bool clearIfExists = true);
    void remove(const Word& file);
}



// Namespace for directory functions
namespace Directory
{
	char* home();
    bool exists(const Word& path);
    void create(const Word& path, bool deleteIfExists);
    void remove(const Word& path);
	void recurseRemove(const Word& path);
	bool current(Word& path);
	Word makePath(const Word& base, const Word& append);
	void change(const Word& newDir);
}



#endif
