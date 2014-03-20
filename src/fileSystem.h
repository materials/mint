/* fileSystem.h -- Deal with files and directories
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
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
