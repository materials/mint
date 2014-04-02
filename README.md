Materials Interface (Mint)
====

Installation
----

**Method 1** (*preferred*)
- Copy code from "update.sh" to your local machine
- Type on your local machine "chmod u+x update.sh"
- Run on your local machine "./update.sh"
- Edit the Makefile with your settings and run "make"
- In the future, just run "./update.sh" anytime that you want to update to a new version

**Method 2**
- Download all files in the repository using the "Download ZIP" link at https://github.com/materials/mint
- Unzip the files on your local machine
- Edit the Makefile with your settings and run "make"
- In the future, just run "./update.sh" anytime that you want to update to a new version

Getting Started
----

Once installed, run "mint -h" for help getting started. In general, a single call to mint will have the form
	mint _input_files_ -function_1 args_1 -function_2 args_2 ...
where _input_files_ is a list of files (structures, settings, etc.) -function_1 and -function_2 are functions to execute and args_1 and args_2 are arguments to each. Any number of functions can be passed in a single call. A list of functions that are available can be obtained by running "mint -h functions". As a practical example, the call to
	mint structureFile -conventional -symmetry -print screen
would read the structure in structureFile, convert it to a conventional cell, print the symmetry of the conventional cell, and print the structure to stdout.

