/* Copyright 2011-2014 Kyle Michel, Logan Ward
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



#ifndef CIF_H
#define CIF_H



#include "iso.h"
#include "text.h"
#include "list.h"
#include "fileSystem.h"



// Class for CIF files
class CIF
{
	
	// Tags
	enum Tag {tag_unknown, atom_site_label, atom_site_type_symbol, atom_site_fract_x, atom_site_fract_y, \
		atom_site_fract_z, atom_site_occupancy, cell_length_a, cell_length_b, cell_length_c, cell_angle_alpha, \
		cell_angle_beta, cell_angle_gamma, space_group_IT_number, space_group_symop_operation_xyz, \
		symmetry_Int_Tables_number, symmetry_space_group_name_H_M, symmetry_equiv_pos_as_xyz, chemical_formula_sum, \
		cell_formula_units_Z};
	
	// Functions
	static bool checkTag(Tag tag, OList<Word>& curLine, Vector3D& lengths, Vector3D& angles, \
		Word& spaceGroupName, OList<Element>& elements, List<double>& composition, int& Z);
	static Words getLine(const OList<Word>& line);
	static Tag getTag(const Word& word);
	static void getElementAndNumber(Word& element, Word& number, const Word& word);
	
public:
	
	// Read structure
	static ISO read(const Text& content, double clusterTol = 0.1);
	static ISO read(const Word& file, double clusterTol = 0.1)		{ return read(Read::text(file), clusterTol); }
	
	// Write structure
    static void write(const Word& file, const ISO& iso, double tol);

	// Check if file is in correct format
	static bool isFormat(const Text& content);
	static bool isFormat(const Word& file)		{ return isFormat(Read::text(file)); }
};



#endif
