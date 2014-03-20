/* bonds.cpp -- Get data about all bonds in a structure
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#include "bonds.h"



/* OList<Bond> Bonds::find(const ISO& iso)
 *
 * Get list of bonds in the structure
 */

OList<Bond> Bonds::find(const ISO& iso)
{
	
	// Loop over pairs of elements in the structure
	int i, j, k, m;
	bool found;
	double min;
	double max;
	double dis;
	OList<Bond> res;
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = i; j < iso.atoms().length(); ++j)
		{
			
			// Look for elements to see if a bond exists
			min = -1;
			
			// Al O
		    if (((iso.atoms()[i][0].element().symbol() == "Al") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "Al")))
		    {
				min = 1.4;
				max = 1.8;
		    }

		    // B Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "B")) || \
		        ((iso.atoms()[i][0].element().symbol() == "B") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")))
		    {
				min = 1.5;
				max = 1.9;
		    }

		    // Be H
		    else if (((iso.atoms()[i][0].element().symbol() == "Be") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "Be")))
		    {
				min = 1.1;
				max = 1.7;
		    }

		    // B H
		    else if (((iso.atoms()[i][0].element().symbol() == "B") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "B")))
		    {
				min = 1.0;
				max = 1.5;
		    }

		    // Br O
		    else if (((iso.atoms()[i][0].element().symbol() == "Br") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "Br")))
		    {
				min = 1.6;
				max = 2.0;
		    }

		    // C C
		    else if ((iso.atoms()[i][0].element().symbol() == "C") && \
				(iso.atoms()[j][0].element().symbol() == "C"))
		    {
				min = 1.0;
				max = 1.8;
		    }

		    // C H
		    else if (((iso.atoms()[i][0].element().symbol() == "C") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "C")))
		    {
				min = 0.8;
				max = 1.4;
		    }

		    // C O
		    else if (((iso.atoms()[i][0].element().symbol() == "C") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "C")))
		    {
				min = 0.9;
				max = 1.7;
		    }

		    // Cu Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Cu") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "Cu")))
		    {
				min = 1.8;
				max = 2.2;
		    }

		    // F Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")))
		    {
				min = 1.2;
				max = 1.9;
		    }

		    // H As
		    else if (((iso.atoms()[i][0].element().symbol() == "As") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "As")))
		    {
				min = 1.3;
				max = 1.7;
		    }

		    // Ca H
		    else if (((iso.atoms()[i][0].element().symbol() == "Ca") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "Ca")))
		    {
				min = 1.7;
				max = 2.3;
		    }

		    // H F
		    else if (((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "F")))
		    {
				min = 0.7;
				max = 1.3;
		    }

		    // H K
		    else if (((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "K")) || \
		        ((iso.atoms()[i][0].element().symbol() == "K") && \
				(iso.atoms()[j][0].element().symbol() == "H")))
		    {
				min = 1.9;
				max = 2.5;
		    }

		    // H Zn
		    else if (((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "Zn")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Zn") && \
				(iso.atoms()[j][0].element().symbol() == "H")))
		    {
				min = 1.3;
				max = 1.8;
		    }

		    // K F
		    else if (((iso.atoms()[i][0].element().symbol() == "K") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "K")))
		    {
				min = 1.8;
				max = 2.4;
		    }

		    // Li H
		    else if (((iso.atoms()[i][0].element().symbol() == "Li") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "Li")))
		    {
				min = 1.3;
				max = 1.8;
		    }

		    // Li O
		    else if (((iso.atoms()[i][0].element().symbol() == "Li") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "Li")))
		    {
				min = 1.3;
				max = 1.9;
		    }

		    // Na Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Na") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "Na")))
		    {
				min = 2.1;
				max = 2.9;
		    }

		    // Al N
		    else if (((iso.atoms()[i][0].element().symbol() == "Al") && \
				(iso.atoms()[j][0].element().symbol() == "N")) || \
		        ((iso.atoms()[i][0].element().symbol() == "N") && \
				(iso.atoms()[j][0].element().symbol() == "Al")))
		    {
				min = 1.4;
				max = 2.0;
		    }

		    // N O
		    else if (((iso.atoms()[i][0].element().symbol() == "N") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "N")))
		    {
				min = 0.8;
				max = 1.7;
		    }

		    // O O
		    else if ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "O"))
		    {
				min = 0.9;
				max = 1.7;
		    }

		    // P Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "P")) || \
		        ((iso.atoms()[i][0].element().symbol() == "P") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")))
		    {
				min = 1.7;
				max = 2.4;
		    }

		    // P O
		    else if (((iso.atoms()[i][0].element().symbol() == "P") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "P")))
		    {
				min = 1.2;
				max = 1.8;
		    }

		    // S H
		    else if (((iso.atoms()[i][0].element().symbol() == "S") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "S")))
		    {
				min = 1.1;
				max = 1.6;
		    }

		    // Si H
		    else if (((iso.atoms()[i][0].element().symbol() == "Si") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "Si")))
		    {
				min = 1.2;
				max = 1.8;
		    }

		    // Si Si
		    else if ((iso.atoms()[i][0].element().symbol() == "Si") && \
				(iso.atoms()[j][0].element().symbol() == "Si"))
		    {
				min = 2.0;
				max = 2.5;
		    }

		    // Zn S
		    else if (((iso.atoms()[i][0].element().symbol() == "Zn") && \
				(iso.atoms()[j][0].element().symbol() == "S")) || \
		        ((iso.atoms()[i][0].element().symbol() == "S") && \
				(iso.atoms()[j][0].element().symbol() == "Zn")))
		    {
				min = 1.8;
				max = 2.3;
		    }

		    // Cl Cl
		    else if ((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "Cl"))
		    {
				min = 1.7;
				max = 2.2;
		    }

		    // Al Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Al") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "Al")))
		    {
				min = 1.8;
				max = 2.4;
		    }

		    // Al S
		    else if (((iso.atoms()[i][0].element().symbol() == "Al") && \
				(iso.atoms()[j][0].element().symbol() == "S")) || \
		        ((iso.atoms()[i][0].element().symbol() == "S") && \
				(iso.atoms()[j][0].element().symbol() == "Al")))
		    {
				min = 1.8;
				max = 2.3;
		    }

		    // Be O
		    else if (((iso.atoms()[i][0].element().symbol() == "Be") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "Be")))
		    {
				min = 1.1;
				max = 1.6;
		    }

		    // B N
		    else if (((iso.atoms()[i][0].element().symbol() == "B") && \
				(iso.atoms()[j][0].element().symbol() == "N")) || \
		        ((iso.atoms()[i][0].element().symbol() == "N") && \
				(iso.atoms()[j][0].element().symbol() == "B")))
		    {
				min = 1.1;
				max = 1.9;
		    }

		    // B S
		    else if (((iso.atoms()[i][0].element().symbol() == "B") && \
				(iso.atoms()[j][0].element().symbol() == "S")) || \
		        ((iso.atoms()[i][0].element().symbol() == "S") && \
				(iso.atoms()[j][0].element().symbol() == "B")))
		    {
				min = 1.4;
				max = 1.9;
		    }

		    // C Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "C")) || \
		        ((iso.atoms()[i][0].element().symbol() == "C") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")))
		    {
				min = 1.4;
				max = 2.0;
		    }

		    // C I
		    else if (((iso.atoms()[i][0].element().symbol() == "C") && \
				(iso.atoms()[j][0].element().symbol() == "I")) || \
		        ((iso.atoms()[i][0].element().symbol() == "I") && \
				(iso.atoms()[j][0].element().symbol() == "C")))
		    {
				min = 1.9;
				max = 2.4;
		    }

		    // Ge Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Ge") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "Ge")))
		    {
				min = 1.9;
				max = 2.4;
		    }

		    // H P
		    else if (((iso.atoms()[i][0].element().symbol() == "P") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "P")))
		    {
				min = 1.1;
				max = 2.1;
		    }

		    // Cu Cu
		    else if ((iso.atoms()[i][0].element().symbol() == "Cu") && \
				(iso.atoms()[j][0].element().symbol() == "Cu"))
		    {
				min = 2.0;
				max = 2.4;
		    }

		    // As F
		    else if (((iso.atoms()[i][0].element().symbol() == "As") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "As")))
		    {
				min = 1.5;
				max = 2.0;
		    }

		    // F F
		    else if ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "F"))
		    {
				min = 1.2;
				max = 1.6;
		    }

		    // Si F
		    else if (((iso.atoms()[i][0].element().symbol() == "Si") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "Si")))
		    {
				min = 1.4;
				max = 1.8;
		    }

		    // H Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")))
		    {
				min = 1.1;
				max = 1.5;
		    }

		    // H Ge
		    else if (((iso.atoms()[i][0].element().symbol() == "Ge") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "Ge")))
		    {
				min = 1.3;
				max = 1.7;
		    }

		    // H Mg
		    else if (((iso.atoms()[i][0].element().symbol() == "Mg") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "Mg")))
		    {
				min = 1.6;
				max = 2.1;
		    }

		    // I I
		    else if ((iso.atoms()[i][0].element().symbol() == "I") && \
				(iso.atoms()[j][0].element().symbol() == "I"))
		    {
				min = 2.4;
				max = 2.8;
		    }

		    // Li Br
		    else if (((iso.atoms()[i][0].element().symbol() == "Li") && \
				(iso.atoms()[j][0].element().symbol() == "Br")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Br") && \
				(iso.atoms()[j][0].element().symbol() == "Li")))
		    {
				min = 2.0;
				max = 2.4;
		    }

		    // Mg Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Mg") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "Mg")))
		    {
				min = 2.0;
				max = 2.4;
		    }

		    // Mg O
		    else if (((iso.atoms()[i][0].element().symbol() == "Mg") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "Mg")))
		    {
				min = 1.5;
				max = 1.9;
		    }

		    // Na F
		    else if (((iso.atoms()[i][0].element().symbol() == "Na") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "Na")))
		    {
				min = 1.7;
				max = 2.1;
		    }

		    // N F
		    else if (((iso.atoms()[i][0].element().symbol() == "N") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "N")))
		    {
				min = 1.1;
				max = 1.7;
		    }

		    // N S
		    else if (((iso.atoms()[i][0].element().symbol() == "N") && \
				(iso.atoms()[j][0].element().symbol() == "S")) || \
		        ((iso.atoms()[i][0].element().symbol() == "S") && \
				(iso.atoms()[j][0].element().symbol() == "N")))
		    {
				min = 1.2;
				max = 1.7;
		    }

		    // P P
		    else if ((iso.atoms()[i][0].element().symbol() == "P") && \
				(iso.atoms()[j][0].element().symbol() == "P"))
		    {
				min = 1.7;
				max = 2.3;
		    }

		    // Se O
		    else if (((iso.atoms()[i][0].element().symbol() == "Se") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "Se")))
		    {
				min = 1.4;
				max = 1.8;
		    }

		    // Si C
		    else if (((iso.atoms()[i][0].element().symbol() == "Si") && \
				(iso.atoms()[j][0].element().symbol() == "C")) || \
		        ((iso.atoms()[i][0].element().symbol() == "C") && \
				(iso.atoms()[j][0].element().symbol() == "Si")))
		    {
				min = 1.6;
				max = 2.0;
		    }

		    // Si N
		    else if (((iso.atoms()[i][0].element().symbol() == "Si") && \
				(iso.atoms()[j][0].element().symbol() == "N")) || \
		        ((iso.atoms()[i][0].element().symbol() == "N") && \
				(iso.atoms()[j][0].element().symbol() == "Si")))
		    {
				min = 1.4;
				max = 1.8;
		    }

		    // Ti Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Ti") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "Ti")))
		    {
				min = 1.9;
				max = 2.3;
		    }

		    // Al F
		    else if (((iso.atoms()[i][0].element().symbol() == "Al") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "Al")))
		    {
				min = 1.4;
				max = 1.8;
		    }

		    // As As
		    else if ((iso.atoms()[i][0].element().symbol() == "As") && \
				(iso.atoms()[j][0].element().symbol() == "As"))
		    {
				min = 1.9;
				max = 2.3;
		    }

		    // Be Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Be") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "Be")))
		    {
				min = 1.6;
				max = 2.0;
		    }

		    // B O
		    else if (((iso.atoms()[i][0].element().symbol() == "B") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "B")))
		    {
				min = 1.0;
				max = 1.4;
		    }

		    // C F
		    else if (((iso.atoms()[i][0].element().symbol() == "C") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "C")))
		    {
				min = 1.1;
				max = 1.5;
		    }

		    // Br Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Br") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "Br")))
		    {
				min = 1.9;
				max = 2.3;
		    }

		    // Cl O
		    else if (((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")))
		    {
				min = 1.2;
				max = 1.9;
		    }

		    // C S
		    else if (((iso.atoms()[i][0].element().symbol() == "C") && \
				(iso.atoms()[j][0].element().symbol() == "S")) || \
		        ((iso.atoms()[i][0].element().symbol() == "S") && \
				(iso.atoms()[j][0].element().symbol() == "C")))
		    {
				min = 1.5;
				max = 2.0;
		    }

		    // Cu F
		    else if (((iso.atoms()[i][0].element().symbol() == "Cu") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "Cu")))
		    {
				min = 1.5;
				max = 1.9;
		    }

		    // F Br
		    else if (((iso.atoms()[i][0].element().symbol() == "Br") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "Br")))
		    {
				min = 1.5;
				max = 1.9;
		    }

		    // F O
		    else if (((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "F")))
		    {
				min = 1.2;
				max = 1.8;
		    }

		    // H Br
		    else if (((iso.atoms()[i][0].element().symbol() == "Br") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "Br")))
		    {
				min = 1.2;
				max = 1.6;
		    }

		    // Cu H
		    else if (((iso.atoms()[i][0].element().symbol() == "Cu") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "Cu")))
		    {
				min = 1.2;
				max = 1.6;
		    }

		    // H H 
		    else if ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "H"))
		    {
				min = 0.6;
				max = 1.1;
		    }

		    // H N
		    else if (((iso.atoms()[i][0].element().symbol() == "N") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "N")))
		    {
				min = 0.7;
				max = 1.2;
		    }

		    // Se H
		    else if (((iso.atoms()[i][0].element().symbol() == "Se") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "Se")))
		    {
				min = 1.2;
				max = 1.7;
		    }

		    // K Br
		    else if (((iso.atoms()[i][0].element().symbol() == "Br") && \
				(iso.atoms()[j][0].element().symbol() == "K")) || \
		        ((iso.atoms()[i][0].element().symbol() == "K") && \
				(iso.atoms()[j][0].element().symbol() == "Br")))
		    {
				min = 2.6;
				max = 3.0;
		    }

		    // Li Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Li") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "Li")))
		    {
				min = 1.8;
				max = 2.4;
		    }

		    // Li Li
		    else if ((iso.atoms()[i][0].element().symbol() == "Li") && \
				(iso.atoms()[j][0].element().symbol() == "Li"))
		    {
				min = 2.5;
				max = 2.9;
		    }

		    // Mg F
		    else if (((iso.atoms()[i][0].element().symbol() == "Mg") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "Mg")))
		    {
				min = 1.5;
				max = 1.9;
		    }

		    // Mg S
		    else if (((iso.atoms()[i][0].element().symbol() == "Mg") && \
				(iso.atoms()[j][0].element().symbol() == "S")) || \
		        ((iso.atoms()[i][0].element().symbol() == "S") && \
				(iso.atoms()[j][0].element().symbol() == "Mg")))
		    {
				min = 1.9;
				max = 2.4;
		    }

		    // Na H
		    else if (((iso.atoms()[i][0].element().symbol() == "Na") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "Na")))
		    {
				min = 1.6;
				max = 2.1;
		    }

		    // Na O
		    else if (((iso.atoms()[i][0].element().symbol() == "Na") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "Na")))
		    {
				min = 1.7;
				max = 2.1;
		    }

		    // O S
		    else if (((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "S")) || \
		        ((iso.atoms()[i][0].element().symbol() == "S") && \
				(iso.atoms()[j][0].element().symbol() == "O")))
		    {
				min = 1.2;
				max = 1.6;
		    }

		    // P S
		    else if (((iso.atoms()[i][0].element().symbol() == "P") && \
				(iso.atoms()[j][0].element().symbol() == "S")) || \
		        ((iso.atoms()[i][0].element().symbol() == "S") && \
				(iso.atoms()[j][0].element().symbol() == "P")))
		    {
				min = 1.7;
				max = 2.1;
		    }

		    // Se Se
		    else if ((iso.atoms()[i][0].element().symbol() == "Se") && \
				(iso.atoms()[j][0].element().symbol() == "Se"))
		    {
				min = 1.9;
				max = 2.3;
		    }

		    // Si Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Si") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "Si")))
		    {
				min = 1.8;
				max = 2.2;
		    }

		    // S S
		    else if ((iso.atoms()[i][0].element().symbol() == "S") && \
				(iso.atoms()[j][0].element().symbol() == "S"))
		    {
				min = 1.7;
				max = 2.2;
		    }

		    // Zn C
		    else if (((iso.atoms()[i][0].element().symbol() == "Zn") && \
				(iso.atoms()[j][0].element().symbol() == "C")) || \
		        ((iso.atoms()[i][0].element().symbol() == "C") && \
				(iso.atoms()[j][0].element().symbol() == "Zn")))
		    {
				min = 1.7;
				max = 2.1;
		    }

		    // Al H
		    else if (((iso.atoms()[i][0].element().symbol() == "Al") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "Al")))
		    {
				min = 1.4;
				max = 1.9;
		    }

		    // B B
		    else if ((iso.atoms()[i][0].element().symbol() == "B") && \
				(iso.atoms()[j][0].element().symbol() == "B"))
		    {
				min = 1.4;
				max = 1.9;
		    }

		    // Be F
		    else if (((iso.atoms()[i][0].element().symbol() == "Be") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "Be")))
		    {
				min = 1.2;
				max = 1.6;
		    }

		    // B F
		    else if (((iso.atoms()[i][0].element().symbol() == "B") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "B")))
		    {
				min = 1.1;
				max = 1.5;
		    }

		    // Br Br
		    else if ((iso.atoms()[i][0].element().symbol() == "Br") && \
				(iso.atoms()[j][0].element().symbol() == "Br"))
		    {
				min = 2.1;
				max = 2.5;
		    }

		    // Br C
		    else if (((iso.atoms()[i][0].element().symbol() == "Br") && \
				(iso.atoms()[j][0].element().symbol() == "C")) || \
		        ((iso.atoms()[i][0].element().symbol() == "C") && \
				(iso.atoms()[j][0].element().symbol() == "Br")))
		    {
				min = 1.6;
				max = 2.1;
		    }

		    // Ge C
		    else if (((iso.atoms()[i][0].element().symbol() == "Ge") && \
				(iso.atoms()[j][0].element().symbol() == "C")) || \
		        ((iso.atoms()[i][0].element().symbol() == "C") && \
				(iso.atoms()[j][0].element().symbol() == "Ge")))
		    {
				min = 1.8;
				max = 2.2;
		    }

		    // Cl Ca
		    else if (((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "Ca")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Ca") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")))
		    {
				min = 2.2;
				max = 2.6;
		    }

		    // C N
		    else if (((iso.atoms()[i][0].element().symbol() == "C") && \
				(iso.atoms()[j][0].element().symbol() == "N")) || \
		        ((iso.atoms()[i][0].element().symbol() == "N") && \
				(iso.atoms()[j][0].element().symbol() == "C")))
		    {
				min = 1.0;
				max = 1.6;
		    }

		    // C Se
		    else if (((iso.atoms()[i][0].element().symbol() == "Se") && \
				(iso.atoms()[j][0].element().symbol() == "C")) || \
		        ((iso.atoms()[i][0].element().symbol() == "C") && \
				(iso.atoms()[j][0].element().symbol() == "Se")))
		    {
				min = 1.5;
				max = 2.1;
		    }

		    // Cu O
		    else if (((iso.atoms()[i][0].element().symbol() == "Cu") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "Cu")))
		    {
				min = 1.5;
				max = 1.9;
		    }

		    // F Ca
		    else if (((iso.atoms()[i][0].element().symbol() == "Ca") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "Ca")))
		    {
				min = 1.8;
				max = 2.2;
		    }

		    // F P
		    else if (((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "P")) || \
		        ((iso.atoms()[i][0].element().symbol() == "P") && \
				(iso.atoms()[j][0].element().symbol() == "F")))
		    {
				min = 1.3;
				max = 1.7;
		    }

		    // H Ar
		    else if (((iso.atoms()[i][0].element().symbol() == "Ar") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "Ar")))
		    {
				min = 1.1;
				max = 1.5;
		    }

		    // He He
		    else if ((iso.atoms()[i][0].element().symbol() == "He") && \
				(iso.atoms()[j][0].element().symbol() == "He"))
		    {
				min = 0.8;
				max = 1.2;
		    }

		    // H I
		    else if (((iso.atoms()[i][0].element().symbol() == "I") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "I")))
		    {
				min = 1.4;
				max = 1.8;
		    }

		    // H O
		    else if (((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "H")) || \
		        ((iso.atoms()[i][0].element().symbol() == "H") && \
				(iso.atoms()[j][0].element().symbol() == "O")))
		    {
				min = 0.8;
				max = 1.2;
		    }

		    // K Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "K")) || \
		        ((iso.atoms()[i][0].element().symbol() == "K") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")))
		    {
				min = 2.5;
				max = 2.9;
		    }

		    // Li F
		    else if (((iso.atoms()[i][0].element().symbol() == "Li") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "Li")))
		    {
				min = 1.4;
				max = 1.8;
		    }

		    // Na Li
		    else if (((iso.atoms()[i][0].element().symbol() == "Na") && \
				(iso.atoms()[j][0].element().symbol() == "Li")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Li") && \
				(iso.atoms()[j][0].element().symbol() == "Na")))
		    {
				min = 2.7;
				max = 3.1;
		    }

		    // Na Br
		    else if (((iso.atoms()[i][0].element().symbol() == "Na") && \
				(iso.atoms()[j][0].element().symbol() == "Br")) || \
		        ((iso.atoms()[i][0].element().symbol() == "Br") && \
				(iso.atoms()[j][0].element().symbol() == "Na")))
		    {
				min = 2.3;
				max = 2.7;
		    }

		    // Na K
		    else if (((iso.atoms()[i][0].element().symbol() == "Na") && \
				(iso.atoms()[j][0].element().symbol() == "K")) || \
		        ((iso.atoms()[i][0].element().symbol() == "K") && \
				(iso.atoms()[j][0].element().symbol() == "Na")))
		    {
				min = 3.4;
				max = 3.8;
		    }

		    // N Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "N")) || \
		        ((iso.atoms()[i][0].element().symbol() == "N") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")))
		    {
				min = 1.4;
				max = 2.2;
		    }

		    // N N
		    else if ((iso.atoms()[i][0].element().symbol() == "N") && \
				(iso.atoms()[j][0].element().symbol() == "N"))
		    {
				min = 0.8;
				max = 2.4;
		    }

		    // Ca O
		    else if (((iso.atoms()[i][0].element().symbol() == "Ca") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "Ca")))
		    {
				min = 1.8;
				max = 2.2;
		    }

		    // Si O
		    else if (((iso.atoms()[i][0].element().symbol() == "Si") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "Si")))
		    {
				min = 1.3;
				max = 1.7;
		    }

		    // P N
		    else if (((iso.atoms()[i][0].element().symbol() == "P") && \
				(iso.atoms()[j][0].element().symbol() == "N")) || \
		        ((iso.atoms()[i][0].element().symbol() == "N") && \
				(iso.atoms()[j][0].element().symbol() == "P")))
		    {
				min = 1.3;
				max = 1.7;
		    }

		    // S Cl
		    else if (((iso.atoms()[i][0].element().symbol() == "Cl") && \
				(iso.atoms()[j][0].element().symbol() == "S")) || \
		        ((iso.atoms()[i][0].element().symbol() == "S") && \
				(iso.atoms()[j][0].element().symbol() == "Cl")))
		    {
				min = 1.8;
				max = 2.3;
		    }

		    // S F
		    else if (((iso.atoms()[i][0].element().symbol() == "S") && \
				(iso.atoms()[j][0].element().symbol() == "F")) || \
		        ((iso.atoms()[i][0].element().symbol() == "F") && \
				(iso.atoms()[j][0].element().symbol() == "S")))
		    {
				min = 1.4;
				max = 1.8;
		    }

		    // Si S
		    else if (((iso.atoms()[i][0].element().symbol() == "Si") && \
				(iso.atoms()[j][0].element().symbol() == "S")) || \
		        ((iso.atoms()[i][0].element().symbol() == "S") && \
				(iso.atoms()[j][0].element().symbol() == "Si")))
		    {
				min = 1.7;
				max = 2.1;
		    }

		    // S Se
		    else if (((iso.atoms()[i][0].element().symbol() == "Se") && \
				(iso.atoms()[j][0].element().symbol() == "S")) || \
		        ((iso.atoms()[i][0].element().symbol() == "S") && \
				(iso.atoms()[j][0].element().symbol() == "Se")))
		    {
				min = 1.8;
				max = 2.2;
		    }

		    // Nb O
		    else if (((iso.atoms()[i][0].element().symbol() == "Nb") && \
				(iso.atoms()[j][0].element().symbol() == "O")) || \
		        ((iso.atoms()[i][0].element().symbol() == "O") && \
				(iso.atoms()[j][0].element().symbol() == "Nb")))
		    {
				min = 1.6;
				max = 2.1;
		    }
			
			// Did not find elements
			if (min < 0)
				continue;
			
			// Loop over pairs of atoms of current element
			found = false;
			for (k = 0; k < iso.atoms()[i].length(); ++k)
			{
				for (m = 0; m < iso.atoms()[j].length(); ++m)
				{
					
					// Get distance of atom with itself
					if (iso.atoms()[i][k].atomNumber() == iso.atoms()[j][m].atomNumber())
						dis = iso.basis().secondDistance(iso.atoms()[i][k].fractional(), FRACTIONAL, \
							iso.atoms()[j][m].fractional(), FRACTIONAL);
					
					// Get distance between atoms
					else
						dis = iso.basis().distance(iso.atoms()[i][k].fractional(), FRACTIONAL, \
							iso.atoms()[j][m].fractional(), FRACTIONAL);
					
					// Found bond
					if ((dis >= min) && (dis <= max))
					{
						found = true;
						break;
					}
				}
				if (found)
					break;
			}
			
			// Found pair
			if (found)
			{
				res.add();
				res.last().set(iso.atoms()[i][0].element(), iso.atoms()[j][0].element(), min, max);
			}
		}
	}
	
	// Return bonds
	return res;
}

