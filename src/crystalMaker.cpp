/* crystalMaker.cpp -- Structure for CrystalMaker files
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */



#include "crystalMaker.h"
#include "bonds.h"
#include "language.h"
#include "num.h"
#include "output.h"



/* void CrystalMaker::write(const Word& file, const ISO& iso)
 *
 * Write structure to crystalmaker file
 */

void CrystalMaker::write(const Word& file, const ISO& iso)
{
	
	// Precision for printing numbers
	int prec = 14;
	
	// Setup output if file was set
	int origStream = Output::streamID();
	PrintMethod origMethod = Output::method();
	if (file.length() > 0)
	{
		
		// Open file for writing if needed
		if (file != "stdout")
			Output::setStream(Output::addStream(file));

		// Set output method
		Output::method(STANDARD);
	}
	
	// Write cell
	int i;
	Output::newline();
    Output::print("CELL ");
    for (i = 0; i < 3; ++i)
	{
		Output::print(iso.basis().lengths()[i], 8);
		Output::print(" ");
	}
	for (i = 0; i < 3; ++i)
	{
		Output::print(Num<double>::toDegrees(iso.basis().angles()[i]), 8);
		Output::print(" ");
	}
	
	// Write lattice type
	Output::newline();
	Output::print("LATC P");
    
    // Write symmetry operations
	Output::newline();
    Output::print("SYMM +x +y +z");
    
    // Write plot range
	Output::newline();
    Output::print("XYZR 0.0 1.0 0.0 1.0 0.0 1.0");
    
    // Write model, 1 for ball and stick, 3 for polyhedral
	Output::newline();
    Output::print("MODL 1");
    
    // Write background color
    Output::newline();
	Output::print("BKCL 1.0 1.0 1.0");
	
	// Get bonds
	OList<Bond> bonds = Bonds::find(iso);
	
	// Print bonds
	for (i = 0; i < bonds.length(); ++i)
	{
		Output::newline();
		Output::print("BMAX ");
		Output::print(bonds[i].elem1().symbol());
		Output::print(" ");
		Output::print(bonds[i].elem2().symbol());
		Output::print(" ");
		Output::print(bonds[i].max());
	}
	
	// Write atom types
    Output::newline();
	Output::print("TYPE");
    for (i = 0; i < iso.atoms().length(); ++i)
    {
		Output::newline();
        if (iso.atoms()[i][0].element().number() == 1)
            Output::print("H  0.10 1.00000 0.74243 0.74243");
        else if (iso.atoms()[i][0].element().number() == 2)
            Output::print("He 1.00 0.20000 0.20000 0.20000");
        else if (iso.atoms()[i][0].element().number() == 3)
            Output::print("Li 0.90 0.32311 0.53387 0.69078");
        else if (iso.atoms()[i][0].element().number() == 4)
            Output::print("Be 0.41 0.61572 0.99997 0.61050");
        else if (iso.atoms()[i][0].element().number() == 5)
            Output::print("B  0.25 0.53341 0.53341 0.71707");
        else if (iso.atoms()[i][0].element().number() == 6)
            Output::print("C  0.77 0.06577 0.02538 0.00287");
        else if (iso.atoms()[i][0].element().number() == 7)
            Output::print("N  0.30 0.50660 0.68658 1.00000");
        else if (iso.atoms()[i][0].element().number() == 8)
            Output::print("O  1.21 1.00000 0.00000 0.00000");
        else if (iso.atoms()[i][0].element().number() == 9)
            Output::print("F  1.19 0.00000 0.80000 0.00000");
        else if (iso.atoms()[i][0].element().number() == 10)
            Output::print("Ne 1.00 0.20000 0.20000 0.20000");
        else if (iso.atoms()[i][0].element().number() == 11)
            Output::print("Na 1.16 0.78946 0.77423 0.00002");
        else if (iso.atoms()[i][0].element().number() == 12)
            Output::print("Mg 0.86 1.00000 0.78787 0.01648");
        else if (iso.atoms()[i][0].element().number() == 13)
            Output::print("Al 0.53 0.09751 0.67741 1.00000");
        else if (iso.atoms()[i][0].element().number() == 14)
            Output::print("Si 0.40 0.00000 0.00000 1.00000");
        else if (iso.atoms()[i][0].element().number() == 15)
            Output::print("P  0.31 0.38148 0.38148 0.38148");
        else if (iso.atoms()[i][0].element().number() == 16)
            Output::print("S  0.26 1.00000 0.98070 0.00000");
        else if (iso.atoms()[i][0].element().number() == 17)
            Output::print("Cl 1.67 0.07999 1.00000 0.05585");
        else if (iso.atoms()[i][0].element().number() == 18)
            Output::print("Ar 1.00 0.20000 0.20000 0.20000");
        else if (iso.atoms()[i][0].element().number() == 19)
            Output::print("K  1.52 0.45776 0.00000 1.00000");
        else if (iso.atoms()[i][0].element().number() == 20)
            Output::print("Ca 1.14 0.36245 0.61630 0.77632");
        else if (iso.atoms()[i][0].element().number() == 21)
            Output::print("Sc 0.89 0.98804 0.41819 1.00000");
        else if (iso.atoms()[i][0].element().number() == 22)
            Output::print("Ti 0.75 0.47236 0.79393 1.00000");
        else if (iso.atoms()[i][0].element().number() == 23)
            Output::print("V  1.00 0.20000 0.20000 0.20000");
        else if (iso.atoms()[i][0].element().number() == 24)
            Output::print("Cr 0.76 1.00000 0.00000 0.60000");
        else if (iso.atoms()[i][0].element().number() == 25)
            Output::print("Mn 0.81 1.00000 0.00000 0.60000");
        else if (iso.atoms()[i][0].element().number() == 26)
            Output::print("Fe 0.69 0.71051 0.44662 0.00136");
        else if (iso.atoms()[i][0].element().number() == 27)
            Output::print("Co 0.54 0.00000 0.00000 0.68666");
        else if (iso.atoms()[i][0].element().number() == 28)
            Output::print("Ni 0.70 0.86664 0.86664 0.86664");
        else if (iso.atoms()[i][0].element().number() == 29)
            Output::print("Cu 0.71 0.00000 0.00000 1.00000");
        else if (iso.atoms()[i][0].element().number() == 30)
            Output::print("Zn 0.74 0.56123 0.56445 0.50799");
        else if (iso.atoms()[i][0].element().number() == 31)
            Output::print("Ga 0.76 0.62292 0.89293 0.45486");
        else if (iso.atoms()[i][0].element().number() == 32)
            Output::print("Ge 0.53 0.49557 0.43499 0.65193");
        else if (iso.atoms()[i][0].element().number() == 33)
            Output::print("As 0.72 0.28579 0.16544 0.31698");
        else if (iso.atoms()[i][0].element().number() == 34)
            Output::print("Se 0.56 0.60420 0.93874 0.06122");
        else if (iso.atoms()[i][0].element().number() == 35)
            Output::print("Br 1.82 0.49644 0.19332 0.01076");
        else if (iso.atoms()[i][0].element().number() == 36)
            Output::print("Kr 1.00  0.20000 0.20000 0.20000");
        else if (iso.atoms()[i][0].element().number() == 37)
            Output::print("Rb 1.66 1.00000 0.00000 0.60000");
        else if (iso.atoms()[i][0].element().number() == 38)
            Output::print("Sr 1.32 0.00000 1.00000 0.15256");
        else if (iso.atoms()[i][0].element().number() == 39)
            Output::print("Y  1.04 0.40258 0.59738 0.55813");
        else if (iso.atoms()[i][0].element().number() == 40)
            Output::print("Zr 0.86 0.00000 1.00000 0.00000");
        else if (iso.atoms()[i][0].element().number() == 41)
            Output::print("Nb 0.78 0.29992 0.70007 0.46458");
        else if (iso.atoms()[i][0].element().number() == 42)
            Output::print("Mo 0.79 0.70584 0.52601 0.68925");
        else if (iso.atoms()[i][0].element().number() == 43)
            Output::print("Tc 0.79 0.80574 0.68698 0.79477");
        else if (iso.atoms()[i][0].element().number() == 44)
            Output::print("Ru 0.82 0.81183 0.72113 0.68089");
        else if (iso.atoms()[i][0].element().number() == 45)
            Output::print("Rh 0.81 0.80748 0.82205 0.67068");
        else if (iso.atoms()[i][0].element().number() == 46)
            Output::print("Pd 0.78 0.75978 0.76817 0.72453");
        else if (iso.atoms()[i][0].element().number() == 47)
            Output::print("Ag 1.29 0.72032 0.73631 0.74339");
        else if (iso.atoms()[i][0].element().number() == 48)
            Output::print("Cd 0.92 0.55711 0.50755 1.00000");
        else if (iso.atoms()[i][0].element().number() == 49)
            Output::print("In 0.94 0.84378 0.50401 0.73483");
        else if (iso.atoms()[i][0].element().number() == 50)
            Output::print("Sn 0.69 0.60763 0.56051 0.72926");
        else if (iso.atoms()[i][0].element().number() == 51)
            Output::print("Sb 0.90 0.84627 0.51498 0.31315");
        else if (iso.atoms()[i][0].element().number() == 52)
            Output::print("Te 1.11 0.67958 0.63586 0.32038");
        else if (iso.atoms()[i][0].element().number() == 53)
            Output::print("I  2.06 0.33724 0.00272 0.50657");
        else if (iso.atoms()[i][0].element().number() == 54)
            Output::print("Xe 0.62 0.60661 0.63217 0.97304");
        else if (iso.atoms()[i][0].element().number() == 55)
            Output::print("Cs 1.81 1.00000 0.00000 0.60000");
        else if (iso.atoms()[i][0].element().number() == 56)
            Output::print("Ba 1.49 0.74342 0.39631 0.45338");
        else if (iso.atoms()[i][0].element().number() == 57)
            Output::print("La 1.36 0.35340 0.77057 0.28736");
        else if (iso.atoms()[i][0].element().number() == 58)
            Output::print("Ce 1.15 0.82054 0.99071 0.02373");
        else if (iso.atoms()[i][0].element().number() == 59)
            Output::print("Pr 1.32 0.99129 0.88559 0.02315");
        else if (iso.atoms()[i][0].element().number() == 60)
            Output::print("Nd 1.30 0.98700 0.55560 0.02744");
        else if (iso.atoms()[i][0].element().number() == 61)
            Output::print("Pm 1.28 0.49999 0.49999 0.49999");
        else if (iso.atoms()[i][0].element().number() == 62)
            Output::print("Sm 1.10 0.99042 0.02402 0.49194");
        else if (iso.atoms()[i][0].element().number() == 63)
            Output::print("Eu 1.31 0.98366 0.03078 0.83615");
        else if (iso.atoms()[i][0].element().number() == 64)
            Output::print("Gd 1.08 0.75325 0.01444 1.00000");
        else if (iso.atoms()[i][0].element().number() == 65)
            Output::print("Tb 1.18 0.44314 0.01662 0.99782");
        else if (iso.atoms()[i][0].element().number() == 66)
            Output::print("Dy 1.05 0.19390 0.02373 0.99071");
        else if (iso.atoms()[i][0].element().number() == 67)
            Output::print("Ho 1.04 0.02837 0.25875 0.98607");
        else if (iso.atoms()[i][0].element().number() == 68)
            Output::print("Er 1.03 0.28687 0.45071 0.23043");
        else if (iso.atoms()[i][0].element().number() == 69)
            Output::print("Tm 1.02 0.15903 0.79509 0.98584");
        else if (iso.atoms()[i][0].element().number() == 70)
            Output::print("Yb 1.13 0.15322 0.99164 0.95836");
        else if (iso.atoms()[i][0].element().number() == 71)
            Output::print("Lu 1.00 0.15096 0.99390 0.71031");
        else if (iso.atoms()[i][0].element().number() == 72)
            Output::print("Hf 0.85 0.70703 0.70552 0.35090");
        else if (iso.atoms()[i][0].element().number() == 73)
            Output::print("Ta 0.78 0.71951 0.60693 0.33840");
        else if (iso.atoms()[i][0].element().number() == 74)
            Output::print("W  0.74 0.55615 0.54257 0.50178");
        else if (iso.atoms()[i][0].element().number() == 75)
            Output::print("Re 0.77 0.70294 0.69400 0.55789");
        else if (iso.atoms()[i][0].element().number() == 76)
            Output::print("Os 0.77 0.78703 0.69511 0.47378");
        else if (iso.atoms()[i][0].element().number() == 77)
            Output::print("Ir 0.77 0.78975 0.81032 0.45048");
        else if (iso.atoms()[i][0].element().number() == 78)
            Output::print("Pt 0.74 0.76088 0.76088 0.76088");
        else if (iso.atoms()[i][0].element().number() == 79)
            Output::print("Au 1.51 0.99628 0.70149 0.22106");
        else if (iso.atoms()[i][0].element().number() == 80)
            Output::print("Hg 0.83 0.82939 0.72125 0.79823");
        else if (iso.atoms()[i][0].element().number() == 81)
            Output::print("Tl 1.03 0.58798 0.53854 0.42649");
        else if (iso.atoms()[i][0].element().number() == 82)
            Output::print("Pb 1.49 0.32386 0.32592 0.35729");
        else if (iso.atoms()[i][0].element().number() == 83)
            Output::print("Bi 1.17 0.62942 0.42309 0.73683");
        else if (iso.atoms()[i][0].element().number() == 84)
            Output::print("Po 1.08 0.49999 0.49999 0.49999");
        else if (iso.atoms()[i][0].element().number() == 85)
            Output::print("At 0.76 0.49999 0.49999 0.49999");
        else if (iso.atoms()[i][0].element().number() == 86)
            Output::print("Rn 1.00 0.20000 0.20000 0.20000");
        else if (iso.atoms()[i][0].element().number() == 87)
            Output::print("Fr 1.94 0.49999 0.49999 0.49999");
        else if (iso.atoms()[i][0].element().number() == 88)
            Output::print("Ra 1.62 0.42959 0.66658 0.34786");
        else if (iso.atoms()[i][0].element().number() == 89)
            Output::print("Ac 1.26 0.39344 0.62100 0.45034");
        else if (iso.atoms()[i][0].element().number() == 90)
            Output::print("Th 1.19 0.14893 0.99596 0.47105");
        else if (iso.atoms()[i][0].element().number() == 91)
            Output::print("Pa 1.09 0.16100 0.98386 0.20855");
        else if (iso.atoms()[i][0].element().number() == 92)
            Output::print("U  0.87 0.47773 0.63362 0.66714");
        else if (iso.atoms()[i][0].element().number() == 93)
            Output::print("Np 1.00 0.20000 0.20000 0.20000");
        else if (iso.atoms()[i][0].element().number() == 94)
            Output::print("Pu 1.00 0.36564 0.36431 0.67912");
        else if (iso.atoms()[i][0].element().number() == 95)
            Output::print("Am 1.12 0.49999 0.49999 0.49999");
        else if (iso.atoms()[i][0].element().number() == 96)
            Output::print("Cm 1.11 0.49999 0.49999 0.49999");
        else if (iso.atoms()[i][0].element().number() == 97)
            Output::print("Bk 1.00 0.20000 0.20000 0.20000");
        else
        {
			Output::print(iso.atoms()[i][0].element().symbol());
			Output::print("  ");
            Output::print("1.0 0.20000 0.20000 0.20000");
        }
    }
	
	// Setup coordinates
	int j, k;
	Output message;
	message.addLines(iso.numAtoms());
	for (i = 0; i < iso.atoms().length(); ++i)
	{
		for (j = 0; j < iso.atoms()[i].length(); ++j)
		{
			message.addLine();
			message.add(iso.atoms()[i][j].element().symbol());
			message.add(iso.atoms()[i][j].element().symbol() + \
				Language::numberToWord(iso.atoms()[i][j].atomNumber() + 1));
			for (k = 0; k < 3; ++k)
				message.add(iso.atoms()[i][j].fractional()[k], prec);
		}
	}
	
	// Make alignment
	List<PrintAlign> align(5);
	align[0] = align[1] = LEFT;
	align[2] = align[3] = align[4] = RIGHT;
	
	// Print coordinates
	Output::newline();
	Output::print("ATOM");
	Output::newline();
	Output::print(message, align);
	
	// Reset output if file was set
	if (file.length() > 0)
	{
		if (file != "stdout")
			Output::removeStream(Output::streamID());
		Output::setStream(origStream);
		Output::method(origMethod);
	}
}
