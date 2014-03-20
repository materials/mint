/* output.h -- Print a formatted message
 *
 * Copyright (C) 2011-2013 by Northwestern University, All Rights Reserved
 * 
 * Contact: Kyle Michel (kylemichel@gmail.com)
 */


#ifndef OUTPUT_H
#define OUTPUT_H



#include "num.h"
#include "list.h"
#include "text.h"
#include <iostream>
#include <fstream>
using namespace std;



// Define properties
enum PrintType {ORDINARY, WARNING, ERROR};
enum PrintMethod {STANDARD, RESTRICTED};
enum PrintAlign {LEFT, RIGHT};



// Class to store information about a stream
class StreamInfo
{
	
	// Variables
	bool _streamIsSet;
	ostream* _stream;
	ofstream _outfile;
	Word _outfileName;
	int _numLinesBefore;
	int _numLinesAfter;
	int _havePrinted;
	bool _hold;
	PrintMethod _method;
	PrintType _type;
	bool _showTermination;
	bool _ordinaryOn;
	bool _warningsOn;
	bool _errorsOn;
	bool _print;
	int _maxLevel;
	int _spacesPerLevel;
	int _currentLevel;
	List<bool> _prevOrdinaryOn;
	List<bool> _prevWarningsOn;
	List<bool> _prevErrorsOn;
	
	// Functions
	void updatePrint();
	
public:
	
	// Constructor and destructor
	StreamInfo();
	~StreamInfo();
	
	// Startup function
	void setup(ostream* input, int numBlankLinesAtStart = 0, int numBlankLinesAtEnd = 1);
	void setup(const Word& file, int numBlankLinesAtStart = 0, int numBlankLinesAtEnd = 1);
	
	// Functions
	void quietOn(bool warningsOff, bool errorsOff);
	void quietOff();
	
	// Set functions
	void havePrinted(int input)			{     _havePrinted = input; }
	void hold(bool input)				{            _hold = input; }
	void method(PrintMethod input)		{          _method = input; updatePrint(); }
	void type(PrintType input)			{            _type = input; updatePrint(); }
	void showTermination(bool input)	{ _showTermination = input; }
	void ordinaryOn(bool input)			{      _ordinaryOn = input; updatePrint(); }
	void warningsOn(bool input)			{      _warningsOn = input; updatePrint(); }
	void errorsOn(bool input)			{        _errorsOn = input; updatePrint(); }
	void maxLevel(int input)			{        _maxLevel = input; }
	void spacesPerLevel(int input)		{  _spacesPerLevel = input; }
	void currentLevel(int input)		{    _currentLevel = input; updatePrint(); }
	
	// Access functions
	bool print() const				{ return _print; }
	bool hold() const				{ return _hold; }
	bool ordinaryOn() const			{ return _ordinaryOn; }
	bool warningsOn() const			{ return _warningsOn; }
	bool errorsOn() const			{ return _errorsOn; }
	bool showTermination() const	{ return _showTermination; }
	bool streamIsSet() const		{ return _streamIsSet; }
	int havePrinted() const			{ return _havePrinted; }
	int numLinesBefore() const		{ return _numLinesBefore; }
	int currentLevel() const		{ return _currentLevel; }
	int maxLevel() const			{ return _maxLevel; }
	int spacesPerLevel() const		{ return _spacesPerLevel; }
	PrintMethod method() const		{ return _method; }
	PrintType type() const			{ return _type; }
	ostream& stream() const			{ return *_stream; }
};



// Class for output
class Output
{
	
	// Variables
	bool _addSpaces;
	Text _data;
	
	// Functions
	static double checkZero(double in, int prec);
	static int roundToOne(double in);
	
	// Static helper variables
	static char _arg[10];
	static char _buffer[50];
	
	// Static variables to control run time output
	static List<int> _ids;
	static OList<StreamInfo> _streams;
	static int _primary;
	static int _curStream;
	
public:
	
	// Constructor
	Output() { _addSpaces = true; }
	
	// Initialize
	static void initialize();
	
	// Functions
	void clear()								{ _data.clear(); }
	void addSpaces(bool input)					{ _addSpaces = input; }
	void add(int input);
	void add(double input, int prec = 8);
	void addSci(double input, int prec = 8);
	void add(const char* input)					{ _data.addWord(input); }
	void add(const Word& input)					{ _data.addWord(input); }
	void add(const Words& input)				{ _data.addWords(input); }
	void addTab();
	void addLine()								{ _data.addLine(); }
	void addLine(const char* line)				{ _data.addLine(line); }
	void addLines(int numToAdd)					{ _data.addLines(numToAdd); }
	void addWords(int numToAdd)					{ _data.addWords(numToAdd); }
	void addWords(int line, int numToAdd)		{ _data.addWords(line, numToAdd); }
	int numLines() const						{ return _data.length(); }
   
	// Add and remove streams
	static int addStream(ostream* input, int numBlankLinesAtStart = 0, int numBlankLinesAtEnd = 1);
	static int addStream(const Word& input, int numBlankLinesAtStart = 0, int numBlankLinesAtEnd = 1);
	static void removeStream(int inID);
	
	// Change between streams
	static void setPrimary() { _curStream = _primary; }
	static void setPrimary(int inID);
	static void setStream(int inID);
	
	// Direct access to streams
	static ostream& stream() { return _streams[_curStream].stream(); }
	
	// Access information about streams
	static bool ordinaryOn() 		{ return _streams[_curStream].ordinaryOn(); }
	static bool warningsOn() 		{ return _streams[_curStream].warningsOn(); }
	static bool errorsOn() 			{ return _streams[_curStream].errorsOn(); }
	static int streamID() 			{ return _ids[_curStream]; }
	static int primaryStreamID() 	{ return _ids[_primary]; }
	static PrintMethod method()		{ return _streams[_curStream].method(); }
	static int level() 				{ return _streams[_curStream].currentLevel(); }
	static List<bool> onOff();
	
	// Set properties of streams
	static void maxLevel(int input) 			{ _streams[_curStream].maxLevel(input); }
	static void spacesPerLevel(int input) 		{ _streams[_curStream].spacesPerLevel(input); }
    static void ordinaryOn(bool input) 			{ _streams[_curStream].ordinaryOn(input); }
    static void warningsOn(bool input) 			{ _streams[_curStream].warningsOn(input); }
	static void errorsOn(bool input) 			{ _streams[_curStream].errorsOn(input); }
	static void method(PrintMethod input);
	static void level(int input) 				{ _streams[_curStream].currentLevel(input); }
	static void onOff(const List<bool>& input)	{ ordinaryOn(input[0]); warningsOn(input[1]); errorsOn(input[2]); }
	static void quietOn(bool warningsOff = false, bool errorsOff = false)
		{ _streams[_curStream].quietOn(warningsOff, errorsOff); }
	static void quietOff()						{ _streams[_curStream].quietOff(); }
	
	// Static print functions
	static void increase() { _streams[_curStream].currentLevel(_streams[_curStream].currentLevel() + 1); }
    static void decrease() { _streams[_curStream].currentLevel(_streams[_curStream].currentLevel() - 1); }
    static void newline(PrintType input = ORDINARY);
	static void tab();
	static void print(char message);
	static void print(const char* messsage);
	static void print(const Word& message)	{ print(message.array()); }
    static void print(int message);
	static void print(unsigned long int message);
    static void print(double message, int prec = -1);
	static void printSci(double message, int prec = -1);
	static void print(const Words& message, bool useComma = true, bool useAnd = true);
    static void print(const List<int>& message, bool useComma = true, bool useAnd = true);
    static void print(const List<double>& message, int prec, bool useComma = true, bool useAnd = true);
	static void print(const Vector3D& message, int prec, bool useComma = true);
	static void print(const Output& message, PrintAlign align, PrintType input = ORDINARY);
	static void print(const Output& message, const List<PrintAlign>& align, PrintType input = ORDINARY);
	
	// Static print functions with padding
	static void printPadded(int message, int width, PrintAlign align);
	static void printPadded(double message, int width, PrintAlign align, int prec = -1);
	static void printPaddedSci(double message, int width, PrintAlign align, int prec = -1);

	// Static exit functions
	static void turnOffExitFun() { _streams[_curStream].showTermination(false); }
	static void quit();
	
	// Friends
	friend class StreamInfo;
	friend class OutputHelper;
};



#endif
