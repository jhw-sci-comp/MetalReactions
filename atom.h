#ifndef ATOM_T
#define ATOM_T

#include <string>
#include <cstring>
#include <iostream>


/* represents an Atom of an element */
class Atom {
public:
	Atom();
	Atom(std::string, std::string, int);
	Atom(const Atom&);
	~Atom();
	std::string getElementName() const;
	std::string getElementSymbol() const;
	int getValenceElectrons() const;
	std::string getType() const;
	float getElectronegativity() const;
	unsigned int getAtomicNumber() const;
	void setElementName(std::string);
	void setElementSymbol(std::string);
	void setValenceElectrons(int);
	void setType(std::string);
	void setAtomicNumber(unsigned int);
	Atom& operator=(const Atom&);


private:
	std::string elementSymbol = "";
	std::string elementName = "";
	int valenceElectrons = 0;
	std::string type = "";
	float electronegativity = 0.0;
	unsigned int atomicNumber = 0;
};

#endif
