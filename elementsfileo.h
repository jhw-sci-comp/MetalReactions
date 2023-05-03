#ifndef ELEMENTSFILEO_H
#define ELEMENTSFILEO_H

#include <fstream>
#include <list>
#include <vector>
#include <algorithm>
#include "Atom.h"
#include "elementsfileexception.h"

class ElementsFileO {
public:
	ElementsFileO();
	bool selectAtom(unsigned int, Atom&);
	bool selectAtom(std::string, Atom&);
	void selectListOfAtoms(std::list<std::string>);
	void selectListOfAtoms(std::string);
	std::vector<Atom> getSelectedAtoms() const;
	void deleteSelectedAtoms();


private:
	std::vector<Atom> selectedAtoms;
};


#endif