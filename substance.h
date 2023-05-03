#ifndef SUBSTANCE_H
#define SUBSTANCE_H

#include <vector>
#include "atom.h"


class Substance {
public:
	Substance();
	Substance(std::vector<Atom>, std::vector<int>);
	Substance(const Substance&);
	~Substance();
	std::vector<Atom> getAtoms() const;
	std::vector<int> getNumberOfAtoms() const;
	void addAtoms(std::vector<Atom>);
	void addAtoms(Atom);
	void addNumberOfAtoms(std::vector<int>);
	void addNumberOfAtoms(int);
	Substance& operator=(const Substance&);

private:
	std::vector<Atom> atoms;
	std::vector<int> numberOfAtoms;
};


bool operator== (const Substance&, const Substance&); //used to find substance in vector

#endif
