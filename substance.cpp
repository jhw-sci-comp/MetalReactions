#include "substance.h"


Substance::Substance() {

}

Substance::Substance(std::vector<Atom> atoms_, std::vector<int> numberOfAtoms_) {
	for (auto a : atoms_) {
		this->atoms.push_back(a);
	}

	for (auto n : numberOfAtoms_) {
		this->numberOfAtoms.push_back(n);
	}
}


Substance::Substance(const Substance& sub) {
	this->atoms.clear();
	this->numberOfAtoms.clear();


	for (unsigned int i = 0; i < sub.getAtoms().size(); i++) {
		this->atoms.push_back(sub.getAtoms().at(i));		
		this->numberOfAtoms.push_back(sub.getNumberOfAtoms().at(i));

	}

}



Substance::~Substance() {
	this->atoms.clear();
	this->numberOfAtoms.clear();
}


std::vector<Atom> Substance::getAtoms() const {
	return this->atoms;
}


std::vector<int> Substance::getNumberOfAtoms() const {
	return this->numberOfAtoms;
}


void Substance::addAtoms(std::vector<Atom> atoms_) {
	for (auto a : atoms_) {
		atoms.push_back(a);
	}
}


void Substance::addAtoms(Atom atom_) {
	atoms.push_back(atom_);
}


void Substance::addNumberOfAtoms(std::vector<int> num_) {
	for (auto n : num_) {
		numberOfAtoms.push_back(n);
	}
}


void Substance::addNumberOfAtoms(int num_) {
	numberOfAtoms.push_back(num_);
}


bool operator== (const Substance& lSub, const Substance& rSub) {
	unsigned int comp_sub = 0;

	if (lSub.getAtoms().size() == rSub.getAtoms().size()) {
		for (unsigned int i = 0; i < lSub.getAtoms().size(); i++) {
			if ((lSub.getAtoms().at(i).getElementSymbol() == rSub.getAtoms().at(i).getElementSymbol())
				&& (lSub.getNumberOfAtoms().at(i) == rSub.getNumberOfAtoms().at(i))) {
				comp_sub += 1;
			}
		}

		if (comp_sub == lSub.getAtoms().size()) {
			return true;
		}
		else {
			return false;
		}
	}
	else {
		return false;
	}

}


Substance& Substance::operator=(const Substance& sub) {
	this->atoms.clear();
	this->numberOfAtoms.clear();

	for (unsigned int i = 0; i < sub.getAtoms().size(); i++) {
		this->atoms.push_back(sub.getAtoms().at(i));
		this->numberOfAtoms.push_back(sub.getNumberOfAtoms().at(i));
	}

	return *this;
}
