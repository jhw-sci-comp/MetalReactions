#include "atom.h"


Atom::Atom(): elementSymbol{ "" }, elementName{ "" }, valenceElectrons{ 0 }
{}


Atom::Atom(std::string elementSymbol_, std::string elementName_, int valenceElectrons_) :
	elementSymbol{ elementSymbol_ }, elementName{ elementName_ }, valenceElectrons{ valenceElectrons_ }
{}


Atom::Atom(const Atom& atom) {
	this->elementSymbol = atom.getElementSymbol();
	this->elementName = atom.getElementName();
	this->valenceElectrons = atom.getValenceElectrons();
	this->electronegativity = atom.getElectronegativity();
	this->type = atom.getType();
	this->atomicNumber = atom.getAtomicNumber();
}



Atom::~Atom() {
	this->elementSymbol = "";
	this->elementName = "";
	this->valenceElectrons = 0;
	this->electronegativity = 0.0;
	this->atomicNumber = 0;
}


std::string Atom::getElementName() const {
	return this->elementName;
}


std::string Atom::getElementSymbol() const {
	return this->elementSymbol;
}


int Atom::getValenceElectrons() const {
	return this->valenceElectrons;
}


std::string Atom::getType() const {
	return type;
}


float Atom::getElectronegativity() const {
	return this->electronegativity;
}


unsigned int Atom::getAtomicNumber() const {
	return this->atomicNumber;
}


void Atom::setElementName(std::string elementName_) {
	this->elementName.assign(elementName_);

}


void Atom::setElementSymbol(std::string elementSymbol_) {
	this->elementSymbol.assign(elementSymbol_);
}


void Atom::setValenceElectrons(int valenceElectrons_) {
	this->valenceElectrons = valenceElectrons_;
}


void Atom::setType(std::string type_) {
	this->type = type_;
}


void Atom::setAtomicNumber(unsigned int atomicNumber_) {
	this->atomicNumber = atomicNumber_;
}


Atom& Atom::operator=(const Atom& atom) {
	this->elementSymbol = atom.getElementSymbol();
	this->elementName = atom.getElementName();
	this->valenceElectrons = atom.getValenceElectrons();
	this->electronegativity = atom.getElectronegativity();
	this->type = atom.getType();
	this->atomicNumber = atom.getAtomicNumber();

	return *this;
}