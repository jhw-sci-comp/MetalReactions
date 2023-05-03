#include "elementsfileo.h"

ElementsFileO::ElementsFileO() {}


bool ElementsFileO::selectAtom(unsigned int atomicNumber, Atom& selectedAtom) {
	std::string str;
	std::fstream elemFile;
	bool elemFound = false;
	std::size_t separator_pos, start_pos = 0;
	std::vector<std::string> substrings;
	
	std::exception e;
	
	elemFile.open("elements.dat");
	if (!elemFile) {		
		throw e;
	}
	

	while (getline(elemFile, str) && (elemFound == false)) {

		separator_pos = str.find_first_of(",", start_pos);
		while (separator_pos != std::string::npos) {
			substrings.push_back(str.substr(start_pos, separator_pos - start_pos));

			start_pos = separator_pos + 1;
			separator_pos = str.find_first_of(",", start_pos);

		}
		substrings.push_back(str.substr(start_pos, separator_pos - start_pos));

		start_pos = 0;


		if (std::to_string(atomicNumber).compare(substrings[0]) == 0) {
			elemFound = true;

			selectedAtom.setAtomicNumber(std::stoi(substrings[0]));
			selectedAtom.setElementSymbol(substrings[1]);
			selectedAtom.setElementName(substrings[2]);
			selectedAtom.setValenceElectrons(std::stoi(substrings[3]));
			selectedAtom.setType(substrings[4]);
		}

		substrings.clear();
	}

	elemFile.close();
	

	return elemFound;
}


bool ElementsFileO::selectAtom(std::string elemStr, Atom& selectedAtom) {
	std::string str, elemNameStr, elemStr_temp;
	std::fstream elemFile;
	bool elemFound = false;
	std::size_t separator_pos, start_pos = 0;
	std::vector<std::string> substrings;

	
	elemFile.open("elements.dat");

	while (getline(elemFile, str) && (elemFound == false)) {

		separator_pos = str.find_first_of(",", start_pos);
		while (separator_pos != std::string::npos) {
			substrings.push_back(str.substr(start_pos, separator_pos - start_pos));

			start_pos = separator_pos + 1;
			separator_pos = str.find_first_of(",", start_pos);
		}
		substrings.push_back(str.substr(start_pos, separator_pos - start_pos));

		start_pos = 0;

		elemNameStr.assign(substrings[2]);
		std::transform(elemNameStr.begin(), elemNameStr.end(), elemNameStr.begin(), ::tolower);
		elemStr_temp.assign(elemStr);
		std::transform(elemStr_temp.begin(), elemStr_temp.end(), elemStr_temp.begin(), ::tolower);


		if ((elemStr.compare(substrings[1]) == 0) || (elemStr_temp.compare(elemNameStr) == 0)) {
			elemFound = true;

			selectedAtom.setAtomicNumber(std::stoi(substrings[0]));
			selectedAtom.setElementSymbol(substrings[1]);
			selectedAtom.setElementName(substrings[2]);
			selectedAtom.setValenceElectrons(std::stoi(substrings[3]));
			selectedAtom.setType(substrings[4]);
		}

		substrings.clear();
	}

		elemFile.close();
	

	return elemFound;
}



void ElementsFileO::selectListOfAtoms(std::list<std::string> inputList) {
	std::string str;
	bool isInteger = true;
	Atom atom_temp;


	for (std::list<std::string>::iterator it = inputList.begin(); it != inputList.end(); it++) {		
		str = *it;

		for (unsigned int i = 0; (i < str.size()) && (isInteger == true); i++) {
			if (isdigit(str[i]) == false) {
				isInteger = false;
			}
		}

		if (isInteger) {
			selectAtom(std::stoi(str), atom_temp);
		}
		else {
			selectAtom(str, atom_temp);
		}

		selectedAtoms.push_back(atom_temp);

		isInteger = true;
	}
	
}


void ElementsFileO::selectListOfAtoms(std::string selectType) {
	std::string str;
	std::fstream elemFile;
	std::size_t separator_pos, start_pos = 0;
	std::vector<std::string> substrings;
	Atom atom_temp;

	ElementsFileException e;
	
			
	elemFile.open("elements.dat");
	if (!elemFile) {
		throw e;
	}

	while (getline(elemFile, str)) {
		separator_pos = str.find_first_of(",", start_pos);
		while (separator_pos != std::string::npos) {
			substrings.push_back(str.substr(start_pos, separator_pos - start_pos));

			start_pos = separator_pos + 1;
			separator_pos = str.find_first_of(",", start_pos);
		}
		substrings.push_back(str.substr(start_pos, separator_pos - start_pos));

		start_pos = 0;

		if (selectType.compare(substrings.at(4)) == 0) {

			atom_temp.setAtomicNumber(std::stoi(substrings[0]));
			atom_temp.setElementSymbol(substrings[1]);
			atom_temp.setElementName(substrings[2]);
			atom_temp.setValenceElectrons(std::stoi(substrings[3]));
			atom_temp.setType(substrings[4]);

			selectedAtoms.push_back(atom_temp);
		}

		substrings.clear();
	}

	elemFile.close();
	
}


std::vector<Atom> ElementsFileO::getSelectedAtoms() const {
	return selectedAtoms;
}


void ElementsFileO::deleteSelectedAtoms() {
	selectedAtoms.clear();
}