#include "reactionui.h"

ReactionUI::ReactionUI() {}


ReactionUI::~ReactionUI() {
	delete pReaction;
}


void ReactionUI::displaySelection() {
	std::string inputStr;
	std::list<string> reactants_strings;
	std::list<string> atomic_numbers_strings;
	std::vector<Substance> reactants;
	int inputChoice;
	Atom atom;
	Substance sub;
	std::list<string>::iterator it;
	AtomicNumberInputException atomic_number_except;
	ReactionInputException rie;

	try {
		elementsFileO.selectListOfAtoms("metal");

		cout << "Select metal from list by entering the atomic number (Z):" << endl << endl;

		cout << "Z\t" << "Element name" << endl;
		cout << "--------------------" << endl;

		for (auto a : elementsFileO.getSelectedAtoms()) {
			cout << a.getAtomicNumber() << "\t" << a.getElementName() << " (" << a.getElementSymbol() << ")" << endl;
			atomic_numbers_strings.push_back(to_string(a.getAtomicNumber()));
		}
		cout << endl;

		elementsFileO.deleteSelectedAtoms();
	

		cout << ">";
		cin >> inputStr;

		it = find(atomic_numbers_strings.begin(), atomic_numbers_strings.end(), inputStr);
		if(it == atomic_numbers_strings.end()) {
			throw atomic_number_except;
		}

		reactants_strings.push_back(inputStr);
		elementsFileO.selectListOfAtoms(reactants_strings);
		
		sub.addAtoms(elementsFileO.getSelectedAtoms());
		sub.addNumberOfAtoms(1);
		
		reactants.push_back(sub);


		cout << endl << endl;
		cout << "Select substance from list from list by entering the number:" << endl << endl;

		cout << "No.\t" << "Substance" << endl;
		cout << "--------------------" << endl;
		cout << "1\t" << "Oxygen" << endl;
		cout << "2\t" << "HCl (Hydrochloric acid)" << endl;
		cout << endl;

		cout << ">";
		cin >> inputChoice;

		if ((inputChoice != 1) && (inputChoice != 2)) {
			throw rie;
		}


		switch (inputChoice) {
		case 1:
			pReaction = new MetalOxygenReaction(reactants);
			pReaction->calculateProducts();
			pReaction->balanceEquation();
			break;

		case 2:
			pReaction = new MetalHydrochloricAcidReaction(reactants);
			pReaction->calculateProducts();
			pReaction->balanceEquation();
			break;
		}

	}
	catch (ElementsFileException e) {
		cout << "ERROR: Could not find file 'elements.dat'." << endl;
		throw e;
	}
	catch (AtomicNumberInputException ae) {
		cout << "ERROR: Invalid input for atomic number." << endl;
		throw ae;
	}
	catch (ReactionInputException re) {
		cout << "ERROR: Invalid input for reaction choice." << endl;
		throw re;
	}

}


void ReactionUI::displayChemEq() {
	unsigned int x_index = 0;

	cout << endl;

	for (unsigned int i = 0; i < pReaction->getReactants().size(); i++) {
		if (pReaction->getX().at(x_index) != 1) {
			cout << pReaction->getX().at(x_index) << " ";
		}

		for (unsigned int j = 0; j < pReaction->getReactants().at(i).getAtoms().size(); j++) {
			cout << pReaction->getReactants().at(i).getAtoms().at(j).getElementSymbol();
			if (pReaction->getReactants().at(i).getNumberOfAtoms().at(j) != 1) {
				cout << "_{" << pReaction->getReactants().at(i).getNumberOfAtoms().at(j) << "}";
			}
		}

		if (i != pReaction->getReactants().size() - 1) {
			cout << " + ";
		}
		else {
			cout << " -> ";
		}

		x_index++;
	}

	for (unsigned int i = 0; i < pReaction->getProducts().size(); i++) {
		if ((x_index < pReaction->getX().size()) && (pReaction->getX().at(x_index) != 1)) {
			cout << pReaction->getX().at(x_index) << " ";
		}
		else if ((x_index >= pReaction->getX().size()) && (pReaction->getX0() != 1)) {
			cout << pReaction->getX0() << " ";
		}

		for (unsigned int j = 0; j < pReaction->getProducts().at(i).getAtoms().size(); j++) {
			cout << pReaction->getProducts().at(i).getAtoms().at(j).getElementSymbol();
			if (pReaction->getProducts().at(i).getNumberOfAtoms().at(j) != 1) {
				cout << "_{" << pReaction->getProducts().at(i).getNumberOfAtoms().at(j) << "}";
			}
		}

		if (i != pReaction->getProducts().size() - 1) {
			cout << " + ";
		}

		x_index++;
	}
	cout << endl;
}


