#include "metaloxygenreaction.h"

#include <iostream>


MetalOxygenReaction::MetalOxygenReaction() {}


MetalOxygenReaction::MetalOxygenReaction(std::vector<Substance> reactants_) {
	Atom oxygen;
	Substance o2;
	ElementsFileO elementsFileO;

	for (auto r : reactants_) {
		this->reactants.push_back(r);
	}

	elementsFileO.selectAtom("Oxygen", oxygen);
	o2.addAtoms(oxygen);
	o2.addNumberOfAtoms(2);

	reactants.push_back(o2);
}


MetalOxygenReaction::~MetalOxygenReaction() {
	ionicCharges.clear();
}


std::vector<int> MetalOxygenReaction::getIonicCharges() {
	return ionicCharges;
}


void MetalOxygenReaction::calculateProducts() {
	int lcm_valencies = 0, charge1, charge2;
	Substance sub_temp;


	if (this->reactants.size() == 2) {
		charge1 = reactants.at(0).getAtoms().at(0).getValenceElectrons() > 4 ? -(8 - reactants.at(0).getAtoms().at(0).getValenceElectrons()) : reactants.at(0).getAtoms().at(0).getValenceElectrons();
		charge2 = reactants.at(1).getAtoms().at(0).getValenceElectrons() > 4 ? -(8 - reactants.at(1).getAtoms().at(0).getValenceElectrons()) : reactants.at(1).getAtoms().at(0).getValenceElectrons();

		lcm_valencies = std::lcm(std::abs(charge1), std::abs(charge2));
		sub_temp.addAtoms(reactants.at(0).getAtoms().at(0));
		sub_temp.addNumberOfAtoms(lcm_valencies / std::abs(charge1));

		sub_temp.addAtoms(reactants.at(1).getAtoms().at(0));
		sub_temp.addNumberOfAtoms(lcm_valencies / std::abs(charge2));
		products.push_back(sub_temp);

		ionicCharges.push_back(charge1);
		ionicCharges.push_back(charge2);

	}
	else {
		//TODO: EXCEPTIONS!!!!!!!!!!!!!
	}
}


