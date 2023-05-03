#include "metalhydrochloricacidreaction.h"


MetalHydrochloricAcidReaction::MetalHydrochloricAcidReaction() {}


MetalHydrochloricAcidReaction::MetalHydrochloricAcidReaction(std::vector<Substance> reactants_) {
	Atom hydrogen, chlorine;
	Substance hCl;
	ElementsFileO elementsFileO;

	for (auto r : reactants_) {
		this->reactants.push_back(r);
	}

	elementsFileO.selectAtom("Hydrogen", hydrogen);
	hCl.addAtoms(hydrogen);
	hCl.addNumberOfAtoms(1);
	elementsFileO.selectAtom("Chlorine", chlorine);
	hCl.addAtoms(chlorine);
	hCl.addNumberOfAtoms(1);

	reactants.push_back(hCl);
}


MetalHydrochloricAcidReaction::~MetalHydrochloricAcidReaction() {
	ionicCharges.clear();
}


std::vector<int> MetalHydrochloricAcidReaction::getIonicCharges() {
	return ionicCharges;
}


void MetalHydrochloricAcidReaction::calculateProducts() {
	int lcm_valencies = 0, charge1, charge2;
	Substance metalChloride, h2;
	Atom hydrogen;
	ElementsFileO elementsFileO;


	if (this->reactants.size() == 2) {
		charge1 = reactants.at(0).getAtoms().at(0).getValenceElectrons() > 4 ? -(8 - reactants.at(0).getAtoms().at(0).getValenceElectrons()) : reactants.at(0).getAtoms().at(0).getValenceElectrons();
		charge2 = reactants.at(1).getAtoms().at(1).getValenceElectrons() > 4 ? -(8 - reactants.at(1).getAtoms().at(1).getValenceElectrons()) : reactants.at(1).getAtoms().at(1).getValenceElectrons();

		lcm_valencies = std::lcm(std::abs(charge1), std::abs(charge2));
		metalChloride.addAtoms(reactants.at(0).getAtoms().at(0));
		metalChloride.addNumberOfAtoms(lcm_valencies / std::abs(charge1));

		metalChloride.addAtoms(reactants.at(1).getAtoms().at(1));
		metalChloride.addNumberOfAtoms(lcm_valencies / std::abs(charge2));
		products.push_back(metalChloride);

		ionicCharges.push_back(charge1);
		ionicCharges.push_back(charge2);

		elementsFileO.selectAtom("Hydrogen", hydrogen);
		h2.addAtoms(hydrogen);
		h2.addNumberOfAtoms(2);

		products.push_back(h2);

	}
	else {
		//TODO: EXCEPTIONS!!!!!!!!!!!!!
	}
}