#ifndef METALOXYGENREACTION_H
#define METALOXYGENREACTION_H

#include <numeric>
#include "chemicalreaction.h"
#include "elementsfileo.h"

class MetalOxygenReaction : public ChemicalReaction {
public:
	MetalOxygenReaction();
	MetalOxygenReaction(std::vector<Substance>);
	~MetalOxygenReaction();
	void calculateProducts();
	std::vector<int> getIonicCharges();


private:
	std::vector<int> ionicCharges;
};

#endif
