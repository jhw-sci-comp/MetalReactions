#ifndef METALHYDROCHLORICACIDREACTION_H
#define METALHYDROCHLORICACIDREACTION_H

#include <numeric>
#include "chemicalreaction.h"
#include "elementsfileo.h"

class MetalHydrochloricAcidReaction : public ChemicalReaction {
public:
	MetalHydrochloricAcidReaction();
	MetalHydrochloricAcidReaction(std::vector<Substance>);
	~MetalHydrochloricAcidReaction();
	void calculateProducts();
	std::vector<int> getIonicCharges();


private:
	std::vector<int> ionicCharges;
};

#endif
