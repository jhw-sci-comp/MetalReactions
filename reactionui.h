#ifndef REACTIONUI_H
#define REACTIONUI_H

#include <iostream>
#include "elementsfileo.h"
#include "metaloxygenreaction.h"
#include "metalhydrochloricacidreaction.h"
#include "substance.h"
#include "atomicnumberinputexception.h"
#include "reactioninputexception.h"


using namespace std;

class ReactionUI {
public:
	ReactionUI();
	~ReactionUI();
	void displaySelection();
	void displayChemEq();

private:
	ElementsFileO elementsFileO;
	ChemicalReaction* pReaction = nullptr;
};

#endif
