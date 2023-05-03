#include <iostream>
#include "metaloxygenreaction.h"
#include "reactionui.h"

using namespace std;

int main()
{
	try {
		ReactionUI ui;
		ui.displaySelection();
		ui.displayChemEq();
	}
	catch (exception e) {
		cout << "Program terminated with exception." << endl;
	}

    return 0;
}
