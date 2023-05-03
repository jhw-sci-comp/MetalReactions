#include "chemicalreaction.h"


ChemicalReaction::ChemicalReaction() {}

ChemicalReaction::~ChemicalReaction() {
	reactants.clear();
	products.clear();
}



void ChemicalReaction::balanceEquation() {
	std::vector<std::string> elemSymbSet;
	Matrix<int>* reactionCoefficientMatrix;
	std::vector<std::string>::iterator iter;
	std::vector<Substance>::iterator iter_sub;
	unsigned int index;

	unsigned int indexCombi;

	std::vector<Substance> productVectorSubs;
	std::vector<int>* productAtomsVector;

	unsigned int start = 0;
	unsigned int lengthCombi;
	std::vector<Substance> combination; 
	std::vector<std::vector<Substance>> combinationsVec;
	unsigned int fixedElements = 0;
	unsigned int level = 0;

	int det_reactionCoefficientMatrix = 0;

	std::vector<std::pair<int, int>> bpg_edges;
	std::vector<std::pair<int, int>> edges_temp;
	BalancedBipartiteGraph* pBBPG;
	std::vector<int> rowRearrangement;

	for (auto s : reactants) {
		for (auto a : s.getAtoms()) {
			iter = std::find(elemSymbSet.begin(), elemSymbSet.end(), a.getElementSymbol());

			if (iter == elemSymbSet.end()) {
				elemSymbSet.push_back(a.getElementSymbol());
			}
		}
	}

	for (auto s : products) {
		for (auto a : s.getAtoms()) {
			iter = std::find(elemSymbSet.begin(), elemSymbSet.end(), a.getElementSymbol());

			if (iter == elemSymbSet.end()) {
				elemSymbSet.push_back(a.getElementSymbol());
			}
		}
	}

	reactionCoefficientMatrix = new Matrix<int>(elemSymbSet.size(), elemSymbSet.size());

	productAtomsVector = new std::vector<int>(elemSymbSet.size(), 0);

	// assemble matrix
	for (unsigned int k = 0; k < reactants.size(); k++) {
		for (unsigned int l = 0; l < reactants.at(k).getAtoms().size(); l++) {
			iter = std::find(elemSymbSet.begin(), elemSymbSet.end(), reactants.at(k).getAtoms().at(l).getElementSymbol());

			index = std::distance(elemSymbSet.begin(), iter);

			reactionCoefficientMatrix->setEntry(index, k, reactants.at(k).getNumberOfAtoms().at(l));

			bpg_edges.push_back(std::make_pair(index, k));
		}
	}



	if (reactants.size() < elemSymbSet.size()) {		
		lengthCombi = elemSymbSet.size() - reactants.size();
		
		findCombinations(start, lengthCombi, combination, combinationsVec, fixedElements, level);
		
		
		for (unsigned int i = 0; (i < combinationsVec.size()) && (det_reactionCoefficientMatrix == 0); i++) {
			for (unsigned int l = 0; l < lengthCombi; l++) {
				reactionCoefficientMatrix->setCol(elemSymbSet.size() - lengthCombi + l, 0);
			}

			edges_temp.clear();

			for (unsigned int j = 0; j < combinationsVec.at(i).size(); j++) {
				for (unsigned int k = 0; k < combinationsVec.at(i).at(j).getAtoms().size(); k++) {

					iter = std::find(elemSymbSet.begin(), elemSymbSet.end(), combinationsVec.at(i).at(j).getAtoms().at(k).getElementSymbol());

					index = std::distance(elemSymbSet.begin(), iter);

					reactionCoefficientMatrix->setEntry(index, elemSymbSet.size() - combinationsVec.at(i).size() + j, -combinationsVec.at(i).at(j).getNumberOfAtoms().at(k));

					edges_temp.push_back(std::make_pair(index, elemSymbSet.size() - combinationsVec.at(i).size() + j));
				}

			}
			
			det_reactionCoefficientMatrix = det(*reactionCoefficientMatrix);
			indexCombi = i;
			
			if (det_reactionCoefficientMatrix != 0) {
				for (auto p : edges_temp) {
					bpg_edges.push_back(p);
				}

				combination = combinationsVec.at(i);
			}
		}

		
		//determine remaining substances for left side
		for (auto s : products) {
			iter_sub = std::find(combinationsVec.at(indexCombi).begin(), combinationsVec.at(indexCombi).end(), s);

			if (iter_sub == combinationsVec.at(indexCombi).end()) {
				productVectorSubs.push_back(s);
			}
		}


		//assemble right side of system
		for (auto s : productVectorSubs) {
			for (unsigned int i = 0; i < s.getAtoms().size(); i++) {
				iter = std::find(elemSymbSet.begin(), elemSymbSet.end(), s.getAtoms().at(i).getElementSymbol());

				index = std::distance(elemSymbSet.begin(), iter);

				productAtomsVector->at(index) += s.getNumberOfAtoms().at(i);
			}
		}

	}
	else {
		//assemble right side of system		
		for (auto s : products) {
			for (unsigned int i = 0; i < s.getAtoms().size(); i++) {
				iter = std::find(elemSymbSet.begin(), elemSymbSet.end(), s.getAtoms().at(i).getElementSymbol());

				index = std::distance(elemSymbSet.begin(), iter);

				productAtomsVector->at(index) += s.getNumberOfAtoms().at(i);
			}
		}
	}

	//create bipartite graph from matrix and determine maximum matching which is a perfect matching
	pBBPG = new BalancedBipartiteGraph(elemSymbSet.size(), bpg_edges);
	rowRearrangement = pBBPG->findMaxMatch();

	//rearrange matrix according to bipartite graph
	std::vector<unsigned int> newIndex;
	std::vector<unsigned int>::iterator it_newIndex;
	unsigned int product_temp;

	for (unsigned int i = 0; i < elemSymbSet.size(); i++) {
		newIndex.push_back(i);
	}



	for (unsigned int k = 0; k < rowRearrangement.size(); k++) {
		reactionCoefficientMatrix->swapRows(newIndex.at(k), rowRearrangement.at(k));

		product_temp = productAtomsVector->at(newIndex.at(k));
		productAtomsVector->at(newIndex.at(k)) = productAtomsVector->at(rowRearrangement.at(k));
		productAtomsVector->at(rowRearrangement.at(k)) = product_temp;

		it_newIndex = std::find(newIndex.begin(), newIndex.end(), rowRearrangement.at(k));
		index = std::distance(newIndex.begin(), it_newIndex);
		newIndex.at(k) = rowRearrangement.at(k);
		newIndex.at(index) = k;
	}


	//solving system of linear equations
	x = solveSLE_ICRB(*reactionCoefficientMatrix, *productAtomsVector, x_0);

	//rearrange entries of product vector
	std::vector<Substance> products_temp(products);

	for (auto s : combination) {
		iter_sub = std::find(products_temp.begin(), products_temp.end(), s);
		products_temp.erase(iter_sub);
	}

	products.clear();

	for (auto s : combination) {
		products.push_back(s);
	}

	for (auto s : products_temp) {
		products.push_back(s);
	}

	delete pBBPG;
	delete reactionCoefficientMatrix;
	delete productAtomsVector;
}





std::vector<Substance> ChemicalReaction::getReactants() {
	return reactants;
}


std::vector<Substance> ChemicalReaction::getProducts() {
	return products;
}


void ChemicalReaction::addReactants(std::vector<Substance> reactants_) {
	for (auto r : reactants_) {
		reactants.push_back(r);
	}
}


void ChemicalReaction::addReactants(Substance reactant_) {
	reactants.push_back(reactant_);
}


void ChemicalReaction::addProducts(std::vector<Substance> products_) {
	for (auto p : products_) {
		products.push_back(p);
	}
}


void ChemicalReaction::addProducts(Substance product_) {
	products.push_back(product_);
}



//find all possible combinaations of product substance of length lengthCombi
void ChemicalReaction::findCombinations(unsigned int startElement, unsigned int lengthCombi, std::vector<Substance>& combination, std::vector<std::vector<Substance>>& combinationsVec, unsigned int& fixedElements, unsigned int& level) {

	for (unsigned int i = startElement; i < products.size() - lengthCombi + level + 1; i++) {
		combination.push_back(products.at(i));
		fixedElements++;

		if ((lengthCombi - fixedElements) > 0) {
			level++;
			findCombinations(i + 1, lengthCombi, combination, combinationsVec, fixedElements, level);
		}
		else if (i < products.size() - 1) {
			combinationsVec.push_back(combination);			
			combination.pop_back();
			fixedElements--;
		}
		else {
			combinationsVec.push_back(combination);								
		}
	}
	for (unsigned int k = 0; k < level; k++) {
		combination.pop_back();
		fixedElements--;
	}

	level--;

}



std::vector<double> ChemicalReaction::solveSLE_ICRB(const Matrix<int>& rcMatrix, std::vector<int> prodVector, double& x_0_var) {
	unsigned int n = prodVector.size();
	std::vector<double> x_vec(n, 0);
	std::vector<double> r_vec(n, 0);
	double s = 0.0;
	int stop = 0;
	unsigned int i_max = 100;
	double tol = 0.0001;
	double r = 0.0;



	for (unsigned int k = 0; (k < i_max) && (stop == 0); k++) {
		for (unsigned int i = 0; i < n; i++) {
			s = 0.0;
			r_vec = x_vec;

			for (unsigned int j = 0; j < n; j++) {
				if (j != i) {
					s = s + static_cast<double>(rcMatrix.getEntry(i, j) * x_vec.at(j));
				}
			}
			x_vec.at(i) = (static_cast<double>(prodVector.at(i)) - s) / static_cast<double>(rcMatrix.getEntry(i, i));

			if (std::fabs((r_vec.at(i) - x_vec.at(i)) / r_vec.at(i)) <= tol) {
				stop = 1;
			}
		}
	}

	
	r_vec = x_vec;
	stop = 0;

	s = 0.0;
	for (auto e : r_vec) {
		s += e;
	}

	for (unsigned int k = 1; (k < i_max) && (stop == 0); k++) {
		x_0_var = static_cast<double>(k);
		r = std::fabs((k * s - std::round(k * s)) / s);

		if (r <= tol) {
			stop = 1;
		}
	}	

	for (unsigned int i = 0; i < x_vec.size(); i++) {
		x_vec.at(i) = x_0_var * x_vec.at(i);
	}

	return x_vec;
}


std::vector<double> ChemicalReaction::getX() const {
	return x;
}


double ChemicalReaction::getX0() const {
	return x_0;
}
