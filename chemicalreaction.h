#ifndef CHEMICALREACTION_H
#define CHEMICALREACTION_H


#include <math.h>
#include "substance.h"
#include "Matrix.t"
#include "balancedbipartitegraph.h"


class ChemicalReaction {
public:
	ChemicalReaction();
	virtual ~ChemicalReaction();
	virtual void calculateProducts() = 0;
	void balanceEquation();
	std::vector<Substance> getReactants();
	std::vector<Substance> getProducts();
	void addReactants(std::vector<Substance>);
	void addReactants(Substance);
	void addProducts(std::vector<Substance>);
	void addProducts(Substance);
	std::vector<double> getX() const;
	double getX0() const;

protected:
	std::vector<Substance> reactants;
	std::vector<Substance> products;
	std::vector<double> x;
	double x_0 = 0.0;


private:
	void findCombinations(unsigned int, unsigned int, std::vector<Substance>&, std::vector<std::vector<Substance>>&, unsigned int&, unsigned int&);
	std::vector<double> solveSLE_ICRB(const Matrix<int>&, std::vector<int>, double&);
};


#endif
